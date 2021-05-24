// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// smeargol distributes count data across the Gene Ontology DAG provided and
// prints the GO terms, their roots and depths and distributed counts in a
// tsv table to stdout. It logs gene identifiers that do not have GO term
// annotations to stderr. The graph analysis assumes Ensembl gene identifiers
// and Gene Ontology graph structure.
//
// The input counts file is a tab-delimited file with the first column being
// Ensembl gene ID (ENSG00000000000) and remaining columns being count data.
// The first row is expected to be labelled with the first column being Geneid
// and the remaining columns holding the names of the samples.
//
// The Gene Ontology is required to be in Owl format. The file can be
// obtained from http://current.geneontology.org/ontology/go.owl.
//
// The ENSG to GO mapping is expected to be in RDF N-Triples or N-Quads in
// the form:
//
//  <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .
//
// for each GO term to Ensembl gene annotation.
//
// All input files are expected to be gzip compressed and the output is
// written uncompressed to standard output and standard error.
package main

import (
	"flag"
	"fmt"
	"log"
	"math/big"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/traverse"
	"gonum.org/v1/gonum/mat"

	"github.com/kortschak/gogo"
)

func main() {
	var (
		in       = flag.String("in", "", "specify the counts input (.tsv.gz - required)")
		ontopath = flag.String("ontology", "", "specify the GO file (.owl.gz - required)")
		mappath  = flag.String("map", "", "specify the ENSG to GO mapping (.nt.gz/.nq.gz - required)")
		lean     = flag.Bool("lean", true, "only load relevant parts of ontology")
		cut      = flag.Float64("cut", 1, "minimum valid singular value")
		frac     = flag.Float64("frac", 0.75, "include singular values up to this cumulative fraction")
		debug    = flag.Bool("debug", false, "output binary assignments - only small sets")
		help     = flag.Bool("help", false, "print help text")
	)
	flag.Parse()

	if *help {
		flag.Usage()
		fmt.Fprintf(os.Stderr, `
%s distributes count data across the gene ontology DAG provided and
prints the GO terms, their roots and depths and distributed counts in a
tsv table to stdout. It logs gene identifiers that do not have GO term
annotations to stderr. The graph analysis assumes Ensembl gene identifiers
and Gene Ontology graph structure.

The input counts file is a tab-delimited file with the first column being
Ensembl gene ID (ENSG00000000000) and remaining columns being count data.
The first row is expected to be labelled with the first column being Geneid
and the remaining columns holding the names of the samples.

The Gene Ontology is required to be in Owl format. The file can be
obtained from http://current.geneontology.org/ontology/go.owl.

The ENSG to GO mapping is expected to be in RDF N-Triples or N-Quads in
the form:

 <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .

for each GO term to Ensembl gene annotation.

All input files are expected to be gzip compressed and the output is
written uncompressed to standard output and standard error.

Copyright ©2020 Dan Kortschak. All rights reserved.

`, filepath.Base(os.Args[0]))
		os.Exit(0)
	}

	if *in == "" || *ontopath == "" || *mappath == "" {
		flag.Usage()
		os.Exit(2)
	}

	log.Println(os.Args)
	for _, d := range []string{
		"matrices",
		"plots",
	} {
		err := os.Mkdir(d, 0o755)
		if err != nil {
			log.Fatal(err)
		}
	}

	log.Println("[loading count data]")
	data, err := mappingCounts(*in)
	if err != nil {
		log.Fatalf("failed to load count data: %v", err)
	}

	if *lean {
		log.Println("[loading lean ontology]")
	} else {
		log.Println("[loading ontology]")
	}
	ontology, err := ontologyGraph(*ontopath, *lean)
	if err != nil {
		log.Fatalf("failed to load ontology: %v", err)
	}

	log.Println("[loading gene to ontology mappings]")
	err = connectGeneIDsTo(ontology, *mappath, data.counts)
	if err != nil {
		log.Fatalf("failed to connect gene IDs to ontology: %v", err)
	}

	log.Println("[smearing counts]")
	roots := ontology.Roots(false)
	sort.Slice(roots, func(i, j int) bool { return roots[i].Value < roots[j].Value })
	ontoData := distributeCounts(ontology, roots, data)

	log.Println("[writing smeared count matrices]")

	var dw *debugWriter
	if *debug {
		dw = newDebugWriter(ontology, ontoData, data)
	}

	var wg sync.WaitGroup
	for i := range ontoData {
		i := i
		wg.Add(1)
		go func() {
			defer wg.Done()
			lastD := -1
			var goTerms []string
			walkDownSubClassesFrom(roots[i], ontology, func(r, t rdf.Term, d int) {
				dw.record(i, d, r, t)

				if lastD == -1 || d == lastD {
					goTerms = append(goTerms, t.Value)
					lastD = d
					return
				}
				lastD = d

				// Write out matrices for this depth. Note that d is now
				// referring to the next level.
				writeCountData(r.Value, d-1, goTerms, data, ontoData[i], *cut, *frac)
				goTerms = goTerms[:0]
				goTerms = append(goTerms, t.Value)
			})

			// Write out last depth.
			writeCountData(roots[i].Value, lastD, goTerms, data, ontoData[i], *cut, *frac)
		}()
	}
	wg.Wait()
	dw.flush()
}

// writeCountData writes out a matrix of gene expression data summed according
// to the bit vector data collected during the walk of the GO DAG. It also
// performs an SVD of the matrix, plotting the singular values and obtaining
// an optimal truncation for each GO level/aspect.
func writeCountData(root string, depth int, goTerms []string, data *countData, ontoData map[string]ontoCounts, cut, frac float64) error {
	if len(goTerms) == 0 || len(data.geneIDs) == 0 {
		return nil
	}
	root = strip(root, "<obo:", ">")

	sort.Strings(goTerms)
	m := mat.NewDense(len(data.geneIDs), len(goTerms), nil) // Assume all samples have same genes.
	for sample, name := range data.names {
		for col, term := range goTerms {
			counts, ok := ontoData[term]
			if !ok {
				continue
			}
			for row, geneID := range data.geneIDs {
				if counts.vector[sample].Bit(row) == 0 {
					continue
				}
				m.Set(row, col, data.counts[geneID][sample])
			}
		}

		path := fmt.Sprintf("%s_%s_%03d", name, root, depth)
		err := optimalTruncation(path, m, cut, frac)
		if err != nil {
			log.Println(err)
		}
		err = writeMatrix(path, data.geneIDs, goTerms, m)
		if err != nil {
			return err
		}

		m.Zero()
	}

	return nil
}

func writeMatrix(path string, rows, cols []string, data *mat.Dense) (err error) {
	f, err := os.Create(filepath.Join("matrices", path+".tsv"))
	if err != nil {
		return err
	}
	defer func() {
		err = f.Close()
	}()

	_, err = f.Write([]byte{'\t'})
	if err != nil {
		return err
	}
	_, err = f.WriteString(strings.Join(stripSlice(cols, "<obo:", ">"), "\t"))
	if err != nil {
		return err
	}
	_, err = f.Write([]byte{'\n'})
	if err != nil {
		return err
	}

	for r, id := range rows {
		_, err = f.WriteString(id)
		if err != nil {
			return err
		}
		for c := range cols {
			_, err = fmt.Fprintf(f, "\t%v", data.At(r, c))
			if err != nil {
				return err
			}
		}
		_, err = f.Write([]byte{'\n'})
		if err != nil {
			return err
		}
	}
	return nil
}

func stripSlice(s []string, prefix, suffix string) []string {
	n := make([]string, len(s))
	for i, e := range s {
		n[i] = strip(e, prefix, suffix)
	}
	return n
}

func strip(s, prefix, suffix string) string {
	return strings.TrimSuffix(strings.TrimPrefix(s, prefix), suffix)
}

type ontoCounts struct {
	vector []big.Int
}

// distributeCounts performs a breadth-first traversal from each of the leaf-most
// terms associated with each of the genes held by data, for each of the sample
// in data independently.
// Each gene ontology aspect is analysed separately since the aspects are not
// connected. The analyses are performed in parallel; the length of the
// returned slice will be the same as the number of roots passed in.
func distributeCounts(g *gogo.Graph, roots []rdf.Term, data *countData) []map[string]ontoCounts {
	ontoData := make([]map[string]ontoCounts, len(roots))
	for i := range ontoData {
		ontoData[i] = make(map[string]ontoCounts)
	}

	dfs := make([]traverse.DepthFirst, len(roots))
	for i := range dfs {
		dfs[i] = traverse.DepthFirst{Traverse: isSubClassOfGO}
	}
	for geneid, counts := range data.counts {
		var wg sync.WaitGroup
		for i, aspect := range leafiestFor(geneid, g, roots) {
			i := i
			aspect := aspect
			wg.Add(1)
			go func() {
				defer wg.Done()
				dfs[i].Reset()
				for _, l := range aspect {
					updateOntoData(ontoData[i], l, geneid, counts, data)
					dfs[i].Walk(g, l, func(n graph.Node) bool {
						t := n.(rdf.Term)
						updateOntoData(ontoData[i], t, geneid, counts, data)
						return false
					})
				}
			}()
		}
		wg.Wait()
	}

	return ontoData
}

func updateOntoData(ontoData map[string]ontoCounts, t rdf.Term, geneid string, counts []float64, data *countData) {
	dst, ok := ontoData[t.Value]
	if !ok {
		vector := make([]big.Int, len(counts))
		for _, v := range vector {
			v.SetBit(&v, len(data.geneIdx), 0)
		}
		dst.vector = vector
		ontoData[t.Value] = dst
	}
	for j, c := range counts {
		if c == 0 {
			continue
		}
		dst.vector[j].SetBit(&dst.vector[j], data.geneIdx[geneid], 1)
	}
}
