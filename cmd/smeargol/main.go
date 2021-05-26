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
// written uncompressed to a matrix and a plot directory. A summary
// document is written to the specified out file in JSON format corresponding
// to the following Go structs.
//
//  type SummaryDoc struct {
//  	// Roots is the set of roots in the Gene Ontology.
//  	Roots []string
//
//  	// Summaries contains the summaries of a smeargol
//  	// analysis.
//  	Summaries [][]*Summary
//  }
//
//  type Summary struct {
//  	// Name is the name of the sample.
//  	Name string
//
//  	// Root is the root GO term for the summary.
//  	Root string
//
//  	// Depth is the distance from the root.
//  	Depth int
//
//  	// Rows and Cols are the dimensions of the matrix
//  	// describing the GO level. Rows corresponds to the
//  	// number of genes and Cols corresponds to the number
//  	// of GO terms in the level.
//  	Rows, Cols int
//
//  	// OptimalRank and FractionalRank are the calculated
//  	// ranks of the summary matrix. OptimalRank is
//  	// calculated according to the method of Matan Gavish
//  	// and David L. Donoho https://arxiv.org/abs/1305.5870.
//  	// FractionalRank is the rank calculated using the
//  	// user-provided fraction parameters.
//  	OptimalRank, FractionalRank int
//
//  	// Sigma is the complete set of singular values.
//  	Sigma []float64
//  }
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/mat"
)

func main() {
	var (
		in       = flag.String("in", "", "specify the counts input (.tsv.gz - required)")
		out      = flag.String("out", "", "specify the summary output file")
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
written uncompressed to a matrix and a plot directory. A summary
document is written to the specified out file in JSON format corresponding
to the following Go structs.

  type SummaryDoc struct {
  	// Roots is the set of roots in the Gene Ontology.
  	Roots []string

  	// Summaries contains the summaries of a smeargol
  	// analysis.
  	Summaries [][]*Summary
  }

  type Summary struct {
  	// Name is the name of the sample.
  	Name string

  	// Root is the root GO term for the summary.
  	Root string

  	// Depth is the distance from the root.
  	Depth int

  	// Rows and Cols are the dimensions of the matrix
  	// describing the GO level. Rows corresponds to the
  	// number of genes and Cols corresponds to the number
  	// of GO terms in the level.
  	Rows, Cols int

  	// OptimalRank and FractionalRank are the calculated
  	// ranks of the summary matrix. OptimalRank is
  	// calculated according to the method of Matan Gavish
  	// and David L. Donoho https://arxiv.org/abs/1305.5870.
  	// FractionalRank is the rank calculated using the
  	// user-provided fraction parameters.
  	OptimalRank, FractionalRank int

  	// Sigma is the complete set of singular values.
  	Sigma []float64
  }

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

	summaries := make([][]*Summary, len(ontoData))
	var wg sync.WaitGroup
	for k := range ontoData {
		k := k
		wg.Add(1)
		go func() {
			defer wg.Done()
			lastD := -1
			var goTerms []string
			walkDownSubClassesFrom(roots[k], ontology, func(r, t rdf.Term, d int) {
				dw.record(k, d, r, t)

				if lastD == -1 || d == lastD {
					goTerms = append(goTerms, t.Value)
					lastD = d
					return
				}
				lastD = d

				// Write out matrices for this depth. Note that d is now
				// referring to the next level.
				s, err := writeCountData(r.Value, d-1, goTerms, data, ontoData[k], *cut, *frac)
				if err != nil {
					log.Println(err)
				}
				summaries[k] = append(summaries[k], s...)
				goTerms = goTerms[:0]
				goTerms = append(goTerms, t.Value)
			})

			// Write out last depth.
			s, err := writeCountData(roots[k].Value, lastD, goTerms, data, ontoData[k], *cut, *frac)
			if err != nil {
				log.Println(err)
			}
			summaries[k] = append(summaries[k], s...)
			sort.Slice(summaries[k], func(i, j int) bool {
				s := summaries[k]
				switch {
				case s[i].Depth < s[j].Depth:
					return true
				case s[i].Depth > s[j].Depth:
					return false
				default:
					return s[i].Name < s[j].Name
				}
			})
		}()
	}
	wg.Wait()

	dw.flush()
	if *out != "" {
		rootNames := make([]string, len(roots))
		for i, r := range roots {
			rootNames[i] = "GO:" + strip(r.Value, "<obo:GO_", ">")
		}
		b, err := json.MarshalIndent(SummaryDoc{rootNames, summaries}, "", "\t")
		if err != nil {
			log.Fatal(err)
		}
		err = ioutil.WriteFile(*out, b, 0o644)
		if err != nil {
			log.Fatal(err)
		}
	}
}

// writeCountData writes out a matrix of gene expression data summed according
// to the bit vector data collected during the walk of the GO DAG. It also
// performs an SVD of the matrix, plotting the singular values and obtaining
// an optimal truncation for each GO level/aspect.
func writeCountData(root string, depth int, goTerms []string, data *countData, ontoData map[string]ontoCounts, cut, frac float64) ([]*Summary, error) {
	if len(goTerms) == 0 || len(data.geneIDs) == 0 {
		return nil, nil
	}
	root = strip(root, "<obo:", ">")

	sort.Strings(goTerms)
	m := mat.NewDense(len(data.geneIDs), len(goTerms), nil) // Assume all samples have same genes.
	var summaries []*Summary
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
		s, err := optimalTruncation(path, m, cut, frac)
		s.Name = name
		s.Root = root
		s.Depth = depth
		summaries = append(summaries, s)
		if err != nil {
			log.Println(err)
		}
		err = writeMatrix(path, data.geneIDs, goTerms, m)
		if err != nil {
			return summaries, err
		}

		m.Zero()
	}

	return summaries, nil
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
