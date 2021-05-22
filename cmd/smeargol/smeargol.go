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
	"compress/gzip"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/traverse"

	"github.com/kortschak/gogo"

	"github.com/kortschak/smeargol/internal/owl"
)

func main() {
	var (
		in       = flag.String("in", "", "specify the counts input (.tsv.gz - required)")
		ontopath = flag.String("ontology", "", "specify the GO file (.owl.gz - required)")
		mappath  = flag.String("map", "", "specify the ENSG to GO mapping (.nt.gz/.nq.gz - required)")
		lean     = flag.Bool("lean", true, "only load relevant parts of ontology")
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
	ontoData := distributeCounts(ontology, roots, data)

	log.Println("[printing smeared counts]")
	fmt.Printf("go_term\tgo_root\tgo_aspect\tdepth\t%s\n", strings.Join(data.names, "\t"))
	for i := range ontoData {
		walkDownSubClassesFrom(roots[i], ontology, func(r, t rdf.Term, d int) {
			counts, ok := ontoData[i][t.Value]
			if !ok {
				return
			}
			fmt.Printf("%s\t%s\t%s\t%d", t.Value, r.Value, nameSpaceOf(t, ontology), d)
			for _, v := range counts {
				fmt.Printf("\t%v", v)
			}
			fmt.Println()
		})
	}
}

// countData holds count data for a set of named samples each with a collection
// of features.
type countData struct {
	// names is the names of the samples.
	names []string

	// counts holds a set of counts for
	// features keyed by the map key.
	// The length of each []float64
	// must match the length of names
	// and indexing into the []float64
	// reflects indexing into names.
	counts map[string][]float64
}

// mappingCounts returns the count data held in the file at path.
func mappingCounts(path string) (*countData, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}

	c := csv.NewReader(r)
	c.Comma = '\t'
	c.Comment = '#'

	labels, err := c.Read()
	if err != nil {
		if err == io.EOF {
			return nil, io.ErrUnexpectedEOF
		}
		return nil, err
	}
	if labels[0] != "Geneid" {
		return nil, fmt.Errorf(`unexpected first column name: %q != "Geneid"`, labels[0])
	}
	samples := labels[1:]

	data := make(map[string][]float64, len(samples))

	c.ReuseRecord = true
	for {
		counts, err := c.Read()
		if err != nil {
			if err != io.EOF {
				return nil, err
			}
			break
		}
		geneid := counts[0]
		for i, f := range counts[1:] {
			v, err := strconv.ParseFloat(f, 64)
			if err != nil {
				return nil, fmt.Errorf("error parsing value for %q in sample %q: %v", geneid, samples[i], err)
			}
			data[geneid] = append(data[geneid], v)
		}
	}

	return &countData{
		names:  samples,
		counts: data,
	}, nil
}

// ontologyGraph returns the graph for the ontology stored in an OBO in OWL
// file. The namespaces are not expanded to full IRI namespaces.
func ontologyGraph(path string, lean bool) (*gogo.Graph, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}

	g := gogo.NewGraph()
	dec, err := owl.NewDecoder(r)
	if err != nil {
		return nil, err
	}
	for {
		s, err := dec.UnmarshalLocal()
		if err != nil {
			if err != io.EOF {
				return nil, err
			}
			break
		}

		if lean {
			// Filter on statements that are actually used. The filter
			// can be more restrictive since we only need "GO subclass GO"
			// and "GO hasOBONamespace literal" statements, but this is
			// good enough.
			switch s.Predicate.Value {
			// This list must include all predicates used in traversals
			// except <local:annotates> which comes from connectGeneIDsTo.
			case "<rdfs:subClassOf>", "<oboInOwl:hasOBONamespace>":
			default:
				continue
			}
		}
		g.AddStatement(s)
	}

	return g, nil
}

// connectGeneIDsTo adds the statements in path to the destination graph.
// The statements are expected to have local IRI namespaces and be in the
// following form:
//
//   <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .
//
// Only ENSG identifiers that match the names in the counts map are added
// to the graph.
func connectGeneIDsTo(dst *gogo.Graph, path string, counts map[string][]float64) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return err
	}

	dec := rdf.NewDecoder(r)
	for {
		s, err := dec.Unmarshal()
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}

		// Only keep annotations needed for the given counts.
		id := strings.TrimSuffix(strings.TrimPrefix(s.Object.Value, "<ensembl:"), ">")
		if _, ok := counts[id]; !ok {
			continue
		}

		s.Subject.UID = 0
		s.Predicate.UID = 0
		s.Object.UID = 0
		dst.AddStatement(s)
	}
}

// distributeCounts performs a breadth-first traversal from each of the leaf-most
// terms associated with each of the genes held by data, for each of the sample
// in data independently.
// Each gene ontology aspect is analysed separately since the aspects are not
// connected. The analyses are performed in parallel; the length of the
// returned slice will be the same as the number of roots passed in.
func distributeCounts(g *gogo.Graph, roots []rdf.Term, data *countData) []map[string][]float64 {
	ontoData := make([]map[string][]float64, len(roots))
	for i := range ontoData {
		ontoData[i] = make(map[string][]float64)
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
					if !dfs[i].Visited(l) {
						dst := ontoData[i][l.Value]
						if dst == nil {
							dst = make([]float64, len(counts))
							ontoData[i][l.Value] = dst
						}
						floats.Add(dst, counts)
					}
					dfs[i].Walk(g, l, func(n graph.Node) bool {
						t := n.(rdf.Term)
						dst := ontoData[i][t.Value]
						if dst == nil {
							dst = make([]float64, len(counts))
							ontoData[i][t.Value] = dst
						}
						floats.Add(dst, counts)
						return false
					})
				}
			}()
		}
		wg.Wait()
	}

	return ontoData
}

// isSubClassOfGO is a traverse edge filter. It accepts statements where
//
//  any -- <rdfs:subClassOf> -> <obo:GO_*
//
// for out queries from a term.
func isSubClassOfGO(e graph.Edge) bool {
	return gogo.ConnectedByAny(e, func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<rdfs:subClassOf>" &&
			strings.HasPrefix(s.Object.Value, "<obo:GO_")
	})
}

// leafiestFor return the leaf-most terms for gene from each of the ontology roots.
// The leaf sets are returned separated so that ontology count mutation can be
// performed concurrently without locking.
func leafiestFor(geneid string, g *gogo.Graph, roots []rdf.Term) [][]rdf.Term {
	leafiest := make([][]rdf.Term, len(roots))
	found := make([]bool, len(roots))
	var wg sync.WaitGroup
	for a, r := range roots {
		a := a
		r := r
		wg.Add(1)
		go func() {
			defer wg.Done()
			from, ok := g.TermFor("<ensembl:" + geneid + ">")
			if !ok {
				return
			}

			terms := g.Query(from).In(func(s *rdf.Statement) bool {
				return s.Predicate.Value == "<local:annotates>"
			}).Unique().Result()

			if len(terms) != 0 {
				found[a] = true
			}

			var depths []gogo.Descendant
			for _, q := range terms {
				ok, d := g.IsDescendantOf(r, q)
				if ok {
					depths = append(depths, gogo.Descendant{Term: q, Depth: d})
				}
			}
			sort.Sort(byDepth(depths))

			for i := 0; i < len(depths); i++ {
				a := depths[i]
				for j := i + 1; j < len(depths); {
					ok, _ := g.IsDescendantOf(a.Term, depths[j].Term)
					if ok {
						copy(depths[j:], depths[j+1:])
						depths = depths[:len(depths)-1]
					} else {
						j++
					}
				}
			}
			for _, d := range depths {
				leafiest[a] = append(leafiest[a], d.Term)
			}
		}()
	}
	wg.Wait()
	if !anyTrue(found) {
		log.Printf("no GO term found for %s", geneid)
	}
	return leafiest
}

func anyTrue(t []bool) bool {
	for _, v := range t {
		if v {
			return true
		}
	}
	return false
}

// byDepth sorts gogo.Descendents by depth, leafiest first.
type byDepth []gogo.Descendant

func (d byDepth) Len() int           { return len(d) }
func (d byDepth) Less(i, j int) bool { return d[i].Depth > d[j].Depth }
func (d byDepth) Swap(i, j int)      { d[i], d[j] = d[j], d[i] }

// nameSpaceOf returns the ontology aspect for the term t which is expected
// to be an <obo:GO_*> term.
func nameSpaceOf(t rdf.Term, in *gogo.Graph) string {
	ns := in.Query(t).Out(func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<oboInOwl:hasOBONamespace>"
	}).Result()
	switch len(ns) {
	case 0:
		return "NA"
	case 1:
		text, _, kind, err := ns[0].Parts()
		if err != nil {
			panic(fmt.Errorf("invalid term in graph: %w", err))
		}
		if kind == rdf.Literal {
			return text
		}
		return t.Value
	default:
		return t.Value
	}
}

// walkDownSubClassesFrom performs a breadth-first enumeration of GO subclass
// terms in g starting from r, and calling fn for each term, including r.
func walkDownSubClassesFrom(r rdf.Term, g *gogo.Graph, fn func(root, term rdf.Term, depth int)) {
	bf := traverse.BreadthFirst{Traverse: goIsSubClassOf}
	bf.Walk(reverse{g}, r, func(n graph.Node, d int) bool {
		fn(r, n.(rdf.Term), d)
		return false
	})
}

// goIsSubClassOf is a traverse edge filter. It accepts statements where
//
//  <obo:GO_* <- <rdfs:subClassOf> -- any
//
// for in queries from a term.
func goIsSubClassOf(e graph.Edge) bool {
	return gogo.ConnectedByAny(e, func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<rdfs:subClassOf>" &&
			strings.HasPrefix(s.Subject.Value, "<obo:GO_")
	})
}

// reverse implements the traverse.Graph reversing the direction of edges.
type reverse struct {
	*gogo.Graph
}

func (g reverse) From(id int64) graph.Nodes      { return g.Graph.To(id) }
func (g reverse) Edge(uid, vid int64) graph.Edge { return g.Graph.Edge(vid, uid) }
