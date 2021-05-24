// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strings"
	"sync"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/traverse"

	"github.com/kortschak/gogo"
	"github.com/kortschak/smeargol/internal/owl"
)

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
		id := strip(s.Object.Value, "<ensembl:", ">")
		if _, ok := counts[id]; !ok {
			continue
		}

		s.Subject.UID = 0
		s.Predicate.UID = 0
		s.Object.UID = 0
		dst.AddStatement(s)
	}
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
