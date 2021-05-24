// Copyright ©2021 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"math/big"
	"os"
	"strings"
	"sync"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/encoding"
	"gonum.org/v1/gonum/graph/encoding/dot"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/iterator"

	"github.com/kortschak/gogo"
)

type debugWriter struct {
	ontology *gogo.Graph
	ontoData []map[string]ontoCounts
	data     *countData

	mu     sync.Mutex
	depths map[string]int

	buffers []bytes.Buffer
}

func newDebugWriter(ontology *gogo.Graph, ontoData []map[string]ontoCounts, data *countData) *debugWriter {
	return &debugWriter{
		ontology: ontology,
		depths:   make(map[string]int),
		ontoData: ontoData,
		data:     data,
		buffers:  make([]bytes.Buffer, len(ontoData)),
	}
}

func (d *debugWriter) record(aspect, depth int, root, term rdf.Term) {
	if d == nil {
		return
	}
	d.mu.Lock()
	d.depths[term.Value] = depth
	d.mu.Unlock()
	counts, ok := d.ontoData[aspect][term.Value]
	if !ok {
		return
	}
	fmt.Fprintf(&d.buffers[aspect], "%s\t%s\t%s\t%d", strip(term.Value, "<obo:", ">"), strip(root.Value, "<obo:", ">"), nameSpaceOf(term, d.ontology), depth)
	for _, v := range counts.vector {
		fmt.Fprintf(&d.buffers[aspect], "\t%0*b", len(d.data.geneIDs), &v)
	}
	fmt.Fprintln(&d.buffers[aspect])
}

func (d *debugWriter) flush() {
	if d == nil {
		return
	}
	fmt.Printf("/*\ngo_term\tgo_root\tgo_aspect\tdepth\t%s\n", strings.Join(d.data.names, "\t"))
	for i := range d.buffers {
		_, err := io.Copy(os.Stdout, &d.buffers[i])
		if err != nil {
			log.Println(err)
		}
	}
	fmt.Println("*/")

	g := newDebugGraph(d.ontology, d.depths, d.data, d.ontoData)
	b, err := dot.MarshalMulti(g, "debug", "", "\t")
	if err != nil {
		log.Println(err)
		return
	}
	fmt.Printf("%s\n", b)
}

type debugGraph struct {
	*gogo.Graph

	depths   map[string]int
	data     *countData
	ontoData []map[string]ontoCounts
}

func newDebugGraph(g *gogo.Graph, depths map[string]int, data *countData, ontoData []map[string]ontoCounts) *debugGraph {
	c := gogo.NewGraph()
	it := g.AllStatements()
	for it.Next() {
		s := it.Statement()
		// Reverse annotation edges so that genes are placed under GO terms.
		if s.Predicate.Value == "<local:annotates>" {
			s = &rdf.Statement{Subject: s.Object, Predicate: s.Predicate, Object: s.Subject}
		}
		c.AddStatement(s)
	}
	return &debugGraph{
		Graph:    c,
		depths:   depths,
		data:     data,
		ontoData: ontoData,
	}
}

func (g *debugGraph) DOTAttributers() (graph, node, edge encoding.Attributer) {
	return attr{{Key: "rankdir", Value: "BT"}}, attr{}, attr{}
}

type attr []encoding.Attribute

func (a attr) Attributes() []encoding.Attribute {
	return a
}

func (g *debugGraph) Nodes() graph.Nodes {
	return g.filtered(g.Graph.Nodes())
}

func (g *debugGraph) From(uid int64) graph.Nodes {
	return g.filtered(g.Graph.From(uid))
}

func (g *debugGraph) filtered(it graph.Nodes) graph.Nodes {
	var dotNodes []graph.Node
	for it.Next() {
		term := it.Node().(rdf.Term)
		switch {
		case strings.HasPrefix(term.Value, "<ensembl:"):
			counts, ok := g.data.counts[strip(term.Value, "<ensembl:", ">")]
			if !ok {
				continue
			}
			for _, v := range counts {
				if v != 0 {
					dotNodes = append(dotNodes, geneNode{
						Term:   term,
						counts: counts,
					})
					break
				}
			}
		case strings.HasPrefix(term.Value, "<obo:GO_"):
			for _, aspect := range g.ontoData {
				counts, ok := aspect[term.Value]
				if !ok {
					continue
				}
				for i := range counts.vector {
					v := &counts.vector[i]
					if v.BitLen() != 0 {
						bits := make([]*big.Int, len(counts.vector))
						for j := range counts.vector {
							bits[j] = &counts.vector[j]
						}
						dotNodes = append(dotNodes, &goTermNode{
							Term:  term,
							depth: g.depths[term.Value],
							wid:   len(g.data.geneIDs),
							bits:  bits,
						})
						break
					}
				}
			}
		}
	}
	if len(dotNodes) == 0 {
		return graph.Empty
	}
	return iterator.NewOrderedNodes(dotNodes)
}

func (g *debugGraph) Lines(uid, vid int64) graph.Lines {
	it := g.Graph.Lines(uid, vid)
	lines := make([]graph.Line, 0, it.Len())
	for it.Next() {
		l := it.Line().(*rdf.Statement)
		switch l.Predicate.Value {
		case "<local:annotates>":
			lines = append(lines, dotLine{
				Statement: l,
				attrs: []encoding.Attribute{
					{Key: "label", Value: "annotates"},
					{Key: "dir", Value: "back"}, // Re-reverse the edge direction.
				},
			})
		default:
			lines = append(lines, dotLine{
				Statement: l,
				attrs: []encoding.Attribute{
					{Key: "label", Value: "subclass_of"},
				},
			})
		}
	}
	return iterator.NewOrderedLines(lines)
}

// geneNode implements graph.Node and dot.Node to allow the
// RDF term value to be given to the DOT encoder.
type geneNode struct {
	rdf.Term

	counts []float64
}

func (n geneNode) DOTID() string { return n.Term.Value }
func (n geneNode) Attributes() []encoding.Attribute {
	counts := make([]string, len(n.counts))
	for i, c := range n.counts {
		counts[i] = fmt.Sprintf("%d:%v", i, c)
	}
	return []encoding.Attribute{
		{Key: "label", Value: fmt.Sprintf("%s\n%s", strip(n.Value, "<ensembl:", ">"), strings.Join(counts, "\n"))},
	}
}

// goTermNode implements graph.Node and dot.Node to allow the
// RDF term value to be given to the DOT encoder.
type goTermNode struct {
	rdf.Term
	depth int
	wid   int
	bits  []*big.Int
}

func (n *goTermNode) DOTID() string { return n.Term.Value }
func (n *goTermNode) Attributes() []encoding.Attribute {
	bits := make([]string, len(n.bits))
	for i, b := range n.bits {
		bits[i] = fmt.Sprintf("%d:%0*b", i, n.wid, b)
	}
	return []encoding.Attribute{
		{Key: "label", Value: fmt.Sprintf("GO:%s [%d]\n%s", strip(n.Value, "<obo:GO_", ">"), n.depth, strings.Join(bits, "\n"))},
	}
}

// dotLine implements graph.Line and encoding.Attributer to
// allow the line's RDF term value to be given to the DOT
// encoder and for the nodes to be shimmed to the dotNode
// type.
//
// Because the graph here is directed and we are not performing
// any line reversals, it is safe not to implement the
// ReversedLine method on dotLine; it will never be called.
type dotLine struct {
	*rdf.Statement
	attrs []encoding.Attribute
}

func (l dotLine) From() graph.Node                 { return l.Subject }
func (l dotLine) To() graph.Node                   { return l.Object }
func (l dotLine) Attributes() []encoding.Attribute { return l.attrs }
