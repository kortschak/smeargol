// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math/big"
	"sync"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/traverse"

	"github.com/kortschak/gogo"
)

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
