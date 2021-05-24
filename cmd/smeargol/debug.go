// Copyright ©2021 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"strings"

	"github.com/kortschak/gogo"
	"gonum.org/v1/gonum/graph/formats/rdf"
)

type debugWriter struct {
	ontology *gogo.Graph
	ontoData []map[string]ontoCounts
	data     *countData

	buffers []bytes.Buffer
}

func newDebugWriter(ontology *gogo.Graph, ontoData []map[string]ontoCounts, data *countData) *debugWriter {
	return &debugWriter{
		ontology: ontology,
		ontoData: ontoData,
		data:     data,
		buffers:  make([]bytes.Buffer, len(ontoData)),
	}
}

func (d *debugWriter) record(aspect, depth int, root, term rdf.Term) {
	if d == nil {
		return
	}
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
	fmt.Printf("go_term\tgo_root\tgo_aspect\tdepth\t%s\n", strings.Join(d.data.names, "\t"))
	for i := range d.buffers {
		_, err := io.Copy(os.Stdout, &d.buffers[i])
		if err != nil {
			log.Println(err)
		}
	}
}
