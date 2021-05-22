// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// goglink maps Ensembl ENSG gene identifiers to GO terms based on
// Ensembl database cross-reference data.
package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"

	"gonum.org/v1/gonum/graph/formats/rdf"

	"github.com/kortschak/gogo"
)

func main() {
	var (
		orgPath  = flag.String("org", "", "specify the Ensembl organism data (.nt.gz/.nq.gz - required)")
		xrefPath = flag.String("xref", "", "specify the Ensembl xref data (.nt.gz/.nq.gz - required)")
		help     = flag.Bool("help", false, "print help text")
	)

	flag.Parse()

	if *help {
		flag.Usage()
		fmt.Fprintf(os.Stderr, `
%s maps ENSG identifiers to GO terms based on Ensembl cross-reference
data. It outputs the mapping as RDF triples in the form:

 <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .

for each GO term to Ensembl gene annotation.

Input data can be obtained from ftp://ftp.ensembl.org/pub/current_rdf
in Turtle format. These files must first be converted to N-Triples.

All input files are expected to be gzip compressed and the output is
written uncompressed to standard output and standard error.

Copyright ©2020 Dan Kortschak. All rights reserved.

`, filepath.Base(os.Args[0]))
		os.Exit(0)
	}

	if *orgPath == "" || *xrefPath == "" {
		flag.Usage()
		os.Exit(2)
	}

	g := gogo.NewGraph()
	var dec rdf.Decoder
	for _, path := range []string{*orgPath, *xrefPath} {
		f, err := os.Open(path)
		if err != nil {
			log.Fatal(err)
		}
		r, err := gzip.NewReader(f)
		if err != nil {
			log.Fatal(err)
		}

		dec.Reset(r)
		for {
			s, err := dec.Unmarshal()
			if err != nil {
				if err != io.EOF {
					log.Fatalf("error during decoding: %v", err)
				}
				break
			}

			switch s.Predicate.Value {
			case "<obo:SO_transcribed_from>":
			case "<http://purl.obolibrary.org/obo/SO_transcribed_from>":
				s.Subject.Value = "<transcript:" + strings.TrimPrefix(s.Subject.Value, "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/")
				s.Predicate.Value = "<obo:SO_transcribed_from>"
				s.Object.Value = "<ensembl:" + strings.TrimPrefix(s.Object.Value, "<http://rdf.ebi.ac.uk/resource/ensembl/")
			case "<rdfs:seeAlso>":
			case "<http://www.w3.org/2000/01/rdf-schema#seeAlso>":
				if !strings.HasPrefix(s.Object.Value, "<http://identifiers.org/go/GO:") {
					continue
				}
				s.Subject.Value = "<transcript:" + strings.TrimPrefix(s.Subject.Value, "<http://rdf.ebi.ac.uk/resource/ensembl.transcript/")
				s.Predicate.Value = "<rdfs:seeAlso>"
				s.Object.Value = "<obo:GO_" + strings.TrimPrefix(s.Object.Value, "<http://identifiers.org/go/GO:")
			default:
				continue
			}

			s.Subject.UID = 0
			s.Predicate.UID = 0
			s.Object.UID = 0

			g.AddStatement(s)
		}
		f.Close()
	}

	nodes := g.Nodes()
	for nodes.Next() {
		gene := nodes.Node().(rdf.Term)
		if !strings.HasPrefix(gene.Value, "<ensembl:") {
			continue
		}

		// We are emitting directly, so we need to ensure statement
		// uniqueness. A seen per start node is enough for this. If
		// we were adding to another graph, the deduplication could
		// be handled by the destination graph.
		seen := make(map[int64]bool)

		// Get all GO terms reachable from the ENSG via an ENST
		// since that is how the Ensembl GO annotation work.
		terms := g.Query(gene).In(func(s *rdf.Statement) bool {
			// <transcript:Y> <obo:SO_transcribed_from> <ensembl:X> .
			return s.Predicate.Value == "<obo:SO_transcribed_from>"

		}).Out(func(s *rdf.Statement) bool {
			if seen[s.Object.UID] {
				return false
			}

			// <transcript:Y> <rdfs:seeAlso> <obo:GO_Z> .
			ok := s.Predicate.Value == "<rdfs:seeAlso>" &&
				strings.HasPrefix(s.Object.Value, "<obo:GO_")
			if ok {
				seen[s.Object.UID] = true
			}
			return ok

		}).Result()

		for _, t := range terms {
			fmt.Println(&rdf.Statement{
				Subject:   rdf.Term{Value: t.Value},
				Predicate: rdf.Term{Value: "<local:annotates>"},
				Object:    rdf.Term{Value: gene.Value},
			})
		}
	}
}
