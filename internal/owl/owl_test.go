// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package owl

import (
	"bytes"
	"compress/gzip"
	"encoding/xml"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/pkg/diff"
	"github.com/pkg/diff/write"

	"gonum.org/v1/gonum/graph/formats/rdf"
)

// Test data obtained from http://geneontology.org/docs/download-ontology/
// Gene Ontology Consortium data and data products are licensed under the
// Creative Commons Attribution 4.0 Unported License.
// https://creativecommons.org/licenses/by/4.0/legalcode
//
// N-Triple data sets were derived from the owl files with the following
// python using the rdflib library available from https://github.com/RDFLib/rdflib/.
//
//  #!/usr/bin/python3
//
//  import rdflib
//  import sys
//
//  g = rdflib.Graph()
//  g.load(sys.argv[1])
//
//  print(g.serialize(format='nt').decode("utf-8"))

func TestOwl(t *testing.T) {
	tests, err := filepath.Glob("testdata/goslim_*.owl.gz")
	if err != nil {
		t.Fatalf("failed to get test data paths: %v", err)
	}
	for _, path := range tests {
		name := strings.TrimSuffix(filepath.Base(path), ".gz")

		f, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		r, err := gzip.NewReader(f)
		if err != nil {
			t.Fatal(err)
		}

		var got []*rdf.Statement
		dec, err := NewDecoder(r)
		if err != nil {
			t.Fatal(err)
		}
		for {
			s, err := dec.Unmarshal()
			if err != nil {
				if err != io.EOF {
					t.Errorf("error during decoding: %v", err)
				}
				break
			}
			got = append(got, s)
		}
		f.Close()

		gotCan, err := rdf.URDNA2015(nil, got)
		if err != nil {
			t.Errorf("error during canonicalisation of %q: %v", name, err)
		}

		wantCan, err := canonicalFromNT(path)
		if err != nil {
			t.Errorf("error during golden data canonicalisation for %q: %v", name, err)
		}

		if !equalCanonicalGraphs(gotCan, wantCan) {
			var got, want strings.Builder
			for _, s := range gotCan {
				fmt.Fprintln(&got, s)
			}
			for _, s := range wantCan {
				fmt.Fprintln(&want, s)
			}
			var buf bytes.Buffer
			err := diff.Text("got", "want", got.String(), want.String(), &buf, write.TerminalColor())
			if err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			t.Errorf("unexpected canonical graph for %q:\n%s", name, &buf)
		}

		gotNamespaces := dec.Namespaces()
		wantNamespaces := namespacesFor(name)
		if !reflect.DeepEqual(gotNamespaces, wantNamespaces) {
			var got, want strings.Builder
			for _, n := range gotNamespaces {
				fmt.Fprintf(&got, "%+v\n", n)
			}
			for _, n := range wantNamespaces {
				fmt.Fprintf(&want, "%+v\n", n)
			}
			var buf bytes.Buffer
			err := diff.Text("got", "want", got.String(), want.String(), &buf, write.TerminalColor())
			if err != nil {
				t.Errorf("unexpected error: %v", err)
			}
			t.Errorf("unexpected namespaces returned for %q:\n%s", name, &buf)
		}
	}
}

func canonicalFromNT(path string) ([]*rdf.Statement, error) {
	path = strings.TrimSuffix(path, ".owl.gz") + ".nt.gz"
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}

	var statements []*rdf.Statement
	dec := rdf.NewDecoder(r)
	for {
		s, err := dec.Unmarshal()
		if err != nil {
			if err != io.EOF {
				return nil, err
			}
			break
		}
		statements = append(statements, s)
	}
	f.Close()

	return rdf.URDNA2015(nil, statements)
}

func equalCanonicalGraphs(a, b []*rdf.Statement) bool {
	if len(a) != len(b) {
		return false
	}
	for i, ai := range a {
		if ai.String() != b[i].String() {
			return false
		}
	}
	return true
}

func namespacesFor(file string) []xml.Attr {
	return []xml.Attr{
		{
			Name:  xml.Name{Space: "xmlns", Local: "terms"},
			Value: "http://www.geneontology.org/formats/oboInOwl#http://purl.org/dc/terms/",
		},
		{
			Name:  xml.Name{Space: "", Local: "xmlns"},
			Value: "http://purl.obolibrary.org/obo/go/subsets/" + file + "#",
		},
		{
			Name:  xml.Name{Space: "xml", Local: "base"},
			Value: "http://purl.obolibrary.org/obo/go/subsets/" + file,
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "oboInOwl"},
			Value: "http://www.geneontology.org/formats/oboInOwl#",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "rdf"},
			Value: "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "rdfs"},
			Value: "http://www.w3.org/2000/01/rdf-schema#",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "xml"},
			Value: "http://www.w3.org/XML/1998/namespace",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "go"},
			Value: "http://purl.obolibrary.org/obo/go#",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "xsd"},
			Value: "http://www.w3.org/2001/XMLSchema#",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "obo"},
			Value: "http://purl.obolibrary.org/obo/",
		},
		{
			Name:  xml.Name{Space: "xmlns", Local: "owl"},
			Value: "http://www.w3.org/2002/07/owl#",
		},
	}
}
