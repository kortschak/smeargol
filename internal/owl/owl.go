// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package owl

import (
	"encoding/xml"
	"fmt"
	"io"
	"sort"
	"strings"

	"gonum.org/v1/gonum/graph/formats/rdf"
)

// Decoder is a Gene Ontology OBO in OWL decoder. rdf.Statements returned
// by calls to the Unmarshal and UnmarshalLocal methods have their Terms'
// UID fields set so that unique terms will have unique IDs and so can be
// used directly in a graph.Multi, or in a graph.Graph if all predicate
// terms are identical. IDs created by the decoder all exist within a single
// namespace and so Terms can be uniquely identified by their UID. Term
// UIDs are based from 1 to allow RDF-aware client graphs to assign ID if
// no ID has been assigned.
type Decoder struct {
	xml        *xml.Decoder
	namespaces []xml.Attr

	strings store
	ids     map[string]int64

	curr int
	buf  []*rdf.Statement
	seen map[[3]int64]bool
}

// NewDecoder returns a new Decoder that takes input from r.
func NewDecoder(r io.Reader) (*Decoder, error) {
	dec := &Decoder{
		xml:     xml.NewDecoder(r),
		strings: make(store),
		ids:     make(map[string]int64),
		seen:    make(map[[3]int64]bool),
	}
	for dec.namespaces == nil {
		err := dec.fillBuffer()
		if err != nil {
			return nil, err
		}
	}
	return dec, nil
}

// Reset resets the decoder to use the provided io.Reader, retaining
// the existing Term ID mapping. A new XML namespace is obtained from the
// XML stream in r.
func (dec *Decoder) Reset(r io.Reader) error {
	dec.namespaces = nil
	dec.xml = xml.NewDecoder(r)
	for dec.namespaces == nil {
		err := dec.fillBuffer()
		if err != nil {
			return err
		}
	}
	for i := range dec.buf[dec.curr:] {
		dec.buf[dec.curr+i] = nil
	}
	dec.curr = 0
	dec.buf = dec.buf[:0]
	dec.strings = make(store)
	return nil
}

// Namespace returns the namespace collected from the XML stream. The
// value returned by Namespaces in valid after the Decoder is returned
// by NewDecoder or a successful Reset.
func (dec *Decoder) Namespaces() []xml.Attr {
	return dec.namespaces
}

// Unmarshal returns the next unique statement from the input stream.
func (dec *Decoder) Unmarshal() (*rdf.Statement, error) {
	for {
		for len(dec.buf[dec.curr:]) == 0 {
			err := dec.fillBuffer()
			if err != nil {
				return nil, err
			}
		}
		s := dec.buf[dec.curr]
		dec.buf[dec.curr] = nil
		dec.curr++
		if len(dec.buf[dec.curr:]) == 0 {
			dec.curr = 0
			dec.buf = dec.buf[:0]
		}
		s.Subject.Value = dec.strings.intern(s.Subject.Value)
		s.Predicate.Value = dec.strings.intern(s.Predicate.Value)
		s.Object.Value = dec.strings.intern(s.Object.Value)
		s.Subject.UID = dec.idFor(s.Subject.Value)
		s.Object.UID = dec.idFor(s.Object.Value)
		s.Predicate.UID = dec.idFor(s.Predicate.Value)
		triple := [3]int64{s.Subject.UID, s.Predicate.UID, s.Object.UID}
		if !dec.seen[triple] {
			dec.seen[triple] = true
			return s, nil
		}
	}
}

// UnmarshalLocal returns the next unique statement from the input stream, but
// replaces full IRI namespace text with the qualified name prefix obtained
// from the decoders internal namespaces. The namespaces can be obtained by
// using the Namespaces method.
func (dec *Decoder) UnmarshalLocal() (*rdf.Statement, error) {
	s, err := dec.Unmarshal()
	if err != nil {
		return nil, err
	}
	subj, err := dec.compactTerm(s.Subject)
	if err != nil {
		return s, err
	}
	s.Subject = subj
	pred, err := dec.compactTerm(s.Predicate)
	if err != nil {
		return s, err
	}
	s.Predicate = pred
	obj, err := dec.compactTerm(s.Object)
	if err != nil {
		return s, err
	}
	s.Object = obj
	return s, nil
}

func (dec *Decoder) compactTerm(term rdf.Term) (rdf.Term, error) {
	text, qual, kind, err := term.Parts()
	if err != nil {
		return term, err
	}
	uid := term.UID
	switch kind {
	case rdf.IRI:
		new, changed := dec.compactIRI(text)
		if changed {
			term, err := rdf.NewIRITerm(new)
			if err != nil {
				return term, err
			}
			term.UID = uid
			return term, nil
		}
	case rdf.Literal:
		if qual == "" {
			return term, nil
		}
		new, changed := dec.compactIRI(qual)
		if changed {
			term, err := rdf.NewLiteralTerm(text, new)
			if err != nil {
				return term, err
			}
			term.UID = uid
			return term, nil
		}
	}
	return term, nil
}

func (dec *Decoder) compactIRI(iri string) (new string, changed bool) {
	// dec.namespaces is ordered longest to shortest
	// to ensure prefixes are not eagerly chosen.
	for _, ns := range dec.namespaces {
		if strings.HasPrefix(iri, ns.Value) {
			suffix := strings.TrimPrefix(iri, ns.Value)
			if len(suffix) == 0 {
				return iri, false
			}
			return ns.Name.Local + ":" + strings.TrimPrefix(iri, ns.Value), true
		}
	}
	return iri, false
}

func (dec *Decoder) idFor(s string) int64 {
	id, ok := dec.ids[s]
	if ok {
		return id
	}
	id = int64(len(dec.ids)) + 1
	dec.ids[s] = id
	return id
}

func (dec *Decoder) fillBuffer() (err error) {
	defer func() {
		r := recover()
		switch r := r.(type) {
		case nil:
			return
		case error:
			err = r
		default:
			panic(r)
		}
	}()
	tok, err := dec.xml.Token()
	if err != nil {
		if err == io.EOF {
			dec.strings = nil
		}
		return err
	}
	switch tok := tok.(type) {
	case xml.StartElement:
		switch tok.Name.Local {
		case "AnnotationProperty":
			var a annotationProperty
			err = dec.xml.DecodeElement(&a, &tok)
			if err != nil {
				return err
			}
			dec.buf = a.collect(dec.buf)

		case "Axiom":
			var a axiom
			err = dec.xml.DecodeElement(&a, &tok)
			if err != nil {
				return err
			}
			dec.buf = a.collect(dec.buf)

		case "Class":
			var c class
			err = dec.xml.DecodeElement(&c, &tok)
			if err != nil {
				return err
			}
			dec.buf = c.collect(dec.buf)

		case "ObjectProperty":
			var o objectProperty
			err = dec.xml.DecodeElement(&o, &tok)
			if err != nil {
				return err
			}
			dec.buf = o.collect(dec.buf)

		case "Ontology":
			var o ontology
			err = dec.xml.DecodeElement(&o, &tok)
			if err != nil {
				return err
			}
			dec.buf = o.collect(dec.buf)

		case "RDF":
			for _, attr := range tok.Attr {
				if attr.Name.Space == "http://www.w3.org/XML/1998/namespace" {
					attr.Name.Space = "xml"
				}
				dec.namespaces = append(dec.namespaces, attr)
			}
			sort.Sort(byLength(dec.namespaces))

		default:
			panic(fmt.Sprintf("%+v", tok.Name))
		}

	case xml.EndElement:
	case xml.CharData:
	case xml.Comment:
	case xml.Directive:
	case xml.ProcInst:
	}
	return nil
}

// store is a string internment implementation.
type store map[string]string

// intern returns an interned version of the parameter.
func (is store) intern(s string) string {
	if s == "" {
		return ""
	}
	t, ok := is[s]
	if ok {
		return t
	}
	is[s] = s
	return s
}

type byLength []xml.Attr

func (a byLength) Len() int           { return len(a) }
func (a byLength) Less(i, j int) bool { return len(a[i].Value) > len(a[j].Value) }
func (a byLength) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
