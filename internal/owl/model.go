// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package owl

import (
	"crypto/md5"
	"encoding/xml"
	"fmt"
	"hash"
	"reflect"
	"strings"

	"gonum.org/v1/gonum/graph/formats/rdf"
)

// This file contains the logic required to unmarshal the RDF/XML and map it
// to RDF N-Triples. If terms are missed in tests, this is the first place to
// look for fixing it.
//
// For a detailed description of the logic here, see:
// https://www.w3.org/TR/owl2-mapping-to-rdf/.

var (
	rdfType  = mustTerm(rdf.NewIRITerm("http://www.w3.org/1999/02/22-rdf-syntax-ns#type"))
	rdfFirst = mustTerm(rdf.NewIRITerm("http://www.w3.org/1999/02/22-rdf-syntax-ns#first"))
	rdfRest  = mustTerm(rdf.NewIRITerm("http://www.w3.org/1999/02/22-rdf-syntax-ns#rest"))
	rdfNil   = mustTerm(rdf.NewIRITerm("http://www.w3.org/1999/02/22-rdf-syntax-ns#nil"))

	firstName = xml.Name{Space: "http://www.w3.org/1999/02/22-rdf-syntax-ns#", Local: "first"}
)

type annotationProperty struct {
	XMLName xml.Name

	About string      `xml:"about,attr"`
	ID    rdfDataType `xml:"id"`

	Comment         []rdfDataType `xml:"comment"`
	HasDbXref       []rdfDataType `xml:"hasDbXref"`
	HasOBONamespace []rdfDataType `xml:"hasOBONamespace"`
	HasScope        []rdfDataType `xml:"hasScope"`
	IsClassLevel    []rdfDataType `xml:"is_class_level"`
	IsMetadataTag   []rdfDataType `xml:"is_metadata_tag"`
	Label           []rdfDataType `xml:"label"`
	Shorthand       []rdfDataType `xml:"shorthand"`
	SubPropertyOf   []rdfDataType `xml:"subPropertyOf"`
}

func (a annotationProperty) collect(dst []*rdf.Statement) []*rdf.Statement {
	claim, text, qual, kind := a.ID.claim()
	if kind != rdf.Invalid {
		subj := mustTerm(rdf.NewIRITerm(a.About))
		pred := mustTerm(rdf.NewIRITerm(claim))
		var obj rdf.Term
		switch kind {
		case rdf.IRI:
			obj = mustTerm(rdf.NewIRITerm(text))
		case rdf.Literal:
			obj = mustTerm(rdf.NewLiteralTerm(text, qual))
		}
		s := &rdf.Statement{Subject: subj, Predicate: pred, Object: obj}
		dst = append(dst, s)
	}
	dst = labelType(dst, a.About, a.XMLName)
	dst = collect(dst, a.About, a)
	return dst
}

type axiom struct {
	XMLName xml.Name

	Source   rdfDataType `xml:"annotatedSource"`
	Property rdfDataType `xml:"annotatedProperty"`
	Target   rdfDataType `xml:"annotatedTarget"`

	Comment   []rdfDataType `xml:"comment"`
	HasDbXref []rdfDataType `xml:"hasDbXref"`
	Label     []rdfDataType `xml:"label"`
}

func (a axiom) collect(dst []*rdf.Statement) []*rdf.Statement {
	label := blankLabel(md5.New(),
		a.Source.Resource,
		a.Property.Resource,
		a.Target.Resource, a.Target.Text, a.Target.Datatype)
	blank := mustTerm(rdf.NewBlankTerm(label))

	typ := mustTerm(rdf.NewIRITerm(a.XMLName.Space + a.XMLName.Local))

	source := mustTerm(rdf.NewIRITerm(a.Source.Resource))
	sourcePred := mustTerm(rdf.NewIRITerm(a.Source.XMLName.Space + a.Source.XMLName.Local))

	property := mustTerm(rdf.NewIRITerm(a.Property.Resource))
	propertyPred := mustTerm(rdf.NewIRITerm(a.Property.XMLName.Space + a.Property.XMLName.Local))

	_, text, qual, kind := a.Target.claim()
	var target rdf.Term
	switch kind {
	case rdf.Literal:
		target = mustTerm(rdf.NewLiteralTerm(text, qual))
	case rdf.IRI:
		target = mustTerm(rdf.NewIRITerm(text))
	default:
		panic(fmt.Sprintf("unexpected term kind: %s", kind))
	}
	targetPred := mustTerm(rdf.NewIRITerm(a.Target.XMLName.Space + a.Target.XMLName.Local))

	dst = append(dst,
		&rdf.Statement{Subject: blank, Predicate: rdfType, Object: typ},
		&rdf.Statement{Subject: blank, Predicate: sourcePred, Object: source},
		&rdf.Statement{Subject: blank, Predicate: propertyPred, Object: property},
		&rdf.Statement{Subject: blank, Predicate: targetPred, Object: target},
	)

	for _, links := range [][]rdfDataType{
		a.Label,
		a.HasDbXref,
		a.Comment,
	} {
		for _, p := range links {
			pred := mustTerm(rdf.NewIRITerm(p.XMLName.Space + p.XMLName.Local))
			_, text, qual, kind := p.claim()
			var obj rdf.Term
			switch kind {
			case rdf.Literal:
				obj = mustTerm(rdf.NewLiteralTerm(text, qual))
			case rdf.IRI:
				obj = mustTerm(rdf.NewIRITerm(text))
			default:
				panic(fmt.Sprintf("unexpected term kind: %s", kind))
			}
			dst = append(dst, &rdf.Statement{Subject: blank, Predicate: pred, Object: obj})
		}
	}

	return dst
}

type class struct {
	XMLName xml.Name

	About string      `xml:"about,attr"`
	ID    rdfDataType `xml:"id"`

	HasBroadSynonym   []rdfDataType `xml:"hasBroadSynonym"`
	Comment           []rdfDataType `xml:"comment"`
	Consider          []rdfDataType `xml:"consider"`
	CreatedBy         []rdfDataType `xml:"created_by"`
	CreationDate      []rdfDataType `xml:"creation_date"`
	HasDbXref         []rdfDataType `xml:"hasDbXref"`
	Deprecated        []rdfDataType `xml:"deprecated"`
	DisjointWith      []rdfDataType `xml:"disjointWith"`
	HasExactSynonym   []rdfDataType `xml:"hasExactSynonym"`
	HasAlternativeID  []rdfDataType `xml:"hasAlternativeId"`
	HasOBONamespace   []rdfDataType `xml:"hasOBONamespace"`
	IAO0000115        []rdfDataType `xml:"IAO_0000115"`
	IAO0000233        []rdfDataType `xml:"IAO_0000233"`
	IAO0000589        []rdfDataType `xml:"IAO_0000589"`
	IAO0100001        []rdfDataType `xml:"IAO_0100001"`
	InSubSet          []rdfDataType `xml:"inSubset"`
	Label             []rdfDataType `xml:"label"`
	HasNarrowSynonym  []rdfDataType `xml:"hasNarrowSynonym"`
	HasRelatedSynonym []rdfDataType `xml:"hasRelatedSynonym"`
	RO0002161         []rdfDataType `xml:"RO_0002161"`

	EquivalentClass equivalentClasses `xml:"equivalentClass"`
	SubClassOf      subClassOfs       `xml:"subClassOf"`
}

func (c class) collect(dst []*rdf.Statement) []*rdf.Statement {
	claim, text, qual, kind := c.ID.claim()
	if kind != rdf.Invalid {
		subj := mustTerm(rdf.NewIRITerm(c.About))
		pred := mustTerm(rdf.NewIRITerm(claim))
		var obj rdf.Term
		switch kind {
		case rdf.IRI:
			obj = mustTerm(rdf.NewIRITerm(text))
		case rdf.Literal:
			obj = mustTerm(rdf.NewLiteralTerm(text, qual))
		}
		s := &rdf.Statement{Subject: subj, Predicate: pred, Object: obj}
		dst = append(dst, s)
	}

	dst = labelType(dst, c.About, c.XMLName)
	dst = collect(dst, c.About, c)
	dst = c.SubClassOf.collect(dst, c.About)
	dst = c.EquivalentClass.collect(dst, c.About, c.XMLName)
	return dst
}

type subClassOf struct {
	XMLName xml.Name

	Resource    string       `xml:"resource,attr"`
	Restriction restrictions `xml:"Restriction"`
}

func (s subClassOf) claim() (pred, obj, qual string, kind rdf.Kind) {
	if s.Resource == "" {
		return
	}
	return s.XMLName.Space + s.XMLName.Local, s.Resource, "", rdf.IRI
}

type subClassOfs []subClassOf

func (l subClassOfs) collect(dst []*rdf.Statement, about string) []*rdf.Statement {
	if l == nil {
		return dst
	}
	for _, c := range l {
		claim, text, qual, kind := c.claim()
		if kind != rdf.Invalid {
			subj := mustTerm(rdf.NewIRITerm(about))
			pred := mustTerm(rdf.NewIRITerm(claim))
			var obj rdf.Term
			switch kind {
			case rdf.IRI:
				obj = mustTerm(rdf.NewIRITerm(text))
			case rdf.Literal:
				obj = mustTerm(rdf.NewLiteralTerm(text, qual))
			}
			s := &rdf.Statement{Subject: subj, Predicate: pred, Object: obj}
			dst = append(dst, s)
		}
		dst = c.Restriction.collect(dst, about, c.XMLName)
	}
	return dst
}

type equivalentClass struct {
	XMLName xml.Name

	Class []intersectClass `xml:"Class"`
}

type intersectClass struct {
	XMLName xml.Name

	IntersectionOf []intersectionOf `xml:"intersectionOf"`
}

type intersectionOf struct {
	XMLName xml.Name

	ParseType   string        `xml:"parseType,attr"`
	Description []description `xml:"Description"`
	Restriction restrictions  `xml:"Restriction"`
}

type description struct {
	XMLName xml.Name

	About string `xml:"about,attr"`
}

func (e description) claim() (pred, obj, qual string, kind rdf.Kind) {
	if e.About == "" {
		return
	}
	return e.XMLName.Space + e.XMLName.Local, e.About, "", rdf.IRI
}

type equivalentClasses []equivalentClass

func (l equivalentClasses) collect(dst []*rdf.Statement, about string, in xml.Name) []*rdf.Statement {
	if l == nil {
		return dst
	}

	h := md5.New()

	subj := mustTerm(rdf.NewIRITerm(about))
	typ := mustTerm(rdf.NewIRITerm(in.Space + in.Local))
	s := &rdf.Statement{Subject: subj, Predicate: rdfType, Object: typ}
	dst = append(dst, s)

	var ecLabel string
	for _, equivalentClass := range l {
		ecLabel = blankLabel(h,
			equivalentClass.XMLName.Space, equivalentClass.XMLName.Local,
			about, in.Space, in.Local, ecLabel,
		)

		pred := mustTerm(rdf.NewIRITerm(equivalentClass.XMLName.Space + equivalentClass.XMLName.Local))
		ecBlank := mustTerm(rdf.NewBlankTerm(ecLabel))
		s := &rdf.Statement{Subject: subj, Predicate: pred, Object: ecBlank}
		dst = append(dst, s)

		for _, class := range equivalentClass.Class {
			typ := mustTerm(rdf.NewIRITerm(in.Space + in.Local))
			s := &rdf.Statement{Subject: ecBlank, Predicate: rdfType, Object: typ}
			dst = append(dst, s)

			cLabel := blankLabel(h,
				class.XMLName.Space, class.XMLName.Local,
				ecLabel,
			)

			for i, intersectionOf := range class.IntersectionOf {
				iLabel := blankLabel(h,
					intersectionOf.XMLName.Space, intersectionOf.XMLName.Local,
					cLabel, intersectionOf.ParseType,
				)

				iBlank := mustTerm(rdf.NewBlankTerm(iLabel))

				if i == 0 {
					// _:ecBlank <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#Class> .
					// _:ecBlank <http://www.w3.org/2002/07/owl#intersectionOf> _:iBlank .
					pred := mustTerm(rdf.NewIRITerm(intersectionOf.XMLName.Space + intersectionOf.XMLName.Local))
					dst = append(dst,
						&rdf.Statement{Subject: ecBlank, Predicate: rdfType, Object: typ},
						&rdf.Statement{Subject: ecBlank, Predicate: pred, Object: iBlank},
					)
				}

				var dLabel string
				for i, description := range intersectionOf.Description {

					dLabel = blankLabel(h,
						description.XMLName.Space, description.XMLName.Local,
						iLabel, description.About, dLabel,
					)

					_, text, qual, kind := description.claim()
					if kind != rdf.Invalid {
						var obj rdf.Term
						switch kind {
						case rdf.IRI:
							obj = mustTerm(rdf.NewIRITerm(text))
						case rdf.Literal:
							obj = mustTerm(rdf.NewLiteralTerm(text, qual))
						}

						s = &rdf.Statement{Subject: iBlank, Predicate: rdfFirst, Object: obj}
						dst = append(dst, s)
					}

					// _:iBlank <http://www.w3.org/1999/02/22-rdf-syntax-ns#rest> _:nextIBlank .
					iLabel = blankLabel(h, iLabel)
					nextIBlank := mustTerm(rdf.NewBlankTerm(iLabel))
					s = &rdf.Statement{Subject: iBlank, Predicate: rdfRest, Object: nextIBlank}
					dst = append(dst, s)
					iBlank = nextIBlank

					if i < len(intersectionOf.Restriction) {
						dst = intersectionOf.Restriction[i].collect(dst, iLabel, firstName)
					}

					if i < len(intersectionOf.Description)-1 {
						// _:iBlank <http://www.w3.org/1999/02/22-rdf-syntax-ns#rest> _:nextIBlank .
						iLabel = blankLabel(h, iLabel)
						nextIBlank := mustTerm(rdf.NewBlankTerm(iLabel))
						s = &rdf.Statement{Subject: iBlank, Predicate: rdfRest, Object: nextIBlank}
						iBlank = nextIBlank
					} else {
						// _:iBlank <http://www.w3.org/1999/02/22-rdf-syntax-ns#rest> <http://www.w3.org/1999/02/22-rdf-syntax-ns#nil> .
						s = &rdf.Statement{Subject: iBlank, Predicate: rdfRest, Object: rdfNil}
					}
					dst = append(dst, s)
				}
			}
		}
	}
	return dst
}

type restriction struct {
	XMLName xml.Name

	OnProperty     rdfDataType `xml:"onProperty"`
	SomeValuesFrom rdfDataType `xml:"someValuesFrom"`
}

func (r restriction) claim() (pred, obj, qual string, kind rdf.Kind) {
	return "", r.XMLName.Space + r.XMLName.Local, "", rdf.IRI
}

func (r restriction) collect(dst []*rdf.Statement, about string, in xml.Name) []*rdf.Statement {
	h := md5.New()

	// Rule swarmlet: try IRI and then blank node.
	subj, err := rdf.NewIRITerm(about)
	if err != nil {
		subj = mustTerm(rdf.NewBlankTerm(about))
	}

	parent := mustTerm(rdf.NewIRITerm(in.Space + in.Local))
	typ := mustTerm(rdf.NewIRITerm(r.XMLName.Space + r.XMLName.Local))

	label := blankLabel(h,
		r.XMLName.Space, r.XMLName.Local,
		about, in.Space, in.Local,
		r.OnProperty.Resource, r.OnProperty.Text,
		r.SomeValuesFrom.Resource, r.SomeValuesFrom.Text)

	blank := mustTerm(rdf.NewBlankTerm(label))

	dst = append(dst,
		&rdf.Statement{Subject: subj, Predicate: parent, Object: blank},
		&rdf.Statement{Subject: blank, Predicate: rdfType, Object: typ})

	for _, d := range []rdfDataType{
		r.OnProperty,
		r.SomeValuesFrom,
	} {
		if claim, text, qual, kind := d.claim(); kind != rdf.Invalid {
			pred := mustTerm(rdf.NewIRITerm(claim))
			var obj rdf.Term
			switch kind {
			case rdf.IRI:
				obj = mustTerm(rdf.NewIRITerm(text))
			case rdf.Literal:
				obj = mustTerm(rdf.NewLiteralTerm(text, qual))
			}
			s := &rdf.Statement{Subject: blank, Predicate: pred, Object: obj}
			dst = append(dst, s)
		}
	}

	return dst
}

// https://www.w3.org/TR/owl-ref/#Restriction
type restrictions []restriction

func (l restrictions) collect(dst []*rdf.Statement, about string, in xml.Name) []*rdf.Statement {
	if len(l) == 0 {
		return dst
	}

	for _, r := range l {
		dst = r.collect(dst, about, in)
	}

	return dst
}

type objectProperty struct {
	XMLName xml.Name

	About string      `xml:"about,attr"`
	ID    rdfDataType `xml:"id"`

	HasDbXref       []rdfDataType `xml:"hasDbXref"`
	HasOBONamespace []rdfDataType `xml:"hasOBONamespace"`
	InverseOf       []rdfDataType `xml:"inverseOf"`
	Label           []rdfDataType `xml:"label"`
	Shorthand       []rdfDataType `xml:"shorthand"`
	SubPropertyOf   []rdfDataType `xml:"subPropertyOf"`
	Type            []rdfDataType `xml:"type"`

	PropertyChainAxiom *propertyChainAxiom `xml:"propertyChainAxiom"`
}

type propertyChainAxiom struct {
	XMLName xml.Name

	ParseType   string        `xml:"parseType,attr"`
	Description []description `xml:"Description"`
}

func (o objectProperty) collect(dst []*rdf.Statement) []*rdf.Statement {
	claim, text, qual, kind := o.ID.claim()
	if kind != rdf.Invalid {
		subj := mustTerm(rdf.NewIRITerm(o.About))
		pred := mustTerm(rdf.NewIRITerm(claim))
		var obj rdf.Term
		switch kind {
		case rdf.IRI:
			obj = mustTerm(rdf.NewIRITerm(text))
		case rdf.Literal:
			obj = mustTerm(rdf.NewLiteralTerm(text, qual))
		}
		s := &rdf.Statement{Subject: subj, Predicate: pred, Object: obj}
		dst = append(dst, s)
	}
	dst = labelType(dst, o.About, o.XMLName)
	dst = collect(dst, o.About, o)
	dst = o.PropertyChainAxiom.collect(dst, o.About)
	return dst
}

func (c *propertyChainAxiom) collect(dst []*rdf.Statement, about string) []*rdf.Statement {
	if c == nil {
		return dst
	}

	h := md5.New()
	label := blankLabel(h, about, c.XMLName.Space+c.XMLName.Local)

	subj := mustTerm(rdf.NewIRITerm(about))
	pred := mustTerm(rdf.NewIRITerm(c.XMLName.Space + c.XMLName.Local))
	blank := mustTerm(rdf.NewBlankTerm(label))
	s := &rdf.Statement{Subject: subj, Predicate: pred, Object: blank}
	dst = append(dst, s)

	for i, d := range c.Description {
		obj := mustTerm(rdf.NewIRITerm(d.About))
		s := &rdf.Statement{Subject: blank, Predicate: rdfFirst, Object: obj}
		dst = append(dst, s)
		if i < len(c.Description)-1 {
			// _:blank <http://www.w3.org/1999/02/22-rdf-syntax-ns#rest> _:nextBlank .
			label = blankLabel(h, label)
			nextBlank := mustTerm(rdf.NewBlankTerm(label))
			s = &rdf.Statement{Subject: blank, Predicate: rdfRest, Object: nextBlank}
			blank = nextBlank
		} else {
			// _:blank <http://www.w3.org/1999/02/22-rdf-syntax-ns#rest> <http://www.w3.org/1999/02/22-rdf-syntax-ns#nil> .
			s = &rdf.Statement{Subject: blank, Predicate: rdfRest, Object: rdfNil}
		}
		dst = append(dst, s)
	}
	return dst
}

type ontology struct {
	XMLName xml.Name

	About string `xml:"about,attr"`

	DefaultNamespace    []rdfDataType `xml:"default-namespace"`
	Description         []rdfDataType `xml:"description"`
	HasOBOFormatVersion []rdfDataType `xml:"hasOBOFormatVersion"`
	License             []rdfDataType `xml:"license"`
	Title               []rdfDataType `xml:"title"`
	VersionIRI          []rdfDataType `xml:"versionIRI"`
}

func (o ontology) collect(dst []*rdf.Statement) []*rdf.Statement {
	dst = labelType(dst, o.About, o.XMLName)
	dst = collect(dst, o.About, o)
	return dst
}

func labelType(dst []*rdf.Statement, about string, name xml.Name) []*rdf.Statement {
	subj := mustTerm(rdf.NewIRITerm(about))
	typ := mustTerm(rdf.NewIRITerm(name.Space + name.Local))
	s := &rdf.Statement{Subject: subj, Predicate: rdfType, Object: typ}
	dst = append(dst, s)
	return dst
}

type rdfDataType struct {
	XMLName xml.Name

	Resource    string        `xml:"resource,attr"`
	Text        string        `xml:",chardata"`
	Datatype    string        `xml:"datatype,attr"`
	Restriction []restriction `xml:"Restriction"`
}

func (r rdfDataType) claim() (pred, obj, qual string, kind rdf.Kind) {
	switch {
	case strings.TrimSpace(r.Resource) != "":
		return r.XMLName.Space + r.XMLName.Local, r.Resource, "", rdf.IRI
	case strings.TrimSpace(r.Datatype) != "":
		return r.XMLName.Space + r.XMLName.Local, r.Text, r.Datatype, rdf.Literal
	default:
		return
	}
}

func collect(dst []*rdf.Statement, about string, v interface{}) []*rdf.Statement {
	rv := reflect.ValueOf(v)
	for i := 0; i < rv.NumField(); i++ {
		f := rv.Field(i)
		if f.Kind() != reflect.Slice {
			continue
		}
		ft := f.Type()
		if ft.Elem() != reflect.TypeOf(rdfDataType{}) {
			continue
		}
		for _, e := range f.Interface().([]rdfDataType) {
			claim, text, qual, kind := e.claim()
			if kind != rdf.Invalid {
				subj := mustTerm(rdf.NewIRITerm(about))
				pred := mustTerm(rdf.NewIRITerm(claim))
				var obj rdf.Term
				switch kind {
				case rdf.IRI:
					obj = mustTerm(rdf.NewIRITerm(text))
				case rdf.Literal:
					obj = mustTerm(rdf.NewLiteralTerm(text, qual))
				}
				s := &rdf.Statement{Subject: subj, Predicate: pred, Object: obj}
				dst = append(dst, s)
			}
		}
	}
	return dst
}

func blankLabel(h hash.Hash, parts ...string) string {
	h.Reset()
	for _, p := range parts {
		h.Write([]byte(p)) //nolint:errcheck
	}
	return hex(h.Sum(nil))
}

func hex(data []byte) string {
	const digit = "0123456789abcdef"
	buf := make([]byte, 0, len(data)*2)
	for _, b := range data {
		buf = append(buf, digit[b>>4], digit[b&0xf])
	}
	return string(buf)
}

func mustTerm(t rdf.Term, err error) rdf.Term {
	if err != nil {
		panic(err)
	}
	return t
}
