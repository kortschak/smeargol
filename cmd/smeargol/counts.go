// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
)

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

	// geneIDs and geneIdx are mappings
	// between internal gene index and
	// external gene identifier.
	geneIDs []string
	geneIdx map[string]int
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
	geneIdx := make(map[string]int)
	var geneIDs []string

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
		geneIdx[geneid] = len(geneIDs)
		geneIDs = append(geneIDs, geneid)
		for i, f := range counts[1:] {
			v, err := strconv.ParseFloat(f, 64)
			if err != nil {
				return nil, fmt.Errorf("error parsing value for %q in sample %q: %v", geneid, samples[i], err)
			}
			data[geneid] = append(data[geneid], v)
		}
	}

	return &countData{
		names:   samples,
		counts:  data,
		geneIDs: geneIDs,
		geneIdx: geneIdx,
	}, nil
}
