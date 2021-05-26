// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

type SummaryDoc struct {
	// Roots is the set of roots in the Gene Ontology.
	Roots []string

	// Summaries contains the summaries of a smeargol
	// analysis.
	Summaries [][]*Summary
}

type Summary struct {
	// Name is the name of the sample.
	Name string

	// Root is the root GO term for the summary.
	Root string

	// Depth is the distance from the root.
	Depth int

	// Rows and Cols are the dimensions of the matrix
	// describing the GO level. Rows corresponds to the
	// number of genes and Cols corresponds to the number
	// of GO terms in the level.
	Rows, Cols int

	// OptimalRank and FractionalRank are the calculated
	// ranks of the summary matrix. OptimalRank is
	// calculated according to the method of Matan Gavish
	// and David L. Donoho https://arxiv.org/abs/1305.5870.
	// FractionalRank is the rank calculated using the
	// user-provided fraction parameters.
	OptimalRank, FractionalRank int

	// Sigma is the complete set of singular values.
	Sigma []float64
}

// https://arxiv.org/abs/1305.5870
func optimalTruncation(path string, m *mat.Dense, cut, frac float64) (*Summary, error) {
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDThin)
	if !ok {
		return nil, fmt.Errorf("could not factorise %q", path)
	}
	sigma := svd.Values(nil)

	sum := make([]float64, len(sigma))
	floats.CumSum(sum, sigma)
	var rFrac int
	var f float64
	max := sum[len(sum)-1]
	if max != 0 {
		floats.Scale(1/max, sum)
		rFrac = idxAbove(frac, sum)
		switch {
		case rFrac < len(sigma):
			f = sigma[rFrac]
		case len(sigma) != 0:
			f = sigma[0]
		}
	}

	sigmaCut := sigma[:idxBelow(cut, sigma)]

	rows, cols := m.Dims()
	t := tau(rows, cols, sigmaCut)
	rOpt := idxBelow(t, sigmaCut)

	err := plotValues(path, sigmaCut, t, f, rOpt, rFrac)

	return &Summary{Rows: rows, Cols: cols, OptimalRank: rOpt, FractionalRank: rFrac, Sigma: sigma}, err
}

func idxAbove(thresh float64, s []float64) int {
	for i, v := range s {
		if v > thresh {
			return i
		}
	}
	return len(s)
}

func idxBelow(thresh float64, s []float64) int {
	for i, v := range s {
		if v < thresh {
			return i
		}
	}
	return len(s)
}

// https://arxiv.org/abs/1305.5870 Eq. 4.
func tau(rows, cols int, values []float64) float64 {
	if len(values) == 0 {
		return 0
	}
	reverseFloats(values)
	m := stat.Quantile(0.5, 1, values, nil)
	reverseFloats(values)
	return omega(rows, cols) * m
}

func reverseFloats(f []float64) {
	for i, j := 0, len(f)-1; i < j; i, j = i+1, j-1 {
		f[i], f[j] = f[j], f[i]
	}
}

// https://arxiv.org/abs/1305.5870 Eq. 5.
func omega(rows, cols int) float64 {
	beta := float64(rows) / float64(cols)
	beta2 := beta * beta
	return 0.56*beta2*beta - 0.95*beta2 + 1.82*beta + 1.43
}
