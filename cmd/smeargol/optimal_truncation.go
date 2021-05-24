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

// https://arxiv.org/abs/1305.5870
func optimalTruncation(path string, m *mat.Dense, cut, frac float64) error {
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDThin)
	if !ok {
		return fmt.Errorf("could not factorise %q", path)
	}
	sigma := svd.Values(nil)

	sum := make([]float64, len(sigma))
	floats.CumSum(sum, sigma)
	floats.Scale(1/sum[len(sum)-1], sum)
	rFrac := idxAbove(frac, sum)
	var f float64
	switch {
	case rFrac < len(sigma):
		f = sigma[rFrac]
	case len(sigma) != 0:
		f = sigma[0]
	}

	sigma = sigma[:idxBelow(cut, sigma)]
	rows, cols := m.Dims()
	t := tau(rows, cols, sigma)
	rOpt := idxBelow(t, sigma)

	if false {
		fmt.Printf("%s: %d[:%d]\n\t%v\n\t%v\n", path, len(sigma), rOpt, sigma, sigma[:rOpt])
	}

	return plotValues(path, sigma, t, f, rOpt, rFrac)
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
