// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"image/color"
	"math"
	"path/filepath"
	"strconv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// plot values plots the singular values to plots/path.png along with
// the optimal and user specified fraction thresholds.
func plotValues(path string, sigma []float64, tau, frac float64, rOpt, rFrac int) error {
	p := plot.New()
	p.Title.Text = fmt.Sprintf("Singular Values\n%s", path)
	p.Y.Scale = logScale{}
	p.Y.Tick.Marker = logTicks{}
	sigXYs := sliceToXYs(sigma)
	if len(sigXYs) != 0 {
		values, err := plotter.NewLine(sigXYs)
		if err != nil {
			return err
		}

		threshOpt, err := plotter.NewLine(plotter.XYs{{X: 0, Y: tau}, {X: sigXYs[len(sigXYs)-1].X, Y: tau}})
		if err != nil {
			return err
		}
		threshOpt.Color = color.RGBA{B: 255, A: 255}

		threshFrac, err := plotter.NewLine(plotter.XYs{{X: 0, Y: frac}, {X: sigXYs[len(sigXYs)-1].X, Y: frac}})
		if err != nil {
			return err
		}
		threshFrac.Color = color.RGBA{R: 255, A: 255}

		p.Add(values, threshOpt, threshFrac)
	}
	return p.Save(18*vg.Centimeter, 15*vg.Centimeter, filepath.Join("plots", path+".png"))
}

func sliceToXYs(s []float64) plotter.XYs {
	xy := make(plotter.XYs, len(s))
	for i, v := range s {
		if v == 0 {
			return xy[:i]
		}
		xy[i] = plotter.XY{X: float64(i), Y: v}
	}
	return xy
}

type logScale struct{}

func (logScale) Normalize(min, max, x float64) float64 {
	min = math.Max(min, 1e-16)
	max = math.Max(max, 1e-16)
	x = math.Max(x, 1e-16)
	logMin := math.Log(min)
	return (math.Log(x) - logMin) / (math.Log(max) - logMin)
}

type logTicks struct{ powers int }

func (t logTicks) Ticks(min, max float64) []plot.Tick {
	min = math.Max(min, 1e-16)
	max = math.Max(max, 1e-16)
	if t.powers < 1 {
		t.powers = 1
	}

	val := math.Pow10(int(math.Log10(min)))
	max = math.Pow10(int(math.Ceil(math.Log10(max))))
	var ticks []plot.Tick
	for val < max {
		for i := 1; i < 10; i++ {
			if i == 1 {
				ticks = append(ticks, plot.Tick{Value: val, Label: strconv.FormatFloat(val, 'e', 0, 64)})
			}
			if t.powers != 1 {
				break
			}
			ticks = append(ticks, plot.Tick{Value: val * float64(i)})
		}
		val *= math.Pow10(t.powers)
	}
	ticks = append(ticks, plot.Tick{Value: val, Label: strconv.FormatFloat(val, 'e', 0, 64)})

	return ticks
}
