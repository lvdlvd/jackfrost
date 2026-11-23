// Open water (flight snow, -water). On Sentinel-2 open water is dark and the
// luminance gate keeps it free of snow. In summer orthophotos lakes and
// rivers can be bright turquoise, as bright as the fields, so the gate lets
// snow onto them; in mid-winter the larger lakes and the rivers stay open.
// Water areas (GeoJSON Polygons) are rasterized once onto a 5 m grid over the
// ground grid; the bilinearly sampled mask (a soft shoreline of a few metres)
// removes the snow there, leaving the water's own appearance.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"sort"
)

var waterPath = flag.String("water", "", "flight snow: water polygons GeoJSON; no snow inside")

const wcell = 5.0 // metres

type waterMask struct {
	g    *ground
	w, h int
	m    []float32
}

func loadWater(path string, g *ground) *waterMask {
	b, err := os.ReadFile(path)
	if err != nil {
		fatal(err)
	}
	var fc struct {
		Features []struct {
			Geometry struct {
				Type        string          `json:"type"`
				Coordinates json.RawMessage `json:"coordinates"`
			} `json:"geometry"`
		} `json:"features"`
	}
	if err := json.Unmarshal(b, &fc); err != nil {
		fatal(err)
	}
	wm := &waterMask{g: g, w: int(float64(g.w)*cell/wcell) + 1, h: int(float64(g.h)*cell/wcell) + 1}
	wm.m = make([]float32, wm.w*wm.h)
	n := 0
	for _, f := range fc.Features {
		if f.Geometry.Type != "Polygon" {
			continue
		}
		var rings [][][2]float64
		if json.Unmarshal(f.Geometry.Coordinates, &rings) != nil || len(rings) == 0 {
			continue
		}
		wm.fill(rings)
		n++
	}
	fmt.Printf("water: %d polygons\n", n)
	return wm
}

// fill rasterizes one polygon (outer ring and holes, even-odd) by
// scanlines through the cell centres.
func (wm *waterMask) fill(rings [][][2]float64) {
	type edge struct{ x0, y0, x1, y1 float64 }
	var edges []edge
	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, r := range rings {
		for i := 0; i+1 < len(r); i++ {
			e0, n0 := wm.g.toEN(r[i][1], r[i][0])
			e1, n1 := wm.g.toEN(r[i+1][1], r[i+1][0])
			x0, y0 := (e0-wm.g.e0)/wcell, (n0-wm.g.n0)/wcell
			x1, y1 := (e1-wm.g.e0)/wcell, (n1-wm.g.n0)/wcell
			edges = append(edges, edge{x0, y0, x1, y1})
			ymin, ymax = math.Min(ymin, math.Min(y0, y1)), math.Max(ymax, math.Max(y0, y1))
		}
	}
	var xs []float64
	for y := max(0, int(math.Floor(ymin))); y <= min(wm.h-1, int(math.Ceil(ymax))); y++ {
		yc := float64(y)
		xs = xs[:0]
		for _, e := range edges {
			if (e.y0 <= yc) != (e.y1 <= yc) {
				xs = append(xs, e.x0+(yc-e.y0)/(e.y1-e.y0)*(e.x1-e.x0))
			}
		}
		sort.Float64s(xs)
		for i := 0; i+1 < len(xs); i += 2 {
			for x := max(0, int(math.Ceil(xs[i]))); x <= min(wm.w-1, int(math.Floor(xs[i+1]))); x++ {
				wm.m[y*wm.w+x] = 1
			}
		}
	}
}

func (wm *waterMask) at(e, n float32) float32 {
	fx, fy := (float64(e)-wm.g.e0)/wcell, (float64(n)-wm.g.n0)/wcell
	x0, y0 := int(math.Floor(fx)), int(math.Floor(fy))
	tx, ty := float32(fx-float64(x0)), float32(fy-float64(y0))
	c := func(x, y int) float32 {
		if x < 0 || y < 0 || x >= wm.w || y >= wm.h {
			return 0
		}
		return wm.m[y*wm.w+x]
	}
	return (c(x0, y0)*(1-tx)+c(x0+1, y0)*tx)*(1-ty) + (c(x0, y0+1)*(1-tx)+c(x0+1, y0+1)*tx)*ty
}

// clear removes the snow on water and returns the frame's water fraction.
func (wm *waterMask) clear(pcov, pe, pn *img) float64 {
	var sum float64
	for i := range pcov.px {
		if v := wm.at(pe.px[i], pn.px[i]); v > 0 {
			pcov.px[i] *= 1 - v
			sum += float64(v)
		}
	}
	return sum / float64(len(pcov.px))
}
