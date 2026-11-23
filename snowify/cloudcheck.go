// Checks of the cloud renderer's geometry (flight cloud, -motion and
// -probe): they render nothing.
package main

import (
	"fmt"
	"image"
	"math"
	"path/filepath"
	"sort"
)

// motionCheck (-motion A,B) tests that the cloud moves through the picture
// as the geometry says. For the nodes of frame A where the cloud is opaque it
// takes the 3-D point where the ray's transmittance fell to one half, lets
// it drift with the wind to the time of frame B and projects it into B's
// camera: the predicted cloud motion. The same for the ground points of A
// gives the predicted ground motion. Then it measures both, independently of
// the prediction, by block matching: cloud patches on the opacity of A and B,
// ground patches on the clean frames, each searched over a window that holds
// both the cloud and the ground hypotheses. Displacements are in full-frame
// pixels.
func (r *renderer) motionCheck(poses []pose, a, b int) {
	if a < 0 || b < 0 || a >= len(poses) || b >= len(poses) {
		fmt.Printf("motion: frames %d and %d not both in the poses\n", a, b)
		return
	}
	pa, pb := &poses[a], &poses[b]
	ca, cb := r.camAt(*pa), r.camAt(*pb)
	n := r.gw * r.gh
	opac := func(c *camState) ([]float64, []vec3, []vec3) {
		T := make([]float64, n)
		hit := make([]vec3, n)
		gnd := make([]vec3, n)
		for k := 0; k < n; k++ {
			w := c.v.fromCam(r.rays[k])
			wx, wy, wz := w[0]*c.scaleE, w[1], w[2]
			xg, yg, zg, tg, okG := r.ground(c, wx, wy, wz)
			tEnd := math.Inf(1)
			if okG {
				tEnd = tg
				gnd[k] = vec3{xg, yg, zg}
			} else {
				gnd[k] = vec3{math.NaN(), 0, 0}
			}
			t, _, th := r.march(c, wx, wy, wz, tEnd, r.jit[k])
			T[k] = 1 - t
			if th >= 0 {
				hit[k] = vec3{c.cx + th*wx, c.cy + th*wy, c.cz + th*wz}
			} else {
				hit[k] = vec3{math.NaN(), 0, 0}
			}
		}
		return T, hit, gnd
	}
	oa, hitA, gndA := opac(ca)
	ob, _, _ := opac(cb)
	project := func(c *camState, x, y, z float64) (float64, float64, bool) {
		e := (x - c.cx) / c.scaleE
		nn := y - c.cy
		px, py, ok := r.lens.Pixel(c.v.toCam(vec3{e, nn, z - c.cz}))
		return px, py, ok && px >= 0 && py >= 0 && px < float64(r.lens.w) && py < float64(r.lens.h)
	}
	grey := func(idx int) []float64 {
		im, err := decodeRGBA(filepath.Join(*framesDir, fmt.Sprintf(*framePat, idx)))
		must(err)
		return halfGrey(im, r.gw, r.gh)
	}
	ga, gb := grey(a), grey(b)

	const half = 10 // patch half size, nodes
	patchStd := func(f []float64, i, j int) float64 {
		var s, s2, m float64
		for q := -half; q <= half; q++ {
			for p := -half; p <= half; p++ {
				v := f[(j+q)*r.gw+i+p]
				s += v
				s2 += v * v
				m++
			}
		}
		return math.Sqrt(math.Max(0, s2/m-(s/m)*(s/m)))
	}
	inside := func(i, j, m int) bool { return i >= m && j >= m && i < r.gw-m && j < r.gh-m }
	// match searches B for A's patch at (i, j) over the window around
	// (ci, cj) of radius rad nodes and returns the best displacement
	match := func(fa, fb []float64, i, j int, ci, cj float64, rad int) (float64, float64, bool) {
		best, bi, bj := math.Inf(1), 0, 0
		for dj := -rad; dj <= rad; dj++ {
			for di := -rad; di <= rad; di++ {
				ti, tj := int(math.Round(ci))+di, int(math.Round(cj))+dj
				if !inside(ti, tj, half) {
					continue
				}
				var ssd float64
				for q := -half; q <= half && ssd < best; q += 2 {
					for p := -half; p <= half; p += 2 {
						d := fa[(j+q)*r.gw+i+p] - fb[(tj+q)*r.gw+ti+p]
						ssd += d * d
					}
				}
				if ssd < best {
					best, bi, bj = ssd, ti, tj
				}
			}
		}
		if math.IsInf(best, 1) {
			return 0, 0, false
		}
		return float64(bi - i), float64(bj - j), true
	}
	var cloudPred, cloudMeas, cloudErr, groundPred, groundMeas, groundErr [][2]float64
	for j := 40; j < r.gh-40; j += 12 {
		for i := 40; i < r.gw-40; i += 12 {
			k := j*r.gw + i
			// predicted motion of the ground point and of the cloud point
			g := gndA[k]
			if math.IsNaN(g[0]) {
				continue
			}
			gx, gy, okg := project(cb, g[0], g[1], g[2])
			if !okg {
				continue
			}
			pgi, pgj := (gx-1)/2, (gy-1)/2 // node coordinates in B
			if oa[k] > 0.5 && !math.IsNaN(hitA[k][0]) {
				h := hitA[k]
				fx, fy := h[0]+ca.dfx, h[1]+ca.dfy
				cx, cy, okc := project(cb, fx-cb.dfx, fy-cb.dfy, h[2])
				if !okc || patchStd(oa, i, j) < 0.08 || !inside(i, j, half) {
					continue
				}
				pci, pcj := (cx-1)/2, (cy-1)/2
				rad := int(math.Hypot(pci-pgi, pcj-pgj)/2) + 20
				di, dj, ok := match(oa, ob, i, j, (pci+pgi)/2, (pcj+pgj)/2, rad)
				if !ok {
					continue
				}
				cloudPred = append(cloudPred, [2]float64{2 * (pci - float64(i)), 2 * (pcj - float64(j))})
				cloudMeas = append(cloudMeas, [2]float64{2 * di, 2 * dj})
				cloudErr = append(cloudErr, [2]float64{2 * (di - (pci - float64(i))), 2 * (dj - (pcj - float64(j)))})
				groundPred = append(groundPred, [2]float64{2 * (pgi - float64(i)), 2 * (pgj - float64(j))})
				continue
			}
			if oa[k] < 0.05 && ob[k] < 0.05 && inside(i, j, half) && patchStd(ga, i, j) > 0.03 {
				di, dj, ok := match(ga, gb, i, j, pgi, pgj, 20)
				if !ok {
					continue
				}
				groundMeas = append(groundMeas, [2]float64{2 * di, 2 * dj})
				groundErr = append(groundErr, [2]float64{2 * (di - (pgi - float64(i))), 2 * (dj - (pgj - float64(j)))})
			}
		}
	}
	med := func(v [][2]float64) string {
		if len(v) == 0 {
			return "none"
		}
		xs, ys := make([]float64, len(v)), make([]float64, len(v))
		for k := range v {
			xs[k], ys[k] = v[k][0], v[k][1]
		}
		sort.Float64s(xs)
		sort.Float64s(ys)
		return fmt.Sprintf("(%+.1f, %+.1f) px, n %d", xs[len(xs)/2], ys[len(ys)/2], len(v))
	}
	norms := func(v [][2]float64) string {
		if len(v) == 0 {
			return "none"
		}
		d := make([]float64, len(v))
		for k := range v {
			d[k] = math.Hypot(v[k][0], v[k][1])
		}
		sort.Float64s(d)
		return fmt.Sprintf("p50 %.1f px, p90 %.1f px", d[len(d)/2], d[len(d)*9/10])
	}
	fmt.Printf("motion %d -> %d (%.2f s): wind %.2f N %.2f E m/s\n", a, b, pb.T-pa.T, r.windN, r.windE)
	fmt.Printf("  cloud nodes:  predicted %s, measured %s, |measured - predicted| %s\n", med(cloudPred), med(cloudMeas), norms(cloudErr))
	fmt.Printf("  the ground under those nodes would have moved %s\n", med(groundPred))
	fmt.Printf("  clear ground: measured on the clean frames %s, |measured - predicted| %s\n", med(groundMeas), norms(groundErr))
}

// halfGrey is the frame's luminance (0..1) at the node grid.
func halfGrey(im *image.RGBA, gw, gh int) []float64 {
	g := make([]float64, gw*gh)
	b := im.Bounds()
	for j := 0; j < gh; j++ {
		for i := 0; i < gw; i++ {
			var s float64
			for q := 0; q < 2; q++ {
				for p := 0; p < 2; p++ {
					o := im.PixOffset(b.Min.X+2*i+p, b.Min.Y+2*j+q)
					s += 0.299*float64(im.Pix[o]) + 0.587*float64(im.Pix[o+1]) + 0.114*float64(im.Pix[o+2])
				}
			}
			g[j*gw+i] = s / (4 * 255)
		}
	}
	return g
}

// probe (-probe idx,idx,...) prints, per pose, what the ray through the
// principal point meets: the ground height, the transmittance of the cloud
// in front of it, and the height of the point where that transmittance fell
// to one half (NaN where it never did): the height whose motion an optical
// flow measured at the principal point sees.
func (r *renderer) probe(poses []pose, list map[int]bool) {
	fmt.Println("index,t_s,cam_msl,ground_msl,centre_T,cloud_hit_msl")
	for n, p := range poses {
		if !list[n] {
			continue
		}
		c := r.camAt(p)
		w := c.v.fromCam(r.lens.Ray(r.lens.cx, r.lens.cy))
		wx, wy, wz := w[0]*c.scaleE, w[1], w[2]
		_, _, zg, tg, ok := r.ground(c, wx, wy, wz)
		if !ok {
			continue
		}
		T, _, th := r.march(c, wx, wy, wz, tg, 0)
		zh := math.NaN()
		if th >= 0 {
			zh = c.cz + th*wz
		}
		fmt.Printf("%d,%.3f,%.1f,%.1f,%.4f,%.1f\n", n, p.T, c.cz, zg, T, zh)
	}
}
