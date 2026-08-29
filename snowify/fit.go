// Stats-driven parameter fitting, per the 2026-08-29 protocol: the
// statistics drive the iteration against the real winter training data;
// the user's eye is reserved for the bigger gates.
//
// During fitting the snow cover is forced to the REAL winter snow mask, so
// the optimizer matches pure appearance over the same geography — the
// noise-field cover generator (the s knob) stays out of the loop.
package main

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
)

// params are the tunable appearance constants of the snowify pipeline. The
// defaults are the hand-calibrated v0.3 values.
type params struct {
	TexMid, TexWidth float64 // open-land texture gate: below TexMid is field
	LumMid, LumWidth float64 // open-land luminance gate: above LumMid is land
	LatticeDiv       float64 // top-hat DN that makes a fully-kept dark line
	PaintP50         float64 // snow luminance at median shading
	PaintExp         float64 // relief-shading exponent
	PaintMin         float64 // darkest and brightest paint
	PaintMax         float64
	PsfSigma         float64 // boundary blur, px

	// Green-share forest gate (green/(r+g+b) above GreenMid = forest, no
	// snow). Off at the default GreenMid 1. S2's luminance+texture gates
	// separate forest from field on reflectance data; rendered orthos
	// compress that contrast and need the hue axis instead.
	GreenMid, GreenWidth float64
}

// The fit-2 winner (2026-08-29, -fit 600 with widened bounds, objective
// 0.161 against the three 31TDL winter chips), rounded. The v0.3
// hand-calibrated values it replaced: TexMid 100/18, LumMid 260/35,
// LatticeDiv 150, Paint 5600^0.8 [3800,8500], PsfSigma 0.5.
var defaultParams = params{
	TexMid: 256, TexWidth: 25,
	LumMid: 339, LumWidth: 61,
	LatticeDiv: 274,
	PaintP50:   6800, PaintExp: 2.37, PaintMin: 4480, PaintMax: 8400,
	PsfSigma: 0.73,
	GreenMid: 1, GreenWidth: 0.02,
}

// The search space. Scale-like dimensions sample log-uniform.
type dim struct {
	name   string
	lo, hi float64
	log    bool
	get    func(*params) *float64
}

var dims = []dim{
	{"TexMid", 60, 260, false, func(p *params) *float64 { return &p.TexMid }},
	{"TexWidth", 8, 40, true, func(p *params) *float64 { return &p.TexWidth }},
	{"LumMid", 180, 350, false, func(p *params) *float64 { return &p.LumMid }},
	{"LumWidth", 15, 70, true, func(p *params) *float64 { return &p.LumWidth }},
	{"LatticeDiv", 60, 400, true, func(p *params) *float64 { return &p.LatticeDiv }},
	{"PaintP50", 4200, 8500, false, func(p *params) *float64 { return &p.PaintP50 }},
	{"PaintExp", 0.4, 2.4, false, func(p *params) *float64 { return &p.PaintExp }},
	{"PaintMin", 3000, 4500, false, func(p *params) *float64 { return &p.PaintMin }},
	{"PaintMax", 7500, 9000, false, func(p *params) *float64 { return &p.PaintMax }},
	{"PsfSigma", 0.3, 1.2, true, func(p *params) *float64 { return &p.PsfSigma }},
}

// fixed enforces the paint ordering PaintMin < PaintP50 < PaintMax.
func (p params) fixed() params {
	p.PaintMin = math.Min(p.PaintMin, p.PaintP50-200)
	p.PaintMax = math.Max(p.PaintMax, p.PaintP50+500)
	return p
}

func sampleUniform(rng *rand.Rand) params {
	p := defaultParams
	for _, d := range dims {
		if d.log {
			*d.get(&p) = d.lo * math.Exp(rng.Float64()*math.Log(d.hi/d.lo))
		} else {
			*d.get(&p) = d.lo + rng.Float64()*(d.hi-d.lo)
		}
	}
	return p.fixed()
}

// sampleAround perturbs each dimension of best by sd of its range.
func sampleAround(best params, sd float64, rng *rand.Rand) params {
	p := best
	for _, d := range dims {
		v := d.get(&p)
		if d.log {
			*v *= math.Exp(rng.NormFloat64() * sd * math.Log(d.hi/d.lo))
		} else {
			*v += rng.NormFloat64() * sd * (d.hi - d.lo)
		}
		*v = math.Min(math.Max(*v, d.lo), d.hi)
	}
	return p.fixed()
}

// fit random-searches the space (first half uniform, second half refining
// around the best), streaming every improvement.
func fit(ctxs []*chipCtx, iters int) params {
	rng := rand.New(rand.NewPCG(7, 7))
	best := defaultParams
	bestObj := objective(ctxs, best)
	fmt.Printf("fit start obj %.4f %+v\n", bestObj, best)
	for i := 0; i < iters; i++ {
		var p params
		if i < iters/2 {
			p = sampleUniform(rng)
		} else {
			p = sampleAround(best, 0.12, rng)
		}
		if obj := objective(ctxs, p); obj < bestObj {
			bestObj, best = obj, p
			fmt.Printf("iter %d obj %.4f %+v\n", i, obj, p)
		}
	}
	fmt.Printf("fit done obj %.4f %+v\n", bestObj, best)
	return best
}

// objective is the summed squared log-ratio between the snowified and the
// real winter stat vectors, all chips, chips evaluated in parallel.
func objective(ctxs []*chipCtx, p params) float64 {
	var mu sync.Mutex
	var wg sync.WaitGroup
	total := 0.0
	for _, c := range ctxs {
		wg.Add(1)
		go func(c *chipCtx) {
			defer wg.Done()
			s := statVec(c.snowify(p, c.fitCover), c)
			d := 0.0
			for i := range s {
				r := math.Log((s[i] + 1e-3) / (c.ref[i] + 1e-3))
				d += r * r
			}
			mu.Lock()
			total += d
			mu.Unlock()
		}(c)
	}
	wg.Wait()
	return total
}

// ---- the stat vector ---------------------------------------------------

var statNames = []string{
	"lum_p10", "lum_p50", "lum_p90", // in the real snow mask
	"gradmask_p50", "gradmask_p90", // sharpness inside the mask
	"gradchip_p50", "gradchip_p90", // sharpness chip-wide
	"white_frac", // fraction of the mask brighter than 7000 DN
}

const nStats = 8

func statVec(rgb [3]*img, c *chipCtx) [nStats]float64 {
	lum := luminance(rgb[0], rgb[1], rgb[2])
	gm := gradMag(lum)
	var hl, hgm, hgc hist
	white, nmask := 0, 0
	for i, m := range c.mask {
		hgc.add(gm.px[i])
		if m {
			hl.add(lum.px[i])
			hgm.add(gm.px[i])
			nmask++
			if lum.px[i] > 7000 {
				white++
			}
		}
	}
	if nmask == 0 {
		nmask = 1
	}
	return [nStats]float64{
		hl.pct(0.1), hl.pct(0.5), hl.pct(0.9),
		hgm.pct(0.5), hgm.pct(0.9),
		hgc.pct(0.5), hgc.pct(0.9),
		float64(white) / float64(nmask),
	}
}

// hist gives percentiles in O(n) with 4 DN resolution — the fit loop calls
// this thousands of times, sorting would dominate.
type hist struct {
	bins [4096]int32
	n    int
}

func (h *hist) add(v float32) {
	i := int(v) >> 2
	if i < 0 {
		i = 0
	} else if i > 4095 {
		i = 4095
	}
	h.bins[i]++
	h.n++
}

func (h *hist) pct(p float64) float64 {
	target := int(p * float64(h.n))
	acc := 0
	for i, b := range h.bins {
		acc += int(b)
		if acc >= target {
			return float64(i*4 + 2)
		}
	}
	return 4095 * 4
}
