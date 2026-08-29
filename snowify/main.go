// Snowify — classical image-space snow synthesis on Sentinel-2 chips.
//
// JackFrost's first committed tool. The python probes (jackfrost_pull*.py,
// jackfrost_stats2.py on macbook1) established the calibration: over tile
// 31TDL, snow-covered land is ~9x brighter than clear land and carries ~6x
// the gradient energy — snow is not a whiteout, it is bright relief crossed
// by a dark lattice of hedgerows, treelines, streams and roads.
//
// The pipeline: steep texture/luminance gates confine snow to open land
// (forests stay dark), a morphological top-hat lets the dark lattice
// survive, the paint follows the terrain's low-frequency shading, partial
// cover (the s knob) sweeps a smooth noise field so patches keep sharp
// boundaries, and a sub-pixel blur on the coverage map reproduces the
// sensor's point-spread at the snow edges.
//
// Working protocol (user, 2026-08-29): the statistics drive the tuning —
// see fit.go — with the user's eye reserved for bigger gates. Run with
// -fit N to refit the appearance constants against the real winter
// captures; without it, params come from -params (or the defaults).
//
// Self-validation: chips of 31TDL, summer S2B_20250807 snowified and
// panelled against the real winter capture S2A_20240120 of the same
// (MGRS-aligned) tile: summer | s=0.5 | s=1.0 | real-cover | real winter.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"image"
	"image/png"
	"math"
	"os"
)

const chipSize = 1536 // 10 m px -> 15.4 km chips

// Snow is faintly cold: a touch less red, a touch more blue than neutral.
var snowTint = [3]float32{0.97, 1.0, 1.06}

// The three validation chips, chosen by the v0 probe for high winter snow
// cover with a mixed landscape, both scenes cloud-free. 10 m grid coords.
var chips = [][2]int{{8363, 380}, {5701, 1900}, {5321, 3801}} // x, y

var (
	rawDir     = flag.String("raw", "raw", "directory with the pulled S2 band tifs")
	outDir     = flag.String("out", "snowify_out", "output directory")
	summer     = flag.String("summer", "S2B_31TDL_20250807_0_L2A", "summer scene id")
	winter     = flag.String("winter", "S2A_31TDL_20240120_0_L2A", "winter scene id")
	fitIters   = flag.Int("fit", 0, "refit the appearance params with this many search iterations")
	paramsFile = flag.String("params", "", "json file with appearance params (default: built-in)")
)

func main() {
	flag.Parse()
	if err := os.MkdirAll(*outDir, 0755); err != nil {
		fatal(err)
	}

	p := defaultParams
	if *paramsFile != "" {
		b, err := os.ReadFile(*paramsFile)
		if err != nil {
			fatal(err)
		}
		if err := json.Unmarshal(b, &p); err != nil {
			fatal(err)
		}
		fmt.Println("params from", *paramsFile)
	}

	if *framesDir != "" {
		runFrames(p)
		return
	}

	if *autoChips {
		pickChips()
	}

	// Load one band at a time (a full 10 m band is 480 MB as float32) and
	// keep only the chip crops.
	sum := map[string][]*img{}
	wtr := map[string][]*img{}
	for _, band := range []string{"red", "green", "blue"} {
		sum[band] = mustChips("summer", *summer, band)
		wtr[band] = mustChips("winter", *winter, band)
	}
	wtr["swir16"] = mustChips("winter", *winter, "swir16")

	ctxs := make([]*chipCtx, len(chips))
	for n := range chips {
		ctxs[n] = prep(n,
			[3]*img{sum["red"][n], sum["green"][n], sum["blue"][n]},
			[3]*img{wtr["red"][n], wtr["green"][n], wtr["blue"][n]},
			wtr["swir16"][n])
		fmt.Printf("chip%d ready, winter ref %v\n", n, ctxs[n].ref)
	}

	fmt.Printf("objective with current params: %.4f\n", objective(ctxs, p))
	if *fitIters > 0 {
		p = fit(ctxs, *fitIters)
		writeJSON(*outDir+"/fitted.json", p)
	}

	grades := map[string]map[string][nStats]float64{}
	for _, c := range ctxs {
		half := c.snowify(p, c.noiseCover(0.5))
		full := c.snowify(p, c.noiseCover(1.0))
		realcov := c.snowify(p, c.fitCover)

		name := fmt.Sprintf("%s/chip%d_summer_s05_s10_realcov_winter.png", *outDir, c.id)
		writePanel(name, [][3]*img{c.sum, half, full, realcov, c.win})

		g := map[string][nStats]float64{
			"summer":  statVec(c.sum, c),
			"s10":     statVec(full, c),
			"realcov": statVec(realcov, c),
			"winter":  c.ref,
		}
		grades[fmt.Sprintf("chip%d", c.id)] = g
		enc, _ := json.MarshalIndent(g, "", " ")
		fmt.Printf("chip%d stats %v\n%s\n", c.id, statNames, enc)
	}
	writeJSON(*outDir+"/grades.json", grades)
	fmt.Println("JACKFROST-SNOWIFY-DONE")
}

// ---- the per-chip pipeline ---------------------------------------------

// chipCtx holds a chip's inputs and the parameter-independent intermediate
// rasters, so the fit loop only redoes the cheap per-pixel work.
type chipCtx struct {
	id       int
	sum, win [3]*img
	lum      *img // summer luminance
	texture  *img // smoothed local gradient energy: forest vs field
	green    *img // green share g/(r+g+b): the forest hue axis
	tophat   *img // morphological closing minus lum: thin dark features
	shade    *img // low-frequency relief shading
	shMed    float64
	noise    *img   // smooth field the s knob sweeps through
	mask     []bool // real winter snow (NDSI > 0.4)
	fitCover *img   // the real mask as a cover map, for fitting
	ref      [nStats]float64
}

func prep(id int, sum, win [3]*img, wswir *img) *chipCtx {
	c := prepLite(id, sum, 1) // S2 chips are the calibration scale
	c.win = win

	c.mask = make([]bool, len(c.lum.px))
	c.fitCover = newImg(c.lum.w, c.lum.h)
	wg, sw := win[1], wswir
	for i := range c.mask {
		ndsi := (wg.px[i] - sw.px[i]) / max32(wg.px[i]+sw.px[i], 1)
		if ndsi > 0.4 {
			c.mask[i] = true
			c.fitCover.px[i] = 1
		}
	}
	c.ref = statVec(win, c)
	return c
}

// coverP combines the open-land gates, the dark lattice and the cover map
// into the per-pixel snow probability, softened by the sensor point-spread.
func (c *chipCtx) coverP(P params, cover *img) *img {
	p := newImg(c.lum.w, c.lum.h)
	for i := range p.px {
		open := sigf((float32(P.TexMid)-c.texture.px[i])/float32(P.TexWidth)) *
			sigf((c.lum.px[i]-float32(P.LumMid))/float32(P.LumWidth)) *
			sigf((float32(P.GreenMid)-c.green.px[i])/float32(P.GreenWidth))
		lattice := clamp01(c.tophat.px[i] / float32(P.LatticeDiv))
		p.px[i] = open * cover.px[i] * (1 - lattice)
	}
	return gauss(p, P.PsfSigma)
}

// paintImg is the snow luminance: the calibrated level modulated by the
// terrain's own low-frequency shading (relative to its median, so it is
// radiometry-independent).
func (c *chipCtx) paintImg(P params) *img {
	paint := newImg(c.lum.w, c.lum.h)
	for i := range paint.px {
		l := P.PaintP50 * math.Pow(float64(c.shade.px[i])/c.shMed, P.PaintExp)
		paint.px[i] = float32(math.Min(math.Max(l, P.PaintMin), P.PaintMax))
	}
	return paint
}

// snowify applies snow with the given cover map and appearance params.
func (c *chipCtx) snowify(P params, cover *img) [3]*img {
	p := c.coverP(P, cover)
	paint := c.paintImg(P)
	var out [3]*img
	for ch, src := range c.sum {
		out[ch] = newImg(src.w, src.h)
		for i := range out[ch].px {
			out[ch].px[i] = (1-p.px[i])*src.px[i] + p.px[i]*paint.px[i]*snowTint[ch]
		}
	}
	return out
}

// noiseCover is the deployment cover: degree s in [0,1] sweeps a threshold
// through the smooth noise field, so partial cover has patch boundaries.
func (c *chipCtx) noiseCover(s float64) *img {
	return each(c.noise, func(v float32) float32 {
		return sigf((float32(s) - v) / 0.06)
	})
}

func sigf(x float32) float32 {
	return float32(1 / (1 + math.Exp(-float64(x))))
}

// ---- raster loading ----------------------------------------------------

// mustChips reads one band file and returns the chip crops on the 10 m
// grid. S2 band resolutions differ: swir16/scl come as 20 m rasters and are
// upsampled 2x nearest; chip coordinates are halved accordingly.
func mustChips(tag, scene, band string) []*img {
	path := fmt.Sprintf("%s/%s_%s_%s.tif", *rawDir, tag, scene, band)
	full, err := loadBand(path)
	if err != nil {
		fatal(err)
	}
	scale := 10980 / full.w // 1 for 10 m bands, 2 for 20 m
	out := make([]*img, len(chips))
	for i, c := range chips {
		ch := crop(full, c[0]/scale, c[1]/scale, chipSize/scale)
		if scale == 2 {
			ch = upsample2(ch)
		}
		out[i] = ch
	}
	fmt.Println("loaded", path)
	return out
}

// ---- output ------------------------------------------------------------

// writePanel renders RGB panels side by side, fixed scale DN/8000, gamma 1.8
// — the same rendering for all panels so brightness comparison is honest.
func writePanel(path string, panels [][3]*img) {
	w, h := panels[0][0].w, panels[0][0].h
	out := image.NewNRGBA(image.Rect(0, 0, w*len(panels), h))
	for pi, p := range panels {
		for y := 0; y < h; y++ {
			for x := 0; x < w; x++ {
				o := out.PixOffset(pi*w+x, y)
				for c := 0; c < 3; c++ {
					v := math.Min(math.Max(float64(p[c].px[y*w+x])/8000, 0), 1)
					out.Pix[o+c] = uint8(math.Pow(v, 1/1.8) * 255)
				}
				out.Pix[o+3] = 255
			}
		}
	}
	f, err := os.Create(path)
	if err != nil {
		fatal(err)
	}
	defer f.Close()
	if err := png.Encode(f, out); err != nil {
		fatal(err)
	}
	fmt.Println("wrote", path)
}

func writeJSON(path string, v any) {
	f, err := os.Create(path)
	if err != nil {
		fatal(err)
	}
	defer f.Close()
	enc := json.NewEncoder(f)
	enc.SetIndent("", " ")
	if err := enc.Encode(v); err != nil {
		fatal(err)
	}
	fmt.Println("wrote", path)
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "snowify:", err)
	os.Exit(1)
}
