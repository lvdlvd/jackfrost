// Rendered-frame mode: apply snowify to the vnavsim nadir-fisheye renders
// (8-bit PNGs) — the closer-to-ground path, and the input to the
// margin-under-snow matcher eval (snowified eval set vs its clear twin).
//
// One canonical parameter set, two radiometric domains. S2 L2A DNs are
// reflectance-scaled and DIM (summer fields ~550 DN of 8000); rendered
// orthos are display-bright. So the gates run in calibration space:
// linearize the frame (v^1.8), then scale so the frame's median luminance
// lands on the calibration's summer clear-land p50 — fields hit the field
// level, forests fall below the luminance gate, exactly as on S2. (A naive
// fixed linearization put snow on the forests and skipped the fields.)
// The BLEND then happens back in the frame's own display-linear space,
// with the snow's display brightness taken from the same DN/8000 mapping
// the validated S2 panels use.
//
// Scale: the calibration lives at S2's 10 m/px. -gsd rescales the
// physical-length constants (texture integration, lattice radius, relief
// shading, patch noise) to the frame's ground sample distance; the PSF
// blur stays in pixel units — it belongs to the sensor, not the terrain.
// Note the near-nadir GSD of the 480 px fisheye at typical corridor
// altitude is ~10 m/px, so the default is no rescaling at all. The warp-
// smeared frame edges read as high texture and are excluded by the
// ordinary open-land gate.
package main

import (
	"flag"
	"fmt"
	"image"
	"image/png"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sort"
)

var (
	framesDir = flag.String("frames", "", "directory of rendered 8-bit PNG frames to snowify (instead of S2 chips)")
	frameGSD  = flag.Float64("gsd", 10, "frame ground sample distance, m/px")
	degree    = flag.Float64("s", 1.0, "snow degree 0..1 for frame mode")
)

// summerCalLumP50 is the calibration's summer clear-land luminance median
// (jackfrost_stats2.json, chips 447-703) — the anchor frame radiometry is
// normalized to.
const summerCalLumP50 = 550.0

func runFrames(p params) {
	paths, err := filepath.Glob(filepath.Join(*framesDir, "*.png"))
	if err != nil || len(paths) == 0 {
		fatal(fmt.Errorf("no frames in %s", *framesDir))
	}
	sort.Strings(paths)
	k := 10 / *frameGSD
	for n, path := range paths {
		lin, err := loadFramePNG(path)
		if err != nil {
			fatal(err)
		}

		// gates in calibration space: scale display-linear -> DN so the
		// frame's median luminance sits at the summer field level
		var hl hist
		lum := luminance(lin[0], lin[1], lin[2])
		for _, v := range lum.px {
			hl.add(v * 4000) // spread 0..1 over the histogram's DN range
		}
		a := float32(summerCalLumP50 / (hl.pct(0.5) / 4000))
		var dn [3]*img
		for c := range dn {
			dn[c] = each(lin[c], func(v float32) float32 { return v * a })
		}

		ctx := prepLite(n, dn, k)
		pcov := ctx.coverP(p, ctx.noiseCover(*degree))
		paint := ctx.paintImg(p)

		// blend in the frame's display-linear space; the snow's display
		// brightness comes from the panels' validated DN/8000 mapping
		var out [3]*img
		for c := range out {
			out[c] = newImg(lum.w, lum.h)
			for i := range out[c].px {
				snow := paint.px[i] / 8000 * snowTint[c]
				out[c].px[i] = (1-pcov.px[i])*lin[c].px[i] + pcov.px[i]*snow
			}
		}
		dst := filepath.Join(*outDir, filepath.Base(path))
		if err := writeFramePNG(dst, out); err != nil {
			fatal(err)
		}
		fmt.Println("snowified", dst)
	}
	fmt.Println("JACKFROST-FRAMES-DONE")
}

// loadFramePNG reads an 8-bit frame into display-linear [0,1] channels.
func loadFramePNG(path string) ([3]*img, error) {
	f, err := os.Open(path)
	if err != nil {
		return [3]*img{}, err
	}
	defer f.Close()
	src, _, err := image.Decode(f)
	if err != nil {
		return [3]*img{}, fmt.Errorf("%s: %w", path, err)
	}
	b := src.Bounds()
	var rgb [3]*img
	for c := range rgb {
		rgb[c] = newImg(b.Dx(), b.Dy())
	}
	i := 0
	for y := b.Min.Y; y < b.Max.Y; y++ {
		for x := b.Min.X; x < b.Max.X; x++ {
			r, g, bb, _ := src.At(x, y).RGBA() // 16-bit
			for c, v := range [3]uint32{r, g, bb} {
				rgb[c].px[i] = float32(math.Pow(float64(v)/65535, 1.8))
			}
			i++
		}
	}
	return rgb, nil
}

// writeFramePNG writes display-linear [0,1] channels back to an 8-bit PNG.
func writeFramePNG(path string, rgb [3]*img) error {
	w, h := rgb[0].w, rgb[0].h
	out := image.NewNRGBA(image.Rect(0, 0, w, h))
	for i := 0; i < w*h; i++ {
		for c := 0; c < 3; c++ {
			v := math.Min(math.Max(float64(rgb[c].px[i]), 0), 1)
			out.Pix[i*4+c] = uint8(math.Pow(v, 1/1.8) * 255)
		}
		out.Pix[i*4+3] = 255
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return png.Encode(f, out)
}

// prepLite computes the parameter-independent rasters snowify needs, with
// the physical-length constants scaled by k = 10 m / GSD (k=1 for S2).
// The 1.5 px pre-smooth stays in pixel units: it suppresses sensor noise.
func prepLite(id int, sum [3]*img, k float64) *chipCtx {
	c := &chipCtx{id: id, sum: sum}
	c.lum = luminance(sum[0], sum[1], sum[2])
	c.green = newImg(c.lum.w, c.lum.h)
	for i := range c.green.px {
		c.green.px[i] = sum[1].px[i] /
			max32(sum[0].px[i]+sum[1].px[i]+sum[2].px[i], 1e-6)
	}
	c.texture = gauss(gradMag(gauss(c.lum, 1.5)), 4*k)
	r := int(math.Max(2, math.Round(2*k)))
	c.tophat = each2(erode(dilate(c.lum, r), r), c.lum,
		func(closed, l float32) float32 { return closed - l })
	c.shade = gauss(c.lum, 24*k)
	var hs hist
	for _, v := range c.shade.px {
		hs.add(v)
	}
	c.shMed = hs.pct(0.5)

	rng := rand.New(rand.NewPCG(1000+uint64(id), 0))
	c.noise = gauss(randn(c.lum.w, c.lum.h, rng), 24*k)
	normalize(c.noise, 0.5, 0.25)
	return c
}
