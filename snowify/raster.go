// Rasters: a minimal float32 image type, the filters snowify needs, and
// GeoTIFF band loading via github.com/lvdlvd/geotiff.
package main

import (
	"fmt"
	"math"
	"math/rand/v2"
	"os"

	"github.com/lvdlvd/geotiff"
)

type img struct {
	w, h int
	px   []float32
}

func newImg(w, h int) *img { return &img{w, h, make([]float32, w*h)} }

func loadBand(path string) (*img, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	ifd, err := geotiff.OpenTIFF(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	w, h := ifd.WxL()
	px, err := geotiff.ReadSamples(ifd)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return &img{w, h, px}, nil
}

func crop(src *img, x0, y0, size int) *img {
	dst := newImg(size, size)
	for y := 0; y < size; y++ {
		copy(dst.px[y*size:(y+1)*size], src.px[(y0+y)*src.w+x0:])
	}
	return dst
}

// downsample2 halves resolution by 2x2 averaging.
func downsample2(src *img) *img {
	dst := newImg(src.w/2, src.h/2)
	for y := 0; y < dst.h; y++ {
		for x := 0; x < dst.w; x++ {
			dst.px[y*dst.w+x] = (src.px[(2*y)*src.w+2*x] +
				src.px[(2*y)*src.w+2*x+1] +
				src.px[(2*y+1)*src.w+2*x] +
				src.px[(2*y+1)*src.w+2*x+1]) / 4
		}
	}
	return dst
}

func upsample2(src *img) *img {
	dst := newImg(src.w*2, src.h*2)
	for y := 0; y < dst.h; y++ {
		for x := 0; x < dst.w; x++ {
			dst.px[y*dst.w+x] = src.px[(y/2)*src.w+x/2]
		}
	}
	return dst
}

func luminance(r, g, b *img) *img {
	out := newImg(r.w, r.h)
	for i := range out.px {
		out.px[i] = 0.3*r.px[i] + 0.5*g.px[i] + 0.2*b.px[i]
	}
	return out
}

// gauss is a separable gaussian blur with edge clamping.
func gauss(src *img, sigma float64) *img {
	radius := int(math.Ceil(3 * sigma))
	kernel := make([]float32, 2*radius+1)
	var sum float32
	for i := range kernel {
		d := float64(i - radius)
		kernel[i] = float32(math.Exp(-d * d / (2 * sigma * sigma)))
		sum += kernel[i]
	}
	for i := range kernel {
		kernel[i] /= sum
	}
	tmp, dst := newImg(src.w, src.h), newImg(src.w, src.h)
	for y := 0; y < src.h; y++ { // horizontal pass
		for x := 0; x < src.w; x++ {
			var acc float32
			for k, kv := range kernel {
				acc += kv * src.px[y*src.w+clampi(x+k-radius, src.w)]
			}
			tmp.px[y*src.w+x] = acc
		}
	}
	for y := 0; y < src.h; y++ { // vertical pass
		for x := 0; x < src.w; x++ {
			var acc float32
			for k, kv := range kernel {
				acc += kv * tmp.px[clampi(y+k-radius, src.h)*src.w+x]
			}
			dst.px[y*src.w+x] = acc
		}
	}
	return dst
}

// gradMag is the central-difference gradient magnitude (matches np.gradient
// so the calibration numbers stay comparable).
func gradMag(src *img) *img {
	dst := newImg(src.w, src.h)
	for y := 0; y < src.h; y++ {
		for x := 0; x < src.w; x++ {
			gx := src.px[y*src.w+clampi(x+1, src.w)] - src.px[y*src.w+clampi(x-1, src.w)]
			gy := src.px[clampi(y+1, src.h)*src.w+x] - src.px[clampi(y-1, src.h)*src.w+x]
			dst.px[y*src.w+x] = float32(math.Hypot(float64(gx), float64(gy)) / 2)
		}
	}
	return dst
}

// dilate/erode are separable max/min filters of the given radius.
func dilate(src *img, r int) *img { return minmax(src, r, func(a, b float32) bool { return a > b }) }
func erode(src *img, r int) *img  { return minmax(src, r, func(a, b float32) bool { return a < b }) }

func minmax(src *img, r int, better func(a, b float32) bool) *img {
	tmp, dst := newImg(src.w, src.h), newImg(src.w, src.h)
	for y := 0; y < src.h; y++ {
		for x := 0; x < src.w; x++ {
			best := src.px[y*src.w+x]
			for d := -r; d <= r; d++ {
				if v := src.px[y*src.w+clampi(x+d, src.w)]; better(v, best) {
					best = v
				}
			}
			tmp.px[y*src.w+x] = best
		}
	}
	for y := 0; y < src.h; y++ {
		for x := 0; x < src.w; x++ {
			best := tmp.px[y*src.w+x]
			for d := -r; d <= r; d++ {
				if v := tmp.px[clampi(y+d, src.h)*src.w+x]; better(v, best) {
					best = v
				}
			}
			dst.px[y*src.w+x] = best
		}
	}
	return dst
}

func randn(w, h int, rng *rand.Rand) *img {
	out := newImg(w, h)
	for i := range out.px {
		out.px[i] = float32(rng.NormFloat64())
	}
	return out
}

// normalize rescales to the given mean and standard deviation.
func normalize(a *img, mean, std float32) {
	var m, v float64
	for _, x := range a.px {
		m += float64(x)
	}
	m /= float64(len(a.px))
	for _, x := range a.px {
		v += (float64(x) - m) * (float64(x) - m)
	}
	sd := math.Sqrt(v/float64(len(a.px))) + 1e-9
	for i, x := range a.px {
		a.px[i] = mean + std*float32((float64(x)-m)/sd)
	}
}

func each(a *img, f func(float32) float32) *img {
	out := newImg(a.w, a.h)
	for i, v := range a.px {
		out.px[i] = f(v)
	}
	return out
}

func each2(a, b *img, f func(x, y float32) float32) *img {
	out := newImg(a.w, a.h)
	for i := range a.px {
		out.px[i] = f(a.px[i], b.px[i])
	}
	return out
}

func clampi(i, n int) int {
	if i < 0 {
		return 0
	}
	if i >= n {
		return n - 1
	}
	return i
}

func clamp01(v float32) float32 {
	if v < 0 {
		return 0
	}
	if v > 1 {
		return 1
	}
	return v
}

func max32(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}
