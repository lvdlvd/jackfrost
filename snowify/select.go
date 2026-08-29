// Automatic validation-chip selection, so any pulled tile pair works —
// the 31TDL chip coordinates in main.go were picked by exactly this
// procedure (in the v0 python probe) and stay hardcoded for
// reproducibility of the calibration record.
package main

import (
	"flag"
	"fmt"
	"sort"
)

var autoChips = flag.Bool("autochips", false, "pick validation chips from the winter snow/cloud masks instead of the built-in 31TDL coordinates")

// pickChips scans the winter scene at 20 m for chip windows that are snowy
// but mixed (0.15 < snow < 0.85) and near cloud-free in BOTH seasons, and
// keeps the three snowiest non-overlapping ones.
func pickChips() {
	wg2 := downsample2(mustLoad("winter", *winter, "green")) // 10 m -> 20 m
	wsw := mustLoad("winter", *winter, "swir16")
	wscl := mustLoad("winter", *winter, "scl")
	sscl := mustLoad("summer", *summer, "scl")

	n := len(wg2.px)
	snow := make([]bool, n)
	bad := make([]bool, n)
	for i := range snow {
		ndsi := (wg2.px[i] - wsw.px[i]) / max32(wg2.px[i]+wsw.px[i], 1)
		wc, sc := wscl.px[i], sscl.px[i]
		wcloud := wc >= 8 && wc <= 10
		snow[i] = ndsi > 0.4 && wc > 0 && !wcloud
		bad[i] = wcloud || wc == 0 || (sc >= 8 && sc <= 10) || sc == 0
	}

	var nv, ns, nb int
	for i := range snow {
		if wscl.px[i] > 0 {
			nv++
		}
		if snow[i] {
			ns++
		}
		if bad[i] {
			nb++
		}
	}
	fmt.Printf("tile fractions: valid %.2f snow %.2f bad %.2f\n",
		float64(nv)/float64(n), float64(ns)/float64(n), float64(nb)/float64(n))

	cc := chipSize / 2 // chip size on the 20 m grid
	stride := cc / 4
	type cand struct {
		frac float64
		x, y int
	}
	var cands []cand
	for y := 0; y+cc <= wg2.h; y += stride {
		for x := 0; x+cc <= wg2.w; x += stride {
			var ns, nb int
			for yy := y; yy < y+cc; yy++ {
				for xx := x; xx < x+cc; xx++ {
					if snow[yy*wg2.w+xx] {
						ns++
					}
					if bad[yy*wg2.w+xx] {
						nb++
					}
				}
			}
			sf := float64(ns) / float64(cc*cc)
			if float64(nb)/float64(cc*cc) < 0.02 && sf > 0.15 && sf < 0.85 {
				cands = append(cands, cand{sf, x, y})
			}
		}
	}
	if len(cands) == 0 {
		fatal(fmt.Errorf("no usable chip window (snowy-but-mixed, cloud-free both seasons)"))
	}
	sort.Slice(cands, func(i, j int) bool { return cands[i].frac > cands[j].frac })

	chips = chips[:0]
	for _, c := range cands {
		ok := true
		for _, p := range chips {
			if abs(c.x*2-p[0]) < chipSize && abs(c.y*2-p[1]) < chipSize {
				ok = false
				break
			}
		}
		if ok {
			chips = append(chips, [2]int{c.x * 2, c.y * 2}) // back to 10 m
			fmt.Printf("chip at %d,%d winter snow frac %.3f\n", c.x*2, c.y*2, c.frac)
		}
		if len(chips) == 3 {
			break
		}
	}
}

func mustLoad(tag, scene, band string) *img {
	im, err := loadBand(fmt.Sprintf("%s/%s_%s_%s.tif", *rawDir, tag, scene, band))
	if err != nil {
		fatal(err)
	}
	return im
}

func abs(v int) int {
	if v < 0 {
		return -v
	}
	return v
}
