// Flight clouds (-weather cloud): scattered low cumulus over the frames of a
// flight, rendered through the calibrated lens with the full attitude of
// each pose.
//
// Theory of operation.
//
// The world is a local plane: x east, y north, metres from the middle of the
// flight (equirectangular at that latitude; at each camera the east axis is
// rescaled to the camera's latitude, so the geometry around a camera is true
// to well under a metre). Heights are in the DEM's datum, which must be that
// of the poses.
//
// The cloud layer is a volume between a base and a top that follow the
// LOCAL ground: the terrain smoothed with a Gaussian of -cloud-smooth-km,
// plus -cloud-base and -cloud-top metres. Inside that band the clouds are a
// 2-D field D on a 10 m grid, anchored in the world and carried by the wind
// (constant per run): the field coordinate of world point p at flight time t
// is p - wind t. D is a sum of cumulus cells, each a cluster of a few smooth
// puffs (kernel (1 - d^2/r^2)^2), modulated by three octaves of value noise
// for ragged edges, with cell sizes drawn log-normally around -cloud-cell-km.
// Cells are added at random places until the requested share of the field is
// covered as seen from straight above (vertical optical depth 1 or more);
// that is the only use of -cloud-cover. From D follow, per column:
//
//	the top:     base + (top - base) * sqrt((D - thr) / 0.9), clipped: flat
//	             bases, domed tops, the biggest cells reaching the top; raised
//	             and lowered by up to a fifth by value noise at 45 and 17 m
//	             evaluated at the sample (the cauliflower), never above the top;
//	the density: -cloud-sigma (extinction, 1/m) * smoothstep(thr, thr + 0.7, D),
//	             so a cell is opaque in its core and thin at its edge;
//	             the bottom 25 m and top 20 m of each column fade in.
//
// Each output pixel is a ray of the calibrated lens (computed on a half
// resolution grid of nodes and interpolated, because nothing in the cloud is
// finer than the 10 m grid) marched from the camera through the band, front
// to back, with 4 m steps inside cloud and 12 m outside (the entry found to a
// metre by bisection, so thin edges carry no step grain), a per-node jitter
// against banding, and ending at the ground (terrain fixed point on the
// ray). Every in-cloud sample scatters light
// amb + (1 - amb) exp(-0.004 d_sun), d_sun the distance to the cloud top
// along the sun ray (one lookup at the column the sun ray leaves through),
// so sun-facing domes are white and the lee sides and deeper parts grey.
// The three flight cases fall out of the same march: below the base the
// rays never enter the band on the way down (the camera only looks down);
// inside the band, inside a cell, the first metres are cloud and the frame
// is white; above, the rays pass the tops and reach the ground in the gaps.
//
// Shadows: from the ground point of each node a ray is marched toward the
// sun (low-precision almanac at the pose's time, sun.go) through the band
// with the full extinction; the ground keeps amb_shadow + (1 - amb_shadow)
// exp(-tau) of its brightness (sky light only under an opaque cell),
// slightly bluer.
//
// The output pixel is ground * T * shadow + scattered light (premultiplied),
// on top of the input frame; nothing else of the input changes. Not
// modelled: cloud growth and decay during the flight, the camera's
// auto-exposure reacting to bright cloud, haze, light scattered under the
// clouds, cloud shadows on other clouds.
//
// The wind is -wind N,E (m/s), or the mean of the samples of -wind-log
// whose height above the field lies inside the layer: a JSON-lines file with
// h_field_m, wind_n_ms and wind_e_ms per line, e.g. the wind a flight
// simulator flew in.
package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"image"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

var (
	cloudCover   = flag.Float64("cloud-cover", 0.375, "flight cloud: share of the sky covered, seen from straight above (0.375 = 3/8)")
	cloudBase    = flag.Float64("cloud-base", 500, "flight cloud: cloud base above the local ground, m")
	cloudTop     = flag.Float64("cloud-top", 800, "flight cloud: highest cloud tops above the local ground, m")
	cloudSeed    = flag.Int64("cloud-seed", 1, "flight cloud: seed of the cloud field")
	cloudSmooth  = flag.Float64("cloud-smooth-km", 2, "flight cloud: the local ground is the terrain smoothed with this Gaussian sigma, km")
	cloudCellKm  = flag.Float64("cloud-cell-km", 0.6, "flight cloud: median cumulus cell size (diameter of its cluster), km")
	cloudSigma   = flag.Float64("cloud-sigma", 0.06, "flight cloud: extinction in a cell's core, 1/m")
	cloudMargin  = flag.Float64("cloud-margin-km", 9, "flight cloud: world kept around the track, km")
	cloudFieldM  = flag.Float64("cloud-field-m", 10, "flight cloud: cloud field grid, m")
	shadowAmb    = flag.Float64("shadow-amb", 0.42, "flight cloud: share of the light left under an opaque cell")
	shadows      = flag.Bool("shadows", true, "flight cloud: cloud shadows on the ground")
	windSpec     = flag.String("wind", "", "flight cloud: N,E wind the cloud moves with, m/s")
	windLog      = flag.String("wind-log", "", "flight cloud: JSON lines with h_field_m, wind_n_ms, wind_e_ms: the wind is their mean inside the layer (instead of -wind)")
	timeStart    = flag.String("t0", "", "flight cloud: UTC time (RFC 3339) of flight time 0, for the sun, when the poses carry no time_utc")
	cloudStats   = flag.String("cloud-stats", "", "flight cloud: append per-frame cloud statistics to this csv")
	cloudMaskDir = flag.String("cloud-masks", "", "flight cloud: also write the half-resolution cloud opacity here")
	cloudProfile = flag.Bool("profile", false, "flight cloud: print the layer under the flight and the phases, render nothing")
	cloudProbe   = flag.String("probe", "", "flight cloud: idx,idx,...: print what the principal point's ray meets (csv), render nothing")
	cloudMotion  = flag.String("motion", "", "flight cloud: A,B[;A,B]: check the cloud's motion between frames A and B against the geometry, render nothing")
)

// grid is a raster on the world plane with bilinear lookup, clamped at the
// edges.
type grid struct {
	x0, y0, step float64
	nx, ny       int
	v            []float32
}

func newGrid(x0, y0, x1, y1, step float64) *grid {
	nx := int(math.Ceil((x1-x0)/step)) + 1
	ny := int(math.Ceil((y1-y0)/step)) + 1
	return &grid{x0: x0, y0: y0, step: step, nx: nx, ny: ny, v: make([]float32, nx*ny)}
}

func (g *grid) xy(i, j int) (float64, float64) {
	return g.x0 + float64(i)*g.step, g.y0 + float64(j)*g.step
}

func (g *grid) at(x, y float64) float64 {
	gx, gy := (x-g.x0)/g.step, (y-g.y0)/g.step
	ix, iy := int(math.Floor(gx)), int(math.Floor(gy))
	if ix < 0 {
		ix, gx = 0, 0
	} else if ix > g.nx-2 {
		ix, gx = g.nx-2, float64(g.nx-1)
	}
	if iy < 0 {
		iy, gy = 0, 0
	} else if iy > g.ny-2 {
		iy, gy = g.ny-2, float64(g.ny-1)
	}
	fx, fy := gx-float64(ix), gy-float64(iy)
	k := iy*g.nx + ix
	a := float64(g.v[k])*(1-fx) + float64(g.v[k+1])*fx
	b := float64(g.v[k+g.nx])*(1-fx) + float64(g.v[k+g.nx+1])*fx
	return a*(1-fy) + b*fy
}

// world is the local plane.
type world struct {
	lat0, lon0, mN0, mE0 float64
}

func (w *world) xy(lat, lon float64) (float64, float64) {
	return (lon - w.lon0) * w.mE0, (lat - w.lat0) * w.mN0
}

func (w *world) latlon(x, y float64) (float64, float64) {
	return w.lat0 + y/w.mN0, w.lon0 + x/w.mE0
}

// value noise: bilinear-interpolated lattice hash with a smooth fade, in
// world metres.
func hash2(ix, iy int64, seed int64) float64 {
	h := uint64(ix)*0x9E3779B97F4A7C15 ^ uint64(iy)*0xC2B2AE3D27D4EB4F ^ uint64(seed)*0x165667B19E3779F9
	h ^= h >> 29
	h *= 0xBF58476D1CE4E5B9
	h ^= h >> 32
	return float64(h&0xFFFFFF) / float64(0xFFFFFF)
}

func smooth(t float64) float64 { return t * t * (3 - 2*t) }

func vnoise(x, y, cell float64, seed int64) float64 {
	gx, gy := x/cell, y/cell
	ix, iy := math.Floor(gx), math.Floor(gy)
	fx, fy := smooth(gx-ix), smooth(gy-iy)
	i, j := int64(ix), int64(iy)
	v00 := hash2(i, j, seed)
	v10 := hash2(i+1, j, seed)
	v01 := hash2(i, j+1, seed)
	v11 := hash2(i+1, j+1, seed)
	return (v00*(1-fx)+v10*fx)*(1-fy) + (v01*(1-fx)+v11*fx)*fy
}

// cloudNoise is three octaves of vnoise at cell, cell/2.7 and cell/7.1.
func cloudNoise(x, y, cell float64, seed int64) float64 {
	return 0.55*vnoise(x, y, cell, seed) +
		0.30*vnoise(x, y, cell/2.7, seed+101) +
		0.15*vnoise(x, y, cell/7.1, seed+202)
}

func smoothstep(a, b, x float64) float64 {
	t := (x - a) / (b - a)
	if t < 0 {
		t = 0
	}
	if t > 1 {
		t = 1
	}
	return t * t * (3 - 2*t)
}

// cloudShape holds the functions of D that make the column.
type cloudShape struct {
	thr, depth, sigma float64
	seed              int64
}

func (c *cloudShape) hrel(D float64) float64 {
	f := (D - c.thr) / 0.9
	if f <= 0 {
		return 0
	}
	if f > 1 {
		f = 1
	}
	return c.depth * math.Sqrt(f)
}

func (c *cloudShape) dens(D float64) float64 {
	return c.sigma * smoothstep(c.thr, c.thr+0.7, D)
}

// top is hrel with the cauliflower: the dome raised and lowered by up to
// a fifth of its height by value noise at 45 and 17 m, evaluated at the
// field point (fx, fy) because the 10 m grid cannot hold it.
func (c *cloudShape) top(D, fx, fy float64) float64 {
	hr := c.hrel(D)
	if hr <= 0 {
		return 0
	}
	n := 0.6*vnoise(fx, fy, 45, c.seed+51) + 0.4*vnoise(fx, fy, 17, c.seed+53)
	return math.Min(c.depth, hr*(1+0.4*(n-0.5))) // the highest tops stay at the top
}

// rho is the extinction at height zr above the base in a column of field
// value D at field point (fx, fy), 0 outside the cloud; and the column's top.
func (c *cloudShape) rho(D, zr, fx, fy float64) (float64, float64) {
	if D <= c.thr || zr < 0 {
		return 0, 0
	}
	hr := c.top(D, fx, fy)
	if zr > hr {
		return 0, hr
	}
	return c.dens(D) * smoothstep(0, 25, zr) * smoothstep(0, 20, hr-zr), hr
}

type renderer struct {
	w       world
	lens    *Pinhole
	terr    *grid // ground, 30 m
	base    *grid // cloud base, 100 m
	field   *grid // D, field coordinates
	shape   cloudShape
	windE   float64
	windN   float64
	dc      float64 // D where the vertical optical depth reaches 1
	t0Unix  float64
	gw, gh  int // node grid
	rays    []vec3
	jit     []float64
	ctNodes []int // nodes inside the 128 px box at the principal point
}

const amb = 0.5 // share of a cloud sample's light that does not depend on the sun path
// sunExt is the effective extinction for the sun's light inside cloud, far
// below the droplets' (multiple scattering carries light forward): a lit
// top stays white, the lee of a 200 m dome is grey.
const sunExt = 0.004
const cloudWhite = 250 // a sunlit top, in frame counts

func newRenderer(poses []pose, lens *Pinhole, dem *demSet) *renderer {
	r := &renderer{lens: lens}
	r.shape = cloudShape{thr: 0.30, depth: *cloudTop - *cloudBase, sigma: *cloudSigma, seed: *cloudSeed}
	if r.shape.depth <= 50 {
		log.Fatalf("-cloud-top %.0f must be well above -cloud-base %.0f", *cloudTop, *cloudBase)
	}
	ts := poses[0].TimeUTC
	if *timeStart != "" {
		ts = *timeStart
	}
	t0, err := time.Parse(time.RFC3339Nano, ts)
	if err != nil {
		log.Fatalf("the time of pose 0 %q (time_utc, or -t0): %v", ts, err)
	}
	r.t0Unix = float64(t0.UnixNano()) / 1e9
	if *timeStart == "" {
		r.t0Unix -= poses[0].T
	}

	// the plane, centred on the flight
	var sLat, sLon float64
	for _, p := range poses {
		sLat += p.Lat
		sLon += p.Lon
	}
	r.w.lat0, r.w.lon0 = sLat/float64(len(poses)), sLon/float64(len(poses))
	r.w.mN0, r.w.mE0, _ = metresPerDeg(r.w.lat0)
	dx0, dy0, dx1, dy1 := math.Inf(1), math.Inf(1), math.Inf(-1), math.Inf(-1)
	for _, p := range poses {
		x, y := r.w.xy(p.Lat, p.Lon)
		dx0, dy0 = math.Min(dx0, x), math.Min(dy0, y)
		dx1, dy1 = math.Max(dx1, x), math.Max(dy1, y)
	}
	m := *cloudMargin * 1000
	dx0, dy0, dx1, dy1 = dx0-m, dy0-m, dx1+m, dy1+m

	// wind at the layer's height above the field
	r.windN, r.windE = cloudWind()
	tMax := poses[len(poses)-1].T

	// terrain and the layer's base
	r.terr = newGrid(dx0, dy0, dx1, dy1, 30)
	holes := 0
	for j := 0; j < r.terr.ny; j++ {
		for i := 0; i < r.terr.nx; i++ {
			lat, lon := r.w.latlon(r.terr.xy(i, j))
			// adjacent clips need not abut exactly: a node in the seam
			// takes its neighbour's height a few metres away
			for _, d := range [][2]float64{{0, 0}, {1e-4, 0}, {-1e-4, 0}, {0, 1e-4}, {0, -1e-4}, {1e-4, 1e-4}, {-1e-4, -1e-4}} {
				if dem.covers(lat+d[0], lon+d[1]) {
					lat, lon = lat+d[0], lon+d[1]
					break
				}
			}
			if !dem.covers(lat, lon) {
				holes++
			}
			r.terr.v[j*r.terr.nx+i] = float32(dem.elevation(lat, lon))
		}
	}
	if holes > 0 {
		log.Fatalf("%d terrain nodes outside the DEM %s: shrink -cloud-margin-km", holes, *demDir)
	}
	r.base = newGrid(dx0, dy0, dx1, dy1, 100)
	for j := 0; j < r.base.ny; j++ {
		for i := 0; i < r.base.nx; i++ {
			x, y := r.base.xy(i, j)
			r.base.v[j*r.base.nx+i] = float32(r.terr.at(x, y))
		}
	}
	gaussBlur(r.base, *cloudSmooth*1000/r.base.step)
	for k := range r.base.v {
		r.base.v[k] += float32(*cloudBase)
	}

	// the field covers the domain at every flight time
	fx0 := dx0 - math.Max(0, r.windE*tMax)
	fx1 := dx1 - math.Min(0, r.windE*tMax)
	fy0 := dy0 - math.Max(0, r.windN*tMax)
	fy1 := dy1 - math.Min(0, r.windN*tMax)
	r.field = newGrid(fx0, fy0, fx1, fy1, *cloudFieldM)
	r.dc = r.buildField()

	r.setupNodes()
	p0, pl := poses[0], poses[len(poses)-1]
	el0, az0 := sunPosition(r.t0Unix+p0.T, p0.Lat, p0.Lon)
	el1, az1 := sunPosition(r.t0Unix+pl.T, pl.Lat, pl.Lon)
	fmt.Printf("clouds: layer %.0f-%.0f m above the ground smoothed at %.1f km, cells ~%.2f km, cover %.3f from above, sigma %.3f/m, seed %d\n",
		*cloudBase, *cloudTop, *cloudSmooth, *cloudCellKm, *cloudCover, *cloudSigma, *cloudSeed)
	fmt.Printf("  wind %.2f N %.2f E m/s (%.1f m/s to %.0f deg), drift over the flight %.0f m\n",
		r.windN, r.windE, math.Hypot(r.windN, r.windE), math.Mod(deg(math.Atan2(r.windE, r.windN))+360, 360), math.Hypot(r.windN, r.windE)*tMax)
	fmt.Printf("  sun el %.2f az %.2f at the start, el %.2f az %.2f at the end; shadows %v\n", el0, az0, el1, az1, *shadows)
	fmt.Printf("  field %dx%d at %.0f m, D threshold %.2f, opaque from above at D %.3f\n", r.field.nx, r.field.ny, *cloudFieldM, r.shape.thr, r.dc)
	return r
}

// cloudDiagnostics prints the layer under the flight, and runs -probe or
// -motion; true when one of them (or -profile) means nothing is rendered.
func cloudDiagnostics(r *renderer, poses []pose) bool {
	r.profile(poses)
	if *cloudProfile {
		return true
	}
	if *cloudProbe != "" {
		want := map[int]bool{}
		for _, f := range strings.Split(*cloudProbe, ",") {
			k, err := strconv.Atoi(strings.TrimSpace(f))
			must(err)
			want[k] = true
		}
		r.probe(poses, want)
		return true
	}
	if *cloudMotion != "" {
		for _, pair := range strings.Split(*cloudMotion, ";") {
			var a, b int
			if _, err := fmt.Sscanf(pair, "%d,%d", &a, &b); err != nil {
				log.Fatalf("-motion %q: want A,B[;A,B...]", *cloudMotion)
			}
			r.motionCheck(poses, a, b)
		}
		return true
	}
	return false
}

// cloudWind is -wind "N,E" m/s, or the mean -wind-log wind over the samples
// whose height above the field lies inside the layer.
func cloudWind() (n, e float64) {
	if *windSpec != "" {
		parts := strings.Split(*windSpec, ",")
		if len(parts) != 2 {
			log.Fatalf("-wind %q: want N,E in m/s", *windSpec)
		}
		n, err1 := strconv.ParseFloat(strings.TrimSpace(parts[0]), 64)
		e, err2 := strconv.ParseFloat(strings.TrimSpace(parts[1]), 64)
		if err1 != nil || err2 != nil {
			log.Fatalf("-wind %q: want N,E in m/s", *windSpec)
		}
		return n, e
	}
	if *windLog == "" {
		log.Fatal("flight cloud: -wind N,E or -wind-log FILE required")
	}
	f, err := os.Open(*windLog)
	must(err)
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<20)
	var sn, se float64
	cnt := 0
	for sc.Scan() {
		var t struct {
			H  float64 `json:"h_field_m"`
			WN float64 `json:"wind_n_ms"`
			WE float64 `json:"wind_e_ms"`
		}
		if json.Unmarshal(sc.Bytes(), &t) != nil {
			continue
		}
		if t.H >= *cloudBase && t.H <= *cloudTop {
			sn += t.WN
			se += t.WE
			cnt++
		}
	}
	must(sc.Err())
	if cnt == 0 {
		log.Fatalf("%s: no sample between %.0f and %.0f m above the field", *windLog, *cloudBase, *cloudTop)
	}
	fmt.Printf("  wind from %d samples %.0f-%.0f m above the field\n", cnt, *cloudBase, *cloudTop)
	return sn / float64(cnt), se / float64(cnt)
}

// gaussBlur smooths a grid in place, separably, sigma in cells, edges
// clamped.
func gaussBlur(g *grid, sigma float64) {
	if sigma < 0.5 {
		return
	}
	rad := int(3 * sigma)
	k := make([]float64, 2*rad+1)
	var s float64
	for i := range k {
		d := float64(i - rad)
		k[i] = math.Exp(-d * d / (2 * sigma * sigma))
		s += k[i]
	}
	for i := range k {
		k[i] /= s
	}
	tmp := make([]float32, len(g.v))
	clampI := func(i, n int) int { return min(max(i, 0), n-1) }
	for j := 0; j < g.ny; j++ {
		for i := 0; i < g.nx; i++ {
			var a float64
			for q := -rad; q <= rad; q++ {
				a += k[q+rad] * float64(g.v[j*g.nx+clampI(i+q, g.nx)])
			}
			tmp[j*g.nx+i] = float32(a)
		}
	}
	for j := 0; j < g.ny; j++ {
		for i := 0; i < g.nx; i++ {
			var a float64
			for q := -rad; q <= rad; q++ {
				a += k[q+rad] * float64(tmp[clampI(j+q, g.ny)*g.nx+i])
			}
			g.v[j*g.nx+i] = float32(a)
		}
	}
}

// buildField stamps cumulus cells into r.field until -cloud-cover of it is
// opaque from above, and returns the D at which that happens.
func (r *renderer) buildField() float64 {
	sh := &r.shape
	// D where dens * hrel (the vertical optical depth, ignoring the faded
	// base and top) reaches 1
	lo, hi := sh.thr, sh.thr+3
	for i := 0; i < 50; i++ {
		mid := (lo + hi) / 2
		if sh.dens(mid)*sh.hrel(mid) < 1 {
			lo = mid
		} else {
			hi = mid
		}
	}
	dc := (lo + hi) / 2
	g := r.field
	seed := *cloudSeed
	// edge modulation: three octaves of world-anchored value noise
	mod := make([]float32, len(g.v))
	for j := 0; j < g.ny; j++ {
		for i := 0; i < g.nx; i++ {
			x, y := g.xy(i, j)
			n := 0.5*vnoise(x, y, 170, seed+11) + 0.3*vnoise(x, y, 65, seed+23) + 0.2*vnoise(x, y, 27, seed+37)
			mod[j*g.nx+i] = float32(0.45 + 1.1*n)
		}
	}
	rng := rand.New(rand.NewSource(seed))
	total := len(g.v)
	covered := 0
	target := int(*cloudCover * float64(total))
	medianR := *cloudCellKm * 1000 / 2
	nCells := 0
	fdc := float32(dc)
	for covered < target {
		cx := g.x0 + rng.Float64()*float64(g.nx-1)*g.step
		cy := g.y0 + rng.Float64()*float64(g.ny-1)*g.step
		R := medianR * math.Exp(0.45*rng.NormFloat64())
		R = math.Min(math.Max(R, 0.35*medianR), 2.3*medianR)
		nPuff := 3 + int(R/100)
		for q := 0; q < nPuff; q++ {
			a := rng.Float64() * 2 * math.Pi
			d := 0.55 * R * math.Sqrt(rng.Float64())
			px, py := cx+d*math.Cos(a), cy+d*math.Sin(a)
			pr := R * (0.45 + 0.35*rng.Float64())
			i0 := max(0, int((px-pr-g.x0)/g.step))
			i1 := min(g.nx-1, int((px+pr-g.x0)/g.step)+1)
			j0 := max(0, int((py-pr-g.y0)/g.step))
			j1 := min(g.ny-1, int((py+pr-g.y0)/g.step)+1)
			for j := j0; j <= j1; j++ {
				for i := i0; i <= i1; i++ {
					x, y := g.xy(i, j)
					q2 := ((x-px)*(x-px) + (y-py)*(y-py)) / (pr * pr)
					if q2 >= 1 {
						continue
					}
					k := j*g.nx + i
					old := g.v[k]
					g.v[k] += float32((1-q2)*(1-q2)) * mod[k]
					if old < fdc && g.v[k] >= fdc {
						covered++
					}
				}
			}
		}
		nCells++
	}
	fmt.Printf("  %d cumulus cells stamped, %.3f of the field opaque from above\n", nCells, float64(covered)/float64(total))
	return dc
}

// setupNodes computes the camera ray of every node of the half-resolution
// grid, its jitter, and which nodes fall in the 128 px box at the principal
// point.
func (r *renderer) setupNodes() {
	l := r.lens
	r.gw, r.gh = l.w/2, l.h/2
	r.rays = make([]vec3, r.gw*r.gh)
	r.jit = make([]float64, r.gw*r.gh)
	for j := 0; j < r.gh; j++ {
		for i := 0; i < r.gw; i++ {
			k := j*r.gw + i
			r.rays[k] = l.Ray(float64(2*i+1), float64(2*j+1))
			r.jit[k] = hash2(int64(i), int64(j), 977)
		}
	}
	for j := 0; j < r.gh; j++ {
		for i := 0; i < r.gw; i++ {
			if math.Abs(float64(2*i+1)-l.cx) <= 64 && math.Abs(float64(2*j+1)-l.cy) <= 64 {
				r.ctNodes = append(r.ctNodes, j*r.gw+i)
			}
		}
	}
}

type frameStats struct {
	phase                      string
	inCloud                    bool
	camAboveBase               float64
	frameCloud, frameOpacity   float64
	centreCloud, centreOpacity float64
	shadowFrac                 float64
}

// camState is what one pose needs from the world.
type camState struct {
	v              *view
	cx, cy, cz     float64
	scaleE         float64
	dfx, dfy       float64 // world -> field offset at this time
	sx, sy, sz     float64 // unit vector to the sun, ENU
	zbCam, terrCam float64
}

func (r *renderer) camAt(p pose) *camState {
	roll, pitch, yaw := p.attitude()
	c := &camState{v: newView(p.Lat, p.Lon, p.AltMSL, roll, pitch, yaw)}
	c.cx, c.cy = r.w.xy(p.Lat, p.Lon)
	c.cz = p.AltMSL
	c.scaleE = r.w.mE0 / c.v.mPerDegE
	c.dfx, c.dfy = -r.windE*p.T, -r.windN*p.T
	el, az := sunPosition(r.t0Unix+p.T, p.Lat, p.Lon)
	c.sx, c.sy, c.sz = math.Sin(rad(az))*math.Cos(rad(el)), math.Cos(rad(az))*math.Cos(rad(el)), math.Sin(rad(el))
	c.zbCam = r.base.at(c.cx, c.cy)
	c.terrCam = r.terr.at(c.cx, c.cy)
	return c
}

// phase is below / band / above the layer at the camera, and whether the
// camera is inside a cloud.
func (r *renderer) phase(c *camState) (string, bool) {
	zr := c.cz - c.zbCam
	switch {
	case zr < 0:
		return "below", false
	case zr > r.shape.depth:
		return "above", false
	}
	D := r.field.at(c.cx+c.dfx, c.cy+c.dfy)
	rho, _ := r.shape.rho(D, zr, c.cx+c.dfx, c.cy+c.dfy)
	return "band", rho > 0.2*r.shape.sigma
}

// ground is the terrain point on a ray (world-plane direction wx, wy, wz) by
// fixed point from the height under the camera.
func (r *renderer) ground(c *camState, wx, wy, wz float64) (x, y, z, t float64, ok bool) {
	if wz > -0.02 {
		return 0, 0, 0, 0, false
	}
	z = c.terrCam
	for it := 0; it < 8; it++ {
		t = (z - c.cz) / wz
		x, y = c.cx+t*wx, c.cy+t*wy
		z = r.terr.at(x, y)
	}
	return x, y, z, (z - c.cz) / wz, true
}

// march is the front-to-back integration of one ray through the band up to
// tEnd: transmittance and premultiplied scattered light (0..1).
// tHalf is where the transmittance first falls below one half, -1 never.
func (r *renderer) march(c *camState, wx, wy, wz, tEnd, jit float64) (T, L, tHalf float64) {
	tHalf = -1
	sh := &r.shape
	T = 1
	lo, hi := c.zbCam-300, c.zbCam+sh.depth+300
	t0, t1 := 0.0, math.Min(tEnd, 25000)
	switch {
	case wz < -1e-6:
		t0 = math.Max(0, (hi-c.cz)/wz)
		t1 = math.Min(t1, (lo-c.cz)/wz)
	case wz > 1e-6:
		t0 = math.Max(0, (lo-c.cz)/wz)
		t1 = math.Min(t1, (hi-c.cz)/wz)
	default:
		if c.cz < lo || c.cz > hi {
			return 1, 0, -1
		}
	}
	const coarse, fine = 12.0, 4.0
	tanEl := c.sz / math.Hypot(c.sx, c.sy)
	hx, hy := c.sx/math.Hypot(c.sx, c.sy), c.sy/math.Hypot(c.sx, c.sy)
	// at returns the extinction at ray parameter t and what the light needs
	type smp struct{ rho, hr, zb, z, fx, fy float64 }
	at := func(t float64) smp {
		x, y, z := c.cx+t*wx, c.cy+t*wy, c.cz+t*wz
		zb := r.base.at(x, y)
		fx, fy := x+c.dfx, y+c.dfy
		rho, hr := sh.rho(r.field.at(fx, fy), z-zb, fx, fy)
		return smp{rho, hr, zb, z, fx, fy}
	}
	// a jittered start only off the camera: from inside the band the first
	// metres may already be cloud, and skipping a random part of them
	// shows as grain
	t := t0 + jit*fine
	if t0 == 0 {
		t = 0
	}
	prevOut := -1.0 // ray parameter of the last sample outside cloud, -1 none
	for t < t1 && T > 0.004 {
		x, y, z := c.cx+t*wx, c.cy+t*wy, c.cz+t*wz
		zb := r.base.at(x, y)
		zr := z - zb
		if zr > sh.depth+2 {
			if wz >= 0 {
				break
			}
			t += math.Max(coarse, (zr-sh.depth-1)/(-wz)+jit*coarse)
			prevOut = -1
			continue
		}
		if zr < -2 {
			if wz <= 0 {
				break
			}
			t += math.Max(coarse, (-zr-1)/wz+jit*coarse)
			prevOut = -1
			continue
		}
		q := at(t)
		if q.rho <= 0 {
			prevOut = t
			t += coarse
			continue
		}
		if prevOut >= 0 {
			// entered cloud during a coarse step: find the entry to a
			// metre, or the entry depth (and so the opacity of thin
			// edges) would carry the step's quantisation as grain
			a, b := prevOut, t
			for i := 0; i < 4; i++ {
				m := (a + b) / 2
				if at(m).rho > 0 {
					b = m
				} else {
					a = m
				}
			}
			t = b
			q = at(t)
			prevOut = -1
		}
		a := 1 - math.Exp(-q.rho*fine)
		// light: distance to the top along the sun ray, from the column the
		// sun ray leaves through
		off := (q.zb + q.hr - q.z) / tanEl
		hr2 := sh.top(r.field.at(q.fx+hx*off, q.fy+hy*off), q.fx+hx*off, q.fy+hy*off)
		dz := math.Max(0, q.zb+hr2-q.z)
		Ls := amb + (1-amb)*math.Exp(-sunExt*dz/c.sz)
		L += T * a * Ls
		T *= 1 - a
		if T < 0.5 && tHalf < 0 {
			tHalf = t
		}
		t += fine
	}
	return T, L, tHalf
}

// shadow is the share of light a ground point keeps under the layer.
func (r *renderer) shadow(c *camState, xg, yg, zg, jit float64) float64 {
	sh := &r.shape
	zb := r.base.at(xg, yg)
	tA := math.Max(0, (zb-zg)/c.sz)
	tB := (zb + sh.depth - zg) / c.sz
	if tB <= 0 {
		return 1
	}
	const step = 5.0
	tau := 0.0
	for t := tA + jit*step; t < tB && tau < 7; t += step {
		x, y, z := xg+t*c.sx, yg+t*c.sy, zg+t*c.sz
		rho, _ := sh.rho(r.field.at(x+c.dfx, y+c.dfy), z-zb, x+c.dfx, y+c.dfy)
		tau += rho * step
	}
	return *shadowAmb + (1-*shadowAmb)*math.Exp(-tau)
}

// frame renders one pose onto its input frame; the mask (half resolution
// opacity) only with -cloud-masks.
func (r *renderer) frame(p pose, src *image.RGBA) (*image.RGBA, frameStats, *image.Gray) {
	c := r.camAt(p)
	var st frameStats
	st.phase, st.inCloud = r.phase(c)
	st.camAboveBase = c.cz - c.zbCam
	n := r.gw * r.gh
	Tn := make([]float32, n)
	Ln := make([]float32, n)
	Sn := make([]float32, n)
	for k := 0; k < n; k++ {
		w := c.v.fromCam(r.rays[k])
		wx, wy, wz := w[0]*c.scaleE, w[1], w[2]
		xg, yg, zg, tg, okG := r.ground(c, wx, wy, wz)
		tEnd := math.Inf(1)
		if okG {
			tEnd = tg
		}
		T, L, _ := r.march(c, wx, wy, wz, tEnd, r.jit[k])
		S := 1.0
		if okG && T > 0.004 && *shadows {
			S = r.shadow(c, xg, yg, zg, r.jit[k])
		}
		Tn[k], Ln[k], Sn[k] = float32(T), float32(L), float32(S)
	}
	// stats on the nodes
	var nShadeVis, nVis int
	for k := 0; k < n; k++ {
		if Tn[k] < 0.5 {
			st.frameCloud++
		} else {
			nVis++
			if Sn[k] < float32((1+*shadowAmb)/2) {
				nShadeVis++
			}
		}
		st.frameOpacity += 1 - float64(Tn[k])
	}
	st.frameCloud /= float64(n)
	st.frameOpacity /= float64(n)
	if nVis > 0 {
		st.shadowFrac = float64(nShadeVis) / float64(nVis)
	}
	for _, k := range r.ctNodes {
		if Tn[k] < 0.5 {
			st.centreCloud++
		}
		st.centreOpacity += 1 - float64(Tn[k])
	}
	st.centreCloud /= float64(max(len(r.ctNodes), 1))
	st.centreOpacity /= float64(max(len(r.ctNodes), 1))

	// compose at full resolution: node (i, j) sits at pixel edge
	// coordinates (2i+1, 2j+1), so pixel (px, py), centre (px+.5, py+.5),
	// is at node coordinates ((px-.5)/2, (py-.5)/2)
	W, H := r.lens.w, r.lens.h
	dst := image.NewRGBA(image.Rect(0, 0, W, H))
	sb := src.Bounds()
	for py := 0; py < H; py++ {
		gy := math.Max(0, math.Min(float64(r.gh-1), (float64(py)-0.5)/2))
		j0 := min(int(gy), r.gh-2)
		fy := gy - float64(j0)
		for px := 0; px < W; px++ {
			gx := math.Max(0, math.Min(float64(r.gw-1), (float64(px)-0.5)/2))
			i0 := min(int(gx), r.gw-2)
			fx := gx - float64(i0)
			k := j0*r.gw + i0
			bil := func(a []float32) float64 {
				return (float64(a[k])*(1-fx)+float64(a[k+1])*fx)*(1-fy) +
					(float64(a[k+r.gw])*(1-fx)+float64(a[k+r.gw+1])*fx)*fy
			}
			T, L, S := bil(Tn), bil(Ln), bil(Sn)
			so := src.PixOffset(sb.Min.X+px, sb.Min.Y+py)
			sp := src.Pix[so : so+3 : so+3]
			gT := T * S
			bT := T * (S + 0.08*(1-S)) // shadow is lit by the blue sky
			cl := L * cloudWhite
			do := dst.PixOffset(px, py)
			dst.Pix[do+0] = clamp8(float64(sp[0])*gT + cl*0.985)
			dst.Pix[do+1] = clamp8(float64(sp[1])*gT + cl*0.995)
			dst.Pix[do+2] = clamp8(float64(sp[2])*bT + cl*1.02)
			dst.Pix[do+3] = 255
		}
	}
	var mask *image.Gray
	if *cloudMaskDir != "" {
		mask = image.NewGray(image.Rect(0, 0, r.gw, r.gh))
		for k := 0; k < n; k++ {
			mask.Pix[k] = clamp8(255 * (1 - float64(Tn[k])))
		}
	}
	return dst, st, mask
}

func clamp8(v float64) uint8 {
	if v <= 0 {
		return 0
	}
	if v >= 255 {
		return 255
	}
	return uint8(v + 0.5)
}

// profile prints the layer under the flight, the time in each phase and the
// cover realised over the flight's footprint (from above, opaque, sampled
// once a second over the ground within 45 deg of the nadir).
func (r *renderer) profile(poses []pose) {
	dt := 0.0
	if len(poses) > 1 {
		dt = poses[1].T - poses[0].T
	}
	type span struct {
		name     string
		from, to float64
	}
	var spans []span
	dur := map[string]float64{}
	inCloud := 0.0
	var covN, covHit int
	lastPrint := -1e9
	for k, p := range poses {
		c := r.camAt(p)
		ph, inc := r.phase(c)
		dur[ph] += dt
		if inc {
			inCloud += dt
		}
		if len(spans) == 0 || spans[len(spans)-1].name != ph {
			spans = append(spans, span{ph, p.T, p.T})
		}
		spans[len(spans)-1].to = p.T
		if p.T-lastPrint >= 30 || k == len(poses)-1 {
			lastPrint = p.T
			fmt.Printf("  t %6.1f  alt %6.0f  ground %6.0f  base %6.0f  top %6.0f  %s%s\n", p.T, c.cz, c.terrCam, c.zbCam, c.zbCam+r.shape.depth, ph,
				map[bool]string{true: " IN CLOUD", false: ""}[inc])
		}
		if k%15 != 0 {
			continue
		}
		rad := math.Max(100, c.cz-c.terrCam)
		for a := -10; a <= 10; a++ {
			for b := -10; b <= 10; b++ {
				x, y := c.cx+float64(a)/10*rad, c.cy+float64(b)/10*rad
				covN++
				if r.field.at(x+c.dfx, y+c.dfy) >= r.dc {
					covHit++
				}
			}
		}
	}
	for _, s := range spans {
		fmt.Printf("  phase %-5s  t %6.1f - %6.1f s\n", s.name, s.from, s.to)
	}
	fmt.Printf("  time below %.1f s, in the band %.1f s (in cloud %.1f s), above %.1f s\n", dur["below"], dur["band"], inCloud, dur["above"])
	fmt.Printf("  realised cover over the flight's footprint (opaque from above): %.3f\n", float64(covHit)/float64(max(covN, 1)))
}
