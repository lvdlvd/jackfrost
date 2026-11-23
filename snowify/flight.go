// Flight mode (-poses): weather on the frames of a simulated flight. The
// frames are renders of the ground through a calibrated camera, one per pose
// of an aircraft; flight mode puts snow on the ground (this file), a layer of
// cumulus in the air (clouds.go), or both, so that the weather stays where it
// is in the world from one frame to the next:
//
//	snowify -poses poses.json -calib lens.json -dem demdir -frames clear -out winter \
//	        -weather snow -params params/winter_ortho.json [-roads roads.geojson -water water.geojson]
//	snowify -poses poses.json -calib lens.json -dem demdir -frames clear -out cloudy \
//	        -weather cloud -cloud-cover 0.375 -cloud-base 500 -cloud-top 800 -wind-log wind.jsonl
//	snowify ... -weather snow+cloud ...   # both, cloud over snow, in one pass
//
// Inputs. poses.json is {"poses": [...]}, one pose per frame, in frame
// order, each with lat_deg, lon_deg, alt_m_msl and heading_deg, optionally
// the body attitude body_roll_deg, body_pitch_deg, body_yaw_deg (else level
// at the heading), and for clouds t_s (flight time, s) and time_utc
// (RFC 3339, for the sun; the first pose's is enough). The frame of pose n
// is -frames/<-pattern with n>. -calib is the lens (lens.go), -dem a
// directory of GeoTIFF elevation tiles in the frames' height datum (dem.go).
//
// Theory of operation, snow.
//
// Frame mode treats every frame as an independent picture: one -gsd for the
// whole run, and a partial-cover noise field drawn per frame in image
// coordinates. On a 15 Hz flight that makes the snow patches jump from frame
// to frame, and frames from a low aircraft run from about 0.1 to 1 m/px, so
// no single -gsd is right either. Flight mode keeps the snow model as it is
// (open-land gates, top-hat lattice, relief paint, PSF, snow tint, the
// calibrated params, coverP and paintImg unchanged) and changes where its
// inputs come from:
//
//   - each pixel's ground point: the rays of every 8th pixel (lens.go) are
//     marched to the DEM and bisected onto it, and ground east/north/height
//     are interpolated bilinearly to every pixel;
//   - the frame's GSD is the median over those nodes of the square root of
//     the ground area of one pixel; k = 10 m / GSD scales the physical-length
//     constants exactly as frame mode's -gsd does (texture 4k px, top-hat
//     radius 2k px, relief shading 24k px); the PSF and the 1.5 px pre-smooth
//     stay in pixels;
//   - the cover noise is frame mode's recipe (unit normal noise, gaussian of
//     24 cells, normalized to mean 0.5 and std 0.25) on a 10 m ground grid
//     around the route: the calibration's own field at its own 10 m scale,
//     sampled at each pixel's ground point, so a snow patch lies at the same
//     place on the ground in every frame;
//   - the snow degree s depends on the ground height under the pixel: -s-low
//     up to -h-low, rising linearly to -s-high at -h-high;
//   - optionally (-params-high) a second calibrated parameter set is blended
//     in by the frame's median ground height between -ph-low and -ph-high,
//     per frame and smoothly, so there is no jump between frames;
//   - optionally (-roads) the asphalt of plowed roads stays clear, and
//     (-water) open water stays open: see roads.go and water.go.
//
// The large blurs (sigma up to thousands of px) run as three box passes,
// the usual box approximation of a gaussian; the top-hat's dilate and erode
// run as van Herk / Gil-Werman; sigmas under 3 px use the exact gauss.
//
// Both weathers: the snow is computed and encoded to 8 bits exactly as a
// snow-only run writes it, and the cloud is drawn over that, so snow+cloud
// equals a cloud run over the output of a snow run. Frames are independent
// given their pose, so -workers frames run in parallel; outputs are written
// under a temporary name and renamed, so an output that exists is complete
// and a rerun skips it.
package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"image"
	"image/draw"
	"image/png"
	"log"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

var (
	posesPath   = flag.String("poses", "", "flight mode: poses.json of the rendered frames (one pose per frame); see flight.go")
	weatherSpec = flag.String("weather", "snow", "flight mode: snow, cloud or snow+cloud")
	calibPath   = flag.String("calib", "", "flight mode: lens calibration the frames were rendered through")
	demDir      = flag.String("dem", "", "flight mode: directory of DEM GeoTIFFs (the frames' height datum, usually MSL)")
	framePat    = flag.String("pattern", "frame_%05d.png", "flight mode: frame file name for a pose index")
	fromIdx     = flag.Int("from", 0, "flight mode: first pose index")
	toIdx       = flag.Int("to", -1, "flight mode: last pose index (inclusive), -1 = all")
	stepIdx     = flag.Int("step", 1, "flight mode: every n-th pose from -from")
	listIdx     = flag.String("list", "", "flight mode: comma-separated pose indices instead of -from/-to")
	nWorkers    = flag.Int("workers", runtime.NumCPU(), "flight mode: frames in parallel")
	sLow        = flag.Float64("s-low", 0.85, "flight snow: snow degree on ground up to -h-low")
	sHigh       = flag.Float64("s-high", 1.0, "flight snow: snow degree on ground from -h-high")
	hLow        = flag.Float64("h-low", 600, "flight snow: ground height (m) where s starts rising")
	hHigh       = flag.Float64("h-high", 800, "flight snow: ground height (m) where s reaches -s-high")
	paramsHigh  = flag.String("params-high", "", "flight snow: second param set, blended in by the frame's median ground height")
	phLow       = flag.Float64("ph-low", 700, "flight snow: frame median ground height where -params-high starts")
	phHigh      = flag.Float64("ph-high", 900, "flight snow: frame median ground height where -params-high is all")
	calScale    = flag.Bool("calscale", true, "flight snow: compute texture and top-hat in calibration units (see prepFlight); false = frame mode's -gsd scaling")
	calFoot     = flag.Float64("cal-foot", 0.25, "flight snow: sigma of the luminance smoothing for the top-hat, in 10 m pixels")
	gateFoot    = flag.Float64("gate-foot", 0.25, "flight snow: sigma of a 10 m sensor's footprint, in 10 m pixels, for the luminance gate")
	aFixed      = flag.Float64("a-fixed", 0, "flight snow: fixed display-linear to calibration-DN scale for every frame (0 = per-frame median, as frame mode; -1 = the median of the per-frame values over every 150th frame of the flight); per frame, a town or forest frame shifts the gates and the snow line moves between frames")
	smoothMid   = flag.Float64("smooth-mid", 0.10, "flight snow: relative fine texture below which dark land is open land (0 = off); see prepFlight")
	smoothWidth = flag.Float64("smooth-width", 0.02, "flight snow: softness of -smooth-mid")
	debugRaw    = flag.Bool("debug-raw", false, "flight snow: with -debug, also dump the feature rasters as float32")
	debugDir    = flag.String("debug", "", "flight snow: also write the snow probability per frame here")
)

// cell is the ground grid's spacing: the calibration's 10 m.
const cell = 10.0

// node is the image grid step at which rays are cast to the ground.
const node = 8

type pose struct {
	Index    int      `json:"index"`
	T        float64  `json:"t_s"`
	TimeUTC  string   `json:"time_utc"`
	Lat      float64  `json:"lat_deg"`
	Lon      float64  `json:"lon_deg"`
	AltMSL   float64  `json:"alt_m_msl"`
	Heading  float64  `json:"heading_deg"`
	BodyRoll *float64 `json:"body_roll_deg"`
	BodyPit  *float64 `json:"body_pitch_deg"`
	BodyYaw  *float64 `json:"body_yaw_deg"`
}

// attitude is the body roll, pitch and yaw, level at the heading where the
// pose has none.
func (p pose) attitude() (roll, pitch, yaw float64) {
	roll, pitch, yaw = 0, 0, p.Heading
	if p.BodyRoll != nil {
		roll = *p.BodyRoll
	}
	if p.BodyPit != nil {
		pitch = *p.BodyPit
	}
	if p.BodyYaw != nil {
		yaw = *p.BodyYaw
	}
	return
}

// ground is the route's ground grid: heights from the DEM and the cover
// noise, both on a 10 m east/north grid about (lat0, lon0).
type ground struct {
	lat0, lon0     float64
	mLat, mLon     float64 // metres per degree
	e0, n0         float64 // south-west corner, metres
	w, h           int
	elev, noise    []float32
	minElev, maxEl float32
}

// snowSetup is what every frame's snow needs from the flight.
type snowSetup struct {
	pLow   params
	pHigh  *params
	lens   *Pinhole
	rays   []vec3
	nx, ny int
	g      *ground
	rs     *roadSet
	wm     *waterMask
}

func runFlight(pLow params) {
	var doSnow, doCloud bool
	switch *weatherSpec {
	case "snow":
		doSnow = true
	case "cloud":
		doCloud = true
	case "snow+cloud", "cloud+snow":
		doSnow, doCloud = true, true
	default:
		fatal(fmt.Errorf("-weather %q: want snow, cloud or snow+cloud", *weatherSpec))
	}
	if *calibPath == "" || *demDir == "" {
		fatal(fmt.Errorf("flight mode needs -calib and -dem"))
	}

	var pf struct {
		Poses []pose `json:"poses"`
	}
	b, err := os.ReadFile(*posesPath)
	if err != nil {
		fatal(err)
	}
	if err := json.Unmarshal(b, &pf); err != nil {
		fatal(err)
	}
	if len(pf.Poses) == 0 {
		fatal(fmt.Errorf("%s: no poses", *posesPath))
	}
	lens, err := LoadPinhole(*calibPath)
	if err != nil {
		fatal(err)
	}
	dem, err := loadDEM(*demDir)
	if err != nil {
		fatal(err)
	}

	var cr *renderer
	if doCloud {
		cr = newRenderer(pf.Poses, lens, dem)
		if cloudDiagnostics(cr, pf.Poses) {
			return // -profile, -probe, -motion: nothing rendered
		}
	}
	if *framesDir == "" {
		fatal(fmt.Errorf("flight mode needs -frames"))
	}
	if err := os.MkdirAll(*outDir, 0o755); err != nil {
		fatal(err)
	}
	var ss *snowSetup
	if doSnow {
		ss = newSnowSetup(pLow, pf.Poses, lens, dem)
	}

	var todo []int
	if *listIdx != "" {
		for _, s := range strings.Split(*listIdx, ",") {
			n, err := strconv.Atoi(strings.TrimSpace(s))
			if err != nil {
				fatal(err)
			}
			todo = append(todo, n)
		}
	} else {
		last := *toIdx
		if last < 0 || last >= len(pf.Poses) {
			last = len(pf.Poses) - 1
		}
		for n := *fromIdx; n <= last; n += max(*stepIdx, 1) {
			todo = append(todo, n)
		}
	}
	for _, n := range todo {
		if n < 0 || n >= len(pf.Poses) {
			fatal(fmt.Errorf("pose index %d out of range", n))
		}
	}

	var statsW *bufio.Writer
	if cr != nil && *cloudStats != "" {
		f, err := os.OpenFile(*cloudStats, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0o644)
		if err != nil {
			fatal(err)
		}
		defer f.Close()
		statsW = bufio.NewWriter(f)
		if st, _ := f.Stat(); st.Size() == 0 {
			fmt.Fprintln(statsW, "index,t_s,phase,in_cloud,cam_above_base_m,frame_cloud,frame_opacity,centre_cloud,centre_opacity,shadow_frac")
		}
		defer statsW.Flush()
	}
	if cr != nil && *cloudMaskDir != "" {
		if err := os.MkdirAll(*cloudMaskDir, 0o755); err != nil {
			fatal(err)
		}
	}

	jobs := make(chan int)
	var wg sync.WaitGroup
	var mu sync.Mutex
	nDone, nSkip := 0, 0
	start := time.Now()
	for w := 0; w < *nWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for n := range jobs {
				name := fmt.Sprintf(*framePat, n)
				dst := filepath.Join(*outDir, name)
				if _, err := os.Stat(dst); err == nil {
					mu.Lock()
					nSkip++
					mu.Unlock()
					continue
				}
				src, err := decodeRGBA(filepath.Join(*framesDir, name))
				if err != nil {
					fatal(err)
				}
				if src.Bounds().Dx() != lens.w || src.Bounds().Dy() != lens.h {
					fatal(fmt.Errorf("%s: %dx%d, the lens is %dx%d", name, src.Bounds().Dx(), src.Bounds().Dy(), lens.w, lens.h))
				}
				var lines []string
				if ss != nil {
					out, line := ss.frame(n, name, pf.Poses[n], linFromRGBA(src))
					src = encodeRGBA(out)
					lines = append(lines, line)
				}
				var st frameStats
				var mask *image.Gray
				if cr != nil {
					src, st, mask = cr.frame(pf.Poses[n], src)
					lines = append(lines, fmt.Sprintf("cloud %.3f", st.frameCloud))
					if mask != nil {
						writePNG(filepath.Join(*cloudMaskDir, name), mask)
					}
				}
				writePNG(dst, src)
				mu.Lock()
				nDone++
				fmt.Println(dst+":", strings.Join(lines, " | "))
				if statsW != nil {
					fmt.Fprintf(statsW, "%d,%.3f,%s,%v,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f\n", n, pf.Poses[n].T, st.phase, st.inCloud, st.camAboveBase,
						st.frameCloud, st.frameOpacity, st.centreCloud, st.centreOpacity, st.shadowFrac)
					if nDone%200 == 0 {
						statsW.Flush()
					}
				}
				mu.Unlock()
			}
		}()
	}
	for _, n := range todo {
		jobs <- n
	}
	close(jobs)
	wg.Wait()
	fmt.Printf("flight %s: %d frames -> %s (%d skipped as existing), %.1f s\n", *weatherSpec, nDone, *outDir, nSkip, time.Since(start).Seconds())
	fmt.Println("JACKFROST-FLIGHT-DONE")
}

func newSnowSetup(pLow params, poses []pose, lens *Pinhole, dem *demSet) *snowSetup {
	ss := &snowSetup{pLow: pLow, lens: lens}
	if *paramsHigh != "" {
		b, err := os.ReadFile(*paramsHigh)
		if err != nil {
			fatal(err)
		}
		ph := defaultParams
		if err := json.Unmarshal(b, &ph); err != nil {
			fatal(err)
		}
		ss.pHigh = &ph
		fmt.Println("params-high from", *paramsHigh)
	}
	ss.g = buildGround(poses, dem)
	fmt.Printf("ground grid %dx%d at %.0f m, heights %.0f..%.0f m\n", ss.g.w, ss.g.h, cell, ss.g.minElev, ss.g.maxEl)
	if *roadsPath != "" {
		ss.rs = loadRoads(*roadsPath, ss.g)
	}
	if *waterPath != "" {
		ss.wm = loadWater(*waterPath, ss.g)
	}

	// the camera-frame rays of the nodes are the same for every frame
	W, H := lens.Size()
	ss.nx, ss.ny = (W-1+node-1)/node+1, (H-1+node-1)/node+1
	ss.rays = make([]vec3, ss.nx*ss.ny)
	for j := 0; j < ss.ny; j++ {
		for i := 0; i < ss.nx; i++ {
			x, y := nodeAt(i, W), nodeAt(j, H)
			ss.rays[j*ss.nx+i] = lens.Ray(float64(x)+0.5, float64(y)+0.5)
		}
	}

	// -a-fixed -1: the flight's own scale, the median of the per-frame
	// values over every 150th frame of the whole flight (not of -from/-to or
	// -list, so that parts rendered separately agree).
	if *aFixed < 0 {
		var as []float64
		for n := 0; n < len(poses); n += 150 {
			src, err := decodeRGBA(filepath.Join(*framesDir, fmt.Sprintf(*framePat, n)))
			if err != nil {
				fatal(err)
			}
			lin := linFromRGBA(src)
			var hl hist
			for _, v := range luminance(lin[0], lin[1], lin[2]).px {
				hl.add(v * 4000)
			}
			as = append(as, summerCalLumP50/(hl.pct(0.5)/4000))
		}
		sort.Float64s(as)
		*aFixed = as[len(as)/2]
		fmt.Printf("a-fixed %.0f: the median of %d frames (%.0f..%.0f)\n", *aFixed, len(as), as[0], as[len(as)-1])
	}
	return ss
}

func nodeAt(i, size int) int {
	if x := i * node; x < size-1 {
		return x
	}
	return size - 1
}

// frame snowifies one frame (display-linear RGB) and returns it with its
// log line.
func (ss *snowSetup) frame(n int, name string, ps pose, lin [3]*img) ([3]*img, string) {
	g, nx, ny := ss.g, ss.nx, ss.ny
	w, h := lin[0].w, lin[0].h

	// gates in calibration space, as frame mode
	var hl hist
	lum := luminance(lin[0], lin[1], lin[2])
	for _, v := range lum.px {
		hl.add(v * 4000)
	}
	a := float32(summerCalLumP50 / (hl.pct(0.5) / 4000))
	if *aFixed > 0 {
		a = float32(*aFixed) // one scale for the flight: see the flag
	}
	var dn [3]*img
	for c := range dn {
		dn[c] = each(lin[c], func(v float32) float32 { return v * a })
	}

	// the camera in the ground grid's metres; its basis in ENU
	ce, cn := g.toEN(ps.Lat, ps.Lon)
	cu := ps.AltMSL
	if gz := float64(g.at(g.elev, ce, cn)); cu < gz+2 {
		cu = gz + 2 // on the runway the camera can sit under the DEM
	}
	xC, yC, zC := camBasisNadir(ps.Heading)
	if ps.BodyRoll != nil || ps.BodyPit != nil || ps.BodyYaw != nil {
		xC, yC, zC = camBasisBody(ps.attitude())
	}

	// ground points of the nodes
	ne, nn, nh := make([]float64, nx*ny), make([]float64, nx*ny), make([]float64, nx*ny)
	for k, r := range ss.rays {
		d := vec3{
			r[0]*xC[0] + r[1]*yC[0] + r[2]*zC[0],
			r[0]*xC[1] + r[1]*yC[1] + r[2]*zC[1],
			r[0]*xC[2] + r[1]*yC[2] + r[2]*zC[2],
		}
		ne[k], nn[k], nh[k] = g.intersect(ce, cn, cu, d)
	}

	// frame GSD: median ground length of a pixel over the nodes
	var gsds []float64
	for j := 0; j+1 < ny; j++ {
		for i := 0; i+1 < nx; i++ {
			k := j*nx + i
			dx := float64(nodeAt(i+1, w) - nodeAt(i, w))
			dy := float64(nodeAt(j+1, h) - nodeAt(j, h))
			ex, nxx := (ne[k+1]-ne[k])/dx, (nn[k+1]-nn[k])/dx
			ey, nyy := (ne[k+nx]-ne[k])/dy, (nn[k+nx]-nn[k])/dy
			gsds = append(gsds, math.Sqrt(math.Abs(ex*nyy-ey*nxx)))
		}
	}
	sort.Float64s(gsds)
	gsd := math.Max(gsds[len(gsds)/2], 0.1) // taxiing: keep the filters finite
	k := 10 / gsd

	// per-pixel ground point, by bilinear interpolation of the nodes
	pe, pn, phh := newImg(w, h), newImg(w, h), newImg(w, h)
	for y := 0; y < h; y++ {
		j := min(y/node, ny-2)
		y0, y1 := nodeAt(j, h), nodeAt(j+1, h)
		ty := float64(y-y0) / float64(y1-y0)
		for x := 0; x < w; x++ {
			i := min(x/node, nx-2)
			x0, x1 := nodeAt(i, w), nodeAt(i+1, w)
			tx := float64(x-x0) / float64(x1-x0)
			q := j*nx + i
			bl := func(v []float64) float32 {
				return float32((v[q]*(1-tx)+v[q+1]*tx)*(1-ty) + (v[q+nx]*(1-tx)+v[q+nx+1]*tx)*ty)
			}
			pe.px[y*w+x], pn.px[y*w+x], phh.px[y*w+x] = bl(ne), bl(nn), bl(nh)
		}
	}

	// frame parameters: the low set, or the blend by median ground height
	hs := append([]float64(nil), nh...)
	sort.Float64s(hs)
	medH := hs[len(hs)/2]
	p := ss.pLow
	wHigh := 0.0
	if ss.pHigh != nil {
		wHigh = clamp01f((medH - *phLow) / (*phHigh - *phLow))
		p = blendParams(ss.pLow, *ss.pHigh, wHigh)
	}

	ctx := prepFlight(n, dn, k, p)
	cover := newImg(w, h)
	for i := range cover.px {
		v := g.at(g.noise, float64(pe.px[i]), float64(pn.px[i]))
		s := *sLow + (*sHigh-*sLow)*clamp01f((float64(phh.px[i])-*hLow)/(*hHigh-*hLow))
		cover.px[i] = sigf((float32(s) - v) / 0.06)
	}
	pcov := ctx.coverP(p, cover)
	roadFrac := 0.0
	if ss.rs != nil {
		roadFrac = ss.rs.clear(pcov, lin, pe, pn, p.PsfSigma)
	}
	waterFrac := 0.0
	if ss.wm != nil {
		waterFrac = ss.wm.clear(pcov, pe, pn)
	}
	paint := ctx.paintImg(p)

	var out [3]*img
	for c := range out {
		out[c] = newImg(w, h)
		for i := range out[c].px {
			snow := paint.px[i] / 8000 * snowTint[c]
			out[c].px[i] = (1-pcov.px[i])*lin[c].px[i] + pcov.px[i]*snow
		}
	}
	if *debugDir != "" {
		writeGray(filepath.Join(*debugDir, strings.TrimSuffix(name, ".png")+"_pcov.png"), pcov)
		pr := func(a *img) string {
			v := append([]float32(nil), a.px...)
			sort.Slice(v, func(i, j int) bool { return v[i] < v[j] })
			q := func(f float64) float32 { return v[int(f*float64(len(v)-1))] }
			return fmt.Sprintf("p10 %.0f p50 %.0f p90 %.0f p99 %.0f", q(0.1), q(0.5), q(0.9), q(0.99))
		}
		fmt.Printf("debug %s: lum %s | texture %s | tophat %s | shade %s\n", name, pr(ctx.lum), pr(ctx.texture), pr(ctx.tophat), pr(ctx.shade))
		if *debugRaw {
			// raw float32 rasters for offline gate analysis; lat = lat0 +
			// pn/mLat, lon = lon0 + pe/mLon
			base := filepath.Join(*debugDir, strings.TrimSuffix(name, ".png"))
			texFine := gaussAny(gradMag(gauss(ctx.lum, 1.5)), 4*k)
			for _, r := range []struct {
				s string
				a *img
			}{{"lum", ctx.lum}, {"tex", ctx.texture}, {"texfine", texFine}, {"tophat", ctx.tophat}, {"green", ctx.green}, {"pe", pe}, {"pn", pn}} {
				writeRaw(base+"_"+r.s+".f32", r.a)
			}
			fmt.Printf("debugraw %s w %d h %d lat0 %.9f lon0 %.9f mLat %.6f mLon %.6f\n", name, w, h, g.lat0, g.lon0, g.mLat, g.mLon)
		}
	}
	var mc float64
	for _, v := range pcov.px {
		mc += float64(v)
	}
	return out, fmt.Sprintf("snow gsd %.3f k %.1f medH %.0f wHigh %.2f cover %.3f roads %.4f water %.4f", gsd, k, medH, wHigh, mc/float64(len(pcov.px)), roadFrac, waterFrac)
}

// prepFlight is prepLite with the fast large filters. Same rasters, same
// scaling by k, same 1.5 px pre-smooth.
func prepFlight(id int, sum [3]*img, k float64, P params) *chipCtx {
	c := &chipCtx{id: id, sum: sum}
	c.lum = luminance(sum[0], sum[1], sum[2])
	c.green = newImg(c.lum.w, c.lum.h)
	for i := range c.green.px {
		c.green.px[i] = sum[1].px[i] /
			max32(sum[0].px[i]+sum[1].px[i]+sum[2].px[i], 1e-6)
	}
	c.shade = gaussAny(c.lum, 24*k)
	r := int(math.Max(2, math.Round(2*k)))
	if !*calScale {
		c.texture = gaussAny(gradMag(gauss(c.lum, 1.5)), 4*k)
		closed := erodeFast(dilateFast(c.lum, r), r)
		c.tophat = each2(closed, c.lum, func(cl, l float32) float32 { return cl - l })
	} else {
		// features in calibration units: the pre-smooth in metres too,
		// the gradient per 10 m, the top-hat of the luminance as a 10 m
		// sensor would see it (footprint sigma -cal-foot px at 10 m)
		kk := float32(k)
		c.texture = gaussAny(each(gradMag(gaussAny(c.lum, 1.5*k)), func(v float32) float32 { return v * kk }), 4*k)
		lumCal := gaussAny(c.lum, *calFoot*k)
		closed := erodeFast(dilateFast(lumCal, r), r)
		c.tophat = each2(closed, lumCal, func(cl, l float32) float32 { return cl - l })
		// the land-cover gates decide per 10 m footprint, not per pixel:
		// per pixel, in-field texture near LumMid speckles the snow.
		// coverP reads only the gates from c.lum; shade is done above.
		lumGate := gaussAny(c.lum, *gateFoot*k)
		if *smoothMid > 0 {
			// Dark but smooth is a dark crop field, not forest: at these
			// resolutions the luminance gate alone cannot tell them apart
			// (both 170-260 DN), tree crowns can. Relative fine texture,
			// all lengths in metres: gradient of the luminance pre-smoothed
			// 0.75 m, per metre, over the 5 m local mean, averaged over
			// 10 m. Fitted against OpenStreetMap forest polygons on four
			// frames of a lowland flight: at 0.10, 99 % of dark forest
			// stays dark and 55 % of dark non-forest (crops; the rest is
			// tree lines, gardens) is snowed. The gate becomes
			// 1-(1-lumgate)(1-smooth), handed to coverP as the luminance
			// that gives that gate.
			gsd := 10 / k
			pre := gaussAny(c.lum, math.Max(0.7, 0.75/gsd))
			gm := gradMag(pre)
			loc := gaussAny(c.lum, 5/gsd)
			rel := each2(gm, loc, func(g, l float32) float32 { return g / float32(gsd) / (l + 20) })
			tex := gaussAny(rel, 10/gsd)
			for i := range lumGate.px {
				gt := 1 / (1 + math.Exp(-float64(lumGate.px[i]-float32(P.LumMid))/P.LumWidth))
				sm := 1 / (1 + math.Exp(-(*smoothMid-float64(tex.px[i]))/(*smoothWidth)))
				g2 := math.Min(math.Max(1-(1-gt)*(1-sm), 1e-6), 1-1e-6)
				lumGate.px[i] = float32(P.LumMid + P.LumWidth*math.Log(g2/(1-g2)))
			}
		}
		c.lum = lumGate
	}
	var hs hist
	for _, v := range c.shade.px {
		hs.add(v)
	}
	c.shMed = hs.pct(0.5)
	return c
}

func blendParams(a, b params, t float64) params {
	l := func(x, y float64) float64 { return x + (y-x)*t }
	return params{
		TexMid: l(a.TexMid, b.TexMid), TexWidth: l(a.TexWidth, b.TexWidth),
		LumMid: l(a.LumMid, b.LumMid), LumWidth: l(a.LumWidth, b.LumWidth),
		LatticeDiv: l(a.LatticeDiv, b.LatticeDiv),
		PaintP50:   l(a.PaintP50, b.PaintP50), PaintExp: l(a.PaintExp, b.PaintExp),
		PaintMin: l(a.PaintMin, b.PaintMin), PaintMax: l(a.PaintMax, b.PaintMax),
		PsfSigma: l(a.PsfSigma, b.PsfSigma),
		GreenMid: l(a.GreenMid, b.GreenMid), GreenWidth: l(a.GreenWidth, b.GreenWidth),
	}
}

func clamp01f(v float64) float64 { return math.Max(0, math.Min(1, v)) }

// ---- the ground grid ------------------------------------------------------

func buildGround(poses []pose, dem *demSet) *ground {
	minLat, maxLat, minLon, maxLon := 90.0, -90.0, 180.0, -180.0
	for _, p := range poses {
		minLat, maxLat = math.Min(minLat, p.Lat), math.Max(maxLat, p.Lat)
		minLon, maxLon = math.Min(minLon, p.Lon), math.Max(maxLon, p.Lon)
	}
	g := &ground{lat0: (minLat + maxLat) / 2, lon0: (minLon + maxLon) / 2}
	g.mLat = 6371000 * math.Pi / 180
	g.mLon = g.mLat * math.Cos(g.lat0*math.Pi/180)
	const margin = 5000.0 // a frame reaches about 2 km from the aircraft at cruise
	e1, n1 := g.toEN(minLat, minLon)
	e2, n2 := g.toEN(maxLat, maxLon)
	g.e0, g.n0 = e1-margin, n1-margin
	g.w, g.h = int((e2-e1+2*margin)/cell)+1, int((n2-n1+2*margin)/cell)+1

	g.elev = make([]float32, g.w*g.h)
	g.minElev, g.maxEl = 1e9, -1e9
	for j := 0; j < g.h; j++ {
		for i := 0; i < g.w; i++ {
			lat := g.lat0 + (g.n0+float64(j)*cell)/g.mLat
			lon := g.lon0 + (g.e0+float64(i)*cell)/g.mLon
			v := float32(dem.elevation(lat, lon)) // 0 outside the DEM: the grid margin only
			g.elev[j*g.w+i] = v
			g.minElev, g.maxEl = min(g.minElev, v), max(g.maxEl, v)
		}
	}

	// frame mode's noise recipe at the calibration's 10 m: gauss(randn, 24),
	// normalized to mean 0.5 and std 0.25; one fixed seed for the route
	rng := rand.New(rand.NewPCG(1000, 0))
	nz := randn(g.w, g.h, rng)
	nz = gauss(nz, 24)
	normalize(nz, 0.5, 0.25)
	g.noise = nz.px
	return g
}

func (g *ground) toEN(lat, lon float64) (e, n float64) {
	return (lon - g.lon0) * g.mLon, (lat - g.lat0) * g.mLat
}

// at samples a grid bilinearly at ground metres (e, n), clamped at the edge.
func (g *ground) at(v []float32, e, n float64) float32 {
	fx, fy := (e-g.e0)/cell, (n-g.n0)/cell
	x0, y0 := int(math.Floor(fx)), int(math.Floor(fy))
	tx, ty := float32(fx-float64(x0)), float32(fy-float64(y0))
	c := func(x, y int) float32 {
		x = max(0, min(g.w-1, x))
		y = max(0, min(g.h-1, y))
		return v[y*g.w+x]
	}
	return (c(x0, y0)*(1-tx)+c(x0+1, y0)*tx)*(1-ty) + (c(x0, y0+1)*(1-tx)+c(x0+1, y0+1)*tx)*ty
}

// intersect finds the first ground point along the ENU ray d from (e, n, u):
// march in steps of at most 10 m until the ray is below the ground, then
// bisect. A ray that never meets the ground ends where it passes the grid's
// lowest height.
func (g *ground) intersect(e, n, u float64, d vec3) (float64, float64, float64) {
	if d[2] >= -1e-3 {
		d[2] = -1e-3
	}
	above := func(t float64) bool {
		return u+t*d[2] > float64(g.at(g.elev, e+t*d[0], n+t*d[1]))
	}
	t0 := math.Max(0, (u-float64(g.maxEl)-1)/(-d[2]))
	tEnd := (u - float64(g.minElev) + 1) / (-d[2])
	step := 10 / math.Max(math.Hypot(d[0], d[1]), -d[2])
	t1 := t0
	for t1 < tEnd && above(t1) {
		t0 = t1
		t1 += step
	}
	t1 = math.Min(t1, tEnd)
	for i := 0; i < 24; i++ {
		m := (t0 + t1) / 2
		if above(m) {
			t0 = m
		} else {
			t1 = m
		}
	}
	t := (t0 + t1) / 2
	pe, pn := e+t*d[0], n+t*d[1]
	return pe, pn, float64(g.at(g.elev, pe, pn))
}

// ---- fast filters -------------------------------------------------------

// gaussAny is the exact gauss for small sigma and three box passes (widths
// after Kovesi) for large ones; both clamp at the edges.
func gaussAny(src *img, sigma float64) *img {
	if sigma < 3 {
		return gauss(src, sigma)
	}
	const n = 3
	wIdeal := math.Sqrt(12*sigma*sigma/n + 1)
	wl := int(math.Floor(wIdeal))
	if wl%2 == 0 {
		wl--
	}
	m := int(math.Round((12*sigma*sigma - float64(n*wl*wl) - float64(4*n*wl) - 3*n) / float64(-4*wl-4)))
	out := src
	for i := 0; i < n; i++ {
		wd := wl
		if i >= m {
			wd = wl + 2
		}
		out = boxBlur(out, (wd-1)/2)
	}
	return out
}

// boxBlur is a (2r+1)-wide box, horizontal then vertical, edges replicated.
func boxBlur(src *img, r int) *img {
	w, h := src.w, src.h
	tmp, dst := newImg(w, h), newImg(w, h)
	line := make([]float32, max(w, h))
	out := make([]float32, max(w, h))
	for y := 0; y < h; y++ {
		copy(line, src.px[y*w:(y+1)*w])
		box1(line[:w], out[:w], r)
		copy(tmp.px[y*w:(y+1)*w], out[:w])
	}
	for x := 0; x < w; x++ {
		for y := 0; y < h; y++ {
			line[y] = tmp.px[y*w+x]
		}
		box1(line[:h], out[:h], r)
		for y := 0; y < h; y++ {
			dst.px[y*w+x] = out[y]
		}
	}
	return dst
}

func box1(in, out []float32, r int) {
	n := len(in)
	pre := make([]float64, n+2*r+2)
	for i := 0; i < n+2*r+1; i++ {
		pre[i+1] = pre[i] + float64(in[max(0, min(n-1, i-r))])
	}
	inv := 1 / float64(2*r+1)
	for x := 0; x < n; x++ {
		out[x] = float32((pre[x+2*r+1] - pre[x]) * inv)
	}
}

// dilateFast/erodeFast are the separable square max/min of radius r, in O(1)
// per pixel (van Herk / Gil-Werman), edges replicated.
func dilateFast(src *img, r int) *img { return mmFast(src, r, false) }
func erodeFast(src *img, r int) *img  { return mmFast(src, r, true) }

func mmFast(src *img, r int, isMin bool) *img {
	w, h := src.w, src.h
	tmp, dst := newImg(w, h), newImg(w, h)
	line := make([]float32, max(w, h))
	out := make([]float32, max(w, h))
	for y := 0; y < h; y++ {
		copy(line, src.px[y*w:(y+1)*w])
		vhgw(line[:w], out[:w], r, isMin)
		copy(tmp.px[y*w:(y+1)*w], out[:w])
	}
	for x := 0; x < w; x++ {
		for y := 0; y < h; y++ {
			line[y] = tmp.px[y*w+x]
		}
		vhgw(line[:h], out[:h], r, isMin)
		for y := 0; y < h; y++ {
			dst.px[y*w+x] = out[y]
		}
	}
	return dst
}

func vhgw(in, out []float32, r int, isMin bool) {
	n := len(in)
	k := 2*r + 1
	m := n + 2*r
	m += (k - m%k) % k // whole blocks
	ext := make([]float32, m)
	for i := range ext {
		ext[i] = in[max(0, min(n-1, i-r))]
	}
	if isMin {
		for i := range ext {
			ext[i] = -ext[i]
		}
	}
	g, hh := make([]float32, m), make([]float32, m)
	for b := 0; b < m; b += k {
		g[b] = ext[b]
		for i := b + 1; i < b+k; i++ {
			g[i] = max(g[i-1], ext[i])
		}
		hh[b+k-1] = ext[b+k-1]
		for i := b + k - 2; i >= b; i-- {
			hh[i] = max(hh[i+1], ext[i])
		}
	}
	for x := 0; x < n; x++ {
		v := max(hh[x], g[x+2*r])
		if isMin {
			v = -v
		}
		out[x] = v
	}
}

// ---- frame io: 8-bit PNG, display-linear as v^1.8 (frame mode's
// arithmetic, with a lookup table for the decode) ----------------------------

var lin8 = func() (t [256]float32) {
	for i := range t {
		t[i] = float32(math.Pow(float64(i)/255, 1.8))
	}
	return
}()

func decodeRGBA(path string) (*image.RGBA, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	src, err := png.Decode(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	if rgba, ok := src.(*image.RGBA); ok && rgba.Bounds().Min == (image.Point{}) {
		return rgba, nil
	}
	b := src.Bounds()
	rgba := image.NewRGBA(image.Rect(0, 0, b.Dx(), b.Dy()))
	draw.Draw(rgba, rgba.Bounds(), src, b.Min, draw.Src)
	return rgba, nil
}

func linFromRGBA(src *image.RGBA) [3]*img {
	w, h := src.Bounds().Dx(), src.Bounds().Dy()
	var rgb [3]*img
	for c := range rgb {
		rgb[c] = newImg(w, h)
	}
	for y := 0; y < h; y++ {
		row := src.Pix[y*src.Stride:]
		for x := 0; x < w; x++ {
			for c := 0; c < 3; c++ {
				rgb[c].px[y*w+x] = lin8[row[4*x+c]]
			}
		}
	}
	return rgb
}

func encodeRGBA(rgb [3]*img) *image.RGBA {
	w, h := rgb[0].w, rgb[0].h
	out := image.NewRGBA(image.Rect(0, 0, w, h))
	for i := 0; i < w*h; i++ {
		for c := 0; c < 3; c++ {
			v := math.Min(math.Max(float64(rgb[c].px[i]), 0), 1)
			out.Pix[i*4+c] = uint8(math.Pow(v, 1/1.8) * 255)
		}
		out.Pix[i*4+3] = 255
	}
	return out
}

// writePNG writes under a temporary name and renames, so a file that exists
// is complete.
func writePNG(path string, im image.Image) {
	tmp := path + ".tmp"
	f, err := os.Create(tmp)
	if err != nil {
		fatal(err)
	}
	enc := png.Encoder{CompressionLevel: png.BestSpeed}
	if err := enc.Encode(f, im); err != nil {
		f.Close()
		fatal(fmt.Errorf("%s: %w", path, err))
	}
	if err := f.Close(); err != nil {
		fatal(err)
	}
	if err := os.Rename(tmp, path); err != nil {
		fatal(err)
	}
}

func writeGray(path string, a *img) {
	out := image.NewGray(image.Rect(0, 0, a.w, a.h))
	for i, v := range a.px {
		out.Pix[i] = uint8(clamp01(v) * 255)
	}
	writePNG(path, out)
}

func writeRaw(path string, a *img) {
	b := make([]byte, 4*len(a.px))
	for i, v := range a.px {
		u := math.Float32bits(v)
		b[4*i], b[4*i+1], b[4*i+2], b[4*i+3] = byte(u), byte(u>>8), byte(u>>16), byte(u>>24)
	}
	if err := os.WriteFile(path, b, 0644); err != nil {
		fatal(err)
	}
}

// must is for the cloud renderer's setup: any error there ends the run.
func must(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
