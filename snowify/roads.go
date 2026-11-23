// Plowed roads (flight snow, -roads). Snowify keeps dark thin features
// through its top-hat, which at Sentinel-2's 10 m covers the roads (sub-pixel
// there). At 0.4-1 m/px the roads of a summer orthophoto are light grey
// asphalt, brighter than the fields, so the top-hat does not keep them, the
// open-land gates pass them, and they vanish under the snow. In mid-winter
// the important ones are plowed and the minor ones often are not.
//
// Which roads are plowed. Each road class has a probability of being plowed
// (plowProb, -plow to change it): near certain for motorways, rarely for
// service roads. Plowing goes by route, not by the hundred-metre pieces a
// map cuts a road into, so the decision is drawn from a smooth random field
// per class, value noise with a correlation length of -plow-km, sampled at
// the middle of each road piece: neighbouring pieces of a class agree, a
// minor road is plowed or not for kilometres at a time, and the share of
// plowed road of each class is its probability. The field's threshold for a
// probability is its empirical quantile over the ground grid. The decision
// is fixed for the flight (-plow-seed), so a road is plowed in every frame
// or in none.
//
// Where a plowed road is clear: a pixel stays clear (keeps its summer
// appearance, the bare asphalt) when its ground point lies within the road's
// half-width (by class, halfWidth below, plus -road-slack) of the centreline
// AND the pixel is grey (chroma (max-min)/max of the display RGB below
// -road-chroma). The corridor only says where to look; the colour test picks
// the asphalt actually in the picture, so the cleared strip sits where the
// imagery's road is, not where the map draws it. Verges, fields and trees
// inside the corridor are green or brown and get snow as anywhere else. The
// clear mask is softened by the same PSF as the cover.
//
// The input is GeoJSON LineStrings or MultiLineStrings with the
// OpenStreetMap "highway" class as a property, e.g. the paved roads of an
// extract.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"sort"
	"strconv"
	"strings"
)

var (
	roadsPath  = flag.String("roads", "", "flight snow: roads GeoJSON (LineStrings with a \"highway\" property); plowed ones keep their asphalt clear")
	roadSlack  = flag.Float64("road-slack", 1.5, "flight snow: corridor margin beyond the class half-width, m (map position error)")
	roadChroma = flag.Float64("road-chroma", 0.12, "flight snow: max chroma (max-min)/max of a cleared asphalt pixel")
	plowSpec   = flag.String("plow", "", "flight snow: class=probability,... overriding the plowed share of a road class (all=1 plows every road)")
	plowKm     = flag.Float64("plow-km", 3, "flight snow: correlation length of the plowing decision along the road network, km")
	plowSeed   = flag.Uint64("plow-seed", 1, "flight snow: seed of the plowing decision")
)

const bucket = 50.0 // metres

// halfWidth is the paved half-width by highway class, metres: typical
// European carriageways (motorway one direction with shoulder, two-lane main
// roads, village streets, single-lane service roads). Classes not listed:
// 3 m.
var halfWidth = map[string]float32{
	"motorway": 6, "trunk": 4.5, "primary": 4.5, "secondary": 3.5, "tertiary": 3.5,
	"motorway_link": 3.5, "trunk_link": 3.5, "primary_link": 3.5, "secondary_link": 3, "tertiary_link": 3,
	"residential": 3, "unclassified": 3, "living_street": 3, "service": 2.5,
}

// plowProb is the share of each class that is plowed on a winter day: the
// network a road service clears first is the through roads, then the bus
// and access roads, and only some of the streets and service roads. These
// are estimates, not measurements. A link takes its road's value; a class
// not listed takes "other".
var plowProb = map[string]float64{
	"motorway": 1, "trunk": 0.97, "primary": 0.93, "secondary": 0.8, "tertiary": 0.6,
	"unclassified": 0.4, "residential": 0.35, "living_street": 0.2, "service": 0.1,
	"other": 0.2,
}

func plowProbability(class string) float64 {
	if p, ok := plowProb[class]; ok {
		return p
	}
	if p, ok := plowProb[strings.TrimSuffix(class, "_link")]; ok {
		return p
	}
	return plowProb["other"]
}

// parsePlow applies -plow: "all=1" or "class=p,class=p".
func parsePlow(spec string) {
	if spec == "" {
		return
	}
	for _, kv := range strings.Split(spec, ",") {
		k, v, ok := strings.Cut(strings.TrimSpace(kv), "=")
		p, err := strconv.ParseFloat(v, 64)
		if !ok || err != nil || p < 0 || p > 1 {
			fatal(fmt.Errorf("-plow %q: want class=probability,... with probabilities in [0,1]", spec))
		}
		if k == "all" {
			for c := range plowProb {
				plowProb[c] = p
			}
			continue
		}
		plowProb[k] = p
	}
}

type seg struct{ ae, an, be, bn, half float32 }

type roadSet struct {
	g       *ground
	bw, bh  int
	buckets [][]int32
	segs    []seg
}

// plowField is the smooth decision field of one class: value noise, three
// octaves, in ground metres.
func plowField(e, n float64, class string) float64 {
	seed := int64(*plowSeed)
	for _, c := range class {
		seed = seed*131 + int64(c)
	}
	return cloudNoise(e, n, *plowKm*1000, seed)
}

// plowThreshold is the field value below which a fraction p of the ground
// lies: the empirical quantile over random points of the grid. The field's
// distribution does not depend on the seed, so one sample serves every
// class.
func plowThreshold(p float64, qs []float64) float64 {
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	return qs[min(int(p*float64(len(qs))), len(qs)-1)]
}

func loadRoads(path string, g *ground) *roadSet {
	parsePlow(*plowSpec)
	b, err := os.ReadFile(path)
	if err != nil {
		fatal(err)
	}
	var fc struct {
		Features []struct {
			Properties map[string]any `json:"properties"`
			Geometry   struct {
				Type        string          `json:"type"`
				Coordinates json.RawMessage `json:"coordinates"`
			} `json:"geometry"`
		} `json:"features"`
	}
	if err := json.Unmarshal(b, &fc); err != nil {
		fatal(err)
	}
	// the field's quantiles, from 20000 points of the grid
	rng := rand.New(rand.NewPCG(*plowSeed, 7))
	qs := make([]float64, 20000)
	for i := range qs {
		qs[i] = plowField(g.e0+rng.Float64()*float64(g.w)*cell, g.n0+rng.Float64()*float64(g.h)*cell, "quantiles")
	}
	sort.Float64s(qs)

	rs := &roadSet{g: g}
	rs.bw = int(float64(g.w)*cell/bucket) + 1
	rs.bh = int(float64(g.h)*cell/bucket) + 1
	rs.buckets = make([][]int32, rs.bw*rs.bh)
	type tally struct{ all, plowed float64 }
	byClass := map[string]*tally{}
	for _, f := range fc.Features {
		var lines [][][2]float64
		switch f.Geometry.Type {
		case "LineString":
			var l [][2]float64
			if json.Unmarshal(f.Geometry.Coordinates, &l) != nil {
				continue
			}
			lines = append(lines, l)
		case "MultiLineString":
			if json.Unmarshal(f.Geometry.Coordinates, &lines) != nil {
				continue
			}
		}
		class, _ := f.Properties["highway"].(string)
		half, ok := halfWidth[class]
		if !ok {
			half = 3
		}
		half += float32(*roadSlack)
		for _, l := range lines {
			if len(l) < 2 {
				continue
			}
			// the decision at the piece's middle vertex, its length for the tally
			length := 0.0
			for i := 0; i+1 < len(l); i++ {
				ae, an := g.toEN(l[i][1], l[i][0])
				be, bn := g.toEN(l[i+1][1], l[i+1][0])
				length += math.Hypot(be-ae, bn-an)
			}
			me, mn := g.toEN(l[len(l)/2][1], l[len(l)/2][0])
			plowed := plowField(me, mn, class) < plowThreshold(plowProbability(class), qs)
			t := byClass[class]
			if t == nil {
				t = &tally{}
				byClass[class] = t
			}
			t.all += length
			if !plowed {
				continue
			}
			t.plowed += length
			for i := 0; i+1 < len(l); i++ {
				ae, an := g.toEN(l[i][1], l[i][0])
				be, bn := g.toEN(l[i+1][1], l[i+1][0])
				rs.add(seg{float32(ae), float32(an), float32(be), float32(bn), half})
			}
		}
	}
	var classes []string
	for c := range byClass {
		classes = append(classes, c)
	}
	sort.Strings(classes)
	var parts []string
	for _, c := range classes {
		t := byClass[c]
		parts = append(parts, fmt.Sprintf("%s %.0f of %.0f km", c, t.plowed/1000, t.all/1000))
	}
	fmt.Printf("roads: %d plowed segments; plowed by class: %s\n", len(rs.segs), strings.Join(parts, ", "))
	return rs
}

func (rs *roadSet) add(s seg) {
	h := s.half
	x0 := int((min(s.ae, s.be) - h - float32(rs.g.e0)) / bucket)
	x1 := int((max(s.ae, s.be) + h - float32(rs.g.e0)) / bucket)
	y0 := int((min(s.an, s.bn) - h - float32(rs.g.n0)) / bucket)
	y1 := int((max(s.an, s.bn) + h - float32(rs.g.n0)) / bucket)
	if x1 < 0 || y1 < 0 || x0 >= rs.bw || y0 >= rs.bh {
		return
	}
	id := int32(len(rs.segs))
	rs.segs = append(rs.segs, s)
	for y := max(0, y0); y <= min(rs.bh-1, y1); y++ {
		for x := max(0, x0); x <= min(rs.bw-1, x1); x++ {
			rs.buckets[y*rs.bw+x] = append(rs.buckets[y*rs.bw+x], id)
		}
	}
}

// near reports whether ground point (e, n) is within a plowed corridor.
func (rs *roadSet) near(e, n float32) bool {
	x := int((e - float32(rs.g.e0)) / bucket)
	y := int((n - float32(rs.g.n0)) / bucket)
	if x < 0 || y < 0 || x >= rs.bw || y >= rs.bh {
		return false
	}
	for _, id := range rs.buckets[y*rs.bw+x] {
		s := rs.segs[id]
		h2 := s.half * s.half
		dx, dy := s.be-s.ae, s.bn-s.an
		l2 := dx*dx + dy*dy
		t := float32(0)
		if l2 > 0 {
			t = max(0, min(1, ((e-s.ae)*dx+(n-s.an)*dy)/l2))
		}
		px, py := s.ae+t*dx-e, s.an+t*dy-n
		if px*px+py*py <= h2 {
			return true
		}
	}
	return false
}

// clear zeroes the snow probability on plowed asphalt and returns the
// cleared fraction of the frame.
func (rs *roadSet) clear(pcov *img, lin [3]*img, pe, pn *img, psf float64) float64 {
	m := newImg(pcov.w, pcov.h)
	for i := range m.px {
		if !rs.near(pe.px[i], pn.px[i]) {
			continue
		}
		var mx, mn float64 = 0, 1
		for c := 0; c < 3; c++ {
			v := math.Pow(float64(lin[c].px[i]), 1/1.8)
			mx, mn = math.Max(mx, v), math.Min(mn, v)
		}
		if mx > 0 && (mx-mn)/mx < *roadChroma {
			m.px[i] = 1
		}
	}
	m = gauss(m, psf)
	var sum float64
	for i := range pcov.px {
		pcov.px[i] *= 1 - m.px[i]
		sum += float64(m.px[i])
	}
	return sum / float64(len(m.px))
}
