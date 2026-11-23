package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// pixel -> ray -> pixel closes everywhere on the sensor, the corners (beyond
// the polynomial's fold, on the extension) included; the extension is on.
func TestLensRoundTrip(t *testing.T) {
	l, err := LoadPinhole("testdata/lens.json")
	if err != nil {
		t.Fatal(err)
	}
	if l.FoldDeg() == 0 || l.ValidDeg() == 0 {
		t.Fatalf("test lens: fold %.2f deg, valid %.2f deg; want both", l.FoldDeg(), l.ValidDeg())
	}
	for _, p := range [][2]float64{{0.5, 0.5}, {1919.5, 0.5}, {960, 600}, {100, 1100}, {1919.5, 1199.5}, {1500, 300}} {
		d := l.Ray(p[0], p[1])
		qx, qy, ok := l.Pixel(d)
		if !ok || math.Hypot(qx-p[0], qy-p[1]) > 1e-3 {
			t.Errorf("pixel %v -> ray %v -> %.4f,%.4f", p, d, qx, qy)
		}
	}
}

// A level pose's nadir lands on the principal point whatever the heading;
// banked right it moves right, pitched up it moves down.
func TestNadirProjection(t *testing.T) {
	l, err := LoadPinhole("testdata/lens.json")
	if err != nil {
		t.Fatal(err)
	}
	const lat, lon, h, gnd = 46.7, 6.9, 1646.0, 900.0
	nadir := func(roll, pitch, yaw float64) (float64, float64) {
		v := newView(lat, lon, h, roll, pitch, yaw)
		px, py, ok := l.Pixel(v.toCam(v.enu(lat, lon, gnd)))
		if !ok {
			t.Fatalf("nadir not imaged at roll %v pitch %v yaw %v", roll, pitch, yaw)
		}
		return px, py
	}
	for _, yaw := range []float64{0, 137, 232.4} {
		px, py := nadir(0, 0, yaw)
		if math.Hypot(px-l.cx, py-l.cy) > 1e-6 {
			t.Errorf("yaw %v: nadir at %.4f,%.4f, principal point %.4f,%.4f", yaw, px, py, l.cx, l.cy)
		}
	}
	if px, py := nadir(15, 0, 137); px <= l.cx+100 || math.Abs(py-l.cy) > 10 {
		t.Errorf("banked right 15: nadir at %.1f,%.1f, want well right of %.1f", px, py, l.cx)
	}
	if px, py := nadir(0, 6, 137); py <= l.cy+50 || math.Abs(px-l.cx) > 10 {
		t.Errorf("pitched up 6: nadir at %.1f,%.1f, want below %.1f", px, py, l.cy)
	}
	// the snow's camera basis is the cloud's
	x1, y1, z1 := camBasisBody(3, -2, 77)
	v := newView(lat, lon, h, 3, -2, 77)
	if x1 != v.right || y1 != v.down || z1 != v.bore {
		t.Errorf("camBasisBody and newView disagree")
	}
}

// 2026-08-10 12:00 UTC at 46.7N 6.9E is some 22 minutes after local solar
// noon, declination about +15.6 deg: elevation near 90 - 46.7 + 15.6 =
// 58.9 deg less a little, azimuth a little west of south.
func TestSun(t *testing.T) {
	el, az := sunPosition(1786363200, 46.7, 6.9) // 2026-08-10T12:00:00Z
	if el < 57 || el > 59.5 || az < 185 || az > 200 {
		t.Errorf("sun elevation %.2f azimuth %.2f", el, az)
	}
}

// The plowed share of each class's length is its probability, decisions
// are the same on every load, and -plow all=1 plows everything.
func TestPlowShare(t *testing.T) {
	g := &ground{lat0: 47, lon0: 8}
	g.mLat = 6371000 * math.Pi / 180
	g.mLon = g.mLat * math.Cos(g.lat0*math.Pi/180)
	g.e0, g.n0 = -30000, -30000
	g.w, g.h = 6001, 6001
	// 300 m road pieces on a 600 m lattice over 60 km, alternately
	// secondary and residential
	type feat struct {
		Type       string            `json:"type"`
		Properties map[string]string `json:"properties"`
		Geometry   struct {
			Type        string       `json:"type"`
			Coordinates [][2]float64 `json:"coordinates"`
		} `json:"geometry"`
	}
	var fs []feat
	k := 0
	for e := -29000.0; e < 29000; e += 600 {
		for n := -29000.0; n < 29000; n += 600 {
			var f feat
			f.Type = "Feature"
			f.Properties = map[string]string{"highway": []string{"secondary", "residential"}[k%2]}
			f.Geometry.Type = "LineString"
			for _, d := range []float64{0, 150, 300} {
				f.Geometry.Coordinates = append(f.Geometry.Coordinates, [2]float64{g.lon0 + (e+d)/g.mLon, g.lat0 + n/g.mLat})
			}
			fs = append(fs, f)
			k++
		}
	}
	b, err := json.Marshal(map[string]any{"type": "FeatureCollection", "features": fs})
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "roads.geojson")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	share := func(rs *roadSet, class string) float64 {
		// the plowed pieces are those whose middle vertex is in a corridor
		all, plowed := 0, 0
		for _, f := range fs {
			if f.Properties["highway"] != class {
				continue
			}
			all++
			m := f.Geometry.Coordinates[1]
			e, n := g.toEN(m[1], m[0])
			if rs.near(float32(e), float32(n)) {
				plowed++
			}
		}
		return float64(plowed) / float64(all)
	}
	rs := loadRoads(path, g)
	sec, res := share(rs, "secondary"), share(rs, "residential")
	if math.Abs(sec-plowProb["secondary"]) > 0.08 || math.Abs(res-plowProb["residential"]) > 0.08 {
		t.Errorf("plowed share secondary %.3f (want %.2f), residential %.3f (want %.2f)", sec, plowProb["secondary"], res, plowProb["residential"])
	}
	if again := loadRoads(path, g); len(again.segs) != len(rs.segs) {
		t.Errorf("a second load plowed %d segments, the first %d", len(again.segs), len(rs.segs))
	}
	saved := map[string]float64{}
	for c, p := range plowProb {
		saved[c] = p
	}
	defer func() { plowProb = saved }()
	parsePlow("all=1")
	if all := loadRoads(path, g); share(all, "secondary") != 1 || share(all, "residential") != 1 {
		t.Errorf("all=1 left roads unplowed")
	}
}
