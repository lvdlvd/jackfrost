// Elevation for flight mode: every georeferenced GeoTIFF in a directory,
// heights as stored (MSL for the usual global DEMs), bilinear between posts.
//
// The tiepoint ties raster point (I,J) to model point (X,Y), and what raster
// point (0,0) is depends on GTRasterTypeGeoKey: in PixelIsArea (the GeoTIFF
// default) the top-left corner of pixel (0,0), in PixelIsPoint its centre.
// Copernicus GLO-30 is PixelIsPoint; taking its tiepoint as a corner samples
// it half a post (15 m) north-west of the point asked for.
package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/lvdlvd/geotiff"
)

type demTile struct {
	w, h                 int
	dLon, dLat           float64
	originLon, originLat float64 // top-left corner of pixel (0,0)
	elev                 []float32
}

// demSet is the tiles in glob order; where tiles overlap the first wins.
type demSet struct{ tiles []*demTile }

func loadDEM(dir string) (*demSet, error) {
	paths, err := filepath.Glob(filepath.Join(dir, "*.tif"))
	if err != nil {
		return nil, err
	}
	d := &demSet{}
	for _, p := range paths {
		t, err := loadDEMTile(p)
		if err != nil {
			return nil, err
		}
		if t != nil {
			d.tiles = append(d.tiles, t)
		}
	}
	if len(d.tiles) == 0 {
		return nil, fmt.Errorf("no georeferenced .tif in %s", dir)
	}
	return d, nil
}

// loadDEMTile reads one tile; nil, nil for a .tif that is not a
// georeferenced raster.
func loadDEMTile(path string) (*demTile, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	ifd, err := geotiff.OpenTIFF(f)
	if err != nil {
		return nil, nil
	}
	scale, ok := ifd.Value(geotiff.ModelPixelScale).([]float64)
	if !ok || len(scale) < 2 {
		return nil, nil
	}
	tp, ok := ifd.Value(geotiff.ModelTiepoint).([]float64)
	if !ok || len(tp) < 6 {
		return nil, nil
	}
	w, h := ifd.WxL()
	if w == 0 || h == 0 {
		return nil, nil
	}
	ci, cj := tp[0], tp[1] // the raster point tied, in corner coordinates
	if rasterType(ifd) == 2 {
		ci, cj = ci+0.5, cj+0.5
	}
	t := &demTile{w: w, h: h, dLon: scale[0], dLat: scale[1],
		originLon: tp[3] - ci*scale[0], originLat: tp[4] + cj*scale[1]}
	if t.elev, err = geotiff.ReadSamples(ifd); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return t, nil
}

// rasterType is GTRasterTypeGeoKey (1025): 1 PixelIsArea, 2 PixelIsPoint,
// 0 when absent (which means PixelIsArea).
func rasterType(ifd *geotiff.IFD) int {
	keys, ok := ifd.Value(geotiff.GeoKeyDirectory).([]uint16)
	if !ok || len(keys) < 4 {
		return 0
	}
	for i := 4; i+3 < len(keys); i += 4 {
		if keys[i] == 1025 && keys[i+1] == 0 {
			return int(keys[i+3])
		}
	}
	return 0
}

func (t *demTile) contains(lat, lon float64) bool {
	return lon >= t.originLon && lon <= t.originLon+float64(t.w)*t.dLon &&
		lat <= t.originLat && lat >= t.originLat-float64(t.h)*t.dLat
}

func (d *demSet) tile(lat, lon float64) *demTile {
	for _, t := range d.tiles {
		if t.contains(lat, lon) {
			return t
		}
	}
	return nil
}

// covers reports whether a tile holds the point; elevation is 0 both at sea
// level and where none does.
func (d *demSet) covers(lat, lon float64) bool { return d.tile(lat, lon) != nil }

func (d *demSet) elevation(lat, lon float64) float64 {
	if t := d.tile(lat, lon); t != nil {
		return t.sample(lat, lon)
	}
	return 0
}

// sample is bilinear between post centres, clamped at the tile's edge.
func (t *demTile) sample(lat, lon float64) float64 {
	fx := (lon-t.originLon)/t.dLon - 0.5
	fy := (t.originLat-lat)/t.dLat - 0.5
	x0, y0 := int(math.Floor(fx)), int(math.Floor(fy))
	tx, ty := fx-float64(x0), fy-float64(y0)
	at := func(x, y int) float64 {
		x = max(0, min(t.w-1, x))
		y = max(0, min(t.h-1, y))
		return float64(t.elev[y*t.w+x])
	}
	return (at(x0, y0)*(1-tx)+at(x0+1, y0)*tx)*(1-ty) + (at(x0, y0+1)*(1-tx)+at(x0+1, y0+1)*tx)*ty
}
