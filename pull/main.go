// pull fetches a paired summer/winter Sentinel-2 L2A capture of one MGRS
// tile from the public earth-search STAC (sentinel-cogs bucket, no auth):
// the snowiest acceptable winter scene over the bbox, and the clearest
// summer scene of the SAME tile, so the pair is pixel-aligned on the fixed
// MGRS grid. Bands: RGB + B11 (SWIR, for NDSI) + SCL.
//
// Datetimes must be full RFC3339 — the API 400s on bare dates.
package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"net/http"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

const api = "https://earth-search.aws.element84.com/v1/search"

var bands = []string{"red", "green", "blue", "swir16", "scl"}

var (
	bboxStr  = flag.String("bbox", "", "region of interest: w,s,e,n (degrees)")
	outDir   = flag.String("out", "raw", "output directory")
	winters  = flag.Int("winters", 4, "how many past winters to search")
	maxCloud = flag.Float64("maxcloud", 30, "max cloud cover %% for the winter scene")
)

type feature struct {
	ID         string         `json:"id"`
	Properties map[string]any `json:"properties"`
	Assets     map[string]struct {
		Href string `json:"href"`
	} `json:"assets"`
}

func (f feature) num(key string) float64 {
	if v, ok := f.Properties[key].(float64); ok {
		return v
	}
	return -1
}

func (f feature) tile() string {
	parts := strings.Split(f.ID, "_")
	if len(parts) < 2 {
		return ""
	}
	return parts[1]
}

func search(bbox []float64, dateRange string) ([]feature, error) {
	q := map[string]any{
		"collections": []string{"sentinel-2-l2a"},
		"bbox":        bbox,
		"datetime":    dateRange,
		"limit":       100,
	}
	body, _ := json.Marshal(q)
	resp, err := http.Post(api, "application/json", bytes.NewReader(body))
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()
	if resp.StatusCode != 200 {
		msg, _ := io.ReadAll(resp.Body)
		return nil, fmt.Errorf("search %s: %s: %s", dateRange, resp.Status, msg)
	}
	var out struct {
		Features []feature `json:"features"`
	}
	return out.Features, json.NewDecoder(resp.Body).Decode(&out)
}

func download(f feature, tag string) error {
	for _, band := range bands {
		a, ok := f.Assets[band]
		if !ok {
			return fmt.Errorf("%s: no asset %q", f.ID, band)
		}
		dst := fmt.Sprintf("%s/%s_%s_%s.tif", *outDir, tag, f.ID, band)
		if _, err := os.Stat(dst); err == nil {
			fmt.Println("have", dst)
			continue
		}
		fmt.Println("pull", dst)
		resp, err := http.Get(a.Href)
		if err != nil {
			return err
		}
		if resp.StatusCode != 200 {
			resp.Body.Close()
			return fmt.Errorf("%s: %s", a.Href, resp.Status)
		}
		tmp := dst + ".part"
		w, err := os.Create(tmp)
		if err != nil {
			resp.Body.Close()
			return err
		}
		n, err := io.Copy(w, resp.Body)
		resp.Body.Close()
		w.Close()
		if err != nil {
			return err
		}
		if err := os.Rename(tmp, dst); err != nil {
			return err
		}
		fmt.Printf("  %d MB\n", n>>20)
	}
	return nil
}

func main() {
	flag.Parse()
	parts := strings.Split(*bboxStr, ",")
	if len(parts) != 4 {
		fmt.Fprintln(os.Stderr, "pull: -bbox w,s,e,n is required")
		os.Exit(1)
	}
	bbox := make([]float64, 4)
	for i, p := range parts {
		v, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil {
			fmt.Fprintln(os.Stderr, "pull: bad bbox:", err)
			os.Exit(1)
		}
		bbox[i] = v
	}
	if err := os.MkdirAll(*outDir, 0755); err != nil {
		fmt.Fprintln(os.Stderr, "pull:", err)
		os.Exit(1)
	}

	// winter: snowiest acceptable scene over the bbox
	var wc []feature
	year := time.Now().Year()
	for y := 0; y < *winters; y++ {
		rng := fmt.Sprintf("%d-11-15T00:00:00Z/%d-03-15T23:59:59Z", year-1-y, year-y)
		fs, err := search(bbox, rng)
		if err != nil {
			fmt.Fprintln(os.Stderr, "pull:", err)
			os.Exit(1)
		}
		wc = append(wc, fs...)
	}
	wc = filter(wc, func(f feature) bool {
		// partial-coverage tiles (swath edges) make chip selection
		// impossible: the summer and winter slivers may not even overlap
		return f.num("eo:cloud_cover") < *maxCloud &&
			f.num("s2:nodata_pixel_percentage") < 20
	})
	if len(wc) == 0 {
		fmt.Fprintln(os.Stderr, "pull: no acceptable winter scene")
		os.Exit(1)
	}
	sort.Slice(wc, func(i, j int) bool {
		return wc[i].num("s2:snow_ice_percentage") > wc[j].num("s2:snow_ice_percentage")
	})
	for _, f := range wc[:min(5, len(wc))] {
		fmt.Printf("winter candidate %s snow%% %.1f cloud%% %.1f\n",
			f.ID, f.num("s2:snow_ice_percentage"), f.num("eo:cloud_cover"))
	}
	winter := wc[0]

	// summer: clearest scene of the SAME tile
	var sc []feature
	for y := 0; y < 2; y++ {
		rng := fmt.Sprintf("%d-06-01T00:00:00Z/%d-09-10T23:59:59Z", year-y, year-y)
		fs, err := search(bbox, rng)
		if err != nil {
			fmt.Fprintln(os.Stderr, "pull:", err)
			os.Exit(1)
		}
		sc = append(sc, fs...)
	}
	sc = filter(sc, func(f feature) bool {
		return f.tile() == winter.tile() && f.num("s2:nodata_pixel_percentage") < 20
	})
	if len(sc) == 0 {
		fmt.Fprintf(os.Stderr, "pull: no summer scene for tile %s\n", winter.tile())
		os.Exit(1)
	}
	sort.Slice(sc, func(i, j int) bool {
		// clearest sky first, residual (glacier) snow as tie-breaker
		a := sc[i].num("eo:cloud_cover") + sc[i].num("s2:snow_ice_percentage")
		b := sc[j].num("eo:cloud_cover") + sc[j].num("s2:snow_ice_percentage")
		return a < b
	})
	summer := sc[0]

	fmt.Printf("pair: tile %s  summer %s (cloud%% %.1f)  winter %s (snow%% %.1f cloud%% %.1f)\n",
		winter.tile(), summer.ID, summer.num("eo:cloud_cover"),
		winter.ID, winter.num("s2:snow_ice_percentage"), winter.num("eo:cloud_cover"))

	if err := download(summer, "summer"); err != nil {
		fmt.Fprintln(os.Stderr, "pull:", err)
		os.Exit(1)
	}
	if err := download(winter, "winter"); err != nil {
		fmt.Fprintln(os.Stderr, "pull:", err)
		os.Exit(1)
	}
	fmt.Println("JACKFROST-PULL-DONE", summer.ID, winter.ID)
}

func filter(fs []feature, keep func(feature) bool) []feature {
	var out []feature
	for _, f := range fs {
		if keep(f) {
			out = append(out, f)
		}
	}
	return out
}
