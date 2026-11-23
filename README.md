# jackfrost

Simulated winter for aerial and satellite imagery: snow cover calibrated
against real Sentinel-2 winter captures, and a layer of scattered cumulus
with its shadows, both drawn onto the frames of a simulated flight so that
the weather stays where it is on the ground from one frame to the next.

![summer | s=0.5 | s=1.0 | real cover | real winter](doc/chip0_panel.jpg)

*A 15 km validation chip (MGRS 31TDL, 10 m/px). Left to right: the real
summer capture; snowified at s=0.5 and s=1.0 (synthetic patchy cover);
snowified with the cover forced to the real winter snow mask; the real
winter capture of the same tile. Only the rightmost panel is real winter.*

Snow is the extreme appearance shift for visual navigation: it erases the
field mosaic and most photometric texture, while dark linear features —
forest boundaries, hedgerows, streams, plowed roads — stay visible or grow
MORE salient. Clouds are the other one: they hide parts of the ground,
cast shadows that move with the wind, and move against the ground
themselves. jackfrost produces both, so that a vision system can be tested
on the same flight in clear weather, under snow, under cloud, and under
both.

Contents:

- [The tools](#the-tools)
- [Building](#building)
- [How the snow works](#how-the-snow-works)
- [How the clouds work](#how-the-clouds-work)
- [Flight mode: inputs](#flight-mode-inputs)
- [How-to](#how-to)
- [Outputs](#outputs)
- [Option reference](#option-reference)
- [What is not modelled](#what-is-not-modelled)

## The tools

`snowify` has three modes, chosen by the flags given:

| mode | chosen by | what it does |
|---|---|---|
| chip mode | no `-frames`, no `-poses` | snowifies Sentinel-2 summer chips and grades them against the real winter capture of the same tile; with `-fit N` it refits the snow's appearance constants |
| frame mode | `-frames DIR` | snowifies every PNG in a directory, each frame on its own |
| flight mode | `-poses FILE` | snow, clouds or both on the rendered frames of a flight, tied to the ground through the camera and the terrain |

`pull` fetches a paired summer/winter Sentinel-2 capture of any region, for
calibrating the snow on a new landscape.

## Building

jackfrost is Go (1.25 or later) with one dependency,
[github.com/lvdlvd/geotiff](https://github.com/lvdlvd/geotiff), which
`go.mod` expects checked out beside it:

```
git clone https://github.com/lvdlvd/geotiff
git clone https://github.com/lvdlvd/jackfrost
cd jackfrost
go install ./snowify ./pull      # into $(go env GOPATH)/bin
go test ./snowify
```

## How the snow works

`snowify` whitens open land toward a calibrated snow luminance modulated
by the terrain's own relief shading, lets thin dark features survive via a
morphological top-hat, keeps forests and water dark with steep
texture/luminance gates, and sweeps partial cover (the snow degree s)
through a smooth noise field so patches keep sharp boundaries. A sub-pixel
blur on the coverage map stands in for the sensor point-spread at snow
edges.

### Calibration

Paired same-tile Sentinel-2 L2A captures (MGRS 31TDL/31TCL; summer
S2B_20250807 vs winter S2A_20240120, 36% snow): snow-covered land is ~9x
brighter than clear land and carries ~6x the gradient energy — snow is
bright relief crossed by a dark lattice, not a whiteout. The appearance
constants (`defaultParams`) are fitted, not hand-set: `-fit N` random-
searches them to match an 8-statistic vector (luminance and gradient
percentiles in and out of the real snow mask, white fraction) against the
real winter chips, with the snow cover forced to the real winter NDSI mask
so appearance is fitted decoupled from cover placement.

### Landscape priors

Snow appearance is landscape-dependent: open plains bury the field
structure under smooth uniform white, bocage keeps a dark hedgerow
lattice, the high foreland pushes brightness further and turns the
lattice nearly off. Parameters fitted on one landscape do not transfer
(objective inflates ~30x), but per-landscape fits recover near-home
quality. `params/` holds the fitted sets for `-params`:

| file | what |
|---|---|
| `bocage.json` | hedgerow country; the built-in default |
| `lowland.json` | open plains |
| `alpine.json` | the high foreland |
| `winter_ortho.json` | open plains, with the gates set for summer orthophotos: the set for rendered aerial frames |

`frames_params.json` carries frame mode's overrides for rendered aerial
frames.

### From satellite chips to frames

Frame mode normalizes each frame's radiometry into the calibration space
for the gates (median luminance to the summer clear-land level) and blends
in the frame's own display space; `-gsd` rescales the physical-length
constants when the frame resolution differs from Sentinel-2's 10 m/px.

On the frames of a flight that is not enough: frames come 15 a second, so
a noise field drawn in each picture makes the snow patches jump from frame
to frame, and from a low aircraft the resolution runs from 0.1 to 1 m/px
within one flight. Flight mode keeps the snow model and changes where its
inputs come from. Every 8th pixel's ray is cast through the lens onto the
DEM, and from those ground points:

- each frame's ground resolution sets the scale of the physical-length
  constants (texture, top-hat radius, relief shading);
- the partial-cover noise is a 10 m field on the ground, sampled at each
  pixel's ground point, so a snow patch lies at the same place in every
  frame;
- the snow degree rises with the ground height under the pixel, from
  `-s-low` below `-h-low` to `-s-high` above `-h-high`, and a second
  parameter set (`-params-high`) can be blended in by the frame's ground
  height;
- open water (`-water`) stays open: summer orthophotos often show lakes
  and rivers as bright as the fields, which the luminance gate alone would
  snow over;
- the asphalt of plowed roads (`-roads`) stays clear. A road is plowed by
  chance, by class: motorways always, main roads mostly, side streets and
  service roads seldom (table below; `-plow` changes it). The decision is
  drawn from a smooth random field along the network, so a road is plowed
  or not for kilometres at a time (`-plow-km`), and it is fixed for the
  flight. On a plowed road a pixel stays clear where its ground point is on
  the carriageway and the pixel is grey, so the cleared strip follows the
  asphalt in the picture, not the map's centreline.

| class | plowed | | class | plowed |
|---|---|---|---|---|
| motorway | 1.00 | | unclassified | 0.40 |
| trunk | 0.97 | | residential | 0.35 |
| primary | 0.93 | | living_street | 0.20 |
| secondary | 0.80 | | service | 0.10 |
| tertiary | 0.60 | | any other class | 0.20 |

A `_link` road takes its road's share. These are estimates of how a road
service works through a snowfall, not measurements.

The theory of operation is at the top of `snowify/flight.go`.

## How the clouds work

The clouds are a 3-D layer of scattered cumulus, rendered through the lens
of the frames with the full attitude of each pose:

- **the layer** lies between `-cloud-base` and `-cloud-top` metres above
  the local ground (the terrain smoothed over `-cloud-smooth-km`), so it
  follows valleys and plateaus but not every hill;
- **the cells** are clusters of smooth puffs with log-normal sizes around
  `-cloud-cell-km`, ragged at the edge by value noise, stamped at random
  until `-cloud-cover` of the sky is opaque seen from straight above. Flat
  bases, domed and cauliflowered tops, the biggest cells reaching the top
  of the layer; opaque in the core (`-cloud-sigma`), thin at the edge;
- **the wind** carries the whole field: `-wind N,E`, or the mean of a wind
  log inside the layer (`-wind-log`);
- **the light**: each pixel's ray is marched through the layer to the
  ground, front to back; sunlit domes are white, lee sides and deeper parts
  grey, the sun's position computed from the time of each pose;
- **the shadows**: from each ground point a ray toward the sun gives the
  shadow, slightly blue from the sky light (`-shadow-amb`, `-shadows`);
- **the aircraft** may be below the layer (looking down, the camera never
  sees the clouds, only their shadows), inside it (passing through a cell
  turns the frame white), or above it (the ground shows in the gaps).

The cloud is drawn over the input frame and nothing else of it changes.
With `-weather snow+cloud` it is drawn over the snow in the same pass; the
result is the same as a cloud run over the output of a snow run.

The theory of operation is at the top of `snowify/clouds.go`.

## Flight mode: inputs

| input | flag | format |
|---|---|---|
| poses | `-poses` | JSON, below: one pose per frame, in frame order |
| frames | `-frames` | 8-bit RGB PNGs, `-pattern` with the pose index, e.g. `frame_00042.png` |
| lens | `-calib` | the calibration the frames were rendered through, below |
| terrain | `-dem` | a directory of GeoTIFF elevation tiles in the poses' height datum (MSL for the usual global DEMs, e.g. Copernicus GLO-30); PixelIsArea and PixelIsPoint rasters both work |
| roads | `-roads` | GeoJSON LineStrings or MultiLineStrings with an OpenStreetMap `highway` class property; snow only, optional |
| water | `-water` | GeoJSON Polygons (outer ring and holes); snow only, optional |
| wind | `-wind` or `-wind-log` | clouds only, below |

Poses:

```json
{"poses": [
  {"t_s": 0.0, "time_utc": "2026-01-20T10:30:00Z",
   "lat_deg": 47.1000, "lon_deg": 8.2000, "alt_m_msl": 1250.0,
   "heading_deg": 232.4,
   "body_roll_deg": 0.0, "body_pitch_deg": 2.1, "body_yaw_deg": 232.4},
  ...
]}
```

`lat_deg`, `lon_deg`, `alt_m_msl` and `heading_deg` are required. The
body attitude (aerospace convention: yaw, pitch nose up, roll right wing
down) is optional; without it the camera is level at the heading. The
camera looks straight down through the floor: image right is the right
wing, image top the nose. For clouds, `t_s` is the flight time in seconds
(it moves the clouds with the wind) and `time_utc` (RFC 3339; only the
first pose's is read) places the sun; `-t0` gives the time of `t_s = 0`
when the poses carry none.

Lens (Brown-Conrady, OpenCV's convention, pixel centres at integers):

```json
{"image_size": [1920, 1200], "fx": 900, "fy": 900, "cx": 959.5, "cy": 599.5,
 "k1": -0.28, "k2": 0.07, "k3": -0.008, "p1": 0.0004, "p2": -0.0003,
 "valid_theta_deg": 50}
```

`valid_theta_deg`, optional, is the field angle of the outermost
calibration corner; beyond it the model is continued smoothly instead of
evaluated, so a wide lens's corners have rays even where its polynomial
folds back (`snowify/lens.go`). The frames must be `image_size`.

Wind log, for `-wind-log`: JSON lines, the wind the aircraft flew in, with
its height above the field; the clouds move with the mean of the samples
inside the layer:

```
{"h_field_m": 512.3, "wind_n_ms": 3.9, "wind_e_ms": 4.4}
```

## How-to

### Check and refit the snow calibration

```
snowify -out check                     # the 31TDL chips: panels and grades.json
snowify -out refit -fit 600            # refit the appearance constants: refit/fitted.json
```

### Calibrate on a new landscape

```
pull -bbox "7.4,46.1,8.3,46.6" -out raw_alpine    # ends with the pair's scene ids
snowify -raw raw_alpine -summer <summer id> -winter <winter id> -autochips -fit 600 -out alpine
```

`pull` searches the last `-winters` winters for the snowiest acceptable
scene over the box and the clearest summer scene of the same MGRS tile,
from the public earth-search STAC (sentinel-cogs, no account).
`-autochips` picks validation chips from the winter scene's snow and cloud
masks.

### Snow on a directory of frames, each on its own

```
snowify -frames renders/ -out renders_snow/ -s 0.7 -gsd 2 -params frames_params.json
```

### Snow on a flight

```
snowify -poses poses.json -calib lens.json -dem dem/ -frames clear/ -out snow/ \
        -weather snow -params params/winter_ortho.json -a-fixed -1 \
        -roads roads.geojson -water water.geojson -h-low 600 -h-high 800
```

`-a-fixed -1` takes one brightness scale for the whole flight from every
150th frame, so the snow line does not move when a town or a forest fills
the frame. Without `-roads` every road disappears under the snow; without
`-water` bright lakes may be snowed over.

Fewer plowed roads, a quieter network, or every road plowed:

```
snowify ... -plow secondary=0.5,tertiary=0.3,residential=0.1
snowify ... -plow all=1
```

A different but equally likely snowfall: `-plow-seed 2`. The snow's cover
noise is fixed; for less snow lower `-s-low` and `-s-high`.

### Clouds on a flight

```
snowify -poses poses.json -calib lens.json -dem dem/ -frames clear/ -out cloud/ \
        -weather cloud -cloud-cover 0.375 -cloud-base 500 -cloud-top 800 \
        -wind 3.9,4.5 -cloud-stats cloud/stats.csv -cloud-masks cloud/masks
```

Broken cloud low over the ground, small cells, no shadows:

```
snowify ... -weather cloud -cloud-cover 0.6 -cloud-base 250 -cloud-top 450 -cloud-cell-km 0.3 -shadows=false
```

Another sky with the same statistics: `-cloud-seed 2`.

### Snow and clouds together

```
snowify -poses poses.json -calib lens.json -dem dem/ -frames clear/ -out winter_cloud/ \
        -weather snow+cloud -params params/winter_ortho.json -a-fixed -1 \
        -roads roads.geojson -water water.geojson \
        -cloud-cover 0.375 -cloud-base 500 -cloud-top 800 -wind 3.9,4.5
```

### Part of a flight, restarts, and parallel runs

```
snowify ... -from 3000 -to 3999          # poses 3000..3999
snowify ... -from 0 -step 15             # one frame a second at 15 Hz
snowify ... -list 900,5400,11500         # just these
```

An output frame is written under a temporary name and renamed, so a frame
that exists is complete: a stopped run continues where it stopped when run
again, and several machines can share one flight by `-from`/`-to`. A frame
is never rewritten; delete the output to redo it.

### Look before rendering a cloudy flight

```
snowify -poses poses.json -calib lens.json -dem dem/ -weather cloud -wind 3.9,4.5 -profile
```

prints the layer under the flight every 30 s, when the aircraft is below,
in or above it, how long it spends inside cloud, and the cover realised
over the flight's footprint. `-probe 900,5400` prints what the ray through
the principal point of those frames meets; `-motion 900,901` checks that
the clouds move through the frames as the geometry and the wind say
(measured by block matching against the prediction; it needs `-frames`).
None of the three renders anything.

## Outputs

| mode | written to `-out` |
|---|---|
| chip mode | `chip<n>_summer_s05_s10_realcov_winter.png` panels, `grades.json` (the 8 statistics per panel), `fitted.json` with `-fit` |
| frame mode | the snowified frames under their own names |
| flight mode | the frames under `-pattern`, RGB PNG, same size as the input |

Flight mode also writes, on request:

- `-cloud-masks DIR`: per frame a half-resolution grey PNG of the cloud's
  opacity (255 opaque);
- `-cloud-stats FILE`: a CSV row per frame, appended:

  | column | meaning |
  |---|---|
  | `index`, `t_s` | the pose |
  | `phase` | `below`, `band` or `above` the layer at the camera |
  | `in_cloud` | the camera is inside a cell |
  | `cam_above_base_m` | camera height above the cloud base |
  | `frame_cloud` | share of the frame where the cloud's transmittance is under one half |
  | `frame_opacity` | mean opacity over the frame |
  | `centre_cloud`, `centre_opacity` | the same in a 128 px box at the principal point |
  | `shadow_frac` | share of the visible ground in cloud shadow |

- `-debug DIR`: per frame the snow probability as a grey PNG, and with
  `-debug-raw` the snow model's feature rasters as little-endian float32.

Each frame's log line gives, for snow, its ground resolution, ground
height, snow cover and the shares of plowed road and open water; for
clouds, its cloud share.

## Option reference

Defaults in brackets.

### All modes

| option | meaning |
|---|---|
| `-out DIR` | output directory [`snowify_out`] |
| `-params FILE` | snow appearance parameters, JSON (`params/`) [built-in, `bocage.json`] |

### Chip mode

| option | meaning |
|---|---|
| `-raw DIR` | the pulled Sentinel-2 band GeoTIFFs [`raw`] |
| `-summer ID` | summer scene [`S2B_31TDL_20250807_0_L2A`] |
| `-winter ID` | winter scene [`S2A_31TDL_20240120_0_L2A`] |
| `-autochips` | pick the validation chips from the winter scene's snow and cloud masks instead of the built-in 31TDL ones |
| `-fit N` | refit the appearance parameters with N search iterations [0: no fit] |

### Frame mode

| option | meaning |
|---|---|
| `-frames DIR` | the PNGs to snowify; selects frame mode (with `-poses`: the flight's frames) |
| `-s X` | snow degree 0..1 [1] |
| `-gsd M` | the frames' ground resolution, m/px [10] |

### Flight mode, common

| option | meaning |
|---|---|
| `-poses FILE` | the poses, one per frame; selects flight mode |
| `-weather W` | `snow`, `cloud` or `snow+cloud` [`snow`] |
| `-frames DIR` | the input frames |
| `-calib FILE` | the lens the frames were rendered through |
| `-dem DIR` | GeoTIFF elevation tiles |
| `-pattern P` | frame file name with the pose index [`frame_%05d.png`] |
| `-from N`, `-to N` | first and last pose index, inclusive [0, -1 = the last] |
| `-step N` | every N-th pose from `-from` [1] |
| `-list A,B,..` | just these pose indices, instead of `-from`/`-to` |
| `-workers N` | frames in parallel [the CPUs] |

### Flight mode, snow

| option | meaning |
|---|---|
| `-s-low X`, `-s-high X` | snow degree on low ground and on high ground [0.85, 1] |
| `-h-low M`, `-h-high M` | ground heights between which the degree rises from `-s-low` to `-s-high` [600, 800] |
| `-params-high FILE` | a second parameter set, blended in by the frame's median ground height |
| `-ph-low M`, `-ph-high M` | ground heights where the blend starts and is complete [700, 900] |
| `-a-fixed A` | brightness scale into the calibration: 0 per frame, -1 one for the flight from every 150th frame, >0 this value [0] |
| `-roads FILE` | roads; plowed ones keep their asphalt clear |
| `-plow C=P,..` | plowed share by road class, `all=P` for every class [table above] |
| `-plow-km K` | correlation length of plowing along the network, km [3] |
| `-plow-seed N` | seed of the plowing decision [1] |
| `-road-slack M` | corridor margin beyond the carriageway's half-width, for map error, m [1.5] |
| `-road-chroma C` | largest chroma (max-min)/max of a clear asphalt pixel [0.12] |
| `-water FILE` | water polygons; no snow inside |
| `-calscale` | compute texture and top-hat in calibration units; `-calscale=false` scales them as frame mode's `-gsd` does [true] |
| `-cal-foot X` | luminance smoothing for the top-hat, in 10 m pixels [0.25] |
| `-gate-foot X` | footprint of the luminance gate, in 10 m pixels [0.25] |
| `-smooth-mid X` | relative fine texture below which dark land is a dark field, not forest; 0 off [0.10] |
| `-smooth-width X` | softness of `-smooth-mid` [0.02] |
| `-debug DIR` | write the snow probability per frame here |
| `-debug-raw` | with `-debug`, also the feature rasters as float32 |

### Flight mode, clouds

| option | meaning |
|---|---|
| `-cloud-cover X` | share of the sky opaque from straight above [0.375] |
| `-cloud-base M` | cloud base above the local ground, m [500] |
| `-cloud-top M` | highest tops above the local ground, m [800] |
| `-cloud-cell-km K` | median cumulus size, km [0.6] |
| `-cloud-sigma S` | extinction in a cell's core, 1/m [0.06] |
| `-cloud-smooth-km K` | smoothing of the terrain the layer follows, km [2] |
| `-cloud-seed N` | seed of the cloud field [1] |
| `-cloud-margin-km K` | world kept around the track; the DEM must cover it, km [9] |
| `-cloud-field-m M` | grid of the cloud field, m [10] |
| `-wind N,E` | wind the clouds move with, m/s north and east |
| `-wind-log FILE` | instead of `-wind`: the mean of a wind log inside the layer |
| `-t0 TIME` | UTC time of flight time 0, when the poses carry no `time_utc` |
| `-shadows` | cloud shadows on the ground [true] |
| `-shadow-amb X` | share of the light left under an opaque cell [0.42] |
| `-cloud-masks DIR` | write the cloud opacity per frame here |
| `-cloud-stats FILE` | append the per-frame cloud statistics to this CSV |
| `-profile` | print the layer under the flight, render nothing |
| `-probe A,B,..` | print what the principal point's ray meets in these frames, render nothing |
| `-motion A,B[;C,D]` | check the clouds' motion between frames, render nothing |

### pull

| option | meaning |
|---|---|
| `-bbox W,S,E,N` | region of interest, degrees |
| `-out DIR` | output directory [`raw`] |
| `-winters N` | how many past winters to search [4] |
| `-maxcloud P` | largest cloud cover of the winter scene, % [30] |

## What is not modelled

- Snow: falling snow, snow on trees and roofs as distinct from open land,
  drifts, and melt over the flight. Water stays open everywhere it is
  mapped; ice is not drawn.
- Clouds: growth and decay over the flight, haze, the camera's exposure
  reacting to bright cloud, light scattered under the clouds, and cloud
  shadows on other clouds. The wind is one constant vector for the layer.
- Both: the frames are taken as rendered under a clear summer sky; the
  lighting of the rendered ground is not changed beyond the cloud shadows.

GeoTIFF reading via [github.com/lvdlvd/geotiff](https://github.com/lvdlvd/geotiff);
Sentinel-2 data from the public earth-search STAC (sentinel-cogs).
