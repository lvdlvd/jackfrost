// Camera geometry for flight mode: the lens the frames were rendered
// through, and the attitude that points it. A frame's pixels can only be
// tied to the ground if the ray of each pixel is cast exactly as the
// renderer cast it, so this is the renderer's model, not an approximation.
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

type vec3 [3]float64

func rad(d float64) float64 { return d * math.Pi / 180 }
func deg(r float64) float64 { return r * 180 / math.Pi }

func dot(a, b vec3) float64 { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }

// camBasisNadir returns the camera frame C axes (x=right, y=down, z=boresight)
// in ENU (east,north,up) for a nadir camera at the given heading (deg from N).
// Boresight points straight down; image-up (−y) points along the heading so the
// frame reads "forward-up". roll is assumed 0.
func camBasisNadir(headingDeg float64) (xC, yC, zC vec3) {
	h := rad(headingDeg)
	sh, ch := math.Sin(h), math.Cos(h)
	xC = vec3{ch, -sh, 0}  // right
	yC = vec3{-sh, -ch, 0} // down-in-image (opposite heading)
	zC = vec3{0, 0, -1}    // boresight = down
	return
}

// camBasisBody is the camera frame in ENU for a camera bolted to the
// aircraft looking down through the floor: image right = right wing, image
// top = nose, boresight = body down. Attitude is the aerospace sequence yaw
// (true heading), pitch (nose up positive), roll (right wing down positive),
// in degrees. Level flight gives camBasisNadir(yaw).
//
// Signs, for the reader of a frame: the boresight is the belly, so a right
// bank swings it to the LEFT of the ground track (and the point straight
// below the aircraft moves right of the image centre), and pitch up swings
// it FORWARD (the point below moves toward the bottom edge).
func camBasisBody(rollDeg, pitchDeg, yawDeg float64) (xC, yC, zC vec3) {
	sr, cr := math.Sincos(rad(rollDeg))
	sp, cp := math.Sincos(rad(pitchDeg))
	sy, cy := math.Sincos(rad(yawDeg))
	// body axes in north-east-down: the columns of Rz(yaw) Ry(pitch) Rx(roll)
	nose := vec3{cy * cp, sy * cp, -sp}
	wing := vec3{cy*sp*sr - sy*cr, sy*sp*sr + cy*cr, cp * sr}
	down := vec3{cy*sp*cr + sy*sr, sy*sp*cr - cy*sr, cp * cr}
	enu := func(v vec3) vec3 { return vec3{v[1], v[0], -v[2]} }
	n := enu(nose)
	return enu(wing), vec3{-n[0], -n[1], -n[2]}, enu(down)
}

// Pinhole is a rectilinear lens with Brown-Conrady distortion (OpenCV's k1 k2
// k3 radial, p1 p2 tangential), read from a calibration JSON:
//
//	{"image_size": [w, h], "fx": .., "fy": .., "cx": .., "cy": ..,
//	 "k1": .., "k2": .., "k3": .., "p1": .., "p2": .., "valid_theta_deg": ..}
//
// Edge convention: pixel i spans [i, i+1); the file is in OpenCV's
// convention (pixel centres at integers), shifted by +0.5 on load.
//
// Edge extension. The radial polynomial r_d = r (1 + k1 r^2 + k2 r^4 + k3 r^6),
// r = tan(theta), is a fit to chessboard corners and means nothing beyond the
// outermost corner. A wide lens has k1 < 0, so further out the polynomial
// stops increasing and folds back, and pixels beyond the fold have no ray at
// all.
//
// So beyond theta0 (ValidDeg) the model is continued, not evaluated: the
// distorted radius grows linearly in THETA with the slope d(r_d)/d(theta)
// the polynomial has at theta0, i.e. the lens is taken to be equidistant
// from there on at the scale it was last measured to have. Linear in theta
// and not in tan(theta), because that is what such a lens does: a tangent
// in r would nearly double the pixels per degree by the corner. The
// tangential terms are frozen at their theta0 value along each azimuth. The
// result is continuous with a continuous radial derivative, strictly
// monotonic, and invertible everywhere. It is an extrapolation: geometry
// beyond theta0 is plausible, not calibrated.
//
// theta0 comes from the file's optional "valid_theta_deg" (the field angle
// of the outermost calibration corner). Without it: 0.9 of the fold angle
// when the polynomial folds within 80 deg, else no extension (the plain
// model). SetValidDeg overrides.
type Pinhole struct {
	w, h               int
	fx, fy, cx, cy     float64
	k1, k2, k3, p1, p2 float64
	thetaMax           float64
	// extension: theta0 == 0 means none; r0 = tan(theta0); slope =
	// d(r_d)/d(theta) at theta0
	theta0, r0, slope float64
}

func LoadPinhole(path string) (*Pinhole, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var in struct {
		ImageSize [2]int  `json:"image_size"`
		FX        float64 `json:"fx"`
		FY        float64 `json:"fy"`
		CX        float64 `json:"cx"`
		CY        float64 `json:"cy"`
		K1        float64 `json:"k1"`
		K2        float64 `json:"k2"`
		K3        float64 `json:"k3"`
		P1        float64 `json:"p1"`
		P2        float64 `json:"p2"`
		ValidDeg  float64 `json:"valid_theta_deg"`
	}
	if err := json.Unmarshal(b, &in); err != nil {
		return nil, fmt.Errorf("%s: %v", path, err)
	}
	if in.ImageSize[0] <= 0 || in.ImageSize[1] <= 0 || in.FX <= 0 || in.FY <= 0 {
		return nil, fmt.Errorf("%s: image_size, fx, fy must be positive", path)
	}
	l := &Pinhole{w: in.ImageSize[0], h: in.ImageSize[1],
		fx: in.FX, fy: in.FY, cx: in.CX + 0.5, cy: in.CY + 0.5,
		k1: in.K1, k2: in.K2, k3: in.K3, p1: in.P1, p2: in.P2, thetaMax: rad(80)}
	valid := in.ValidDeg
	if fold := l.FoldDeg(); valid == 0 && fold > 0 {
		valid = 0.9 * fold
	}
	if valid == 0 {
		return l, nil
	}
	if err := l.SetValidDeg(valid); err != nil {
		return nil, fmt.Errorf("%s: %v", path, err)
	}
	return l, nil
}

func (l *Pinhole) Size() (int, int) { return l.w, l.h }

// radial is the polynomial r_d(r) and its derivative in r = tan(theta).
func (l *Pinhole) radial(r float64) (rd, slope float64) {
	r2 := r * r
	return r * (1 + r2*(l.k1+r2*(l.k2+r2*l.k3))),
		1 + r2*(3*l.k1+r2*(5*l.k2+r2*7*l.k3))
}

// FoldDeg is the field angle where the radial polynomial stops increasing,
// 0 if it is monotonic out to the 80 deg the model accepts at all.
func (l *Pinhole) FoldDeg() float64 {
	for t := 0.0; t < l.thetaMax; t += rad(0.01) {
		if _, s := l.radial(math.Tan(t)); s <= 0 {
			return deg(t)
		}
	}
	return 0
}

// ValidDeg is the field angle beyond which the model is extended instead of
// evaluated (see the type comment); 0 when there is no extension.
func (l *Pinhole) ValidDeg() float64 { return deg(l.theta0) }

// SetValidDeg sets that angle. It must lie inside the fold.
func (l *Pinhole) SetValidDeg(d float64) error {
	if fold := l.FoldDeg(); d <= 0 || d >= deg(l.thetaMax) || (fold > 0 && d >= fold) {
		return fmt.Errorf("valid field angle %.2f deg: must be positive and inside the polynomial's fold at %.2f deg", d, fold)
	}
	l.theta0 = rad(d)
	l.r0 = math.Tan(l.theta0)
	_, s := l.radial(l.r0)
	l.slope = s * (1 + l.r0*l.r0)
	return nil
}

// distort is the forward model on the tan plane: ideal (x, y) = d.xy / d.z to
// distorted, both in focal lengths. Pixel is this plus the camera matrix.
func (l *Pinhole) distort(x, y float64) (xd, yd float64) {
	r := math.Hypot(x, y)
	ext := l.theta0 > 0 && r > l.r0
	ux, uy := 0.0, 0.0
	if ext {
		ux, uy = x/r, y/r
		x, y = l.r0*ux, l.r0*uy
	}
	r2 := x*x + y*y
	k := 1 + r2*(l.k1+r2*(l.k2+r2*l.k3))
	xd = x*k + 2*l.p1*x*y + l.p2*(r2+2*x*x)
	yd = y*k + l.p1*(r2+2*y*y) + 2*l.p2*x*y
	if ext {
		e := l.slope * (math.Atan(r) - l.theta0)
		xd, yd = xd+e*ux, yd+e*uy
	}
	return xd, yd
}

// Pixel is the forward model, camera-frame direction to pixel. d need not be
// unit length. front is false behind the lens or more than 80 deg off axis;
// the pixel may lie outside the image.
func (l *Pinhole) Pixel(d vec3) (float64, float64, bool) {
	if d[2] <= 0 || math.Atan2(math.Hypot(d[0], d[1]), d[2]) > l.thetaMax {
		return 0, 0, false
	}
	xd, yd := l.distort(d[0]/d[2], d[1]/d[2])
	return l.cx + l.fx*xd, l.cy + l.fy*yd, true
}

// Ray inverts Pixel. The unknown is (a, b) = theta * (cos az, sin az), in
// which the model is smooth through the axis and through theta0. A radial
// bisection gives the start (the radial part is monotonic by construction),
// Newton with a difference Jacobian then accounts for the tangential terms.
// A pixel further out than the model reaches at 80 deg gets the 80 deg ray.
func (l *Pinhole) Ray(px, py float64) vec3 {
	xd, yd := (px-l.cx)/l.fx, (py-l.cy)/l.fy
	rd := math.Hypot(xd, yd)
	if rd < 1e-12 {
		return vec3{0, 0, 1}
	}
	fwd := func(a, b float64) (float64, float64) {
		t := math.Hypot(a, b)
		if t < 1e-12 {
			return 0, 0
		}
		r := math.Tan(math.Min(t, rad(89)))
		return l.distort(r*a/t, r*b/t)
	}
	lo, hi := 0.0, l.thetaMax
	for i := 0; i < 40; i++ {
		mid := (lo + hi) / 2
		if q, _ := fwd(mid, 0); q < rd {
			lo = mid
		} else {
			hi = mid
		}
	}
	t := (lo + hi) / 2
	a, b := t*xd/rd, t*yd/rd
	const h = 1e-6
	for i := 0; i < 10; i++ {
		fx0, fy0 := fwd(a, b)
		ex, ey := fx0-xd, fy0-yd
		if math.Abs(ex) < 1e-12 && math.Abs(ey) < 1e-12 {
			break
		}
		fxa, fya := fwd(a+h, b)
		fxb, fyb := fwd(a, b+h)
		j00, j01 := (fxa-fx0)/h, (fxb-fx0)/h
		j10, j11 := (fya-fy0)/h, (fyb-fy0)/h
		det := j00*j11 - j01*j10
		if det == 0 {
			break
		}
		a -= (j11*ex - j01*ey) / det
		b -= (j00*ey - j10*ex) / det
	}
	t = math.Min(math.Hypot(a, b), l.thetaMax)
	if t < 1e-12 {
		return vec3{0, 0, 1}
	}
	s := math.Sin(t) / math.Hypot(a, b)
	return vec3{s * a, s * b, math.Cos(t)}
}

// view is a camera at a geodetic position with an aircraft attitude, for the
// cloud renderer, which works in true metres around each camera: the camera
// basis is camBasisBody's, and enu includes the earth's curvature.
type view struct {
	lat, lon, h        float64
	right, down, bore  vec3 // camera axes in local ENU
	mPerDegN, mPerDegE float64
	primeR             float64
}

const (
	wgs84A  = 6378137.0
	wgs84E2 = 6.69437999014e-3
)

// metresPerDeg is the WGS84 north and east metres per degree at lat.
func metresPerDeg(lat float64) (n, e, primeR float64) {
	s := math.Sin(rad(lat))
	w := math.Sqrt(1 - wgs84E2*s*s)
	primeR = wgs84A / w
	meridR := wgs84A * (1 - wgs84E2) / (w * w * w)
	return rad(1) * meridR, rad(1) * primeR * math.Cos(rad(lat)), primeR
}

func newView(lat, lon, h, rollDeg, pitchDeg, yawDeg float64) *view {
	right, down, bore := camBasisBody(rollDeg, pitchDeg, yawDeg)
	n, e, primeR := metresPerDeg(lat)
	return &view{lat: lat, lon: lon, h: h, right: right, down: down, bore: bore,
		mPerDegN: n, mPerDegE: e, primeR: primeR}
}

// enu is the offset of a world point from the camera, metres east, north, up.
func (v *view) enu(lat, lon, h float64) vec3 {
	e := (lon - v.lon) * v.mPerDegE * (1 + h/v.primeR)
	n := (lat - v.lat) * v.mPerDegN * (1 + h/v.primeR)
	return vec3{e, n, h - v.h - (e*e+n*n)/(2*v.primeR)}
}

func (v *view) toCam(p vec3) vec3 {
	return vec3{dot(p, v.right), dot(p, v.down), dot(p, v.bore)}
}

func (v *view) fromCam(d vec3) vec3 {
	return vec3{
		d[0]*v.right[0] + d[1]*v.down[0] + d[2]*v.bore[0],
		d[0]*v.right[1] + d[1]*v.down[1] + d[2]*v.bore[1],
		d[0]*v.right[2] + d[1]*v.down[2] + d[2]*v.bore[2],
	}
}
