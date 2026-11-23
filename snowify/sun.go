// The sun for the cloud renderer's light and shadows.
package main

import "math"

// sunPosition is the sun's elevation and azimuth (degrees, azimuth from true
// north through east) at a unix time and place: the US Naval Observatory's
// low-precision almanac (about 0.01 deg in this century, ignoring
// refraction), which is far finer than a cloud shadow needs.
func sunPosition(unix float64, lat, lon float64) (elDeg, azDeg float64) {
	n := unix/86400 + 2440587.5 - 2451545.0 // days from J2000.0
	L := math.Mod(280.460+0.9856474*n, 360)
	g := rad(math.Mod(357.528+0.9856003*n, 360))
	lambda := rad(L + 1.915*math.Sin(g) + 0.020*math.Sin(2*g))
	eps := rad(23.439 - 0.0000004*n)
	ra := math.Atan2(math.Cos(eps)*math.Sin(lambda), math.Cos(lambda))
	dec := math.Asin(math.Sin(eps) * math.Sin(lambda))
	gmst := math.Mod(18.697374558+24.06570982441908*n, 24)
	lha := rad(gmst*15+lon) - ra
	phi := rad(lat)
	el := math.Asin(math.Sin(phi)*math.Sin(dec) + math.Cos(phi)*math.Cos(dec)*math.Cos(lha))
	az := math.Atan2(-math.Sin(lha), math.Tan(dec)*math.Cos(phi)-math.Sin(phi)*math.Cos(lha))
	return deg(el), math.Mod(deg(az)+360, 360)
}
