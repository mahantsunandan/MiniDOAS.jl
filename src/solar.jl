# This file is part of MiniDOAS.jl (Differential Optical Absorption Spectroscopy in Julia).
# Copyright (c) 2026 Sunandan Mahant <sunandanmahant@outlook.com>; BSD-3-Clause, see LICENSE.

# =============================================================================
# Solar position and geometric air mass factors
# =============================================================================

"""
    solar_position(t::DateTime, latitude, longitude; refraction=false)
    solar_position(s::MFCSpectrum; refraction=false)

The position of the sun at the time `t` (UTC) seen from `latitude`, `longitude` (degrees,
north and east positive), with the NOAA solar calculator's algorithm (after Meeus,
*Astronomical Algorithms*; about 0.01° from 1800 to 2100). Returns a named tuple:

  * `zenith`      solar zenith angle (SZA), degrees
  * `elevation`   90° - zenith
  * `azimuth`     degrees clockwise from north
  * `declination` degrees
  * `equation_of_time` minutes

With `refraction=true` zenith and elevation include the NOAA estimate of atmospheric
refraction (the apparent position); by default they are geometric, as usually reported
with DOAS measurements. For an `MFCSpectrum`, its mid-exposure time and its latitude and
longitude are used.
"""
function solar_position(t::DateTime, latitude::Real, longitude::Real; refraction::Bool=false)
    -90 <= latitude <= 90 || throw(ArgumentError("latitude $latitude outside [-90, 90]"))
    jd = datetime2julian(t)
    jc = (jd - 2451545.0) / 36525.0
    l0 = mod(280.46646 + jc * (36000.76983 + jc * 0.0003032), 360.0)    # geometric mean longitude
    m = 357.52911 + jc * (35999.05029 - 0.0001537 * jc)                  # mean anomaly
    e = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)              # eccentricity of the orbit
    c = sind(m) * (1.914602 - jc * (0.004817 + 0.000014 * jc)) + sind(2m) * (0.019993 - 0.000101 * jc) +
        sind(3m) * 0.000289                                              # equation of the centre
    ω = 125.04 - 1934.136 * jc
    app_long = l0 + c - 0.00569 - 0.00478 * sind(ω)                     # apparent longitude
    ε0 = 23 + (26 + (21.448 - jc * (46.815 + jc * (0.00059 - jc * 0.001813))) / 60) / 60
    ε = ε0 + 0.00256 * cosd(ω)                                           # corrected obliquity
    decl = asind(sind(ε) * sind(app_long))
    y = tand(ε / 2)^2
    eqtime = 4 * rad2deg(y * sind(2l0) - 2e * sind(m) + 4e * y * sind(m) * cosd(2l0) -
                         0.5 * y^2 * sind(4l0) - 1.25 * e^2 * sind(2m))
    minutes = (Dates.value(t - DateTime(Date(t))) / 60_000)            # UTC minutes of the day
    tst = mod(minutes + eqtime + 4 * longitude, 1440.0)                  # true solar time
    ha = tst / 4 < 0 ? tst / 4 + 180 : tst / 4 - 180                     # hour angle
    cosz = sind(latitude) * sind(decl) + cosd(latitude) * cosd(decl) * cosd(ha)
    zen = acosd(clamp(cosz, -1.0, 1.0))
    den = cosd(latitude) * sind(zen)
    az = if abs(den) < 1e-12
        latitude > 0 ? 180.0 : 0.0
    else
        a = acosd(clamp((sind(latitude) * cosd(zen) - sind(decl)) / den, -1.0, 1.0))
        ha > 0 ? mod(a + 180, 360) : mod(540 - a, 360)
    end
    if refraction
        el = 90 - zen
        te = tand(el)
        r = el > 85 ? 0.0 :
            el > 5 ? 58.1 / te - 0.07 / te^3 + 0.000086 / te^5 :
            el > -0.575 ? 1735 + el * (-518.2 + el * (103.4 + el * (-12.79 + el * 0.711))) :
            -20.772 / te
        zen -= r / 3600
    end
    (zenith=zen, elevation=90 - zen, azimuth=az, declination=decl, equation_of_time=eqtime)
end

function solar_position(s::MFCSpectrum; refraction::Bool=false)
    (isnan(s.latitude) || isnan(s.longitude)) && throw(ArgumentError("the spectrum $(s.name) has no position"))
    solar_position(measurement_time(s), s.latitude, s.longitude; refraction=refraction)
end

"""
    geometric_amf(elevation)

The geometric air mass factor 1/sin(α) of a tropospheric absorber seen at the viewing
elevation α (degrees) - the single-scattering approximation used for MAX-DOAS (e.g.
Hönninger et al., 2004).
"""
geometric_amf(elevation::Real) = 1 / sind(elevation)

"""
    sza_amf(sza)

The geometric air mass factor 1/cos(θ) for the solar zenith angle θ (degrees): direct-sun
observations, or roughly a stratospheric absorber seen in zenith-scattered light at
moderate θ.
"""
sza_amf(sza::Real) = 1 / cosd(sza)

"""
    geometric_vcd(dscd, elevation)

The tropospheric vertical column from a MAX-DOAS differential slant column (measured at
elevation α against a zenith reference) in the geometric approximation:
dSCD / (1/sin α - 1). Reasonable for α of about 15-30° in a well-mixed boundary layer
(Brinksma et al., 2008); NaN for α ≥ 90° or α ≤ 0°.
"""
function geometric_vcd(dscd::Real, elevation::Real)
    (elevation <= 0 || elevation >= 90) && return NaN
    dscd / (1 / sind(elevation) - 1)
end
