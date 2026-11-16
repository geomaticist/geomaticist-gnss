//! Observation models: geometry (range with Sagnac, azimuth / elevation), the standard atmosphere
//! (Saastamoinen + Niell), Sun and Moon, the degree-2 solid tide, and the satellite PCO offset.
//!
//! The ANTEX part lives in `geomaticist-gnss-products`. Operation order is kept fixed so results
//! reproduce to the last bit.
//!
//! References: Saastamoinen (1972); Niell (1996); IERS Conventions 2010 ch. 7.1.1;
//! Montenbruck & Gill (2000) 3.3.2; ANTEX 1.4.

use geomaticist_gnss_core::consts::{CLIGHT, FE_WGS84, OMGE, RE_WGS84};
use geomaticist_gnss_products::AntennaCal;

use crate::linalg::{cross3, dot3, mul3, norm3, sub3};
use crate::nputil::np_interp;

/// Degrees to radians.
pub const D2R: f64 = std::f64::consts::PI / 180.0;

/// ECEF (m) -> latitude, longitude (rad), height (m) on WGS 84 (six iterations).
#[must_use]
pub fn ecef_to_llh(xyz: [f64; 3]) -> (f64, f64, f64) {
    let [x, y, z] = xyz;
    let (a, f) = (RE_WGS84, FE_WGS84);
    let e2 = f * (2.0 - f);
    let p = x.hypot(y);
    let mut lat = z.atan2(p * (1.0 - e2));
    for _ in 0..6 {
        let n = a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
        let h = p / lat.cos() - n;
        lat = z.atan2(p * (1.0 - e2 * n / (n + h)));
    }
    let n = a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    let h = p / lat.cos() - n;
    (lat, y.atan2(x), h)
}

/// Latitude, longitude (rad), height (m) -> ECEF.
#[must_use]
pub fn llh_to_ecef(lat: f64, lon: f64, h: f64) -> [f64; 3] {
    let (a, f) = (RE_WGS84, FE_WGS84);
    let e2 = f * (2.0 - f);
    let n = a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    [
        (n + h) * lat.cos() * lon.cos(),
        (n + h) * lat.cos() * lon.sin(),
        (n * (1.0 - e2) + h) * lat.sin(),
    ]
}

/// Rows E, N, U: `ENU = R · dXYZ`.
#[must_use]
pub fn enu_matrix(lat: f64, lon: f64) -> [[f64; 3]; 3] {
    let (sl, cl, sp, cp) = (lon.sin(), lon.cos(), lat.sin(), lat.cos());
    [
        [-sl, cl, 0.0],
        [-sp * cl, -sp * sl, cp],
        [cp * cl, cp * sl, sp],
    ]
}

/// Range with the Sagnac term and the unit vector receiver -> satellite (ECEF), for
/// one satellite / receiver pair.
#[inline]
#[must_use]
pub fn geometry(sat: [f64; 3], rec: [f64; 3]) -> (f64, [f64; 3]) {
    let d = sub3(sat, rec);
    let r = norm3(d);
    let e = [d[0] / r, d[1] / r, d[2] / r];
    let sag = OMGE * (sat[0] * rec[1] - sat[1] * rec[0]) / CLIGHT;
    (r + sag, e)
}

/// Azimuth, elevation (rad) and the ENU unit vector from an ECEF unit vector and the station's ENU
/// matrix.
#[inline]
#[must_use]
pub fn azel(e: [f64; 3], r: &[[f64; 3]; 3]) -> (f64, f64, [f64; 3]) {
    // enu = e @ R.T: enu[k] = sum_i e[i] * R[k][i]
    let enu = [
        e[0] * r[0][0] + e[1] * r[0][1] + e[2] * r[0][2],
        e[0] * r[1][0] + e[1] * r[1][1] + e[2] * r[1][2],
        e[0] * r[2][0] + e[1] * r[2][1] + e[2] * r[2][2],
    ];
    let mut az = enu[0].atan2(enu[1]);
    if az < 0.0 {
        az += 2.0 * std::f64::consts::PI;
    }
    let el = enu[2].clamp(-1.0, 1.0).asin();
    (az, el, enu)
}

/// Zenith hydrostatic (+ wet, for humidity > 0) delay (m) of the standard atmosphere (RTKLIB
/// tropmodel).
#[must_use]
pub fn saastamoinen_zenith(lat: f64, h: f64, humidity: f64) -> (f64, f64) {
    if !(-100.0..=1e4).contains(&h) {
        return (0.0, 0.0);
    }
    let hgt = h.max(0.0);
    let pres = 1013.25 * (1.0 - 2.2557e-5 * hgt).powf(5.2568);
    let temp = 15.0 - 6.5e-3 * hgt + 273.16;
    let e = 6.108 * humidity * ((17.15 * temp - 4684.0) / (temp - 38.45)).exp();
    let zhd = 0.002_276_8 * pres / (1.0 - 0.002_66 * (2.0 * lat).cos() - 0.000_28 * hgt / 1e3);
    let zwd = 0.002_277 * (1255.0 / temp + 0.05) * e;
    (zhd, zwd)
}

const NMF_LAT: [f64; 5] = [15.0, 30.0, 45.0, 60.0, 75.0];
const NMF_H_AVG: [[f64; 5]; 3] = [
    [
        1.276_993_4e-3,
        1.268_323_0e-3,
        1.246_539_7e-3,
        1.219_604_9e-3,
        1.204_599_6e-3,
    ],
    [
        2.915_369_5e-3,
        2.915_229_9e-3,
        2.928_844_5e-3,
        2.902_256_5e-3,
        2.902_491_2e-3,
    ],
    [
        62.610_505e-3,
        62.837_393e-3,
        63.721_774e-3,
        63.824_265e-3,
        64.258_455e-3,
    ],
];
const NMF_H_AMP: [[f64; 5]; 3] = [
    [
        0.0,
        1.270_962_6e-5,
        2.652_366_2e-5,
        3.400_045_2e-5,
        4.120_219_1e-5,
    ],
    [
        0.0,
        2.141_497_9e-5,
        3.016_077_9e-5,
        7.256_272_2e-5,
        11.723_375e-5,
    ],
    [
        0.0,
        9.012_840_0e-5,
        4.349_703_7e-5,
        84.795_348e-5,
        170.372_06e-5,
    ],
];
const NMF_W: [[f64; 5]; 3] = [
    [
        5.802_189_7e-4,
        5.679_484_7e-4,
        5.811_801_9e-4,
        5.972_754_2e-4,
        6.164_169_3e-4,
    ],
    [
        1.427_526_8e-3,
        1.513_862_5e-3,
        1.457_275_2e-3,
        1.500_742_8e-3,
        1.759_908_2e-3,
    ],
    [
        4.347_296_1e-2,
        4.672_951_0e-2,
        4.390_893_1e-2,
        4.462_698_2e-2,
        5.473_603_8e-2,
    ],
];
/// Niell's height-correction coefficients.
pub const NMF_HT: (f64, f64, f64) = (2.53e-5, 5.49e-3, 1.14e-3);

/// The continued-fraction mapping function.
#[inline]
#[must_use]
pub fn mapf(el: f64, a: f64, b: f64, c: f64) -> f64 {
    let s = el.sin();
    (1.0 + a / (1.0 + b / (1.0 + c))) / (s + a / (s + b / (s + c)))
}

/// Niell (1996) hydrostatic and wet mapping coefficients at a site and day: `(ah[3], aw[3])`.
#[must_use]
pub fn niell_coeffs(lat: f64, doy: f64) -> ([f64; 3], [f64; 3]) {
    let alat = lat.abs() / D2R;
    let y = (doy - 28.0) / 365.25 + if lat < 0.0 { 0.5 } else { 0.0 };
    let cosy = (2.0 * std::f64::consts::PI * y).cos();
    let mut ah = [0.0; 3];
    let mut aw = [0.0; 3];
    for i in 0..3 {
        ah[i] = np_interp(alat, &NMF_LAT, &NMF_H_AVG[i])
            - np_interp(alat, &NMF_LAT, &NMF_H_AMP[i]) * cosy;
        aw[i] = np_interp(alat, &NMF_LAT, &NMF_W[i]);
    }
    (ah, aw)
}

/// Niell hydrostatic (with the height correction) and wet mapping at one elevation (rad), from the
/// coefficients of [`niell_coeffs`].
#[inline]
#[must_use]
pub fn niell_at(el: f64, ah: &[f64; 3], aw: &[f64; 3], h: f64) -> (f64, f64) {
    let el = el.max(1e-3);
    let dm = (1.0 / el.sin() - mapf(el, NMF_HT.0, NMF_HT.1, NMF_HT.2)) * h / 1e3;
    (
        mapf(el, ah[0], ah[1], ah[2]) + dm,
        mapf(el, aw[0], aw[1], aw[2]),
    )
}

/// Niell mapping functions for elevations (rad). `doy` may be fractional (an integer gives the same
/// arithmetic).
#[must_use]
pub fn niell(el: &[f64], lat: f64, h: f64, doy: f64) -> (Vec<f64>, Vec<f64>) {
    let (ah, aw) = niell_coeffs(lat, doy);
    let mut mh = Vec::with_capacity(el.len());
    let mut mw = Vec::with_capacity(el.len());
    for &e in el {
        let (a, b) = niell_at(e, &ah, &aw, h);
        mh.push(a);
        mw.push(b);
    }
    (mh, mw)
}

/// Greenwich mean sidereal angle (rad) at GPS seconds (UT1 ≈ UTC, GPS − UTC taken as 18 s:
/// a 0.02 mm tide effect before 2017).
#[must_use]
pub fn gmst(t_gps: f64) -> f64 {
    let jd = 2_444_244.5 + (t_gps - 18.0) / 86400.0;
    ((280.460_618_37 + 360.985_647_366_29 * (jd - 2_451_545.0)).rem_euclid(360.0)).to_radians()
}

/// Low-precision Sun and Moon positions, ECEF (m), at one GPS time (Montenbruck & Gill 3.3.2).
#[must_use]
pub fn sun_moon_ecef(t: f64) -> ([f64; 3], [f64; 3]) {
    let jd_tt = 2_444_244.5 + (t + 51.184) / 86400.0;
    let tc = (jd_tt - 2_451_545.0) / 36525.0;
    let eps = 23.439_291_11_f64.to_radians();
    let prec = (1.3972 * tc).to_radians();
    // Sun
    let m = (357.5256 + 35999.049 * tc).to_radians();
    let lam = 282.9400_f64.to_radians()
        + m
        + ((6892.0 * m.sin() + 72.0 * (2.0 * m).sin()) / 3600.0).to_radians()
        + prec;
    let r = (149.619 - 2.499 * m.cos() - 0.021 * (2.0 * m).cos()) * 1e9;
    let sun = [
        r * lam.cos(),
        r * lam.sin() * eps.cos(),
        r * lam.sin() * eps.sin(),
    ];
    // Moon
    let l0 = 218.31617 + 481_267.88088 * tc;
    let l = (134.96292 + 477_198.86753 * tc).to_radians();
    let lp = (357.52543 + 35999.04944 * tc).to_radians();
    let f = (93.27283 + 483_202.01873 * tc).to_radians();
    let d = (297.85027 + 445_267.11135 * tc).to_radians();
    let dl = 22640.0 * l.sin() + 769.0 * (2.0 * l).sin() - 4586.0 * (l - 2.0 * d).sin()
        + 2370.0 * (2.0 * d).sin()
        - 668.0 * lp.sin()
        - 412.0 * (2.0 * f).sin()
        - 212.0 * (2.0 * l - 2.0 * d).sin()
        - 206.0 * (l + lp - 2.0 * d).sin()
        + 192.0 * (l + 2.0 * d).sin()
        - 165.0 * (lp - 2.0 * d).sin()
        + 148.0 * (l - lp).sin()
        - 125.0 * d.sin()
        - 110.0 * (l + lp).sin()
        - 55.0 * (2.0 * f - 2.0 * d).sin();
    let lamm = (l0 + dl / 3600.0).to_radians() + prec;
    let beta = ((18520.0
        * (f + (dl / 3600.0).to_radians()
            + ((412.0 * (2.0 * f).sin() + 541.0 * lp.sin()) / 3600.0).to_radians())
        .sin()
        - 526.0 * (f - 2.0 * d).sin()
        + 44.0 * (l + f - 2.0 * d).sin()
        - 31.0 * (-l + f - 2.0 * d).sin()
        - 25.0 * (-2.0 * l + f).sin()
        - 23.0 * (lp + f - 2.0 * d).sin()
        + 21.0 * (-l + f).sin()
        + 11.0 * (-lp + f - 2.0 * d).sin())
        / 3600.0)
        .to_radians();
    let rm = (385_000.0
        - 20905.0 * l.cos()
        - 3699.0 * (2.0 * d - l).cos()
        - 2956.0 * (2.0 * d).cos()
        - 570.0 * (2.0 * l).cos()
        + 246.0 * (2.0 * l - 2.0 * d).cos()
        - 205.0 * (lp - 2.0 * d).cos()
        - 171.0 * (l + 2.0 * d).cos()
        - 152.0 * (l + lp - 2.0 * d).cos())
        * 1e3;
    let xe = rm * beta.cos() * lamm.cos();
    let ye = rm * beta.cos() * lamm.sin();
    let ze = rm * beta.sin();
    let moon = [
        xe,
        ye * eps.cos() - ze * eps.sin(),
        ye * eps.sin() + ze * eps.cos(),
    ];
    let g = gmst(t);
    let (cg, sg) = (g.cos(), g.sin());
    let rot = |v: [f64; 3]| [cg * v[0] + sg * v[1], -sg * v[0] + cg * v[1], v[2]];
    (rot(sun), rot(moon))
}

/// Sun positions at every time.
#[must_use]
pub fn sun_ecef(t: &[f64]) -> Vec<[f64; 3]> {
    t.iter().map(|&x| sun_moon_ecef(x).0).collect()
}

/// GM of the Earth, Sun and Moon the degree-2 tide uses.
pub const GM_EARTH: f64 = 3.986_004_415e14;
/// Sun.
pub const GM_SUN: f64 = 1.327_124_420_76e20;
/// Moon.
pub const GM_MOON: f64 = 4.9028e12;

/// Solid-earth tide displacement (m, ECEF) at each time: IERS 2010 step 1, degree 2, with the latitude
/// dependence of h2 and l2 (no step-2 corrections, no permanent-tide removal: tide-free, as ITRF).
/// (The full IERS model is [`crate::tides::solid_tide`]; B9.)
#[must_use]
pub fn solid_tide_degree2(xyz: [f64; 3], t: &[f64]) -> Vec<[f64; 3]> {
    let rn = norm3(xyz);
    let ru = [xyz[0] / rn, xyz[1] / rn, xyz[2] / rn];
    let lat = ru[2].asin();
    let p2 = (3.0 * lat.sin().powi(2) - 1.0) / 2.0;
    let h2 = 0.6078 - 0.0006 * p2;
    let l2 = 0.0847 + 0.0002 * p2;
    t.iter()
        .map(|&tt| {
            let (sun, moon) = sun_moon_ecef(tt);
            let mut out = [0.0; 3];
            for (body, gm) in [(sun, GM_SUN), (moon, GM_MOON)] {
                let rb = norm3(body);
                let bu = [body[0] / rb, body[1] / rb, body[2] / rb];
                let k = gm * rn.powi(4) / (GM_EARTH * rb.powi(3));
                let c = dot3(bu, ru);
                let f1 = h2 * (1.5 * c * c - 0.5);
                for i in 0..3 {
                    out[i] += k * (f1 * ru[i] + 3.0 * l2 * c * (bu[i] - c * ru[i]));
                }
            }
            out
        })
        .collect()
}

/// Centre of mass -> ionosphere-free antenna phase centre (m, ECEF), nominal attitude, for one
/// satellite position and Sun position.
#[must_use]
pub fn sat_pco_offset(
    cal: Option<&AntennaCal>,
    sys: u8,
    b1: u8,
    b2: u8,
    sat: [f64; 3],
    sun: [f64; 3],
    a: f64,
    b: f64,
) -> [f64; 3] {
    let Some(cal) = cal else {
        return [0.0; 3];
    };
    let (Some(k1), Some(k2)) = (cal.freq_key(sys, b1), cal.freq_key(sys, b2)) else {
        return [0.0; 3];
    };
    let p1 = cal.pco[&k1];
    let p2 = cal.pco[&k2];
    let off = [
        a * p1[0] - b * p2[0],
        a * p1[1] - b * p2[1],
        a * p1[2] - b * p2[2],
    ];
    let (ex, ey, ez) = crate::attitude::sat_axes(sat, sun);
    let mut out = mul3(ex, off[0]);
    let ey = mul3(ey, off[1]);
    let ez = mul3(ez, off[2]);
    for i in 0..3 {
        out[i] = out[i] + ey[i] + ez[i];
    }
    out
}

/// `ez` of the body frame: `-r / |r|` (shared with the attitude helpers).
#[inline]
#[must_use]
pub fn nadir_unit(sat: [f64; 3]) -> [f64; 3] {
    let n = norm3(sat);
    [-sat[0] / n, -sat[1] / n, -sat[2] / n]
}

/// `(sun - sat) / |sun - sat|`.
#[inline]
#[must_use]
pub fn sun_unit(sat: [f64; 3], sun: [f64; 3]) -> [f64; 3] {
    let es = sub3(sun, sat);
    let n = norm3(es);
    [es[0] / n, es[1] / n, es[2] / n]
}

/// Unit cross product helper.
#[inline]
#[must_use]
pub fn unit_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let c = cross3(a, b);
    let n = norm3(c);
    [c[0] / n, c[1] / n, c[2] / n]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geodesy_round_trip() {
        let xyz = llh_to_ecef(64.978_f64.to_radians(), (-147.499_f64).to_radians(), 319.0);
        let (lat, lon, h) = ecef_to_llh(xyz);
        assert!((lat.to_degrees() - 64.978).abs() < 1e-11);
        assert!((lon.to_degrees() + 147.499).abs() < 1e-11);
        assert!((h - 319.0).abs() < 1e-6);
        let r = enu_matrix(lat, lon);
        let (rng, e) = geometry([15600e3, 7540e3, 20140e3], xyz);
        assert!(rng > 2.0e7 && rng < 2.6e7);
        let (az, el, _) = azel(e, &r);
        assert!((0.0..2.0 * std::f64::consts::PI).contains(&az) && el.abs() <= 1.6);
    }

    #[test]
    fn tropo_and_bodies() {
        let (zhd, zwd) = saastamoinen_zenith(1.1, 300.0, 0.0);
        assert!((zhd - 2.2).abs() < 0.1 && zwd == 0.0);
        let (mh, mw) = niell(&[0.5, 1.0], 1.1, 300.0, 267.0);
        assert!(mh[0] > mh[1] && mw[0] > mw[1] && mh[1] > 1.0);
        let (sun, moon) = sun_moon_ecef(1_411_084_800.0);
        assert!((norm3(sun) / 1.496e11 - 1.0).abs() < 0.03);
        assert!((norm3(moon) / 3.84e8 - 1.0).abs() < 0.08);
        let tide = solid_tide_degree2([4_395_470.63, 617_743.11, 4_565_607.21], &[1_411_084_800.0]);
        assert!(norm3(tide[0]) < 0.5);
    }
}
