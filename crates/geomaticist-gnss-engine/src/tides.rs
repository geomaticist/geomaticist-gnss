//! The full IERS 2010 solid-earth tide (steps 1 and 2), the pole tide and ocean tide loading.
//!
//! After the IERS routines DEHANTTIDEINEL, ST1IDIU, ST1ISEM, ST1L1, STEP2DIU and
//! STEP2LON (Dehant, Mathews, Petit, Bizouard; IERS Conventions Software Licence): a re-implementation,
//! not IERS software, with the engine's own low-precision Sun / Moon. Ocean loading arguments as RTKLIB's
//! `tide_oload` (T. Takasu, BSD-2) and the IERS ARG2 routine.
//!
//! The static engine uses this model under `Compat::tide = TideModel::IersFull`; the original
//! degree-2 step-1 model is [`crate::models::solid_tide_degree2`].

use geomaticist_gnss_core::time::leap_at;

use crate::linalg::{dot3, norm3};
use crate::models::{ecef_to_llh, enu_matrix, sun_moon_ecef};

/// Earth radius of the IERS routines.
pub const RE_IERS: f64 = 6_378_136.6;
/// Sun / Earth mass ratio.
pub const MASS_SUN: f64 = 332_946.0482;
/// Moon / Earth mass ratio.
pub const MASS_MOON: f64 = 0.012_300_037_1;

/// Julian centuries of TT since J2000 (TT = GPS + 51.184 s).
#[must_use]
pub fn tt_centuries(t_gps: f64) -> f64 {
    let jd = 2_444_244.5 + (t_gps + 51.184) / 86400.0;
    (jd - 2_451_545.0) / 36525.0
}

/// Hours of the UTC day (UT1 ≈ UTC).
#[must_use]
pub fn ut_hours(t_gps: f64, leap: i32) -> f64 {
    (t_gps - f64::from(leap)).rem_euclid(86400.0) / 3600.0
}

fn neu_to_xyz(
    dr: f64,
    dn: f64,
    de: f64,
    sinphi: f64,
    cosphi: f64,
    sinla: f64,
    cosla: f64,
) -> [f64; 3] {
    [
        dr * cosla * cosphi - de * sinla - dn * sinphi * cosla,
        dr * sinla * cosphi + de * cosla - dn * sinphi * sinla,
        dr * sinphi + dn * cosphi,
    ]
}

/// Step 1 of the IERS tide for one epoch.
#[must_use]
pub fn step1(xsta: [f64; 3], xsun: [f64; 3], xmon: [f64; 3]) -> [f64; 3] {
    let rsta = norm3(xsta);
    let rsun = norm3(xsun);
    let rmon = norm3(xmon);
    let scsun = dot3(xsun, xsta) / rsta / rsun;
    let scmon = dot3(xmon, xsta) / rsta / rmon;
    let cosphi = xsta[0].hypot(xsta[1]) / rsta;
    let sinphi = xsta[2] / rsta;
    let sinla = xsta[1] / cosphi / rsta;
    let cosla = xsta[0] / cosphi / rsta;
    let (h20, l20, h3, l3) = (0.6078, 0.0847, 0.292, 0.015);
    let h2 = h20 - 0.0006 * (1.0 - 1.5 * cosphi * cosphi);
    let l2 = l20 + 0.0002 * (1.0 - 1.5 * cosphi * cosphi);
    let fac2sun = MASS_SUN * RE_IERS * (RE_IERS / rsun).powi(3);
    let fac2mon = MASS_MOON * RE_IERS * (RE_IERS / rmon).powi(3);
    let fac3sun = fac2sun * (RE_IERS / rsun);
    let fac3mon = fac2mon * (RE_IERS / rmon);
    let mut out = [0.0; 3];
    let usta = [xsta[0] / rsta, xsta[1] / rsta, xsta[2] / rsta];
    for (sc, xb, rb, f2, f3) in [
        (scsun, xsun, rsun, fac2sun, fac3sun),
        (scmon, xmon, rmon, fac2mon, fac3mon),
    ] {
        let p2 = 3.0 * (h2 / 2.0 - l2) * sc * sc - h2 / 2.0;
        let p3 = 2.5 * (h3 - 3.0 * l3) * sc * sc * sc + 1.5 * (l3 - h3) * sc;
        let x2 = 3.0 * l2 * sc;
        let x3 = 1.5 * l3 * (5.0 * sc * sc - 1.0);
        for i in 0..3 {
            out[i] += (f2 * x2) * xb[i] / rb + (f2 * p2) * usta[i];
        }
        for i in 0..3 {
            out[i] += (f3 * x3) * xb[i] / rb + (f3 * p3) * usta[i];
        }
    }
    let cos2phi = cosphi * cosphi - sinphi * sinphi;
    for (xb, rb, f2) in [(xsun, rsun, fac2sun), (xmon, rmon, fac2mon)] {
        let rb2 = rb * rb;
        let q = xb[2] * (xb[0] * sinla - xb[1] * cosla) / rb2;
        let qe = xb[2] * (xb[0] * cosla + xb[1] * sinla) / rb2;
        let (dhi, dli) = (-0.0025, -0.0007);
        let dr = -3.0 * dhi * sinphi * cosphi * f2 * q;
        let dn = -3.0 * dli * cos2phi * f2 * q;
        let de = -3.0 * dli * sinphi * f2 * qe;
        let v = neu_to_xyz(dr, dn, de, sinphi, cosphi, sinla, cosla);
        for i in 0..3 {
            out[i] += v[i];
        }
        let c2l = cosla * cosla - sinla * sinla;
        let s2l = 2.0 * cosla * sinla;
        let u = ((xb[0] * xb[0] - xb[1] * xb[1]) * s2l - 2.0 * xb[0] * xb[1] * c2l) / rb2;
        let w = ((xb[0] * xb[0] - xb[1] * xb[1]) * c2l + 2.0 * xb[0] * xb[1] * s2l) / rb2;
        let (dhi, dli) = (-0.0022, -0.0007);
        let dr = -0.75 * dhi * cosphi * cosphi * f2 * u;
        let dn = 1.5 * dli * sinphi * cosphi * f2 * u;
        let de = -1.5 * dli * cosphi * f2 * w;
        let v = neu_to_xyz(dr, dn, de, sinphi, cosphi, sinla, cosla);
        for i in 0..3 {
            out[i] += v[i];
        }
        let l1 = 0.0012;
        let dn = 3.0 * (-l1 * sinphi * sinphi * f2 * qe);
        let de = 3.0 * (l1 * sinphi * cos2phi * f2 * q);
        let v = neu_to_xyz(0.0, dn, de, sinphi, cosphi, sinla, cosla);
        for i in 0..3 {
            out[i] += v[i];
        }
        let l1 = 0.0024;
        let dn = 3.0 * (-l1 / 2.0 * sinphi * cosphi * f2 * w);
        let de = 3.0 * (-l1 / 2.0 * sinphi * sinphi * cosphi * f2 * u);
        let v = neu_to_xyz(0.0, dn, de, sinphi, cosphi, sinla, cosla);
        for i in 0..3 {
            out[i] += v[i];
        }
    }
    out
}

/// Diurnal terms (multipliers of s, h, p, N', ps; radial in / out of phase; north-east in / out), mm.
const DIU: [[f64; 9]; 31] = [
    [-3.0, 0.0, 2.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
    [-3.0, 2.0, 0.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
    [-2.0, 0.0, 1.0, -1.0, 0.0, -0.02, 0.0, 0.0, 0.0],
    [-2.0, 0.0, 1.0, 0.0, 0.0, -0.08, 0.0, -0.01, 0.01],
    [-2.0, 2.0, -1.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0, -1.0, 0.0, -0.10, 0.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0, 0.0, 0.0, -0.51, 0.0, -0.02, 0.03],
    [-1.0, 2.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
    [0.0, -2.0, 1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
    [0.0, 0.0, -1.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 0.0, 0.06, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 1.0, 0.0, 0.01, 0.0, 0.0, 0.0],
    [0.0, 2.0, -1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
    [1.0, -3.0, 0.0, 0.0, 1.0, -0.06, 0.0, 0.0, 0.0],
    [1.0, -2.0, 0.0, -1.0, 0.0, 0.01, 0.0, 0.0, 0.0],
    [1.0, -2.0, 0.0, 0.0, 0.0, -1.23, -0.07, 0.06, 0.01],
    [1.0, -1.0, 0.0, 0.0, -1.0, 0.02, 0.0, 0.0, 0.0],
    [1.0, -1.0, 0.0, 0.0, 1.0, 0.04, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, -1.0, 0.0, -0.22, 0.01, 0.01, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0, 12.00, -0.80, -0.67, -0.03],
    [1.0, 0.0, 0.0, 1.0, 0.0, 1.73, -0.12, -0.10, 0.0],
    [1.0, 0.0, 0.0, 2.0, 0.0, -0.04, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, -1.0, -0.50, -0.01, 0.03, 0.0],
    [1.0, 1.0, 0.0, 0.0, 1.0, 0.01, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 1.0, -1.0, -0.01, 0.0, 0.0, 0.0],
    [1.0, 2.0, -2.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
    [1.0, 2.0, 0.0, 0.0, 0.0, -0.11, 0.01, 0.01, 0.0],
    [2.0, -2.0, 1.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
    [2.0, 0.0, -1.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
];
/// Long-period terms.
const LON: [[f64; 9]; 5] = [
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.47, 0.23, 0.16, 0.07],
    [0.0, 2.0, 0.0, 0.0, 0.0, -0.20, -0.12, -0.11, -0.05],
    [1.0, 0.0, -1.0, 0.0, 0.0, -0.11, -0.08, -0.09, -0.04],
    [2.0, 0.0, 0.0, 0.0, 0.0, -0.13, -0.11, -0.15, -0.07],
    [2.0, 0.0, 0.0, 1.0, 0.0, -0.05, -0.05, -0.06, -0.03],
];

/// Mean longitudes (degrees): `(s0, s, h, p, N', ps, tau_no_s)`.
#[must_use]
pub fn fund_args(t: f64) -> [f64; 7] {
    let s = 218.316_645_63 + (481_267.881_94 + (-0.001_466_388_9 + 0.000_001_851_39 * t) * t) * t;
    let pr = (1.396_971_278 + (0.000_308_889 + (0.000_000_021 + 0.000_000_007 * t) * t) * t) * t;
    let h = 280.46645
        + (36000.769_748_9 + (0.000_303_222_22 + (0.000_000_020 - 0.000_000_006_54 * t) * t) * t)
            * t;
    let p = 83.353_243_12
        + (4069.013_635_25
            + (-0.010_321_722_22 + (-0.000_012_499_1 + 0.000_000_052_63 * t) * t) * t)
            * t;
    let zns = 234.955_444_99
        + (1934.136_261_97
            + (-0.002_075_611_11 + (-0.000_002_139_44 + 0.000_000_016_50 * t) * t) * t)
            * t;
    let ps = 282.937_340_98
        + (1.719_457_666_67
            + (0.000_456_888_89 + (-0.000_000_017_78 - 0.000_000_003_34 * t) * t) * t)
            * t;
    let tau_no_s = 280.460_618_4 + (36000.770_053_6 + (0.000_387_93 - 0.000_000_025_8 * t) * t) * t;
    let m = |x: f64| x.rem_euclid(360.0);
    [m(s), m(s + pr), m(h), m(p), m(zns), m(ps), tau_no_s]
}

/// Step 2 (frequency-dependent corrections) for one epoch: `t` TT centuries, `fhr` UT hours.
#[must_use]
pub fn step2(xsta: [f64; 3], t: f64, fhr: f64) -> [f64; 3] {
    let [s0, s, h, p, zns, ps, tau0] = fund_args(t);
    let tau = (fhr * 15.0 + tau0 - s0).rem_euclid(360.0);
    let rsta = norm3(xsta);
    let sinphi = xsta[2] / rsta;
    let cosphi = xsta[0].hypot(xsta[1]) / rsta;
    let cosla = xsta[0] / cosphi / rsta;
    let sinla = xsta[1] / cosphi / rsta;
    let zla = xsta[1].atan2(xsta[0]);
    let ang = [s, h, p, zns, ps];
    let (mut a5, mut a6, mut a7, mut a8) = (0.0, 0.0, 0.0, 0.0);
    // sn @ DIU[:, 5]: a sequential dot over the 31 terms
    let mut sn = [0.0; 31];
    let mut cs = [0.0; 31];
    for (k, row) in DIU.iter().enumerate() {
        let mut th = tau;
        for i in 0..5 {
            th += ang[i] * row[i];
        }
        let th = th.to_radians() + zla;
        sn[k] = th.sin();
        cs[k] = th.cos();
    }
    for k in 0..31 {
        a5 += sn[k] * DIU[k][5];
        a6 += cs[k] * DIU[k][6];
        a7 += sn[k] * DIU[k][7];
        a8 += cs[k] * DIU[k][8];
    }
    let (mut b7, mut b8) = (0.0, 0.0);
    for k in 0..31 {
        b7 += cs[k] * DIU[k][7];
        b8 += sn[k] * DIU[k][8];
    }
    let dr = (2.0 * sinphi * cosphi) * (a5 + a6);
    let dn = (cosphi * cosphi - sinphi * sinphi) * (a7 + a8);
    let de = sinphi * (b7 - b8);
    let mut out = neu_to_xyz(dr, dn, de, sinphi, cosphi, sinla, cosla);
    let mut cl = [0.0; 5];
    let mut sl = [0.0; 5];
    for (k, row) in LON.iter().enumerate() {
        let mut th = 0.0;
        for i in 0..5 {
            th += ang[i] * row[i];
        }
        let th = th.to_radians();
        cl[k] = th.cos();
        sl[k] = th.sin();
    }
    let (mut c5, mut c7, mut c6, mut c8) = (0.0, 0.0, 0.0, 0.0);
    for k in 0..5 {
        c5 += cl[k] * LON[k][5];
        c7 += sl[k] * LON[k][7];
        c6 += cl[k] * LON[k][6];
        c8 += sl[k] * LON[k][8];
    }
    let p2 = (3.0 * sinphi * sinphi - 1.0) / 2.0;
    let dr = p2 * (c5 + c7);
    let dn = (2.0 * cosphi * sinphi) * (c6 + c8);
    let v = neu_to_xyz(dr, dn, 0.0, sinphi, cosphi, sinla, cosla);
    for i in 0..3 {
        out[i] = (out[i] + v[i]) / 1000.0;
    }
    out
}

/// IERS 2010 solid-earth tide (conventional tide-free) for one epoch given Sun and Moon (ECEF).
#[must_use]
pub fn solid_tide_bodies(
    xsta: [f64; 3],
    xsun: [f64; 3],
    xmon: [f64; 3],
    t: f64,
    fhr: f64,
) -> [f64; 3] {
    let a = step1(xsta, xsun, xmon);
    let b = step2(xsta, t, fhr);
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// IERS 2010 solid-earth tide at GPS seconds `t` (the leap second of the first epoch for all).
#[must_use]
pub fn solid_tide(xsta: [f64; 3], t: &[f64]) -> Vec<[f64; 3]> {
    let Some(&t0) = t.first() else {
        return Vec::new();
    };
    let leap = leap_at(t0);
    t.iter()
        .map(|&tt| {
            let (sun, moon) = sun_moon_ecef(tt);
            solid_tide_bodies(xsta, sun, moon, tt_centuries(tt), ut_hours(tt, leap))
        })
        .collect()
}

/// Pole tide displacement (m, ECEF) for a session (its mean time `t_mean`), IERS 2010 eq. 7.26 with the
/// 2018 linear mean pole.
#[must_use]
pub fn pole_tide(xsta: [f64; 3], t_mean: f64, xp_arcsec: f64, yp_arcsec: f64) -> [f64; 3] {
    let yr = 2000.0 + (t_mean - 630_763_213.0) / (365.25 * 86400.0);
    let xbar = (55.0 + 1.677 * (yr - 2000.0)) / 1000.0;
    let ybar = (320.5 + 3.460 * (yr - 2000.0)) / 1000.0;
    let (m1, m2) = (xp_arcsec - xbar, -(yp_arcsec - ybar));
    let (lat, lon, _) = ecef_to_llh(xsta);
    let th = std::f64::consts::FRAC_PI_2 - lat;
    let lam = lon;
    let sr = -0.033 * (2.0 * th).sin() * (m1 * lam.cos() + m2 * lam.sin());
    let sth = -0.009 * (2.0 * th).cos() * (m1 * lam.cos() + m2 * lam.sin());
    let slam = 0.009 * th.cos() * (m1 * lam.sin() - m2 * lam.cos());
    let r = enu_matrix(lat, lam);
    crate::linalg::mat3t_vec(&r, [slam, -sth, sr])
}

/// Constituents in BLQ order (M2 S2 N2 K2 K1 O1 P1 Q1 Mf Mm Ssa): speed (rad/s), multipliers of h0, s0,
/// p0, phase (cycles).
const OL_ARGS: [(f64, f64, f64, f64, f64); 11] = [
    (1.40519e-4, 2.0, -2.0, 0.0, 0.0),
    (1.45444e-4, 0.0, 0.0, 0.0, 0.0),
    (1.37880e-4, 2.0, -3.0, 1.0, 0.0),
    (1.45842e-4, 2.0, 0.0, 0.0, 0.0),
    (0.72921e-4, 1.0, 0.0, 0.0, 0.25),
    (0.67598e-4, 1.0, -2.0, 0.0, -0.25),
    (0.72523e-4, -1.0, 0.0, 0.0, -0.25),
    (0.64959e-4, 1.0, -3.0, 1.0, -0.25),
    (0.53234e-5, 0.0, 2.0, 0.0, 0.0),
    (0.26392e-5, 0.0, 1.0, -1.0, 0.0),
    (0.03982e-5, 2.0, 0.0, 0.0, 0.0),
];

/// Ocean tide loading displacement (m, ECEF) at each time from BLQ coefficients (6 rows: amplitudes radial,
/// west, south; phases radial, west, south).
#[must_use]
pub fn ocean_loading(blq: &[[f64; 11]; 6], xsta: [f64; 3], t: &[f64]) -> Vec<[f64; 3]> {
    let Some(&t0) = t.first() else {
        return Vec::new();
    };
    let leap = f64::from(leap_at(t0));
    let (lat, lon, _) = ecef_to_llh(xsta);
    let r = enu_matrix(lat, lon);
    t.iter()
        .map(|&tt| {
            let utc = tt - leap;
            let fday = utc.rem_euclid(86400.0);
            let day0 = utc - fday;
            let days = (day0 + 1831.0 * 86400.0) / 86400.0 + 1.0;
            let tc = (27392.500_528 + 1.000_000_035 * days) / 36525.0;
            let h0 = (279.69668 + 36000.768_930_485 * tc + 3.03e-4 * tc * tc).to_radians();
            let s0 = (270.434_358 + 481_267.883_141_37 * tc - 0.001_133 * tc * tc
                + 1.9e-6 * tc * tc * tc)
                .to_radians();
            let p0 = (334.329_653 + 4069.034_032_957_7 * tc
                - 0.010_325 * tc * tc
                - 1.2e-5 * tc * tc * tc)
                .to_radians();
            let (mut du, mut dw, mut ds) = (0.0, 0.0, 0.0);
            for (k, &(speed, mh, ms, mp, cyc)) in OL_ARGS.iter().enumerate() {
                let arg =
                    speed * fday + mh * h0 + ms * s0 + mp * p0 + 2.0 * std::f64::consts::PI * cyc;
                du += blq[0][k] * (arg - blq[3][k].to_radians()).cos();
                dw += blq[1][k] * (arg - blq[4][k].to_radians()).cos();
                ds += blq[2][k] * (arg - blq[5][k].to_radians()).cos();
            }
            let enu = [-dw, -ds, du];
            // enu @ R: out[i] = sum_k enu[k] * R[k][i]
            [
                enu[0] * r[0][0] + enu[1] * r[1][0] + enu[2] * r[2][0],
                enu[0] * r[0][1] + enu[1] * r[1][1] + enu[2] * r[2][1],
                enu[0] * r[0][2] + enu[1] * r[1][2] + enu[2] * r[2][2],
            ]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tide_magnitudes() {
        let x = [4_395_470.63, 617_743.11, 4_565_607.21];
        let t = [1_411_084_800.0, 1_411_088_400.0];
        let s = solid_tide(x, &t);
        assert!(norm3(s[0]) < 0.5 && norm3(s[0]) > 0.0);
        let d2 = crate::models::solid_tide_degree2(x, &t);
        let diff = norm3(crate::linalg::sub3(s[0], d2[0]));
        assert!(diff < 0.02, "full vs degree-2 differ by {diff} m");
        let p = pole_tide(x, t[0], 0.15, 0.35);
        assert!(norm3(p) < 0.05);
        let blq = [
            [0.001; 11],
            [0.001; 11],
            [0.001; 11],
            [0.0; 11],
            [0.0; 11],
            [0.0; 11],
        ];
        let o = ocean_loading(&blq, x, &t);
        assert!(norm3(o[0]) < 0.02);
        assert!((ut_hours(t[0], 18) - (t[0] - 18.0).rem_euclid(86400.0) / 3600.0).abs() < 1e-12);
    }
}
