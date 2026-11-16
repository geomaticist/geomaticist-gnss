//! Satellite attitude, phase wind-up, nadir angle, Earth shadow, and the relativistic clock and Shapiro
//! terms. Wind-up follows RTKLIB's `model_phw` arrangement (T. Takasu,
//! BSD-2-Clause); the rest is written from Wu et al. (1993) and IERS 2010.

use geomaticist_gnss_core::consts::{CLIGHT, RE_WGS84};

use crate::linalg::{cross3, dot3, norm3, sub3};

/// GM of the Earth the Shapiro term uses.
pub const GM_EARTH_PPP: f64 = 3.986_004_418e14;

/// Nominal attitude body axes (x, y, z): z toward the Earth's centre, y = z × (Sun − satellite),
/// x = y × z.
#[must_use]
pub fn sat_axes(sat: [f64; 3], sun: [f64; 3]) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let n = norm3(sat);
    let ez = [-sat[0] / n, -sat[1] / n, -sat[2] / n];
    let es = sub3(sun, sat);
    let ns = norm3(es);
    let es = [es[0] / ns, es[1] / ns, es[2] / ns];
    let ey = cross3(ez, es);
    let ny = norm3(ey);
    let ey = [ey[0] / ny, ey[1] / ny, ey[2] / ny];
    let ex = cross3(ey, ez);
    (ex, ey, ez)
}

/// Phase wind-up (cycles) along one satellite's arc, continuous (unwrapped) from `prev` (the value at the
/// epoch before, or `None`). `r_enu` is the station's ENU matrix (rows E, N, U).
#[must_use]
pub fn windup(
    sat: &[[f64; 3]],
    rec: [f64; 3],
    sun: &[[f64; 3]],
    r_enu: &[[f64; 3]; 3],
    prev: Option<f64>,
) -> Vec<f64> {
    let exr = r_enu[1];
    let eyr = [-r_enu[0][0], -r_enu[0][1], -r_enu[0][2]];
    let mut out = Vec::with_capacity(sat.len());
    let mut last = prev;
    for (s, su) in sat.iter().zip(sun) {
        let (exs, eys, _) = sat_axes(*s, *su);
        let k = sub3(rec, *s);
        let nk = norm3(k);
        let k = [k[0] / nk, k[1] / nk, k[2] / nk];
        let eks = cross3(k, eys);
        let ekr = cross3(k, eyr);
        let kxs = dot3(k, exs);
        let kxr = dot3(k, exr);
        let mut ds = [0.0; 3];
        let mut dr = [0.0; 3];
        for i in 0..3 {
            ds[i] = exs[i] - k[i] * kxs - eks[i];
            dr[i] = exr[i] - k[i] * kxr + ekr[i];
        }
        let c = dot3(ds, dr) / (norm3(ds) * norm3(dr));
        let mut ph = c.clamp(-1.0, 1.0).acos() / (2.0 * std::f64::consts::PI);
        let sgn = dot3(k, cross3(ds, dr));
        if sgn < 0.0 {
            ph = -ph;
        }
        if ph.is_finite() {
            let v = match last {
                None => ph,
                Some(l) => ph + (l - ph + 0.5).floor(),
            };
            last = Some(v);
            out.push(v);
        } else {
            out.push(f64::NAN);
        }
    }
    out
}

/// Nadir angle (deg) of the receiver seen from the satellite.
#[must_use]
pub fn nadir_angle(sat: [f64; 3], rec: [f64; 3]) -> f64 {
    let n = norm3(sat);
    let ez = [-sat[0] / n, -sat[1] / n, -sat[2] / n];
    let k = sub3(rec, sat);
    let nk = norm3(k);
    let k = [k[0] / nk, k[1] / nk, k[2] / nk];
    dot3(ez, k).clamp(-1.0, 1.0).acos().to_degrees()
}

/// True where the satellite is in the Earth's (cylindrical) shadow.
#[must_use]
pub fn in_shadow(sat: [f64; 3], sun: [f64; 3]) -> bool {
    let n = norm3(sun);
    let su = [sun[0] / n, sun[1] / n, sun[2] / n];
    let along = dot3(sat, su);
    let perp = norm3([
        sat[0] - along * su[0],
        sat[1] - along * su[1],
        sat[2] - along * su[2],
    ]);
    along < 0.0 && perp < RE_WGS84
}

/// Periodic relativistic clock correction (s) for a precise satellite clock: `−2 r·v / c²`.
#[must_use]
pub fn rel_clock(sat: [f64; 3], vel: [f64; 3]) -> f64 {
    -2.0 * dot3(sat, vel) / (CLIGHT * CLIGHT)
}

/// Shapiro (gravitational) path delay (m).
#[must_use]
pub fn shapiro(sat: [f64; 3], rec: [f64; 3]) -> f64 {
    let rs = norm3(sat);
    let rr = norm3(rec);
    let rho = norm3(sub3(sat, rec));
    2.0 * GM_EARTH_PPP / (CLIGHT * CLIGHT) * ((rs + rr + rho) / (rs + rr - rho)).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axes_are_orthonormal() {
        let sat = [15600e3, 7540e3, 20140e3];
        let sun = [1.0e11, 5.0e10, -3.0e10];
        let (ex, ey, ez) = sat_axes(sat, sun);
        assert!(dot3(ex, ey).abs() < 1e-12 && dot3(ey, ez).abs() < 1e-12);
        assert!((norm3(ex) - 1.0).abs() < 1e-12);
        let rec = [4_395_470.63, 617_743.11, 4_565_607.21];
        assert!(nadir_angle(sat, rec) < 15.0);
        assert!(!in_shadow(sat, sun));
        assert!(shapiro(sat, rec) > 0.0 && shapiro(sat, rec) < 0.03);
        let w = windup(
            &[sat, sat],
            rec,
            &[sun, sun],
            &crate::models::enu_matrix(1.1, -2.5),
            None,
        );
        assert_eq!(w[0], w[1]);
    }
}
