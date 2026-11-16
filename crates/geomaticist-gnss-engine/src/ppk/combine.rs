//! The forward and backward passes merged epoch by epoch.
//!
//! The two passes see the same observations at an epoch but reach it with different histories (the data
//! before it, the data after it), so their agreement is a real check. The rules:
//!
//! - only one pass has a position: it is taken;
//! - both fixed and they agree: fixed, the mean of the two weighted by their covariances (a pass may
//!   hold fewer satellites at an epoch). Their errors are largely shared (same phases, same integers),
//!   so the covariance is the smaller of the two, not halved;
//! - both fixed but they disagree: one of them has wrong integers and nothing says which. The epoch is
//!   demoted to float, its σ widened to the disagreement, and it is counted in
//!   [`super::PpkSummary::fb_conflicts`];
//! - one fixed, the other float: the fixed one, provided it lies inside the float one's uncertainty;
//!   otherwise the epoch is float (and counted as a conflict);
//! - both float: the covariance-weighted mean (the passes' float errors are largely independent);
//! - otherwise the better of the two.
//!
//! Taking the better status and checking that two fixes agree is what RTKLIB's combined solution does
//! (T. Takasu, BSD-2-Clause; `postpos.c` `combres`); the float cross-check of a one-sided fix and the
//! treatment of the fixed covariance are this engine's.

use super::options::PpkOptions;
use super::{EpochSolution, FixStatus, PpkSummary, enu_at};
use crate::linalg::{inv3, norm3};

/// Two fixes agree when they differ by less than this many combined standard deviations ...
pub const AGREE_SIGMA: f64 = 4.0;
/// ... or by less than this many metres in any case (east or north; twice that in height).
pub const AGREE_FLOOR_M: f64 = 0.03;

fn agree(a: &EpochSolution, b: &EpochSolution, d: [f64; 3]) -> bool {
    (0..3).all(|i| {
        let s = (a.cov_enu[i][i].max(0.0) + b.cov_enu[i][i].max(0.0)).sqrt();
        let floor = if i == 2 {
            2.0 * AGREE_FLOOR_M
        } else {
            AGREE_FLOOR_M
        };
        d[i].abs() <= (AGREE_SIGMA * s).max(floor)
    })
}

fn add_cov(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

/// Covariance-weighted mean of two float solutions; `None` when a covariance cannot be inverted.
fn weighted(
    f: &EpochSolution,
    b: &EpochSolution,
    d_enu: [f64; 3],
) -> Option<([f64; 3], [[f64; 3]; 3])> {
    let wf = inv3(&f.cov_enu)?;
    let wb = inv3(&b.cov_enu)?;
    let c = inv3(&add_cov(&wf, &wb))?;
    // x = xf + C Wb (xb - xf), in ENU: the gain applied to the difference
    let mut g = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                g[i] += c[i][j] * wb[j][k] * d_enu[k];
            }
        }
    }
    Some((g, c))
}

/// Merges the passes. `summary` receives the agreement statistics.
#[must_use]
pub fn combine(
    fwd: &[EpochSolution],
    bwd: &[EpochSolution],
    _opts: &PpkOptions,
    summary: &mut PpkSummary,
) -> Vec<EpochSolution> {
    let n = fwd.len().min(bwd.len());
    let mut out = Vec::with_capacity(n);
    let mut ss = [0.0; 3];
    for k in 0..n {
        let (f, b) = (&fwd[k], &bwd[k]);
        if !f.has_position() || !b.has_position() {
            out.push(if f.has_position() { *f } else { *b });
            continue;
        }
        // backward minus forward, east / north / up at the forward position
        let d = enu_at(f.xyz, b.xyz);
        let enu_to_xyz = |v: [f64; 3]| -> [f64; 3] {
            let (lat, lon, _) = crate::models::ecef_to_llh(f.xyz);
            let r = crate::models::enu_matrix(lat, lon);
            let dv = crate::linalg::mat3t_vec(&r, v);
            [f.xyz[0] + dv[0], f.xyz[1] + dv[1], f.xyz[2] + dv[2]]
        };
        let mut s = *f;
        s.n_sats = f.n_sats.max(b.n_sats);
        s.ratio = f.ratio.min(b.ratio);
        match (f.status, b.status) {
            (FixStatus::Fixed, FixStatus::Fixed) => {
                summary.both_fixed += 1;
                for i in 0..3 {
                    ss[i] += d[i] * d[i];
                }
                summary.fb_max = summary.fb_max.max(norm3(d));
                if agree(f, b, d) {
                    // weighted by their covariances: a pass that holds fewer satellites at this epoch
                    // counts less; the covariance is the better one's, not smaller
                    let g = weighted(f, b, d).map_or([0.5 * d[0], 0.5 * d[1], 0.5 * d[2]], |x| x.0);
                    s.set_xyz(enu_to_xyz(g));
                    let tf = f.cov_enu[0][0] + f.cov_enu[1][1] + f.cov_enu[2][2];
                    let tb = b.cov_enu[0][0] + b.cov_enu[1][1] + b.cov_enu[2][2];
                    s.cov_enu = if tf <= tb { f.cov_enu } else { b.cov_enu };
                } else {
                    s.set_xyz(enu_to_xyz([0.5 * d[0], 0.5 * d[1], 0.5 * d[2]]));
                    summary.fb_conflicts += 1;
                    s.status = FixStatus::Float;
                    s.cov_enu = add_cov(&f.cov_enu, &b.cov_enu);
                    for i in 0..3 {
                        s.cov_enu[i][i] += 0.25 * d[i] * d[i];
                    }
                }
            }
            (FixStatus::Fixed, FixStatus::Float) | (FixStatus::Float, FixStatus::Fixed) => {
                let (fx, fl) = if f.status == FixStatus::Fixed {
                    (f, b)
                } else {
                    (b, f)
                };
                let inside = (0..3).all(|i| {
                    d[i].abs()
                        <= AGREE_SIGMA * fl.cov_enu[i][i].max(0.0).sqrt() + 2.0 * AGREE_FLOOR_M
                });
                if inside {
                    s = *fx;
                    s.ratio = fx.ratio;
                } else {
                    summary.fb_conflicts += 1;
                    s = *fl;
                    for i in 0..3 {
                        s.cov_enu[i][i] += 0.25 * d[i] * d[i];
                    }
                }
            }
            (FixStatus::Float, FixStatus::Float)
            | (FixStatus::Single, FixStatus::Single)
            | (FixStatus::Ppp, FixStatus::Ppp) => {
                if let Some((g, c)) = weighted(f, b, d) {
                    s.set_xyz(enu_to_xyz(g));
                    s.cov_enu = c;
                    // two float solutions that disagree beyond their σ: say so in the σ
                    if !agree(f, b, d) {
                        for i in 0..3 {
                            s.cov_enu[i][i] += 0.25 * d[i] * d[i];
                        }
                    }
                }
            }
            _ => {
                s = if f.status.better_than(b.status) {
                    *f
                } else {
                    *b
                };
            }
        }
        out.push(s);
    }
    if summary.both_fixed > 0 {
        let nb = summary.both_fixed as f64;
        summary.fb_rms_enu = [
            (ss[0] / nb).sqrt(),
            (ss[1] / nb).sqrt(),
            (ss[2] / nb).sqrt(),
        ];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ppk::llh_deg;

    fn sol(status: FixStatus, east: f64, sig: f64) -> EpochSolution {
        let base = [4_395_470.6, 617_743.1, 4_565_607.2];
        let (lat, lon, _) = crate::models::ecef_to_llh(base);
        let r = crate::models::enu_matrix(lat, lon);
        let d = crate::linalg::mat3t_vec(&r, [east, 0.0, 0.0]);
        let xyz = [base[0] + d[0], base[1] + d[1], base[2] + d[2]];
        let mut s = EpochSolution::none(100.0);
        s.status = status;
        s.xyz = xyz;
        s.llh = llh_deg(xyz);
        s.cov_enu = [
            [sig * sig, 0.0, 0.0],
            [0.0, sig * sig, 0.0],
            [0.0, 0.0, sig * sig],
        ];
        s
    }

    #[test]
    fn rules() {
        let o = PpkOptions::default();
        let mut sm = PpkSummary::default();
        // both fixed, 4 mm apart: fixed, the mean
        let c = combine(
            &[sol(FixStatus::Fixed, 0.0, 0.01)],
            &[sol(FixStatus::Fixed, 0.004, 0.01)],
            &o,
            &mut sm,
        );
        assert_eq!(c[0].status, FixStatus::Fixed);
        let d = enu_at(sol(FixStatus::Fixed, 0.0, 0.01).xyz, c[0].xyz);
        assert!((d[0] - 0.002).abs() < 1e-6, "{d:?}");
        assert!((c[0].cov_enu[0][0] - 1e-4).abs() < 1e-12);
        // both fixed, 20 cm apart: a conflict, float
        let c = combine(
            &[sol(FixStatus::Fixed, 0.0, 0.01)],
            &[sol(FixStatus::Fixed, 0.2, 0.01)],
            &o,
            &mut sm,
        );
        assert_eq!(c[0].status, FixStatus::Float);
        assert_eq!(sm.fb_conflicts, 1);
        assert!(c[0].cov_enu[0][0].sqrt() >= 0.1);
        // fixed beside a float that contains it: fixed; beside one that does not: float
        let c = combine(
            &[sol(FixStatus::Fixed, 0.0, 0.01)],
            &[sol(FixStatus::Float, 0.3, 0.2)],
            &o,
            &mut sm,
        );
        assert_eq!(c[0].status, FixStatus::Fixed);
        let c = combine(
            &[sol(FixStatus::Float, 0.0, 0.02)],
            &[sol(FixStatus::Fixed, 0.5, 0.01)],
            &o,
            &mut sm,
        );
        assert_eq!(c[0].status, FixStatus::Float);
        // two floats: the weighted mean, the tighter one counts more
        let c = combine(
            &[sol(FixStatus::Float, 0.0, 0.1)],
            &[sol(FixStatus::Float, 0.1, 0.2)],
            &o,
            &mut sm,
        );
        let d = enu_at(sol(FixStatus::Float, 0.0, 0.1).xyz, c[0].xyz);
        assert!((d[0] - 0.02).abs() < 1e-6, "{d:?}");
        // one side only
        let c = combine(
            &[EpochSolution::none(1.0)],
            &[sol(FixStatus::Single, 0.0, 3.0)],
            &o,
            &mut sm,
        );
        assert_eq!(c[0].status, FixStatus::Single);
    }
}
