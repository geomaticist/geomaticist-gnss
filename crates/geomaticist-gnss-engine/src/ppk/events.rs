//! Positions at arbitrary times (camera exposures, marker presses) by interpolation between epochs, and
//! averages over stop-and-go occupations.
//!
//! An event lies between two epochs of the track. Two interpolations are offered:
//!
//! - **2-point (linear)**: the two epochs around the event;
//! - **4-point (cubic)**: the two epochs before and the two after, a cubic through the four positions
//!   (Lagrange form, so uneven spacing is handled).
//!
//! The times used are the epochs' receive times ([`EpochSolution::t`]).
//!
//! Rules that protect the result:
//!
//! - an event before the first or after the last epoch with a position gets **no position**: nothing is
//!   ever extrapolated;
//! - when epochs are missing next to the event, the event gets no position unless the number missing is
//!   within the limit the caller sets (`max_missing`); such an event has status [`EventStatus::Gap`], and its
//!   σ grows with the longer span;
//! - the cubic needs four epochs with no more than the allowed gaps between them; where it cannot be
//!   formed (the first and last interval of the data, beside a gap) the linear form is used and
//!   [`EventResult::method`] says so.
//!
//! σ of an event: the epochs' σ combined as fully correlated (neighbouring epochs share their errors, so
//! nothing is gained by averaging them), plus an allowance for the interpolation itself estimated from the
//! curvature of the track around the event.

use geomaticist_gnss_core::event::Event;
use serde::{Deserialize, Serialize};

use super::{EpochSolution, FixStatus, enu_at, llh_deg};
use crate::linalg::mat3t_vec;
use crate::models::{ecef_to_llh, enu_matrix};

/// How a position at an event is interpolated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum EventInterp {
    /// The two epochs around the event, linear.
    Linear2,
    /// Two epochs before and two after, cubic.
    #[default]
    Cubic4,
    /// The two epochs around the event with their Doppler velocities, cubic (Hermite form). Worth it at
    /// low rates (1 Hz), where the receiver's velocity knows what the track does between epochs; at 5 Hz
    /// and above the 4-point cubic is as good or better. Needs Doppler in the rover file; without a
    /// velocity at both epochs the linear form is used.
    Hermite2,
}

impl EventInterp {
    /// `"linear2"` or `"cubic4"`.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            EventInterp::Linear2 => "linear2",
            EventInterp::Cubic4 => "cubic4",
            EventInterp::Hermite2 => "hermite2",
        }
    }
}

/// Noise assumed for a Doppler velocity in the σ of a Hermite result, m/s ...
const HERMITE_VEL_SIGMA: f64 = 0.05;
/// ... plus this much per m/s² of acceleration between the two epochs (the receiver's velocity lags).
const HERMITE_VEL_LAG_S: f64 = 0.03;

/// Whether an event got a position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum EventStatus {
    /// Interpolated between neighbouring epochs.
    Ok,
    /// Interpolated across missing epochs (within the limit).
    Gap,
    /// No position; [`EventResult::why`] says why.
    NoSolution,
}

/// The position at one event.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EventResult {
    /// Position of the event in the input list.
    pub index: usize,
    /// The event's label.
    pub label: String,
    /// The event's time, GPS seconds.
    pub t: f64,
    /// Whether there is a position.
    pub status: EventStatus,
    /// Why there is none (empty otherwise).
    pub why: String,
    /// The worst status of the epochs used.
    pub fix: FixStatus,
    /// The point wanted: the antenna position plus the event's offset (ECEF m).
    pub xyz: [f64; 3],
    /// Its latitude, longitude (degrees) and ellipsoidal height (m).
    pub llh: [f64; 3],
    /// The antenna position at the event (ECEF m).
    pub antenna_xyz: [f64; 3],
    /// Standard deviations east, north, up (m).
    pub sigma_enu: [f64; 3],
    /// Indices of the epochs used, in time order.
    pub epochs_used: Vec<usize>,
    /// Their statuses.
    pub epochs_fix: Vec<FixStatus>,
    /// The interpolation actually used.
    pub method: EventInterp,
    /// Seconds between the two epochs around the event.
    pub gap_s: f64,
    /// Epochs missing between them.
    pub missing: usize,
}

/// The median spacing of sorted times (0 for fewer than two).
#[must_use]
pub fn median_interval(t: &[f64]) -> f64 {
    if t.len() < 2 {
        return 0.0;
    }
    let mut d: Vec<f64> = t
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|v| *v > 0.0)
        .collect();
    if d.is_empty() {
        return 0.0;
    }
    d.sort_by(f64::total_cmp);
    d[d.len() / 2]
}

/// Lagrange weights of the nodes `t` at `x`.
#[must_use]
pub fn lagrange_weights(t: &[f64], x: f64) -> Vec<f64> {
    (0..t.len())
        .map(|i| {
            let mut w = 1.0;
            for j in 0..t.len() {
                if j != i {
                    w *= (x - t[j]) / (t[i] - t[j]);
                }
            }
            w
        })
        .collect()
}

/// A cubic through two positions with their velocities (Hermite form), at `x` between `t0` and `t1`.
/// The interpolation behind [`EventInterp::Hermite2`].
#[must_use]
pub fn hermite(
    t0: f64,
    p0: [f64; 3],
    v0: [f64; 3],
    t1: f64,
    p1: [f64; 3],
    v1: [f64; 3],
    x: f64,
) -> [f64; 3] {
    let h = t1 - t0;
    let s = (x - t0) / h;
    let (s2, s3) = (s * s, s * s * s);
    let (h00, h10, h01, h11) = (
        2.0 * s3 - 3.0 * s2 + 1.0,
        s3 - 2.0 * s2 + s,
        -2.0 * s3 + 3.0 * s2,
        s3 - s2,
    );
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = h00 * p0[i] + h10 * h * v0[i] + h01 * p1[i] + h11 * h * v1[i];
    }
    out
}

fn no_solution(index: usize, ev: &Event, why: &str, method: EventInterp) -> EventResult {
    EventResult {
        index,
        label: ev.label.clone(),
        t: ev.t,
        status: EventStatus::NoSolution,
        why: why.to_owned(),
        fix: FixStatus::None,
        xyz: [f64::NAN; 3],
        llh: [f64::NAN; 3],
        antenna_xyz: [f64::NAN; 3],
        sigma_enu: [f64::NAN; 3],
        epochs_used: Vec::new(),
        epochs_fix: Vec::new(),
        method,
        gap_s: f64::NAN,
        missing: 0,
    }
}

/// Positions of `events` on `track`.
///
/// `max_missing`: how many missing epochs between the two epochs around an event are tolerated.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn interpolate(
    track: &[EpochSolution],
    events: &[Event],
    method: EventInterp,
    max_missing: usize,
) -> Vec<EventResult> {
    // the epochs with a position, in time order
    let valid: Vec<usize> = (0..track.len())
        .filter(|&i| track[i].has_position())
        .collect();
    let tv: Vec<f64> = valid.iter().map(|&i| track[i].t).collect();
    let h = median_interval(&track.iter().map(|e| e.t_tag).collect::<Vec<_>>());
    let missing_in = |span: f64| -> usize {
        if h > 0.0 {
            ((span / h).round() as usize).saturating_sub(1)
        } else {
            0
        }
    };
    let mut out = Vec::with_capacity(events.len());
    for (index, ev) in events.iter().enumerate() {
        let te = ev.t;
        if !te.is_finite() {
            out.push(no_solution(index, ev, "the event has no time", method));
            continue;
        }
        if tv.len() < 2 {
            out.push(no_solution(
                index,
                ev,
                "no positions to interpolate between",
                method,
            ));
            continue;
        }
        const EDGE: f64 = 1e-6;
        if te < tv[0] - EDGE {
            out.push(no_solution(
                index,
                ev,
                "before the first epoch with a position",
                method,
            ));
            continue;
        }
        if te > tv[tv.len() - 1] + EDGE {
            out.push(no_solution(
                index,
                ev,
                "after the last epoch with a position",
                method,
            ));
            continue;
        }
        // k: tv[k] <= te <= tv[k + 1]
        let k = tv.partition_point(|x| *x <= te).clamp(1, tv.len() - 1) - 1;
        let gap = tv[k + 1] - tv[k];
        let missing = missing_in(gap);
        if missing > max_missing {
            let mut r = no_solution(
                index,
                ev,
                &format!(
                    "{missing} epoch{} without a position beside the event (limit {max_missing})",
                    if missing == 1 { "" } else { "s" }
                ),
                method,
            );
            r.gap_s = gap;
            r.missing = missing;
            out.push(r);
            continue;
        }
        // the cubic needs a usable epoch on each outer side
        let outer_ok = k >= 1
            && k + 2 < tv.len()
            && missing_in(tv[k] - tv[k - 1]) <= max_missing
            && missing_in(tv[k + 2] - tv[k + 1]) <= max_missing;
        let has_vel = |q: usize| track[valid[q]].vel_enu.iter().all(|v| v.is_finite());
        let used_method = match method {
            EventInterp::Cubic4 if outer_ok => EventInterp::Cubic4,
            EventInterp::Hermite2 if has_vel(k) && has_vel(k + 1) => EventInterp::Hermite2,
            _ => EventInterp::Linear2,
        };
        let nodes: Vec<usize> = match used_method {
            EventInterp::Cubic4 => vec![k - 1, k, k + 1, k + 2],
            EventInterp::Linear2 | EventInterp::Hermite2 => vec![k, k + 1],
        };
        let origin = track[valid[k]].xyz;
        let (lat, lon, _) = ecef_to_llh(origin);
        let r = enu_matrix(lat, lon);
        let tn: Vec<f64> = nodes.iter().map(|&q| tv[q]).collect();
        let pn: Vec<[f64; 3]> = nodes
            .iter()
            .map(|&q| enu_at(origin, track[valid[q]].xyz))
            .collect();
        let w = lagrange_weights(&tn, te);
        let mut p = [0.0; 3];
        let mut sig = [0.0; 3];
        for (q, &node) in nodes.iter().enumerate() {
            let s = track[valid[node]].sigma_enu();
            for i in 0..3 {
                p[i] += w[q] * pn[q][i];
                sig[i] += w[q].abs() * s[i];
            }
        }
        // allowance for the interpolation itself
        let mut interp = [0.0; 3];
        match used_method {
            EventInterp::Hermite2 => {
                // each velocity is in the east / north / up of its own epoch: bring both to this frame
                let v = |q: usize| -> [f64; 3] {
                    let e = &track[valid[q]];
                    let (la, lo, _) = ecef_to_llh(e.xyz);
                    let own = [
                        f64::from(e.vel_enu[0]),
                        f64::from(e.vel_enu[1]),
                        f64::from(e.vel_enu[2]),
                    ];
                    crate::linalg::mat3_vec(&r, mat3t_vec(&enu_matrix(la, lo), own))
                };
                let (v0, v1) = (v(k), v(k + 1));
                p = hermite(tn[0], pn[0], v0, tn[1], pn[1], v1, te);
                // the velocities' share of the result: gap × (|h10| + |h11|) × their noise
                let s = (te - tn[0]) / gap;
                let lever = gap * ((s * s * s - 2.0 * s * s + s).abs() + (s * s * s - s * s).abs());
                for i in 0..3 {
                    let acc = (v1[i] - v0[i]).abs() / gap;
                    interp[i] = lever * (HERMITE_VEL_SIGMA + HERMITE_VEL_LAG_S * acc);
                }
            }
            EventInterp::Cubic4 => {
                // the two quadratics through three of the four nodes differ by the third derivative
                let w1 = lagrange_weights(&tn[..3], te);
                let w2 = lagrange_weights(&tn[1..], te);
                for i in 0..3 {
                    let q1: f64 = (0..3).map(|q| w1[q] * pn[q][i]).sum();
                    let q2: f64 = (0..3).map(|q| w2[q] * pn[q + 1][i]).sum();
                    interp[i] = 0.25 * (q1 - q2).abs();
                }
            }
            EventInterp::Linear2 => {
                // |a| τ (g - τ) / 2 with the acceleration from a neighbouring epoch, 1 m/s² when there is none
                let tau = te - tv[k];
                let third = if k >= 1 && missing_in(tv[k] - tv[k - 1]) <= max_missing {
                    Some(k - 1)
                } else if k + 2 < tv.len() && missing_in(tv[k + 2] - tv[k + 1]) <= max_missing {
                    Some(k + 2)
                } else {
                    None
                };
                for i in 0..3 {
                    let acc = third.map_or(1.0, |q| {
                        let (ta, tb, tc) = if q < k {
                            (tv[q], tv[k], tv[k + 1])
                        } else {
                            (tv[k], tv[k + 1], tv[q])
                        };
                        let e = |j: usize| enu_at(origin, track[valid[j]].xyz)[i];
                        let (pa, pb, pc) = if q < k {
                            (e(q), e(k), e(k + 1))
                        } else {
                            (e(k), e(k + 1), e(q))
                        };
                        // second divided difference × 2
                        2.0 * ((pc - pb) / (tc - tb) - (pb - pa) / (tb - ta)) / (tc - ta)
                    });
                    interp[i] = 0.5 * acc.abs() * tau * (gap - tau);
                }
            }
        }
        for i in 0..3 {
            sig[i] = sig[i].hypot(interp[i]);
        }
        let d = mat3t_vec(&r, p);
        let antenna = [origin[0] + d[0], origin[1] + d[1], origin[2] + d[2]];
        let xyz = match ev.offset {
            Some(o) => {
                let (la, lo, _) = ecef_to_llh(antenna);
                let dv = mat3t_vec(&enu_matrix(la, lo), o.enu());
                [antenna[0] + dv[0], antenna[1] + dv[1], antenna[2] + dv[2]]
            }
            None => antenna,
        };
        let epochs_used: Vec<usize> = nodes.iter().map(|&q| valid[q]).collect();
        let epochs_fix: Vec<FixStatus> = epochs_used.iter().map(|&i| track[i].status).collect();
        let fix = epochs_fix.iter().fold(FixStatus::Fixed, |a, b| a.worst(*b));
        out.push(EventResult {
            index,
            label: ev.label.clone(),
            t: te,
            status: if missing > 0 {
                EventStatus::Gap
            } else {
                EventStatus::Ok
            },
            why: String::new(),
            fix,
            xyz,
            llh: llh_deg(xyz),
            antenna_xyz: antenna,
            sigma_enu: sig,
            epochs_used,
            epochs_fix,
            method: used_method,
            gap_s: gap,
            missing,
        });
    }
    out
}

pub use geomaticist_gnss_core::event::Occupation;

/// The averaged position of an occupation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OccupationResult {
    /// The point's name.
    pub name: String,
    /// Start, GPS seconds.
    pub t0: f64,
    /// End, GPS seconds.
    pub t1: f64,
    /// The status of the epochs averaged (the best available in the occupation).
    pub fix: FixStatus,
    /// Epochs averaged.
    pub n: usize,
    /// Epochs of the occupation with any position.
    pub n_all: usize,
    /// Mean position, ECEF m (NaN when no epoch has a position).
    pub xyz: [f64; 3],
    /// Latitude, longitude (degrees), ellipsoidal height (m).
    pub llh: [f64; 3],
    /// Standard deviations east, north, up (m): the epochs' RMS σ (not divided by √n: the epochs of one
    /// occupation share most of their errors).
    pub sigma_enu: [f64; 3],
    /// Scatter of the epochs about the mean, east / north / up (m).
    pub scatter_enu: [f64; 3],
}

/// Averages the track over each occupation: the fixed epochs when there are any, else the float ones,
/// else the code-only ones.
#[must_use]
pub fn average_occupations(
    track: &[EpochSolution],
    occupations: &[Occupation],
) -> Vec<OccupationResult> {
    occupations
        .iter()
        .map(|oc| {
            let inside: Vec<&EpochSolution> = track
                .iter()
                .filter(|e| e.has_position() && e.t_tag >= oc.t0 - 1e-6 && e.t_tag <= oc.t1 + 1e-6)
                .collect();
            let fix = inside.iter().fold(FixStatus::None, |a, e| {
                if e.status.better_than(a) { e.status } else { a }
            });
            let sel: Vec<&&EpochSolution> = inside.iter().filter(|e| e.status == fix).collect();
            let mut res = OccupationResult {
                name: oc.name.clone(),
                t0: oc.t0,
                t1: oc.t1,
                fix,
                n: sel.len(),
                n_all: inside.len(),
                xyz: [f64::NAN; 3],
                llh: [f64::NAN; 3],
                sigma_enu: [f64::NAN; 3],
                scatter_enu: [f64::NAN; 3],
            };
            if sel.is_empty() {
                return res;
            }
            let n = sel.len() as f64;
            let mut m = [0.0; 3];
            for e in &sel {
                for i in 0..3 {
                    m[i] += e.xyz[i] / n;
                }
            }
            let (mut s2, mut sc) = ([0.0; 3], [0.0; 3]);
            for e in &sel {
                let d = enu_at(m, e.xyz);
                for i in 0..3 {
                    s2[i] += e.cov_enu[i][i].max(0.0) / n;
                    sc[i] += d[i] * d[i] / n;
                }
            }
            res.xyz = m;
            res.llh = llh_deg(m);
            res.sigma_enu = [s2[0].sqrt(), s2[1].sqrt(), s2[2].sqrt()];
            res.scatter_enu = [sc[0].sqrt(), sc[1].sqrt(), sc[2].sqrt()];
            res
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use super::*;
    use geomaticist_gnss_core::event::EventOffset;

    const BASE: [f64; 3] = [4_395_470.6, 617_743.1, 4_565_607.2];

    /// A track whose east / north / up (m) at time `t` is `f(t)`, epochs every `h` seconds.
    fn track(f: &dyn Fn(f64) -> [f64; 3], h: f64, n: usize) -> Vec<EpochSolution> {
        let (lat, lon, _) = ecef_to_llh(BASE);
        let r = enu_matrix(lat, lon);
        (0..n)
            .map(|k| {
                let t = 1000.0 + h * k as f64;
                let d = mat3t_vec(&r, f(t));
                let mut s = EpochSolution::none(t);
                s.status = FixStatus::Fixed;
                s.set_xyz([BASE[0] + d[0], BASE[1] + d[1], BASE[2] + d[2]]);
                s.cov_enu = [[1e-4, 0.0, 0.0], [0.0, 1e-4, 0.0], [0.0, 0.0, 4e-4]];
                s
            })
            .collect()
    }

    fn err(f: &dyn Fn(f64) -> [f64; 3], r: &EventResult) -> f64 {
        let p = enu_at(BASE, r.antenna_xyz);
        let q = f(r.t);
        ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
    }

    #[test]
    fn straight_line_is_exact_for_both() {
        let f = |t: f64| [15.0 * (t - 1000.0), -3.0 * (t - 1000.0), 100.0];
        let tr = track(&f, 0.2, 50);
        let ev: Vec<Event> = (0..40)
            .map(|k| Event::at(1000.3 + 0.2137 * f64::from(k)))
            .collect();
        for m in [EventInterp::Linear2, EventInterp::Cubic4] {
            for r in interpolate(&tr, &ev, m, 0) {
                assert_eq!(r.status, EventStatus::Ok);
                assert!(err(&f, &r) < 1e-6, "{m:?} {}", err(&f, &r));
            }
        }
    }

    #[test]
    fn a_turn_favours_the_cubic() {
        // a circle of 30 m radius flown at 12 m/s: 4.8 m/s² of centripetal acceleration
        let w = 12.0 / 30.0;
        let f = move |t: f64| [30.0 * (w * t).cos(), 30.0 * (w * t).sin(), 100.0];
        let ev: Vec<Event> = (0..50)
            .map(|k| Event::at(1001.0 + 0.3371 * f64::from(k)))
            .collect();
        let tr = track(&f, 0.2, 200);
        let max = |m| {
            interpolate(&tr, &ev, m, 0)
                .iter()
                .map(|r| err(&f, r))
                .fold(0.0, f64::max)
        };
        let (lin, cub) = (max(EventInterp::Linear2), max(EventInterp::Cubic4));
        // linear: a h² / 8 = 4.8 × 0.04 / 8 = 24 mm at mid-interval; cubic: well under a millimetre
        assert!(lin > 0.015 && lin < 0.025, "{lin}");
        assert!(cub < 2e-4, "{cub}");
        // the σ of the linear result admits the interpolation error
        let r = &interpolate(&tr, &[Event::at(1010.1)], EventInterp::Linear2, 0)[0];
        assert!(r.sigma_enu[0].hypot(r.sigma_enu[1]) > err(&f, r));
    }

    #[test]
    fn never_extrapolates_and_respects_the_gap_limit() {
        let f = |t: f64| [10.0 * (t - 1000.0), 0.0, 50.0];
        let mut tr = track(&f, 0.2, 50);
        // epoch 20 lost
        tr[20].status = FixStatus::None;
        let ev = [
            Event::at(999.9),
            Event::at(1000.0 + 0.2 * 49.0 + 0.05),
            Event::at(1000.0 + 0.2 * 19.5), // beside the lost epoch
            Event::at(1000.0 + 0.2 * 30.5),
        ];
        let r = interpolate(&tr, &ev, EventInterp::Cubic4, 0);
        assert_eq!(r[0].status, EventStatus::NoSolution);
        assert!(r[0].why.contains("before"));
        assert_eq!(r[1].status, EventStatus::NoSolution);
        assert!(r[1].why.contains("after"));
        assert_eq!(r[2].status, EventStatus::NoSolution);
        assert_eq!(r[2].missing, 1);
        assert_eq!(r[3].status, EventStatus::Ok);
        assert_eq!(r[3].epochs_used, vec![29, 30, 31, 32]);
        // with one missing epoch allowed: flagged, a larger σ, still exact on a straight line
        let g = &interpolate(&tr, &ev, EventInterp::Cubic4, 1)[2];
        assert_eq!(g.status, EventStatus::Gap);
        assert_eq!(g.epochs_used, vec![18, 19, 21, 22]);
        assert!(err(&f, g) < 1e-6);
        let lin_gap = &interpolate(&tr, &ev, EventInterp::Linear2, 1)[2];
        let lin_ok = &interpolate(&tr, &ev, EventInterp::Linear2, 1)[3];
        assert!(lin_gap.gap_s > 0.39 && lin_gap.gap_s < 0.41);
        assert!(lin_gap.sigma_enu[0] >= lin_ok.sigma_enu[0]);
        // first interval: the cubic falls back to linear
        let first = &interpolate(&tr, &[Event::at(1000.1)], EventInterp::Cubic4, 0)[0];
        assert_eq!(first.method, EventInterp::Linear2);
        assert_eq!(first.status, EventStatus::Ok);
    }

    #[test]
    fn offset_and_status() {
        let f = |t: f64| [5.0 * (t - 1000.0), 0.0, 80.0];
        let mut tr = track(&f, 1.0, 10);
        tr[5].status = FixStatus::Float;
        let mut e = Event::at(1004.5);
        e.label = "7".into();
        e.offset = Some(EventOffset {
            north: 0.4,
            east: -0.3,
            down: 0.25,
        });
        let r = &interpolate(&tr, &[e], EventInterp::Linear2, 0)[0];
        assert_eq!(r.fix, FixStatus::Float);
        assert_eq!(r.label, "7");
        let d = enu_at(r.antenna_xyz, r.xyz);
        assert!(
            (d[0] + 0.3).abs() < 1e-6 && (d[1] - 0.4).abs() < 1e-6 && (d[2] + 0.25).abs() < 1e-6,
            "{d:?}"
        );
    }

    /// The table of the user documentation: interpolation error on exact tracks (no measurement noise),
    /// maximum over many events, in millimetres. Run with `--nocapture` to print it.
    #[test]
    fn interpolation_table() {
        type Path = fn(f64) -> ([f64; 3], [f64; 3]);
        // position and velocity, east / north / up
        let straight: Path = |t| ([15.0 * t, -3.0 * t, 100.0], [15.0, -3.0, 0.0]);
        // a circle of 30 m radius at 12 m/s: 4.8 m/s² (a tight survey turn)
        let turn: Path = |t| {
            let w = 12.0 / 30.0;
            (
                [30.0 * (w * t).cos(), 30.0 * (w * t).sin(), 100.0],
                [-12.0 * (w * t).sin(), 12.0 * (w * t).cos(), 0.0],
            )
        };
        // speeding up along a line at 2 m/s² while climbing at 1 m/s² (start of a flight line)
        let accel: Path = |t| ([t * t, 0.0, 100.0 + 0.5 * t * t], [2.0 * t, 0.0, t]);
        let mut table = String::new();
        let mut worst = std::collections::BTreeMap::new();
        for (name, f) in [
            ("straight 15 m/s", straight),
            ("turn 4.8 m/s2", turn),
            ("speeding up 2 m/s2", accel),
        ] {
            for hz in [1.0, 5.0] {
                let h = 1.0 / hz;
                let pos = move |t: f64| f(t - 1000.0).0;
                let mut tr = track(&pos, h, (40.0 * hz) as usize);
                // the velocity in the east / north / up of each epoch's own position
                let (blat, blon, _) = ecef_to_llh(BASE);
                let rb = enu_matrix(blat, blon);
                for e in &mut tr {
                    let (lat, lon, _) = ecef_to_llh(e.xyz);
                    let v = crate::linalg::mat3_vec(
                        &enu_matrix(lat, lon),
                        mat3t_vec(&rb, f(e.t - 1000.0).1),
                    );
                    e.vel_enu = [v[0] as f32, v[1] as f32, v[2] as f32];
                }
                let ev: Vec<Event> = (0..200)
                    .map(|k| Event::at(1003.0 + 0.1731 * f64::from(k)))
                    .collect();
                let mut line = format!("{name:<20} {hz:>3} Hz");
                for m in [
                    EventInterp::Linear2,
                    EventInterp::Cubic4,
                    EventInterp::Hermite2,
                ] {
                    let r = interpolate(&tr, &ev, m, 0);
                    assert!(
                        r.iter()
                            .all(|x| x.status == EventStatus::Ok && x.method == m)
                    );
                    let e: Vec<f64> = r.iter().map(|x| err(&pos, x)).collect();
                    let max = e.iter().copied().fold(0.0, f64::max);
                    let rms = (e.iter().map(|x| x * x).sum::<f64>() / e.len() as f64).sqrt();
                    let _ = write!(line, "  {:>9.2} {:>9.2}", 1e3 * rms, 1e3 * max);
                    worst.insert((name, hz as u32, m.name()), max);
                }
                table.push_str(&line);
                table.push('\n');
            }
        }
        println!(
            "track                rate    linear2 rms / max     cubic4 rms / max   hermite2 rms / max  (mm)\n{table}"
        );
        // a straight line at constant speed is exact for all (Hermite to the f32 of its velocities)
        for m in ["linear2", "cubic4"] {
            assert!(worst[&("straight 15 m/s", 1, m)] < 1e-5);
        }
        assert!(worst[&("straight 15 m/s", 1, "hermite2")] < 5e-5);
        // linear: a h² / 8 in a turn — 600 mm at 1 Hz, 24 mm at 5 Hz
        assert!((worst[&("turn 4.8 m/s2", 1, "linear2")] - 0.598).abs() < 0.02);
        assert!((worst[&("turn 4.8 m/s2", 5, "linear2")] - 0.024).abs() < 0.002);
        // cubic: 18 mm at 1 Hz, 0.03 mm at 5 Hz (it goes with the fourth power of the interval)
        assert!(worst[&("turn 4.8 m/s2", 1, "cubic4")] < 0.020);
        assert!(worst[&("turn 4.8 m/s2", 5, "cubic4")] < 5e-5);
        // Hermite with exact velocities: 2 mm at 1 Hz
        assert!(worst[&("turn 4.8 m/s2", 1, "hermite2")] < 2.5e-3);
        // constant acceleration is a quadratic: exact for the cubic and for Hermite, not for linear
        assert!(worst[&("speeding up 2 m/s2", 1, "cubic4")] < 1e-5);
        assert!(worst[&("speeding up 2 m/s2", 1, "hermite2")] < 5e-5);
        assert!(worst[&("speeding up 2 m/s2", 1, "linear2")] > 0.25);
    }

    #[test]
    fn hermite_falls_back_without_velocity_and_prices_its_velocities() {
        let f = |t: f64| [10.0 * (t - 1000.0), 0.0, 50.0];
        let mut tr = track(&f, 1.0, 12);
        let r = &interpolate(&tr, &[Event::at(1004.5)], EventInterp::Hermite2, 0)[0];
        assert_eq!(r.method, EventInterp::Linear2);
        for e in &mut tr {
            e.vel_enu = [10.0, 0.0, 0.0];
        }
        let h = &interpolate(&tr, &[Event::at(1004.5)], EventInterp::Hermite2, 0)[0];
        assert_eq!(h.method, EventInterp::Hermite2);
        assert!(err(&f, h) < 1e-4);
        // 1 s × 0.25 × 5 cm/s of velocity noise on top of the 1 cm of the positions
        assert!(
            (h.sigma_enu[0] - 0.01f64.hypot(0.0125)).abs() < 1e-4,
            "{}",
            h.sigma_enu[0]
        );
    }

    #[test]
    fn hermite_reproduces_a_cubic() {
        let p = |t: f64| [t * t * t - 2.0 * t, 3.0 * t * t, 1.0 - t];
        let v = |t: f64| [3.0 * t * t - 2.0, 6.0 * t, -1.0];
        let x = hermite(0.0, p(0.0), v(0.0), 0.2, p(0.2), v(0.2), 0.13);
        let q = p(0.13);
        assert!((0..3).all(|i| (x[i] - q[i]).abs() < 1e-12));
    }

    #[test]
    fn occupations_average_the_best_epochs() {
        let f = |_t: f64| [1.0, 2.0, 3.0];
        let mut tr = track(&f, 1.0, 20);
        for e in tr.iter_mut().take(8).skip(5) {
            e.status = FixStatus::Float;
        }
        let oc = [
            Occupation {
                name: "P1".into(),
                t0: 1003.0,
                t1: 1010.0,
            },
            Occupation {
                name: "P2".into(),
                t0: 2000.0,
                t1: 2010.0,
            },
        ];
        let r = average_occupations(&tr, &oc);
        assert_eq!(r[0].fix, FixStatus::Fixed);
        assert_eq!((r[0].n, r[0].n_all), (5, 8));
        let d = enu_at(BASE, r[0].xyz);
        assert!((d[0] - 1.0).abs() < 1e-6 && (d[2] - 3.0).abs() < 1e-6);
        assert_eq!(r[1].n, 0);
        assert!(r[1].xyz[0].is_nan());
    }
}
