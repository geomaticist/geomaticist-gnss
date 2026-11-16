//! Cycle slips of one satellite's between-receiver single difference: find them, tell a one-epoch spike
//! from a step, test the hypothesis "no slip" directly (B2, the zero test), repair a step whose integer
//! jumps are certain, else start a new arc.
//!
//! The rolling robust scatter is computed per window directly (21 values: a small sort), which is
//! linear in the series length.

use crate::nputil::{mean, median, py_round_i, std_ddof};
use crate::options::SlipOptions;

/// What `process` found.
#[derive(Debug, Clone, PartialEq)]
pub struct SlipResult {
    /// An arc starts here.
    pub new_arc: Vec<bool>,
    /// Epoch dropped (spike).
    pub drop: Vec<bool>,
    /// Added to the SD ionosphere-free phase (repairs), m.
    pub corr_if: Vec<f64>,
    /// The same repairs on the GF phase, m.
    pub corr_gf: Vec<f64>,
    /// And on the Melbourne-Wübbena, cycles.
    pub corr_mw: Vec<f64>,
    /// (epoch, ΔN1, ΔN2).
    pub repairs: Vec<(usize, i64, i64)>,
    /// Event counts: start, gap, hard, step, repaired, no_slip, spike, wl_link.
    pub events: SlipEvents,
    /// An arc starts here whose wide lane is the previous piece's (ΔNw = 0 across the break).
    pub wl_link: Vec<bool>,
}

/// The event counters of one satellite (and, summed, of a solve).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SlipEvents {
    /// Arcs started at the first usable epoch.
    pub start: usize,
    /// Arcs started after a gap.
    pub gap: usize,
    /// Arcs started at a hard break.
    pub hard: usize,
    /// Arcs started at a step.
    pub step: usize,
    /// Steps repaired.
    pub repaired: usize,
    /// Flags that were no slip.
    pub no_slip: usize,
    /// Spikes dropped.
    pub spike: usize,
    /// Pieces linked for the wide lane.
    pub wl_link: usize,
}

impl SlipEvents {
    /// Sum another satellite's counts into this one.
    pub fn add(&mut self, o: &SlipEvents) {
        self.start += o.start;
        self.gap += o.gap;
        self.hard += o.hard;
        self.step += o.step;
        self.repaired += o.repaired;
        self.no_slip += o.no_slip;
        self.spike += o.spike;
        self.wl_link += o.wl_link;
    }

    /// The counts in a fixed key order, for the notes and dumps.
    #[must_use]
    pub fn pairs(&self) -> [(&'static str, usize); 8] {
        [
            ("start", self.start),
            ("gap", self.gap),
            ("hard", self.hard),
            ("step", self.step),
            ("repaired", self.repaired),
            ("no_slip", self.no_slip),
            ("spike", self.spike),
            ("wl_link", self.wl_link),
        ]
    }
}

/// `max(a, b)`: `b` only when `b > a` (so a NaN first argument stays NaN).
#[inline]
fn py_max(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// Value at `at` of a straight line through `(tt, yy)` and the fit's scatter (NaN under 3 points). The
/// line is parametrised at `tt[-1]` as the design matrix `[1, tt - tt[-1]]`.
#[must_use]
pub fn line(tt: &[f64], yy: &[f64], at: f64) -> (f64, f64) {
    let n = tt.len();
    if n == 1 {
        return (yy[0], f64::NAN);
    }
    let t0 = tt[n - 1];
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0, 0.0, 0.0, 0.0);
    for (&t, &y) in tt.iter().zip(yy) {
        let x = t - t0;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let nf = n as f64;
    let den = nf * sxx - sx * sx;
    let c1 = if den != 0.0 {
        (nf * sxy - sx * sy) / den
    } else {
        0.0
    };
    let c0 = (sy - c1 * sx) / nf;
    let s = if n >= 3 {
        let mut ss = 0.0;
        for (&t, &y) in tt.iter().zip(yy) {
            let r = y - (c0 + c1 * (t - t0));
            ss += r * r;
        }
        (ss / (n - 2) as f64).sqrt()
    } else {
        f64::NAN
    };
    (c0 + c1 * (at - t0), s)
}

/// `1.4826 ×` the median absolute deviation over ±`w` points (NaN where under 5 finite points).
#[must_use]
pub fn rolling_mad(v: &[f64], w: usize) -> Vec<f64> {
    let n = v.len();
    let mut out = vec![f64::NAN; n];
    let mut seg: Vec<f64> = Vec::with_capacity(2 * w + 1);
    let mut dev: Vec<f64> = Vec::with_capacity(2 * w + 1);
    for i in 0..n {
        seg.clear();
        let lo = i.saturating_sub(w);
        let hi = (i + w + 1).min(n);
        seg.extend(v[lo..hi].iter().copied().filter(|x| x.is_finite()));
        if seg.len() >= 5 {
            let m = median(&seg);
            dev.clear();
            dev.extend(seg.iter().map(|x| (x - m).abs()));
            out[i] = 1.4826 * median(&dev);
        }
    }
    out
}

/// Per epoch: the ionosphere's activity (robust scatter of the GF second difference) and the MW noise
/// (robust scatter of its first difference / √2), over ±`w` usable epochs.
#[must_use]
pub fn activity(gf: &[f64], mw: &[f64], idx: &[usize], w: usize) -> (Vec<f64>, Vec<f64>) {
    let ne = gf.len();
    let mut a_gf = vec![f64::NAN; ne];
    let mut a_mw = vec![f64::NAN; ne];
    if idx.len() >= 3 {
        let g: Vec<f64> = idx.iter().map(|&i| gf[i]).collect();
        let mut d2 = vec![f64::NAN; idx.len()];
        for i in 2..idx.len() {
            d2[i] = g[i] - 2.0 * g[i - 1] + g[i - 2];
        }
        let r = rolling_mad(&d2, w);
        for (k, &i) in idx.iter().enumerate() {
            a_gf[i] = r[k];
        }
        let m: Vec<f64> = idx.iter().map(|&i| mw[i]).collect();
        let mut d1 = vec![f64::NAN; idx.len()];
        for i in 1..idx.len() {
            d1[i] = m[i] - m[i - 1];
        }
        let r = rolling_mad(&d1, w);
        for (k, &i) in idx.iter().enumerate() {
            a_mw[i] = r[k] / std::f64::consts::SQRT_2;
        }
    }
    (a_gf, a_mw)
}

/// One satellite's SD series (all `ne` long); `ok` marks usable epochs, `hard` a break that is never
/// repaired.
#[must_use]
pub fn process(
    t: &[f64],
    gf_in: &[f64],
    mw_in: &[f64],
    ok: &[bool],
    lli: &[bool],
    hard: &[bool],
    lam1: f64,
    lam2: f64,
    a_if: f64,
    b_if: f64,
    o: &SlipOptions,
) -> SlipResult {
    let ne = t.len();
    let mut new_arc = vec![false; ne];
    let mut drop = vec![false; ne];
    let mut corr_if = vec![0.0; ne];
    let mut corr_gf = vec![0.0; ne];
    let mut corr_mw = vec![0.0; ne];
    let mut wl_link = vec![false; ne];
    let mut ev = SlipEvents::default();
    let mut repairs = Vec::new();
    let mut gf = gf_in.to_vec();
    let mut mw = mw_in.to_vec();
    let idx: Vec<usize> = (0..ne).filter(|&i| ok[i] && gf[i].is_finite()).collect();
    let mut arc: Vec<usize> = Vec::new();
    let dl = lam1 - lam2;
    let (act_gf, act_mw) = activity(&gf, &mw, &idx, o.activity_window);

    let gfw = o.gf_window;
    let mut i = 0usize;
    while i < idx.len() {
        let k = idx[i];
        if arc.is_empty() {
            new_arc[k] = true;
            ev.start += 1;
            arc = vec![k];
            i += 1;
            continue;
        }
        let prev = arc[arc.len() - 1];
        if hard[k] {
            new_arc[k] = true;
            ev.hard += 1;
            arc = vec![k];
            i += 1;
            continue;
        }
        if k - prev > o.max_gap + 1 {
            new_arc[k] = true;
            ev.gap += 1;
            arc = vec![k];
            i += 1;
            continue;
        }
        // stats()
        let e_b: Vec<usize> = arc[arc.len().saturating_sub(30)..].to_vec();
        let m_b: Vec<f64> = e_b
            .iter()
            .map(|&e| mw[e])
            .filter(|x| x.is_finite())
            .collect();
        let g: Vec<usize> = arc[arc.len().saturating_sub(gfw)..].to_vec();
        let (mu, s_eff, n_b) = if m_b.is_empty() {
            (f64::NAN, o.mw_sigma_min, 0usize)
        } else {
            let mu = mean(&m_b);
            let s = if m_b.len() >= 3 {
                std_ddof(&m_b, 1)
            } else {
                f64::NAN
            };
            let mut s_eff = py_max(
                if s.is_finite() {
                    s
                } else {
                    o.mw_floor / o.mw_k
                },
                o.mw_sigma_min,
            );
            let am = act_mw[arc[arc.len() - 1]];
            if am.is_finite() {
                s_eff = py_max(s_eff, am);
            }
            (mu, s_eff, m_b.len())
        };
        let tg: Vec<f64> = g.iter().map(|&e| t[e]).collect();
        let gg: Vec<f64> = g.iter().map(|&e| gf[e]).collect();
        let test = |kk: usize, gf: &[f64], mw: &[f64]| -> (bool, bool) {
            let (p, _sg) = line(&tg, &gg, t[kk]);
            let thr = if act_gf[kk].is_finite() {
                py_max(o.slip_gf, o.gf_k * act_gf[kk])
            } else {
                o.slip_gf
            };
            let f_mw = mw[kk].is_finite()
                && mu.is_finite()
                && (mw[kk] - mu).abs() > py_max(o.mw_k * s_eff, o.mw_floor);
            (f_mw, (gf[kk] - p).abs() > thr)
        };
        let (f_mw, f_gf) = test(k, &gf, &mw);
        if !(f_mw || f_gf || lli[k]) {
            arc.push(k);
            i += 1;
            continue;
        }
        // spike: the next epoch agrees with the arc as it was
        if o.spikes && !lli[k] && i + 1 < idx.len() {
            let k2 = idx[i + 1];
            if !hard[k2] && !lli[k2] && k2 - prev <= o.max_gap + 1 {
                let (m2, g2) = test(k2, &gf, &mw);
                if !(m2 || g2) {
                    drop[k] = true;
                    ev.spike += 1;
                    i += 1;
                    continue;
                }
            }
        }
        // a step: measure it from the epochs after
        let mut fixed = false;
        let mut link = false;
        if (o.repair || o.zero_test) && n_b >= 3 {
            let n_after = o
                .repair_after
                .max(if o.zero_test { o.zero_after } else { 0 });
            let mut after_all: Vec<usize> = vec![k];
            let mut first_lli: Option<usize> = None;
            let hi = (i + n_after).min(idx.len());
            for &k3 in &idx[(i + 1).min(hi)..hi] {
                if hard[k3] || k3 - after_all[after_all.len() - 1] > o.max_gap + 1 {
                    break;
                }
                if lli[k3] && first_lli.is_none() {
                    first_lli = Some(after_all.len());
                }
                after_all.push(k3);
            }
            let e_bw: Vec<usize> = e_b.iter().copied().filter(|&e| mw[e].is_finite()).collect();
            let mw_bw: Vec<f64> = e_bw.iter().map(|&e| mw[e]).collect();
            if o.zero_test {
                let aft: Vec<usize> = after_all[..after_all.len().min(o.zero_after)].to_vec();
                let aw: Vec<f64> = aft
                    .iter()
                    .map(|&e| mw[e])
                    .filter(|x| x.is_finite())
                    .collect();
                let (pb, sgb) = line(&tg, &gg, t[k]);
                let ta: Vec<f64> = aft.iter().rev().map(|&e| t[e]).collect();
                let ga: Vec<f64> = aft.iter().rev().map(|&e| gf[e]).collect();
                let (pa, sga) = line(&ta, &ga, t[k]);
                let sgb = if sgb.is_finite() { sgb } else { o.slip_gf };
                let sga_ = if sga.is_finite() { sga } else { sgb };
                let d_g = pa - pb;
                let s_g = py_max(sgb.hypot(sga_), o.zero_gf_floor);
                let gf_zero = d_g.abs() <= o.zero_gf_k * s_g
                    && d_g.abs() <= dl.abs() / 2.0
                    && (!sga.is_finite() || sga <= o.zero_fit_max);
                let mut mw_zero = false;
                if aw.len() >= o.zero_mw_n && mw_bw.len() >= 3 {
                    let d_w = mean(&aw) - mean(&mw_bw);
                    let s_w = (py_max(std_ddof(&aw, 1), o.mw_sigma_min) / (aw.len() as f64).sqrt())
                        .hypot(s_eff / (mw_bw.len() as f64).sqrt());
                    mw_zero = d_w.abs() <= o.zero_mw && s_w <= o.zero_mw_sigma;
                    link = o.wl_link && d_w.abs() <= o.wl_link_mw && s_w <= o.zero_mw_sigma;
                }
                if gf_zero && mw_zero {
                    ev.no_slip += 1;
                    arc.push(k);
                    fixed = true;
                }
            }
            let cut = first_lli.unwrap_or(after_all.len()).min(after_all.len());
            let after: Vec<usize> = after_all[..cut.min(o.repair_after)].to_vec();
            let after_w: Vec<f64> = after
                .iter()
                .map(|&e| mw[e])
                .filter(|x| x.is_finite())
                .collect();
            if !fixed && o.repair && after.len() >= 2 && after_w.len() >= 2 && mw_bw.len() >= 3 {
                let d_w = mean(&after_w) - mean(&mw_bw);
                let s_after = std_ddof(&after_w, 1);
                let s_w = (py_max(s_after, o.mw_sigma_min) / (after_w.len() as f64).sqrt())
                    .hypot(s_eff / (mw_bw.len() as f64).sqrt());
                let (pb, sgb) = line(&tg, &gg, t[k]);
                let ta: Vec<f64> = after.iter().rev().map(|&e| t[e]).collect();
                let ga: Vec<f64> = after.iter().rev().map(|&e| gf[e]).collect();
                let (pa, sga) = line(&ta, &ga, t[k]);
                let sgb = if sgb.is_finite() { sgb } else { o.slip_gf };
                let sga = if sga.is_finite() { sga } else { o.slip_gf };
                let d_g = pa - pb;
                let s_g = sgb.hypot(sga * if after.len() < 3 { 1.5 } else { 1.0 });
                let nw = py_round_i(d_w);
                let n1f = (d_g - lam2 * nw as f64) / dl;
                let n1 = py_round_i(n1f);
                let s_n1 = s_g / dl.abs();
                let ok_w = (d_w - nw as f64).abs() <= o.repair_frac
                    && s_w <= o.repair_sigma * 2.0
                    && (nw.abs() as f64) < 1e5;
                let ok_1 = (n1f - n1 as f64).abs() <= o.repair_frac && s_n1 <= o.repair_sigma;
                if ok_w && ok_1 {
                    let n2 = n1 - nw;
                    if n1 == 0 && n2 == 0 {
                        ev.no_slip += 1;
                    } else {
                        let dgf = lam1 * n1 as f64 - lam2 * n2 as f64;
                        let dif = a_if * lam1 * n1 as f64 - b_if * lam2 * n2 as f64;
                        for e in k..ne {
                            gf[e] -= dgf;
                            mw[e] -= nw as f64;
                            corr_gf[e] -= dgf;
                            corr_mw[e] -= nw as f64;
                            corr_if[e] -= dif;
                        }
                        repairs.push((k, n1, n2));
                        ev.repaired += 1;
                    }
                    arc.push(k);
                    fixed = true;
                }
            }
        }
        if !fixed {
            new_arc[k] = true;
            ev.step += 1;
            arc = vec![k];
            if link {
                wl_link[k] = true;
                ev.wl_link += 1;
            }
        }
        i += 1;
    }
    SlipResult {
        new_arc,
        drop,
        corr_if,
        corr_gf,
        corr_mw,
        repairs,
        events: ev,
        wl_link,
    }
}

/// The likeliest single step in a residual series `r` (σ per point `sig`): (index where the new level
/// starts, size, whole-series statistic, local statistic, local size), or `None`.
#[must_use]
pub fn find_step(
    r: &[f64],
    sig: &[f64],
    min_seg: usize,
    window: usize,
) -> Option<(usize, f64, f64, f64, f64)> {
    let n = r.len();
    if n < 2 * min_seg {
        return None;
    }
    let mut cw = vec![0.0; n + 1];
    let mut cr = vec![0.0; n + 1];
    for i in 0..n {
        let w = 1.0 / (sig[i] * sig[i]);
        cw[i + 1] = cw[i] + w;
        cr[i + 1] = cr[i] + w * r[i];
    }
    let mut best_k = 0usize;
    let mut best_z = f64::NEG_INFINITY;
    let mut best_d = 0.0;
    let mut first = true;
    for p in min_seg..=n - min_seg {
        let (wl, wr) = (cw[p], cw[n] - cw[p]);
        let (ml, mr) = (cr[p] / wl, (cr[n] - cr[p]) / wr);
        let z = (mr - ml).abs() / (1.0 / wl + 1.0 / wr).sqrt();
        if first || z > best_z {
            first = false;
            best_z = z;
            best_k = p;
            best_d = mr - ml;
        }
    }
    let q = best_k;
    let lo = q.saturating_sub(window);
    let hi = (q + window).min(n);
    let (wl2, wr2) = (cw[q] - cw[lo], cw[hi] - cw[q]);
    let d_loc = (cr[hi] - cr[q]) / wr2 - (cr[q] - cr[lo]) / wl2;
    let z_loc = d_loc.abs() / (1.0 / wl2 + 1.0 / wr2).sqrt();
    Some((q, best_d, best_z, z_loc, d_loc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_and_mad() {
        let (p, s) = line(&[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0], 3.0);
        assert!((p - 7.0).abs() < 1e-12 && s.abs() < 1e-12);
        let v: Vec<f64> = (0..30)
            .map(|i| if i == 15 { 100.0 } else { (i % 3) as f64 })
            .collect();
        let r = rolling_mad(&v, 10);
        assert!(r[15] < 2.0 && r[0].is_finite());
        assert!(rolling_mad(&[1.0, 2.0, 3.0], 10)[0].is_nan());
    }

    #[test]
    fn a_clean_series_is_one_arc_and_a_step_is_found() {
        let ne = 60;
        let t: Vec<f64> = (0..ne).map(|i| 30.0 * i as f64).collect();
        let gf = vec![0.1; ne];
        let mw = vec![0.3; ne];
        let ok = vec![true; ne];
        let lli = vec![false; ne];
        let hard = vec![false; ne];
        let o = SlipOptions::default();
        let r = process(
            &t, &gf, &mw, &ok, &lli, &hard, 0.19, 0.244, 2.5457, 1.5457, &o,
        );
        assert_eq!(r.events.start, 1);
        assert_eq!(r.new_arc.iter().filter(|&&b| b).count(), 1);
        let res: Vec<f64> = (0..40).map(|i| if i >= 20 { 0.05 } else { 0.0 }).collect();
        let sig = vec![0.003; 40];
        let (q, d, z, zl, _) = find_step(&res, &sig, 4, 6).unwrap();
        assert_eq!(q, 20);
        assert!((d - 0.05).abs() < 1e-12 && z > 8.0 && zl > 6.0);
    }
}
