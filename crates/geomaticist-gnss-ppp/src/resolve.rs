//! PPP-AR on the float solution, with its own `partial_fix` that has `greedy_above`
//! (the engine's `partial_fix` lacks that path: above `greedy_above` ambiguities the one farthest from
//! its integer, in its own σ, is dropped instead of trying every one-out subset).
//!
//! With the analysis centre's phase biases (OSB) and integer-recovering clocks a satellite's ambiguities
//! are integers but for a receiver bias shared by every satellite of a group (one system, one signal
//! pair): the wide lane is fixed per group from the Melbourne-Wübbena means (`ambiguity::wide_lane`),
//! the narrow lane from the float ionosphere-free ambiguities as differences against one arc of the group
//! by LAMBDA with the ratio test, and the fixed differences constrain the float solution.
//!
//! The wide-lane σ rule is `PppOptions::wl_estimator` (`Ar1`, the default: one independent MW value
//! per `ar_wl_corr_s` of span; `Legacy`: one per ten epochs).

use std::collections::{BTreeMap, BTreeSet};

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::text::round_to;
use geomaticist_gnss_engine::ambiguity::{lambda_fix, wide_lane};
use geomaticist_gnss_engine::arr::{B2, F2, I2};
use geomaticist_gnss_engine::dump::{jfx, jsl};
use geomaticist_gnss_engine::linalg::{DMat, dot, norm3};
use geomaticist_gnss_engine::nputil::{argmax, mean, median, std_ddof};
use serde_json::{Map, Value};

use crate::options::{PppOptions, WlEstimator, sys_rank};

/// What `_resolve` returns: the record (`ok`, `why`, `wl`, `n_fixed` …) and, when fixed, the position
/// change and the unscaled fixed covariance `(dx, Qxx)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    /// The AR record as a JSON object, without `dx` / `Qxx`.
    pub meta: Map<String, Value>,
    /// `(dx, Qxx)` when `ok`.
    pub fix: Option<([f64; 3], [[f64; 3]; 3])>,
}

fn not_ok(why: String) -> Resolved {
    let mut m = Map::new();
    m.insert("ok".into(), Value::Bool(false));
    m.insert("why".into(), Value::String(why));
    Resolved { meta: m, fix: None }
}

/// The outcome of [`partial_fix`].
#[derive(Debug, Clone, PartialEq)]
pub struct Fix {
    /// Whether a subset passed.
    pub ok: bool,
    /// Indices (into the input) that are fixed.
    pub keep: Vec<usize>,
    /// Their integers.
    pub fixed: Vec<f64>,
    /// The ratio of the accepted (or last) search.
    pub ratio: f64,
    /// Indices dropped, in order.
    pub dropped: Vec<usize>,
    /// Searches run.
    pub tried: usize,
}

/// LAMBDA with the ratio test; on failure drop one ambiguity at a time until a subset of at least
/// `min_fix` passes.
/// `ratio_min` a number; above `greedy_above` ambiguities the one farthest from its integer (in its own
/// σ when the ratio fails) is dropped instead of searching every one-out subset. `None` when a search
/// finds the covariance not positive definite.
#[must_use]
pub fn partial_fix(
    a: &[f64],
    q: &DMat,
    ratio_min: f64,
    min_fix: usize,
    max_offset: f64,
    greedy_above: Option<usize>,
) -> Option<Fix> {
    let n = a.len();
    let mut idx: Vec<usize> = (0..n).collect();
    let mut dropped: Vec<usize> = Vec::new();
    let mut tried = 0usize;
    if n < min_fix {
        return Some(Fix {
            ok: false,
            keep: Vec::new(),
            fixed: Vec::new(),
            ratio: 0.0,
            dropped,
            tried,
        });
    }
    let run = |sel: &[usize], tried: &mut usize| -> Option<(Vec<f64>, Vec<f64>, f64)> {
        *tried += 1;
        let asub: Vec<f64> = sel.iter().map(|&i| a[i]).collect();
        let (f, s) = lambda_fix(&asub, &q.sub(sel), 2).ok()?;
        let r = if s[0] > 0.0 {
            s[1] / s[0]
        } else {
            f64::INFINITY
        };
        Some((f[0].clone(), s, r))
    };
    let (mut f, mut _s, mut r) = run(&idx, &mut tried)?;
    loop {
        let mut off: Vec<f64> = idx
            .iter()
            .zip(&f)
            .map(|(&i, fv)| (a[i] - fv).abs())
            .collect();
        let off_max = off.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if r >= ratio_min && (max_offset <= 0.0 || off_max <= max_offset) {
            return Some(Fix {
                ok: true,
                keep: idx,
                fixed: f,
                ratio: r,
                dropped,
                tried,
            });
        }
        if idx.len() - 1 < min_fix {
            return Some(Fix {
                ok: false,
                keep: idx,
                fixed: f,
                ratio: r,
                dropped,
                tried,
            });
        }
        let best: (f64, usize, Vec<usize>, Vec<f64>, Vec<f64>);
        if r >= ratio_min || greedy_above.is_some_and(|g| idx.len() > g) {
            if r < ratio_min {
                for (k, &i) in idx.iter().enumerate() {
                    off[k] /= q.at(i, i).max(1e-12).sqrt();
                }
            }
            let qi = argmax(&off);
            let sel: Vec<usize> = idx
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != qi)
                .map(|(_, &v)| v)
                .collect();
            let (f2, s2, r2) = run(&sel, &mut tried)?;
            best = (r2, qi, sel, f2, s2);
        } else {
            let mut b: Option<(f64, usize, Vec<usize>, Vec<f64>, Vec<f64>)> = None;
            for qi in 0..idx.len() {
                let sel: Vec<usize> = idx
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| *k != qi)
                    .map(|(_, &v)| v)
                    .collect();
                let (f2, s2, r2) = run(&sel, &mut tried)?;
                if b.as_ref().is_none_or(|bb| r2 > bb.0) {
                    b = Some((r2, qi, sel, f2, s2));
                }
            }
            best = b?;
        }
        let (r2, qi, sel, f2, s2) = best;
        dropped.push(idx[qi]);
        idx = sel;
        r = r2;
        f = f2;
        _s = s2;
    }
}

/// PPP-AR on the float solution (`grp`: the per-satellite group strings).
#[must_use]
pub fn resolve(
    opts: &PppOptions,
    t: &[f64],
    arc: &I2,
    use_p: &B2,
    mw: &F2,
    osb4: &[Option<[f64; 4]>],
    lam1s: &[f64],
    lam2s: &[f64],
    sysc: &[char],
    amb0: &F2,
    sol: &[f64],
    ninv: &DMat,
    sigma0: f64,
    i_a: usize,
    namb: usize,
    grp: &[String],
) -> Resolved {
    let (ne, ns) = (arc.ne, arc.ns);
    // cells of each arc with phase used, row-major
    let mut cells: Vec<Vec<(usize, usize)>> = vec![Vec::new(); namb];
    for k in 0..ne {
        for j in 0..ns {
            let q = arc.at(k, j);
            if q >= 0 && (q as usize) < namb && use_p.at(k, j) {
                cells[q as usize].push((k, j));
            }
        }
    }
    // per candidate arc: (q, j, mw mean, mw σ)
    let mut info: Vec<(usize, usize, f64, f64)> = Vec::new();
    for q in 0..namb {
        let kk = &cells[q];
        if kk.is_empty() {
            continue;
        }
        let j = kk[0].1;
        let Some([vc1, vc2, vl1, vl2]) = osb4[j] else {
            continue;
        };
        if !opts.ar_systems.contains(sysc[j]) {
            continue;
        }
        let kmin = kk.iter().map(|c| c.0).min().unwrap_or(0);
        let kmax = kk.iter().map(|c| c.0).max().unwrap_or(0);
        if t[kmax] - t[kmin] < opts.ar_min_arc_s {
            continue;
        }
        let (f1, f2) = (CLIGHT / lam1s[j], CLIGHT / lam2s[j]);
        let lw = CLIGHT / (f1 - f2);
        let bias = ((f1 * vl1 - f2 * vl2) / (f1 - f2) - (f1 * vc1 + f2 * vc2) / (f1 + f2)) / lw;
        let w: Vec<f64> = kk
            .iter()
            .map(|&(k, _)| mw.at(k, j) - bias)
            .filter(|v| v.is_finite())
            .collect();
        if w.len() < 20 {
            continue;
        }
        let med = median(&w);
        let w: Vec<f64> = w.into_iter().filter(|v| (v - med).abs() < 1.5).collect();
        if w.len() < 20 {
            continue;
        }
        // the MW errors are correlated over minutes (code multipath)
        let span = t[kmax] - t[kmin];
        let n_eff = match opts.wl_estimator {
            // ar_wl_corr_s: one independent value per correlation time (B4)
            WlEstimator::Ar1 => (span / opts.ar_wl_corr_s).max(1.0),
            // the rule before B4: one in ten epochs, whatever the rate
            WlEstimator::Legacy => (w.len() as f64 / 10.0).max(1.0),
        };
        let sig = (std_ddof(&w, 0) / n_eff.sqrt()).max(0.02);
        info.push((q, j, mean(&w), sig));
    }
    if info.len() < opts.ar_min_fix + 1 {
        return not_ok(format!(
            "only {} arc(s) long enough with phase biases",
            info.len()
        ));
    }
    // 1. wide lane, per group
    let mut nwl: BTreeMap<usize, i64> = BTreeMap::new();
    let (mut wl_tried, mut wl_fixed) = (0usize, 0usize);
    let groups: BTreeSet<&String> = info.iter().map(|i| &grp[i.1]).collect();
    for s in groups {
        let sel: Vec<usize> = (0..info.len()).filter(|&k| &grp[info[k].1] == s).collect();
        let mws: Vec<f64> = sel.iter().map(|&k| info[k].2).collect();
        let sgs: Vec<f64> = sel.iter().map(|&k| info[k].3).collect();
        let (n, ok, _b, _off) = wide_lane(&mws, &sgs, opts.ar_wl_frac, opts.ar_wl_sig);
        wl_tried += sel.len();
        for (i, &k) in sel.iter().enumerate() {
            if ok[i] {
                nwl.insert(info[k].0, n[i]);
                wl_fixed += 1;
            }
        }
    }
    let wl_json = Value::Array(vec![Value::from(wl_fixed), Value::from(wl_tried)]);
    // 2. narrow lane: the float N1 of each wide-lane-fixed arc
    let mut arcs: Vec<usize> = Vec::new();
    let mut n1: Vec<f64> = Vec::new();
    let mut scale: Vec<f64> = Vec::new();
    let mut const_m: Vec<f64> = Vec::new();
    let mut sys_of: Vec<&String> = Vec::new();
    for &(q, j, _m, _s) in &info {
        let Some(&nw) = nwl.get(&q) else { continue };
        let (f1, f2) = (CLIGHT / lam1s[j], CLIGHT / lam2s[j]);
        let lam_nl = CLIGHT / (f1 + f2);
        let bl2 = CLIGHT * f2 / (f1 * f1 - f2 * f2);
        let k0 = (0..ne).find(|&k| arc.at(k, j) == q as i64).unwrap_or(0);
        let a0 = amb0.at(k0, j);
        let b = a0 + sol[i_a + q]; // the arc's whole ionosphere-free ambiguity (m)
        arcs.push(q);
        n1.push((b - bl2 * nw as f64) / lam_nl);
        scale.push(lam_nl);
        const_m.push(bl2 * nw as f64 - a0);
        sys_of.push(&grp[j]);
    }
    // differences against the longest-observed arc of each group (one per group: its receiver bias)
    let mut d_rows: Vec<(usize, usize)> = Vec::new(); // (i, ref)
    let mut ref_of: BTreeMap<&String, usize> = BTreeMap::new();
    let sys_set: BTreeSet<&String> = sys_of.iter().copied().collect();
    for s in sys_set {
        let idx: Vec<usize> = (0..sys_of.len()).filter(|&i| sys_of[i] == s).collect();
        if idx.len() < 2 {
            continue;
        }
        let mut rf = idx[0];
        let mut best = 0usize;
        for &i in &idx {
            let c = cells[arcs[i]].len();
            if c > best || i == idx[0] {
                best = c;
                rf = i;
            }
        }
        ref_of.insert(s, rf);
        for &i in &idx {
            if i != rf {
                d_rows.push((i, rf));
            }
        }
    }
    let fail_wl = |why: String, extra: &[(&str, Value)]| -> Resolved {
        let mut r = not_ok(why);
        r.meta.insert("wl".into(), wl_json.clone());
        for (k, v) in extra {
            r.meta.insert((*k).into(), v.clone());
        }
        r
    };
    if d_rows.len() < opts.ar_min_fix {
        return fail_wl(
            format!("wide lane fixed on {wl_fixed} of {wl_tried} arcs: too few to go on"),
            &[],
        );
    }
    // float differences (cycles) and their covariance; parameters -> N1: 1 / λ_NL each
    let na = arcs.len();
    let s02 = sigma0 * sigma0;
    let mut qn = DMat::zeros(na);
    for i in 0..na {
        for j in 0..na {
            let v = ((1.0 / scale[i]) * ninv.at(i_a + arcs[i], i_a + arcs[j])) * (1.0 / scale[j]);
            *qn.at_mut(i, j) = v * s02
                + if i == j {
                    opts.ar_nl_floor * opts.ar_nl_floor
                } else {
                    0.0
                };
        }
    }
    let nr = d_rows.len();
    let a: Vec<f64> = d_rows.iter().map(|&(i, r)| n1[i] - n1[r]).collect();
    let mut qa = DMat::zeros(nr);
    for (r, &(i, rf)) in d_rows.iter().enumerate() {
        for (s, &(j, rs)) in d_rows.iter().enumerate() {
            let m_j = qn.at(i, j) - qn.at(rf, j);
            let m_r = qn.at(i, rs) - qn.at(rf, rs);
            *qa.at_mut(r, s) = m_j - m_r;
        }
    }
    let Some(fx) = partial_fix(
        &a,
        &qa,
        opts.ar_ratio,
        opts.ar_min_fix,
        opts.ar_nl_offset,
        Some(12),
    ) else {
        return fail_wl(
            "narrow lane: the covariance is not positive definite".into(),
            &[("n_cand", Value::from(na))],
        );
    };
    let n_cand = a.len();
    if !fx.ok {
        return fail_wl(
            format!(
                "narrow lane: no set passes the ratio test (ratio {:.1})",
                fx.ratio
            ),
            &[("n_cand", Value::from(n_cand)), ("ratio", jfx(fx.ratio))],
        );
    }
    // 3. constrain: rows G (in the ambiguity parameters, metres) and their targets
    let nf = fx.keep.len();
    let mut gi: Vec<(usize, usize)> = Vec::with_capacity(nf); // (col +1, col -1)
    let mut tgt = vec![0.0; nf];
    let mut ahat = vec![0.0; nf];
    for (r, (&row, &nfix)) in fx.keep.iter().zip(&fx.fixed).enumerate() {
        let (i, k) = d_rows[row];
        gi.push((i_a + arcs[i], i_a + arcs[k]));
        tgt[r] = scale[i] * nfix + const_m[i] - const_m[k];
        ahat[r] = sol[i_a + arcs[i]] - sol[i_a + arcs[k]];
    }
    let mut qaa = DMat::zeros(nf);
    for r in 0..nf {
        for s in 0..nf {
            let (pi, mi) = gi[r];
            let (pj, mj) = gi[s];
            let gn_j = ninv.at(pi, pj) - ninv.at(mi, pj);
            let gn_m = ninv.at(pi, mj) - ninv.at(mi, mj);
            *qaa.at_mut(r, s) = gn_j - gn_m;
        }
    }
    let mut qxa = vec![[0.0; 3]; nf]; // qxa[r][i] = Qxa[i][r]
    for r in 0..nf {
        let (p, m) = gi[r];
        for i in 0..3 {
            qxa[r][i] = ninv.at(i, p) - ninv.at(i, m);
        }
    }
    let rhs: Vec<f64> = (0..nf).map(|r| ahat[r] - tgt[r]).collect();
    let Some(kv) = qaa.solve(&rhs) else {
        return fail_wl(
            "narrow lane: the constraint covariance is singular".into(),
            &[("n_cand", Value::from(n_cand)), ("ratio", jfx(fx.ratio))],
        );
    };
    let mut dx = [0.0; 3];
    for i in 0..3 {
        let col: Vec<f64> = (0..nf).map(|r| qxa[r][i]).collect();
        dx[i] = -dot(&col, &kv);
    }
    let cols: Vec<Vec<f64>> = (0..3)
        .map(|i| (0..nf).map(|r| qxa[r][i]).collect())
        .collect();
    let Some(sol_cols) = qaa.solve_many(&cols) else {
        return fail_wl(
            "narrow lane: the constraint covariance is singular".into(),
            &[("n_cand", Value::from(n_cand)), ("ratio", jfx(fx.ratio))],
        );
    };
    let mut qxx = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            qxx[i][j] = ninv.at(i, j) - dot(&cols[i], &sol_cols[j]);
        }
    }
    let mut qsym = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            qsym[i][j] = (qxx[i][j] + qxx[j][i]) / 2.0;
        }
    }
    let chi = dot(&rhs, &kv) / s02.max(1e-12) / nf as f64;
    if norm3(dx) > 0.30 {
        return fail_wl(
            format!("the fix moves the position {:.2} m: not taken", norm3(dx)),
            &[("n_cand", Value::from(n_cand)), ("ratio", jfx(fx.ratio))],
        );
    }
    let mut systems: Vec<char> = ref_of
        .keys()
        .filter_map(|g| g.chars().next())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    systems.sort_by_key(|&c| sys_rank(c));
    let mut m = Map::new();
    m.insert("ok".into(), Value::Bool(true));
    m.insert("n_fixed".into(), Value::from(nf));
    m.insert("n_cand".into(), Value::from(n_cand));
    m.insert("ratio".into(), jfx(fx.ratio));
    m.insert("wl".into(), wl_json);
    m.insert("chi".into(), jfx(round_to(chi, 2)));
    m.insert(
        "systems".into(),
        Value::String(systems.into_iter().collect()),
    );
    m.insert(
        "groups".into(),
        jsl(&ref_of.keys().map(|s| s.as_str()).collect::<Vec<_>>()),
    );
    m.insert("dropped".into(), Value::from(fx.dropped.len()));
    Resolved {
        meta: m,
        fix: Some((dx, qsym)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partial_fix_accepts_an_easy_set_and_drops_a_biased_one() {
        let a = [3.02, -1.98, 7.01, 0.99, 5.0];
        let q = DMat::scaled(&DMat::eye(5), 1e-4);
        let f = partial_fix(&a, &q, 3.0, 4, 0.3, Some(12)).unwrap();
        assert!(f.ok && f.keep.len() == 5);
        assert_eq!(f.fixed, vec![3.0, -2.0, 7.0, 1.0, 5.0]);
        // one ambiguity half a cycle off: the one-out search drops it
        let a2 = [3.02, -1.98, 7.01, 0.99, 5.5];
        let f2 = partial_fix(&a2, &q, 3.0, 4, 0.3, None).unwrap();
        assert!(f2.ok && f2.dropped == vec![4], "{f2:?}");
        // the greedy path (above 2 ambiguities) drops the farthest in its own σ first
        let mut q2 = DMat::scaled(&DMat::eye(5), 1e-4);
        *q2.at_mut(4, 4) = 1e-2;
        let f3 = partial_fix(&a2, &q2, 3.0, 4, 0.3, Some(2)).unwrap();
        assert!(f3.ok && f3.dropped == vec![4], "{f3:?}");
    }
}
