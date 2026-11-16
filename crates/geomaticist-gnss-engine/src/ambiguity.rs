//! Integer ambiguities: the LAMBDA method (decorrelation and search), the ratio test, partial fixing, and
//! the wide-lane fix from Melbourne-Wübbena averages.
//!
//! LAMBDA after Teunissen (1995), de Jonge & Tiberius (1996), the search as MLAMBDA
//! (Chang, Yang & Zhou 2005); the same algorithm RTKLIB uses, written here from the papers.
//! Rounding to an integer is half to even everywhere.

use crate::linalg::DMat;
use crate::nputil::{argmax, np_sum};

/// The ambiguity covariance is not positive definite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NotPositiveDefinite;

/// `Q = Lᵀ diag(D) L` with `L` unit lower triangular (from the last row up).
pub fn ldl(q: &DMat) -> Result<(DMat, Vec<f64>), NotPositiveDefinite> {
    let n = q.n;
    let mut a = q.clone();
    let mut l = DMat::zeros(n);
    let mut d = vec![0.0; n];
    for i in (0..n).rev() {
        d[i] = a.at(i, i);
        if d[i] <= 0.0 {
            return Err(NotPositiveDefinite);
        }
        let s = d[i].sqrt();
        for k in 0..=i {
            *l.at_mut(i, k) = a.at(i, k) / s;
        }
        for j in 0..i {
            let lij = l.at(i, j);
            for k in 0..=j {
                let v = l.at(i, k) * lij;
                *a.at_mut(j, k) -= v;
            }
        }
        let lii = l.at(i, i);
        for k in 0..=i {
            *l.at_mut(i, k) /= lii;
        }
    }
    Ok((l, d))
}

fn gauss(l: &mut DMat, z: &mut DMat, i: usize, j: usize) {
    let n = l.n;
    let mu = l.at(i, j).round_ties_even();
    if mu != 0.0 {
        for r in i..n {
            let v = mu * l.at(r, i);
            *l.at_mut(r, j) -= v;
        }
        for r in 0..n {
            let v = mu * z.at(r, i);
            *z.at_mut(r, j) -= v;
        }
    }
}

fn perm(l: &mut DMat, d: &mut [f64], j: usize, delta: f64, z: &mut DMat) {
    let n = l.n;
    let eta = d[j] / delta;
    let lam = d[j + 1] * l.at(j + 1, j) / delta;
    d[j] = eta * d[j + 1];
    d[j + 1] = delta;
    let a0: Vec<f64> = (0..j).map(|k| l.at(j, k)).collect();
    let a1: Vec<f64> = (0..j).map(|k| l.at(j + 1, k)).collect();
    let ljj = l.at(j + 1, j);
    for k in 0..j {
        *l.at_mut(j, k) = -ljj * a0[k] + a1[k];
        *l.at_mut(j + 1, k) = eta * a0[k] + lam * a1[k];
    }
    *l.at_mut(j + 1, j) = lam;
    for r in j + 2..n {
        let (x, y) = (l.at(r, j), l.at(r, j + 1));
        *l.at_mut(r, j) = y;
        *l.at_mut(r, j + 1) = x;
    }
    for r in 0..n {
        let (x, y) = (z.at(r, j), z.at(r, j + 1));
        *z.at_mut(r, j) = y;
        *z.at_mut(r, j + 1) = x;
    }
}

/// Decorrelate in place: returns `Z` with `Zᵀ Q Z = Lᵀ D L` (new `L`, `D`).
pub fn reduction(l: &mut DMat, d: &mut [f64]) -> DMat {
    let n = d.len();
    let mut z = DMat::eye(n);
    if n < 2 {
        return z;
    }
    let mut j = n as i64 - 2;
    let mut k = n as i64 - 2;
    while j >= 0 {
        let ju = j as usize;
        if j <= k {
            for i in ju + 1..n {
                gauss(l, &mut z, i, ju);
            }
        }
        let delta = d[ju] + l.at(ju + 1, ju).powi(2) * d[ju + 1];
        if delta + 1e-6 < d[ju + 1] {
            perm(l, d, ju, delta, &mut z);
            k = j;
            j = n as i64 - 2;
        } else {
            j -= 1;
        }
    }
    z
}

#[inline]
fn sgn(x: f64) -> f64 {
    if x <= 0.0 { -1.0 } else { 1.0 }
}

/// The `m` integer vectors nearest `zs` in the metric `(Lᵀ D L)⁻¹`: `(vectors (m of n), squared norms)`
/// sorted ascending.
#[must_use]
pub fn search(l: &DMat, d: &[f64], zs: &[f64], m: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n = d.len();
    let loop_max = 100_000;
    let mut s_mat = DMat::zeros(n);
    let mut dist = vec![0.0; n];
    let mut zb = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut step = vec![0.0; n];
    let mut zn: Vec<Vec<f64>> = vec![vec![0.0; n]; m];
    let mut s = vec![f64::INFINITY; m];
    let (mut nn, mut imax, mut maxdist) = (0usize, 0usize, f64::INFINITY);
    let mut k = n - 1;
    zb[k] = zs[k];
    z[k] = zb[k].round_ties_even();
    let mut y = zb[k] - z[k];
    step[k] = sgn(y);
    for _ in 0..loop_max {
        let newdist = dist[k] + y * y / d[k];
        if newdist < maxdist {
            if k != 0 {
                k -= 1;
                dist[k] = newdist;
                for c in 0..=k {
                    let v = s_mat.at(k + 1, c) + (z[k + 1] - zb[k + 1]) * l.at(k + 1, c);
                    *s_mat.at_mut(k, c) = v;
                }
                zb[k] = zs[k] + s_mat.at(k, k);
                z[k] = zb[k].round_ties_even();
                y = zb[k] - z[k];
                step[k] = sgn(y);
            } else {
                if nn < m {
                    if nn == 0 || newdist > s[imax] {
                        imax = nn;
                    }
                    zn[nn].clone_from(&z);
                    s[nn] = newdist;
                    nn += 1;
                } else {
                    if newdist < s[imax] {
                        zn[imax].clone_from(&z);
                        s[imax] = newdist;
                        imax = argmax(&s);
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0];
                y = zb[0] - z[0];
                step[0] = -step[0] - sgn(step[0]);
            }
        } else {
            if k == n - 1 {
                break;
            }
            k += 1;
            z[k] += step[k];
            y = zb[k] - z[k];
            step[k] = -step[k] - sgn(step[k]);
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| s[a].total_cmp(&s[b]));
    (
        order.iter().map(|&i| zn[i].clone()).collect(),
        order.iter().map(|&i| s[i]).collect(),
    )
}

/// Integer least squares: the best `m` integer vectors for float `a` with covariance `q`, and their
/// squared distances, ascending.
pub fn lambda_fix(
    a: &[f64],
    q: &DMat,
    m: usize,
) -> Result<(Vec<Vec<f64>>, Vec<f64>), NotPositiveDefinite> {
    let n = a.len();
    if n == 0 {
        return Ok((vec![Vec::new(); m], vec![0.0; m]));
    }
    let base: Vec<f64> = a.iter().map(|x| x.round_ties_even()).collect();
    let (mut l, mut d) = ldl(q)?;
    let z = reduction(&mut l, &mut d);
    let diff: Vec<f64> = a.iter().zip(&base).map(|(x, b)| x - b).collect();
    let zs: Vec<f64> = (0..n)
        .map(|j| (0..n).map(|i| z.at(i, j) * diff[i]).sum())
        .collect();
    let (e, s) = search(&l, &d, &zs, m);
    // back: a = Z⁻ᵀ z, i.e. solve Zᵀ F = E for each vector
    let zt = DMat {
        n,
        data: (0..n * n).map(|k| z.at(k % n, k / n)).collect(),
    };
    let cols = zt.solve_many(&e).ok_or(NotPositiveDefinite)?;
    let out: Vec<Vec<f64>> = cols
        .iter()
        .map(|c| {
            c.iter()
                .zip(&base)
                .map(|(f, b)| f.round_ties_even() + b)
                .collect()
        })
        .collect();
    Ok((out, s))
}

/// Bootstrapping success rate of the decorrelated ambiguities (a lower bound).
pub fn success_rate(q: &DMat) -> Result<f64, NotPositiveDefinite> {
    if q.n == 0 {
        return Ok(1.0);
    }
    let (mut l, mut d) = ldl(q)?;
    reduction(&mut l, &mut d);
    let mut p = 1.0;
    for &dd in &d {
        p *= libm::erf(1.0 / (2.0 * (2.0 * dd).sqrt()));
    }
    Ok(p)
}

/// The outcome of a partial fix.
#[derive(Debug, Clone, PartialEq)]
pub struct FixResult {
    /// Whether a subset passed.
    pub ok: bool,
    /// Indices (into the input) that are fixed.
    pub keep: Vec<usize>,
    /// Their integers.
    pub fixed: Vec<f64>,
    /// The ratio of the accepted (or last) search.
    pub ratio: f64,
    /// The two squared norms.
    pub s: (f64, f64),
    /// Indices dropped, in order.
    pub dropped: Vec<usize>,
    /// Searches run.
    pub tried: usize,
    /// `par_fix`: the fixed integer combinations `C a = c` (rows), else `None`.
    pub c_mat: Option<Vec<Vec<f64>>>,
    /// `par_fix`: the right-hand side `c`.
    pub c_vec: Option<Vec<f64>>,
    /// `par_fix`: bootstrapped success rate of the subset.
    pub p: f64,
}

impl FixResult {
    fn none(
        idx: Vec<usize>,
        fixed: Vec<f64>,
        ratio: f64,
        s: (f64, f64),
        dropped: Vec<usize>,
        tried: usize,
    ) -> Self {
        Self {
            ok: false,
            keep: idx,
            fixed,
            ratio,
            s,
            dropped,
            tried,
            c_mat: None,
            c_vec: None,
            p: f64::NAN,
        }
    }
}

fn run_subset(
    a: &[f64],
    q: &DMat,
    sel: &[usize],
) -> Result<(Vec<f64>, Vec<f64>, f64), NotPositiveDefinite> {
    let asub: Vec<f64> = sel.iter().map(|&i| a[i]).collect();
    let qsub = q.sub(sel);
    let (f, s) = lambda_fix(&asub, &qsub, 2)?;
    let r = if s[0] > 0.0 {
        s[1] / s[0]
    } else {
        f64::INFINITY
    };
    Ok((f[0].clone(), s, r))
}

/// LAMBDA with the ratio test on all of `a`; on failure drop, one at a time, the ambiguity whose removal
/// gives the highest ratio, until a subset of at least `min_fix` passes or too few are left. `thr(n)` gives
/// the threshold for a set of `n`.
pub fn partial_fix(
    a: &[f64],
    q: &DMat,
    thr: &dyn Fn(usize) -> f64,
    min_fix: usize,
    max_offset: f64,
) -> Result<FixResult, NotPositiveDefinite> {
    let n = a.len();
    let mut idx: Vec<usize> = (0..n).collect();
    let mut dropped: Vec<usize> = Vec::new();
    let mut tried = 0usize;
    if n < min_fix {
        return Ok(FixResult::none(
            Vec::new(),
            Vec::new(),
            0.0,
            (f64::NAN, f64::NAN),
            dropped,
            0,
        ));
    }
    let (mut f, mut s, mut r) = run_subset(a, q, &idx)?;
    tried += 1;
    loop {
        let off: Vec<f64> = idx
            .iter()
            .zip(&f)
            .map(|(&i, fv)| (a[i] - fv).abs())
            .collect();
        let off_max = off.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if r >= thr(idx.len()) && (max_offset <= 0.0 || off_max <= max_offset) {
            return Ok(FixResult {
                ok: true,
                keep: idx,
                fixed: f,
                ratio: r,
                s: (s[0], s[1]),
                dropped,
                tried,
                c_mat: None,
                c_vec: None,
                p: f64::NAN,
            });
        }
        if idx.len() - 1 < min_fix {
            return Ok(FixResult::none(idx, f, r, (s[0], s[1]), dropped, tried));
        }
        let best;
        if r >= thr(idx.len()) {
            // the set passes the ratio but an offset is too large: drop the farthest from its integer
            let qi = argmax(&off);
            let sel: Vec<usize> = idx
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != qi)
                .map(|(_, &v)| v)
                .collect();
            let (f2, s2, r2) = run_subset(a, q, &sel)?;
            tried += 1;
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
                let (f2, s2, r2) = run_subset(a, q, &sel)?;
                tried += 1;
                if b.as_ref().is_none_or(|bb| r2 > bb.0) {
                    b = Some((r2, qi, sel, f2, s2));
                }
            }
            best = b.ok_or(NotPositiveDefinite)?;
        }
        let (r2, qi, sel, f2, s2) = best;
        dropped.push(idx[qi]);
        idx = sel;
        r = r2;
        f = f2;
        s = s2;
    }
}

/// Partial ambiguity resolution in the decorrelated space: the largest bootstrapping-ordered subset with
/// success rate ≥ `p0`, searched on its own and accepted by the ratio test. What is fixed is integer
/// combinations `C a = c`.
pub fn par_fix(
    a: &[f64],
    q: &DMat,
    p0: f64,
    thr: &dyn Fn(usize) -> f64,
    min_fix: usize,
) -> Result<FixResult, NotPositiveDefinite> {
    let n = a.len();
    if n < min_fix {
        return Ok(FixResult::none(
            Vec::new(),
            Vec::new(),
            0.0,
            (f64::NAN, f64::NAN),
            Vec::new(),
            0,
        ));
    }
    let base: Vec<f64> = a.iter().map(|x| x.round_ties_even()).collect();
    let (mut l, mut d) = ldl(q)?;
    let z = reduction(&mut l, &mut d);
    let diff: Vec<f64> = a.iter().zip(&base).map(|(x, b)| x - b).collect();
    let zs: Vec<f64> = (0..n)
        .map(|j| (0..n).map(|i| z.at(i, j) * diff[i]).sum())
        .collect();
    let mut ps = Vec::with_capacity(n);
    let mut acc = 1.0;
    for &dd in d.iter().rev() {
        acc *= libm::erf(1.0 / (2.0 * (2.0 * dd).sqrt()));
        ps.push(acc);
    }
    let mut m = ps.iter().filter(|&&p| p >= p0).count();
    let mut tried = 0usize;
    while m >= min_fix {
        let sub: Vec<usize> = (n - m..n).collect();
        let lsub = l.sub(&sub);
        let dsub: Vec<f64> = sub.iter().map(|&i| d[i]).collect();
        let zsub: Vec<f64> = sub.iter().map(|&i| zs[i]).collect();
        let (e, s) = search(&lsub, &dsub, &zsub, 2);
        tried += 1;
        let r = if s[0] > 0.0 {
            s[1] / s[0]
        } else {
            f64::INFINITY
        };
        if r >= thr(m) {
            // C = Z[:, sub].T  (m rows of n)
            let c: Vec<Vec<f64>> = sub
                .iter()
                .map(|&j| (0..n).map(|i| z.at(i, j)).collect())
                .collect();
            let fixed: Vec<f64> = e[0].iter().map(|x| x.round_ties_even()).collect();
            let cvec: Vec<f64> = c
                .iter()
                .zip(&fixed)
                .map(|(row, f)| f + row.iter().zip(&base).map(|(x, b)| x * b).sum::<f64>())
                .collect();
            return Ok(FixResult {
                ok: true,
                keep: sub,
                fixed,
                ratio: r,
                s: (s[0], s[1]),
                dropped: Vec::new(),
                tried,
                c_mat: Some(c),
                c_vec: Some(cvec),
                p: ps[m - 1],
            });
        }
        m -= 1;
    }
    Ok(FixResult::none(
        Vec::new(),
        Vec::new(),
        0.0,
        (f64::NAN, f64::NAN),
        Vec::new(),
        tried,
    ))
}

/// Wide-lane integers of SD arcs from their Melbourne-Wübbena means and σ: the arcs share a receiver bias
/// `b` (the weighted circular mean of the fractional parts of the arcs good enough to vote), arc `i` is fixed
/// to `round(mw_i - b)` when its offset is under `frac_max` and its σ under `sig_max`. Returns
/// `(N (0 where not fixed), fixed, b, offsets)`.
#[must_use]
pub fn wide_lane(
    mw: &[f64],
    sig: &[f64],
    frac_max: f64,
    sig_max: f64,
) -> (Vec<i64>, Vec<bool>, f64, Vec<f64>) {
    let n = mw.len();
    let good: Vec<bool> = (0..n)
        .map(|i| mw[i].is_finite() && sig[i].is_finite() && sig[i] <= sig_max)
        .collect();
    if !good.iter().any(|&g| g) {
        return (vec![0; n], vec![false; n], 0.0, vec![f64::NAN; n]);
    }
    let gi: Vec<usize> = (0..n).filter(|&i| good[i]).collect();
    let w: Vec<f64> = gi.iter().map(|&i| 1.0 / sig[i].max(0.02).powi(2)).collect();
    let re: Vec<f64> = gi
        .iter()
        .zip(&w)
        .map(|(&i, wi)| wi * (2.0 * std::f64::consts::PI * mw[i]).cos())
        .collect();
    let im: Vec<f64> = gi
        .iter()
        .zip(&w)
        .map(|(&i, wi)| wi * (2.0 * std::f64::consts::PI * mw[i]).sin())
        .collect();
    let mut b = np_sum(&im).atan2(np_sum(&re)) / (2.0 * std::f64::consts::PI);
    for _ in 0..2 {
        let mut num = Vec::new();
        let mut den = Vec::new();
        for (k, &i) in gi.iter().enumerate() {
            let v = mw[i] - b;
            let off = v - v.round_ties_even();
            if off.abs() <= frac_max {
                num.push(w[k] * off);
                den.push(w[k]);
            }
        }
        if !num.is_empty() {
            b += np_sum(&num) / np_sum(&den);
        }
    }
    let mut nw = vec![0i64; n];
    let mut off = vec![f64::NAN; n];
    let mut fixed = vec![false; n];
    for i in 0..n {
        if mw[i].is_finite() {
            let nn = (mw[i] - b).round_ties_even();
            nw[i] = nn as i64;
            off[i] = mw[i] - b - nn;
            fixed[i] = good[i] && off[i].abs() <= frac_max;
        }
    }
    let n_out: Vec<i64> = (0..n).map(|i| if fixed[i] { nw[i] } else { 0 }).collect();
    (n_out, fixed, b, off)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spd(n: usize, seed: u64) -> DMat {
        // a deterministic SPD matrix: A Aᵀ + n I, scaled
        let mut s = seed;
        let mut rnd = || {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let a: Vec<f64> = (0..n * n).map(|_| rnd()).collect();
        let mut q = DMat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                let v: f64 = (0..n).map(|k| a[i * n + k] * a[j * n + k]).sum::<f64>()
                    + if i == j { n as f64 } else { 0.0 };
                *q.at_mut(i, j) = v * 0.02;
            }
        }
        q
    }

    #[test]
    fn ldl_reproduces_q_and_z_is_unimodular() {
        let q = spd(5, 7);
        let (l, d) = ldl(&q).unwrap();
        for i in 0..5 {
            for j in 0..5 {
                let v: f64 = (0..5).map(|k| l.at(k, i) * d[k] * l.at(k, j)).sum();
                assert!((v - q.at(i, j)).abs() < 1e-12);
            }
        }
        let (mut l2, mut d2) = ldl(&q).unwrap();
        let z = reduction(&mut l2, &mut d2);
        for v in &z.data {
            assert_eq!(v.fract(), 0.0);
        }
        let a: Vec<f64> = (0..5).map(|i| 3.3 + i as f64 * 0.9).collect();
        let (f, s) = lambda_fix(&a, &q, 2).unwrap();
        assert_eq!(f.len(), 2);
        assert!(s[0] <= s[1]);
        let r = success_rate(&q).unwrap();
        assert!(r > 0.0 && r <= 1.0);
        let fr = partial_fix(&a, &q, &|_| 3.0, 2, 0.0).unwrap();
        assert!(fr.tried >= 1);
        let pr = par_fix(&a, &q, 0.5, &|_| 1.0, 2).unwrap();
        assert!(pr.ok);
    }

    #[test]
    fn wide_lane_bias_and_integers() {
        let mw = [3.1, -2.9, 10.12, f64::NAN, 5.6];
        let sig = [0.05, 0.05, 0.05, 0.05, 0.05];
        let (n, f, b, off) = wide_lane(&mw, &sig, 0.25, 0.15);
        assert!((b - 0.1).abs() < 0.03, "{b}");
        assert_eq!(&n[..3], &[3, -3, 10]);
        assert!(!f[3] && !f[4] && off[3].is_nan());
    }
}
