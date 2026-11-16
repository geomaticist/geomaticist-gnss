//! Integer fixing and the checks: wide lanes from the Melbourne-Wübbena means, the narrow lane by LAMBDA
//! on double differences (partial fixing, group drop, conditional least squares), the ionosphere-weighted
//! L1 / L2 solution and its fix, the iono-free / iono-weighted and with / without GLONASS cross-checks,
//! the code-only trust test and the realistic σ.
//!
//! Every refusal rule and the `realistic_cov` coefficients are fixed values, kept exactly.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use geomaticist_gnss_core::Sat;
use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::text::round_to;
use serde_json::{Map, Value};

use crate::ambiguity::{self, FixResult};
use crate::arr::{B2, F2, I2};
use crate::dump::{jf, jfl, jfx, jround, jround0, jsl};
use crate::linalg::{DMat, dot, mat3_vec, norm3, rotate_cov3, solve3, sub3};
use crate::models::{ecef_to_llh, enu_matrix};
use crate::nputil::{argmax, mean, median, np_interp, np_sum, var_ddof};
use crate::options::{Options, StaticResult, fix_group};
use crate::slips::SlipEvents;
use crate::solve::{
    Block, Ctx, EngineError, Est, ROVER, RoverBlock, Rows, Scratch, assemble_blocks, assemble_rows,
    block_in, code_check_columns, code_sat_sigma, quality_cov, rounded4, rover_block_in,
    sat_report_entry,
};
use crate::stages::{
    NarrowLaneStage, SigmaStage, Stages, ValidationStage, WideLaneStage, WideLanes,
};

/// A fixed solution: the full parameter vector, its covariance (σ0²-scaled), σ0 and the constraints.
#[derive(Debug, Clone)]
pub struct Sol {
    /// Parameters.
    pub x: Vec<f64>,
    /// Covariance.
    pub cov: DMat,
    /// σ0.
    pub s0: f64,
    /// The fixed combinations `C x = c` (rows).
    pub c_mat: Vec<Vec<f64>>,
    /// Their values.
    pub c_vec: Vec<f64>,
}

/// What `_fix_ambiguities` returns.
#[derive(Debug, Clone)]
pub struct FixOut {
    /// The fixed solution when one passed.
    pub sol: Option<Sol>,
    /// The fix info dict.
    pub info: Map<String, Value>,
    /// The wide lanes.
    pub wl: Option<WideLanes>,
    /// Columns kept float (GEOs, code-biased).
    pub keep_float: Vec<bool>,
}

fn nl_pair(info: &Map<String, Value>) -> [usize; 2] {
    info.get("nl")
        .and_then(Value::as_array)
        .map(|a| {
            [
                a[0].as_u64().unwrap_or(0) as usize,
                a[1].as_u64().unwrap_or(0) as usize,
            ]
        })
        .unwrap_or([0, 0])
}

fn pick(info: &Map<String, Value>, keys: &[&str]) -> Value {
    let mut m = Map::new();
    for k in keys {
        m.insert(
            (*k).to_owned(),
            info.get(*k).cloned().unwrap_or(Value::Null),
        );
    }
    Value::Object(m)
}

fn why_of(info: &Map<String, Value>) -> String {
    info.get("why")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_owned()
}

fn status_of(info: &Map<String, Value>) -> String {
    info.get("status")
        .and_then(Value::as_str)
        .unwrap_or("")
        .to_owned()
}

fn ratio_or0(info: &Map<String, Value>) -> f64 {
    info.get("ratio").and_then(Value::as_f64).unwrap_or(0.0)
}

/// Each SD arc's wide-lane integer from its Melbourne-Wübbena mean.
pub fn wide_lanes(
    opts: &Options,
    arc: &I2,
    n_amb: usize,
    use_: &B2,
    mw: &F2,
    t: &[f64],
    grp_col: &[String],
    links: &[(usize, usize)],
    keep_float: &[bool],
    float_arcs: Option<&[bool]>,
) -> Option<WideLanes> {
    let (ne, ns) = (arc.ne, arc.ns);
    let mut root: Vec<usize> = (0..n_amb).collect();
    fn find(root: &mut [usize], mut a: usize) -> usize {
        while root[a] != a {
            root[a] = root[root[a]];
            a = root[a];
        }
        a
    }
    for &(a_, b_) in links {
        let (ra, rb) = (find(&mut root, a_), find(&mut root, b_));
        if ra != rb {
            root[ra.max(rb)] = ra.min(rb);
        }
    }
    let grp: Vec<usize> = (0..n_amb).map(|a| find(&mut root, a)).collect();
    let mut any = false;
    let mut cnt_a = vec![0.0; n_amb];
    let mut cnt: Vec<f64> = vec![0.0; n_amb];
    let mut s1 = vec![0.0; n_amb];
    let mut s2 = vec![0.0; n_amb];
    let mut t_lo = vec![f64::INFINITY; n_amb];
    let mut t_hi = vec![f64::NEG_INFINITY; n_amb];
    let mut g_lo = vec![f64::INFINITY; n_amb];
    let mut g_hi = vec![f64::NEG_INFINITY; n_amb];
    let mut jj = vec![-1i64; n_amb];
    for e in 0..ne {
        for j in 0..ns {
            let a = arc.at(e, j);
            if use_.at(e, j) && a >= 0 && mw.at(e, j).is_finite() {
                any = true;
                let a = a as usize;
                let g = grp[a];
                let m = mw.at(e, j);
                cnt_a[a] += 1.0;
                cnt[g] += 1.0;
                s1[g] += m;
                s2[g] += m * m;
                t_lo[a] = t_lo[a].min(t[e]);
                t_hi[a] = t_hi[a].max(t[e]);
                g_lo[g] = g_lo[g].min(t[e]);
                g_hi[g] = g_hi[g].max(t[e]);
                jj[a] = j as i64;
            }
        }
    }
    if !any {
        return None;
    }
    let mut mean_g = vec![0.0; n_amb];
    let mut sd = vec![0.0; n_amb];
    let mut span = vec![0.0; n_amb];
    let mut g_span = vec![0.0; n_amb];
    for a in 0..n_amb {
        let c = cnt[a].max(1.0);
        mean_g[a] = s1[a] / c;
        sd[a] = ((s2[a] / c - mean_g[a] * mean_g[a]).max(0.0) * cnt[a] / (cnt[a] - 1.0).max(1.0))
            .sqrt();
        span[a] = if cnt_a[a] > 0.0 {
            t_hi[a] - t_lo[a]
        } else {
            0.0
        };
        g_span[a] = if cnt[a] > 0.0 { g_hi[a] - g_lo[a] } else { 0.0 };
    }
    let sig: Vec<f64> = (0..n_amb)
        .map(|a| {
            let n_eff = (g_span[a] / opts.wl_corr_s + 1.0).clamp(1.0, cnt[a].max(1.0));
            sd[a].max(0.1) / n_eff.sqrt()
        })
        .collect();
    let sys_a: Vec<String> = jj
        .iter()
        .map(|&j| {
            if j >= 0 {
                grp_col[j as usize].clone()
            } else {
                String::new()
            }
        })
        .collect();
    let mut have_g: Vec<bool> = (0..n_amb)
        .map(|a| cnt[a] >= 2.0 && grp[a] == a && jj[a] >= 0)
        .collect();
    for a in 0..n_amb {
        if keep_float[jj[a].max(0) as usize] {
            have_g[a] = false;
        }
    }
    if let Some(fa) = float_arcs
        && fa.iter().any(|&f| f)
    {
        let mut gsum = vec![0.0; n_amb];
        for a in 0..n_amb {
            gsum[grp[a]] += f64::from(u8::from(fa[a]));
        }
        for a in 0..n_amb {
            if gsum[a] != 0.0 {
                have_g[a] = false;
            }
        }
    }
    let mut nw_g = vec![0i64; n_amb];
    let mut ok_g = vec![false; n_amb];
    let mut cand_g = vec![false; n_amb];
    let mut bias: BTreeMap<String, f64> = BTreeMap::new();
    let mut bias_full = vec![0.0; n_amb];
    let (mut n_fix, mut n_cand) = (0usize, 0usize);
    let groups: BTreeSet<String> = (0..n_amb)
        .filter(|&a| have_g[a])
        .map(|a| sys_a[a].clone())
        .collect();
    for s in groups {
        let first = s.chars().next().unwrap_or(' ');
        if !opts.fix_systems.contains(first) {
            continue;
        }
        let sel: Vec<usize> = (0..n_amb).filter(|&a| have_g[a] && sys_a[a] == s).collect();
        let mws: Vec<f64> = sel.iter().map(|&a| mean_g[a]).collect();
        let sgs: Vec<f64> = sel.iter().map(|&a| sig[a]).collect();
        let (n_, f_, b_, _off) = ambiguity::wide_lane(&mws, &sgs, opts.wl_frac, opts.wl_sigma);
        for (k, &a) in sel.iter().enumerate() {
            nw_g[a] = n_[k];
            ok_g[a] = f_[k];
            cand_g[a] = true;
        }
        bias.insert(s.clone(), round_to(b_, 3));
        for a in 0..n_amb {
            if sys_a[a] == s {
                bias_full[a] = b_;
            }
        }
        n_cand += sel.len();
        n_fix += f_.iter().filter(|&&f| f).count();
    }
    let groups_n: BTreeSet<usize> = (0..n_amb)
        .filter(|&a| cnt_a[a] > 0.0)
        .map(|a| grp[a])
        .collect();
    Some(WideLanes {
        nw: (0..n_amb).map(|a| nw_g[grp[a]]).collect(),
        ok: (0..n_amb).map(|a| ok_g[grp[a]] && cnt_a[a] > 0.0).collect(),
        span,
        cnt: cnt_a.clone(),
        jj,
        sys: sys_a,
        bias,
        n: [n_fix, n_cand],
        groups: groups_n.len(),
        mw: (0..n_amb).map(|a| mean_g[grp[a]] - bias_full[a]).collect(),
        sig: (0..n_amb).map(|a| sig[grp[a]]).collect(),
        cand: (0..n_amb)
            .map(|a| cand_g[grp[a]] && cnt_a[a] > 0.0)
            .collect(),
    })
}

/// The ratio a set of `n` fixed DDs (out of `n_cand` candidates) must pass.
#[must_use]
pub fn nl_threshold(opts: &Options, n: usize, n_cand: usize) -> f64 {
    let xs: Vec<f64> = opts.nl_ratio_curve.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = opts.nl_ratio_curve.iter().map(|p| p.1).collect();
    let r = np_interp(n as f64, &xs, &ys);
    opts.nl_ratio
        .max(r * (n_cand as f64 / n.max(1) as f64).powf(opts.nl_small_scale))
}

/// The candidate description handed to the narrow-lane search.
#[derive(Debug, Clone)]
pub struct Cands {
    /// Parameter index of each candidate.
    pub pidx: Vec<usize>,
    /// Wavelength (m per cycle).
    pub lam: Vec<f64>,
    /// Offset removed (m).
    pub off: Vec<f64>,
    /// Group of each.
    pub sys: Vec<String>,
    /// Epoch count (reference choice).
    pub cnt: Vec<f64>,
    /// Names.
    pub names: Vec<String>,
    /// Share of `nl_sigma`.
    pub w_c: Vec<f64>,
}

impl Cands {
    fn select(&self, keep: &[bool]) -> Cands {
        let f = |v: &Vec<usize>| {
            v.iter()
                .zip(keep)
                .filter(|(_, k)| **k)
                .map(|(x, _)| *x)
                .collect()
        };
        let g = |v: &Vec<f64>| {
            v.iter()
                .zip(keep)
                .filter(|(_, k)| **k)
                .map(|(x, _)| *x)
                .collect()
        };
        let h = |v: &Vec<String>| {
            v.iter()
                .zip(keep)
                .filter(|(_, k)| **k)
                .map(|(x, _)| x.clone())
                .collect()
        };
        Cands {
            pidx: f(&self.pidx),
            lam: g(&self.lam),
            off: g(&self.off),
            sys: h(&self.sys),
            cnt: g(&self.cnt),
            names: h(&self.names),
            w_c: g(&self.w_c),
        }
    }
}

/// Conditional least squares: `x̌ = x̂ + Q Cᵀ S⁻¹ (c − C x̂)`, `Q̌ = Q − Q Cᵀ S⁻¹ C Q`, and `vᵀ gain`.
pub fn conditional(
    x: &[f64],
    q: &DMat,
    c_mat: &[Vec<f64>],
    c_vec: &[f64],
) -> Option<(Vec<f64>, DMat, f64)> {
    let p = x.len();
    let (x_fix, qc, s_mat, dv) = conditional_core(x, q, c_mat, c_vec)?;
    // Q_fix = Q - QC S⁻¹ QCᵀ
    let cols: Vec<Vec<f64>> = (0..p).map(|i| qc[i].clone()).collect(); // rows of QC
    let sinv_qct = s_mat.solve_many(&cols)?; // sinv_qct[i] = S⁻¹ QC[i, :]ᵀ  (nf)
    let mut q_fix = q.clone();
    for i in 0..p {
        for j in 0..p {
            let s = dot(&qc[i], &sinv_qct[j]);
            *q_fix.at_mut(i, j) -= s;
        }
    }
    Some((x_fix, q_fix, dv))
}

/// `conditional` without the fixed covariance: `(x̌, vᵀPv increase)` — the same numbers as
/// `conditional`'s, for the callers that never read `Q_fix`.
pub fn conditional_x(
    x: &[f64],
    q: &DMat,
    c_mat: &[Vec<f64>],
    c_vec: &[f64],
) -> Option<(Vec<f64>, f64)> {
    let (x_fix, _, _, dv) = conditional_core(x, q, c_mat, c_vec)?;
    Some((x_fix, dv))
}

/// The shared part: `x̌`, `Q Cᵀ`, `S = C Q Cᵀ` and the vᵀPv increase.
#[allow(clippy::type_complexity)]
fn conditional_core(
    x: &[f64],
    q: &DMat,
    c_mat: &[Vec<f64>],
    c_vec: &[f64],
) -> Option<(Vec<f64>, Vec<Vec<f64>>, DMat, f64)> {
    let p = x.len();
    let nf = c_vec.len();
    // QC = Q Cᵀ  (P × nf), using the sparsity of C
    let mut qc = vec![vec![0.0; nf]; p];
    for k in 0..nf {
        let nz: Vec<(usize, f64)> = c_mat[k]
            .iter()
            .enumerate()
            .filter(|(_, v)| **v != 0.0)
            .map(|(j, v)| (j, *v))
            .collect();
        for i in 0..p {
            let mut s = 0.0;
            for &(j, cv) in &nz {
                s += q.at(i, j) * cv;
            }
            qc[i][k] = s;
        }
    }
    let mut s_mat = DMat::zeros(nf);
    for a in 0..nf {
        for b in 0..nf {
            let mut s = 0.0;
            for (j, cv) in c_mat[a].iter().enumerate() {
                if *cv != 0.0 {
                    s += cv * qc[j][b];
                }
            }
            *s_mat.at_mut(a, b) = s;
        }
    }
    let v: Vec<f64> = (0..nf).map(|k| c_vec[k] - dot(&c_mat[k], x)).collect();
    let gain = s_mat.solve(&v)?;
    let mut x_fix = x.to_vec();
    for i in 0..p {
        x_fix[i] += dot(&qc[i], &gain);
    }
    Some((x_fix, qc, s_mat, dot(&v, &gain)))
}

/// Integer DD ambiguities on a float solution.
pub fn dd_fix_one(
    opts: &Options,
    x: &[f64],
    q: &DMat,
    s0: f64,
    vtpv: f64,
    nrow: usize,
    n_sys: usize,
    c: &Cands,
    info: &mut Map<String, Value>,
) -> Option<Sol> {
    let p = x.len();
    let nc = c.pidx.len();
    let beta: Vec<f64> = (0..nc)
        .map(|i| (x[c.pidx[i]] - c.off[i]) / c.lam[i])
        .collect();
    let mut qb = DMat::zeros(nc);
    for i in 0..nc {
        for j in 0..nc {
            *qb.at_mut(i, j) =
                q.at(c.pidx[i], c.pidx[j]) / (c.lam[i] * c.lam[j]) * s0 * s0 * opts.fix_qscale;
        }
        *qb.at_mut(i, i) += opts.nl_floor * opts.nl_floor;
    }
    let mut rows: Vec<(usize, usize)> = Vec::new();
    let mut refs: Vec<usize> = Vec::new();
    let groups: BTreeSet<String> = c.sys.iter().cloned().collect();
    for s in groups {
        let loc: Vec<usize> = (0..nc).filter(|&i| c.sys[i] == s).collect();
        if loc.len() < 2 {
            continue;
        }
        let cnts: Vec<f64> = loc.iter().map(|&i| c.cnt[i]).collect();
        let r = loc[argmax(&cnts)];
        for &qi in &loc {
            if qi == r {
                continue;
            }
            let var = qb.at(qi, qi) + qb.at(r, r) - 2.0 * qb.at(qi, r);
            let lim = opts.nl_sigma * c.w_c[qi].min(c.w_c[r]);
            if var <= lim * lim {
                rows.push((qi, r));
            } else if var <= opts.nl_sigma * opts.nl_sigma {
                let e = info
                    .entry("low_out")
                    .or_insert_with(|| Value::Array(Vec::new()));
                if let Value::Array(a) = e {
                    a.push(Value::String(c.names[qi].clone()));
                }
            }
        }
        refs.push(r);
    }
    info.insert(
        "nl".into(),
        Value::Array(vec![Value::from(0), Value::from(rows.len())]),
    );
    if rows.len() < opts.nl_min {
        info.insert(
            "why".into(),
            Value::String(format!("{} narrow-lane candidates", rows.len())),
        );
        return None;
    }
    let nr = rows.len();
    let d: Vec<f64> = rows.iter().map(|&(qi, r)| beta[qi] - beta[r]).collect();
    let mut qd = DMat::zeros(nr);
    for a in 0..nr {
        let (qa, ra) = rows[a];
        for b in 0..nr {
            let (qb_, rb) = rows[b];
            *qd.at_mut(a, b) = qb.at(qa, qb_) - qb.at(qa, rb) - qb.at(ra, qb_) + qb.at(ra, rb);
        }
    }
    let n_cand = nr;
    let thr = |n: usize| nl_threshold(opts, n, n_cand);
    let fr: FixResult;
    if nr <= opts.nl_greedy {
        match ambiguity::partial_fix(&d, &qd, &thr, opts.nl_min, opts.nl_offset) {
            Ok(f) => fr = f,
            Err(_) => {
                info.insert("why".into(), Value::String("covariance".into()));
                return None;
            }
        }
        info.insert("method".into(), Value::String("partial".into()));
    } else {
        match ambiguity::par_fix(&d, &qd, opts.nl_p0, &thr, opts.nl_min) {
            Ok(f) => fr = f,
            Err(_) => {
                info.insert("why".into(), Value::String("covariance".into()));
                return None;
            }
        }
        info.insert("method".into(), Value::String("par".into()));
    }
    if nr <= 40 {
        match ambiguity::success_rate(&qd) {
            Ok(sr) => info.insert("success".into(), jround(sr, 4)),
            Err(_) => {
                info.insert("why".into(), Value::String("covariance".into()));
                return None;
            }
        };
    } else {
        info.insert("success".into(), Value::Null);
    }
    info.insert(
        "ratio".into(),
        if fr.ratio.is_finite() {
            jround(fr.ratio, 2)
        } else {
            Value::Null
        },
    );
    info.insert(
        "dropped".into(),
        jsl(&fr
            .dropped
            .iter()
            .map(|&i| c.names[rows[i].0].clone())
            .collect::<Vec<_>>()),
    );
    if opts.debug {
        let mut li = Map::new();
        li.insert(
            "names".into(),
            jsl(&rows
                .iter()
                .map(|r| c.names[r.0].clone())
                .collect::<Vec<_>>()),
        );
        li.insert(
            "refs".into(),
            jsl(&rows
                .iter()
                .map(|r| c.names[r.1].clone())
                .collect::<Vec<_>>()),
        );
        li.insert("d".into(), jfl(&d));
        li.insert("Qd".into(), crate::dump::jmat(&qd));
        li.insert(
            "keep".into(),
            Value::Array(fr.keep.iter().map(|&k| Value::from(k)).collect()),
        );
        li.insert(
            "fixed".into(),
            if fr.ok && fr.c_mat.is_none() {
                Value::Array(fr.fixed.iter().map(|v| Value::from(*v as i64)).collect())
            } else {
                Value::Null
            },
        );
        li.insert("ok".into(), Value::Bool(fr.ok));
        info.insert("lambda_in".into(), Value::Object(li));
        if fr.c_mat.is_none() {
            let dd: Vec<Value> = fr
                .keep
                .iter()
                .enumerate()
                .map(|(i, &k)| {
                    Value::Array(vec![
                        Value::String(c.names[rows[k].0].clone()),
                        jround(d[k], 2),
                        if fr.ok {
                            Value::from(fr.fixed[i] as i64)
                        } else {
                            Value::Null
                        },
                        jround(qd.at(k, k).sqrt(), 3),
                    ])
                })
                .collect();
            info.insert("dd".into(), Value::Array(dd));
        }
    }
    if !fr.ok {
        let why = if fr.ratio != 0.0 && !fr.ratio.is_nan() {
            if fr.ratio.is_infinite() {
                "ratio inf".to_owned()
            } else {
                format!("ratio {:.2}", fr.ratio)
            }
        } else {
            "no subset".to_owned()
        };
        info.insert("why".into(), Value::String(why));
        return None;
    }
    // the fixed combinations Cd d = cd  ->  C x = c
    let (cd_mat, cd_vec): (Vec<Vec<f64>>, Vec<f64>) = match (&fr.c_mat, &fr.c_vec) {
        (Some(cm), Some(cv)) => (cm.clone(), cv.clone()),
        _ => {
            let mut cm = vec![vec![0.0; nr]; fr.keep.len()];
            for (i, &k) in fr.keep.iter().enumerate() {
                cm[i][k] = 1.0;
            }
            (cm, fr.fixed.clone())
        }
    };
    let nf = cd_vec.len();
    // CB = Cd D / lam  (nf × nc)
    let mut cb = vec![vec![0.0; nc]; nf];
    for k in 0..nf {
        for (r, &(qi, ri)) in rows.iter().enumerate() {
            let v = cd_mat[k][r];
            if v != 0.0 {
                cb[k][qi] += v;
                cb[k][ri] -= v;
            }
        }
        for i in 0..nc {
            cb[k][i] /= c.lam[i];
        }
    }
    let mut c_full = vec![vec![0.0; p]; nf];
    let mut c_vec = vec![0.0; nf];
    for k in 0..nf {
        for i in 0..nc {
            c_full[k][c.pidx[i]] = cb[k][i];
        }
        c_vec[k] = cd_vec[k] + dot(&cb[k], &c.off);
    }
    let (x_fix, q_fix, dvtpv) = conditional(x, q, &c_full, &c_vec)?;
    let dof = (nrow as i64 - p as i64 + n_sys as i64).max(1);
    let s0_fix = ((vtpv + dvtpv) / (dof as f64 + nf as f64)).sqrt();
    let dpos = [x_fix[0] - x[0], x_fix[1] - x[1], x_fix[2] - x[2]];
    let mut q3 = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            q3[i][j] = q.at(i, j) * s0 * s0;
        }
    }
    let shift_z =
        solve3(&q3, dpos).map_or(f64::NAN, |g| crate::linalg::dot3(dpos, g).max(0.0).sqrt());
    info.insert(
        "nl".into(),
        Value::Array(vec![Value::from(nf), Value::from(nr)]),
    );
    info.insert("s0_float".into(), jround(s0, 3));
    info.insert("s0_fixed".into(), jround(s0_fix, 3));
    info.insert(
        "refs".into(),
        jsl(&refs.iter().map(|&r| c.names[r].clone()).collect::<Vec<_>>()),
    );
    info.insert("chi".into(), jround(fr.s.0 / nf as f64, 2));
    info.insert("p".into(), jfx(fr.p));
    info.insert("shift_mm".into(), jround(1000.0 * norm3(dpos), 1));
    info.insert("shift_z".into(), jround(shift_z, 1));
    if s0_fix > opts.fix_s0_max * s0 {
        info.insert(
            "why".into(),
            Value::String(format!("σ0 {s0:.2} → {s0_fix:.2}")),
        );
        info.insert("status".into(), Value::String("refused".into()));
        return None;
    }
    info.insert("status".into(), Value::String("fixed".into()));
    Some(Sol {
        x: x_fix,
        cov: q_fix.scaled(s0_fix * s0_fix),
        s0: s0_fix,
        c_mat: c_full,
        c_vec,
    })
}

/// `_dd_fix_one` on every candidate, then without each group in turn when nothing passes.
pub fn dd_fix(
    opts: &Options,
    x: &[f64],
    q: &DMat,
    s0: f64,
    vtpv: f64,
    nrow: usize,
    n_sys: usize,
    c: &Cands,
    info: &mut Map<String, Value>,
) -> Option<Sol> {
    let res = dd_fix_one(opts, x, q, s0, vtpv, nrow, n_sys, c, info);
    let base_g: Vec<String> = c
        .sys
        .iter()
        .map(|g| g.split('~').next().unwrap_or("").to_owned())
        .collect();
    let groups: BTreeSet<String> = base_g.iter().cloned().collect();
    let why = why_of(info);
    let retry = res.is_none()
        && opts.nl_drop_group
        && groups.len() >= 2
        && (why == "no subset" || why.starts_with("ratio") || status_of(info) == "refused");
    if !retry {
        return res;
    }
    let mut best: Option<(String, Map<String, Value>, Sol)> = None;
    for g in groups {
        let keep: Vec<bool> = base_g.iter().map(|b| *b != g).collect();
        let sub = c.select(&keep);
        let mut inf = Map::new();
        inf.insert("status".into(), Value::String("float".into()));
        inf.insert("why".into(), Value::String(String::new()));
        if let Some(r) = dd_fix_one(opts, x, q, s0, vtpv, nrow, n_sys, &sub, &mut inf) {
            let better = match &best {
                None => true,
                Some((_, bi, _)) => {
                    (nl_pair(&inf)[0], ratio_or0(&inf)) > (nl_pair(bi)[0], ratio_or0(bi))
                }
            };
            if better {
                best = Some((g, inf, r));
            }
        }
    }
    let (g, inf, r) = best?;
    let first = pick(info, &["why", "ratio", "nl", "status"]);
    let first: Map<String, Value> = first
        .as_object()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter(|(k, _)| info.contains_key(k))
        .collect();
    for (k, v) in inf {
        info.insert(k, v);
    }
    info.insert("why".into(), Value::String(String::new()));
    info.insert("dropped_group".into(), Value::String(g));
    info.insert("all_groups".into(), Value::Object(first));
    Some(r)
}

/// GLONASS: the common constant of the SD ambiguities from the code.
#[must_use]
pub fn glo_common(x_amb: &[f64], x_code: &[f64]) -> Option<f64> {
    let d: Vec<f64> = x_amb
        .iter()
        .zip(x_code)
        .map(|(a, b)| a - b)
        .filter(|v| v.is_finite())
        .collect();
    (d.len() >= 2).then(|| median(&d))
}

/// Wide lane, then narrow lane, on the float solution.
pub fn fix_ambiguities(
    ctx: &Ctx,
    opts: &Options,
    x: &[f64],
    q: &DMat,
    s0: f64,
    est: &Est,
    arc: &I2,
    n_amb: usize,
    use_: &B2,
    mw: &F2,
    links: &[(usize, usize)],
    biased: &BTreeMap<usize, f64>,
    pc: Option<&[f64]>,
    float_arcs: Option<&[bool]>,
    w_arc: Option<&[f64]>,
) -> FixOut {
    let mut info = Map::new();
    info.insert("status".into(), Value::String("float".into()));
    info.insert(
        "wl".into(),
        Value::Array(vec![Value::from(0), Value::from(0)]),
    );
    info.insert(
        "nl".into(),
        Value::Array(vec![Value::from(0), Value::from(0)]),
    );
    info.insert("ratio".into(), Value::Null);
    info.insert("why".into(), Value::String(String::new()));
    let mut keep_float: Vec<bool> = if opts.fix_geo {
        vec![false; ctx.ns]
    } else {
        ctx.cols
            .iter()
            .map(|c| crate::options::bds_geo(*c))
            .collect()
    };
    for &j in biased.keys() {
        keep_float[j] = true;
    }
    let grp: Vec<String> = ctx.cols.iter().map(|c| fix_group(*c, opts)).collect();
    let Some(w) = wide_lanes(
        opts,
        arc,
        n_amb,
        use_,
        mw,
        &ctx.t,
        &grp,
        links,
        &keep_float,
        float_arcs,
    ) else {
        info.insert("why".into(), Value::String("no MW".into()));
        return FixOut {
            sol: None,
            info,
            wl: None,
            keep_float,
        };
    };
    info.insert(
        "wl".into(),
        Value::Array(vec![Value::from(w.n[0]), Value::from(w.n[1])]),
    );
    info.insert(
        "wl_bias".into(),
        Value::Object(w.bias.iter().map(|(k, v)| (k.clone(), jf(*v))).collect()),
    );
    let mut lam_n = vec![0.0; n_amb];
    let mut k_w = vec![0.0; n_amb];
    for a in 0..n_amb {
        if w.jj[a] >= 0 {
            let j = w.jj[a] as usize;
            let (f1, f2) = (CLIGHT / ctx.lam1[j], CLIGHT / ctx.lam2[j]);
            lam_n[a] = CLIGHT / (f1 + f2);
            k_w[a] = ctx.bif[j] * ctx.lam2[j];
        }
    }
    let mut cand: Vec<usize> = (0..n_amb)
        .filter(|&a| w.ok[a] && w.span[a] >= opts.fix_min_s)
        .collect();
    let mut off: Vec<f64> = cand.iter().map(|&a| k_w[a] * w.nw[a] as f64).collect();
    let is_r: Vec<bool> = cand
        .iter()
        .map(|&a| ctx.cols[w.jj[a] as usize].sys.letter() == 'R')
        .collect();
    if is_r.iter().any(|&r| r) {
        let c_g = pc.and_then(|pc| {
            let xa: Vec<f64> = cand
                .iter()
                .zip(&is_r)
                .filter(|(_, r)| **r)
                .map(|(&a, _)| x[ctx.i_amb + a])
                .collect();
            let xc: Vec<f64> = cand
                .iter()
                .zip(&is_r)
                .filter(|(_, r)| **r)
                .map(|(&a, _)| pc[a])
                .collect();
            glo_common(&xa, &xc)
        });
        match c_g {
            None => {
                cand = cand
                    .iter()
                    .zip(&is_r)
                    .filter(|(_, r)| !**r)
                    .map(|(&a, _)| a)
                    .collect();
                off = off
                    .iter()
                    .zip(&is_r)
                    .filter(|(_, r)| !**r)
                    .map(|(&o, _)| o)
                    .collect();
                info.insert(
                    "glo".into(),
                    Value::String(
                        "no GLONASS arcs with code to place their common constant".into(),
                    ),
                );
            }
            Some(cg) => {
                for (k, r) in is_r.iter().enumerate() {
                    if *r {
                        off[k] += cg;
                    }
                }
                info.insert("glo_common_m".into(), jround(cg, 3));
            }
        }
    }
    let cands = Cands {
        pidx: cand.iter().map(|&a| ctx.i_amb + a).collect(),
        lam: cand.iter().map(|&a| lam_n[a]).collect(),
        off,
        sys: cand.iter().map(|&a| w.sys[a].clone()).collect(),
        cnt: cand.iter().map(|&a| w.cnt[a]).collect(),
        names: cand
            .iter()
            .map(|&a| ctx.cols[w.jj[a] as usize].to_string())
            .collect(),
        w_c: cand
            .iter()
            .map(|&a| w_arc.map_or(1.0, |wa| wa[a]))
            .collect(),
    };
    let sol = dd_fix(
        opts,
        x,
        q,
        s0,
        est.acc.vtpv,
        est.acc.nrow,
        ctx.n_sys,
        &cands,
        &mut info,
    );
    FixOut {
        sol,
        info,
        wl: Some(w),
        keep_float,
    }
}

/// Each arc's share of `nl_sigma` from the rover's elevation at its first epoch.
#[must_use]
pub fn arc_weight(opts: &Options, arc: &I2, n_amb: usize, use_: &B2, el: &F2) -> Option<Vec<f64>> {
    if opts.ar_el <= 0.0 || n_amb == 0 {
        return None;
    }
    let mut w = vec![1.0; n_amb];
    let mut seen = vec![false; n_amb];
    let s_ar = opts.ar_el.to_radians().sin();
    for e in 0..arc.ne {
        for j in 0..arc.ns {
            let a = arc.at(e, j);
            if use_.at(e, j) && a >= 0 && !seen[a as usize] {
                seen[a as usize] = true;
                w[a as usize] = (el.at(e, j).max(1e-3).sin() / s_ar).min(1.0);
            }
        }
    }
    Some(w)
}

/// σ_I (m) of the SD slant ionosphere at each node.
#[must_use]
pub fn iono_prior(opts: &Options, length_m: f64, spread_m: f64) -> f64 {
    (opts.iono_sigma.powi(2)
        + (opts.iono_ppm * 1e-6 * length_m).powi(2)
        + (opts.iono_data * spread_m).powi(2))
    .sqrt()
}

/// Whether two fixed solutions disagree: `(apart, [dH, dU] mm)`.
#[must_use]
pub fn apart(
    opts: &Options,
    ax: &[f64],
    acov: &DMat,
    bx: &[f64],
    bcov: &DMat,
    r: &[[f64; 3]; 3],
) -> (bool, Value) {
    let de = mat3_vec(r, [ax[0] - bx[0], ax[1] - bx[1], ax[2] - bx[2]]);
    let mut cs = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cs[i][j] = acov.at(i, j) + bcov.at(i, j);
        }
    }
    let cc = rotate_cov3(r, &cs);
    let (dh, du) = (de[0].hypot(de[1]), de[2].abs());
    let sh = (cc[0][0] + cc[1][1]).max(0.0).sqrt();
    let su = cc[2][2].max(0.0).sqrt();
    let ap = dh > opts.iono_agree.max(opts.iono_agree_k * sh)
        || du > (2.0 * opts.iono_agree).max(opts.iono_agree_k * su);
    (
        ap,
        Value::Array(vec![jround(1000.0 * dh, 1), jround(1000.0 * de[2], 1)]),
    )
}

/// Validation of an ionosphere-weighted fix. Returns the fix or `None` when refused.
pub fn iono_check(
    opts: &Options,
    iw: Sol,
    info: &mut Map<String, Value>,
    x: &[f64],
    q: &DMat,
    s0: f64,
    fx: &FixOut,
    r: &[[f64; 3]; 3],
) -> Option<Sol> {
    let nl = nl_pair(info);
    if nl[0] < opts.iono_nl_min || (nl[0] as f64) < opts.iono_share_min * nl[1] as f64 {
        info.insert("status".into(), Value::String("refused".into()));
        info.insert(
            "why".into(),
            Value::String(format!(
                "{} of {} DDs fixed (iono-weighted: at least {} and {:.0} %)",
                nl[0],
                nl[1],
                opts.iono_nl_min,
                100.0 * opts.iono_share_min
            )),
        );
        return None;
    }
    let mut c3 = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c3[i][j] = iw.cov.at(i, j);
        }
    }
    let cf = rotate_cov3(r, &c3);
    let shf = (cf[0][0] + cf[1][1]).max(0.0).sqrt();
    info.insert(
        "sigma_fixed_mm".into(),
        Value::Array(vec![
            jround(1000.0 * shf, 1),
            jround(1000.0 * cf[2][2].max(0.0).sqrt(), 1),
        ]),
    );
    if opts.iono_sh_max > 0.0 && shf > opts.iono_sh_max {
        info.insert("status".into(), Value::String("refused".into()));
        info.insert(
            "why".into(),
            Value::String(format!(
                "fixed σH {:.0} mm: the ionosphere holds it, not the integers",
                1000.0 * shf
            )),
        );
        return None;
    }
    let d = [iw.x[0] - x[0], iw.x[1] - x[1], iw.x[2] - x[2]];
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = q.at(i, j) * s0 * s0;
        }
    }
    let z = solve3(&c, d).map_or(f64::NAN, |g| crate::linalg::dot3(d, g).max(0.0).sqrt());
    info.insert("z_float".into(), jround(z, 2));
    if opts.iono_z_max > 0.0 && z > opts.iono_z_max {
        info.insert("status".into(), Value::String("refused".into()));
        info.insert(
            "why".into(),
            Value::String(format!("{z:.1} σ from the iono-free float")),
        );
        return None;
    }
    if let Some(fs) = &fx.sol {
        let (ap, mm) = apart(opts, &iw.x, &iw.cov, &fs.x, &fs.cov, r);
        info.insert("agree_mm".into(), mm.clone());
        if ap {
            info.insert("conflict".into(), mm);
            if ["ionofree", "float", "arbiter"].contains(&opts.iono_conflict.as_str()) {
                info.insert("status".into(), Value::String("refused".into()));
                info.insert(
                    "why".into(),
                    Value::String("disagrees with the iono-free fix".into()),
                );
                return None;
            }
        }
    }
    Some(iw)
}

/// The third arbiter between an iono-free and an iono-weighted fix that disagree.
fn arbiter(
    opts: &Options,
    iw: &Sol,
    iw_info: &Map<String, Value>,
    fx: &Sol,
    fx_info: &Map<String, Value>,
    w3: Option<&Sol>,
    w3_info: &Map<String, Value>,
    r: &[[f64; 3]; 3],
) -> (String, Map<String, Value>) {
    let growth = |inf: &Map<String, Value>| -> f64 {
        let a = inf.get("s0_fixed").and_then(Value::as_f64);
        let b = inf.get("s0_float").and_then(Value::as_f64);
        match (a, b) {
            (Some(a), Some(b)) if a != 0.0 && b != 0.0 => a / b,
            _ => f64::INFINITY,
        }
    };
    let (gw, gi) = (growth(iw_info), growth(fx_info));
    let mut out = Map::new();
    out.insert(
        "s0_growth".into(),
        Value::Array(vec![jround(gw, 3), jround(gi, 3)]),
    );
    out.insert(
        "w3".into(),
        pick(w3_info, &["status", "nl", "ratio", "why"]),
    );
    let Some(w3) = w3 else {
        out.insert(
            "why".into(),
            Value::String("the wider-prior search did not fix".into()),
        );
        return ("float".into(), out);
    };
    let (apart_w, mm_w) = apart(opts, &w3.x, &w3.cov, &iw.x, &iw.cov, r);
    out.insert("w3_to_weighted_mm".into(), mm_w);
    let (apart_i, mm_i) = apart(opts, &w3.x, &w3.cov, &fx.x, &fx.cov, r);
    out.insert("w3_to_ionofree_mm".into(), mm_i);
    if !apart_w && apart_i && gw <= gi {
        return ("weighted".into(), out);
    }
    if !apart_i && apart_w && gi <= gw {
        return ("ionofree".into(), out);
    }
    out.insert(
        "why".into(),
        Value::String("no side both held by the wider prior and with the smaller σ0 growth".into()),
    );
    ("float".into(), out)
}

/// The ionosphere-weighted L1 / L2 solution and its integer fix.
struct IonoWeighted<'a> {
    est: &'a Est,
    w: &'a WideLanes,
    keep_float: &'a [bool],
    fixable: Option<&'a str>,
    prior_k: f64,
    arc: &'a I2,
    n_amb: usize,
    w_arc: Option<&'a [f64]>,
    pc_arc: &'a [f64],
}

/// The per-cell ionosphere node layout of the weighted solution.
struct IonoLayout {
    lw: B2,
    lf: B2,
    arcs_w: Vec<usize>,
    kn: Vec<Vec<usize>>,
    uu: F2,
    wk: Vec<i64>,
    n_i: usize,
    i_i: usize,
    i_w: usize,
    p2: usize,
    gam: Vec<f64>,
    r21: Vec<f64>,
    l1o: F2,
    l2o: F2,
    var1: F2,
    nw0: Vec<f64>,
    sw: Vec<f64>,
    amb_c: Vec<f64>,
    w_rw: Vec<(usize, usize, f64)>,
}

/// What the builds of one design share: the normal matrix without the rover's rows and columns and
/// without the priors. Every element of it is the same sum of the same numbers in every build while the
/// design holds — the troposphere, ambiguity, ionosphere and wide-lane design entries, the weights and the
/// groups do not depend on the parameters; only the rover's `−e` columns and the residuals do — and the
/// rover's corrections (`Ctx::geom` re-evaluates them when the rover moves over 1 cm; `at_b` names that
/// state) are those it was built with.
struct IonoCache {
    n_fixed: DMat,
    at_b: [f64; 3],
}

impl IonoWeighted<'_> {
    fn build(
        &self,
        ctx: &mut Ctx,
        lay: &IonoLayout,
        xx: &[f64],
        s_i: f64,
        cache: &mut Option<IonoCache>,
    ) -> (DMat, Vec<f64>, crate::stages::Acc) {
        let t_b = Instant::now();
        let (r, reused) = self.build_inner(ctx, lay, xx, s_i, cache);
        if std::env::var_os("GEOM_ENGINE_TIMING").is_some() {
            eprintln!(
                "timing:   iono build P2={} {:.3} s{}",
                lay.p2,
                t_b.elapsed().as_secs_f64(),
                if reused { " (rover rows only)" } else { "" }
            );
        }
        r
    }

    /// The normal equations at `xx`: a full assembly when the design changed (or first), else the cached
    /// matrix plus the rover rows / columns, the right-hand side and vᵀPv — the same numbers either way.
    #[allow(clippy::too_many_lines)]
    fn build_inner(
        &self,
        ctx: &mut Ctx,
        lay: &IonoLayout,
        xx: &[f64],
        s_i: f64,
        cache: &mut Option<IonoCache>,
    ) -> ((DMat, Vec<f64>, crate::stages::Acc), bool) {
        let (rho_b, eb) = ctx.geom(xx);
        let p2 = lay.p2;
        let mut bv = vec![0.0; p2];
        let mut acc = crate::stages::Acc::default();
        let mif = ctx.model_sd(xx, &rho_b);
        let (ne, ns) = (ctx.ne, ctx.ns);
        let mut res_l1 = F2::nans(ne, ns);
        let mut res_l2 = F2::nans(ne, ns);
        let mut res_if = F2::nans(ne, ns);
        let mut res_cd = F2::nans(ne, ns);
        if ns > 0 {
            // per epoch (independent cells; a few hundred thousand of them: not worth a parallel section
            // in a shared pool)
            let ctx_r: &Ctx = ctx;
            res_l1
                .data
                .chunks_mut(ns)
                .zip(res_l2.data.chunks_mut(ns))
                .zip(res_if.data.chunks_mut(ns))
                .zip(res_cd.data.chunks_mut(ns))
                .enumerate()
                .for_each(|(e, (((r1, r2), rif), rcd))| {
                    let ctx = ctx_r;
                    for j in 0..ns {
                        let agf = ctx.cb.antgf.at(e, j) - ctx.ca.antgf.at(e, j);
                        let m1 = mif.at(e, j) - ctx.bif[j] * agf;
                        let m2 = m1 - agf;
                        let a = self.arc.at(e, j);
                        let amb = if a >= 0 {
                            xx[ctx.i_amb + a as usize]
                        } else {
                            0.0
                        };
                        let (nwv, iv) = if lay.lw.at(e, j) {
                            let k = lay.kn[e][j];
                            (
                                xx[lay.i_w + lay.wk[a as usize] as usize],
                                xx[lay.i_i + k] * (1.0 - lay.uu.at(e, j))
                                    + xx[lay.i_i + (k + 1).min(lay.n_i - 1)] * lay.uu.at(e, j),
                            )
                        } else {
                            (0.0, 0.0)
                        };
                        r1[j] = lay.l1o.at(e, j) - m1 - amb + iv;
                        r2[j] = lay.l2o.at(e, j) - m2 - lay.r21[j] * amb
                            + ctx.lam2[j] * nwv
                            + lay.gam[j] * iv;
                        rif[j] = ctx.obs_ph.at(e, j) - mif.at(e, j) - amb;
                        rcd[j] = ctx.obs_cd.at(e, j) - mif.at(e, j);
                    }
                });
        }
        let ntr_ = 3 + ctx.ntr;
        let wf = &self.est.wf_ph;
        let wf_c = &self.est.wf_cd;
        let ctx_r: &Ctx = ctx;
        let form = ctx_r.opts.compat.vtpv;
        // one group's design for one kind into the scratch (`asd`, `var`, `y`): its parameters, or
        // `None` with under two rows
        let design =
            |e: usize, cj_all: &[usize], kind: usize, sc: &mut Scratch| -> Option<Vec<usize>> {
                let ctx = ctx_r;
                let Scratch {
                    sel, asd, var, y, ..
                } = sc;
                sel.clear();
                match kind {
                    0 | 1 => sel.extend(cj_all.iter().copied().filter(|&j| lay.lw.at(e, j))),
                    2 => sel.extend(cj_all.iter().copied().filter(|&j| lay.lf.at(e, j))),
                    _ => {
                        ctx.var_cd.as_ref()?;
                        sel.extend(cj_all.iter().copied().filter(|&j| wf_c.at(e, j) > 0.0));
                    }
                }
                let nr = sel.len();
                if nr < 2 {
                    return None;
                }
                let amb_cols = kind != 3;
                let ion_cols = kind <= 1;
                let nw_cols = kind == 1;
                let u = ntr_
                    + if amb_cols { nr } else { 0 }
                    + if ion_cols { 2 * nr } else { 0 }
                    + if nw_cols { nr } else { 0 };
                asd.reset(nr, u);
                let mut g: Vec<usize> = vec![0, 1, 2];
                let mut c0 = 3;
                for (r, &j) in sel.iter().enumerate() {
                    let ee = eb.at(e, j);
                    asd[(r, 0)] = -ee[0];
                    asd[(r, 1)] = -ee[1];
                    asd[(r, 2)] = -ee[2];
                }
                if ctx.nn > 0 {
                    let k = ctx.kk[e];
                    for (r, &j) in sel.iter().enumerate() {
                        asd[(r, c0)] = ctx.cb.map.at(e, j) * (1.0 - ctx.ww[e]);
                        asd[(r, c0 + 1)] = ctx.cb.map.at(e, j) * ctx.ww[e];
                    }
                    g.push(ctx.i_tb + k);
                    g.push(ctx.i_tb + k + 1);
                    c0 += 2;
                    if ctx.n_ta > 0 {
                        for (r, &j) in sel.iter().enumerate() {
                            asd[(r, c0)] = -ctx.ca.map.at(e, j) * (1.0 - ctx.ww[e]);
                            asd[(r, c0 + 1)] = -ctx.ca.map.at(e, j) * ctx.ww[e];
                        }
                        g.push(ctx.i_ta + k);
                        g.push(ctx.i_ta + k + 1);
                        c0 += 2;
                    }
                }
                if amb_cols {
                    for (r, &j) in sel.iter().enumerate() {
                        asd[(r, c0 + r)] = if kind == 1 { lay.r21[j] } else { 1.0 };
                        g.push(ctx.i_amb + self.arc.at(e, j) as usize);
                    }
                    c0 += nr;
                }
                if ion_cols {
                    for (r, &j) in sel.iter().enumerate() {
                        let g_ = if kind == 1 { lay.gam[j] } else { 1.0 };
                        asd[(r, c0 + 2 * r)] = -g_ * (1.0 - lay.uu.at(e, j));
                        asd[(r, c0 + 2 * r + 1)] = -g_ * lay.uu.at(e, j);
                        g.push(lay.i_i + lay.kn[e][j]);
                        g.push(lay.i_i + lay.kn[e][j] + 1);
                    }
                    c0 += 2 * nr;
                }
                if nw_cols {
                    for (r, &j) in sel.iter().enumerate() {
                        asd[(r, c0 + r)] = -ctx.lam2[j];
                        g.push(lay.i_w + lay.wk[self.arc.at(e, j) as usize] as usize);
                    }
                }
                var.clear();
                var.extend(sel.iter().map(|&j| match kind {
                    3 => ctx.var_cd.as_ref().map_or(f64::NAN, |v| v.at(e, j)) / wf_c.at(e, j),
                    2 => ctx.var_ph.at(e, j) / wf.at(e, j),
                    _ => lay.var1.at(e, j) / wf.at(e, j),
                }));
                let res = match kind {
                    0 => &res_l1,
                    1 => &res_l2,
                    2 => &res_if,
                    _ => &res_cd,
                };
                y.clear();
                y.extend(sel.iter().map(|&j| res.at(e, j)));
                Some(g)
            };
        let reused = cache.as_ref().is_some_and(|c| c.at_b == ctx_r.at_b);
        let mut n = match cache.as_ref().filter(|c| c.at_b == ctx_r.at_b) {
            Some(c) => c.n_fixed.clone(),
            None => DMat::zeros(p2),
        };
        if reused {
            let group =
                |e: usize, cj_all: &[usize], sc: &mut Scratch, out: &mut Vec<RoverBlock>| {
                    for kind in 0..4 {
                        if let Some(g) = design(e, cj_all, kind, sc) {
                            out.push(rover_block_in(sc, g, form));
                        }
                    }
                };
            assemble_rows(&ctx_r.groups_all, &mut n, &mut bv, &mut acc, group);
        } else {
            let group = |e: usize, cj_all: &[usize], sc: &mut Scratch, out: &mut Vec<Block>| {
                for kind in 0..4 {
                    if let Some(g) = design(e, cj_all, kind, sc) {
                        out.push(block_in(sc, g, None, form));
                    }
                }
            };
            assemble_blocks(&ctx_r.groups_all, &mut n, &mut bv, &mut acc, group);
            let mut n_fixed = n.clone();
            for i in 0..ROVER {
                for j in 0..p2 {
                    *n_fixed.at_mut(i, j) = 0.0;
                    *n_fixed.at_mut(j, i) = 0.0;
                }
            }
            *cache = Some(IonoCache {
                n_fixed,
                at_b: ctx_r.at_b,
            });
        }
        let w = 1.0 / (ctx.opts.amb_prior * ctx.opts.amb_prior);
        for a in 0..self.n_amb {
            let i = ctx.i_amb + a;
            n.add(i, i, w);
            bv[i] += w * (lay.amb_c[a] - xx[i]);
        }
        let wi = 1.0 / (s_i * s_i);
        for i in lay.i_i..lay.i_w {
            n.add(i, i, wi);
            bv[i] -= wi * xx[i];
        }
        for (k, i) in (lay.i_w..p2).enumerate() {
            let ww_ = 1.0 / (lay.sw[k] * lay.sw[k]);
            n.add(i, i, ww_);
            bv[i] += ww_ * (lay.nw0[k] - xx[i]);
        }
        ctx.trop_priors(&mut n, &mut bv, xx);
        for &(a_, b_, w_) in &lay.w_rw {
            Ctx::tie(&mut n, &mut bv, xx, a_, b_, w_);
        }
        ((n, bv, acc), reused)
    }

    #[allow(clippy::too_many_lines)]
    fn run(&self, ctx: &mut Ctx) -> Result<(Option<Sol>, Map<String, Value>), EngineError> {
        let opts = ctx.opts.clone();
        let mut info = Map::new();
        info.insert("status".into(), Value::String("float".into()));
        info.insert("mode".into(), Value::String("iono-weighted".into()));
        info.insert(
            "nl".into(),
            Value::Array(vec![Value::from(0), Value::from(0)]),
        );
        info.insert("ratio".into(), Value::Null);
        info.insert("why".into(), Value::String(String::new()));
        let w = self.w;
        let x0 = &self.est.x;
        let (ne, ns, n_amb) = (ctx.ne, ctx.ns, self.n_amb);
        let fixsys = self.fixable.unwrap_or(&opts.fix_systems);
        let fixsys_col: Vec<bool> = ctx.sysc.iter().map(|s| fixsys.contains(*s)).collect();
        let okc: Vec<bool> = (0..n_amb)
            .map(|a| {
                w.jj[a] >= 0
                    && fixsys_col[w.jj[a].max(0) as usize]
                    && !self.keep_float[w.jj[a].max(0) as usize]
            })
            .collect();
        let a_fix: Vec<bool> = (0..n_amb).map(|a| w.ok[a] && okc[a]).collect();
        let a_flt: Vec<bool> = (0..n_amb)
            .map(|a| {
                opts.iono_wl_float
                    && w.cand[a]
                    && !w.ok[a]
                    && okc[a]
                    && w.mw[a].is_finite()
                    && w.sig[a].is_finite()
            })
            .collect();
        let a_l12: Vec<bool> = (0..n_amb).map(|a| a_fix[a] || a_flt[a]).collect();
        let mut lw = B2::falses(ne, ns);
        let mut lf = B2::falses(ne, ns);
        for e in 0..ne {
            for j in 0..ns {
                let u = self.est.use_.at(e, j) && self.est.wf_ph.at(e, j) > 0.0;
                let a = self.arc.at(e, j);
                let l = u && a >= 0 && a_l12[a as usize];
                lw.set(e, j, l);
                lf.set(e, j, u && !l);
            }
        }
        let mut arcs_set: BTreeSet<usize> = BTreeSet::new();
        let mut tlo = vec![f64::INFINITY; n_amb];
        let mut thi = vec![f64::NEG_INFINITY; n_amb];
        for e in 0..ne {
            for j in 0..ns {
                if lw.at(e, j) {
                    let a = self.arc.at(e, j) as usize;
                    arcs_set.insert(a);
                    tlo[a] = tlo[a].min(ctx.t[e]);
                    thi[a] = thi[a].max(ctx.t[e]);
                }
            }
        }
        let arcs_w: Vec<usize> = arcs_set.into_iter().collect();
        if arcs_w.len() <= opts.nl_min {
            info.insert(
                "why".into(),
                Value::String(format!("{} wide-lane arcs", arcs_w.len())),
            );
            return Ok((None, info));
        }
        let spans: Vec<f64> = arcs_w.iter().map(|&a| thi[a] - tlo[a]).collect();
        let total = np_sum(&spans);
        let stp = opts
            .iono_step_s
            .max(total / (opts.iono_max_nodes as f64 - arcs_w.len() as f64).max(1.0));
        let mut kcount = vec![0usize; n_amb];
        for &a in &arcs_w {
            kcount[a] = ((((thi[a] - tlo[a]) / stp) - 1e-9).ceil() as i64).max(1) as usize;
        }
        let mut n0 = vec![0usize; n_amb];
        let mut n_i = 0usize;
        for &a in &arcs_w {
            n0[a] = n_i;
            n_i += kcount[a] + 1;
        }
        let i_i = ctx.i_amb + n_amb;
        let i_w = i_i + n_i;
        let mut wk = vec![-1i64; n_amb];
        for (k, &a) in arcs_w.iter().enumerate() {
            wk[a] = k as i64;
        }
        let p2 = i_w + arcs_w.len();
        let mut kn = vec![vec![0usize; ns]; ne];
        let mut uu = F2::zeros(ne, ns);
        for e in 0..ne {
            for j in 0..ns {
                if lw.at(e, j) {
                    let a = self.arc.at(e, j) as usize;
                    let u_ = (ctx.t[e] - tlo[a]) / stp;
                    let k_ = (u_.floor() as i64).min(kcount[a] as i64 - 1).max(0) as usize;
                    kn[e][j] = n0[a] + k_;
                    uu.set(e, j, u_ - k_ as f64);
                }
            }
        }
        let lb = norm3(sub3([x0[0], x0[1], x0[2]], ctx.base_xyz));
        let gam: Vec<f64> = (0..ns)
            .map(|j| (ctx.lam2[j] / ctx.lam1[j]).powi(2))
            .collect();
        let r21: Vec<f64> = (0..ns).map(|j| ctx.lam2[j] / ctx.lam1[j]).collect();
        let gfo = ctx.sdgf.zip_map(&ctx.corr_gf, |a, b| a + b);
        let mut sv = Vec::new();
        for &a in &arcs_w {
            let j = w.jj[a] as usize;
            let g_: Vec<f64> = (0..ne)
                .filter(|&e| lw.at(e, j) && self.arc.at(e, j) == a as i64)
                .map(|e| gfo.at(e, j) / (gam[j] - 1.0))
                .collect();
            if g_.len() >= 20 {
                let d: Vec<f64> = g_.windows(2).map(|p| p[1] - p[0]).collect();
                sv.push((var_ddof(&g_, 0) - 0.5 * var_ddof(&d, 0)).max(0.0));
            }
        }
        let s_raw = if sv.is_empty() {
            0.0
        } else {
            median(&sv).sqrt()
        };
        let mut s_i = iono_prior(&opts, lb, s_raw) * self.prior_k;
        let l1o = ctx.obs_ph.zip_map(&gfo, |p, g| p - 0.0 * g); // placeholder shape
        let mut l1o = l1o;
        for e in 0..ne {
            for j in 0..ns {
                l1o.set(e, j, ctx.obs_ph.at(e, j) - ctx.bif[j] * gfo.at(e, j));
            }
        }
        let l2o = l1o.zip_map(&gfo, |a, b| a - b);
        let var1 = ctx.var_ph.map(|v| v / 9.0);
        let nw0: Vec<f64> = arcs_w
            .iter()
            .map(|&a| if a_fix[a] { w.nw[a] as f64 } else { w.mw[a] })
            .collect();
        let sw: Vec<f64> = arcs_w
            .iter()
            .map(|&a| {
                if a_fix[a] {
                    1e-3
                } else {
                    w.sig[a].max(0.02) * opts.iono_wl_scale
                }
            })
            .collect();
        let mut x = vec![0.0; p2];
        x[..ctx.i_amb + n_amb].copy_from_slice(&x0[..ctx.i_amb + n_amb]);
        x[i_w..].copy_from_slice(&nw0);
        for &a in &arcs_w {
            let j = w.jj[a] as usize;
            let lam_n = ctx.lam1[j] * ctx.lam2[j] / (ctx.lam1[j] + ctx.lam2[j]);
            x[ctx.i_amb + a] = ctx.lam1[j]
                * (x0[ctx.i_amb + a] - ctx.bif[j] * ctx.lam2[j] * x[i_w + wk[a] as usize])
                / lam_n;
        }
        let amb_c = x[ctx.i_amb..ctx.i_amb + n_amb].to_vec();
        let mut w_rw = Vec::new();
        if opts.iono_rw > 0.0 {
            for &a in &arcs_w {
                for k in 0..kcount[a] {
                    w_rw.push((
                        i_i + n0[a] + k,
                        i_i + n0[a] + k + 1,
                        1.0 / (opts.iono_rw * opts.iono_rw * stp),
                    ));
                }
            }
        }
        let lay = IonoLayout {
            lw,
            lf,
            arcs_w: arcs_w.clone(),
            kn,
            uu,
            wk,
            n_i,
            i_i,
            i_w,
            p2,
            gam,
            r21,
            l1o,
            l2o,
            var1,
            nw0,
            sw,
            amb_c,
            w_rw,
        };
        let s_i0 = s_i;
        let mut cache: Option<IonoCache> = None;
        let mut acc2 = crate::stages::Acc::default();
        let mut q2 = DMat::zeros(p2);
        for rnd in 0..=opts.iono_vce_rounds {
            for it in 0..4 {
                let (nm, bm, _) = self.build(ctx, &lay, &x, s_i, &mut cache);
                let dx = nm
                    .solve(&bm)
                    .ok_or_else(|| EngineError("singular iono-weighted normal equations".into()))?;
                for i in 0..p2 {
                    x[i] += dx[i];
                }
                if dx[..3].iter().map(|v| v.abs()).fold(0.0, f64::max) < 1e-5 && it > 0 {
                    break;
                }
            }
            let (n2, _bm, a2) = self.build(ctx, &lay, &x, s_i, &mut cache);
            acc2 = a2;
            let t_i = Instant::now();
            q2 = n2
                .inverse()
                .ok_or_else(|| EngineError("singular iono-weighted normal matrix".into()))?;
            if std::env::var_os("GEOM_ENGINE_TIMING").is_some() {
                eprintln!(
                    "timing:   iono inverse {:.3} s",
                    t_i.elapsed().as_secs_f64()
                );
            }
            let red: f64 = (i_i..i_w).map(|i| 1.0 - q2.at(i, i) / (s_i * s_i)).sum();
            let s_hat = if red > 0.0 {
                ((i_i..i_w).map(|i| x[i] * x[i]).sum::<f64>() / red).sqrt()
            } else {
                0.0
            };
            info.insert(
                "iono_vce".into(),
                Value::Array(vec![jround(1000.0 * s_hat, 1), jround(red, 1)]),
            );
            if rnd == opts.iono_vce_rounds
                || red < opts.iono_vce_red
                || s_hat <= opts.iono_vce_k * s_i
            {
                break;
            }
            s_i = s_hat * opts.iono_vce_k;
        }
        let n_sys_ = ctx.n_sys;
        let dof = (acc2.nrow as i64 - p2 as i64 + n_sys_ as i64).max(1);
        let mut s0w = (acc2.vtpv / dof as f64).sqrt();
        if s_i > s_i0 {
            info.insert("iono_sigma0_mm".into(), jround(1000.0 * s_i0, 1));
        }
        info.insert("iono_sigma_mm".into(), jround(1000.0 * s_i, 1));
        info.insert("iono_spread_mm".into(), jround(1000.0 * s_raw, 1));
        info.insert("iono_step_s".into(), jround0(stp));
        info.insert("s0_float".into(), jround(s0w, 3));
        info.insert("weighted_arcs".into(), Value::from(arcs_w.len()));
        info.insert(
            "wl_float".into(),
            Value::from(arcs_w.iter().filter(|&&a| a_flt[a]).count()),
        );
        let cand: Vec<usize> = arcs_w
            .iter()
            .copied()
            .filter(|&a| w.span[a] >= opts.fix_min_s)
            .collect();
        let cw: Vec<usize> = cand.iter().copied().filter(|&a| a_flt[a]).collect();
        let mut cands = Cands {
            pidx: cand.iter().map(|&a| ctx.i_amb + a).collect(),
            lam: cand.iter().map(|&a| ctx.lam1[w.jj[a] as usize]).collect(),
            off: vec![0.0; cand.len()],
            sys: cand.iter().map(|&a| w.sys[a].clone()).collect(),
            cnt: cand.iter().map(|&a| w.cnt[a]).collect(),
            names: cand
                .iter()
                .map(|&a| ctx.cols[w.jj[a] as usize].to_string())
                .collect(),
            w_c: cand
                .iter()
                .map(|&a| self.w_arc.map_or(1.0, |wa| wa[a]))
                .collect(),
        };
        let cw_groups: BTreeSet<String> = cw.iter().map(|&a| w.sys[a].clone()).collect();
        for g in cw_groups {
            let mut mem: Vec<usize> = cw.iter().copied().filter(|&a| w.sys[a] == g).collect();
            let fx_: Vec<usize> = cand
                .iter()
                .copied()
                .filter(|&a| a_fix[a] && w.sys[a] == g)
                .collect();
            if !fx_.is_empty() {
                let cnts: Vec<f64> = fx_.iter().map(|&a| w.cnt[a]).collect();
                mem.push(fx_[argmax(&cnts)]);
            }
            if mem.len() < 2 {
                continue;
            }
            for a in mem {
                cands.pidx.push(i_w + lay.wk[a] as usize);
                cands.lam.push(1.0);
                cands.off.push(0.0);
                cands.sys.push(format!("{g}~w"));
                cands.cnt.push(w.cnt[a] + if a_fix[a] { 1e9 } else { 0.0 });
                cands
                    .names
                    .push(format!("{}~w", ctx.cols[w.jj[a] as usize]));
                cands.w_c.push(self.w_arc.map_or(1.0, |wa| wa[a]));
            }
        }
        let glo_a: Vec<usize> = (0..cand.len())
            .filter(|&i| ctx.cols[w.jj[cand[i]] as usize].sys.letter() == 'R')
            .collect();
        let (c_bif, c_lam1, c_lam2, c_iamb) = (
            ctx.bif.clone(),
            ctx.lam1.clone(),
            ctx.lam2.clone(),
            ctx.i_amb,
        );
        let search = |xx: &[f64],
                      qq: &DMat,
                      s0_: f64,
                      acc_: &crate::stages::Acc,
                      inf: &mut Map<String, Value>|
         -> Option<Sol> {
            let mut c = cands.clone();
            if !glo_a.is_empty() {
                let ag: Vec<usize> = glo_a.iter().map(|&i| cand[i]).collect();
                let a1_code: Vec<f64> = ag
                    .iter()
                    .map(|&a| {
                        let j = w.jj[a] as usize;
                        (self.pc_arc[a] - c_bif[j] * c_lam2[j] * xx[i_w + lay.wk[a] as usize])
                            * (1.0 + c_lam1[j] / c_lam2[j])
                    })
                    .collect();
                let xa: Vec<f64> = ag.iter().map(|&a| xx[c_iamb + a]).collect();
                match glo_common(&xa, &a1_code) {
                    None => {
                        let keep: Vec<bool> = c
                            .sys
                            .iter()
                            .map(|g| g.split('~').next().unwrap_or("") != "R")
                            .collect();
                        c = c.select(&keep);
                    }
                    Some(c1) => {
                        for &i in &glo_a {
                            c.off[i] = c1;
                        }
                        inf.insert("glo_common_m".into(), jround(c1, 3));
                    }
                }
            }
            dd_fix(&opts, xx, qq, s0_, acc_.vtpv, acc_.nrow, n_sys_, &c, inf)
        };
        let sys_l: Vec<String> = ctx.cols.iter().map(|c| fix_group(*c, &opts)).collect();
        let post_sigma = |xf: &[f64], ctx: &Ctx| -> Option<f64> {
            let mut ss = 0.0;
            let mut n_ = 0usize;
            for (e, cj_all) in &ctx.groups_all {
                let gs: BTreeSet<String> = cj_all.iter().map(|&j| sys_l[j].clone()).collect();
                for g in gs {
                    let sel: Vec<usize> = cj_all
                        .iter()
                        .copied()
                        .filter(|&j| lay.lw.at(*e, j) && sys_l[j] == g)
                        .collect();
                    if sel.len() >= 2 {
                        let v: Vec<f64> = sel
                            .iter()
                            .map(|&j| {
                                let k = lay.kn[*e][j];
                                xf[i_i + k] * (1.0 - lay.uu.at(*e, j))
                                    + xf[i_i + (k + 1).min(n_i - 1)] * lay.uu.at(*e, j)
                            })
                            .collect();
                        let m = mean(&v);
                        let sq: Vec<f64> = v.iter().map(|x| (x - m) * (x - m)).collect();
                        ss += np_sum(&sq) * sel.len() as f64 / (sel.len() - 1) as f64;
                        n_ += sel.len();
                    }
                }
            }
            (n_ >= 20).then(|| (ss / n_ as f64).sqrt())
        };
        let mut res = search(&x, &q2, s0w, &acc2, &mut info);
        info.insert(
            "wl".into(),
            Value::Array(vec![Value::from(w.n[0]), Value::from(w.n[1])]),
        );
        let mut rounds = 0usize;
        while let Some(r) = &res {
            let (xf, _dv) = conditional_x(&x, &q2, &r.c_mat, &r.c_vec)
                .ok_or_else(|| EngineError("conditional LS failed".into()))?;
            let s_post = post_sigma(&xf, ctx);
            info.insert(
                "iono_post_mm".into(),
                s_post.map_or(Value::Null, |s| jround(1000.0 * s, 1)),
            );
            let Some(sp) = s_post else { break };
            if sp <= opts.iono_refix_k * s_i || rounds >= opts.iono_refix_rounds {
                break;
            }
            rounds += 1;
            // at_prior(s_post, x)
            s_i = sp;
            for _ in 0..2 {
                let (nm, bm, _) = self.build(ctx, &lay, &x, s_i, &mut cache);
                let dx = nm
                    .solve(&bm)
                    .ok_or_else(|| EngineError("singular iono-weighted normal equations".into()))?;
                for i in 0..p2 {
                    x[i] += dx[i];
                }
            }
            let (nm, _bm, a_) = self.build(ctx, &lay, &x, s_i, &mut cache);
            q2 = nm
                .inverse()
                .ok_or_else(|| EngineError("singular iono-weighted normal matrix".into()))?;
            s0w = (a_.vtpv / (a_.nrow as i64 - p2 as i64 + n_sys_ as i64).max(1) as f64).sqrt();
            acc2 = a_;
            let mut inf2 = Map::new();
            inf2.insert("status".into(), Value::String("float".into()));
            inf2.insert("mode".into(), Value::String("iono-weighted".into()));
            inf2.insert(
                "nl".into(),
                Value::Array(vec![Value::from(0), Value::from(0)]),
            );
            inf2.insert("ratio".into(), Value::Null);
            inf2.insert("why".into(), Value::String(String::new()));
            res = search(&x, &q2, s0w, &acc2, &mut inf2);
            for k_ in [
                "status",
                "nl",
                "ratio",
                "why",
                "dropped",
                "method",
                "success",
                "s0_fixed",
                "refs",
                "chi",
                "p",
                "shift_mm",
                "shift_z",
                "dropped_group",
                "all_groups",
            ] {
                match inf2.get(k_) {
                    Some(v) => {
                        info.insert(k_.into(), v.clone());
                    }
                    None => {
                        info.remove(k_);
                    }
                }
            }
            info.insert("s0_float".into(), jround(s0w, 3));
            info.insert("iono_refix".into(), Value::from(rounds));
            info.insert("iono_sigma_mm".into(), jround(1000.0 * s_i, 1));
            if res.is_none() {
                let w0 = why_of(&info);
                info.insert(
                    "why".into(),
                    Value::String(format!("{w0} (ionosphere over its prior)")),
                );
            }
        }
        let Some(mut r) = res else {
            return Ok((None, info));
        };
        if opts.iono_alt {
            let mut xx = x.clone();
            let s_save = s_i;
            let s_alt = s_i * 100.0;
            for _ in 0..2 {
                let (nm, bm, _) = self.build(ctx, &lay, &xx, s_alt, &mut cache);
                let dx = nm
                    .solve(&bm)
                    .ok_or_else(|| EngineError("singular iono-weighted normal equations".into()))?;
                for i in 0..p2 {
                    xx[i] += dx[i];
                }
            }
            let (nm, _bm, _a) = self.build(ctx, &lay, &xx, s_alt, &mut cache);
            let qq = nm
                .inverse()
                .ok_or_else(|| EngineError("singular iono-weighted normal matrix".into()))?;
            let (xf, _) = conditional_x(&xx, &qq, &r.c_mat, &r.c_vec)
                .ok_or_else(|| EngineError("conditional LS failed".into()))?;
            info.insert(
                "vector_loose".into(),
                jfl(&sub3([xf[0], xf[1], xf[2]], ctx.base_xyz)),
            );
            s_i = s_save;
        }
        if (opts.iono_final_k - 1.0).abs() > 1e-9 {
            let s_fin = s_i * opts.iono_final_k;
            let mut xx = x.clone();
            for _ in 0..2 {
                let (nm, bm, _) = self.build(ctx, &lay, &xx, s_fin, &mut cache);
                let dx = nm
                    .solve(&bm)
                    .ok_or_else(|| EngineError("singular iono-weighted normal equations".into()))?;
                for i in 0..p2 {
                    xx[i] += dx[i];
                }
            }
            let (nm, _bm, a_) = self.build(ctx, &lay, &xx, s_fin, &mut cache);
            let qq = nm
                .inverse()
                .ok_or_else(|| EngineError("singular iono-weighted normal matrix".into()))?;
            let (xf, qf, dv) = conditional(&xx, &qq, &r.c_mat, &r.c_vec)
                .ok_or_else(|| EngineError("conditional LS failed".into()))?;
            let s0f = ((a_.vtpv + dv)
                / ((a_.nrow as i64 - p2 as i64 + n_sys_ as i64).max(1) as f64
                    + r.c_vec.len() as f64))
                .sqrt();
            r.x = xf;
            r.cov = qf.scaled(s0f * s0f);
            r.s0 = s0f;
            info.insert("iono_final_mm".into(), jround(1000.0 * s_fin, 1));
        }
        let ix: Vec<f64> = r.x[i_i..i_w].to_vec();
        let rms = if ix.is_empty() {
            0.0
        } else {
            (np_sum(&ix.iter().map(|v| v * v).collect::<Vec<_>>()) / ix.len() as f64).sqrt()
        };
        info.insert("iono_rms_mm".into(), jround(1000.0 * rms, 1));
        let _ = &lay.arcs_w;
        Ok((Some(r), info))
    }
}

/// The second half of `solve`: fixing, checks, trust, report, σ. Returns the result.
#[allow(clippy::too_many_lines)]
pub fn finish(
    ctx: &mut Ctx,
    stages: &mut Stages,
    est: Est,
    q: DMat,
    arc: &I2,
    n_amb: usize,
    links: &[(usize, usize, f64)],
    biased: &BTreeMap<usize, f64>,
    wl_links: &[(usize, usize)],
    vec_float: [f64; 3],
    cov_float: [[f64; 3]; 3],
    glo_why: &str,
    n_split: usize,
    n_gross_code: usize,
    eph_age: &BTreeMap<Sat, f64>,
    sat_slips: &BTreeMap<Sat, (usize, usize)>,
    slip_ev: SlipEvents,
    slip_ev_gf: &BTreeMap<String, usize>,
    mut notes: Vec<String>,
    t_start: Instant,
) -> Result<StaticResult, EngineError> {
    let opts = ctx.opts.clone();
    let (ne, ns) = (ctx.ne, ctx.ns);
    let x = est.x.clone();
    let s0 = est.s0;
    let base_xyz = ctx.base_xyz;
    let t = ctx.t.clone();
    let cov_float_m = DMat {
        n: 3,
        data: cov_float.iter().flatten().copied().collect(),
    };
    let mut fix_info: Map<String, Value> = Map::new();
    fix_info.insert("status".into(), Value::String("float".into()));
    let mut x_out = x.clone();
    let mut cov_out = cov_float_m.clone();
    let mut s0_out = s0;
    // each arc's mean SD iono-free phase minus code
    let mut pc_arc = vec![f64::NAN; n_amb];
    if ctx.var_cd.is_some() && n_amb > 0 {
        let mut n_pc = vec![0usize; n_amb];
        let mut s_pc = vec![0.0; n_amb];
        for e in 0..ne {
            for j in 0..ns {
                let a = arc.at(e, j);
                if est.use_.at(e, j)
                    && est.wf_ph.at(e, j) > 0.0
                    && est.wf_cd.at(e, j) > 0.0
                    && a >= 0
                    && ctx.obs_ph.at(e, j).is_finite()
                    && ctx.obs_cd.at(e, j).is_finite()
                {
                    n_pc[a as usize] += 1;
                    s_pc[a as usize] += ctx.obs_ph.at(e, j) - ctx.obs_cd.at(e, j);
                }
            }
        }
        for a in 0..n_amb {
            if n_pc[a] >= 10 {
                pc_arc[a] = s_pc[a] / n_pc[a] as f64;
            }
        }
    }
    let use_w = est.use_.zip_map(&est.wf_ph, |u, w| u && w > 0.0);
    let w_arc = arc_weight(&opts, arc, n_amb, &use_w, &ctx.cb.el);
    let mut half_arcs: Option<Vec<bool>> = None;
    if ctx.half_sd.any() && n_amb > 0 {
        let mut ha = vec![false; n_amb];
        let mut hs: BTreeSet<Sat> = BTreeSet::new();
        for e in 0..ne {
            for j in 0..ns {
                let a = arc.at(e, j);
                if ctx.half_sd.at(e, j) && est.use_.at(e, j) && a >= 0 {
                    ha[a as usize] = true;
                }
            }
        }
        for e in 0..ne {
            for j in 0..ns {
                let a = arc.at(e, j);
                if ctx.half_sd.at(e, j) && est.use_.at(e, j) && a >= 0 && ha[a as usize] {
                    hs.insert(ctx.cols[j]);
                }
            }
        }
        if ha.iter().any(|&h| h) {
            notes.push(format!(
                "half-cycle flag (LLI 2) on {}: {} arcs kept float",
                hs.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(" "),
                ha.iter().filter(|&&h| h).count()
            ));
        }
        half_arcs = Some(ha);
    }
    let mw_all = ctx.sdmw.zip_map(&ctx.corr_mw, |a, b| a + b);
    let (lat0, lon0, _) = ecef_to_llh(base_xyz);
    let re0 = enu_matrix(lat0, lon0);
    if opts.fix && n_amb > 0 && links.is_empty() {
        let run_fix = |ctx: &Ctx, o: &Options| -> FixOut {
            fix_ambiguities(
                ctx,
                o,
                &x,
                &q,
                s0,
                &est,
                arc,
                n_amb,
                &use_w,
                &mw_all,
                wl_links,
                biased,
                Some(&pc_arc),
                half_arcs.as_deref(),
                w_arc.as_deref(),
            )
        };
        let mut fx = run_fix(ctx, &opts);
        fix_info = fx.info.clone();
        if let Some(s) = &fx.sol {
            x_out = s.x.clone();
            cov_out = s.cov.clone();
            s0_out = s.s0;
        }
        stages.wide_lane = Some(WideLaneStage {
            wl: fx.wl.clone(),
            keep_float: fx.keep_float.clone(),
            half_arcs: half_arcs.clone(),
            w_arc: w_arc.clone(),
            pc_arc: pc_arc.clone(),
        });
        let mut nl_info = fx.info.clone();
        nl_info.insert(
            "s0".into(),
            fx.sol.as_ref().map_or(Value::Null, |s| jf(s.s0)),
        );
        stages.narrow_lane = Some(NarrowLaneStage {
            x: fx.sol.as_ref().map(|s| s.x.clone()),
            cov: fx.sol.as_ref().map(|s| s.cov.clone()),
            info: nl_info,
        });
        let l_km = norm3(sub3([x[0], x[1], x[2]], base_xyz)) / 1000.0;
        let run_iw = |ctx: &mut Ctx,
                      fx: &FixOut,
                      fixable: Option<&str>,
                      prior_k: f64|
         -> Result<(Option<Sol>, Map<String, Value>), EngineError> {
            match &fx.wl {
                None => {
                    let mut info = Map::new();
                    info.insert("status".into(), Value::String("float".into()));
                    info.insert("mode".into(), Value::String("iono-weighted".into()));
                    info.insert(
                        "nl".into(),
                        Value::Array(vec![Value::from(0), Value::from(0)]),
                    );
                    info.insert("ratio".into(), Value::Null);
                    info.insert("why".into(), Value::String("no wide lanes".into()));
                    Ok((None, info))
                }
                Some(w) => IonoWeighted {
                    est: &est,
                    w,
                    keep_float: &fx.keep_float,
                    fixable,
                    prior_k,
                    arc,
                    n_amb,
                    w_arc: w_arc.as_deref(),
                    pc_arc: &pc_arc,
                }
                .run(ctx),
            }
        };
        // conflicted: (iw, fixable)
        let mut conflicted: Option<(Sol, Option<String>)> = None;
        let fs_no_r = opts.fix_systems.replace('R', "");
        let weighted_against = |ctx: &mut Ctx,
                                fx: &FixOut,
                                conflicted: &mut Option<(Sol, Option<String>)>|
         -> Result<(Option<Sol>, Map<String, Value>), EngineError> {
            *conflicted = None;
            let (iw0, mut iw_info) = run_iw(ctx, fx, None, 1.0)?;
            let mut iw = None;
            if let Some(s) = iw0 {
                iw_info.insert(
                    "vector".into(),
                    jfl(&sub3([s.x[0], s.x[1], s.x[2]], base_xyz)),
                );
                let pre = s.clone();
                iw = iono_check(&opts, s, &mut iw_info, &x, &q, s0, fx, &re0);
                if iw.is_none() && iw_info.contains_key("conflict") {
                    *conflicted = Some((pre, None));
                }
            }
            let has_r_cand = fx
                .wl
                .as_ref()
                .is_some_and(|w| (0..w.cand.len()).any(|a| w.cand[a] && w.sys[a] == "R"));
            if iw.is_none() && opts.fix_systems.contains('R') && has_r_cand {
                let (iw2_0, mut inf2) = run_iw(ctx, fx, Some(&fs_no_r), 1.0)?;
                let pre2 = iw2_0.clone();
                let mut iw2 = None;
                if let Some(s) = iw2_0 {
                    inf2.insert(
                        "vector".into(),
                        jfl(&sub3([s.x[0], s.x[1], s.x[2]], base_xyz)),
                    );
                    iw2 = iono_check(&opts, s, &mut inf2, &x, &q, s0, fx, &re0);
                }
                if iw2.is_some() || inf2.contains_key("conflict") {
                    *conflicted = None;
                    if iw2.is_none()
                        && let Some(p2) = pre2
                    {
                        *conflicted = Some((p2, Some(fs_no_r.clone())));
                    }
                    inf2.insert(
                        "glo_first".into(),
                        pick(&iw_info, &["status", "nl", "ratio", "why", "conflict"]),
                    );
                    return Ok((iw2, inf2));
                }
            }
            Ok((iw, iw_info))
        };
        if opts.iono == "weighted" || (opts.iono == "auto" && l_km <= opts.iono_km) {
            let (mut iw, mut iw_info) = weighted_against(ctx, &fx, &mut conflicted)?;
            if iw.is_none()
                && iw_info.contains_key("conflict")
                && fx.sol.is_some()
                && fx.info.contains_key("glo_common_m")
            {
                let mut o2 = opts.clone();
                o2.fix_systems = fs_no_r.clone();
                let fx2 = run_fix(ctx, &o2);
                let first = pick(&fix_info, &["status", "nl", "ratio"]);
                fx = fx2;
                fix_info = fx.info.clone();
                let mut gc = Map::new();
                gc.insert("ionofree".into(), first);
                gc.insert(
                    "weighted".into(),
                    pick(&iw_info, &["status", "nl", "ratio", "conflict"]),
                );
                fix_info.insert("glo_conflict".into(), Value::Object(gc));
                match &fx.sol {
                    Some(s) => {
                        x_out = s.x.clone();
                        cov_out = s.cov.clone();
                        s0_out = s.s0;
                    }
                    None => {
                        x_out = x.clone();
                        cov_out = cov_float_m.clone();
                        s0_out = s0;
                    }
                }
                let r = weighted_against(ctx, &fx, &mut conflicted)?;
                iw = r.0;
                iw_info = r.1;
            }
            let keys = ["status", "ratio", "nl", "why"];
            if let Some(s) = &fx.sol {
                fix_info.insert(
                    "vector".into(),
                    jfl(&sub3([s.x[0], s.x[1], s.x[2]], base_xyz)),
                );
            }
            if let Some(s) = &iw {
                iw_info.insert(
                    "ionofree".into(),
                    pick(&fix_info, &["status", "ratio", "nl", "why", "vector"]),
                );
                if let Some(gc) = fix_info.get("glo_conflict") {
                    iw_info.insert("glo_conflict".into(), gc.clone());
                }
                iw_info.insert(
                    "wl_bias".into(),
                    fix_info.get("wl_bias").cloned().unwrap_or(Value::Null),
                );
                fix_info = iw_info;
                x_out = s.x.clone();
                cov_out = s.cov.clone();
                s0_out = s.s0;
            } else {
                let weighted = pick(
                    &iw_info,
                    &[
                        "status",
                        "ratio",
                        "nl",
                        "why",
                        "z_float",
                        "conflict",
                        "sigma_fixed_mm",
                        "vector",
                    ],
                );
                let mut verdict: Option<String> = None;
                if iw_info.contains_key("conflict")
                    && opts.iono_conflict == "arbiter"
                    && let Some(fs) = &fx.sol
                    && let Some((ciw, cfix)) = conflicted.clone()
                {
                    let (w3, w3_info) = run_iw(ctx, &fx, cfix.as_deref(), opts.arb_prior_k)?;
                    let (v, mut arb) = arbiter(
                        &opts,
                        &ciw,
                        &iw_info,
                        fs,
                        &fix_info,
                        w3.as_ref(),
                        &w3_info,
                        &re0,
                    );
                    arb.insert("verdict".into(), Value::String(v.clone()));
                    verdict = Some(v.clone());
                    if v == "weighted" {
                        iw_info.insert("status".into(), Value::String("fixed".into()));
                        iw_info.insert("why".into(), Value::String(String::new()));
                        iw_info.insert("arbiter".into(), Value::Object(arb));
                        iw_info.insert(
                            "ionofree".into(),
                            pick(&fix_info, &["status", "ratio", "nl", "why", "vector"]),
                        );
                        iw_info.insert(
                            "wl_bias".into(),
                            fix_info.get("wl_bias").cloned().unwrap_or(Value::Null),
                        );
                        fix_info = iw_info.clone();
                        x_out = ciw.x.clone();
                        cov_out = ciw.cov.clone();
                        s0_out = ciw.s0;
                        iw = Some(ciw);
                    } else {
                        fix_info.insert("arbiter".into(), Value::Object(arb));
                    }
                }
                if iw.is_some() {
                } else if iw_info.contains_key("conflict")
                    && ["float", "arbiter"].contains(&opts.iono_conflict.as_str())
                    && verdict.as_deref() != Some("ionofree")
                    && fx.sol.is_some()
                {
                    let mut f2 = Map::new();
                    f2.insert("status".into(), Value::String("float".into()));
                    for k in ["nl", "wl", "ratio", "wl_bias"] {
                        f2.insert(k.into(), fix_info.get(k).cloned().unwrap_or(Value::Null));
                    }
                    f2.insert(
                        "why".into(),
                        Value::String("the iono-weighted and iono-free fixes disagree".into()),
                    );
                    f2.insert(
                        "ionofree".into(),
                        pick(&fix_info, &["status", "ratio", "nl", "why", "vector"]),
                    );
                    if let Some(a) = fix_info.get("arbiter") {
                        f2.insert("arbiter".into(), a.clone());
                    }
                    fix_info = f2;
                    x_out = x.clone();
                    cov_out = cov_float_m.clone();
                    s0_out = s0;
                }
                if iw.is_none() {
                    fix_info.insert("weighted".into(), weighted);
                }
                let _ = keys;
            }
        }
        fix_info
            .entry("mode")
            .or_insert_with(|| Value::String("ionofree".into()));
        let refs_have_r = fix_info
            .get("refs")
            .and_then(Value::as_array)
            .is_some_and(|a| {
                a.iter()
                    .any(|r| r.as_str().is_some_and(|s| s.starts_with('R')))
            });
        if opts.glo_check != "off" && status_of(&fix_info) == "fixed" && refs_have_r {
            let (alt, alt_info): (Option<Sol>, Map<String, Value>) =
                if fix_info.get("mode").and_then(Value::as_str) == Some("iono-weighted") {
                    let (a0, mut ai) = run_iw(ctx, &fx, Some(&fs_no_r), 1.0)?;
                    let a = a0.and_then(|s| iono_check(&opts, s, &mut ai, &x, &q, s0, &fx, &re0));
                    (a, ai)
                } else {
                    let mut o3 = opts.clone();
                    o3.fix_systems = fs_no_r.clone();
                    let f3 = run_fix(ctx, &o3);
                    (f3.sol, f3.info)
                };
            let mut chk = Map::new();
            chk.insert(
                "status".into(),
                Value::String(if alt.is_some() { "fixed" } else { "float" }.into()),
            );
            for k in ["nl", "ratio", "why"] {
                chk.insert(k.into(), alt_info.get(k).cloned().unwrap_or(Value::Null));
            }
            let mut verdict = "kept";
            if let Some(a) = &alt {
                let (ap, mm) = apart(&opts, &x_out, &cov_out, &a.x, &a.cov, &re0);
                chk.insert("apart_mm".into(), mm);
                if ap {
                    verdict = if opts.glo_check == "agree" {
                        "without"
                    } else {
                        "float"
                    };
                }
            } else if opts.glo_check == "confirm" {
                verdict = "float";
            }
            chk.insert("verdict".into(), Value::String(verdict.into()));
            let glo_with = pick(&fix_info, &["status", "nl", "ratio", "refs"]);
            if verdict == "without"
                && let Some(a) = &alt
            {
                let mut ai = alt_info.clone();
                ai.insert("status".into(), Value::String("fixed".into()));
                ai.insert(
                    "mode".into(),
                    fix_info.get("mode").cloned().unwrap_or(Value::Null),
                );
                let wl = match alt_info.get("wl") {
                    Some(v) if !v.is_null() => v.clone(),
                    _ => fix_info.get("wl").cloned().unwrap_or(Value::Null),
                };
                ai.insert("wl".into(), wl);
                ai.insert(
                    "wl_bias".into(),
                    fix_info.get("wl_bias").cloned().unwrap_or(Value::Null),
                );
                ai.insert("glo_with".into(), glo_with);
                fix_info = ai;
                x_out = a.x.clone();
                cov_out = a.cov.clone();
                s0_out = a.s0;
            } else if verdict == "float" {
                let mut f2 = Map::new();
                f2.insert("status".into(), Value::String("float".into()));
                for k in ["mode", "nl", "wl", "ratio", "wl_bias"] {
                    f2.insert(k.into(), fix_info.get(k).cloned().unwrap_or(Value::Null));
                }
                f2.insert(
                    "why".into(),
                    Value::String(
                        if alt.is_some() {
                            "the fixes with and without GLONASS disagree"
                        } else {
                            "the fix needs GLONASS and none without it confirms it"
                        }
                        .into(),
                    ),
                );
                f2.insert("glo_with".into(), glo_with);
                fix_info = f2;
                x_out = x.clone();
                cov_out = cov_float_m.clone();
                s0_out = s0;
            }
            fix_info.insert("glo_check".into(), Value::Object(chk));
        }
    }
    if opts.fix && ctx.sysc.contains(&'R') {
        fix_info.insert("glo".into(), Value::String(glo_why.to_owned()));
        notes.push(glo_why.to_owned());
    }
    if !biased.is_empty() {
        fix_info.insert(
            "code_bias".into(),
            Value::Object(
                biased
                    .iter()
                    .map(|(j, v)| (ctx.cols[*j].to_string(), jf(*v)))
                    .collect(),
            ),
        );
    }
    // can it be trusted?
    let mut failed = String::new();
    let (lat_, lon_, _) = ecef_to_llh(base_xyz);
    let re = enu_matrix(lat_, lon_);
    if let Some(var_cd) = ctx.var_cd.clone() {
        let (rb_, eb_) = ctx.geom(&x_out);
        let m = ctx.model_sd(&x_out, &rb_);
        let rc = ctx.obs_cd.zip_map(&m, |a, b| a - b);
        let have_cd: Vec<bool> = (0..ns)
            .map(|j| (0..ne).any(|e| est.wf_cd.at(e, j) > 0.0))
            .collect();
        let in_c = code_check_columns(&ctx.sysc, &have_cd, biased, &opts);
        let mut n3 = DMat::zeros(3);
        let mut b3 = vec![0.0; 3];
        let mut acc3 = crate::stages::Acc::default();
        let mut used_c: Vec<(usize, Vec<usize>)> = Vec::new();
        for (e, cj_all) in &ctx.groups_all {
            let sel: Vec<usize> = cj_all
                .iter()
                .copied()
                .filter(|&j| est.wf_cd.at(*e, j) > 0.0 && rc.at(*e, j).is_finite() && in_c[j])
                .collect();
            if sel.len() >= 2 {
                let asd = Rows::from_fn(sel.len(), 3, |r, c| -eb_.at(*e, sel[r])[c]);
                let var: Vec<f64> = sel
                    .iter()
                    .map(|&j| var_cd.at(*e, j) / est.wf_cd.at(*e, j))
                    .collect();
                let y: Vec<f64> = sel.iter().map(|&j| rc.at(*e, j)).collect();
                Ctx::normal_add(
                    &mut n3,
                    &mut b3,
                    &asd,
                    &[0, 1, 2],
                    &var,
                    &y,
                    &mut acc3,
                    None,
                    opts.compat.vtpv,
                );
                used_c.push((*e, sel));
            }
        }
        if let Some(d_code) = n3.solve(&b3)
            && let Some(n3i) = n3.inverse()
        {
            let d = [d_code[0], d_code[1], d_code[2]];
            let s0c = ((acc3.vtpv - dot(&d_code, &b3)).max(0.0)
                / (acc3.nrow as i64 - 3).max(1) as f64)
                .sqrt();
            let mut n3m = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    n3m[i][j] = n3i.at(i, j);
                }
            }
            let ch = rotate_cov3(&re, &n3m);
            let mut shc = (ch[0][0] * s0c * s0c + ch[1][1] * s0c * s0c)
                .max(0.0)
                .sqrt();
            let sh_sat = code_sat_sigma(
                &used_c,
                &rc,
                &eb_,
                d,
                &ctx.sysc,
                ns,
                &re,
                opts.code_sat_floor,
            );
            if let Some(ss) = sh_sat {
                fix_info.insert("code_sigma_h_epochs_m".into(), jround(shc, 3));
                shc = shc.max(ss);
            }
            let de = mat3_vec(&re, d);
            let dh = de[0].hypot(de[1]);
            fix_info.insert(
                "code_minus_m".into(),
                Value::Array(de.iter().map(|v| jround(*v, 3)).collect()),
            );
            fix_info.insert("code_sigma_h_m".into(), jround(shc, 3));
            if dh > opts.fail_code.max(4.0 * shc) {
                failed = format!("{dh:.2} m from the code-only solution (σ {shc:.2} m)");
            }
        }
    }
    if status_of(&fix_info) == "fixed" && failed.contains("code-only") {
        fix_info.insert("code_check".into(), Value::String(failed.clone()));
        failed = String::new();
    }
    if status_of(&fix_info) != "fixed" {
        let mut c3 = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                c3[i][j] = cov_out.at(i, j);
            }
        }
        let ce = rotate_cov3(&re, &c3);
        let sh = (ce[0][0] + ce[1][1]).max(0.0).sqrt();
        if sh > opts.fail_sigma && failed.is_empty() {
            failed = format!("float horizontal σ {:.0} mm", 1000.0 * sh);
        }
    }
    if !failed.is_empty() {
        fix_info.insert("failed".into(), Value::String(failed.clone()));
        notes.push(format!("failed: {failed}"));
    }
    let vec = sub3([x_out[0], x_out[1], x_out[2]], base_xyz);
    stages.validation = Some(ValidationStage {
        x: x_out.clone(),
        cov: cov_out.clone(),
        fix: fix_info.clone(),
        failed: failed.clone(),
        s0: s0_out,
        vector: vec,
        notes: notes.clone(),
    });
    let rejected = est.use_.zip_map(&est.wf_ph, |u, w| u && w == 0.0);
    let n_rej_total = rejected.count();
    let n_rej_code = est
        .code_ok
        .zip_map(&est.wf_cd, |c, w| c && w == 0.0)
        .count()
        + n_gross_code;
    let mut used: BTreeMap<String, Vec<Sat>> = BTreeMap::new();
    for j in 0..ns {
        if (0..ne).any(|e| est.use_.at(e, j) && !rejected.at(e, j)) {
            used.entry(ctx.sysc[j].to_string())
                .or_default()
                .push(ctx.cols[j]);
        }
    }
    let mut report: BTreeMap<Sat, Value> = BTreeMap::new();
    for (j, sat) in ctx.cols.iter().enumerate() {
        if !est.use_.any_col(j) {
            continue;
        }
        let u: Vec<usize> = (0..ne)
            .filter(|&e| {
                est.use_.at(e, j) && est.wf_ph.at(e, j) > 0.0 && est.cf_ph.at(e, j).is_finite()
            })
            .collect();
        let rms_mm = if u.is_empty() {
            None
        } else {
            let sq: Vec<f64> = u.iter().map(|&e| est.cf_ph.at(e, j).powi(2)).collect();
            Some(1000.0 * (np_sum(&sq) / sq.len() as f64).sqrt())
        };
        let down = (0..ne)
            .filter(|&e| est.use_.at(e, j) && est.wf_ph.at(e, j) > 0.0 && est.wf_ph.at(e, j) < 1.0)
            .count();
        let zero = (0..ne).filter(|&e| rejected.at(e, j)).count();
        let arcs: BTreeSet<i64> = (0..ne)
            .filter(|&e| est.use_.at(e, j))
            .map(|e| arc.at(e, j))
            .collect();
        let (rep, spk) = sat_slips.get(sat).copied().unwrap_or((0, 0));
        report.insert(
            *sat,
            sat_report_entry(
                u.len(),
                rms_mm,
                down,
                zero,
                arcs.len(),
                rep,
                spk,
                eph_age.get(sat).copied(),
            ),
        );
    }
    let mut trop: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    if ctx.nn > 0 {
        trop.insert(
            "rover_zwd".into(),
            rounded4(&x_out[ctx.i_tb..ctx.i_tb + ctx.nn]),
        );
        if ctx.n_ta > 0 {
            trop.insert(
                "base_zwd".into(),
                rounded4(&x_out[ctx.i_ta..ctx.i_ta + ctx.n_ta]),
            );
        }
    }
    let mut cov_formal = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cov_formal[i][j] = cov_out.at(i, j);
        }
    }
    let mut fdist = 0.0;
    if status_of(&fix_info) == "fixed" {
        let cf = rotate_cov3(&re, &cov_float);
        let sfh = (cf[0][0] + cf[1][1]).max(0.0).sqrt();
        let d = mat3_vec(&re, sub3(vec_float, vec));
        let dfh = d[0].hypot(d[1]);
        fdist = if sfh > 0.0 { dfh / sfh } else { 0.0 };
    }
    let hours = (t[ne - 1] - t[0]) / 3600.0;
    let n_sat: usize = used.values().map(Vec::len).sum();
    let cov_real = quality_cov(
        &cov_formal,
        &re,
        norm3(vec),
        hours,
        &mut fix_info,
        n_amb,
        n_sat,
        &opts,
        &mut notes,
        mat3_vec(&re, vec)[2],
        n_rej_total,
        est.acc.n_ph,
        s0_out,
        fdist,
    );
    stages.sigma = Some(SigmaStage {
        cov_formal,
        cov: cov_real,
        sigma_q: fix_info.get("sigma_q").and_then(Value::as_f64),
        frag: fix_info.get("frag").and_then(Value::as_f64),
        rej: fix_info.get("rej").and_then(Value::as_f64),
        partial: fix_info.get("partial").and_then(Value::as_f64),
        fdist,
        s0: s0_out,
        hours,
        length_m: norm3(vec),
        n_rejected: n_rej_total,
        n_dd: est.acc.n_ph,
        notes: notes.clone(),
    });
    let mut slip_events: BTreeMap<String, usize> = if opts.slip_mode == "gf" {
        slip_ev_gf.clone()
    } else {
        slip_ev
            .pairs()
            .iter()
            .map(|(k, v)| ((*k).to_owned(), *v))
            .collect()
    };
    slip_events.insert("split".into(), n_split);
    let fixed = status_of(&fix_info) == "fixed";
    Ok(StaticResult {
        base_xyz,
        rover_xyz: [x_out[0], x_out[1], x_out[2]],
        vector: vec,
        cov: cov_real,
        cov_formal,
        sigma0: s0_out,
        rms_phase: (est.acc.ss_ph / est.acc.n_ph.max(1) as f64).sqrt(),
        rms_code: if est.acc.n_cd > 0 {
            (est.acc.ss_cd / est.acc.n_cd as f64).sqrt()
        } else {
            f64::NAN
        },
        n_epochs: ne,
        n_dd: est.acc.n_ph,
        n_amb,
        n_rejected: n_rej_total,
        sats: used,
        pairs: ctx.cols.iter().map(|c| (*c, ctx.pairs[c])).collect(),
        trop,
        seconds: t_start.elapsed().as_secs_f64(),
        n_rejected_code: n_rej_code,
        notes,
        sat_report: report,
        slip_events,
        fixed,
        fix: fix_info,
        vector_float: Some(vec_float),
        cov_float: Some(cov_float),
        failed: !failed.is_empty(),
    })
}
