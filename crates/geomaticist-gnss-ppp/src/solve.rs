//! The batch ionosphere-free PPP, split into the
//! seven stages: observations, arcs, float, wide lane, narrow lane, validation, σ.
//!
//! Unknowns: the ARP (ECEF), a receiver clock per epoch (eliminated epoch by epoch), a code bias per
//! system against the reference, a code bias per GLONASS satellite, the zenith wet delay as
//! piecewise-linear nodes tied by a random walk, optional N / E gradients, a float ionosphere-free
//! ambiguity per arc. Data: ionosphere-free code and phase on the clock grid, elevation weighted.
//!
//! The normal matrix is assembled per epoch (`Nl = AᵀWA − g gᵀ / Σw`), solved by
//! Cholesky each iteration and inverted once after the last iteration for
//! the covariance, the ZTD σ and PPP-AR.

use std::collections::{BTreeMap, BTreeSet};

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::obs::{combine, if_coeffs, pick_pair};
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};
use geomaticist_gnss_engine::arr::{B2, F2, I2, Vec3Arr};
use geomaticist_gnss_engine::attitude::{in_shadow, nadir_angle, rel_clock, shapiro, windup};
use geomaticist_gnss_engine::linalg::{DMat, add3, dot, norm3, sub3};
use geomaticist_gnss_engine::models::{
    azel, ecef_to_llh, enu_matrix, geometry, sat_pco_offset, sun_moon_ecef,
};
use geomaticist_gnss_engine::nputil::{
    mean, median, nanmedian, np_interp, np_round, np_sum, py_round_i,
};
use geomaticist_gnss_engine::orbits::Orbits;
use geomaticist_gnss_engine::slips::{self, find_step};
use geomaticist_gnss_engine::solve::{receiver_cal, spp};
use geomaticist_gnss_engine::tides::{ocean_loading, pole_tide, solid_tide};
use geomaticist_gnss_engine::tropo::{Tropo, TropoInputs, troposphere};
use geomaticist_gnss_products::bias::{Erp, pole_at};
use geomaticist_gnss_products::{AntennaCal, Antex, receiver_correction, sat_pcv_if};
use rayon::prelude::*;
use serde_json::{Map, Value};

use crate::options::{Antenna, PppOptions, PppResult, obs_name, realistic_cov, sys_name, sys_rank};
use crate::resolve::resolve;
use crate::stages::{
    ArcsStage, FloatStage, NarrowLaneStage, ObservationsStage, SigmaStage, Stages, ValidationStage,
    WideLaneStage,
};

/// A failed solve (in `solve` / `spp`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PppError(pub String);

impl std::fmt::Display for PppError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for PppError {}

/// Satellite OSB map as the products crate reads it: (satellite, observable) -> metres.
pub type Osb = BTreeMap<(String, String), f64>;

/// Arc number per epoch (-1 where not used) from the arc starts.
#[must_use]
pub fn arcs_from(new_arc: &[bool], use_: &[bool]) -> Vec<i64> {
    let mut arc = vec![-1i64; use_.len()];
    let mut cur = -1i64;
    for k in 0..use_.len() {
        if new_arc[k] {
            cur += 1;
        }
        if use_[k] && cur >= 0 {
            arc[k] = cur;
        }
    }
    arc
}

/// `np.searchsorted(a, v, side="right")`: first index with `a[i] > v`.
fn searchsorted_right(a: &[f64], v: f64) -> usize {
    a.partition_point(|&x| x <= v)
}

/// The station's displacement (m, ECEF) per epoch: solid tide, pole tide, ocean loading (B6, when `otl`).
pub fn displacements(
    opts: &PppOptions,
    x0: [f64; 3],
    t: &[f64],
    erp: Option<&Erp>,
    blq: Option<&[[f64; 11]; 6]>,
    notes: &mut Vec<String>,
    name: &str,
) -> Vec<[f64; 3]> {
    let ne = t.len();
    let mut disp = vec![[0.0; 3]; ne];
    if opts.tides {
        for (d, s) in disp.iter_mut().zip(solid_tide(x0, t)) {
            *d = add3(*d, s);
        }
    }
    if opts.pole_tide
        && let Some(erp) = erp
    {
        let tm = mean(t);
        if let Some((xp, yp)) = pole_at(erp, tm) {
            let p = pole_tide(x0, tm, xp, yp);
            for d in &mut disp {
                *d = add3(*d, p);
            }
        }
    }
    if opts.otl {
        match blq {
            Some(b) => {
                for (d, o) in disp.iter_mut().zip(ocean_loading(b, x0, t)) {
                    *d = add3(*d, o);
                }
            }
            None => notes.push(format!(
                "{name}: ocean tide loading asked for (otl) but no BLQ coefficients given: not applied"
            )),
        }
    }
    disp
}

struct SatCol {
    sat: Sat,
    pr: (u8, u8),
    ph: Vec<f64>,
    cd: Vec<f64>,
    gf: Vec<f64>,
    mw: Vec<f64>,
    lli: Vec<bool>,
    xyz: Vec<[f64; 3]>,
    clk: Vec<f64>,
    lamif: f64,
    ok: Vec<bool>,
    pcvs: Vec<f64>,
    sh: Vec<bool>,
    a: f64,
    b: f64,
    lam1: f64,
    lam2: f64,
    osb4: Option<[f64; 4]>,
    sig: String,
}

/// One epoch's design rows after the clock elimination's inputs: the local column list, `A` (m × nl,
/// row-major), `y`, `w` and the phase / code column indexes.
struct EpochRows {
    loc: Vec<usize>,
    a: Vec<f64>,
    nl: usize,
    y: Vec<f64>,
    w: Vec<f64>,
    jp: Vec<usize>,
    jc: Vec<usize>,
    /// `Aᵀ (A w)` (nl × nl, row-major), `Σ_rows (A w)`, `(A w)ᵀ y`, `Σ w`: the epoch's block before the clock
    /// elimination, computed in parallel over the epochs and scattered into `N` in epoch order.
    blk: Vec<f64>,
    g: Vec<f64>,
    bl: Vec<f64>,
    sw: f64,
}

/// A satellite after the sequential part of the preparation (signals, OSB, antenna lookup); the
/// numerics (orbits, relativity, PCO, PCV, shadow) run in parallel over these.
struct Pre {
    sat: Sat,
    pr: (u8, u8),
    c: geomaticist_gnss_core::Combo,
    ph: Vec<f64>,
    cd: Vec<f64>,
    a: f64,
    b: f64,
    nm: [Option<String>; 4],
    v: [Option<f64>; 4],
    cal: Option<AntennaCal>,
    com: bool,
}

/// Static PPP of one station's observations (already cut to the session).
///
/// `orb` holds the precise orbits and clocks (`Orbits` over an `OrbitSource::Precise`); `osb` the
/// analysis centre's biases (empty: none); `erp` the pole; `blq` the station's ocean loading
/// coefficients, applied only when `opts.otl` (B6); `x0` the start position (else a code solution);
/// `tropo` the GPT3 / VMF3 grids the caller has.
pub fn solve(
    st: &StationObs,
    orb: &Orbits,
    antex: Option<&mut Antex>,
    opts: &PppOptions,
    fcn: &BTreeMap<Sat, i32>,
    osb: &Osb,
    erp: Option<&Erp>,
    x0: Option<[f64; 3]>,
    blq: Option<&[[f64; 11]; 6]>,
    tropo: &TropoInputs,
) -> Result<(PppResult, Stages), PppError> {
    let (res, stages) = solve_inner(st, orb, antex, opts, fcn, osb, erp, x0, blq, tropo, false)?;
    let res = res.ok_or_else(|| PppError("no result".into()))?;
    Ok((res, stages))
}

/// The first two stages of [`solve`] alone: the observations on the epoch grid with every correction that
/// does not depend on the estimate (satellite positions and clocks, biases, tides, wind-up, satellite
/// antennas, the elevation mask at `x0`) and the arcs (cycle slips, a-priori ambiguities, what may be
/// used). The kinematic solution ([`crate::kinematic`]) starts from these, so that it shares every model
/// with the static one.
pub fn prepare(
    st: &StationObs,
    orb: &Orbits,
    antex: Option<&mut Antex>,
    opts: &PppOptions,
    fcn: &BTreeMap<Sat, i32>,
    osb: &Osb,
    erp: Option<&Erp>,
    x0: Option<[f64; 3]>,
    blq: Option<&[[f64; 11]; 6]>,
    tropo: &TropoInputs,
) -> Result<Stages, PppError> {
    Ok(solve_inner(st, orb, antex, opts, fcn, osb, erp, x0, blq, tropo, true)?.1)
}

fn solve_inner(
    st: &StationObs,
    orb: &Orbits,
    mut antex: Option<&mut Antex>,
    opts: &PppOptions,
    fcn: &BTreeMap<Sat, i32>,
    osb: &Osb,
    erp: Option<&Erp>,
    x0: Option<[f64; 3]>,
    blq: Option<&[[f64; 11]; 6]>,
    tropo: &TropoInputs,
    stop_after_arcs: bool,
) -> Result<(Option<PppResult>, Stages), PppError> {
    let started = std::time::Instant::now();
    let osb_on = opts.osb && !osb.is_empty();
    let mut notes: Vec<String> = Vec::new();
    let mut stages = Stages::default();
    let name = st.name.as_str();
    // ---- epochs on the clock grid
    let tt = &st.t;
    let iv = opts.interval;
    let ei: Vec<usize> = (0..tt.len())
        .filter(|&i| ((tt[i] + 0.5 * iv).rem_euclid(iv) - 0.5 * iv).abs() < 0.5)
        .collect();
    if ei.len() < 20 {
        return Err(PppError(format!(
            "{name}: only {} epochs on the {iv:.0} s grid",
            ei.len()
        )));
    }
    let t_obs: Vec<f64> = ei.iter().map(|&i| tt[i]).collect();
    let t: Vec<f64> = t_obs.iter().map(|&x| np_round(x / iv, 0) * iv).collect();
    let ne = t.len();
    // ---- satellites
    let mut cand: Vec<Sat> = Vec::new();
    for sat in st.sats.keys() {
        let s = sat.sys.letter();
        if !opts.systems.contains(s) || opts.exclude.contains(sat) || !orb.has(*sat) {
            continue;
        }
        if sat.sys == GnssSystem::Glonass && !fcn.contains_key(sat) {
            continue;
        }
        if sat.sys == GnssSystem::BeiDou && (sat.prn <= 5 || sat.prn >= 59) {
            continue; // BeiDou GEO
        }
        cand.push(*sat);
    }
    let x0 = match x0 {
        Some(x) => x,
        None => {
            let d: Vec<f64> = tt.windows(2).map(|w| w[1] - w[0]).collect();
            let med = if d.is_empty() { f64::NAN } else { median(&d) };
            let every = py_round_i(300.0 / med.max(1.0)).max(1) as usize;
            spp(st, orb, &opts.systems, fcn, every).map_err(|e| PppError(e.0))?
        }
    };
    let (lat0, lon0, _h0) = ecef_to_llh(x0);
    let r0 = enu_matrix(lat0, lon0);
    // ---- station displacements (tides) at the epochs
    let disp = displacements(opts, x0, &t, erp, blq, &mut notes, name);
    let sun: Vec<[f64; 3]> = t.iter().map(|&x| sun_moon_ecef(x).0).collect();
    // ---- per satellite: combinations, satellite state, corrections that don't depend on the position
    let mut pre: Vec<Pre> = Vec::new();
    let (mut osb_hits, mut osb_miss) = (0usize, 0usize);
    for &sat in &cand {
        let Some(pr) = pick_pair(sat, &[st]) else {
            continue;
        };
        let Ok(c) = combine(st, sat, pr, fcn.get(&sat).copied(), None) else {
            continue;
        };
        let mut ph: Vec<f64> = ei.iter().map(|&i| c.phase_if[i]).collect();
        let mut cd: Vec<f64> = ei.iter().map(|&i| c.code_if[i]).collect();
        if !ph.iter().any(|v| v.is_finite()) {
            continue;
        }
        let (a, b) = if_coeffs(c.f1, c.f2);
        let mut nm: [Option<String>; 4] = [None, None, None, None];
        let mut v: [Option<f64>; 4] = [None; 4];
        if osb_on {
            let mut k = 0;
            for kind in ['C', 'L'] {
                for bb in [pr.0, pr.1] {
                    nm[k] = obs_name(st, sat, st.best_code(sat, kind, bb));
                    v[k] = nm[k]
                        .as_ref()
                        .and_then(|n| osb.get(&(sat.to_string(), n.clone())).copied());
                    k += 1;
                }
            }
            if let (Some(v0), Some(v1)) = (v[0], v[1]) {
                let corr = a * v0 - b * v1;
                for x in &mut cd {
                    *x -= corr;
                }
                osb_hits += 1;
            } else {
                osb_miss += 1;
            }
            if let (Some(v2), Some(v3)) = (v[2], v[3]) {
                let corr = a * v2 - b * v3;
                for x in &mut ph {
                    *x -= corr;
                }
            }
        }
        let com = orb.source(sat).is_some_and(|src| src.center_of_mass());
        let cal = antex
            .as_deref_mut()
            .and_then(|atx| atx.satellite(&sat.to_string(), t[0]));
        if com && cal.is_none() && antex.is_some() {
            notes.push(format!("{sat}: no satellite antenna in ANTEX — left out"));
            continue;
        }
        pre.push(Pre {
            sat,
            pr,
            c,
            ph,
            cd,
            a,
            b,
            nm,
            v,
            cal,
            com,
        });
    }
    // the numerics per satellite, in parallel (each satellite's columns are independent; the order of
    // `cols_data` is the candidates' order, so the result does not depend on the thread count)
    let prepared: Vec<Option<SatCol>> = pre
        .into_par_iter()
        .map(|p| {
            let Pre {
                sat,
                pr,
                c,
                ph,
                cd,
                a,
                b,
                nm,
                v,
                cal,
                com,
            } = p;
            let sys = sat.sys.letter();
            let tr: Vec<f64> = ei
                .iter()
                .zip(&t_obs)
                .map(|(&i, &to)| {
                    let c1 = c.code1[i];
                    to - (if c1.is_finite() { c1 } else { 2.2e7 }) / CLIGHT
                })
                .collect();
            let (_x, clk1, o1, _) = orb.position(sat, &tr);
            let tr2: Vec<f64> = tr.iter().zip(&clk1).map(|(x, y)| x - y).collect();
            let (mut xyz, clk, o2, _) = orb.position(sat, &tr2);
            let trp: Vec<f64> = tr2.iter().map(|x| x + 0.5).collect();
            let trm: Vec<f64> = tr2.iter().map(|x| x - 0.5).collect();
            let (xp, _, _, _) = orb.position(sat, &trp);
            let (xm, _, _, _) = orb.position(sat, &trm);
            let clk: Vec<f64> = (0..ne)
                .map(|e| clk[e] + rel_clock(xyz[e], sub3(xp[e], xm[e])))
                .collect();
            if com && cal.is_some() {
                let cal_r = cal.as_ref();
                for e in 0..ne {
                    let off = sat_pco_offset(cal_r, sys as u8, pr.0, pr.1, xyz[e], sun[e], a, b);
                    xyz[e] = add3(xyz[e], off);
                }
            }
            let ok: Vec<bool> = (0..ne)
                .map(|e| ph[e].is_finite() && cd[e].is_finite() && o1[e] && o2[e])
                .collect();
            if ok.iter().filter(|&&o| o).count() < 20 {
                return None;
            }
            let mut pcvs = vec![0.0; ne];
            if opts.sat_pcv && cal.is_some() {
                let nadir: Vec<f64> = xyz.iter().map(|&p| nadir_angle(p, x0)).collect();
                pcvs = sat_pcv_if(cal.as_ref(), sys as u8, pr.0, pr.1, a, b, &nadir);
            }
            let mut sh = vec![false; ne];
            if opts.shadow {
                let s: Vec<bool> = (0..ne)
                    .map(|e| in_shadow(xyz[e], sun[e]) && ok[e])
                    .collect();
                for k in 0..ne {
                    if s[k] {
                        for e in 0..ne {
                            if t[e] >= t[k] && t[e] <= t[k] + opts.shadow_after_s {
                                sh[e] = true;
                            }
                        }
                    }
                }
            }
            let none = "None".to_owned();
            let sig = format!(
                "{sys} {} {}",
                nm[2].as_ref().unwrap_or(&none),
                nm[3].as_ref().unwrap_or(&none)
            );
            let osb4 = match (v[0], v[1], v[2], v[3]) {
                (Some(a0), Some(a1), Some(a2), Some(a3)) => Some([a0, a1, a2, a3]),
                _ => None,
            };
            Some(SatCol {
                sat,
                pr,
                ph,
                cd,
                gf: ei.iter().map(|&i| c.gf[i]).collect(),
                mw: ei.iter().map(|&i| c.mw[i]).collect(),
                lli: ei.iter().map(|&i| c.lli[i]).collect(),
                xyz,
                clk,
                lamif: CLIGHT / (c.f1 + c.f2),
                ok,
                pcvs,
                sh,
                a,
                b,
                lam1: CLIGHT / c.f1,
                lam2: CLIGHT / c.f2,
                osb4,
                sig,
            })
        })
        .collect();
    let cols_data: Vec<SatCol> = prepared.into_iter().flatten().collect();
    if cols_data.len() < 5 {
        return Err(PppError(format!(
            "{name}: only {} satellites with precise orbits and dual-frequency data",
            cols_data.len()
        )));
    }
    let ns = cols_data.len();
    let cols: Vec<Sat> = cols_data.iter().map(|c| c.sat).collect();
    let pairs: BTreeMap<Sat, (u8, u8)> = cols_data.iter().map(|c| (c.sat, c.pr)).collect();
    let mut p_phase = F2::nans(ne, ns);
    let mut p_code = F2::nans(ne, ns);
    let mut p_gf = F2::nans(ne, ns);
    let mut p_mw = F2::nans(ne, ns);
    let mut p_lli = B2::falses(ne, ns);
    let mut sxyz = Vec3Arr::zeros(ne, ns);
    let mut sclk = F2::zeros(ne, ns);
    let mut lamif = F2::zeros(ne, ns);
    let mut ok_raw = B2::falses(ne, ns);
    let mut pcvs = F2::zeros(ne, ns);
    let mut shadow = B2::falses(ne, ns);
    for (j, c) in cols_data.iter().enumerate() {
        for e in 0..ne {
            p_phase.set(e, j, c.ph[e]);
            p_code.set(e, j, c.cd[e]);
            p_gf.set(e, j, c.gf[e]);
            p_mw.set(e, j, c.mw[e]);
            p_lli.set(e, j, c.lli[e]);
            sxyz.set(e, j, c.xyz[e]);
            sclk.set(e, j, c.clk[e]);
            lamif.set(e, j, c.lamif);
            ok_raw.set(e, j, c.ok[e]);
            pcvs.set(e, j, c.pcvs[e]);
            shadow.set(e, j, c.sh[e]);
        }
    }
    let a_if: Vec<f64> = cols_data.iter().map(|c| c.a).collect();
    let b_if: Vec<f64> = cols_data.iter().map(|c| c.b).collect();
    let lam1s: Vec<f64> = cols_data.iter().map(|c| c.lam1).collect();
    let lam2s: Vec<f64> = cols_data.iter().map(|c| c.lam2).collect();
    let osb4: Vec<Option<[f64; 4]>> = cols_data.iter().map(|c| c.osb4).collect();
    let sigs: Vec<String> = cols_data.iter().map(|c| c.sig.clone()).collect();
    drop(cols_data);
    if osb_on && osb_miss > 0 {
        notes.push(format!(
            "code biases (OSB) missing for {osb_miss} of {} satellites' signals",
            osb_hits + osb_miss
        ));
    }
    let sysc: Vec<char> = cols.iter().map(|s| s.sys.letter()).collect();
    let mut ok = ok_raw.zip_map(&shadow, |o, s| o && !s);
    let n_shadow = ok_raw.zip_map(&shadow, |o, s| o && s).count();
    if n_shadow > 0 {
        notes.push(format!(
            "{n_shadow} satellite-epochs left out in or just after the Earth's shadow (yaw unknown)"
        ));
    }
    // ---- elevation mask and receiver-independent geometry at x0
    let mut az0 = F2::zeros(ne, ns);
    let mut el0 = F2::zeros(ne, ns);
    let elmask = opts.elmask.to_radians();
    for e in 0..ne {
        let pos = add3(x0, disp[e]);
        for j in 0..ns {
            let (_rho, ee) = geometry(sxyz.at(e, j), pos);
            let (a, l, _) = azel(ee, &r0);
            az0.set(e, j, a);
            el0.set(e, j, l);
            if !(l >= elmask) {
                ok.set(e, j, false);
            }
        }
    }
    // ---- wind-up (continuous per satellite; the ambiguity takes its constant)
    let mut wu = F2::zeros(ne, ns);
    if opts.windup {
        for j in 0..ns {
            let m: Vec<usize> = (0..ne)
                .filter(|&e| sxyz.at(e, j)[0].is_finite() && ok.at(e, j))
                .collect();
            if m.is_empty() {
                continue;
            }
            let sp: Vec<[f64; 3]> = m.iter().map(|&e| sxyz.at(e, j)).collect();
            let su: Vec<[f64; 3]> = m.iter().map(|&e| sun[e]).collect();
            let w = windup(&sp, x0, &su, &r0, None);
            for (k, &e) in m.iter().enumerate() {
                wu.set(e, j, if w[k].is_finite() { w[k] } else { 0.0 });
            }
        }
    }
    if opts.debug || stop_after_arcs {
        stages.observations = Some(ObservationsStage {
            t: t.clone(),
            t_obs: t_obs.clone(),
            cols: cols.clone(),
            phase: p_phase.clone(),
            code: p_code.clone(),
            gf: p_gf.clone(),
            mw: p_mw.clone(),
            lli: p_lli.clone(),
            sxyz: sxyz.clone(),
            sclk: sclk.clone(),
            lamif: lamif.clone(),
            ok: ok.clone(),
            ok_raw: ok_raw.clone(),
            shadow: shadow.clone(),
            pcvs: pcvs.clone(),
            az: az0,
            el: el0,
            windup: wu.clone(),
            disp: disp.clone(),
            station: name.to_owned(),
            x0,
            pairs: pairs.clone(),
            signals: sigs.clone(),
            a_if: a_if.clone(),
            b_if: b_if.clone(),
            lam1: lam1s.clone(),
            lam2: lam2s.clone(),
            osb: osb4.iter().map(Option::is_some).collect(),
            osb_hits,
            osb_miss,
            notes: notes.clone(),
        });
    }
    // ---- arcs: slips (MW + GF on each satellite's own series) and phase-minus-code jumps
    let mut arc = I2::filled(ne, ns, -1);
    // per satellite (in parallel: each column's slips are its own): the local arc numbers, or None for a
    // satellite with fewer than 4 usable epochs
    let per_col: Vec<Option<Vec<i64>>> = (0..ns)
        .into_par_iter()
        .map(|j| {
            let mut u = ok.col(j);
            if u.iter().filter(|&&x| x).count() < 4 {
                return None;
            }
            let pc: Vec<f64> = (0..ne)
                .map(|e| p_phase.at(e, j) - p_code.at(e, j))
                .collect();
            let mut hard = vec![false; ne];
            let have: Vec<usize> = (0..ne).filter(|&e| u[e]).collect();
            for q in 1..have.len() {
                if (pc[have[q]] - pc[have[q - 1]]).abs() > opts.slip_pc {
                    hard[have[q]] = true;
                }
            }
            let lli_u: Vec<bool> = (0..ne).map(|e| p_lli.at(e, j) && u[e]).collect();
            let sr = slips::process(
                &t,
                &p_gf.col(j),
                &p_mw.col(j),
                &u,
                &lli_u,
                &hard,
                lam1s[j],
                lam2s[j],
                a_if[j],
                b_if[j],
                &opts.slips,
            );
            for e in 0..ne {
                u[e] &= !sr.drop[e];
            }
            let mut a = arcs_from(&sr.new_arc, &u);
            // drop short arcs
            let qs: BTreeSet<i64> = a.iter().copied().filter(|&q| q >= 0).collect();
            for q in qs {
                let ks: Vec<usize> = (0..ne).filter(|&e| a[e] == q).collect();
                if t[ks[ks.len() - 1]] - t[ks[0]] < opts.min_arc_s {
                    for &e in &ks {
                        a[e] = -1;
                    }
                }
            }
            Some(a)
        })
        .collect();
    let mut n_arcs = 0i64;
    for (j, col) in per_col.into_iter().enumerate() {
        let Some(mut a) = col else {
            for e in 0..ne {
                ok.set(e, j, false); // no arc: its few epochs are out
            }
            continue;
        };
        for e in 0..ne {
            ok.set(e, j, a[e] >= 0);
            if a[e] >= 0 {
                a[e] += n_arcs;
            }
        }
        if let Some(mx) = a.iter().copied().filter(|&q| q >= 0).max() {
            n_arcs = n_arcs.max(mx + 1);
        }
        for e in 0..ne {
            arc.set(e, j, a[e]);
        }
    }
    // renumber arcs densely
    let used_arcs: BTreeSet<i64> = arc.data.iter().copied().filter(|&q| q >= 0).collect();
    let remap: BTreeMap<i64, i64> = used_arcs
        .iter()
        .enumerate()
        .map(|(i, &q)| (q, i as i64))
        .collect();
    for v in &mut arc.data {
        *v = remap.get(v).copied().unwrap_or(-1);
    }
    let mut namb = used_arcs.len();
    // a-priori ambiguity per arc: its median phase minus code
    let mut amb0 = F2::zeros(ne, ns);
    {
        let mut per: Vec<Vec<f64>> = vec![Vec::new(); namb];
        for e in 0..ne {
            for j in 0..ns {
                let q = arc.at(e, j);
                if q >= 0 {
                    per[q as usize].push(p_phase.at(e, j) - p_code.at(e, j));
                }
            }
        }
        let meds: Vec<f64> = per.iter().map(|v| nanmedian(v)).collect();
        for e in 0..ne {
            for j in 0..ns {
                let q = arc.at(e, j);
                if q >= 0 {
                    amb0.set(e, j, meds[q as usize]);
                }
            }
        }
    }
    let mut code_ok = ok.clone();
    if sysc.contains(&'R') && !opts.glo_code {
        for j in 0..ns {
            if sysc[j] == 'R' {
                for e in 0..ne {
                    code_ok.set(e, j, false);
                }
            }
        }
    }
    let mut phase_ok = ok.clone();
    // ---- receiver antenna: each signal calibrated, copied from the same frequency, or substituted
    let cal_r = receiver_cal(antex.as_deref_mut(), st, &mut notes);
    if antex.is_some() && cal_r.is_none() && !st.antenna.is_empty() {
        notes.push("receiver antenna not modelled: heights may be off by several cm".into());
    }
    let mut antenna: Antenna = BTreeMap::new();
    let mut bad = vec![false; ns];
    if let Some(cal) = &cal_r {
        for (j, sat) in cols.iter().enumerate() {
            let s = sat.sys.letter();
            for b in [pairs[sat].0, pairs[sat].1] {
                let (stt, key) = cal.freq_status(s as u8, b);
                let text = if matches!(stt, "same frequency" | "substituted") {
                    format!("{stt} ({})", key.unwrap_or_default())
                } else {
                    stt.to_owned()
                };
                let k = format!("{s}{b:02}");
                let entry = antenna.entry(s).or_default();
                match entry.iter_mut().find(|(kk, _)| *kk == k) {
                    Some(e) => e.1 = text,
                    None => entry.push((k, text)),
                }
                bad[j] |= matches!(stt, "substituted" | "none");
            }
        }
    }
    let mut off_sys: Vec<char> = Vec::new();
    let good_sys: BTreeSet<char> = (0..ns)
        .filter(|&j| !bad[j] && ok.any_col(j))
        .map(|j| sysc[j])
        .collect();
    let mut bad_sys: Vec<char> = (0..ns)
        .filter(|&j| bad[j])
        .map(|j| sysc[j])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    bad_sys.sort_by_key(|&s| sys_rank(s));
    for s in bad_sys {
        let cal_name = cal_r.as_ref().map_or("", |c| c.name.trim());
        let what = antenna
            .get(&s)
            .map(|v| {
                v.iter()
                    .filter(|(_, x)| x.starts_with("substituted") || x.starts_with("none"))
                    .map(|(k, x)| format!("{k} {x}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            })
            .unwrap_or_default();
        let others = good_sys.iter().any(|&g| g != s && g != 'R');
        if opts.uncal == "drop" && others {
            for j in 0..ns {
                if bad[j] && sysc[j] == s {
                    for e in 0..ne {
                        phase_ok.set(e, j, false);
                        code_ok.set(e, j, false);
                    }
                }
            }
            notes.push(format!(
                "{} left out: the antenna {cal_name} has no calibration of its own for {what} (option uncal: 'offset' estimates a vertical offset for it)",
                sys_name(s)
            ));
        } else if opts.uncal == "offset" && others {
            off_sys.push(s);
            notes.push(format!(
                "{}: antenna {cal_name} not calibrated for {what} — a vertical offset of its own estimated",
                sys_name(s)
            ));
        } else {
            notes.push(format!(
                "⚠ {}: antenna {cal_name} not calibrated for {what} — another frequency's values used: heights may be off by about a centimetre",
                sys_name(s)
            ));
        }
    }
    if opts.debug || stop_after_arcs {
        stages.arcs = Some(ArcsStage {
            arc: arc.clone(),
            amb0: amb0.clone(),
            ok: ok.clone(),
            phase_ok: phase_ok.clone(),
            code_ok: code_ok.clone(),
            n_arcs: used_arcs.len(),
            antenna: antenna.clone(),
            off_sys: off_sys.clone(),
            notes: notes.clone(),
        });
    }
    if stop_after_arcs {
        return Ok((None, stages));
    }
    // ---- parameters
    let present: BTreeSet<char> = (0..ns)
        .filter(|&j| phase_ok.any_col(j))
        .map(|j| sysc[j])
        .collect();
    let systems_present: Vec<char> = crate::options::SYS_ORDER
        .chars()
        .filter(|s| present.contains(s))
        .collect();
    let ref_sys = systems_present
        .iter()
        .copied()
        .find(|&s| s != 'R')
        .unwrap_or('R');
    let isb_sys: Vec<char> = systems_present
        .iter()
        .copied()
        .filter(|&s| s != ref_sys && s != 'R')
        .collect();
    let glo_sats: Vec<usize> = if ref_sys != 'R' {
        (0..ns)
            .filter(|&j| sysc[j] == 'R' && phase_ok.any_col(j))
            .collect()
    } else {
        Vec::new()
    };
    let off_sys: Vec<char> = off_sys
        .into_iter()
        .filter(|&s| systems_present.contains(&s) && s != ref_sys)
        .collect();
    let hours = (t[ne - 1] - t[0]) / 3600.0;
    let use_grad =
        opts.gradients == "on" || (opts.gradients == "auto" && hours >= opts.gradient_min_h);
    let zstep = opts.ztd_step_h * 3600.0;
    let nz = (((t[ne - 1] - t[0]) / zstep).ceil() as usize + 1).max(2);
    let zt: Vec<f64> = (0..nz).map(|q| t[0] + q as f64 * zstep).collect();
    let gstep = opts.gradient_step_h * 3600.0;
    let ng = if use_grad {
        (((t[ne - 1] - t[0]) / gstep).ceil() as usize + 1).max(2)
    } else {
        0
    };
    let gt: Vec<f64> = (0..ng).map(|q| t[0] + q as f64 * gstep).collect();
    let i_i = 3usize;
    let i_g = i_i + isb_sys.len();
    let i_o = i_g + glo_sats.len();
    let i_z = i_o + off_sys.len();
    let i_t = i_z + nz;
    let i_a = i_t + 2 * ng;
    let mut npar = i_a + namb;
    let isb_col: BTreeMap<char, usize> = isb_sys
        .iter()
        .enumerate()
        .map(|(k, &s)| (s, i_i + k))
        .collect();
    let glo_col: BTreeMap<usize, usize> = glo_sats
        .iter()
        .enumerate()
        .map(|(k, &j)| (j, i_g + k))
        .collect();
    let off_cols: Vec<(char, usize)> = off_sys
        .iter()
        .enumerate()
        .map(|(k, &s)| (s, i_o + k))
        .collect();
    // interpolation weights of the troposphere / gradient nodes per epoch
    let kz: Vec<usize> = t
        .iter()
        .map(|&x| (searchsorted_right(&zt, x) as i64 - 1).clamp(0, nz as i64 - 2) as usize)
        .collect();
    let wz: Vec<f64> = (0..ne).map(|e| (t[e] - zt[kz[e]]) / zstep).collect();
    let (kg, wg): (Vec<usize>, Vec<f64>) = if ng > 0 {
        let kg: Vec<usize> = t
            .iter()
            .map(|&x| (searchsorted_right(&gt, x) as i64 - 1).clamp(0, ng as i64 - 2) as usize)
            .collect();
        let wg = (0..ne).map(|e| (t[e] - gt[kg[e]]) / gstep).collect();
        (kg, wg)
    } else {
        (vec![0; ne], vec![0.0; ne])
    };
    let efact: Vec<f64> = sysc.iter().map(|&s| opts.efact_of(s)).collect();
    let fif: Vec<f64> = (0..ns)
        .map(|j| (a_if[j] * a_if[j] + b_if[j] * b_if[j]).sqrt())
        .collect();
    let mut x = x0;
    let mut n_rej = 0usize;
    let mut n_split = 0usize;
    let mut sigma0 = 1.0;
    let mut rej_phase = B2::falses(ne, ns);
    let mut rej_code = B2::falses(ne, ns);
    let mut sol: Vec<f64> = Vec::new();
    let mut nmat = DMat::zeros(0);
    let mut resp = F2::nans(ne, ns);
    let mut resc = F2::nans(ne, ns);
    let mut clk = vec![f64::NAN; ne];
    let mut use_p = B2::falses(ne, ns);
    let mut use_c = B2::falses(ne, ns);
    let mut nobs = 0usize;
    let mut trop: Option<Tropo> = None;
    let mut varp = F2::zeros(ne, ns);
    let mut arc_max: i64 = namb as i64 - 1;
    for it in 0..(opts.max_iter + opts.reject_rounds) {
        let (lat, lon, h) = ecef_to_llh(x);
        let r = enu_matrix(lat, lon);
        let mut rho = F2::zeros(ne, ns);
        let mut e_los = Vec3Arr::zeros(ne, ns);
        let mut az = F2::zeros(ne, ns);
        let mut el = F2::zeros(ne, ns);
        let geo: Vec<Vec<(f64, [f64; 3], f64, f64)>> = (0..ne)
            .into_par_iter()
            .map(|e| {
                let pos = add3(x, disp[e]);
                (0..ns)
                    .map(|j| {
                        let (rr, ee) = geometry(sxyz.at(e, j), pos);
                        let (a, l, _) = azel(ee, &r);
                        (rr, ee, a, l)
                    })
                    .collect()
            })
            .collect();
        for (e, row) in geo.into_iter().enumerate() {
            for (j, (rr, ee, a, l)) in row.into_iter().enumerate() {
                rho.set(e, j, rr);
                e_los.set(e, j, ee);
                az.set(e, j, a);
                el.set(e, j, l);
            }
        }
        let tp = troposphere(&opts.trop_model, lat, lon, h, &t, &el, tropo);
        let mut ant = F2::zeros(ne, ns);
        if let Some(cal) = &cal_r {
            for (j, sat) in cols.iter().enumerate() {
                let (b1, b2) = pairs[sat];
                let azc = az.col(j);
                let elc = el.col(j);
                let sys = sat.sys.letter() as u8;
                let c1 = receiver_correction(Some(cal), sys, b1, &azc, &elc, false);
                let c2 = receiver_correction(Some(cal), sys, b2, &azc, &elc, false);
                for e in 0..ne {
                    ant.set(e, j, a_if[j] * c1[e] - b_if[j] * c2[e]);
                }
            }
        }
        let grad_on = opts.trop_grad_apriori && tp.grad != (0.0, 0.0);
        let mut yc = F2::zeros(ne, ns);
        let mut yp = F2::zeros(ne, ns);
        let mut mg = F2::zeros(ne, ns);
        let rows_model: Vec<Vec<(f64, f64, f64, f64)>> = (0..ne)
            .into_par_iter()
            .map(|e| {
                (0..ns)
                    .map(|j| {
                        let l = el.at(e, j);
                        let mgv = 1.0 / (l.sin() * l.tan() + 0.0032);
                        let mut model = rho.at(e, j) - CLIGHT * sclk.at(e, j)
                            + tp.zhd[e] * tp.mh.at(e, j)
                            + ant.at(e, j)
                            + pcvs.at(e, j);
                        if grad_on {
                            let a = az.at(e, j);
                            model += mgv * (tp.grad.0 * a.cos() + tp.grad.1 * a.sin());
                        }
                        if opts.shapiro {
                            model += shapiro(sxyz.at(e, j), x);
                        }
                        let ycv = p_code.at(e, j) - model;
                        let ypv =
                            p_phase.at(e, j) - amb0.at(e, j) - model - lamif.at(e, j) * wu.at(e, j);
                        let s_el = l.max(3f64.to_radians()).sin();
                        let v = (efact[j] * fif[j] * opts.err_phase).powi(2)
                            * (1.0 + 1.0 / (s_el * s_el));
                        (mgv, ycv, ypv, v)
                    })
                    .collect()
            })
            .collect();
        for (e, row) in rows_model.into_iter().enumerate() {
            for (j, (mgv, ycv, ypv, v)) in row.into_iter().enumerate() {
                mg.set(e, j, mgv);
                yc.set(e, j, ycv);
                yp.set(e, j, ypv);
                varp.set(e, j, v);
            }
        }
        let cr2 = opts.code_ratio * opts.code_ratio;
        for e in 0..ne {
            for j in 0..ns {
                use_p.set(
                    e,
                    j,
                    phase_ok.at(e, j) && !rej_phase.at(e, j) && yp.at(e, j).is_finite(),
                );
                use_c.set(
                    e,
                    j,
                    code_ok.at(e, j) && !rej_code.at(e, j) && yc.at(e, j).is_finite(),
                );
            }
        }
        nmat = DMat::zeros(npar);
        let mut bvec = vec![0.0; npar];
        // the design rows and the block of every epoch, in parallel; scattered into N in epoch order below
        let rows: Vec<Option<EpochRows>> = (0..ne)
            .into_par_iter()
            .map(|k| {
                let jp: Vec<usize> = (0..ns).filter(|&j| use_p.at(k, j)).collect();
                let jc: Vec<usize> = (0..ns).filter(|&j| use_c.at(k, j)).collect();
                let m = jp.len() + jc.len();
                if m < 5 || jp.len() < 4 {
                    return None;
                }
                // local columns: xyz, isb, glo, two ztd nodes, gradients, ambiguities of jp
                let colz = [i_z + kz[k], i_z + kz[k] + 1];
                let mut loc: Vec<usize> = vec![0, 1, 2];
                loc.extend(isb_sys.iter().map(|s| isb_col[s]));
                loc.extend(
                    glo_sats
                        .iter()
                        .filter(|j| jc.contains(j))
                        .map(|j| glo_col[j]),
                );
                loc.extend(off_cols.iter().map(|(_, c)| *c));
                loc.extend_from_slice(&colz);
                let g0 = i_t + 2 * kg[k];
                if ng > 0 {
                    loc.extend_from_slice(&[g0, g0 + 1, g0 + 2, g0 + 3]);
                }
                loc.extend(jp.iter().map(|&j| i_a + arc.at(k, j) as usize));
                let nl = loc.len();
                let mut pos_of = vec![usize::MAX; npar];
                for (i, &c) in loc.iter().enumerate() {
                    pos_of[c] = i;
                }
                let mut a = vec![0.0; m * nl];
                let mut y = vec![0.0; m];
                let mut w = vec![0.0; m];
                let js: Vec<usize> = jp.iter().chain(jc.iter()).copied().collect();
                let npk = jp.len();
                for (rr, &j) in js.iter().enumerate() {
                    let row = &mut a[rr * nl..(rr + 1) * nl];
                    let ee = e_los.at(k, j);
                    row[0] = -ee[0];
                    row[1] = -ee[1];
                    row[2] = -ee[2];
                    let mwk = tp.mw.at(k, j);
                    row[pos_of[colz[0]]] = mwk * (1.0 - wz[k]);
                    row[pos_of[colz[1]]] = mwk * wz[k];
                    if ng > 0 {
                        let gn = mg.at(k, j) * az.at(k, j).cos();
                        let ge = mg.at(k, j) * az.at(k, j).sin();
                        row[pos_of[g0]] = gn * (1.0 - wg[k]);
                        row[pos_of[g0 + 1]] = ge * (1.0 - wg[k]);
                        row[pos_of[g0 + 2]] = gn * wg[k];
                        row[pos_of[g0 + 3]] = ge * wg[k];
                    }
                    if rr < npk {
                        row[pos_of[i_a + arc.at(k, j) as usize]] = 1.0;
                    }
                    for &(s, c) in &off_cols {
                        // a system's own vertical antenna offset: -Δu·sin(el)
                        if sysc[j] == s {
                            row[pos_of[c]] = -el.at(k, j).sin();
                        }
                    }
                    if rr >= npk {
                        let s = sysc[j];
                        if let Some(&c) = isb_col.get(&s) {
                            row[pos_of[c]] = 1.0;
                        } else if let Some(&c) = glo_col.get(&j) {
                            row[pos_of[c]] = 1.0;
                        }
                        y[rr] = yc.at(k, j);
                        w[rr] = 1.0 / (varp.at(k, j) * cr2);
                    } else {
                        y[rr] = yp.at(k, j);
                        w[rr] = 1.0 / varp.at(k, j);
                    }
                }
                // Aᵀ (A w), Σ (A w), (A w)ᵀ y, Σ w: each element summed over the rows in order, as the
                // NumPy matmuls do (the clock elimination g gᵀ / Σw is applied at the scatter)
                let sw = np_sum(&w);
                let aw: Vec<f64> = (0..m * nl).map(|q| a[q] * w[q / nl]).collect();
                let mut g = vec![0.0; nl];
                let mut bl = vec![0.0; nl];
                let mut blk = vec![0.0; nl * nl];
                for rr in 0..m {
                    let arow = &a[rr * nl..(rr + 1) * nl];
                    let awrow = &aw[rr * nl..(rr + 1) * nl];
                    for ((gc, bc), &awc) in g.iter_mut().zip(bl.iter_mut()).zip(awrow) {
                        *gc += awc;
                        *bc += awc * y[rr];
                    }
                    for i in 0..nl {
                        let ai = arow[i];
                        for (o, &awc) in blk[i * nl..(i + 1) * nl].iter_mut().zip(awrow) {
                            *o += ai * awc;
                        }
                    }
                }
                Some(EpochRows {
                    loc,
                    a,
                    nl,
                    y,
                    w,
                    jp,
                    jc,
                    blk,
                    g,
                    bl,
                    sw,
                })
            })
            .collect();
        // eliminate the epoch clock (a column of ones) and scatter, epoch by epoch
        for rc in rows.iter().flatten() {
            let nl = rc.nl;
            let wy = dot(&rc.w, &rc.y);
            for i in 0..nl {
                let li = rc.loc[i];
                bvec[li] += rc.bl[i] - rc.g[i] * wy / rc.sw;
                for jj in 0..nl {
                    nmat.add(
                        li,
                        rc.loc[jj],
                        rc.blk[i * nl + jj] - rc.g[i] * rc.g[jj] / rc.sw,
                    );
                }
            }
        }
        // priors and constraints
        let mut pri = vec![0.0; npar];
        pri[i_z] += 1.0 / (opts.ztd_sigma * opts.ztd_sigma);
        bvec[i_z] += opts.ztd_prior / (opts.ztd_sigma * opts.ztd_sigma);
        for q in 0..nz - 1 {
            let v = (opts.ztd_rw * opts.ztd_rw) * zstep;
            let (i0, i1) = (i_z + q, i_z + q + 1);
            nmat.add(i0, i0, 1.0 / v);
            nmat.add(i1, i1, 1.0 / v);
            nmat.add(i0, i1, -1.0 / v);
            nmat.add(i1, i0, -1.0 / v);
        }
        for q in 0..ng {
            for c in 0..2 {
                pri[i_t + 2 * q + c] += 1.0 / (opts.gradient_sigma * opts.gradient_sigma);
            }
        }
        for q in 0..ng.saturating_sub(1) {
            let v = opts.gradient_rw * opts.gradient_rw * (gstep / 3600.0);
            for c in 0..2 {
                let (i0, i1) = (i_t + 2 * q + c, i_t + 2 * (q + 1) + c);
                nmat.add(i0, i0, 1.0 / v);
                nmat.add(i1, i1, 1.0 / v);
                nmat.add(i0, i1, -1.0 / v);
                nmat.add(i1, i0, -1.0 / v);
            }
        }
        for &c in glo_col.values() {
            pri[c] += 1.0 / (10.0 * 10.0);
        }
        for &(_, c) in &off_cols {
            pri[c] += 1.0 / (opts.uncal_sigma * opts.uncal_sigma);
        }
        for i in 0..npar {
            nmat.add(i, i, pri[i]);
        }
        // ambiguities of arcs without data this round (all rejected): pin them
        for i in i_a..npar {
            if nmat.at(i, i).abs() < 1e-12 {
                *nmat.at_mut(i, i) = 1.0;
            }
        }
        sol = nmat
            .solve(&bvec)
            .ok_or_else(|| PppError(format!("{name}: the PPP normal equations are singular")))?;
        // residuals, clocks
        let mut vtpv = 0.0;
        nobs = 0;
        resp = F2::nans(ne, ns);
        resc = F2::nans(ne, ns);
        clk = vec![f64::NAN; ne];
        for (k, rc) in rows.iter().enumerate() {
            let Some(rc) = rc else { continue };
            let m = rc.y.len();
            let xs: Vec<f64> = rc.loc.iter().map(|&c| sol[c]).collect();
            let mut v: Vec<f64> = (0..m)
                .map(|rr| rc.y[rr] - dot(&rc.a[rr * rc.nl..(rr + 1) * rc.nl], &xs))
                .collect();
            let ck = dot(&rc.w, &v) / np_sum(&rc.w);
            for vv in &mut v {
                *vv -= ck;
            }
            clk[k] = ck;
            let wvv: Vec<f64> = (0..m).map(|rr| rc.w[rr] * (v[rr] * v[rr])).collect();
            vtpv += np_sum(&wvv);
            nobs += m;
            for (rr, &j) in rc.jp.iter().enumerate() {
                resp.set(k, j, v[rr]);
            }
            for (rr, &j) in rc.jc.iter().enumerate() {
                resc.set(k, j, v[rc.jp.len() + rr]);
            }
        }
        let n_clk = clk.iter().filter(|c| c.is_finite()).count();
        let dof = (nobs as i64 - npar as i64 - n_clk as i64).max(1) as f64;
        sigma0 = (vtpv / dof).sqrt();
        let dx = [sol[0], sol[1], sol[2]];
        x = add3(x, dx);
        trop = Some(tp);
        if norm3(dx) > 1e-4 && it < opts.max_iter {
            continue;
        }
        // outliers: the standardized residuals, robust scale, worst first
        let zp = resp.zip_map(&varp, |r, v| r.abs() / v.sqrt());
        let zc = resc.zip_map(&varp, |r, v| r.abs() / (v * cr2).sqrt());
        let scale = |z: &F2| -> (f64, f64) {
            let fin: Vec<f64> = z.data.iter().copied().filter(|v| v.is_finite()).collect();
            if fin.is_empty() {
                (1.0, 0.0)
            } else {
                let mx = fin.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                (nanmedian(&z.data) / 0.6745, 0.5 * mx)
            }
        };
        let (sp, hp) = scale(&zp);
        let (sc, hc) = scale(&zc);
        let lim_p = opts.reject_k * sp.max(1e-9);
        let lim_c = opts.reject_k * sc.max(1e-9);
        let thr_p = lim_p.max(hp);
        let thr_c = lim_c.max(hc);
        let n2z = |v: f64| if v.is_nan() { 0.0 } else { v };
        let bad_p = zp.map(|z| n2z(z) > thr_p);
        let bad_c = zc.zip_map(&resc, |z, r| {
            n2z(z) > thr_c || n2z(r).abs() > opts.code_gross
        });
        // a step in an arc's phase residuals splits the arc
        let mut per: Vec<Vec<(usize, usize)>> = vec![Vec::new(); namb];
        for k in 0..ne {
            for j in 0..ns {
                let q = arc.at(k, j);
                if q >= 0 && (q as usize) < namb && use_p.at(k, j) {
                    per[q as usize].push((k, j));
                }
            }
        }
        let mut split = 0usize;
        for q in 0..namb {
            let kk = &per[q];
            if kk.len() < 2 * 8 {
                continue;
            }
            let rr: Vec<f64> = kk.iter().map(|&(k, j)| resp.at(k, j)).collect();
            let ss: Vec<f64> = kk
                .iter()
                .map(|&(k, j)| varp.at(k, j).sqrt() * sigma0.max(1.0))
                .collect();
            let Some((qi, size, z, zl, _dl)) = find_step(&rr, &ss, 6, 8) else {
                continue;
            };
            if z > opts.step_k && zl > opts.step_k_local && size.abs() > opts.step_min {
                let j = kk[0].1;
                let k0 = kk[qi].0;
                arc_max += 1;
                let newq = arc_max;
                for e in k0..ne {
                    if arc.at(e, j) == q as i64 {
                        arc.set(e, j, newq);
                    }
                }
                split += 1;
            }
        }
        if split > 0 {
            namb = (arc_max + 1) as usize;
            npar = i_a + namb;
            n_split += split;
            continue;
        }
        if !bad_p.any() && !bad_c.any() {
            break;
        }
        rej_phase = rej_phase.or(&bad_p);
        rej_code = rej_code.or(&bad_c);
        n_rej = rej_phase.count() + rej_code.count();
    }
    let trop = trop.ok_or_else(|| PppError(format!("{name}: no iteration ran")))?;
    if n_split > 0 {
        notes.push(format!(
            "{n_split} arc(s) split at a step in their phase residuals"
        ));
    }
    let ninv = nmat
        .inverse()
        .ok_or_else(|| PppError(format!("{name}: the PPP normal equations are singular")))?;
    let s02 = sigma0 * sigma0;
    let mut cov = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cov[i][j] = ninv.at(i, j) * s02;
        }
    }
    if opts.debug {
        stages.float = Some(FloatStage {
            x,
            n: nmat.clone(),
            ninv: ninv.clone(),
            sol: sol.clone(),
            resp: resp.clone(),
            resc: resc.clone(),
            clk: clk.clone(),
            arc: arc.clone(),
            use_p: use_p.clone(),
            use_c: use_c.clone(),
            rej_phase: rej_phase.clone(),
            rej_code: rej_code.clone(),
            amb0: amb0.clone(),
            sigma0,
            cov,
            npar,
            i_i,
            i_g,
            i_o,
            i_z,
            i_t,
            i_a,
            nz,
            ng,
            namb,
            isb_sys: isb_sys.clone(),
            ref_sys,
            glo_sats: glo_sats.clone(),
            off_sys: off_sys.clone(),
            n_rej,
            n_split,
            nobs,
            notes: notes.clone(),
        });
    }
    // ---- PPP-AR
    let mut ar: Map<String, Value> = Map::new();
    if opts.ar && osb_on {
        let out = resolve(
            opts, &t, &arc, &use_p, &p_mw, &osb4, &lam1s, &lam2s, &sysc, &amb0, &sol, &ninv,
            sigma0, i_a, namb, &sigs,
        );
        ar = out.meta;
        if let Some((dx, qxx)) = out.fix {
            ar.insert(
                "xyz_float".into(),
                Value::Array(
                    x.iter()
                        .map(|&v| geomaticist_gnss_engine::dump::jfx(v))
                        .collect(),
                ),
            );
            ar.insert(
                "cov_float".into(),
                Value::Array(
                    cov.iter()
                        .map(|r| {
                            Value::Array(
                                r.iter()
                                    .map(|&v| geomaticist_gnss_engine::dump::jfx(v))
                                    .collect(),
                            )
                        })
                        .collect(),
                ),
            );
            x = add3(x, dx);
            for i in 0..3 {
                for j in 0..3 {
                    cov[i][j] = qxx[i][j] * s02;
                }
            }
            notes.push(format!(
                "ambiguities fixed: {} of {} (ratio {:.1})",
                ar["n_fixed"],
                ar["n_cand"],
                ar["ratio"].as_f64().unwrap_or(f64::NAN)
            ));
        } else if let Some(why) = ar.get("why").and_then(Value::as_str) {
            notes.push(format!("ambiguities not fixed: {why}"));
        }
    } else if opts.ar {
        ar.insert("ok".into(), Value::Bool(false));
        ar.insert(
            "why".into(),
            Value::String("no phase biases (OSB) in the products".into()),
        );
    }
    let fixed = ar.get("ok").and_then(Value::as_bool).unwrap_or(false);
    if opts.debug {
        let wl = ar.get("wl").and_then(Value::as_array).map(|w| {
            [
                w[0].as_u64().unwrap_or(0) as usize,
                w[1].as_u64().unwrap_or(0) as usize,
            ]
        });
        stages.wide_lane = Some(WideLaneStage {
            ar: opts.ar,
            wl,
            why: ar.get("why").and_then(Value::as_str).map(str::to_owned),
            ok: ar.get("ok").and_then(Value::as_bool),
        });
        stages.narrow_lane = Some(NarrowLaneStage { ar: ar.clone() });
        stages.validation = Some(ValidationStage {
            x,
            cov,
            fixed,
            sigma0,
            n_rejected: n_rej,
            n_split,
            notes: notes.clone(),
        });
    }
    // ---- results
    let mut used: BTreeMap<char, Vec<Sat>> = BTreeMap::new();
    for (j, s) in cols.iter().enumerate() {
        if use_p.any_col(j) {
            used.entry(s.sys.letter()).or_default().push(*s);
        }
    }
    for v in used.values_mut() {
        v.sort();
    }
    if !trop.note.is_empty() {
        notes.push(trop.note.clone());
    }
    let mut ztd = Vec::with_capacity(nz);
    for q in 0..nz {
        let zq = np_interp(zt[q], &t, &trop.zhd);
        ztd.push((
            zt[q],
            zq + sol[i_z + q],
            ninv.at(i_z + q, i_z + q).sqrt() * sigma0,
        ));
    }
    let gradients: Vec<(f64, f64, f64)> = (0..ng)
        .map(|q| (gt[q], sol[i_t + 2 * q], sol[i_t + 2 * q + 1]))
        .collect();
    let isb: BTreeMap<char, f64> = isb_sys.iter().map(|&s| (s, sol[isb_col[&s]])).collect();
    let offsets: BTreeMap<char, (f64, f64)> = off_cols
        .iter()
        .map(|&(s, c)| (s, (sol[c], ninv.at(c, c).sqrt() * sigma0)))
        .collect();
    let rms_of = |res: &F2, use_: &B2| -> f64 {
        // np.sqrt(np.nanmean(rp ** 2)) on rp = res[use]
        let mut sq = Vec::new();
        let mut n = 0usize;
        for e in 0..ne {
            for j in 0..ns {
                if use_.at(e, j) {
                    let v = res.at(e, j);
                    if v.is_nan() {
                        sq.push(0.0);
                    } else {
                        sq.push(v * v);
                        n += 1;
                    }
                }
            }
        }
        if sq.is_empty() {
            0.0
        } else {
            (np_sum(&sq) / n as f64).sqrt()
        }
    };
    let res = PppResult {
        xyz: x,
        cov_formal: cov,
        sigma0,
        t0: t[0],
        t1: t[ne - 1],
        hours,
        n_epochs: ne,
        n_obs_phase: use_p.count(),
        n_obs_code: use_c.count(),
        n_amb: namb,
        n_rejected: n_rej,
        sats: used,
        rms_phase: rms_of(&resp, &use_p),
        rms_code: rms_of(&resc, &use_c),
        ztd,
        isb,
        gradients,
        notes: notes.clone(),
        seconds: started.elapsed().as_secs_f64(),
        products: String::new(),
        osb_used: osb_hits > 0,
        antenna,
        offsets,
        fixed,
        ar,
        troposphere: trop.used.to_owned(),
    };
    if opts.debug {
        let enu = res.enu_cov();
        stages.sigma = Some(SigmaStage {
            cov_formal: cov,
            cov_enu: enu,
            cov_realistic_enu: realistic_cov(&enu, hours, res.fixed),
            sigma0,
            hours,
            rms_phase: res.rms_phase,
            rms_code: res.rms_code,
            n_obs_phase: res.n_obs_phase,
            n_obs_code: res.n_obs_code,
            n_amb: res.n_amb,
            n_rejected: res.n_rejected,
            notes: notes.clone(),
        });
    }
    Ok((Some(res), stages))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arcs_from_numbers_runs() {
        let new_arc = [true, false, false, true, false, false];
        let use_ = [true, true, false, true, true, false];
        assert_eq!(arcs_from(&new_arc, &use_), vec![0, 0, -1, 1, 1, -1]);
        assert_eq!(arcs_from(&[false, true], &[true, true]), vec![-1, 0]);
        assert_eq!(searchsorted_right(&[0.0, 10.0, 20.0], 10.0), 2);
        assert_eq!(searchsorted_right(&[0.0, 10.0, 20.0], 9.0), 1);
    }
}
