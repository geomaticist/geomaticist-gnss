//! Kinematic precise point positioning: one moving receiver, precise orbits and clocks, no base.
//!
//! The observations, corrections and arcs are those of the static solution ([`crate::solve::prepare`]:
//! ionosphere-free phase and code, satellite positions and clocks with relativity, satellite antenna
//! offsets and variations, code and phase biases, solid Earth and pole tides, ocean loading, wind-up, the
//! cycle-slip arcs), so the two share every model. What differs is the estimate: instead of one position
//! for the session, a filter gives one position per epoch.
//!
//! States: the position, with no memory from one epoch to the next (set to a code position with a loose
//! prior before each update, so no assumption is made about the motion); the receiver clock, also new at
//! every epoch; one clock offset per additional system for its code; the zenith wet delay as a random
//! walk; one float ambiguity per arc. The same filter runs forward and backward and the two tracks are
//! merged ([`geomaticist_gnss_engine::ppk::combine`]).
//!
//! The result is a [`PpkResult`], the type of the differential kinematic solution: the track, the events
//! and the text outputs are used in the same way. Epochs have the status
//! [`FixStatus::Ppp`] (quality flag 6).
//!
//! What to expect: a float PPP has to *converge*. The ambiguities start from the code (metres) and are
//! separated from the position only as the satellites move across the sky, which takes tens of minutes.
//! The forward track is poor at its start and the backward track at its end; the merged track is best in
//! the middle. A session of a few minutes stays at the level of its code. The receiver's antenna is
//! taken at the place of the first epoch for the terms that change slowly with position (tides, the
//! elevation mask, wind-up), which holds for tracks of some tens of kilometres.

use std::collections::BTreeMap;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::event::{Event, Occupation};
use geomaticist_gnss_core::{Sat, StationObs};
use geomaticist_gnss_engine::arr::F2;
use geomaticist_gnss_engine::attitude::shapiro;
use geomaticist_gnss_engine::linalg::{DMat, add3, norm3, rotate_cov3, sub3};
use geomaticist_gnss_engine::models::{azel, ecef_to_llh, enu_matrix, geometry};
use geomaticist_gnss_engine::orbits::Orbits;
use geomaticist_gnss_engine::ppk::{
    self, Direction, EpochSolution, EventInterp, EventStatus, FixStatus, PpkError, PpkOptions,
    PpkResult, Progress, Stage,
};
use geomaticist_gnss_engine::solve::receiver_cal;
use geomaticist_gnss_engine::tropo::{TropoInputs, troposphere};
use geomaticist_gnss_products::bias::Erp;
use geomaticist_gnss_products::{AntennaCal, Antex, receiver_correction};
use serde::{Deserialize, Serialize};

use crate::options::PppOptions;
use crate::solve::Osb;
use crate::stages::{ArcsStage, ObservationsStage};

/// Settings of [`solve`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct KinematicOptions {
    /// The models and the screening of the static solution. Its `interval` is replaced by the data's own
    /// interval and its `min_arc_s` by [`KinematicOptions::min_arc_s`].
    pub ppp: PppOptions,
    /// Which track to return.
    pub direction: Direction,
    /// Also return the forward and backward tracks.
    pub keep_passes: bool,
    /// Standard deviation of the position before each epoch's update, m.
    pub pos_sigma0: f64,
    /// Standard deviation of the receiver clock around its code value before each update, m.
    pub clk_sigma: f64,
    /// Standard deviation of a new ambiguity, m.
    pub amb_sigma0: f64,
    /// Standard deviation of a system's code offset when first seen, m, and its random walk, m per √s.
    pub isb_sigma0: f64,
    /// Random walk of a system's code offset, m per √s.
    pub isb_rw: f64,
    /// How long code errors last, s; at shorter intervals the code is down-weighted.
    pub code_corr_s: f64,
    /// An observation that misfits by more than this many standard deviations is left out (code) or
    /// restarts its ambiguity (phase).
    pub outlier_k: f64,
    /// Arcs shorter than this (s) are not used.
    pub min_arc_s: f64,
    /// Fewest satellites for a solution.
    pub min_sats: usize,
    /// How a position at an event is interpolated.
    pub event_interp: EventInterp,
    /// Missing epochs next to an event that may be skipped.
    pub event_max_missing: usize,
    /// The filter's covariance is scaled by this factor squared for the σ reported.
    pub sigma_scale: f64,
}

impl Default for KinematicOptions {
    fn default() -> Self {
        Self {
            ppp: PppOptions::default(),
            direction: Direction::Combined,
            keep_passes: false,
            pos_sigma0: 100.0,
            clk_sigma: 100.0,
            amb_sigma0: 30.0,
            isb_sigma0: 30.0,
            isb_rw: 1e-4,
            code_corr_s: 5.0,
            outlier_k: 5.0,
            min_arc_s: 30.0,
            min_sats: 5,
            event_interp: EventInterp::Cubic4,
            event_max_missing: 0,
            sigma_scale: 3.0,
        }
    }
}

impl KinematicOptions {
    /// The options with the fields of a JSON object replaced; an unknown field is an error.
    pub fn with_json(&self, json: &str) -> Result<Self, String> {
        let over: serde_json::Map<String, serde_json::Value> =
            serde_json::from_str(json).map_err(|e| e.to_string())?;
        let mut v = serde_json::to_value(self).map_err(|e| e.to_string())?;
        for (k, val) in over {
            if v.get(&k).is_none() {
                return Err(format!("unknown option {k:?}"));
            }
            v[k] = val;
        }
        serde_json::from_value(v).map_err(|e| e.to_string())
    }
}

/// Everything [`solve`] needs.
#[derive(Debug)]
pub struct KinematicInput<'a> {
    /// The moving receiver.
    pub obs: &'a StationObs,
    /// Precise orbits and clocks.
    pub orbits: &'a Orbits,
    /// Satellite and receiver antenna calibrations.
    pub antex: Option<&'a mut Antex>,
    /// GLONASS frequency channels.
    pub fcn: &'a BTreeMap<Sat, i32>,
    /// The analysis centre's biases (empty: none).
    pub osb: &'a Osb,
    /// Earth rotation parameters, for the pole tide.
    pub erp: Option<&'a Erp>,
    /// Ocean loading coefficients, applied when `ppp.otl`.
    pub blq: Option<&'a [[f64; 11]; 6]>,
    /// Troposphere grids.
    pub tropo: &'a TropoInputs<'a>,
    /// A position near the track (default: a code solution).
    pub x0: Option<[f64; 3]>,
    /// Times to position.
    pub events: &'a [Event],
    /// Stop-and-go occupations to average.
    pub occupations: &'a [Occupation],
}

/// What is the same for both passes.
struct Model<'a> {
    o: &'a ObservationsStage,
    a: &'a ArcsStage,
    cal: Option<AntennaCal>,
    opts: &'a KinematicOptions,
    tropo: &'a TropoInputs<'a>,
    /// Per satellite: noise factor (system factor × ionosphere-free amplification).
    fact: Vec<f64>,
    /// Per satellite: index of its system's code offset state, `None` for the reference system.
    isb_of: Vec<Option<usize>>,
    n_isb: usize,
    /// Per arc: first and last epoch.
    span: Vec<(usize, usize)>,
    /// Largest number of arcs alive at once.
    cap: usize,
    dt: f64,
}

/// One satellite's terms at an epoch for a position.
#[derive(Clone, Copy)]
struct Term {
    j: usize,
    /// Code less the model, m (NaN: no code).
    yc: f64,
    /// Phase less the model, the a-priori ambiguity and the wind-up, m (NaN: no phase).
    yp: f64,
    e: [f64; 3],
    mw: f64,
    var: f64,
}

impl Model<'_> {
    /// The observed-minus-computed terms of epoch `k` at antenna-less-tide position `x`.
    fn terms(&self, k: usize, x: [f64; 3]) -> Vec<Term> {
        let (o, a) = (self.o, self.a);
        let ns = o.cols.len();
        let js: Vec<usize> = (0..ns)
            .filter(|&j| a.phase_ok.at(k, j) || a.code_ok.at(k, j))
            .collect();
        if js.is_empty() {
            return Vec::new();
        }
        let (lat, lon, h) = ecef_to_llh(x);
        let r = enu_matrix(lat, lon);
        let pos = add3(x, o.disp[k]);
        let geo: Vec<(f64, [f64; 3], f64, f64)> = js
            .iter()
            .map(|&j| {
                let (rho, e) = geometry(o.sxyz.at(k, j), pos);
                let (az, el, _) = azel(e, &r);
                (rho, e, az, el)
            })
            .collect();
        let mut el = F2::zeros(1, js.len());
        for (q, g) in geo.iter().enumerate() {
            el.set(0, q, g.3);
        }
        let tp = troposphere(
            &self.opts.ppp.trop_model,
            lat,
            lon,
            h,
            &[o.t[k]],
            &el,
            self.tropo,
        );
        js.iter()
            .enumerate()
            .map(|(q, &j)| {
                let (rho, e, az, l) = geo[q];
                let sat = o.cols[j];
                let mut ant = 0.0;
                if let Some(cal) = &self.cal {
                    let (b1, b2) = o.pairs[&sat];
                    let sys = sat.sys.letter() as u8;
                    let c1 = receiver_correction(Some(cal), sys, b1, &[az], &[l], false)[0];
                    let c2 = receiver_correction(Some(cal), sys, b2, &[az], &[l], false)[0];
                    ant = o.a_if[j] * c1 - o.b_if[j] * c2;
                }
                let mut model = rho - CLIGHT * o.sclk.at(k, j)
                    + tp.zhd[0] * tp.mh.at(0, q)
                    + ant
                    + o.pcvs.at(k, j);
                if self.opts.ppp.shapiro {
                    model += shapiro(o.sxyz.at(k, j), x);
                }
                let yc = if a.code_ok.at(k, j) {
                    o.code.at(k, j) - model
                } else {
                    f64::NAN
                };
                let yp = if a.phase_ok.at(k, j) && a.arc.at(k, j) >= 0 {
                    o.phase.at(k, j)
                        - a.amb0.at(k, j)
                        - model
                        - o.lamif.at(k, j) * o.windup.at(k, j)
                } else {
                    f64::NAN
                };
                let s = l.max(3f64.to_radians()).sin();
                let var = (self.fact[j] * self.opts.ppp.err_phase).powi(2) * (1.0 + 1.0 / (s * s));
                Term {
                    j,
                    yc,
                    yp,
                    e,
                    mw: tp.mw.at(0, q),
                    var,
                }
            })
            .collect()
    }

    /// A code-only position from `x`: the position and the clock of every code group (the reference
    /// system, then each system with its own offset), m. `ztd` is the wet delay to take out.
    fn code_position(
        &self,
        k: usize,
        mut x: [f64; 3],
        ztd: f64,
    ) -> Option<([f64; 3], Vec<f64>, Vec<Term>)> {
        let ng = 1 + self.n_isb;
        let mut clk = vec![f64::NAN; ng];
        for it in 0..8 {
            let tm = self.terms(k, x);
            let rows: Vec<&Term> = tm.iter().filter(|t| t.yc.is_finite()).collect();
            let mut on = vec![false; ng];
            for t in &rows {
                on[self.isb_of[t.j].map_or(0, |i| 1 + i)] = true;
            }
            let cols: Vec<usize> = (0..ng).filter(|g| on[*g]).collect();
            let nu = 3 + cols.len();
            if rows.len() < nu + 1 {
                return None;
            }
            // clocks start at the median of their group
            for &g in &cols {
                if !clk[g].is_finite() {
                    let mut v: Vec<f64> = rows
                        .iter()
                        .filter(|t| self.isb_of[t.j].map_or(0, |i| 1 + i) == g)
                        .map(|t| t.yc - t.mw * ztd)
                        .collect();
                    v.sort_by(f64::total_cmp);
                    clk[g] = v[v.len() / 2];
                }
            }
            let mut n = DMat::zeros(nu);
            let mut b = vec![0.0; nu];
            for t in &rows {
                let g = self.isb_of[t.j].map_or(0, |i| 1 + i);
                let c = cols.iter().position(|x| *x == g).unwrap_or(0);
                let mut a = vec![0.0; nu];
                a[0] = -t.e[0];
                a[1] = -t.e[1];
                a[2] = -t.e[2];
                a[3 + c] = 1.0;
                let v = t.yc - t.mw * ztd - clk[g];
                let w = 1.0 / t.var;
                for i in 0..nu {
                    b[i] += w * a[i] * v;
                    for q in 0..nu {
                        n.add(i, q, w * a[i] * a[q]);
                    }
                }
            }
            let dx = n.solve(&b)?;
            x = [x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]];
            for (c, &g) in cols.iter().enumerate() {
                clk[g] += dx[3 + c];
            }
            let step = norm3([dx[0], dx[1], dx[2]]);
            if !step.is_finite() || step > 1.0e7 {
                return None;
            }
            if step < 0.01 || it == 7 {
                let tm = self.terms(k, x);
                return Some((x, clk, tm));
            }
        }
        None
    }
}

/// The filter of one pass.
struct Filter<'a> {
    m: &'a Model<'a>,
    /// States: position (3), clock, wet delay, code offsets, then the ambiguity slots.
    n: usize,
    x: Vec<f64>,
    p: Vec<f64>,
    /// Slot of each arc (`usize::MAX`: not alive).
    slot_of: Vec<usize>,
    free: Vec<usize>,
    isb_set: Vec<bool>,
    x_prev: [f64; 3],
    t_prev: f64,
}

const I_CLK: usize = 3;
const I_ZTD: usize = 4;
const I_ISB: usize = 5;

impl<'a> Filter<'a> {
    fn new(m: &'a Model<'a>) -> Self {
        let n = I_ISB + m.n_isb + m.cap;
        let mut f = Self {
            m,
            n,
            x: vec![0.0; n],
            p: vec![0.0; n * n],
            slot_of: vec![usize::MAX; m.span.len()],
            free: (I_ISB + m.n_isb..n).rev().collect(),
            isb_set: vec![false; m.n_isb],
            x_prev: m.o.x0,
            t_prev: f64::NAN,
        };
        f.x[I_ZTD] = m.opts.ppp.ztd_prior;
        f.p[I_ZTD * n + I_ZTD] = m.opts.ppp.ztd_sigma * m.opts.ppp.ztd_sigma;
        f
    }

    fn reset_state(&mut self, i: usize, value: f64, var: f64) {
        let n = self.n;
        for q in 0..n {
            self.p[i * n + q] = 0.0;
            self.p[q * n + i] = 0.0;
        }
        self.x[i] = value;
        self.p[i * n + i] = var;
    }

    /// One scalar update with a row of at most three single states beside the position.
    fn update(x: &mut [f64], p: &mut [f64], n: usize, e: [f64; 3], idx: &[usize], y: f64, r: f64) {
        let mut ph = vec![0.0; n];
        for (i, v) in ph.iter_mut().enumerate() {
            let row = &p[i * n..(i + 1) * n];
            let mut s = -(row[0] * e[0] + row[1] * e[1] + row[2] * e[2]);
            for &c in idx {
                s += row[c];
            }
            *v = s;
        }
        let mut s = r - (ph[0] * e[0] + ph[1] * e[1] + ph[2] * e[2]);
        let mut pred = -(x[0] * e[0] + x[1] * e[1] + x[2] * e[2]);
        for &c in idx {
            s += ph[c];
            pred += x[c];
        }
        if !(s > 0.0) {
            return;
        }
        let nu = (y - pred) / s;
        for i in 0..n {
            x[i] += ph[i] * nu;
            let k = ph[i] / s;
            if k != 0.0 {
                let row = &mut p[i * n..(i + 1) * n];
                for (q, v) in row.iter_mut().enumerate() {
                    *v -= k * ph[q];
                }
            }
        }
    }

    /// Processes epoch `k`.
    fn step(&mut self, k: usize, reverse: bool) -> EpochSolution {
        let m = self.m;
        let (o, a, opt) = (m.o, m.a, m.opts);
        let n = self.n;
        let t = o.t_obs[k];
        let mut sol = EpochSolution::none(t);
        let dt = if self.t_prev.is_finite() {
            (t - self.t_prev).abs()
        } else {
            0.0
        };
        self.t_prev = t;
        // arcs that have ended (in this pass's direction) give their slot back
        for q in 0..m.span.len() {
            let s = self.slot_of[q];
            if s != usize::MAX {
                let over = if reverse {
                    k < m.span[q].0
                } else {
                    k > m.span[q].1
                };
                if over {
                    self.slot_of[q] = usize::MAX;
                    self.free.push(s);
                }
            }
        }
        self.p[I_ZTD * n + I_ZTD] += opt.ppp.ztd_rw * opt.ppp.ztd_rw * dt;
        for i in 0..m.n_isb {
            if self.isb_set[i] {
                self.p[(I_ISB + i) * n + I_ISB + i] += opt.isb_rw * opt.isb_rw * dt;
            }
        }
        let Some((x_lin, clk, terms)) = m.code_position(k, self.x_prev, self.x[I_ZTD]) else {
            return sol;
        };
        self.x_prev = x_lin;
        let (lat, lon, _) = ecef_to_llh(x_lin);
        let renu = enu_matrix(lat, lon);
        // position and clock are new at every epoch
        for i in 0..3 {
            self.reset_state(i, 0.0, opt.pos_sigma0 * opt.pos_sigma0);
        }
        let c0 = clk.iter().copied().find(|c| c.is_finite()).unwrap_or(0.0);
        let clk_ref = if clk[0].is_finite() { clk[0] } else { c0 };
        self.reset_state(I_CLK, clk_ref, opt.clk_sigma * opt.clk_sigma);
        for i in 0..m.n_isb {
            if !self.isb_set[i] && clk[1 + i].is_finite() {
                self.reset_state(
                    I_ISB + i,
                    clk[1 + i] - clk_ref,
                    opt.isb_sigma0 * opt.isb_sigma0,
                );
                self.isb_set[i] = true;
            }
        }
        // new arcs
        for tm in &terms {
            if !tm.yp.is_finite() {
                continue;
            }
            let q = a.arc.at(k, tm.j) as usize;
            if self.slot_of[q] == usize::MAX {
                let Some(s) = self.free.pop() else { continue };
                self.slot_of[q] = s;
                // the phase was reduced by its arc's phase-minus-code: what is left is near the clock
                let v = tm.yp - clk_ref - tm.mw * self.x[I_ZTD];
                self.reset_state(s, v, opt.amb_sigma0 * opt.amb_sigma0);
            }
        }
        let code_var = opt.ppp.code_ratio * opt.ppp.code_ratio * (opt.code_corr_s / m.dt).max(1.0);
        // rows: (term, is phase, single-state indices)
        let mut rows: Vec<(usize, bool, Vec<usize>)> = Vec::new();
        for (q, tm) in terms.iter().enumerate() {
            if tm.yp.is_finite() {
                let s = self.slot_of[a.arc.at(k, tm.j) as usize];
                if s != usize::MAX {
                    rows.push((q, true, vec![I_CLK, s]));
                }
            }
            if tm.yc.is_finite() {
                let mut idx = vec![I_CLK];
                match m.isb_of[tm.j] {
                    Some(i) if self.isb_set[i] => idx.push(I_ISB + i),
                    Some(_) => continue,
                    None => {}
                }
                rows.push((q, false, idx));
            }
        }
        let n_phase = rows.iter().filter(|r| r.1).count();
        if n_phase < opt.min_sats.max(4) {
            return sol;
        }
        // update on a copy, look at the misfits, and repeat without what does not fit
        let mut skip = vec![false; rows.len()];
        let (mut xs, mut ps) = (self.x.clone(), self.p.clone());
        for _attempt in 0..4 {
            xs.clone_from(&self.x);
            ps.clone_from(&self.p);
            let run = |xs: &mut Vec<f64>, ps: &mut Vec<f64>, skip: &[bool]| {
                // the code first: it fixes the scale of the position before the phases come
                for phase in [false, true] {
                    for (r, (q, is_p, idx)) in rows.iter().enumerate() {
                        if *is_p != phase || skip[r] {
                            continue;
                        }
                        let tm = &terms[*q];
                        // the wet delay enters with its mapping: fold it into the observation's row
                        let (y, var) = if *is_p {
                            (tm.yp, tm.var)
                        } else {
                            (tm.yc, tm.var * code_var)
                        };
                        Self::update_ztd(xs, ps, n, tm.e, idx, tm.mw, y, var);
                    }
                }
            };
            run(&mut xs, &mut ps, &skip);
            let mut worst: Option<(usize, f64)> = None;
            for (r, (q, is_p, idx)) in rows.iter().enumerate() {
                if skip[r] {
                    continue;
                }
                let tm = &terms[*q];
                let mut pred =
                    -(xs[0] * tm.e[0] + xs[1] * tm.e[1] + xs[2] * tm.e[2]) + tm.mw * xs[I_ZTD];
                for &c in idx {
                    pred += xs[c];
                }
                let (y, sd) = if *is_p {
                    (tm.yp, tm.var.sqrt())
                } else {
                    (tm.yc, (tm.var * code_var).sqrt())
                };
                let z = (y - pred).abs() / sd;
                let floor = if *is_p { 0.05 } else { 0.0 };
                if z > opt.outlier_k && (y - pred).abs() > floor && worst.is_none_or(|w| z > w.1) {
                    worst = Some((r, z));
                }
            }
            let Some((r, _)) = worst else { break };
            if rows[r].1 {
                // a phase that does not fit: its ambiguity starts again, from the updated states
                let tm = &terms[rows[r].0];
                let s = rows[r].2[1];
                let v = tm.yp - self.x[I_CLK] - tm.mw * self.x[I_ZTD];
                self.reset_state(s, v, opt.amb_sigma0 * opt.amb_sigma0);
            } else {
                skip[r] = true;
            }
        }
        self.x = xs;
        self.p = ps;
        let x = [
            x_lin[0] + self.x[0],
            x_lin[1] + self.x[1],
            x_lin[2] + self.x[2],
        ];
        self.x_prev = x;
        let cov = [
            [self.p[0], self.p[1], self.p[2]],
            [self.p[n], self.p[n + 1], self.p[n + 2]],
            [self.p[2 * n], self.p[2 * n + 1], self.p[2 * n + 2]],
        ];
        let s2 = opt.sigma_scale * opt.sigma_scale;
        let mut c = rotate_cov3(&renu, &cov);
        for row in &mut c {
            for v in row.iter_mut() {
                *v *= s2;
            }
        }
        sol.status = FixStatus::Ppp;
        sol.t = t - self.x[I_CLK] / CLIGHT;
        sol.cov_enu = c;
        sol.n_sats = n_phase as u16;
        sol.set_xyz(x);
        sol
    }

    /// [`Filter::update`] with the wet delay's mapping as the partial of its state.
    fn update_ztd(
        x: &mut [f64],
        p: &mut [f64],
        n: usize,
        e: [f64; 3],
        idx: &[usize],
        mw: f64,
        y: f64,
        r: f64,
    ) {
        // scale the wet-delay state's row and column by mw for this update, then back
        let z = I_ZTD;
        if mw == 0.0 {
            Self::update(x, p, n, e, idx, y, r);
            return;
        }
        for q in 0..n {
            p[z * n + q] *= mw;
        }
        for q in 0..n {
            p[q * n + z] *= mw;
        }
        x[z] *= mw;
        let mut all = idx.to_vec();
        all.push(z);
        Self::update(x, p, n, e, &all, y, r);
        for q in 0..n {
            p[z * n + q] /= mw;
        }
        for q in 0..n {
            p[q * n + z] /= mw;
        }
        x[z] /= mw;
    }
}

/// Kinematic PPP of `input.obs`. `progress` is called a few times per second; returning `false` stops the
/// run with [`PpkError::Cancelled`].
pub fn solve(
    input: &mut KinematicInput<'_>,
    opts: &KinematicOptions,
    progress: &mut dyn FnMut(Progress) -> bool,
) -> Result<PpkResult, PpkError> {
    let st = input.obs;
    let mut res = PpkResult {
        base_arp: [f64::NAN; 3],
        base_marker: [f64::NAN; 3],
        ..PpkResult::default()
    };
    if st.t.len() < 20 {
        return Err(PpkError::Input("fewer than 20 epochs".into()));
    }
    let dt = ppk::events::median_interval(&st.t);
    let mut popts = opts.ppp.clone();
    popts.interval = dt;
    popts.min_arc_s = opts.min_arc_s;
    if !progress(Progress {
        stage: Stage::Prepare,
        done: 0,
        total: st.t.len(),
    }) {
        return Err(PpkError::Cancelled);
    }
    let stages = crate::solve::prepare(
        st,
        input.orbits,
        input.antex.as_deref_mut(),
        &popts,
        input.fcn,
        input.osb,
        input.erp,
        input.x0,
        input.blq,
        input.tropo,
    )
    .map_err(|e| PpkError::Input(e.0))?;
    let (Some(o), Some(a)) = (stages.observations.as_ref(), stages.arcs.as_ref()) else {
        return Err(PpkError::Input("no observations".into()));
    };
    res.notes.clone_from(&a.notes);
    let mut cal_notes = Vec::new();
    let cal = receiver_cal(input.antex.as_deref_mut(), st, &mut cal_notes);
    if cal.is_none() {
        res.notes.push(
            "receiver antenna without a calibration: the positions are those of the antenna's phase centre"
                .into(),
        );
    }
    let ns = o.cols.len();
    let ne = o.t.len();
    // code groups: the first system present that is not GLONASS is the reference; GLONASS code is not used
    let mut order: Vec<char> = Vec::new();
    for sat in &o.cols {
        let s = sat.sys.letter();
        if s != 'R' && !order.contains(&s) {
            order.push(s);
        }
    }
    let ref_sys = crate::options::SYS_ORDER
        .chars()
        .find(|s| order.contains(s))
        .unwrap_or('G');
    let others: Vec<char> = order.iter().copied().filter(|s| *s != ref_sys).collect();
    let isb_of: Vec<Option<usize>> = o
        .cols
        .iter()
        .map(|sat| others.iter().position(|s| *s == sat.sys.letter()))
        .collect();
    let mut arcs = a.clone();
    for j in 0..ns {
        if o.cols[j].sys.letter() == 'R' {
            for k in 0..ne {
                arcs.code_ok.set(k, j, false);
            }
        }
    }
    let n_arcs = arcs
        .arc
        .data
        .iter()
        .copied()
        .max()
        .map_or(0, |m| (m + 1).max(0) as usize);
    let mut span = vec![(usize::MAX, 0usize); n_arcs];
    for k in 0..ne {
        for j in 0..ns {
            let q = arcs.arc.at(k, j);
            if q >= 0 && arcs.phase_ok.at(k, j) {
                let s = &mut span[q as usize];
                s.0 = s.0.min(k);
                s.1 = s.1.max(k);
            }
        }
    }
    let mut alive = vec![0i32; ne + 1];
    for s in span.iter().filter(|s| s.0 != usize::MAX) {
        alive[s.0] += 1;
        alive[s.1 + 1] -= 1;
    }
    let mut cap = 0;
    let mut cur = 0;
    for v in &alive {
        cur += v;
        cap = cap.max(cur);
    }
    let fact: Vec<f64> = (0..ns)
        .map(|j| {
            opts.ppp.efact_of(o.cols[j].sys.letter())
                * (o.a_if[j] * o.a_if[j] + o.b_if[j] * o.b_if[j]).sqrt()
        })
        .collect();
    let model = Model {
        o,
        a: &arcs,
        cal,
        opts,
        tropo: input.tropo,
        fact,
        isb_of,
        n_isb: others.len(),
        span,
        cap: cap.max(0) as usize + 2,
        dt: dt.max(1e-3),
    };
    let run = |reverse: bool, progress: &mut dyn FnMut(Progress) -> bool| {
        let mut f = Filter::new(&model);
        let mut out = vec![EpochSolution::none(0.0); ne];
        let stage = if reverse {
            Stage::Backward
        } else {
            Stage::Forward
        };
        for q in 0..ne {
            let k = if reverse { ne - 1 - q } else { q };
            out[k] = f.step(k, reverse);
            if q % 200 == 0
                && !progress(Progress {
                    stage,
                    done: q,
                    total: ne,
                })
            {
                return Err(PpkError::Cancelled);
            }
        }
        Ok(out)
    };
    let (fwd, bwd) = match opts.direction {
        Direction::Forward => (run(false, progress)?, Vec::new()),
        Direction::Backward => (Vec::new(), run(true, progress)?),
        Direction::Combined => {
            let (f, b) = std::thread::scope(|s| {
                let h = s.spawn(|| run(true, &mut |_| true));
                let f = run(false, progress);
                (f, h.join())
            });
            let b = b.map_err(|_| PpkError::Input("the backward pass failed".into()))?;
            (f?, b?)
        }
    };
    res.epochs = match opts.direction {
        Direction::Forward => fwd.clone(),
        Direction::Backward => bwd.clone(),
        Direction::Combined => {
            ppk::combine::combine(&fwd, &bwd, &PpkOptions::default(), &mut res.summary)
        }
    };
    let sm = &mut res.summary;
    sm.n_epochs = ne;
    sm.n_float = res.epochs.iter().filter(|e| e.has_position()).count();
    sm.n_none = ne - sm.n_float;
    sm.interval_s = dt;
    sm.slips = n_arcs;
    let mut sig: Vec<String> = o.signals.clone();
    sig.sort();
    sig.dedup();
    sm.signals = sig;
    if sm.n_float == 0 {
        res.warnings.push("no epoch could be positioned".into());
    }
    let far = res
        .epochs
        .iter()
        .filter(|e| e.has_position())
        .map(|e| norm3(sub3(e.xyz, o.x0)))
        .fold(0.0, f64::max);
    if far > 50_000.0 {
        res.warnings.push(format!(
            "the track goes {:.0} km from its start: tides, wind-up and the elevation mask are taken at the \
             start, which costs accuracy at that distance",
            far / 1000.0
        ));
    }
    let hours = (o.t[ne - 1] - o.t[0]) / 3600.0;
    if hours < 0.5 {
        res.warnings.push(format!(
            "the session is {:.0} minutes long: a float precise point solution needs tens of minutes to \
             converge, expect decimetres or more",
            hours * 60.0
        ));
    }
    res.events = ppk::events::interpolate(
        &res.epochs,
        input.events,
        opts.event_interp,
        opts.event_max_missing,
    );
    for e in &res.events {
        match e.status {
            EventStatus::Ok => res.summary.events_ok += 1,
            EventStatus::Gap => res.summary.events_gap += 1,
            EventStatus::NoSolution => res.summary.events_failed += 1,
        }
    }
    if res.summary.events_failed > 0 {
        res.warnings.push(format!(
            "{} of {} events have no position (outside the data, beside a gap, or no solution there)",
            res.summary.events_failed,
            res.events.len()
        ));
    }
    res.occupations = ppk::events::average_occupations(&res.epochs, input.occupations);
    res.notes.extend(cal_notes);
    if opts.keep_passes {
        res.forward = fwd;
        res.backward = bwd;
    }
    Ok(res)
}
