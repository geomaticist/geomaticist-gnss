//! The kinematic estimator: one call per epoch.
//!
//! States: the rover antenna position and one single-difference ambiguity (cycles) per satellite and
//! frequency. The position carries nothing from one epoch to the next: before each update it is set to a
//! code position with a loose prior, so no assumption is made about how the receiver moves. The
//! ambiguities are constants with a very small random walk and are restarted at a cycle slip.
//!
//! Each epoch ([`PpkFilter::step`]):
//!
//! 1. restart the ambiguities the preparation flagged, and those out of the data for too long;
//! 2. refine the code position by double-difference code (decimetres instead of metres);
//! 3. form double differences of phase and code against the highest satellite of each signal group, with
//!    their full covariance (the reference's noise is in every difference of its group);
//! 4. update the float states; a phase that still misfits afterwards has slipped unnoticed: its ambiguity
//!    is restarted and the update repeated; a code that misfits is left out;
//! 5. search the integers of the double-difference ambiguities afresh from the float states (LAMBDA with
//!    the ratio test); an integer that keeps coming out the same over `ar_confirm_s` is *confirmed*;
//! 6. with at least `ar_min_dd` confirmed double differences: compute the fixed position from this
//!    epoch's phases alone with those integers taken as known, check that every phase fits it, and check
//!    that the float states do not contradict the integers; give new satellites their integers from the
//!    fixed position;
//! 7. report the fixed position when everything passes, the float one otherwise.
//!
//! The wrong-fix safeguards, in the order they act: the ratio test on the first search; an elevation mask
//! and a minimum time in lock for candidates; a least number of ambiguities; confirmation (the same
//! integer several times over a time span, not once); the fit of every phase to the fixed position at
//! every epoch, with a limit per phase and one for their RMS; the audit against the float states at every
//! epoch (a set that a fresh search would reject by the same ratio is no longer reported as fixed, and is
//! given up when that lasts); a limit on the move from the float to the fixed position; and, outside this
//! module, the comparison of the forward and the backward pass ([`super::combine`]).
//!
//! The float states are never overwritten by a fix ("continuous" resolution): a wrong set of integers
//! cannot contaminate later epochs. The fixed position has no memory: it depends on the confirmed integers
//! and the current phases only, so two passes that hold the same integers give the same position.
//!
//! The same filter runs the backward pass: it is fed the epochs in reverse order.
//!
//! The integer search itself is the engine's [`crate::ambiguity::lambda_fix`]. The ratio test, the
//! minimum lock and the elevation mask for fixing are common practice (compare RTKLIB by T. Takasu,
//! BSD-2-Clause, whose `rtkpos.c` this code does not copy); the confirmed-integer bookkeeping, the
//! memoryless fixed position, the audit and the residual checks are this engine's own.

use geomaticist_gnss_products::receiver_correction;

use super::options::{ArMode, PpkOptions};
use super::prepare::{PreparedEpoch, Setup, slant_tropo};
use super::{EpochSolution, FixStatus};
use crate::ambiguity::lambda_fix;
use crate::linalg::{DMat, mat3_vec, mat3t_vec, norm3, rotate_cov3, sub3};
use crate::models::{azel, ecef_to_llh, enu_matrix, geometry};

/// An integer counts as confirmed after it has come out this many times (and over `ar_confirm_s`).
const CONFIRM_COUNT: u16 = 3;
/// A phase whose integer is read off the fixed position must lie this close to an integer, cycles.
const INSTANT_MAX_FRAC: f64 = 0.25;

/// One observation (satellite, frequency) of the epoch being processed.
#[derive(Debug, Clone, Copy)]
struct Ob {
    key: usize,
    /// Index into the epoch's satellites.
    s: usize,
    /// Phase less the rover's modelled range, m (what is left: clock difference + λN).
    y: f64,
    /// Code less the rover's modelled range, m (NaN without code).
    yc: f64,
    /// Variance of the single-difference phase, m².
    var: f64,
    el: f64,
    half: bool,
}

/// A double-difference row.
#[derive(Debug, Clone, Copy)]
struct Row {
    /// Partials of the position.
    h: [f64; 3],
    /// Compact state index and partial of the satellite's ambiguity (phase rows).
    a: Option<(usize, f64)>,
    /// The same for the reference's ambiguity.
    b: Option<(usize, f64)>,
    z: f64,
    /// Variance of the satellite's single difference; the reference's is shared by its group.
    var: f64,
    var_ref: f64,
    /// Rows with the same `block` share the reference's noise.
    block: usize,
    /// Index of the observation in the epoch's list.
    ob: usize,
    phase: bool,
}

/// The estimator. Create one per pass.
#[derive(Debug)]
pub struct PpkFilter<'a> {
    setup: &'a Setup<'a>,
    opts: PpkOptions,
    reverse: bool,
    nk: usize,
    /// Ambiguities, cycles.
    x: Vec<f64>,
    /// Their covariance, `nk × nk`, row-major.
    p: Vec<f64>,
    alive: Vec<bool>,
    /// Tag of the last epoch each key was observed at.
    t_seen: Vec<f64>,
    /// Tag at which each key was last restarted.
    t_reset: Vec<f64>,
    /// Backward pass: the slip flag of the epoch processed last (the later one in time).
    slip_later: Vec<bool>,
    /// The integers found, as single-difference labels (defined up to one constant per group).
    label: Vec<Option<i64>>,
    /// Tag at which each label was first found with its present value.
    label_t0: Vec<f64>,
    /// How often each label has come out again since.
    label_n: Vec<u16>,
    /// Tag since which the independent search has been contradicting a confirmed label (NaN: it is not).
    dispute_t0: Vec<f64>,
    t_prev: f64,
    last_ratio: f64,
    /// Variance factor of the code this epoch: its errors last several seconds, so at high rates each
    /// epoch brings less than its noise suggests.
    code_factor: f64,
    /// Tag since which the float states have been contradicting the confirmed integers (NaN: they do not).
    held_dispute_t0: f64,
    /// Number of ambiguity restarts made by the filter's own residual check.
    pub resets: usize,
}

fn chol_in_place(a: &mut [f64], m: usize) -> bool {
    for i in 0..m {
        for j in 0..=i {
            let mut s = a[i * m + j];
            for k in 0..j {
                s -= a[i * m + k] * a[j * m + k];
            }
            if i == j {
                if s <= 0.0 || !s.is_finite() {
                    return false;
                }
                a[i * m + i] = s.sqrt();
            } else {
                a[i * m + j] = s / a[j * m + j];
            }
        }
    }
    true
}

impl<'a> PpkFilter<'a> {
    /// A filter with no ambiguities yet. `reverse`: the epochs will be fed last to first.
    #[must_use]
    pub fn new(setup: &'a Setup<'a>, opts: &PpkOptions, reverse: bool) -> Self {
        let nk = setup.keys.len();
        Self {
            setup,
            opts: opts.clone(),
            reverse,
            nk,
            x: vec![0.0; nk],
            p: vec![0.0; nk * nk],
            alive: vec![false; nk],
            t_seen: vec![f64::NAN; nk],
            t_reset: vec![f64::NAN; nk],
            slip_later: vec![false; nk],
            label: vec![None; nk],
            label_t0: vec![f64::NAN; nk],
            label_n: vec![0; nk],
            dispute_t0: vec![f64::NAN; nk],
            t_prev: f64::NAN,
            last_ratio: 0.0,
            code_factor: 1.0,
            held_dispute_t0: f64::NAN,
            resets: 0,
        }
    }

    /// Forgets every ambiguity (a restart of the filter).
    pub fn reset(&mut self) {
        self.alive.fill(false);
        self.p.fill(0.0);
        self.label.fill(None);
        self.dispute_t0.fill(f64::NAN);
        self.held_dispute_t0 = f64::NAN;
    }

    fn kill(&mut self, k: usize) {
        let nk = self.nk;
        self.alive[k] = false;
        self.label[k] = None;
        self.dispute_t0[k] = f64::NAN;
        for j in 0..nk {
            self.p[k * nk + j] = 0.0;
            self.p[j * nk + k] = 0.0;
        }
    }

    fn start(&mut self, k: usize, value: f64, t: f64) {
        self.kill(k);
        self.alive[k] = true;
        self.x[k] = value;
        self.p[k * self.nk + k] = self.opts.amb_sigma0 * self.opts.amb_sigma0;
        self.t_reset[k] = t;
    }

    /// A value for a new ambiguity of `ob`: from the clock difference the group's other ambiguities imply,
    /// else from phase minus code.
    fn start_value(&self, obs: &[Ob], q: usize) -> Option<f64> {
        let ob = &obs[q];
        let info = &self.setup.keys[ob.key];
        let mut c: Vec<f64> = obs
            .iter()
            .enumerate()
            .filter(|(i, o)| {
                *i != q && self.alive[o.key] && self.setup.keys[o.key].group == info.group
            })
            .map(|(_, o)| o.y - self.setup.keys[o.key].lam * self.x[o.key])
            .collect();
        if c.is_empty() {
            return ob.yc.is_finite().then(|| (ob.y - ob.yc) / info.lam);
        }
        c.sort_by(f64::total_cmp);
        Some((ob.y - c[c.len() / 2]) / info.lam)
    }

    /// The rover's modelled range terms of every satellite of the epoch at `x0`: range, line of sight,
    /// azimuth, elevation, troposphere.
    fn geometry(&self, pe: &PreparedEpoch, x0: [f64; 3]) -> Vec<(f64, [f64; 3], f64, f64, f64)> {
        let (lat, lon, h) = ecef_to_llh(x0);
        let r = enu_matrix(lat, lon);
        let mut g: Vec<(f64, [f64; 3], f64, f64, f64)> = pe
            .sats
            .iter()
            .map(|s| {
                let (rho, e) = geometry(s.xyz, x0);
                let (az, el, _) = azel(e, &r);
                (rho, e, az, el, 0.0)
            })
            .collect();
        let el: Vec<f64> = g.iter().map(|x| x.3).collect();
        let tr = slant_tropo(
            &self.opts.trop_model,
            self.opts.humidity,
            self.setup.tropo,
            pe.t_tag,
            (lat, lon, h),
            &el,
        );
        for (x, t) in g.iter_mut().zip(tr) {
            x.4 = t;
        }
        g
    }

    fn observations(&self, pe: &PreparedEpoch, g: &[(f64, [f64; 3], f64, f64, f64)]) -> Vec<Ob> {
        let o = &self.opts;
        let mask = o.elmask.to_radians();
        let mut obs = Vec::with_capacity(pe.sats.len() * 2);
        for (s, ps) in pe.sats.iter().enumerate() {
            let (rho, _, az, el, trop) = g[s];
            if el < mask {
                continue;
            }
            let sat = self.setup.plan.sats[usize::from(ps.slot)].sat;
            let sin = el.sin();
            let f = if sat.sys == geomaticist_gnss_core::GnssSystem::Glonass {
                o.efact_glonass
            } else {
                1.0
            };
            let var = 2.0
                * f
                * f
                * (o.err_phase * o.err_phase + o.err_phase_el * o.err_phase_el / (sin * sin));
            for b in ps.bands() {
                let key = usize::from(b.key);
                let ant = if self.setup.rover_cal.is_some() {
                    receiver_correction(
                        self.setup.rover_cal.as_ref(),
                        sat.sys.letter() as u8,
                        self.setup.keys[key].band,
                        &[az],
                        &[el],
                        o.pcv_azimuth,
                    )[0]
                } else {
                    0.0
                };
                let model = rho + trop + ant;
                obs.push(Ob {
                    key,
                    s,
                    y: b.phase - model,
                    yc: b.code - model,
                    var,
                    el,
                    half: b.half,
                });
            }
        }
        obs
    }

    /// The double-difference rows of the epoch. `cix[key]` is the compact state index of a key.
    fn rows(
        &self,
        obs: &[Ob],
        g: &[(f64, [f64; 3], f64, f64, f64)],
        cix: &[usize],
        code_off: &[bool],
        phase_rows: bool,
    ) -> Vec<Row> {
        let ng = self.setup.plan.groups.len();
        let cr2 = self.opts.code_ratio * self.opts.code_ratio * self.code_factor;
        let mut rows = Vec::with_capacity(2 * obs.len());
        for grp in 0..ng {
            let members: Vec<usize> = (0..obs.len())
                .filter(|&q| self.setup.keys[obs[q].key].group == grp)
                .collect();
            if members.len() < 2 {
                continue;
            }
            if phase_rows {
                let r = members
                    .iter()
                    .copied()
                    .max_by(|a, b| obs[*a].el.total_cmp(&obs[*b].el))
                    .unwrap_or(members[0]);
                let (er, lr) = (g[obs[r].s].1, self.setup.keys[obs[r].key].lam);
                for &q in members.iter().filter(|q| **q != r) {
                    let e = g[obs[q].s].1;
                    rows.push(Row {
                        h: [-(e[0] - er[0]), -(e[1] - er[1]), -(e[2] - er[2])],
                        a: Some((cix[obs[q].key], self.setup.keys[obs[q].key].lam)),
                        b: Some((cix[obs[r].key], -lr)),
                        z: obs[q].y - obs[r].y,
                        var: obs[q].var,
                        var_ref: obs[r].var,
                        block: 2 * grp,
                        ob: q,
                        phase: true,
                    });
                }
            }
            let coded: Vec<usize> = members
                .iter()
                .copied()
                .filter(|q| obs[*q].yc.is_finite() && !code_off[*q])
                .collect();
            if coded.len() < 2 {
                continue;
            }
            let r = coded
                .iter()
                .copied()
                .max_by(|a, b| obs[*a].el.total_cmp(&obs[*b].el))
                .unwrap_or(coded[0]);
            let er = g[obs[r].s].1;
            for &q in coded.iter().filter(|q| **q != r) {
                let e = g[obs[q].s].1;
                rows.push(Row {
                    h: [-(e[0] - er[0]), -(e[1] - er[1]), -(e[2] - er[2])],
                    a: None,
                    b: None,
                    z: obs[q].yc - obs[r].yc,
                    var: obs[q].var * cr2,
                    var_ref: obs[r].var * cr2,
                    block: 2 * grp + 1,
                    ob: q,
                    phase: false,
                });
            }
        }
        rows
    }

    /// A position from double-difference code alone, starting at `x0`.
    fn code_refine(&self, pe: &PreparedEpoch, x0: [f64; 3]) -> [f64; 3] {
        let mut x = x0;
        for _ in 0..3 {
            let g = self.geometry(pe, x);
            let obs = self.observations(pe, &g);
            let rows = self.rows(&obs, &g, &[], &vec![false; obs.len()], false);
            if rows.len() < 4 {
                return x0;
            }
            let mut n = [[0.0; 3]; 3];
            let mut b = [0.0; 3];
            for r in &rows {
                let w = 1.0 / (r.var + r.var_ref);
                for i in 0..3 {
                    b[i] += w * r.h[i] * r.z;
                    for j in 0..3 {
                        n[i][j] += w * r.h[i] * r.h[j];
                    }
                }
            }
            let Some(dx) = crate::linalg::solve3(&n, b) else {
                return x0;
            };
            if !dx.iter().all(|v| v.is_finite()) || norm3(dx) > 1.0e4 {
                return x0;
            }
            x = [x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]];
            if norm3(dx) < 0.05 {
                break;
            }
        }
        x
    }

    fn single(&self, pe: &PreparedEpoch) -> EpochSolution {
        let mut sol = EpochSolution::none(pe.t_tag);
        sol.t = pe.t;
        sol.age = pe.age;
        if let Some(x) = pe.spp {
            let (lat, lon, _) = ecef_to_llh(x);
            let r = enu_matrix(lat, lon);
            let s2 = self.opts.sigma_single * self.opts.sigma_single;
            let mut c = rotate_cov3(&r, &pe.spp_q);
            for row in &mut c {
                for v in row.iter_mut() {
                    *v *= s2;
                }
            }
            sol.status = FixStatus::Single;
            sol.cov_enu = c;
            sol.n_sats = pe.spp_n;
            sol.vel_enu = vel_enu(pe, &r);
            sol.set_xyz(self.to_mark(x, &r));
        }
        sol
    }

    /// The rover's mark below its antenna reference point (the antenna height taken along the vertical).
    fn to_mark(&self, x: [f64; 3], r: &[[f64; 3]; 3]) -> [f64; 3] {
        let d = self.setup.rover_delta_hen;
        if d == [0.0; 3] {
            return x;
        }
        sub3(x, mat3t_vec(r, [d[1], d[2], d[0]]))
    }

    /// Processes one epoch and returns its solution.
    #[allow(clippy::too_many_lines)]
    pub fn step(&mut self, pe: &PreparedEpoch) -> EpochSolution {
        let nk = self.nk;
        let t = pe.t_tag;
        let o = self.opts.clone();
        // --- time update
        let dt = if self.t_prev.is_finite() {
            (t - self.t_prev).abs()
        } else {
            0.0
        };
        self.t_prev = t;
        let qn = o.amb_rw * o.amb_rw * dt;
        if dt > 0.0 {
            self.code_factor = (o.code_corr_s / dt).max(1.0);
        }
        for k in 0..nk {
            if !self.alive[k] {
                continue;
            }
            if (t - self.t_seen[k]).abs() > self.setup.outage_s {
                self.kill(k);
            } else {
                self.p[k * nk + k] += qn;
            }
        }
        let Some(spp) = pe.spp else {
            return self.single(pe);
        };
        if pe.sats.is_empty() {
            return self.single(pe);
        }

        // --- geometry at a decimetre-level position
        let x0 = self.code_refine(pe, spp);
        let g = self.geometry(pe, x0);
        let obs = self.observations(pe, &g);
        let (lat, lon, _) = ecef_to_llh(x0);
        let renu = enu_matrix(lat, lon);

        // --- restarts asked for by the preparation
        for ps in &pe.sats {
            for b in ps.bands() {
                let k = usize::from(b.key);
                let slip = if self.reverse {
                    let s = self.slip_later[k];
                    self.slip_later[k] = b.slip;
                    s
                } else {
                    b.slip
                };
                if slip && self.alive[k] {
                    self.kill(k);
                }
            }
        }
        // keys seen now; new ones start after the surviving ones are known
        for q in 0..obs.len() {
            let k = obs[q].key;
            self.t_seen[k] = t;
        }
        for q in 0..obs.len() {
            let k = obs[q].key;
            if !self.alive[k] {
                match self.start_value(&obs, q) {
                    Some(v) => self.start(k, v, t),
                    None => self.t_seen[k] = f64::NAN,
                }
            }
        }
        let obs: Vec<Ob> = obs.into_iter().filter(|ob| self.alive[ob.key]).collect();

        // --- float update, repeated when a phase or a code misfits
        let mut code_off = vec![false; obs.len()];
        let mut upd: Option<Update> = None;
        for _attempt in 0..6 {
            let Some(u) = self.update(&obs, &g, &code_off) else {
                break;
            };
            // the worst misfit after the update
            let mut worst: Option<(usize, f64, bool)> = None;
            for (r, row) in u.rows.iter().enumerate() {
                let s = (row.var + row.var_ref).sqrt();
                let v = u.res[r].abs();
                let lim = if row.phase {
                    (o.outlier_k * s).max(0.03)
                } else {
                    (o.outlier_k * s).min(o.code_gross)
                };
                if v > lim && worst.is_none_or(|w| v / s > w.1) {
                    worst = Some((r, v / s, row.phase));
                }
            }
            match worst {
                None => {
                    upd = Some(u);
                    break;
                }
                Some((r, _, true)) => {
                    // an unnoticed slip: restart that ambiguity (the reference's, when most of its group misfit)
                    let row = u.rows[r];
                    let grp_rows: Vec<usize> = (0..u.rows.len())
                        .filter(|&i| u.rows[i].phase && u.rows[i].block == row.block)
                        .collect();
                    let n_bad = grp_rows
                        .iter()
                        .filter(|&&i| {
                            let s = (u.rows[i].var + u.rows[i].var_ref).sqrt();
                            u.res[i].abs() > (o.outlier_k * s).max(0.03)
                        })
                        .count();
                    let q = if n_bad * 2 > grp_rows.len() && grp_rows.len() > 2 {
                        u.ref_ob[row.block / 2]
                    } else {
                        row.ob
                    };
                    let k = obs[q].key;
                    self.kill(k);
                    self.alive[k] = false;
                    if let Some(v) = self.start_value(&obs, q) {
                        self.start(k, v, t);
                    } else {
                        // cannot restart without a reference: give up on this epoch's phase for it
                        self.alive[k] = true;
                        self.x[k] = 0.0;
                        self.p[k * nk + k] = o.amb_sigma0 * o.amb_sigma0 * 1e6;
                        self.t_reset[k] = t;
                    }
                    self.resets += 1;
                    upd = Some(u);
                }
                Some((r, _, false)) => {
                    code_off[u.rows[r].ob] = true;
                    upd = Some(u);
                }
            }
        }
        let Some(u) = upd else {
            return self.single(pe);
        };
        // accept the float states
        for (c, &k) in u.keys.iter().enumerate() {
            self.x[k] = u.xc[3 + c];
            for (d, &k2) in u.keys.iter().enumerate() {
                self.p[k * nk + k2] = u.pc[(3 + c) * u.na + 3 + d];
            }
        }
        let mut sats: Vec<usize> = u
            .rows
            .iter()
            .filter(|r| r.phase)
            .map(|r| obs[r.ob].s)
            .collect();
        sats.extend(
            u.ref_ob
                .iter()
                .filter(|q| **q != usize::MAX)
                .map(|q| obs[*q].s),
        );
        sats.sort_unstable();
        sats.dedup();
        let n_sats = sats.len();
        let pxx = |pc: &[f64], na: usize| -> [[f64; 3]; 3] {
            [
                [pc[0], pc[1], pc[2]],
                [pc[na], pc[na + 1], pc[na + 2]],
                [pc[2 * na], pc[2 * na + 1], pc[2 * na + 2]],
            ]
        };
        let cov_float = pxx(&u.pc, u.na);
        let sig_float = (cov_float[0][0] + cov_float[1][1] + cov_float[2][2])
            .max(0.0)
            .sqrt();
        if n_sats < o.min_sats.max(4) || sig_float > 10.0 {
            return self.single(pe);
        }
        let x_float = [x0[0] + u.xc[0], x0[1] + u.xc[1], x0[2] + u.xc[2]];
        let mut sol = EpochSolution::none(t);
        sol.t = pe.t;
        sol.age = pe.age;
        sol.n_sats = n_sats as u16;
        sol.vel_enu = vel_enu(pe, &renu);
        sol.status = FixStatus::Float;
        let mut xyz = x_float;
        let fs = o.sigma_float_scale * o.sigma_float_scale;
        let mut cov = rotate_cov3(&renu, &cov_float);
        for row in &mut cov {
            for v in row.iter_mut() {
                *v *= fs;
            }
        }

        // --- integers
        if o.ar == ArMode::Continuous {
            match self.fix(&u, &obs, &g, t, sig_float) {
                Some(f) => {
                    sol.ratio = f.ratio.min(999.9) as f32;
                    if f.confirmed {
                        sol.status = FixStatus::Fixed;
                        xyz = [x0[0] + f.dx[0], x0[1] + f.dx[1], x0[2] + f.dx[2]];
                        let l = norm3(sub3(xyz, self.setup.base_arp));
                        let s2 = o.sigma_fixed_scale * o.sigma_fixed_scale;
                        let ppm2 = (o.sigma_ppm * 1e-6 * l).powi(2);
                        cov = rotate_cov3(&renu, &f.cov);
                        for (i, row) in cov.iter_mut().enumerate() {
                            for v in row.iter_mut() {
                                *v *= s2;
                            }
                            let floor = if i == 2 {
                                o.sigma_fixed_v
                            } else {
                                o.sigma_fixed_h
                            };
                            row[i] += floor * floor + ppm2;
                        }
                    }
                }
                None => sol.ratio = self.last_ratio as f32,
            }
        }
        sol.cov_enu = cov;
        sol.set_xyz(self.to_mark(xyz, &renu));
        sol
    }
}

/// The result of one float update, in compact state order: position (3), then the keys in `keys`.
#[derive(Debug, Clone)]
struct Update {
    na: usize,
    keys: Vec<usize>,
    /// Compact index of each key (`usize::MAX`: not a state).
    cix: Vec<usize>,
    xc: Vec<f64>,
    pc: Vec<f64>,
    rows: Vec<Row>,
    /// Residuals after the update.
    res: Vec<f64>,
    /// Per group: the observation that is the phase reference (`usize::MAX`: none).
    ref_ob: Vec<usize>,
}

/// A validated set of integers.
#[derive(Debug, Clone)]
struct Fix {
    ratio: f64,
    /// Fixed position less the linearisation point.
    dx: [f64; 3],
    cov: [[f64; 3]; 3],
    confirmed: bool,
}

impl PpkFilter<'_> {
    fn update(
        &self,
        obs: &[Ob],
        g: &[(f64, [f64; 3], f64, f64, f64)],
        code_off: &[bool],
    ) -> Option<Update> {
        let nk = self.nk;
        let keys: Vec<usize> = (0..nk).filter(|&k| self.alive[k]).collect();
        let na = 3 + keys.len();
        let mut cix = vec![usize::MAX; nk];
        for (c, &k) in keys.iter().enumerate() {
            cix[k] = 3 + c;
        }
        let rows = self.rows(obs, g, &cix, code_off, true);
        let m = rows.len();
        if m < 3 {
            return None;
        }
        // compact prior
        let mut xc = vec![0.0; na];
        let mut pc = vec![0.0; na * na];
        let s0 = self.opts.pos_sigma0 * self.opts.pos_sigma0;
        for i in 0..3 {
            pc[i * na + i] = s0;
        }
        for (c, &k) in keys.iter().enumerate() {
            xc[3 + c] = self.x[k];
            for (d, &k2) in keys.iter().enumerate() {
                pc[(3 + c) * na + 3 + d] = self.p[k * nk + k2];
            }
        }
        // G[r] = P hᵀ of row r
        let mut gm = vec![0.0; m * na];
        let mut nu = vec![0.0; m];
        for (r, row) in rows.iter().enumerate() {
            let gr = &mut gm[r * na..(r + 1) * na];
            for i in 0..na {
                let pi = &pc[i * na..(i + 1) * na];
                let mut v = pi[0] * row.h[0] + pi[1] * row.h[1] + pi[2] * row.h[2];
                if let Some((ia, ca)) = row.a {
                    v += pi[ia] * ca;
                }
                if let Some((ib, cb)) = row.b {
                    v += pi[ib] * cb;
                }
                gr[i] = v;
            }
            let mut pred = 0.0;
            if let Some((ia, ca)) = row.a {
                pred += xc[ia] * ca;
            }
            if let Some((ib, cb)) = row.b {
                pred += xc[ib] * cb;
            }
            nu[r] = row.z - pred;
        }
        // S = H G + R (lower triangle)
        let mut s = vec![0.0; m * m];
        for r in 0..m {
            let row = &rows[r];
            for c in 0..=r {
                let gc = &gm[c * na..(c + 1) * na];
                let mut v = row.h[0] * gc[0] + row.h[1] * gc[1] + row.h[2] * gc[2];
                if let Some((ia, ca)) = row.a {
                    v += gc[ia] * ca;
                }
                if let Some((ib, cb)) = row.b {
                    v += gc[ib] * cb;
                }
                if rows[c].block == row.block {
                    v += row.var_ref;
                }
                if c == r {
                    v += row.var;
                }
                s[r * m + c] = v;
            }
        }
        if !chol_in_place(&mut s, m) {
            return None;
        }
        // W = L⁻¹ G, y = L⁻¹ ν
        let mut y = nu.clone();
        for r in 0..m {
            let (done, cur) = gm.split_at_mut(r * na);
            let wr = &mut cur[..na];
            let mut yr = y[r];
            for c in 0..r {
                let l = s[r * m + c];
                if l != 0.0 {
                    let wc = &done[c * na..(c + 1) * na];
                    for i in 0..na {
                        wr[i] -= l * wc[i];
                    }
                    yr -= l * y[c];
                }
            }
            let d = 1.0 / s[r * m + r];
            for v in wr.iter_mut() {
                *v *= d;
            }
            y[r] = yr * d;
        }
        // x += Wᵀ y, P -= Wᵀ W
        for r in 0..m {
            let wr = &gm[r * na..(r + 1) * na];
            let yr = y[r];
            for i in 0..na {
                xc[i] += wr[i] * yr;
                let wi = wr[i];
                if wi != 0.0 {
                    let pi = &mut pc[i * na..(i + 1) * na];
                    for j in 0..na {
                        pi[j] -= wi * wr[j];
                    }
                }
            }
        }
        // residuals after the update
        let res: Vec<f64> = rows
            .iter()
            .map(|row| {
                let mut pred = row.h[0] * xc[0] + row.h[1] * xc[1] + row.h[2] * xc[2];
                if let Some((ia, ca)) = row.a {
                    pred += xc[ia] * ca;
                }
                if let Some((ib, cb)) = row.b {
                    pred += xc[ib] * cb;
                }
                row.z - pred
            })
            .collect();
        // the phase reference of each group: the observation whose key is `b` of the group's rows
        let mut ref_ob = vec![usize::MAX; self.setup.plan.groups.len()];
        for row in rows.iter().filter(|r| r.phase) {
            if let Some((ib, _)) = row.b
                && let Some(q) = obs.iter().position(|ob| cix[ob.key] == ib)
            {
                ref_ob[row.block / 2] = q;
            }
        }
        Some(Update {
            na,
            keys,
            cix,
            xc,
            pc,
            rows,
            res,
            ref_ob,
        })
    }

    /// Whether the integer of key `k` has come out the same often and long enough to be relied on.
    fn confirmed(&self, k: usize, t: f64) -> bool {
        self.label[k].is_some()
            && self.label_n[k] >= CONFIRM_COUNT
            && (t - self.label_t0[k]).abs() >= self.opts.ar_confirm_s - 1e-9
    }

    /// Stores the integer found for key `k`: the same as before counts as one more confirmation, another
    /// one starts over.
    fn set_label(&mut self, k: usize, v: i64, t: f64) {
        match self.label[k] {
            Some(l) if l == v => self.label_n[k] = self.label_n[k].saturating_add(1),
            _ => {
                self.label[k] = Some(v);
                self.label_t0[k] = t;
                self.label_n[k] = 1;
            }
        }
        self.dispute_t0[k] = f64::NAN;
    }

    /// LAMBDA with the ratio test on float double differences `a` with covariance `q`; when the full set
    /// fails, the least certain fifth is left out and the rest searched again, up to `ar_partial_rounds`
    /// times and never below `min_n`. Returns the indices kept, their integers and the ratio; the best
    /// ratio seen is left in `last_ratio`.
    fn search(&mut self, a: &[f64], q: &DMat, min_n: usize) -> Option<(Vec<usize>, Vec<f64>, f64)> {
        let mut set: Vec<usize> = (0..a.len()).collect();
        let mut rounds = 0;
        loop {
            let nd = set.len();
            if nd < min_n.max(1) {
                return None;
            }
            let asub: Vec<f64> = set.iter().map(|&i| a[i]).collect();
            let qsub = q.sub(&set);
            let found = lambda_fix(&asub, &qsub, 2).ok();
            let ratio = found.as_ref().map_or(0.0, |(_, s)| {
                if s[0] > 0.0 {
                    (s[1] / s[0]).min(999.9)
                } else {
                    999.9
                }
            });
            self.last_ratio = self.last_ratio.max(ratio);
            if let Some((cands, _)) = found
                && ratio >= self.opts.ar_ratio
            {
                return Some((set, cands[0].clone(), ratio));
            }
            rounds += 1;
            let drop = (nd / 5).max(1);
            if rounds > self.opts.ar_partial_rounds || nd < min_n.max(1) + drop {
                return None;
            }
            set.sort_by(|x, y| q.at(*x, *x).total_cmp(&q.at(*y, *y)));
            set.truncate(nd - drop);
            set.sort_unstable();
        }
    }

    /// The position from the phases of the keys in `held` (indices into `obs`) with their integers taken
    /// as known: per signal group one clock, and the three position components, by weighted least squares
    /// on this epoch's phases alone. Returns the position (less the linearisation point), its covariance,
    /// the clock of each group (NaN where the group has no held key), the residual of each held key and the
    /// number of double differences behind the position.
    #[allow(clippy::type_complexity)]
    fn fixed_position(
        &self,
        obs: &[Ob],
        g: &[(f64, [f64; 3], f64, f64, f64)],
        held: &[usize],
    ) -> Option<([f64; 3], [[f64; 3]; 3], Vec<f64>, Vec<f64>, usize)> {
        let ng = self.setup.plan.groups.len();
        // phase less the integer; the first held key of a group is taken out to keep the numbers small
        let mut first = vec![f64::NAN; ng];
        let mut rows: Vec<(usize, [f64; 3], f64, f64)> = Vec::with_capacity(held.len());
        for &q in held {
            let ob = &obs[q];
            let info = &self.setup.keys[ob.key];
            let r = ob.y - info.lam * self.label[ob.key].unwrap_or(0) as f64;
            if !first[info.group].is_finite() {
                first[info.group] = r;
            }
            let e = g[ob.s].1;
            rows.push((
                info.group,
                [-e[0], -e[1], -e[2]],
                r - first[info.group],
                1.0 / ob.var,
            ));
        }
        let mut sw = vec![0.0; ng];
        let mut swh = vec![[0.0; 3]; ng];
        let mut swr = vec![0.0; ng];
        let mut cnt = vec![0usize; ng];
        let mut n = [[0.0; 3]; 3];
        let mut b = [0.0; 3];
        for (grp, h, r, w) in &rows {
            sw[*grp] += w;
            swr[*grp] += w * r;
            cnt[*grp] += 1;
            for i in 0..3 {
                swh[*grp][i] += w * h[i];
                b[i] += w * h[i] * r;
                for j in 0..3 {
                    n[i][j] += w * h[i] * h[j];
                }
            }
        }
        for grp in 0..ng {
            if sw[grp] > 0.0 {
                for i in 0..3 {
                    b[i] -= swh[grp][i] * swr[grp] / sw[grp];
                    for j in 0..3 {
                        n[i][j] -= swh[grp][i] * swh[grp][j] / sw[grp];
                    }
                }
            }
        }
        let nd: usize = cnt.iter().map(|c| c.saturating_sub(1)).sum();
        if nd < 3 {
            return None;
        }
        let cov = crate::linalg::inv3(&n)?;
        let dx = mat3_vec(&cov, b);
        if !dx.iter().all(|v| v.is_finite()) {
            return None;
        }
        let mut clk = vec![f64::NAN; ng];
        for grp in 0..ng {
            if sw[grp] > 0.0 {
                let hd = swh[grp][0] * dx[0] + swh[grp][1] * dx[1] + swh[grp][2] * dx[2];
                clk[grp] = (swr[grp] - hd) / sw[grp];
            }
        }
        let res: Vec<f64> = rows
            .iter()
            .map(|(grp, h, r, _)| r - clk[*grp] - (h[0] * dx[0] + h[1] * dx[1] + h[2] * dx[2]))
            .collect();
        // the clocks with the offset that was taken out put back
        for grp in 0..ng {
            clk[grp] += first[grp];
        }
        Some((dx, cov, clk, res, nd))
    }

    /// The integers, after a float update.
    ///
    /// Every key carries the integer found for it so far (a *label*: its single-difference integer up to
    /// one constant per group) with the time it was first found and how often it has come out again. An
    /// integer is *confirmed* once it has come out the same at least three times over `ar_confirm_s`.
    ///
    /// 1. All candidate double differences are searched from the float states alone (LAMBDA, ratio
    ///    test). This independent search starts the labels when there are too few confirmed ones, and
    ///    afterwards audits them: a confirmed integer the search keeps contradicting for `ar_confirm_s`
    ///    is replaced (and has to be confirmed again).
    /// 2. With the confirmed integers taken as known, the position follows from this epoch's phases
    ///    alone, with no memory: this is the fixed solution. The phases must fit it (each within three
    ///    times `ar_max_res`, their RMS within `ar_max_res`); a satellite that does not fit loses its
    ///    label. The epoch is fixed when at least `ar_min_dd` double differences are left, and the float
    ///    states do not contradict the confirmed integers (were those searched now, no other set may
    ///    beat them by the ratio a fix needs; when that lasts `ar_confirm_s`, all labels are given up).
    /// 3. The other candidates (new or restarted satellites) get their integers from the fixed position:
    ///    with the position and the group's clock known, a phase gives its integer directly. What lies
    ///    within a quarter cycle of an integer gets a label and starts its own confirmation.
    ///
    /// The float states are not touched: the labels only decide what is reported.
    #[allow(clippy::too_many_lines)]
    fn fix(
        &mut self,
        u: &Update,
        obs: &[Ob],
        g: &[(f64, [f64; 3], f64, f64, f64)],
        t: f64,
        sig_float: f64,
    ) -> Option<Fix> {
        let o = self.opts.clone();
        self.last_ratio = 0.0;
        let na = u.na;
        let armask = o.ar_elmask.to_radians();
        let ng = self.setup.plan.groups.len();
        // the keys that may carry an integer this epoch
        let usable = |s: &Self, q: usize| -> bool {
            let ob = &obs[q];
            s.setup.keys[ob.key].fixable && u.cix[ob.key] != usize::MAX && !ob.half
        };
        let candidate = |s: &Self, q: usize| -> bool {
            usable(s, q)
                && (s.label[obs[q].key].is_some()
                    || (obs[q].el >= armask
                        && (t - s.t_reset[obs[q].key]).abs() >= o.ar_min_lock_s))
        };
        // pairs per group; the reference is a confirmed key when there is one
        let mut pairs: Vec<(usize, usize)> = Vec::new(); // (observation, reference observation)
        for grp in 0..ng {
            let mut c: Vec<usize> = (0..obs.len())
                .filter(|&q| self.setup.keys[obs[q].key].group == grp && candidate(self, q))
                .collect();
            if c.len() < 2 {
                continue;
            }
            let rank = |q: usize| -> (u8, f64) {
                let k = obs[q].key;
                let r = if self.confirmed(k, t) {
                    2
                } else {
                    u8::from(self.label[k].is_some())
                };
                (r, obs[q].el)
            };
            c.sort_by(|a, b| {
                let (ra, rb) = (rank(*a), rank(*b));
                rb.0.cmp(&ra.0).then(rb.1.total_cmp(&ra.1))
            });
            for &q in &c[1..] {
                pairs.push((q, c[0]));
            }
        }
        let np = pairs.len();
        // float double differences of all pairs and their covariance
        let ia: Vec<usize> = pairs.iter().map(|p| u.cix[obs[p.0].key]).collect();
        let ib: Vec<usize> = pairs.iter().map(|p| u.cix[obs[p.1].key]).collect();
        let a: Vec<f64> = (0..np).map(|i| u.xc[ia[i]] - u.xc[ib[i]]).collect();
        let mut q = DMat::zeros(np);
        for i in 0..np {
            for j in 0..=i {
                let v =
                    u.pc[ia[i] * na + ia[j]] - u.pc[ia[i] * na + ib[j]] - u.pc[ib[i] * na + ia[j]]
                        + u.pc[ib[i] * na + ib[j]];
                *q.at_mut(i, j) = v;
                *q.at_mut(j, i) = v;
            }
        }
        let keys_of: Vec<(usize, usize)> =
            pairs.iter().map(|p| (obs[p.0].key, obs[p.1].key)).collect();
        let both_confirmed = |s: &Self, i: usize| {
            let (kq, kr) = keys_of[i];
            s.confirmed(kq, t) && s.confirmed(kr, t)
        };
        let n_conf = (0..np).filter(|&i| both_confirmed(self, i)).count();
        let bootstrap = n_conf < o.ar_min_dd;

        // 1. the independent search: starts the labels, then audits them
        let mut ratio = 0.0;
        if np >= o.ar_min_dd.max(1)
            && let Some((kept, ints, r)) = self.search(&a, &q, o.ar_min_dd)
        {
            ratio = r;
            for (&i, z) in kept.iter().zip(&ints) {
                let z = z.round() as i64;
                let (kq, kr) = keys_of[i];
                match (self.label[kq], self.label[kr]) {
                    (Some(lq), Some(lr)) if lq - lr == z => {
                        self.set_label(kq, lq, t);
                        if self.label_n[kr] < CONFIRM_COUNT {
                            self.set_label(kr, lr, t);
                        }
                    }
                    (Some(_), Some(lr)) => {
                        if both_confirmed(self, i) {
                            // a confirmed integer is only given up when contradicted for a while
                            if !self.dispute_t0[kq].is_finite() {
                                self.dispute_t0[kq] = t;
                            } else if (t - self.dispute_t0[kq]).abs() >= o.ar_confirm_s {
                                self.set_label(kq, z + lr, t);
                            }
                        } else if bootstrap {
                            self.set_label(kq, z + lr, t);
                        }
                    }
                    (None, Some(lr)) if bootstrap => self.set_label(kq, z + lr, t),
                    (None, None) if bootstrap => {
                        self.set_label(kr, 0, t);
                        self.set_label(kq, z, t);
                    }
                    _ => {}
                }
            }
        }
        ratio = ratio.max(self.last_ratio);
        let float_only = |ratio: f64| {
            (ratio > 0.0).then_some(Fix {
                ratio,
                dx: [0.0; 3],
                cov: [[0.0; 3]; 3],
                confirmed: false,
            })
        };

        // 2. the fixed position from the confirmed integers and this epoch's phases
        let mut solved: Option<([f64; 3], [[f64; 3]; 3], Vec<f64>)> = None;
        for _ in 0..6 {
            let held: Vec<usize> = (0..obs.len())
                .filter(|&q| usable(self, q) && self.confirmed(obs[q].key, t))
                .collect();
            let Some((dx, cov, clk, res, nd)) = self.fixed_position(obs, g, &held) else {
                break;
            };
            if nd < o.ar_min_dd {
                break;
            }
            let (worst, vmax) = res
                .iter()
                .enumerate()
                .map(|(i, v)| (i, v.abs()))
                .fold((0, 0.0), |m, x| if x.1 > m.1 { x } else { m });
            if vmax > 3.0 * o.ar_max_res {
                let k = obs[held[worst]].key;
                self.label[k] = None;
                self.dispute_t0[k] = f64::NAN;
                continue;
            }
            let rms = (res.iter().map(|v| v * v).sum::<f64>() / res.len().max(1) as f64).sqrt();
            let shift = norm3([dx[0] - u.xc[0], dx[1] - u.xc[1], dx[2] - u.xc[2]]);
            if rms <= o.ar_max_res
                && shift <= o.ar_max_shift_sigma * o.sigma_float_scale * sig_float + 0.05
            {
                solved = Some((dx, cov, clk));
            }
            break;
        }
        let Some((dx, cov, clk)) = solved else {
            return float_only(ratio);
        };

        // the audit of the held integers against the float states: were they searched now, would
        // another set beat them by the ratio a fix needs? Then they are no longer believed: the epoch
        // is not fixed, and when it lasts they are all given up
        let hs: Vec<usize> = (0..np).filter(|&i| both_confirmed(self, i)).collect();
        if !hs.is_empty() {
            let d: Vec<f64> = hs
                .iter()
                .map(|&i| {
                    let (kq, kr) = keys_of[i];
                    a[i] - (self.label[kq].unwrap_or(0) - self.label[kr].unwrap_or(0)) as f64
                })
                .collect();
            let qh = q.sub(&hs);
            let a_h: Vec<f64> = hs.iter().map(|&i| a[i]).collect();
            let beaten = qh.solve(&d).is_some_and(|w| {
                let t_held: f64 = d.iter().zip(&w).map(|(x, y)| x * y).sum();
                lambda_fix(&a_h, &qh, 1).is_ok_and(|(_, s)| {
                    s[0] > 0.0 && t_held > s[0] * (1.0 + 1e-9) && t_held / s[0] >= o.ar_ratio
                })
            });
            if beaten {
                if !self.held_dispute_t0.is_finite() {
                    self.held_dispute_t0 = t;
                } else if (t - self.held_dispute_t0).abs() >= o.ar_confirm_s {
                    self.label.fill(None);
                    self.dispute_t0.fill(f64::NAN);
                    self.held_dispute_t0 = f64::NAN;
                }
                return float_only(ratio);
            }
        }
        self.held_dispute_t0 = f64::NAN;

        // 3. the integers of the other candidates, from the fixed position
        let fresh: Vec<usize> = (0..obs.len())
            .filter(|&q| candidate(self, q) && !self.confirmed(obs[q].key, t))
            .collect();
        let pred = |q: usize| -> f64 {
            let e = g[obs[q].s].1;
            obs[q].y + e[0] * dx[0] + e[1] * dx[1] + e[2] * dx[2]
        };
        for grp in 0..ng {
            let mine: Vec<usize> = fresh
                .iter()
                .copied()
                .filter(|&q| self.setup.keys[obs[q].key].group == grp)
                .collect();
            if mine.is_empty() {
                continue;
            }
            let lam = self.setup.keys[obs[mine[0]].key].lam;
            if clk[grp].is_finite() {
                for &q in &mine {
                    let n = (pred(q) - clk[grp]) / lam;
                    if (n - n.round()).abs() <= INSTANT_MAX_FRAC {
                        self.set_label(obs[q].key, n.round() as i64, t);
                    }
                }
            } else if mine.len() >= 2 {
                // a group with no confirmed key yet: integers between its own satellites
                let r = mine
                    .iter()
                    .copied()
                    .max_by(|x, y| obs[*x].el.total_cmp(&obs[*y].el))
                    .unwrap_or(mine[0]);
                let base = self.label[obs[r].key].unwrap_or(0);
                let mut any = false;
                for &q in mine.iter().filter(|q| **q != r) {
                    let n = (pred(q) - pred(r)) / lam;
                    if (n - n.round()).abs() <= INSTANT_MAX_FRAC {
                        self.set_label(obs[q].key, n.round() as i64 + base, t);
                        any = true;
                    }
                }
                if any {
                    self.set_label(obs[r].key, base, t);
                }
            }
        }
        Some(Fix {
            ratio,
            dx,
            cov,
            confirmed: true,
        })
    }
}

fn vel_enu(pe: &PreparedEpoch, r: &[[f64; 3]; 3]) -> [f32; 3] {
    if !pe.vel[0].is_finite() {
        return [f32::NAN; 3];
    }
    let v = mat3_vec(
        r,
        [
            f64::from(pe.vel[0]),
            f64::from(pe.vel[1]),
            f64::from(pe.vel[2]),
        ],
    );
    [v[0] as f32, v[1] as f32, v[2] as f32]
}
