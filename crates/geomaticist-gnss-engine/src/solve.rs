//! The static baseline solver: common epochs, satellites and pairs, the per-station preparation,
//! corrections, slips and arcs, the float solution, then (in [`crate::fixing`]) the integer fixing and
//! the checks, the trust test and the realistic σ.
//!
//! Operation order is kept fixed so the stage dumps agree to the tolerance tiers.

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::obs::{
    BandSignals, Combo, combine, common_signals, if_coeffs, pick_pair,
};
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};
use geomaticist_gnss_products::{AntennaCal, Antex, receiver_correction};
use serde_json::{Map, Value};

use crate::arr::{B2, F2, I2, Vec3Arr};
use crate::dump::{jf, jfl, jround};
use crate::linalg::{DMat, add3, dot, mat3_vec, norm3, rotate_cov3, sub3};
use crate::models::{azel, ecef_to_llh, enu_matrix, geometry, sat_pco_offset, sun_ecef};
use crate::nputil::{
    argmax, mad_scatter, mean, median, nanmedian, np_round, np_sum, py_round_i, sample_rows,
    searchsorted,
};
use crate::options::{
    Options, StaticResult, TideModel, VtpvForm, antex_name, glo_fixable, realistic_cov,
    sigma_quality,
};
use crate::orbits::{OrbitSource, Orbits};
use crate::slips::{self, SlipEvents, find_step};
use crate::stages::{Acc, ArcsStage, FloatStage, ObservationsStage, Stages, StationEvents};
use crate::tropo::{TropoInputs, troposphere};

/// A failed solve (in `solve` / `spp`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EngineError(pub String);

impl std::fmt::Display for EngineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for EngineError {}

/// Index of the nearest of the sorted times `ts` to each of `t`.
#[must_use]
pub fn nearest(ts: &[f64], t: &[f64]) -> Vec<usize> {
    if ts.len() < 2 {
        return vec![0; t.len()];
    }
    t.iter()
        .map(|&x| {
            let j = searchsorted(ts, x).clamp(1, ts.len() - 1);
            if (ts[j - 1] - x).abs() <= (ts[j] - x).abs() {
                j - 1
            } else {
                j
            }
        })
        .collect()
}

/// One station on the common epochs.
#[derive(Debug, Clone)]
pub struct Station {
    /// Ionosphere-free phase (m).
    pub phase: F2,
    /// Ionosphere-free code (m).
    pub code: F2,
    /// Phase + code + orbit usable.
    pub ok: B2,
    /// LLI or hard break since the previous common epoch.
    pub slip: B2,
    /// Geometry-free L1 − L2 (m).
    pub gf: F2,
    /// Satellite positions.
    pub sat_xyz: Vec3Arr,
    /// Satellite clocks (s).
    pub sat_clk: F2,
    /// Tide per epoch.
    pub tide: Vec<[f64; 3]>,
    /// Preparation events.
    pub events: StationEvents,
    /// Melbourne-Wübbena (cycles).
    pub mw: F2,
    /// Loss of lock only.
    pub lli: B2,
    /// Whole-millisecond phase jump or power cycle.
    pub hard: B2,
    /// A whole-millisecond phase step put back here.
    pub ms: B2,
    /// Half-cycle flag at this epoch.
    pub half: B2,
}

/// The range / direction / troposphere / antenna terms of one station.
#[derive(Debug, Clone)]
pub struct Corrections {
    /// Range with Sagnac (m).
    pub rho: F2,
    /// Unit vector receiver -> satellite.
    pub e: Vec3Arr,
    /// Azimuth (rad).
    pub az: F2,
    /// Elevation (rad).
    pub el: F2,
    /// Hydrostatic delay × mapping (m).
    pub dry: F2,
    /// Wet mapping.
    pub map: F2,
    /// Ionosphere-free antenna correction (m).
    pub ant: F2,
    /// Its geometry-free part (band 1 − band 2).
    pub antgf: F2,
}

fn sat_pairs_of(
    st: &StationObs,
    cols: &[Sat],
    pairs: &BTreeMap<Sat, (u8, u8)>,
    fcn: &BTreeMap<Sat, i32>,
    signals: Option<&BTreeMap<Sat, BandSignals>>,
) -> BTreeMap<usize, Combo> {
    let mut combos = BTreeMap::new();
    for (j, sat) in cols.iter().enumerate() {
        if !st.sats.contains_key(sat) {
            continue;
        }
        let sg = signals.and_then(|s| s.get(sat));
        if let Ok(c) = combine(st, *sat, pairs[sat], fcn.get(sat).copied(), sg) {
            combos.insert(j, c);
        }
    }
    combos
}

/// One station on the common epochs.
pub fn prepare(
    st: &StationObs,
    t: &[f64],
    cols: &[Sat],
    pairs: &BTreeMap<Sat, (u8, u8)>,
    fcn: &BTreeMap<Sat, i32>,
    orb: &Orbits,
    sun: &[[f64; 3]],
    opts: &Options,
    pos: [f64; 3],
    sat_cals: &BTreeMap<Sat, Option<AntennaCal>>,
    signals: Option<&BTreeMap<Sat, BandSignals>>,
) -> Station {
    let (ne, ns) = (t.len(), cols.len());
    let ii = nearest(&st.t, t);
    let tol = opts.epoch_tol.max(1e-3);
    let here: Vec<bool> = ii
        .iter()
        .zip(t)
        .map(|(&i, &x)| st.t.get(i).is_some_and(|v| (v - x).abs() < tol))
        .collect();
    let t_own: Vec<f64> = ii
        .iter()
        .zip(t)
        .zip(&here)
        .map(|((&i, &x), &h)| if h { st.t[i] } else { x })
        .collect();
    let mut phase = F2::nans(ne, ns);
    let mut code = F2::nans(ne, ns);
    let mut c1s = F2::nans(ne, ns);
    let mut slip = B2::falses(ne, ns);
    let mut gf = F2::nans(ne, ns);
    let mut mw = F2::nans(ne, ns);
    let mut ok = B2::falses(ne, ns);
    let mut half = B2::falses(ne, ns);
    let pf: Vec<bool> = ii
        .iter()
        .zip(&here)
        .map(|(&i, &h)| h && st.power_fail.get(i).copied().unwrap_or(false))
        .collect();
    let combos = sat_pairs_of(st, cols, pairs, fcn, signals);
    for (&j, c) in &combos {
        // loss of lock between common epochs is carried to the next common epoch
        let mut cs = vec![0i64; c.lli.len() + 1];
        for (k, &l) in c.lli.iter().enumerate() {
            cs[k + 1] = cs[k] + i64::from(l);
        }
        // cs[k+1] is the cumulative sum at k
        for e in 0..ne {
            let i = ii[e];
            let carried = if e == 0 {
                cs[i + 1] > 0
            } else {
                let p = ii[e - 1];
                cs[i + 1] - cs[p + 1] > 0
            };
            slip.set(e, j, (carried || c.lli[i]) && here[e]);
            half.set(e, j, c.half[i] && here[e]);
            if here[e] {
                phase.set(e, j, c.phase_if[i]);
                gf.set(e, j, c.gf[i]);
                mw.set(e, j, c.mw[i]);
                code.set(e, j, c.code_if[i]);
                c1s.set(e, j, c.code1[i]);
            }
            ok.set(e, j, phase.at(e, j).is_finite() && c1s.at(e, j).is_finite());
        }
    }
    let sat_at = |j: usize, prange: &[f64]| -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        let sat = cols[j];
        let tr: Vec<f64> = t_own
            .iter()
            .zip(prange)
            .map(|(&x, &p)| x - (if p.is_finite() { p } else { 0.0 }) / CLIGHT)
            .collect();
        let (_x, clk1, o1, _) = orb.position(sat, &tr);
        let tr2: Vec<f64> = tr.iter().zip(&clk1).map(|(a, b)| a - b).collect();
        let (mut xyz, clk, o2, com) = orb.position(sat, &tr2);
        if com && sat_cals.contains_key(&sat) {
            let c = &combos[&j];
            let (a, b) = if_coeffs(c.f1, c.f2);
            let cal = sat_cals[&sat].as_ref();
            let (b1, b2) = pairs[&sat];
            for e in 0..ne {
                let off = sat_pco_offset(cal, sat.sys.letter() as u8, b1, b2, xyz[e], sun[e], a, b);
                xyz[e] = add3(xyz[e], off);
            }
        }
        let okv: Vec<bool> = o1.iter().zip(&o2).map(|(a, b)| *a && *b).collect();
        (xyz, clk, okv)
    };
    let mut sxyz = Vec3Arr::zeros(ne, ns);
    let mut sclk = F2::zeros(ne, ns);
    for &j in combos.keys() {
        let (xyz, clk, o) = sat_at(j, &c1s.col(j));
        for e in 0..ne {
            sxyz.set(e, j, xyz[e]);
            sclk.set(e, j, clk[e]);
            if !o[e] {
                ok.set(e, j, false);
            }
        }
    }
    // code sanity against the epoch's receiver clock (median over satellites)
    let mut events = StationEvents::default();
    let mut rho = F2::zeros(ne, ns);
    for e in 0..ne {
        for j in 0..ns {
            rho.set(e, j, geometry(sxyz.at(e, j), pos).0);
        }
    }
    let mut dev = F2::nans(ne, ns);
    for e in 0..ne {
        for j in 0..ns {
            if ok.at(e, j) {
                dev.set(e, j, c1s.at(e, j) - rho.at(e, j) + CLIGHT * sclk.at(e, j));
            }
        }
    }
    let n_fin: Vec<usize> = (0..ne)
        .map(|e| dev.row(e).iter().filter(|v| v.is_finite()).count())
        .collect();
    let rows: Vec<bool> = n_fin.iter().map(|&n| n >= 3).collect();
    let clk_r: Vec<f64> = (0..ne)
        .map(|e| {
            if rows[e] {
                nanmedian(dev.row(e))
            } else {
                f64::NAN
            }
        })
        .collect();
    let mut bad = B2::falses(ne, ns);
    for e in 0..ne {
        for j in 0..ns {
            let d = dev.at(e, j) - clk_r[e];
            bad.set(
                e,
                j,
                ok.at(e, j) && rows[e] && d.is_finite() && d.abs() > opts.code_gross_tx,
            );
        }
    }
    for j in 0..ns {
        if !bad.any_col(j) {
            continue;
        }
        let pr: Vec<f64> = (0..ne)
            .map(|e| {
                if bad.at(e, j) {
                    rho.at(e, j) - CLIGHT * sclk.at(e, j) + clk_r[e]
                } else {
                    c1s.at(e, j)
                }
            })
            .collect();
        let (xyz, clk, _o) = sat_at(j, &pr);
        for e in 0..ne {
            if bad.at(e, j) {
                sxyz.set(e, j, xyz[e]);
                sclk.set(e, j, clk[e]);
                code.set(e, j, f64::NAN);
                mw.set(e, j, f64::NAN);
                events.code_ms += 1;
            }
        }
    }
    for e in 0..ne {
        if !(rows[e] || n_fin[e] > 0) {
            for j in 0..ns {
                ok.set(e, j, false);
            }
        }
    }
    // phase minus code: whole-millisecond phase jumps
    let lli = slip.clone();
    let mut hard = B2::falses(ne, ns);
    let mut msk = B2::falses(ne, ns);
    let ms = CLIGHT * 1e-3;
    for j in 0..ns {
        let have: Vec<usize> = (0..ne)
            .filter(|&e| (phase.at(e, j) - code.at(e, j)).is_finite())
            .collect();
        if have.len() > 1 {
            let pc: Vec<f64> = have
                .iter()
                .map(|&e| phase.at(e, j) - code.at(e, j))
                .collect();
            let dpc: Vec<f64> = pc.windows(2).map(|w| w[1] - w[0]).collect();
            let mut jump: Vec<bool> = dpc.iter().map(|d| d.abs() > opts.slip_pc).collect();
            if opts.ms_repair {
                let (f1, f2) = (combos[&j].f1, combos[&j].f2);
                for q in 0..dpc.len() {
                    let k_ms = np_round(dpc[q] / ms, 0);
                    let whole = jump[q] && k_ms != 0.0 && (dpc[q] - k_ms * ms).abs() < opts.ms_tol;
                    if whole {
                        let e0 = have[q + 1];
                        for e in e0..ne {
                            let v = phase.at(e, j) - k_ms * ms;
                            phase.set(e, j, v);
                            let m = mw.at(e, j) - k_ms * ms * (f1 - f2) / CLIGHT;
                            mw.set(e, j, m);
                        }
                        msk.set(e0, j, true);
                        events.ms_repaired += 1;
                        jump[q] = false;
                    }
                }
            }
            for q in 0..dpc.len() {
                if jump[q] {
                    hard.set(have[q + 1], j, true);
                    events.pc_jump += 1;
                }
            }
        }
    }
    for e in 0..ne {
        if pf[e] {
            for j in 0..ns {
                hard.set(e, j, true);
            }
        }
    }
    let slip = lli.or(&hard);
    Station {
        phase,
        code,
        ok,
        slip,
        gf,
        sat_xyz: sxyz,
        sat_clk: sclk,
        tide: vec![[0.0; 3]; ne],
        events,
        mw,
        lli,
        hard,
        ms: msk,
        half,
    }
}

/// The station antenna's calibration (its own entry, else the same antenna with no radome).
pub fn receiver_cal(
    antex: Option<&mut Antex>,
    st: &StationObs,
    notes: &mut Vec<String>,
) -> Option<AntennaCal> {
    let atx = antex?;
    if st.antenna.is_empty() {
        return None;
    }
    let (entry, why) = antex_name(&st.antenna, &atx.types());
    let cal = if entry.is_empty() {
        None
    } else {
        atx.receiver(&entry)
    };
    // an entry with frequency rows but no ZEN1 / ZEN2 / DZEN gives zero PCV; the engine says so instead of staying silent
    if let Some(c) = &cal
        && c.zen.is_empty()
        && !c.noazi.is_empty()
    {
        notes.push(format!(
            "{}: ANTEX entry {:?} has no ZEN1 / ZEN2 / DZEN line: its PCV is taken as zero",
            st.name,
            c.name.trim()
        ));
    }
    if cal.is_none() {
        notes.push(format!(
            "{}: antenna {:?} not in ANTEX — no PCO / PCV",
            st.name,
            st.antenna.trim()
        ));
    } else if !why.is_empty() {
        notes.push(format!("{}: {why}", st.name));
    }
    cal
}

/// Whether system `sys`'s bands have their own receiver calibration at both stations, or the antennas
/// are the same type (B7: no entry at all disqualifies).
#[must_use]
pub fn system_calibrated(
    cals: &[Option<&AntennaCal>],
    antennas: &[&str],
    sys: char,
    bands: &[u8],
) -> bool {
    let norm = |a: &str| a.split_whitespace().collect::<Vec<_>>().join(" ");
    let set: BTreeSet<String> = antennas.iter().map(|a| norm(a)).collect();
    if set.len() == 1 && !antennas[0].is_empty() {
        return true;
    }
    for cal in cals {
        let Some(cal) = cal else { return false };
        if bands
            .iter()
            .any(|b| !cal.pco.contains_key(&format!("{sys}{b:02}")))
        {
            return false;
        }
    }
    true
}

/// Range, LOS, az / el, hydrostatic delay, wet mapping and the antenna corrections at positions
/// `pos_ep` (one per epoch).
pub fn corrections(
    sta: &Station,
    pos_ep: &[[f64; 3]],
    cols: &[Sat],
    pairs: &BTreeMap<Sat, (u8, u8)>,
    fcn: &BTreeMap<Sat, i32>,
    cal: Option<&AntennaCal>,
    t: &[f64],
    opts: &Options,
    tropo: &TropoInputs,
) -> Corrections {
    let (ne, ns) = (sta.sat_xyz.ne, sta.sat_xyz.ns);
    let mut rho = F2::zeros(ne, ns);
    let mut e = Vec3Arr::zeros(ne, ns);
    let mut az = F2::zeros(ne, ns);
    let mut el = F2::zeros(ne, ns);
    let (lat, lon, h) = ecef_to_llh(pos_ep[0]);
    let r = enu_matrix(lat, lon);
    for ep in 0..ne {
        for j in 0..ns {
            let (rr, ee) = geometry(sta.sat_xyz.at(ep, j), pos_ep[ep]);
            rho.set(ep, j, rr);
            e.set(ep, j, ee);
            let (a, l, _) = azel(ee, &r);
            az.set(ep, j, a);
            el.set(ep, j, l);
        }
    }
    let trop = troposphere(&opts.trop_model, lat, lon, h, t, &el, tropo);
    let mut dry = F2::zeros(ne, ns);
    for ep in 0..ne {
        for j in 0..ns {
            dry.set(ep, j, trop.zhd[ep] * trop.mh.at(ep, j));
        }
    }
    let mut ant = F2::zeros(ne, ns);
    let mut antgf = F2::zeros(ne, ns);
    if let Some(cal) = cal {
        for (j, sat) in cols.iter().enumerate() {
            let (b1, b2) = pairs[sat];
            let f = fcn.get(sat).copied();
            let (Ok(f1), Ok(f2)) = (freq(*sat, b1, f), freq(*sat, b2, f)) else {
                continue;
            };
            let (a, b) = if_coeffs(f1, f2);
            let azc = az.col(j);
            let elc = el.col(j);
            let sys = sat.sys.letter() as u8;
            let c1 = receiver_correction(Some(cal), sys, b1, &azc, &elc, opts.pcv_azimuth);
            let c2 = receiver_correction(Some(cal), sys, b2, &azc, &elc, opts.pcv_azimuth);
            for ep in 0..ne {
                ant.set(ep, j, a * c1[ep] - b * c2[ep]);
                antgf.set(ep, j, c1[ep] - c2[ep]);
            }
        }
    }
    Corrections {
        rho,
        e,
        az,
        el,
        dry,
        map: trop.mw,
        ant,
        antgf,
    }
}

/// Phase variance per station and satellite.
#[must_use]
pub fn variance(el: &F2, sysc: &[char], opts: &Options) -> F2 {
    let mut out = F2::zeros(el.ne, el.ns);
    for j in 0..el.ns {
        let fact = opts.efact_of(sysc[j]) * 3.0;
        for e in 0..el.ne {
            let s = el.at(e, j).max(1.0f64.to_radians()).sin();
            out.set(
                e,
                j,
                fact * fact * (opts.err_a * opts.err_a + opts.err_b * opts.err_b / (s * s)),
            );
        }
    }
    out
}

/// Code-only position (ECEF) of a static receiver: median of epoch solutions.
pub fn spp(
    st: &StationObs,
    orb: &Orbits,
    systems: &str,
    fcn: &BTreeMap<Sat, i32>,
    every: usize,
) -> Result<[f64; 3], EngineError> {
    let out: Vec<[f64; 3]> = spp_track(st, orb, systems, fcn, every)
        .into_iter()
        .map(|x| x.1)
        .collect();
    if out.is_empty() {
        return Err(EngineError(format!(
            "{}: no single-point solution",
            st.name
        )));
    }
    Ok([
        median(&out.iter().map(|x| x[0]).collect::<Vec<_>>()),
        median(&out.iter().map(|x| x[1]).collect::<Vec<_>>()),
        median(&out.iter().map(|x| x[2]).collect::<Vec<_>>()),
    ])
}

/// A single-point solution moved more than this from the session's median (m) counts as a move.
pub const MOVING_M: f64 = 25.0;
/// The share of sampled epochs that must lie farther than [`MOVING_M`] for the receiver to have moved.
pub const MOVING_SHARE: f64 = 0.2;

/// How far a receiver's code-only positions wander over a session.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wander {
    /// Epochs sampled with a solution.
    pub n: usize,
    /// Share of them farther than [`MOVING_M`] from the median position.
    pub share_far: f64,
    /// The largest distance from the median, m.
    pub max_m: f64,
}

impl Wander {
    /// The receiver moved: a kinematic log, which a static baseline can't use.
    #[must_use]
    pub fn moved(&self) -> bool {
        self.n >= 10 && self.share_far > MOVING_SHARE
    }
}

/// Code-only positions of up to ~150 epochs spread over the session, and how far they lie from their median
/// (`None`: too few solutions to tell). A static receiver's stay within a few metres; a rover carried between
/// points spreads over the distance it walked. Cheap (well under a second), so a moving log is refused before
/// the solver spends minutes on it.
#[must_use]
pub fn wander(
    st: &StationObs,
    orb: &Orbits,
    systems: &str,
    fcn: &BTreeMap<Sat, i32>,
) -> Option<Wander> {
    let every = (st.t.len() / 150).max(1);
    let out: Vec<[f64; 3]> = spp_track(st, orb, systems, fcn, every)
        .into_iter()
        .map(|x| x.1)
        .collect();
    if out.len() < 10 {
        return None;
    }
    let m = [
        median(&out.iter().map(|x| x[0]).collect::<Vec<_>>()),
        median(&out.iter().map(|x| x[1]).collect::<Vec<_>>()),
        median(&out.iter().map(|x| x[2]).collect::<Vec<_>>()),
    ];
    let d: Vec<f64> = out.iter().map(|x| norm3(sub3(*x, m))).collect();
    let far = d.iter().filter(|v| **v > MOVING_M).count();
    Some(Wander {
        n: d.len(),
        share_far: far as f64 / d.len() as f64,
        max_m: d.iter().copied().fold(0.0, f64::max),
    })
}

/// The epoch solutions behind [`spp`]: (epoch index, ECEF) of every `every`-th epoch that solved.
#[must_use]
pub fn spp_track(
    st: &StationObs,
    orb: &Orbits,
    systems: &str,
    fcn: &BTreeMap<Sat, i32>,
    every: usize,
) -> Vec<(usize, [f64; 3])> {
    let sats: Vec<Sat> = st
        .sats
        .keys()
        .copied()
        .filter(|s| {
            systems.contains(s.sys.letter())
                && orb.has(*s)
                && (s.sys != GnssSystem::Glonass || fcn.contains_key(s))
        })
        .collect();
    let idx: Vec<usize> = (0..st.t.len()).step_by(every.max(1)).collect();
    let t: Vec<f64> = idx.iter().map(|&i| st.t[i]).collect();
    let n = idx.len();
    let mut p = vec![vec![f64::NAN; sats.len()]; n];
    let mut xs = vec![vec![[0.0; 3]; sats.len()]; n];
    for (j, &s) in sats.iter().enumerate() {
        let Some(pr) = pick_pair(s, &[st]) else {
            continue;
        };
        let Ok(c) = combine(st, s, pr, fcn.get(&s).copied(), None) else {
            continue;
        };
        let tr: Vec<f64> = idx
            .iter()
            .zip(&t)
            .map(|(&i, &tt)| {
                tt - (if c.code1[i].is_finite() {
                    c.code1[i]
                } else {
                    0.0
                }) / CLIGHT
            })
            .collect();
        let (_x, clk, o, _) = orb.position(s, &tr);
        let tr2: Vec<f64> = tr.iter().zip(&clk).map(|(a, b)| a - b).collect();
        let (xyz, clk2, o2, _) = orb.position(s, &tr2);
        for (k, &i) in idx.iter().enumerate() {
            if o[k] && o2[k] && c.code_if[i].is_finite() {
                p[k][j] = c.code_if[i] + CLIGHT * clk2[k];
            }
            xs[k][j] = xyz[k];
        }
    }
    let sysi: Vec<usize> = sats
        .iter()
        .map(|s| "GRECJ".find(s.sys.letter()).unwrap_or(0))
        .collect();
    let mut out: Vec<(usize, [f64; 3])> = Vec::new();
    let x0 = if norm3(st.approx) > 6e6 {
        st.approx
    } else {
        [0.0; 3]
    };
    for k in 0..n {
        let m: Vec<usize> = (0..sats.len()).filter(|&j| p[k][j].is_finite()).collect();
        if m.len() < 6 {
            continue;
        }
        let mut xk = [0.0; 8];
        xk[..3].copy_from_slice(&x0);
        let mut v = vec![0.0; m.len()];
        for _ in 0..8 {
            let mut hrows: Vec<Vec<f64>> = Vec::with_capacity(m.len());
            for (r, &j) in m.iter().enumerate() {
                let (rho, e) = geometry(xs[k][j], [xk[0], xk[1], xk[2]]);
                let mut row = vec![0.0; 8];
                row[0] = -e[0];
                row[1] = -e[1];
                row[2] = -e[2];
                row[3 + sysi[j]] = 1.0;
                v[r] = p[k][j] - (rho + xk[3 + sysi[j]]);
                hrows.push(row);
            }
            let used: Vec<usize> = (3..8)
                .filter(|&c| hrows.iter().any(|r| r[c] != 0.0))
                .collect();
            let cols_u: Vec<usize> = (0..3).chain(used.iter().copied()).collect();
            // least squares by normal equations
            let nu = cols_u.len();
            let mut nmat = DMat::zeros(nu);
            let mut bvec = vec![0.0; nu];
            for (r, row) in hrows.iter().enumerate() {
                for (a, &ca) in cols_u.iter().enumerate() {
                    bvec[a] += row[ca] * v[r];
                    for (b, &cb) in cols_u.iter().enumerate() {
                        nmat.add(a, b, row[ca] * row[cb]);
                    }
                }
            }
            let Some(dx) = nmat.solve(&bvec) else { break };
            for (a, &ca) in cols_u.iter().enumerate() {
                xk[ca] += dx[a];
            }
            if norm3([dx[0], dx[1], dx[2]]) < 1e-3 {
                break;
            }
        }
        let rms = (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
        if rms < 30.0 {
            out.push((idx[k], [xk[0], xk[1], xk[2]]));
        }
    }
    out
}

/// Phase observations to drop: steps over `phase_gross` in an arc's SD phase less the model, with the
/// epoch-to-epoch change common to all satellites out.
#[must_use]
pub fn phase_steps(r: &F2, ok: &B2, arc: &I2, t: &[f64], opts: &Options) -> B2 {
    let (ne, ns) = (r.ne, r.ns);
    let mut bad = B2::falses(ne, ns);
    if ne < 2 {
        return bad;
    }
    let rr = r.zip_map(ok, |v, o| if o { v } else { f64::NAN });
    let mut cm = vec![0.0; ne - 1];
    for e in 0..ne - 1 {
        let mut d = Vec::new();
        for j in 0..ns {
            if ok.at(e + 1, j) && ok.at(e, j) && arc.at(e + 1, j) == arc.at(e, j) {
                let v = rr.at(e + 1, j) - rr.at(e, j);
                if v.is_finite() {
                    d.push(v);
                }
            }
        }
        if d.len() >= 2 {
            cm[e] = median(&d);
        }
    }
    let mut cum = vec![0.0; ne];
    for e in 1..ne {
        cum[e] = cum[e - 1] + cm[e - 1];
    }
    for j in 0..ns {
        let mut ids: Vec<i64> = (0..ne)
            .filter(|&e| ok.at(e, j))
            .map(|e| arc.at(e, j))
            .collect();
        ids.sort_unstable();
        ids.dedup();
        for a in ids {
            let ii: Vec<usize> = (0..ne)
                .filter(|&e| ok.at(e, j) && arc.at(e, j) == a)
                .collect();
            if ii.len() < 2 {
                continue;
            }
            let vals: Vec<f64> = ii.iter().map(|&e| rr.at(e, j) - cum[e]).collect();
            let step: Vec<bool> = (1..ii.len())
                .map(|k| {
                    (vals[k] - vals[k - 1]).abs()
                        > opts.phase_gross + opts.phase_gross_rate * (t[ii[k]] - t[ii[k - 1]])
                })
                .collect();
            if !step.iter().any(|&s| s) {
                continue;
            }
            let mut seg = vec![0usize; ii.len()];
            for k in 1..ii.len() {
                seg[k] = seg[k - 1] + usize::from(step[k - 1]);
            }
            let nseg = seg[ii.len() - 1] + 1;
            let mut counts = vec![0usize; nseg];
            for &s in &seg {
                counts[s] += 1;
            }
            let keep = counts
                .iter()
                .enumerate()
                .fold(0usize, |b, (i, &c)| if c > counts[b] { i } else { b });
            for (k, &e) in ii.iter().enumerate() {
                if seg[k] != keep {
                    bad.set(e, j, true);
                }
            }
        }
    }
    bad
}

/// Satellites whose SD code is biased against the others'. `{column: b (m)}`.
#[must_use]
pub fn code_biases(
    res_cd: &F2,
    ok: &B2,
    t: &[f64],
    sysc: &[char],
    opts: &Options,
) -> BTreeMap<usize, f64> {
    let mut out = BTreeMap::new();
    if opts.code_bias_m <= 0.0 {
        return out;
    }
    let (ne, _ns) = (res_cd.ne, res_cd.ns);
    let systems: BTreeSet<char> = sysc.iter().copied().collect();
    for s in systems {
        let k: Vec<usize> = (0..sysc.len()).filter(|&j| sysc[j] == s).collect();
        if k.len() < 4 {
            continue;
        }
        let v = |e: usize, j: usize| {
            if ok.at(e, j) && res_cd.at(e, j).is_finite() {
                res_cd.at(e, j)
            } else {
                f64::NAN
            }
        };
        let rows: Vec<usize> = (0..ne)
            .filter(|&e| k.iter().filter(|&&j| v(e, j).is_finite()).count() >= 3)
            .collect();
        if rows.len() < opts.code_bias_n {
            continue;
        }
        let meds: Vec<f64> = rows
            .iter()
            .map(|&e| nanmedian(&k.iter().map(|&j| v(e, j)).collect::<Vec<_>>()))
            .collect();
        let mut b = vec![f64::NAN; k.len()];
        let mut se = vec![f64::NAN; k.len()];
        for (i, &j) in k.iter().enumerate() {
            let mut vals = Vec::new();
            let mut tt = Vec::new();
            for (r, &e) in rows.iter().enumerate() {
                let x = v(e, j) - meds[r];
                if x.is_finite() {
                    vals.push(x);
                    tt.push(t[e]);
                }
            }
            if vals.len() < opts.code_bias_n {
                continue;
            }
            b[i] = mean(&vals);
            let n_eff = ((tt[tt.len() - 1] - tt[0]) / opts.wl_corr_s + 1.0).min(vals.len() as f64);
            se[i] = crate::nputil::std_ddof(&vals, 0) / n_eff.sqrt();
        }
        let have: Vec<bool> = b.iter().map(|x| x.is_finite()).collect();
        let n_have = have.iter().filter(|&&h| h).count();
        if n_have < 4 {
            continue;
        }
        let mut flag = vec![false; k.len()];
        for _round in 0..2 {
            let good: Vec<f64> = (0..k.len())
                .filter(|&i| have[i] && !flag[i])
                .map(|i| b[i])
                .collect();
            let med = median(&good);
            let spread = mad_scatter(&good).max(0.3);
            let mut new = vec![false; k.len()];
            for i in 0..k.len() {
                let off = (if have[i] { b[i] } else { med } - med).abs();
                let z =
                    off / (spread * spread + (if have[i] { se[i] } else { 1.0 }).powi(2)).sqrt();
                new[i] = have[i] && off > opts.code_bias_m && z > opts.code_bias_k;
            }
            if 2 * new.iter().filter(|&&x| x).count() >= n_have {
                new = vec![false; k.len()];
            }
            if new == flag {
                break;
            }
            flag = new;
        }
        let rest: Vec<f64> = (0..k.len())
            .filter(|&i| have[i] && !flag[i])
            .map(|i| b[i])
            .collect();
        let m_rest = median(&rest);
        for i in 0..k.len() {
            if flag[i] {
                out.insert(
                    k[i],
                    geomaticist_gnss_core::text::round_to(b[i] - m_rest, 2),
                );
            }
        }
    }
    out
}

/// The columns the code-only check uses.
#[must_use]
pub fn code_check_columns(
    sysc: &[char],
    have_cd: &[bool],
    biased: &BTreeMap<usize, f64>,
    opts: &Options,
) -> Vec<bool> {
    let mut in_c = vec![true; sysc.len()];
    for &j in biased.keys() {
        in_c[j] = false;
    }
    if opts.code_check_glo > 0 {
        let n = (0..sysc.len())
            .filter(|&j| have_cd[j] && in_c[j] && sysc[j] != 'R')
            .count();
        if n >= opts.code_check_glo {
            for j in 0..sysc.len() {
                if sysc[j] == 'R' {
                    in_c[j] = false;
                }
            }
        }
    }
    in_c
}

/// Horizontal σ (m) of a code-only position counting each satellite once.
#[must_use]
pub fn code_sat_sigma(
    used: &[(usize, Vec<usize>)],
    rc: &F2,
    eb: &Vec3Arr,
    d: [f64; 3],
    sysc: &[char],
    ns: usize,
    re: &[[f64; 3]; 3],
    floor: f64,
) -> Option<f64> {
    let mut sm = vec![0.0; ns];
    let mut sn = vec![0.0; ns];
    let mut sd = vec![[0.0; 3]; ns];
    for (e, sel) in used {
        let v: Vec<f64> = sel
            .iter()
            .map(|&j| rc.at(*e, j) + crate::linalg::dot3(eb.at(*e, j), d))
            .collect();
        let mv = mean(&v);
        for (k, &j) in sel.iter().enumerate() {
            sm[j] += v[k] - mv;
            sn[j] += 1.0;
            let ee = eb.at(*e, j);
            for i in 0..3 {
                sd[j][i] += ee[i];
            }
        }
    }
    let have: Vec<usize> = (0..ns).filter(|&j| sn[j] >= 10.0).collect();
    let systems: BTreeSet<char> = have.iter().map(|&j| sysc[j]).collect();
    let dof = have.len() as i64 - 3 - systems.len() as i64;
    if dof < 2 {
        return None;
    }
    let mut m: Vec<f64> = have.iter().map(|&j| sm[j] / sn[j]).collect();
    let mut em: Vec<[f64; 3]> = have
        .iter()
        .map(|&j| [sd[j][0] / sn[j], sd[j][1] / sn[j], sd[j][2] / sn[j]])
        .collect();
    for s in systems {
        let k: Vec<usize> = (0..have.len()).filter(|&i| sysc[have[i]] == s).collect();
        let mk = mean(&k.iter().map(|&i| m[i]).collect::<Vec<_>>());
        let mut me = [0.0; 3];
        for c in 0..3 {
            me[c] = mean(&k.iter().map(|&i| em[i][c]).collect::<Vec<_>>());
        }
        for &i in &k {
            m[i] -= mk;
            for c in 0..3 {
                em[i][c] -= me[c];
            }
        }
    }
    let mut ete = [[0.0; 3]; 3];
    for row in &em {
        for a in 0..3 {
            for b in 0..3 {
                ete[a][b] += row[a] * row[b];
            }
        }
    }
    let inv = crate::linalg::inv3(&ete)?;
    let scale = (dot(&m, &m) / dof as f64).max(floor * floor);
    let mut c = inv;
    for row in &mut c {
        for v in row.iter_mut() {
            *v *= scale;
        }
    }
    let ch = rotate_cov3(re, &c);
    Some((ch[0][0] + ch[1][1]).max(0.0).sqrt())
}

/// A group's SD design rows, row-major `nr × u`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Rows {
    /// Rows (SD observations of the group).
    pub nr: usize,
    /// Columns (the group's parameters).
    pub u: usize,
    /// Elements, row-major.
    pub a: Vec<f64>,
}

impl Rows {
    /// All zeros.
    #[must_use]
    pub fn zeros(nr: usize, u: usize) -> Self {
        Self {
            nr,
            u,
            a: vec![0.0; nr * u],
        }
    }

    /// Reshaped to `nr × u`, all zeros, keeping the allocation.
    pub fn reset(&mut self, nr: usize, u: usize) {
        self.nr = nr;
        self.u = u;
        self.a.clear();
        self.a.resize(nr * u, 0.0);
    }

    /// Filled from `f(r, c)`.
    #[must_use]
    pub fn from_fn(nr: usize, u: usize, f: impl Fn(usize, usize) -> f64) -> Self {
        let mut a = Vec::with_capacity(nr * u);
        for r in 0..nr {
            for c in 0..u {
                a.push(f(r, c));
            }
        }
        Self { nr, u, a }
    }

    /// Row `r`.
    #[inline]
    #[must_use]
    pub fn row(&self, r: usize) -> &[f64] {
        &self.a[r * self.u..(r + 1) * self.u]
    }
}

impl std::ops::Index<(usize, usize)> for Rows {
    type Output = f64;

    #[inline]
    fn index(&self, (r, c): (usize, usize)) -> &f64 {
        &self.a[r * self.u + c]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Rows {
    #[inline]
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut f64 {
        &mut self.a[r * self.u + c]
    }
}

/// One epoch / system / kind's contribution to the normal equations: the dense block on the parameters `g`
/// (`AᵀWA − (AᵀW1)(AᵀW1)ᵀ/1ᵀW1` and the matching right-hand side), with the σ0 and RMS accumulators.
#[derive(Debug, Clone)]
pub struct Block {
    /// Parameter indices.
    pub g: Vec<usize>,
    /// `u × u` block, row-major.
    pub m: Vec<f64>,
    /// Right-hand side.
    pub b: Vec<f64>,
    /// `vᵀPv` of the group.
    pub vtpv: f64,
    /// DDs in the group (`n − 1`).
    pub nrow: usize,
    /// `yᵀy` (for the RMS of `key`).
    pub ss: f64,
    /// `'p'` phase, `'c'` code, `None` for the ionosphere-weighted solution and the code check.
    pub key: Option<char>,
}

/// Per-task scratch of the assembly: one group's design (`asd`, `var`, `y`, the row selection) and the
/// block's temporaries (`w`, `ca`), reused group after group. Without it every block cost about nine
/// allocations — several hundred thousand per build, and the heap lock is what eight solves at once then
/// wait on.
#[derive(Debug, Clone, Default)]
#[allow(missing_docs)]
pub struct Scratch {
    pub asd: Rows,
    pub asd_full: Rows,
    pub sel: Vec<usize>,
    pub var: Vec<f64>,
    pub y: Vec<f64>,
    pub tmp: Vec<f64>,
    pub tmp2: Vec<f64>,
    pub w: Vec<f64>,
    pub ca: Vec<f64>,
}

/// The weights `1 / var` into `w` and their sum.
#[inline]
fn block_weights(var: &[f64], w: &mut Vec<f64>) -> f64 {
    w.clear();
    w.extend(var.iter().map(|v| 1.0 / v));
    np_sum(w)
}

/// `cA[c] = Σ_r Asd[r][c] w[r]` into `ca` (NumPy's axis-0 reduction: sequential over the rows).
#[inline]
fn block_ca(asd: &Rows, w: &[f64], u: usize, ca: &mut Vec<f64>) {
    ca.clear();
    ca.resize(u, 0.0);
    for (row, &wr) in asd.a.chunks(u).zip(w) {
        for (c, &a) in row.iter().enumerate() {
            ca[c] += a * wr;
        }
    }
}

/// The group's vᵀPv. `Original`: `yᵀWy − wy² / sw`, floored at 0 (the original expression: at a large SD
/// receiver clock the two terms cancel to their last bits). `Stable`: the
/// same quantity as `Σ w (y − ȳ_w)²`.
#[inline]
fn block_vtpv(ysd: &[f64], w: &[f64], wy: f64, sw: f64, form: VtpvForm) -> f64 {
    match form {
        VtpvForm::Original => {
            let ywy: f64 = (0..ysd.len()).map(|r| ysd[r] * (w[r] * ysd[r])).sum();
            (ywy - wy * wy / sw).max(0.0)
        }
        VtpvForm::Stable => {
            let ybar = wy / sw;
            (0..ysd.len())
                .map(|r| {
                    let d = ysd[r] - ybar;
                    w[r] * (d * d)
                })
                .sum()
        }
    }
}

/// The block of one group, in a fixed operation order.
#[must_use]
pub fn block(
    asd: &Rows,
    g: Vec<usize>,
    var: &[f64],
    ysd: &[f64],
    key: Option<char>,
    form: VtpvForm,
) -> Block {
    let (mut w, mut ca) = (Vec::new(), Vec::new());
    block_with(asd, g, var, ysd, key, form, &mut w, &mut ca)
}

/// [`block`] on the design held in `sc` (`asd`, `var`, `y`), with its temporaries.
#[must_use]
pub fn block_in(sc: &mut Scratch, g: Vec<usize>, key: Option<char>, form: VtpvForm) -> Block {
    block_with(&sc.asd, g, &sc.var, &sc.y, key, form, &mut sc.w, &mut sc.ca)
}

#[allow(clippy::too_many_arguments)]
fn block_with(
    asd: &Rows,
    g: Vec<usize>,
    var: &[f64],
    ysd: &[f64],
    key: Option<char>,
    form: VtpvForm,
    w: &mut Vec<f64>,
    ca: &mut Vec<f64>,
) -> Block {
    let nr = ysd.len();
    let u = g.len();
    let sw = block_weights(var, w);
    block_ca(asd, w, u, ca);
    let wy: f64 = (0..nr).map(|r| w[r] * ysd[r]).sum();
    // m[c1][c2] = Σ_r Asd[r][c1] (Asd[r][c2] w[r]) − cA[c1] cA[c2] / sw and b[c1] = Σ_r Asd[r][c1] w[r] y[r]
    // − cA[c1] (wy / sw): accumulated row by row (every element's sum runs over the rows in order, from 0)
    let mut m = vec![0.0; u * u];
    let mut b = vec![0.0; u];
    for r in 0..nr {
        let (row, wr, yr) = (asd.row(r), w[r], ysd[r]);
        for c1 in 0..u {
            let a1 = row[c1];
            let mrow = &mut m[c1 * u..(c1 + 1) * u];
            for c2 in 0..u {
                mrow[c2] += a1 * (row[c2] * wr);
            }
            b[c1] += a1 * wr * yr;
        }
    }
    for c1 in 0..u {
        for c2 in 0..u {
            m[c1 * u + c2] -= ca[c1] * ca[c2] / sw;
        }
        b[c1] -= ca[c1] * (wy / sw);
    }
    let ss = if key.is_some() {
        ysd.iter().map(|y| y * y).sum()
    } else {
        0.0
    };
    Block {
        g,
        m,
        b,
        vtpv: block_vtpv(ysd, w, wy, sw, form),
        nrow: nr - 1,
        ss,
        key,
    }
}

/// Adds a block into the normal equations (sequential, in group order: a fixed sum order).
pub fn scatter(n: &mut DMat, bvec: &mut [f64], acc: &mut Acc, blk: &Block) {
    let u = blk.g.len();
    for c1 in 0..u {
        for c2 in 0..u {
            n.add(blk.g[c1], blk.g[c2], blk.m[c1 * u + c2]);
        }
        bvec[blk.g[c1]] += blk.b[c1];
    }
    acc.vtpv += blk.vtpv;
    acc.nrow += blk.nrow;
    match blk.key {
        Some('p') => {
            acc.ss_ph += blk.ss;
            acc.n_ph += blk.nrow;
        }
        Some(_) => {
            acc.ss_cd += blk.ss;
            acc.n_cd += blk.nrow;
        }
        None => {}
    }
}

/// The part of a group's block that changes between the builds of one design: the rows and columns of the
/// first `ROVER` parameters (the rover's coordinates, whose design entries `−e` follow the rover), the
/// right-hand side and vᵀPv. The rest of the block (`AᵀWA − …` on the troposphere, ambiguity, ionosphere
/// and wide-lane columns) is the same number in every build and is kept (`IonoWeighted`'s cache).
#[derive(Debug, Clone)]
pub struct RoverBlock {
    /// Parameter indices.
    pub g: Vec<usize>,
    /// `m[c1][c2]` for `c1 < ROVER`, row-major `ROVER × u`.
    pub rows: Vec<f64>,
    /// `m[c1][c2]` for `c1 ≥ ROVER`, `c2 < ROVER`, row-major `(u − ROVER) × ROVER`.
    pub cols: Vec<f64>,
    /// Right-hand side.
    pub b: Vec<f64>,
    /// `vᵀPv` of the group.
    pub vtpv: f64,
    /// DDs in the group (`n − 1`).
    pub nrow: usize,
}

/// The leading parameters a `RoverBlock` recomputes: the rover's X, Y, Z.
pub const ROVER: usize = 3;

/// The rover part of a group's block: the same numbers `block` gives for those elements.
#[must_use]
pub fn rover_block(
    asd: &Rows,
    g: Vec<usize>,
    var: &[f64],
    ysd: &[f64],
    form: VtpvForm,
) -> RoverBlock {
    let (mut w, mut ca) = (Vec::new(), Vec::new());
    rover_block_with(asd, g, var, ysd, form, &mut w, &mut ca)
}

/// [`rover_block`] on the design held in `sc`.
#[must_use]
pub fn rover_block_in(sc: &mut Scratch, g: Vec<usize>, form: VtpvForm) -> RoverBlock {
    rover_block_with(&sc.asd, g, &sc.var, &sc.y, form, &mut sc.w, &mut sc.ca)
}

fn rover_block_with(
    asd: &Rows,
    g: Vec<usize>,
    var: &[f64],
    ysd: &[f64],
    form: VtpvForm,
    w: &mut Vec<f64>,
    ca: &mut Vec<f64>,
) -> RoverBlock {
    let nr = ysd.len();
    let u = g.len();
    let sw = block_weights(var, w);
    block_ca(asd, w, u, ca);
    let wy: f64 = (0..nr).map(|r| w[r] * ysd[r]).sum();
    let mut rows = vec![0.0; ROVER * u];
    let mut cols = vec![0.0; (u - ROVER) * ROVER];
    let mut b = vec![0.0; u];
    for r in 0..nr {
        let (row, wr, yr) = (asd.row(r), w[r], ysd[r]);
        for c1 in 0..ROVER {
            let a1 = row[c1];
            let mrow = &mut rows[c1 * u..(c1 + 1) * u];
            for c2 in 0..u {
                mrow[c2] += a1 * (row[c2] * wr);
            }
        }
        for c1 in ROVER..u {
            let a1 = row[c1];
            let mrow = &mut cols[(c1 - ROVER) * ROVER..(c1 - ROVER + 1) * ROVER];
            for c2 in 0..ROVER {
                mrow[c2] += a1 * (row[c2] * wr);
            }
        }
        for c1 in 0..u {
            b[c1] += row[c1] * wr * yr;
        }
    }
    for c1 in 0..ROVER {
        for c2 in 0..u {
            rows[c1 * u + c2] -= ca[c1] * ca[c2] / sw;
        }
    }
    for c1 in ROVER..u {
        for c2 in 0..ROVER {
            cols[(c1 - ROVER) * ROVER + c2] -= ca[c1] * ca[c2] / sw;
        }
    }
    for c1 in 0..u {
        b[c1] -= ca[c1] * (wy / sw);
    }
    RoverBlock {
        g,
        rows,
        cols,
        b,
        vtpv: block_vtpv(ysd, w, wy, sw, form),
        nrow: nr - 1,
    }
}

/// Adds a rover block into normal equations whose other elements are already there.
pub fn scatter_rover(n: &mut DMat, bvec: &mut [f64], acc: &mut Acc, blk: &RoverBlock) {
    let u = blk.g.len();
    for c1 in 0..ROVER {
        for c2 in 0..u {
            n.add(blk.g[c1], blk.g[c2], blk.rows[c1 * u + c2]);
        }
    }
    for c1 in ROVER..u {
        for c2 in 0..ROVER {
            n.add(blk.g[c1], blk.g[c2], blk.cols[(c1 - ROVER) * ROVER + c2]);
        }
    }
    for c1 in 0..u {
        bvec[blk.g[c1]] += blk.b[c1];
    }
    acc.vtpv += blk.vtpv;
    acc.nrow += blk.nrow;
}

/// A group's block seen as rows of the normal equations, for the parallel scatter: row `c1` belongs to
/// parameter `g()[c1]`.
pub trait RowBlock: Send + Sync {
    /// Parameter indices.
    fn g(&self) -> &[usize];
    /// Adds row `c1` into `row` (the whole matrix row of parameter `g()[c1]`) and into its right-hand
    /// side element.
    fn add_row(&self, c1: usize, row: &mut [f64], b: &mut f64);
    /// Adds the block's accumulators.
    fn add_acc(&self, acc: &mut Acc);
}

impl RowBlock for Block {
    fn g(&self) -> &[usize] {
        &self.g
    }

    fn add_row(&self, c1: usize, row: &mut [f64], b: &mut f64) {
        let u = self.g.len();
        for c2 in 0..u {
            row[self.g[c2]] += self.m[c1 * u + c2];
        }
        *b += self.b[c1];
    }

    fn add_acc(&self, acc: &mut Acc) {
        acc.vtpv += self.vtpv;
        acc.nrow += self.nrow;
        match self.key {
            Some('p') => {
                acc.ss_ph += self.ss;
                acc.n_ph += self.nrow;
            }
            Some(_) => {
                acc.ss_cd += self.ss;
                acc.n_cd += self.nrow;
            }
            None => {}
        }
    }
}

impl RowBlock for RoverBlock {
    fn g(&self) -> &[usize] {
        &self.g
    }

    fn add_row(&self, c1: usize, row: &mut [f64], b: &mut f64) {
        let u = self.g.len();
        if c1 < ROVER {
            for c2 in 0..u {
                row[self.g[c2]] += self.rows[c1 * u + c2];
            }
        } else {
            for c2 in 0..ROVER {
                row[self.g[c2]] += self.cols[(c1 - ROVER) * ROVER + c2];
            }
        }
        *b += self.b[c1];
    }

    fn add_acc(&self, acc: &mut Acc) {
        acc.vtpv += self.vtpv;
        acc.nrow += self.nrow;
    }
}

/// Groups per rayon task of the block computation. The blocks are pure functions of their groups and are
/// scattered in group order afterwards, so neither the chunking nor the thread count touches a number.
pub const GROUP_CHUNK: usize = 64;

/// Groups whose blocks are alive at once: a round's blocks are computed, then scattered, so a build's
/// worth of them (~1 GB on a long multi-system pair) never exists at once.
pub const GROUP_ROUND: usize = 16 * GROUP_CHUNK;

/// Whether the assembly forks onto the rayon pool. Not when the engine itself runs as a task of the
/// caller's pool (`rayon::current_thread_index()` is `Some`: the app runs several solves at once on one
/// pool, and forks inside a saturated pool only make every build wait on the other solves' threads — the
/// a six-session local job ran 2–3× slower that way). `GEOM_ENGINE_PAR=par` forks anyway, `GEOM_ENGINE_PAR=seq` never
/// does (measurements). The numbers do not depend on the choice.
fn fork_here() -> bool {
    static CHOICE: std::sync::OnceLock<Option<bool>> = std::sync::OnceLock::new();
    let choice = CHOICE.get_or_init(|| match std::env::var("GEOM_ENGINE_PAR").as_deref() {
        Ok("seq") => Some(false),
        Ok("par") => Some(true),
        _ => None,
    });
    choice.unwrap_or_else(|| rayon::current_thread_index().is_none())
}

/// Assembles the normal equations with the numbers of a sequential assembly in group order. Per round of
/// `GROUP_ROUND` groups the blocks are computed in parallel (rayon, `GROUP_CHUNK` groups per task), then
/// scattered by the calling thread in group order, so every element of `N`, `b` and the accumulators
/// receives its additions in the order a sequential loop over the groups would give it.
///
/// The scatter is sequential on purpose: it is a small part of a build (adds into `N`), and the engine
/// runs inside a shared rayon pool when a job processes several baselines at once — one coarse parallel
/// section per round, with the calling thread doing the rest, costs nothing when the pool is busy, where
/// a fork-join per matrix row made each build wait on whichever thread was descheduled.
pub fn assemble_rows<T, F>(
    groups: &[(usize, Vec<usize>)],
    n: &mut DMat,
    bvec: &mut [f64],
    acc: &mut Acc,
    f: F,
) where
    T: RowBlock,
    F: Fn(usize, &[usize], &mut Scratch, &mut Vec<T>) + Sync,
{
    use rayon::prelude::*;
    let p = n.n;
    if !fork_here() {
        // one thread: block by block, scattered at once (no round's worth of blocks in memory)
        let mut sc = Scratch::default();
        let mut out: Vec<T> = Vec::new();
        for (e, cj) in groups {
            out.clear();
            f(*e, cj, &mut sc, &mut out);
            for blk in &out {
                for (c1, &gi) in blk.g().iter().enumerate() {
                    blk.add_row(c1, &mut n.data[gi * p..(gi + 1) * p], &mut bvec[gi]);
                }
                blk.add_acc(acc);
            }
        }
        return;
    }
    for round in groups.chunks(GROUP_ROUND) {
        let chunks: Vec<Vec<T>> = round
            .par_chunks(GROUP_CHUNK)
            .map(|chunk| {
                let mut out = Vec::with_capacity(chunk.len() * 2);
                let mut sc = Scratch::default();
                for (e, cj) in chunk {
                    f(*e, cj, &mut sc, &mut out);
                }
                out
            })
            .collect();
        for blk in chunks.iter().flatten() {
            for (c1, &gi) in blk.g().iter().enumerate() {
                blk.add_row(c1, &mut n.data[gi * p..(gi + 1) * p], &mut bvec[gi]);
            }
            blk.add_acc(acc);
        }
    }
}

/// `assemble_rows` for the full blocks.
pub fn assemble_blocks<F>(
    groups: &[(usize, Vec<usize>)],
    n: &mut DMat,
    bvec: &mut [f64],
    acc: &mut Acc,
    f: F,
) where
    F: Fn(usize, &[usize], &mut Scratch, &mut Vec<Block>) + Sync,
{
    assemble_rows(groups, n, bvec, acc, f);
}

/// The float estimate of one run.
#[derive(Debug, Clone)]
#[allow(missing_docs)]
pub struct Est {
    pub x: Vec<f64>,
    pub n: DMat,
    pub s0: f64,
    pub res_ph: F2,
    pub res_cd: F2,
    pub acc: Acc,
    pub wf_ph: F2,
    pub wf_cd: F2,
    pub cf_ph: F2,
    pub use_: B2,
    pub code_ok: B2,
    pub amb0: Vec<f64>,
    pub p: usize,
}

/// Everything the solve's closures share. The rover-dependent corrections (`cb`) and `at_b` are what
/// `geom` re-evaluates when the rover moves.
#[allow(missing_docs)]
pub struct Ctx<'a> {
    pub opts: Options,
    pub t: Vec<f64>,
    pub ne: usize,
    pub ns: usize,
    pub cols: Vec<Sat>,
    pub sysc: Vec<char>,
    pub pairs: BTreeMap<Sat, (u8, u8)>,
    pub fcn: BTreeMap<Sat, i32>,
    pub base_xyz: [f64; 3],
    pub x_r: [f64; 3],
    pub a: Station,
    pub b: Station,
    pub ca: Corrections,
    pub cb: Corrections,
    pub at_b: [f64; 3],
    pub cal_b: Option<AntennaCal>,
    pub tropo: &'a TropoInputs<'a>,
    pub lam1: Vec<f64>,
    pub lam2: Vec<f64>,
    pub aif: Vec<f64>,
    pub bif: Vec<f64>,
    pub valid: B2,
    pub var_ph: F2,
    pub var_cd: Option<F2>,
    pub obs_ph: F2,
    pub obs_cd: F2,
    pub sdgf: F2,
    pub corr_gf: F2,
    pub sdmw: F2,
    pub corr_mw: F2,
    pub half_sd: B2,
    pub nn: usize,
    pub i_tb: usize,
    pub i_ta: usize,
    pub n_ta: usize,
    pub i_amb: usize,
    pub kk: Vec<usize>,
    pub ww: Vec<f64>,
    pub ntr: usize,
    pub step: f64,
    pub w_rel: f64,
    pub groups_all: Vec<(usize, Vec<usize>)>,
    pub n_sys: usize,
    pub rej_cd0: B2,
}

impl Ctx<'_> {
    /// The wet delay at each epoch from the nodes.
    pub fn wet(&self, xx: &[f64], which: char) -> Vec<f64> {
        if self.nn == 0 || (which == 'A' && self.n_ta == 0) {
            return vec![0.0; self.ne];
        }
        let i0 = if which == 'B' { self.i_tb } else { self.i_ta };
        (0..self.ne)
            .map(|e| {
                xx[i0 + self.kk[e]] * (1.0 - self.ww[e]) + xx[i0 + self.kk[e] + 1] * self.ww[e]
            })
            .collect()
    }

    /// The modelled SD range (m).
    pub fn model_sd(&self, xx: &[f64], rho_b: &F2) -> F2 {
        let wb = self.wet(xx, 'B');
        let wa = self.wet(xx, 'A');
        let mut m = F2::zeros(self.ne, self.ns);
        for e in 0..self.ne {
            for j in 0..self.ns {
                let t1 = rho_b.at(e, j) - self.ca.rho.at(e, j);
                let t2 = CLIGHT * (self.b.sat_clk.at(e, j) - self.a.sat_clk.at(e, j));
                let t3 = t1 - t2;
                let t4 = self.cb.dry.at(e, j) - self.ca.dry.at(e, j);
                let t5 = t3 + t4;
                let t6 = self.cb.ant.at(e, j) - self.ca.ant.at(e, j);
                let mm = t5 + t6;
                let mm = mm + self.cb.map.at(e, j) * wb[e] - self.ca.map.at(e, j) * wa[e];
                m.set(e, j, mm);
            }
        }
        m
    }

    /// Range and LOS from the rover at `xx[:3]`, re-evaluating the rover's troposphere and antenna terms
    /// when it moved over 1 cm.
    pub fn geom(&mut self, xx: &[f64]) -> (F2, Vec3Arr) {
        let x = [xx[0], xx[1], xx[2]];
        let pos: Vec<[f64; 3]> = self.b.tide.iter().map(|tide| add3(x, *tide)).collect();
        if norm3(sub3(x, self.at_b)) > 0.01 {
            let c = corrections(
                &self.b,
                &pos,
                &self.cols,
                &self.pairs,
                &self.fcn,
                self.cal_b.as_ref(),
                &self.t,
                &self.opts,
                self.tropo,
            );
            self.cb.dry = c.dry;
            self.cb.map = c.map;
            self.cb.ant = c.ant;
            self.cb.antgf = c.antgf;
            self.at_b = x;
        }
        let mut rho = F2::zeros(self.ne, self.ns);
        let mut e = Vec3Arr::zeros(self.ne, self.ns);
        for ep in 0..self.ne {
            for j in 0..self.ns {
                let (r, u) = geometry(self.b.sat_xyz.at(ep, j), pos[ep]);
                rho.set(ep, j, r);
                e.set(ep, j, u);
            }
        }
        (rho, e)
    }

    /// Each epoch / system's SD residuals less their weighted mean over `use`.
    pub fn clock_free(&self, res: &F2, var_eff: &F2, use_: &B2, also: Option<&B2>) -> F2 {
        let mut out = F2::nans(self.ne, self.ns);
        for (e, cj0) in &self.groups_all {
            let cj: Vec<usize> = cj0.iter().copied().filter(|&j| use_.at(*e, j)).collect();
            if cj.len() >= 2 {
                let w: Vec<f64> = cj.iter().map(|&j| 1.0 / var_eff.at(*e, j)).collect();
                let wr: Vec<f64> = cj
                    .iter()
                    .zip(&w)
                    .map(|(&j, wi)| wi * res.at(*e, j))
                    .collect();
                let m = np_sum(&wr) / np_sum(&w);
                for &j in &cj {
                    out.set(*e, j, res.at(*e, j) - m);
                }
                if let Some(al) = also {
                    for &j in cj0 {
                        if al.at(*e, j) {
                            out.set(*e, j, res.at(*e, j) - m);
                        }
                    }
                }
            }
        }
        out
    }

    /// A tie between two parameters.
    pub fn tie(n: &mut DMat, b: &mut [f64], xx: &[f64], a_: usize, b_: usize, wr: f64) {
        let d = xx[b_] - xx[a_];
        n.add(a_, a_, wr);
        n.add(b_, b_, wr);
        n.add(a_, b_, -wr);
        n.add(b_, a_, -wr);
        b[a_] += wr * d;
        b[b_] -= wr * d;
    }

    /// The troposphere's prior, random walk and the two stations' tie; the rover held.
    pub fn trop_priors(&self, n: &mut DMat, b: &mut [f64], xx: &[f64]) {
        let o = &self.opts;
        if o.rover_sigma > 0.0 {
            let wp = 1.0 / (o.rover_sigma * o.rover_sigma);
            for i in 0..3 {
                n.add(i, i, wp);
                b[i] += wp * (self.x_r[i] - xx[i]);
            }
        }
        if self.nn > 0 {
            for (i0, nt) in [(self.i_tb, self.nn), (self.i_ta, self.n_ta)] {
                if nt == 0 {
                    continue;
                }
                let w0 = 1.0 / (o.trop_sigma * o.trop_sigma);
                let p0 = if o.trop == "station" {
                    o.trop_init
                } else {
                    0.0
                };
                n.add(i0, i0, w0);
                b[i0] += w0 * (p0 - xx[i0]);
                for k in 0..nt - 1 {
                    Self::tie(
                        n,
                        b,
                        xx,
                        i0 + k,
                        i0 + k + 1,
                        1.0 / (o.trop_rw * o.trop_rw * self.step),
                    );
                }
            }
            if self.n_ta > 0 && self.w_rel > 0.0 {
                for k in 0..self.nn {
                    Self::tie(n, b, xx, self.i_ta + k, self.i_tb + k, self.w_rel);
                }
            }
        }
    }

    /// One epoch / system / signal group's DDs, as the SDs with their common clock eliminated.
    pub fn normal_add(
        n: &mut DMat,
        b: &mut [f64],
        asd: &Rows,
        g: &[usize],
        var: &[f64],
        ysd: &[f64],
        acc: &mut Acc,
        key: Option<char>,
        form: VtpvForm,
    ) {
        scatter(n, b, acc, &block(asd, g.to_vec(), var, ysd, key, form));
    }

    /// IGG-III weights.
    pub fn igg(&self, z: f64) -> f64 {
        let (k0, k1) = (self.opts.igg_k0, self.opts.igg_k1);
        if z > k0 && z <= k1 {
            (k0 / z) * ((k1 - z) / (k1 - k0)).powi(2)
        } else if z > k1 {
            0.0
        } else {
            1.0
        }
    }

    /// `|cf| / sqrt(var)` over the robust scale (never under 0.2).
    pub fn robust_z(&self, cf: &F2, var: &F2, wf: &F2) -> F2 {
        let mut zz: Vec<f64> = Vec::new();
        for i in 0..cf.data.len() {
            if wf.data[i] > 0.0 && cf.data[i].is_finite() {
                zz.push(cf.data[i].abs() / var.data[i].sqrt());
            }
        }
        let s = if zz.is_empty() {
            1.0
        } else {
            1.4826 * median(&zz)
        };
        let scale = s.max(0.2);
        cf.zip_map(var, |c, v| c.abs() / v.sqrt() / scale)
    }

    /// Hard rejection: the worst satellite per epoch / system over `outlier_k`.
    pub fn screen(&self, cf: &F2, var: &F2, wf: &mut F2, scale: f64) -> usize {
        let mut worst: Vec<(f64, usize, usize)> = Vec::new();
        for (e, cj_all) in &self.groups_all {
            let cj: Vec<usize> = cj_all
                .iter()
                .copied()
                .filter(|&j| wf.at(*e, j) > 0.0)
                .collect();
            if cj.len() < 3 {
                continue;
            }
            let z: Vec<f64> = cj
                .iter()
                .map(|&j| cf.at(*e, j).abs() / var.at(*e, j).sqrt() / scale)
                .collect();
            let k = argmax(&z);
            if z[k] > self.opts.outlier_k {
                worst.push((z[k], *e, cj[k]));
            }
        }
        if worst.is_empty() {
            return 0;
        }
        let zmax = worst.iter().map(|w| w.0).fold(f64::NEG_INFINITY, f64::max);
        let cut: Vec<&(f64, usize, usize)> = worst
            .iter()
            .filter(|w| w.0 > self.opts.outlier_k.max(zmax / 2.0))
            .collect();
        for w in &cut {
            wf.set(w.1, w.2, 0.0);
        }
        cut.len()
    }

    /// Batch LS with these arcs from rover position `pos0`: robust rounds.
    pub fn estimate(
        &mut self,
        arc: &I2,
        n_amb: usize,
        links: &[(usize, usize, f64)],
        pos0: [f64; 3],
        warm: Option<&Est>,
    ) -> Result<Est, EngineError> {
        let o = self.opts.clone();
        let use_ = self.valid.zip_map(arc, |v, a| v && a >= 0);
        let p = self.i_amb + n_amb;
        let mut x = vec![0.0; p];
        x[..3].copy_from_slice(&pos0);
        if self.nn > 0 {
            let v = if o.trop == "station" {
                o.trop_init
            } else {
                0.0
            };
            for k in 0..self.nn {
                x[self.i_tb + k] = v;
            }
            for k in 0..self.n_ta {
                x[self.i_ta + k] = o.trop_init;
            }
        }
        let (rho_b, _eb) = self.geom(&x);
        // starting ambiguities: arc mean of the SD residual with each epoch's median removed
        let m0 = self.model_sd(&x, &rho_b);
        let mut r0 = F2::nans(self.ne, self.ns);
        for e in 0..self.ne {
            for j in 0..self.ns {
                if use_.at(e, j) {
                    r0.set(e, j, self.obs_ph.at(e, j) - m0.at(e, j));
                }
            }
        }
        let systems: BTreeSet<char> = self.sysc.iter().copied().collect();
        for s in systems {
            let cj: Vec<usize> = (0..self.ns).filter(|&j| self.sysc[j] == s).collect();
            for e in 0..self.ne {
                let sub: Vec<f64> = cj
                    .iter()
                    .map(|&j| if use_.at(e, j) { r0.at(e, j) } else { f64::NAN })
                    .collect();
                let med = if sub.iter().any(|v| v.is_finite()) {
                    nanmedian(&sub)
                } else {
                    0.0
                };
                for &j in &cj {
                    let v = r0.at(e, j) - med;
                    r0.set(e, j, v);
                }
            }
        }
        let mut amb0 = vec![0.0; n_amb];
        if n_amb > 0 {
            let mut sums = vec![0.0; n_amb];
            let mut cnt = vec![0usize; n_amb];
            for e in 0..self.ne {
                for j in 0..self.ns {
                    if use_.at(e, j) {
                        let a = arc.at(e, j) as usize;
                        sums[a] += r0.at(e, j);
                        cnt[a] += 1;
                    }
                }
            }
            for a in 0..n_amb {
                amb0[a] = sums[a] / (cnt[a].max(1) as f64);
            }
        }
        x[self.i_amb..].copy_from_slice(&amb0);
        let mut wf_ph = use_.map(|u| if u { 1.0 } else { 0.0 });
        if let Some(w) = warm {
            wf_ph = use_.zip_map(&w.wf_ph, |u, v| if u { v } else { 0.0 });
        }
        let code_ok = if self.var_cd.is_some() {
            let mut c = B2::falses(self.ne, self.ns);
            for i in 0..c.data.len() {
                c.data[i] =
                    use_.data[i] && self.obs_cd.data[i].is_finite() && !self.rej_cd0.data[i];
            }
            c
        } else {
            B2::falses(self.ne, self.ns)
        };
        let mut wf_cd = code_ok.map(|u| if u { 1.0 } else { 0.0 });
        if let Some(w) = warm {
            wf_cd = code_ok.zip_map(&w.wf_cd, |u, v| if u { v } else { 0.0 });
        }
        let mut s0 = 1.0;
        let mut last_pos: Option<[f64; 3]> = None;
        let (mut n_mat, mut res_ph, mut res_cd, mut acc);
        n_mat = DMat::zeros(p);
        res_ph = F2::nans(self.ne, self.ns);
        res_cd = F2::nans(self.ne, self.ns);
        acc = Acc::default();
        for rnd in 0..=o.max_reject_rounds {
            let its = if rnd == 0 { 6 } else { 2 };
            for it in 0..its {
                let (rb, eb) = self.geom(&x);
                let (n2, b2, _rp, _rc, _acc) = self.assemble(
                    &x, &rb, &eb, arc, n_amb, links, &amb0, &wf_ph, &wf_cd, &use_,
                );
                let dx = n2
                    .solve(&b2)
                    .ok_or_else(|| EngineError("singular normal equations".into()))?;
                for i in 0..p {
                    x[i] += dx[i];
                }
                let mx = dx[..3].iter().map(|v| v.abs()).fold(0.0, f64::max);
                if mx < 1e-5 && it > 0 {
                    break;
                }
            }
            let (rb, eb) = self.geom(&x);
            let (n2, _b2, rp, rc, a2) = self.assemble(
                &x, &rb, &eb, arc, n_amb, links, &amb0, &wf_ph, &wf_cd, &use_,
            );
            n_mat = n2;
            res_ph = rp;
            res_cd = rc;
            acc = a2;
            let dof = (acc.nrow as i64 - p as i64 + self.n_sys as i64).max(1);
            s0 = (acc.vtpv / dof as f64).sqrt();
            if rnd == o.max_reject_rounds {
                break;
            }
            if o.robust == "igg3"
                && let Some(lp) = last_pos
                && (0..3).map(|i| (x[i] - lp[i]).abs()).fold(0.0, f64::max) < 2e-4
            {
                break;
            }
            last_pos = Some([x[0], x[1], x[2]]);
            let ve_ph = self.var_ph.zip_map(&wf_ph, |v, w| v / w.max(1e-12));
            let cf_ph = self.clock_free(&res_ph, &ve_ph, &wf_ph.map(|w| w > 0.0), None);
            let cf_cd = self.var_cd.as_ref().map(|vc| {
                let ve = vc.zip_map(&wf_cd, |v, w| v / w.max(1e-12));
                self.clock_free(&res_cd, &ve, &wf_cd.map(|w| w > 0.0), None)
            });
            if o.robust == "igg3" {
                let mut changed = 0.0f64;
                let var_cd = self.var_cd.clone();
                let mut sets: Vec<(&F2, &F2, &mut F2, &B2)> =
                    vec![(&cf_ph, &self.var_ph, &mut wf_ph, &use_)];
                if let (Some(cfc), Some(vc)) = (&cf_cd, &var_cd) {
                    sets.push((cfc, vc, &mut wf_cd, &code_ok));
                }
                for (cf, var, wf, base_use) in sets {
                    let z = self.robust_z(cf, var, wf);
                    let mut new = wf.clone();
                    for i in 0..new.data.len() {
                        let m = base_use.data[i] && cf.data[i].is_finite();
                        let v = if m { self.igg(z.data[i]) } else { wf.data[i] };
                        let v = if base_use.data[i] { v } else { 0.0 };
                        changed = changed.max((v - wf.data[i]).abs());
                        new.data[i] = v;
                    }
                    *wf = new;
                }
                if changed < 0.05 {
                    break;
                }
            } else {
                let k_ph = self.screen(&cf_ph, &self.var_ph.clone(), &mut wf_ph, s0);
                let mut k_cd = 0;
                if let (Some(cfc), Some(vc)) = (&cf_cd, &self.var_cd.clone()) {
                    let mut zz = Vec::new();
                    for i in 0..cfc.data.len() {
                        if wf_cd.data[i] > 0.0 && cfc.data[i].is_finite() {
                            zz.push(cfc.data[i].abs() / vc.data[i].sqrt());
                        }
                    }
                    let scale = if zz.is_empty() {
                        1.0
                    } else {
                        1.4826 * median(&zz)
                    };
                    k_cd = self.screen(cfc, vc, &mut wf_cd, scale);
                }
                if k_ph == 0 && k_cd == 0 {
                    break;
                }
            }
        }
        let ve_ph = self.var_ph.zip_map(&wf_ph, |v, w| v / w.max(1e-12));
        let cf_ph = self.clock_free(&res_ph, &ve_ph, &wf_ph.map(|w| w > 0.0), Some(&use_));
        Ok(Est {
            x,
            n: n_mat,
            s0,
            res_ph,
            res_cd,
            acc,
            wf_ph,
            wf_cd,
            cf_ph,
            use_,
            code_ok,
            amb0,
            p,
        })
    }

    /// The normal equations at `xx`.
    #[allow(clippy::type_complexity)]
    pub fn assemble(
        &self,
        xx: &[f64],
        rho_b: &F2,
        eb: &Vec3Arr,
        arc: &I2,
        n_amb: usize,
        links: &[(usize, usize, f64)],
        amb0: &[f64],
        wf_ph: &F2,
        wf_cd: &F2,
        use_: &B2,
    ) -> (DMat, Vec<f64>, F2, F2, Acc) {
        let p = self.i_amb + n_amb;
        let mut n = DMat::zeros(p);
        let mut bvec = vec![0.0; p];
        let msd = self.model_sd(xx, rho_b);
        let mut res_ph = F2::nans(self.ne, self.ns);
        let mut res_cd = F2::nans(self.ne, self.ns);
        for e in 0..self.ne {
            for j in 0..self.ns {
                let a = arc.at(e, j);
                let amb = if a >= 0 {
                    xx[self.i_amb + a as usize]
                } else {
                    0.0
                };
                res_ph.set(e, j, self.obs_ph.at(e, j) - msd.at(e, j) - amb);
                res_cd.set(e, j, self.obs_cd.at(e, j) - msd.at(e, j));
            }
        }
        let mut acc = Acc::default();
        let form = self.opts.compat.vtpv;
        let group = |e: usize, cj_all: &[usize], sc: &mut Scratch, out: &mut Vec<Block>| {
            let cj0: Vec<usize> = cj_all.iter().copied().filter(|&j| use_.at(e, j)).collect();
            let n0 = cj0.len();
            if n0 < 2 {
                return;
            }
            let u = 3 + self.ntr + n0;
            sc.asd_full.reset(n0, u);
            let asd = &mut sc.asd_full;
            let mut g: Vec<usize> = vec![0, 1, 2];
            let mut c0 = 3;
            for (r, &j) in cj0.iter().enumerate() {
                let ee = eb.at(e, j);
                asd[(r, 0)] = -ee[0];
                asd[(r, 1)] = -ee[1];
                asd[(r, 2)] = -ee[2];
            }
            if self.nn > 0 {
                let k = self.kk[e];
                for (r, &j) in cj0.iter().enumerate() {
                    asd[(r, c0)] = self.cb.map.at(e, j) * (1.0 - self.ww[e]);
                    asd[(r, c0 + 1)] = self.cb.map.at(e, j) * self.ww[e];
                }
                g.push(self.i_tb + k);
                g.push(self.i_tb + k + 1);
                c0 += 2;
                if self.n_ta > 0 {
                    for (r, &j) in cj0.iter().enumerate() {
                        asd[(r, c0)] = -self.ca.map.at(e, j) * (1.0 - self.ww[e]);
                        asd[(r, c0 + 1)] = -self.ca.map.at(e, j) * self.ww[e];
                    }
                    g.push(self.i_ta + k);
                    g.push(self.i_ta + k + 1);
                    c0 += 2;
                }
            }
            for r in 0..n0 {
                asd[(r, c0 + r)] = 1.0;
            }
            for &j in &cj0 {
                g.push(self.i_amb + arc.at(e, j) as usize);
            }
            // the rows with a weight, the reference row (highest elevation) subtracted from y
            // (`sel` / `var` / `y` / `asd` of the scratch are the block's design)
            // (the phase design keeps the common columns and each kept row's own ambiguity column —
            // `keep_cols = 0..c0 ++ c0 + sel[k]`, an identity block; the code design the common ones)
            let sub = |sc: &mut Scratch,
                       keep: &dyn Fn(usize) -> bool,
                       amb: bool,
                       var_src: &F2,
                       wf: &F2,
                       res: &F2|
             -> bool {
                sc.sel.clear();
                sc.sel.extend((0..n0).filter(|&r| keep(r)));
                let ns = sc.sel.len();
                if ns < 2 {
                    return false;
                }
                let asd_full = &sc.asd_full;
                sc.asd.reset(ns, c0 + if amb { ns } else { 0 });
                for (r, &sr) in sc.sel.iter().enumerate() {
                    for c in 0..c0 {
                        sc.asd[(r, c)] = asd_full[(sr, c)];
                    }
                    if amb {
                        sc.asd[(r, c0 + r)] = asd_full[(sr, c0 + sr)];
                    }
                }
                sc.var.clear();
                sc.var.extend(
                    sc.sel
                        .iter()
                        .map(|&r| var_src.at(e, cj0[r]) / wf.at(e, cj0[r])),
                );
                sc.tmp.clear();
                sc.tmp.extend(sc.sel.iter().map(|&r| res.at(e, cj0[r])));
                sc.tmp2.clear();
                sc.tmp2
                    .extend(sc.sel.iter().map(|&r| self.cb.el.at(e, cj0[r])));
                let rp = argmax(&sc.tmp2);
                let y_ref = sc.tmp[rp];
                sc.y.clear();
                sc.y.extend(sc.tmp.iter().map(|v| v - y_ref));
                true
            };
            if sub(
                sc,
                &|r| wf_ph.at(e, cj0[r]) > 0.0,
                true,
                &self.var_ph,
                wf_ph,
                &res_ph,
            ) {
                let mut g_sub: Vec<usize> = g[..c0].to_vec();
                g_sub.extend(sc.sel.iter().map(|&r| g[c0 + r]));
                out.push(block_in(sc, g_sub, Some('p'), form));
            }
            if let Some(vc) = &self.var_cd
                && sub(
                    sc,
                    &|r| wf_cd.at(e, cj0[r]) > 0.0,
                    false,
                    vc,
                    wf_cd,
                    &res_cd,
                )
            {
                let g_sub: Vec<usize> = g[..c0].to_vec();
                out.push(block_in(sc, g_sub, Some('c'), form));
            }
        };
        assemble_blocks(&self.groups_all, &mut n, &mut bvec, &mut acc, group);
        if n_amb > 0 {
            let w = 1.0 / (self.opts.amb_prior * self.opts.amb_prior);
            for a in 0..n_amb {
                let i = self.i_amb + a;
                n.add(i, i, w);
                bvec[i] += w * (amb0[a] - xx[i]);
            }
        }
        self.trop_priors(&mut n, &mut bvec, xx);
        for &(a_, b_, dt_) in links {
            Self::tie(
                &mut n,
                &mut bvec,
                xx,
                self.i_amb + a_,
                self.i_amb + b_,
                1.0 / (self.opts.amb_rw * self.opts.amb_rw * dt_),
            );
        }
        (n, bvec, res_ph, res_cd, acc)
    }
}

/// `StaticResult.cov`: `realistic_cov` with the data quality on top.
pub fn quality_cov(
    cov_formal: &[[f64; 3]; 3],
    re: &[[f64; 3]; 3],
    length_m: f64,
    hours: f64,
    fix_info: &mut Map<String, Value>,
    n_amb: usize,
    n_sat: usize,
    opts: &Options,
    notes: &mut Vec<String>,
    dh_m: f64,
    n_rej: usize,
    n_obs: usize,
    s0: f64,
    fdist: f64,
) -> [[f64; 3]; 3] {
    let mut is_fixed = fix_info.get("status").and_then(Value::as_str) == Some("fixed");
    let wl_n: [usize; 2] = fix_info
        .get("wl")
        .and_then(Value::as_array)
        .map(|a| {
            [
                a[0].as_u64().unwrap_or(0) as usize,
                a[1].as_u64().unwrap_or(0) as usize,
            ]
        })
        .unwrap_or([0, 0]);
    let nl0 = fix_info
        .get("nl")
        .and_then(Value::as_array)
        .and_then(|a| a.first())
        .and_then(Value::as_u64)
        .unwrap_or(0) as usize;
    let unfixed = opts.fix && !is_fixed;
    let mut partial = false;
    if is_fixed && hours < 6.0 && wl_n[1] > 0 && (nl0 as f64) < opts.fix_share_min * wl_n[1] as f64
    {
        fix_info.insert("partial".into(), jround(nl0 as f64 / wl_n[1] as f64, 2));
        notes.push(format!(
            "fixed {nl0} DDs of {} wide-lane candidates: the float σ",
            wl_n[1]
        ));
        is_fixed = false;
        partial = true;
    }
    let p = opts.sigma_params(is_fixed);
    let frag = if n_sat > 0 {
        n_amb as f64 / n_sat as f64 / (1.0 + 0.1 * hours)
    } else {
        1.0
    };
    let rej = n_rej as f64 / ((n_obs + n_rej).max(1)) as f64;
    let n_dropped = fix_info
        .get("dropped")
        .and_then(Value::as_array)
        .map_or(0, Vec::len);
    let drop = n_dropped as f64 / n_sat.max(1) as f64;
    fix_info.insert("frag".into(), jround(frag, 2));
    fix_info.insert("rej".into(), jround(rej, 4));
    let ratio = fix_info.get("ratio").and_then(Value::as_f64);
    let q = sigma_quality(
        &p, is_fixed, frag, rej, s0, ratio, fdist, unfixed, drop, partial,
    );
    fix_info.insert("sigma_q".into(), jround(q, 3));
    if q >= 1.5 {
        notes.push(format!(
            "data quality (arcs {frag:.1} a satellite, {:.1} % of the phase rejected, σ0 {s0:.2}, {n_dropped} satellites dropped from the fix): σ × {q:.1}",
            100.0 * rej
        ));
    }
    realistic_cov(cov_formal, re, length_m, hours, is_fixed, opts, dh_m, q)
}

/// The signals of one station as the dump writes them: `{sat: {band: [phase, code]}}`.
fn signals_json(sig: &BTreeMap<Sat, BandSignals>) -> Value {
    let mut m = Map::new();
    for (s, bands) in sig {
        let mut bm = Map::new();
        for (b, (ph, cd)) in bands {
            bm.insert(
                b.to_string(),
                Value::Array(vec![
                    ph.map_or(Value::Null, |c| Value::String(c.as_str().to_owned())),
                    cd.map_or(Value::Null, |c| Value::String(c.as_str().to_owned())),
                ]),
            );
        }
        m.insert(s.to_string(), Value::Object(bm));
    }
    Value::Object(m)
}

/// DD iono-free static baseline (fixed where the integers pass, else float) from `base` (known ARP
/// `base_xyz`) to `rover`. Returns the result and every stage.
pub fn solve(
    base: &StationObs,
    rover: &StationObs,
    orb: &Orbits,
    base_xyz: [f64; 3],
    opts_in: &Options,
    mut antex: Option<&mut Antex>,
    fcn: &BTreeMap<Sat, i32>,
    rover_xyz0: Option<[f64; 3]>,
    tropo: &TropoInputs,
) -> Result<(StaticResult, Stages), EngineError> {
    let t_start = Instant::now();
    let mut opts = opts_in.clone();
    let mut notes: Vec<String> = Vec::new();
    let mut stages = Stages::default();
    let timing = std::env::var_os("GEOM_ENGINE_TIMING").is_some();
    let mut t_lap = Instant::now();
    let mut lap = |what: &str| {
        if timing {
            eprintln!("timing: {what:<14} {:7.3} s", t_lap.elapsed().as_secs_f64());
            t_lap = Instant::now();
        }
    };
    let (glo_ok, glo_why) = glo_fixable(&opts, &base.receiver, &rover.receiver);
    opts.fix_systems = opts.fix_systems.replace('R', "") + if glo_ok { "R" } else { "" };

    // common epochs: the rover's tags with a base epoch within epoch_tol
    let rb: Vec<f64> = base.t.iter().map(|&x| np_round(x, 3)).collect();
    let rr: Vec<f64> = rover.t.iter().map(|&x| np_round(x, 3)).collect();
    let bset: BTreeSet<u64> = rb.iter().map(|x| x.to_bits()).collect();
    let mut t: Vec<f64> = rr
        .iter()
        .copied()
        .filter(|x| bset.contains(&x.to_bits()))
        .collect();
    t.dedup();
    if opts.epoch_tol > 1e-3 && base.t.len() > 1 && !rover.t.is_empty() {
        let j = nearest(&base.t, &rover.t);
        t = rover
            .t
            .iter()
            .zip(&j)
            .filter(|(x, jj)| (base.t[**jj] - **x).abs() < opts.epoch_tol)
            .map(|(&x, _)| np_round(x, 3))
            .collect();
    }
    if opts.interval > 0.0 {
        t.retain(|&x| ((x + 1e-3).rem_euclid(opts.interval)).abs() < 2e-3);
    }
    if t.len() < 10 {
        return Err(EngineError(
            "the two sessions share fewer than 10 epochs".into(),
        ));
    }
    let ne = t.len();

    // satellites and their signal pairs
    let mut cols: Vec<Sat> = Vec::new();
    let mut pairs: BTreeMap<Sat, (u8, u8)> = BTreeMap::new();
    let mut no_orbit: BTreeSet<Sat> = BTreeSet::new();
    let common: BTreeSet<Sat> = base
        .sats
        .keys()
        .filter(|s| rover.sats.contains_key(*s))
        .copied()
        .collect();
    for sat in common {
        if !opts.systems.contains(sat.sys.letter()) || opts.exclude.contains(&sat) {
            continue;
        }
        if !orb.has(sat) {
            no_orbit.insert(sat);
            continue;
        }
        if sat.sys == GnssSystem::Glonass && !fcn.contains_key(&sat) {
            continue;
        }
        let Some(pr) = pick_pair(sat, &[base, rover]) else {
            continue;
        };
        cols.push(sat);
        pairs.insert(sat, pr);
    }
    let col_sys: BTreeSet<char> = cols.iter().map(|c| c.sys.letter()).collect();
    let gone: Vec<char> = no_orbit
        .iter()
        .map(|s| s.sys.letter())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .filter(|s| !col_sys.contains(s))
        .collect();
    if !gone.is_empty() {
        notes.push(format!(
            "{} left out: no orbit for any of their satellites (broadcast or precise)",
            gone.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    let part: Vec<String> = no_orbit
        .iter()
        .filter(|s| !gone.contains(&s.sys.letter()))
        .map(ToString::to_string)
        .collect();
    if !part.is_empty() {
        notes.push(format!(
            "no orbit (broadcast or precise) for {}: left out",
            part.join(" ")
        ));
    }
    if !cols.is_empty() {
        let (mut acc_gone, mut acc_part) = (Vec::new(), Vec::new());
        for sat in cols.clone() {
            let (n_bad, n_ok, worst) = orb.accuracy_out(sat, t[0], t[ne - 1]);
            if n_bad > 0 && n_ok == 0 {
                acc_gone.push(format!(
                    "{sat} ({} m)",
                    geomaticist_gnss_core::text::format_g(worst)
                ));
                cols.retain(|c| *c != sat);
                pairs.remove(&sat);
            } else if n_bad > 0 {
                acc_part.push(format!(
                    "{sat} ({n_bad} at {} m)",
                    geomaticist_gnss_core::text::format_g(worst)
                ));
            }
        }
        if !acc_gone.is_empty() {
            notes.push(format!(
                "left out, its broadcast accuracy says unusable (URA index 15 / SISA NAPA): {}",
                acc_gone.join(", ")
            ));
        }
        if !acc_part.is_empty() {
            notes.push(format!(
                "ephemerides whose broadcast accuracy says unusable not used: {}",
                acc_part.join(", ")
            ));
        }
    }
    let ts: Vec<f64> = sample_rows(ne, 240).iter().map(|&i| t[i]).collect();
    if !cols.is_empty() {
        let mut stale = Vec::new();
        for &sat in &cols {
            let (frac, age) = orb.orbit_gap(sat, &ts);
            if frac > 0.0 {
                let share = if frac >= 0.995 {
                    "all".to_owned()
                } else {
                    format!("{:.0} %", 100.0 * frac)
                };
                stale.push(format!(
                    "{sat} ({share} of the time, nearest {:.0} min away)",
                    age / 60.0
                ));
            }
        }
        if !stale.is_empty() {
            notes.push(format!(
                "no usable broadcast ephemeris for {}: those epochs left out",
                stale.join(", ")
            ));
        }
    }
    let mut eph_age: BTreeMap<Sat, f64> = BTreeMap::new();
    if !cols.is_empty() {
        for &sat in &cols {
            if let Some(a) = orb.age_used(sat, &ts)
                && a.iter().any(|v| v.is_finite())
            {
                let mut best = 0usize;
                let mut bv = f64::NEG_INFINITY;
                for (i, v) in a.iter().enumerate() {
                    if v.is_finite() && v.abs() > bv {
                        bv = v.abs();
                        best = i;
                    }
                }
                eph_age.insert(
                    sat,
                    geomaticist_gnss_core::text::round_to(a[best] / 60.0, 1),
                );
            }
        }
        let mut worst: BTreeMap<char, (Sat, f64)> = BTreeMap::new();
        for (&sat, &a) in &eph_age {
            let s = sat.sys.letter();
            match worst.get(&s) {
                Some((_, w)) if a.abs() <= w.abs() => {}
                _ => {
                    worst.insert(s, (sat, a));
                }
            }
        }
        let mut old: Vec<(Sat, f64)> = worst.values().copied().collect();
        old.sort_by_key(|a| a.0);
        let old: Vec<String> = old
            .iter()
            .filter(|(_, a)| a.abs() > 60.0)
            .map(|(s, a)| format!("{s} {a:+.0} min"))
            .collect();
        if !old.is_empty() {
            notes.push(format!(
                "broadcast ephemeris age used (t − toe) up to {}",
                old.join(", ")
            ));
        }
    }
    if cols.is_empty() {
        return Err(EngineError("no dual-frequency satellites in common".into()));
    }
    // a receiver that moved (a rover's kinematic log): refused here, in well under a second, instead of minutes
    // of arcs split at every move and thousands of ambiguities that end in a failed float
    for (st, role) in [(rover, "rover"), (base, "base")] {
        if let Some(w) = wander(st, orb, &opts.systems, fcn)
            && w.moved()
        {
            return Err(EngineError(format!(
                "the {role} ({}) moved during the session: its single-point positions spread over {:.0} m ({:.0} % of the sampled epochs more than {MOVING_M:.0} m from their median). A moving (kinematic) log can't be processed as a static baseline; process only the times its antenna stood on a point",
                st.name,
                w.max_m,
                100.0 * w.share_far
            )));
        }
    }
    lap("setup:orbits");
    let cal_a = receiver_cal(antex.as_deref_mut(), base, &mut notes);
    let cal_b = receiver_cal(antex.as_deref_mut(), rover, &mut notes);
    lap("setup:antex");
    if opts.cal_rule == "own" && antex.is_some() {
        let mut out_sys = Vec::new();
        for s in col_sys.iter().filter(|s| opts.cal_systems.contains(**s)) {
            if !cols.iter().any(|c| c.sys.letter() == *s) {
                continue;
            }
            let bands: BTreeSet<u8> = cols
                .iter()
                .filter(|c| c.sys.letter() == *s)
                .flat_map(|c| [pairs[c].0, pairs[c].1])
                .collect();
            let bands: Vec<u8> = bands.into_iter().collect();
            if !system_calibrated(
                &[cal_a.as_ref(), cal_b.as_ref()],
                &[&base.antenna, &rover.antenna],
                *s,
                &bands,
            ) {
                out_sys.push(*s);
            }
        }
        if !out_sys.is_empty() {
            cols.retain(|c| !out_sys.contains(&c.sys.letter()));
            notes.push(format!(
                "{} left out: no receiver calibration of their frequencies at both ends (cal_rule own)",
                out_sys.iter().map(ToString::to_string).collect::<Vec<_>>().join(", ")
            ));
        }
        if cols.is_empty() {
            return Err(EngineError(
                "no calibrated dual-frequency satellites in common".into(),
            ));
        }
    }
    let ns = cols.len();
    let sysc: Vec<char> = cols.iter().map(|c| c.sys.letter()).collect();
    let sun = sun_ecef(&t);
    let mut x_r = rover_xyz0.unwrap_or(rover.approx);
    if norm3(x_r) < 6e6 || norm3(sub3(x_r, base_xyz)) > 2e6 {
        x_r = spp(rover, orb, &opts.systems, fcn, 20)?;
        notes.push("rover start from a single-point solution".into());
    }
    let mut sig_a: BTreeMap<Sat, BandSignals> = BTreeMap::new();
    let mut sig_b: BTreeMap<Sat, BandSignals> = BTreeMap::new();
    if opts.same_signals {
        for &sat in &cols {
            let cs = common_signals(sat, pairs[&sat], &[base, rover], 0.5);
            sig_a.insert(sat, cs[0].clone());
            sig_b.insert(sat, cs[1].clone());
        }
    }
    let ms_on = opts.ms_repair && (t[ne - 1] - t[0]) <= opts.ms_max_h * 3600.0 + 1.0;
    let mut o_prep = opts.clone();
    o_prep.ms_repair = ms_on;
    // satellite antenna calibrations at t[0]: only precise (centre-of-mass) orbits need them (reading the
    // satellite blocks of igs20.atx costs ~0.1 s)
    let mut sat_cals: BTreeMap<Sat, Option<AntennaCal>> = BTreeMap::new();
    let any_precise = orb.sources.iter().any(OrbitSource::center_of_mass);
    if let Some(atx) = antex
        && any_precise
    {
        for &sat in &cols {
            sat_cals.insert(sat, atx.satellite(&sat.to_string(), t[0]));
        }
    }
    lap("setup:rest");
    let mut a = prepare(
        base,
        &t,
        &cols,
        &pairs,
        fcn,
        orb,
        &sun,
        &o_prep,
        base_xyz,
        &sat_cals,
        opts.same_signals.then_some(&sig_a),
    );
    let mut b = prepare(
        rover,
        &t,
        &cols,
        &pairs,
        fcn,
        orb,
        &sun,
        &o_prep,
        x_r,
        &sat_cals,
        opts.same_signals.then_some(&sig_b),
    );
    for (st, s) in [(base, &a), (rover, &b)] {
        if s.events.code_ms > 0 || s.events.pc_jump > 0 || s.events.ms_repaired > 0 {
            notes.push(format!(
                "{}: {} codes off the receiver clock (dropped), {} whole-millisecond phase steps put back, {} other phase-code jumps (new arcs)",
                st.name, s.events.code_ms, s.events.ms_repaired, s.events.pc_jump
            ));
        }
    }
    let tide = |xyz: [f64; 3], t: &[f64]| -> Vec<[f64; 3]> {
        match opts.compat.tide {
            TideModel::Degree2Step1 => crate::models::solid_tide_degree2(xyz, t),
            TideModel::IersFull => crate::tides::solid_tide(xyz, t),
        }
    };
    if opts.tides {
        a.tide = tide(base_xyz, &t);
        b.tide = tide(x_r, &t);
    }
    let pos_a: Vec<[f64; 3]> = a.tide.iter().map(|td| add3(base_xyz, *td)).collect();
    lap("prepare");
    let ca = corrections(
        &a,
        &pos_a,
        &cols,
        &pairs,
        fcn,
        cal_a.as_ref(),
        &t,
        &opts,
        tropo,
    );
    let pos_b: Vec<[f64; 3]> = b.tide.iter().map(|td| add3(x_r, *td)).collect();
    let cb = corrections(
        &b,
        &pos_b,
        &cols,
        &pairs,
        fcn,
        cal_b.as_ref(),
        &t,
        &opts,
        tropo,
    );
    let mask = opts.elmask.to_radians();
    let mut valid = B2::falses(ne, ns);
    for e in 0..ne {
        for j in 0..ns {
            valid.set(
                e,
                j,
                a.ok.at(e, j) && b.ok.at(e, j) && cb.el.at(e, j) >= mask && ca.el.at(e, j) >= mask,
            );
        }
    }
    let sdgf = b.gf.zip_map(&a.gf, |x, y| x - y);
    let sdmw = b.mw.zip_map(&a.mw, |x, y| x - y);
    let (mut lam1, mut lam2, mut aif, mut bif) =
        (vec![0.0; ns], vec![0.0; ns], vec![0.0; ns], vec![0.0; ns]);
    for (j, sat) in cols.iter().enumerate() {
        let f = fcn.get(sat).copied();
        let f1 = freq(*sat, pairs[sat].0, f).map_err(|e| EngineError(e.to_string()))?;
        let f2 = freq(*sat, pairs[sat].1, f).map_err(|e| EngineError(e.to_string()))?;
        lam1[j] = CLIGHT / f1;
        lam2[j] = CLIGHT / f2;
        let (x, y) = if_coeffs(f1, f2);
        aif[j] = x;
        bif[j] = y;
    }
    let lli_sd = a.lli.or(&b.lli);
    let hard_sd = a.hard.or(&b.hard);
    stages.observations = Some(ObservationsStage {
        t: t.clone(),
        cols: cols.clone(),
        ph_a: a.phase.clone(),
        ph_b: b.phase.clone(),
        cd_a: a.code.clone(),
        cd_b: b.code.clone(),
        gf_a: a.gf.clone(),
        gf_b: b.gf.clone(),
        mw_a: a.mw.clone(),
        mw_b: b.mw.clone(),
        ok_a: a.ok.clone(),
        ok_b: b.ok.clone(),
        valid: valid.clone(),
        lli: lli_sd.clone(),
        hard: hard_sd.clone(),
        half_a: a.half.clone(),
        half_b: b.half.clone(),
        ms_a: a.ms.clone(),
        ms_b: b.ms.clone(),
        el_a: ca.el.clone(),
        el_b: cb.el.clone(),
        az_a: ca.az.clone(),
        az_b: cb.az.clone(),
        rho_a: ca.rho.clone(),
        rho_b: cb.rho.clone(),
        sat_a: a.sat_xyz.clone(),
        sat_b: b.sat_xyz.clone(),
        clk_a: a.sat_clk.clone(),
        clk_b: b.sat_clk.clone(),
        tide_a: a.tide.clone(),
        tide_b: b.tide.clone(),
        dry_a: ca.dry.clone(),
        dry_b: cb.dry.clone(),
        map_a: ca.map.clone(),
        map_b: cb.map.clone(),
        ant_a: ca.ant.clone(),
        ant_b: cb.ant.clone(),
        lam1: lam1.clone(),
        lam2: lam2.clone(),
        aif: aif.clone(),
        bif: bif.clone(),
        base: base.name.clone(),
        rover: rover.name.clone(),
        base_xyz,
        rover_xyz0: x_r,
        pairs: pairs.clone(),
        signals_base: signals_json(&sig_a),
        signals_rover: signals_json(&sig_b),
        eph_age_min: eph_age.clone(),
        events_base: a.events,
        events_rover: b.events,
        notes: notes.clone(),
    });

    lap("corrections");
    // slips
    let mut corr = F2::zeros(ne, ns);
    let mut corr_mw = F2::zeros(ne, ns);
    let mut corr_gf = F2::zeros(ne, ns);
    let mut starts = B2::falses(ne, ns);
    let mut wl_cont = B2::falses(ne, ns);
    let mut slip_ev = SlipEvents::default();
    let mut slip_ev_gf: BTreeMap<String, usize> = BTreeMap::new();
    let mut sat_slips: BTreeMap<Sat, (usize, usize)> = BTreeMap::new();
    if opts.slip_mode == "gf" {
        let mut slip = a.slip.or(&b.slip);
        for j in 0..ns {
            let have: Vec<usize> = (0..ne)
                .filter(|&e| sdgf.at(e, j).is_finite() && a.ok.at(e, j) && b.ok.at(e, j))
                .collect();
            for k in 1..have.len() {
                if (sdgf.at(have[k], j) - sdgf.at(have[k - 1], j)).abs() > opts.slip_gf {
                    slip.set(have[k], j, true);
                }
            }
        }
        let ll = a.slip.or(&b.slip);
        let (mut n_start, mut n_lli, mut n_gf, mut n_gap) = (0, 0, 0, 0);
        for j in 0..ns {
            let ep: Vec<usize> = (0..ne).filter(|&e| valid.at(e, j)).collect();
            if ep.is_empty() {
                continue;
            }
            starts.set(ep[0], j, true);
            for k in 1..ep.len() {
                let gap = ep[k] - ep[k - 1] > opts.max_gap + 1;
                let sl = (ep[k - 1] + 1..=ep[k]).any(|e| slip.at(e, j));
                let l = (ep[k - 1] + 1..=ep[k]).any(|e| ll.at(e, j));
                if gap || sl {
                    starts.set(ep[k], j, true);
                }
                n_gap += usize::from(gap);
                n_lli += usize::from(!gap && l);
                n_gf += usize::from(!gap && sl && !l);
            }
            n_start += 1;
        }
        slip_ev_gf.insert("start".into(), n_start);
        slip_ev_gf.insert("LLI".into(), n_lli);
        slip_ev_gf.insert("GF".into(), n_gf);
        slip_ev_gf.insert("gap".into(), n_gap);
    } else {
        for j in 0..ns {
            let r = slips::process(
                &t,
                &sdgf.col(j),
                &sdmw.col(j),
                &valid.col(j),
                &lli_sd.col(j),
                &hard_sd.col(j),
                lam1[j],
                lam2[j],
                aif[j],
                bif[j],
                &opts.slips,
            );
            for e in 0..ne {
                starts.set(e, j, r.new_arc[e]);
                wl_cont.set(e, j, r.wl_link[e]);
                if r.drop[e] {
                    valid.set(e, j, false);
                }
                corr.set(e, j, r.corr_if[e]);
                corr_mw.set(e, j, r.corr_mw[e]);
                corr_gf.set(e, j, r.corr_gf[e]);
            }
            slip_ev.add(&r.events);
            sat_slips.insert(cols[j], (r.repairs.len(), r.events.spike));
        }
    }
    lap("slips");
    // arcs
    let ms_sd = if ms_on && opts.ms_split {
        a.ms.or(&b.ms)
    } else {
        B2::falses(ne, ns)
    };
    starts = starts.or(&ms_sd);
    let half_sd = a.half.or(&b.half);
    if half_sd.any() {
        for j in 0..ns {
            if !half_sd.any_col(j) {
                continue;
            }
            let ep: Vec<usize> = (0..ne).filter(|&e| valid.at(e, j)).collect();
            for k in 1..ep.len() {
                if half_sd.at(ep[k], j) != half_sd.at(ep[k - 1], j) {
                    starts.set(ep[k], j, true);
                }
            }
        }
    }
    let mut arc = I2::filled(ne, ns, -1);
    let mut n_arc = 0usize;
    for j in 0..ns {
        let ep: Vec<usize> = (0..ne).filter(|&e| valid.at(e, j)).collect();
        if ep.is_empty() {
            continue;
        }
        let mut id = n_arc as i64 - 1;
        for k in 0..ep.len() {
            let st = k == 0 || starts.at(ep[k], j) || ep[k] - ep[k - 1] > opts.max_gap + 1;
            if st {
                id += 1;
            }
            arc.set(ep[k], j, id);
        }
        n_arc = id as usize + 1;
    }
    let mut counts = vec![0usize; n_arc];
    for &v in &arc.data {
        if v >= 0 {
            counts[v as usize] += 1;
        }
    }
    let short: Vec<bool> = counts.iter().map(|&c| c < opts.min_arc).collect();
    let mut remap = vec![-1i64; n_arc];
    let mut n_amb = 0usize;
    for i in 0..n_arc {
        if !short[i] {
            remap[i] = n_amb as i64;
            n_amb += 1;
        }
    }
    for v in &mut arc.data {
        if *v >= 0 {
            *v = remap[*v as usize];
        }
    }
    let mut wl_links: Vec<(usize, usize)> = Vec::new();
    for e0 in 0..ne {
        for j in 0..ns {
            if (ms_sd.at(e0, j) || wl_cont.at(e0, j)) && arc.at(e0, j) >= 0 {
                let prev = (0..e0).rev().find(|&e| arc.at(e, j) >= 0);
                if let Some(pe) = prev
                    && arc.at(pe, j) != arc.at(e0, j)
                {
                    wl_links.push((arc.at(pe, j) as usize, arc.at(e0, j) as usize));
                }
            }
        }
    }
    let n_short = short.iter().filter(|&&s| s).count();
    let ev_text = if opts.slip_mode == "gf" {
        ["start", "LLI", "GF", "gap"]
            .iter()
            .map(|k| format!("{k} {}", slip_ev_gf.get(*k).copied().unwrap_or(0)))
            .collect::<Vec<_>>()
            .join(", ")
    } else {
        slip_ev
            .pairs()
            .iter()
            .map(|(k, v)| format!("{k} {v}"))
            .collect::<Vec<_>>()
            .join(", ")
    };
    notes.push(format!(
        "arcs: {n_arc} ({ev_text}); {n_short} shorter than {} epochs dropped{}",
        opts.min_arc,
        if wl_links.is_empty() {
            String::new()
        } else {
            format!("; {} pieces linked for the wide lane", wl_links.len())
        }
    ));
    valid = valid.zip_map(&arc, |v, a| v && a >= 0);
    stages.arcs = Some(ArcsStage {
        arc: arc.clone(),
        starts: starts.clone(),
        valid: valid.clone(),
        corr: corr.clone(),
        corr_mw: corr_mw.clone(),
        corr_gf: corr_gf.clone(),
        half: half_sd.clone(),
        ms: ms_sd.clone(),
        n_amb,
        n_arc_raw: n_arc,
        wl_links: wl_links.clone(),
        slip_events: slip_ev,
        sat_slips: sat_slips.clone(),
        notes: notes.clone(),
    });
    // ambiguities random-walking in pieces (amb_rw > 0)
    let mut links: Vec<(usize, usize, f64)> = Vec::new();
    if opts.amb_rw > 0.0 && n_amb > 0 {
        let mut seg = I2::filled(ne, ns, -1);
        let mut nid = 0usize;
        for a_id in 0..n_amb {
            let mut cells: Vec<(usize, usize)> = Vec::new();
            for e in 0..ne {
                for j in 0..ns {
                    if arc.at(e, j) == a_id as i64 {
                        cells.push((e, j));
                    }
                }
            }
            let t0 = t[cells[0].0];
            let pieces: Vec<i64> = cells
                .iter()
                .map(|(e, _)| ((t[*e] - t0) / opts.amb_seg_s).floor() as i64)
                .collect();
            let mut uniq = pieces.clone();
            uniq.sort_unstable();
            uniq.dedup();
            for (k, (e, j)) in cells.iter().enumerate() {
                let inv = uniq.iter().position(|&u| u == pieces[k]).unwrap_or(0);
                seg.set(*e, *j, (nid + inv) as i64);
            }
            for k in 0..uniq.len().saturating_sub(1) {
                links.push((
                    nid + k,
                    nid + k + 1,
                    opts.amb_seg_s * (uniq[k + 1] - uniq[k]) as f64,
                ));
            }
            nid += uniq.len();
        }
        arc = seg;
        n_amb = nid;
    }

    // parameters other than the ambiguities
    let step = opts.trop_step_h * 3600.0;
    let mut nn = if opts.trop == "none" {
        0
    } else {
        ((t[ne - 1] - t[0]) / step - 1e-9).ceil() as usize + 1
    };
    if opts.trop != "none" {
        nn = nn.max(2);
    }
    let i_tb = 3;
    let i_ta = 3 + nn;
    let n_ta = if opts.trop == "station" { nn } else { 0 };
    let i_amb = 3 + nn + n_ta;
    let kk: Vec<usize> = t
        .iter()
        .map(|&x| (((x - t[0]) / step).floor() as usize).min(nn.saturating_sub(2)))
        .collect();
    let ww: Vec<f64> = if nn > 0 {
        t.iter()
            .zip(&kk)
            .map(|(&x, &k)| (x - (t[0] + k as f64 * step)) / step)
            .collect()
    } else {
        vec![0.0; ne]
    };
    let ntr = if nn > 0 { 2 } else { 0 } + if n_ta > 0 { 2 } else { 0 };
    let mut var_ph =
        variance(&ca.el, &sysc, &opts).zip_map(&variance(&cb.el, &sysc, &opts), |x, y| x + y);
    if half_sd.any() {
        for e in 0..ne {
            for j in 0..ns {
                let hv = opts.half_var
                    * (aif[j] * aif[j] + bif[j] * bif[j])
                    * (f64::from(u8::from(a.half.at(e, j))) + f64::from(u8::from(b.half.at(e, j))));
                let v = var_ph.at(e, j) + hv;
                var_ph.set(e, j, v);
            }
        }
    }
    let mut w_rel = 0.0;
    if opts.trop_rel > 0.0 && n_ta > 0 {
        let (_, _, ha) = ecef_to_llh(base_xyz);
        let (_, _, hb) = ecef_to_llh(x_r);
        let s_rel = opts
            .trop_rel
            .hypot(opts.trop_rel_ppm * 1e-6 * norm3(sub3(x_r, base_xyz)))
            .hypot(opts.trop_init * (hb - ha).abs() / 2000.0);
        w_rel = 1.0 / (s_rel * s_rel);
    }
    let var_cd =
        (opts.code_ratio > 0.0).then(|| var_ph.map(|v| v * opts.code_ratio * opts.code_ratio));
    let obs_ph = b
        .phase
        .zip_map(&a.phase, |x, y| x - y)
        .zip_map(&corr, |x, c| x + c);
    let obs_cd = b.code.zip_map(&a.code, |x, y| x - y);
    let systems: BTreeSet<char> = sysc.iter().copied().collect();
    let mut groups_all: Vec<(usize, Vec<usize>)> = Vec::new();
    for e in 0..ne {
        for s in &systems {
            let cj: Vec<usize> = (0..ns)
                .filter(|&j| valid.at(e, j) && sysc[j] == *s)
                .collect();
            if cj.len() >= 2 {
                groups_all.push((e, cj));
            }
        }
    }
    if groups_all.is_empty() {
        return Err(EngineError(
            "no double differences (too few common satellites above the mask)".into(),
        ));
    }
    let n_sys = systems.len();
    let mut ctx = Ctx {
        opts: opts.clone(),
        t: t.clone(),
        ne,
        ns,
        cols: cols.clone(),
        sysc: sysc.clone(),
        pairs: pairs.clone(),
        fcn: fcn.clone(),
        base_xyz,
        x_r,
        a,
        b,
        ca,
        cb,
        at_b: x_r,
        cal_b,
        tropo,
        lam1,
        lam2,
        aif,
        bif,
        valid,
        var_ph,
        var_cd,
        obs_ph,
        obs_cd,
        sdgf,
        corr_gf,
        sdmw,
        corr_mw,
        half_sd,
        nn,
        i_tb,
        i_ta,
        n_ta,
        i_amb,
        kk,
        ww,
        ntr,
        step,
        w_rel,
        groups_all,
        n_sys,
        rej_cd0: B2::falses(ne, ns),
    };
    // gross code errors first
    let mut n_gross_code = 0usize;
    let mut x0 = vec![0.0; 3 + nn + n_ta];
    x0[..3].copy_from_slice(&x_r);
    if nn > 0 {
        let v = if opts.trop == "station" {
            opts.trop_init
        } else {
            0.0
        };
        for k in 0..nn {
            x0[i_tb + k] = v;
        }
        for k in 0..n_ta {
            x0[i_ta + k] = opts.trop_init;
        }
    }
    let rho_b0 = ctx.cb.rho.clone();
    if opts.phase_gross > 0.0 && n_amb > 0 {
        let m = ctx.model_sd(&x0, &rho_b0);
        let r = ctx.obs_ph.zip_map(&m, |x, y| x - y);
        let ok = ctx.valid.zip_map(&arc, |v, a| v && a >= 0);
        let bad = phase_steps(&r, &ok, &arc, &t, &opts);
        if bad.any() {
            let bad_cols: Vec<String> = (0..ns)
                .filter(|&j| bad.any_col(j))
                .map(|j| cols[j].to_string())
                .collect();
            for i in 0..bad.data.len() {
                if bad.data[i] {
                    ctx.valid.data[i] = false;
                    arc.data[i] = -1;
                }
            }
            notes.push(format!(
                "{} phase observations after a step of over {:.0} m dropped ({})",
                bad.count(),
                opts.phase_gross,
                bad_cols.join(", ")
            ));
        }
    }
    if ctx.var_cd.is_some() {
        let m = ctx.model_sd(&x0, &rho_b0);
        let rc = ctx.obs_cd.zip_map(&m, |x, y| x - y);
        for (e, cj) in &ctx.groups_all {
            let cc: Vec<usize> = cj
                .iter()
                .copied()
                .filter(|&j| rc.at(*e, j).is_finite())
                .collect();
            if cc.len() >= 3 {
                let vals: Vec<f64> = cc.iter().map(|&j| rc.at(*e, j)).collect();
                let med = median(&vals);
                for (k, &j) in cc.iter().enumerate() {
                    if (vals[k] - med).abs() > opts.code_gross {
                        ctx.rej_cd0.set(*e, j, true);
                        n_gross_code += 1;
                    }
                }
            }
        }
    }
    lap("arcs+setup");
    let mut est = ctx.estimate(&arc, n_amb, &links, x_r, None)?;
    lap("estimate");
    let mut n_split = 0usize;
    if opts.step_k > 0.0 && links.is_empty() {
        for _rnd in 0..opts.step_rounds {
            let (wf, cf) = (&est.wf_ph, &est.cf_ph);
            let mut zz = Vec::new();
            for i in 0..cf.data.len() {
                if wf.data[i] > 0.0 && cf.data[i].is_finite() {
                    zz.push(cf.data[i].abs() / ctx.var_ph.data[i].sqrt());
                }
            }
            let s_rob = if zz.is_empty() {
                1.0
            } else {
                1.4826 * median(&zz)
            };
            let mut splits: Vec<(i64, usize, usize)> = Vec::new();
            // arcs in id order, each with its cells in row-major order
            let mut cells: Vec<(i64, usize, usize)> = Vec::new();
            for e in 0..ne {
                for j in 0..ns {
                    if est.use_.at(e, j) && cf.at(e, j).is_finite() && arc.at(e, j) >= 0 {
                        cells.push((arc.at(e, j), e, j));
                    }
                }
            }
            cells.sort_by_key(|c| c.0);
            let mut k = 0;
            while k < cells.len() {
                let a_id = cells[k].0;
                let mut hi = k;
                while hi < cells.len() && cells[hi].0 == a_id {
                    hi += 1;
                }
                let grp = &cells[k..hi];
                let r: Vec<f64> = grp.iter().map(|c| cf.at(c.1, c.2)).collect();
                let sig: Vec<f64> = grp
                    .iter()
                    .map(|c| ctx.var_ph.at(c.1, c.2).sqrt() * s_rob)
                    .collect();
                if let Some((q, _dlt, zz, zl, dl)) = find_step(&r, &sig, opts.min_arc.max(4), 6)
                    && zz > opts.step_k
                    && zl > opts.step_k_local
                    && dl.abs() > opts.step_min
                {
                    splits.push((a_id, grp[q].1, grp[0].2));
                }
                k = hi;
            }
            if splits.is_empty() {
                break;
            }
            for (a_id, e0, j) in &splits {
                for e in *e0..ne {
                    if arc.at(e, *j) == *a_id {
                        arc.set(e, *j, n_amb as i64);
                    }
                }
                n_amb += 1;
            }
            n_split += splits.len();
            let pos = [est.x[0], est.x[1], est.x[2]];
            est = ctx.estimate(&arc, n_amb, &links, pos, Some(&est))?;
        }
        if n_split > 0 {
            notes.push(format!("{n_split} arcs split where their residuals step"));
        }
    }
    lap("step splits");
    let biased = if ctx.var_cd.is_some() {
        code_biases(&est.res_cd, &est.code_ok, &t, &sysc, &opts)
    } else {
        BTreeMap::new()
    };
    if !biased.is_empty() {
        for &j in biased.keys() {
            for e in 0..ne {
                ctx.rej_cd0.set(e, j, true);
            }
        }
        let pos = [est.x[0], est.x[1], est.x[2]];
        est = ctx.estimate(&arc, n_amb, &links, pos, Some(&est))?;
        notes.push(format!(
            "code biased against the other satellites (SD iono-free): {}; their codes left out, their wide lanes float",
            biased.iter().map(|(j, v)| format!("{} {v:+.1} m", cols[*j])).collect::<Vec<_>>().join(", ")
        ));
    }
    let q = est
        .n
        .inverse()
        .ok_or_else(|| EngineError("singular normal matrix".into()))?;
    let s0 = est.s0;
    let x = est.x.clone();
    let vec_float = sub3([x[0], x[1], x[2]], base_xyz);
    let mut cov_float = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cov_float[i][j] = q.at(i, j) * s0 * s0;
        }
    }
    stages.float = Some(FloatStage {
        x: x.clone(),
        n: est.n.clone(),
        q: q.clone(),
        arc: arc.clone(),
        use_: est.use_.clone(),
        wf_ph: est.wf_ph.clone(),
        wf_cd: est.wf_cd.clone(),
        res_ph: est.res_ph.clone(),
        res_cd: est.res_cd.clone(),
        cf_ph: est.cf_ph.clone(),
        var_ph: ctx.var_ph.clone(),
        var_cd: ctx.var_cd.clone(),
        s0,
        vector: vec_float,
        cov: cov_float,
        n_amb,
        i_amb,
        i_tb,
        i_ta,
        nn,
        n_ta,
        p: est.p,
        acc: est.acc,
        biased: biased.clone(),
        notes: notes.clone(),
    });

    // fixing, validation, sigma
    lap("float+inverse");
    let out = crate::fixing::finish(
        &mut ctx,
        &mut stages,
        est,
        q,
        &arc,
        n_amb,
        &links,
        &biased,
        &wl_links,
        vec_float,
        cov_float,
        &glo_why,
        n_split,
        n_gross_code,
        &eph_age,
        &sat_slips,
        slip_ev,
        &slip_ev_gf,
        notes,
        t_start,
    )?;
    lap("fixing+checks");
    Ok((out, stages))
}

/// The arcs per satellite, repairs, spikes and ephemeris age for the report.
pub fn sat_report_entry(
    n: usize,
    rms_mm: Option<f64>,
    down: usize,
    zero: usize,
    arcs: usize,
    repaired: usize,
    spikes: usize,
    eph_age_min: Option<f64>,
) -> Value {
    let mut m = Map::new();
    m.insert("n".into(), Value::from(n));
    m.insert(
        "rms_mm".into(),
        rms_mm.map_or(Value::Null, |v| jround(v, 1)),
    );
    m.insert("down".into(), Value::from(down));
    m.insert("zero".into(), Value::from(zero));
    m.insert("arcs".into(), Value::from(arcs));
    m.insert("repaired".into(), Value::from(repaired));
    m.insert("spikes".into(), Value::from(spikes));
    m.insert("eph_age_min".into(), eph_age_min.map_or(Value::Null, jf));
    Value::Object(m)
}

/// Each value rounded to 4 decimals, for the troposphere report.
#[must_use]
pub fn rounded4(v: &[f64]) -> Vec<f64> {
    v.iter()
        .map(|x| geomaticist_gnss_core::text::round_to(*x, 4))
        .collect()
}

/// Distance `|a - b|` between 3-vectors.
#[must_use]
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3(sub3(a, b))
}

/// `R v` helper re-exported for the fixing module.
#[must_use]
pub fn enu_of(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    mat3_vec(r, v)
}

/// The arrays of a float list for JSON.
#[must_use]
pub fn jvec(v: &[f64]) -> Value {
    jfl(v)
}

/// `x` rounded half to even, for the dumps.
#[must_use]
pub fn round_i(x: f64) -> i64 {
    py_round_i(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Satellites standing still (enough for a code-only position).
    #[derive(Debug)]
    struct FixedSats(Vec<(Sat, [f64; 3])>);

    impl crate::orbits::ExternalOrbit for FixedSats {
        fn has(&self, sat: Sat) -> bool {
            self.0.iter().any(|s| s.0 == sat)
        }
        fn center_of_mass(&self) -> bool {
            false
        }
        fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
            let p = self.0.iter().find(|s| s.0 == sat).map_or([0.0; 3], |s| s.1);
            (vec![p; t.len()], vec![0.0; t.len()], vec![true; t.len()])
        }
    }

    /// A receiver whose position at epoch k is `at(k)`, 300 one-second epochs, 8 satellites; GPS L1/L2 code
    /// and phase computed from the geometry (no errors).
    fn walker(at: &dyn Fn(f64) -> [f64; 3]) -> (StationObs, Orbits) {
        use geomaticist_gnss_core::obs::{SatSignals, Signal};
        use geomaticist_gnss_core::signal::{ObsCode, code_key};
        let x0 = at(0.0);
        let r = enu_matrix(60f64.to_radians(), (-150f64).to_radians());
        let mut sats = Vec::new();
        for k in 0..8 {
            let az = f64::from(k) * 45f64.to_radians();
            let el = if k % 2 == 0 { 30f64 } else { 60f64 }.to_radians();
            let enu = [el.cos() * az.sin(), el.cos() * az.cos(), el.sin()];
            // ENU -> ECEF: rows of r are e, n, u
            let d: Vec<f64> = (0..3)
                .map(|i| r[0][i] * enu[0] + r[1][i] * enu[1] + r[2][i] * enu[2])
                .collect();
            let sat: Sat = format!("G{:02}", k + 1)
                .parse()
                .unwrap_or_else(|_| unreachable!());
            sats.push((
                sat,
                [
                    x0[0] + 2.2e7 * d[0],
                    x0[1] + 2.2e7 * d[1],
                    x0[2] + 2.2e7 * d[2],
                ],
            ));
        }
        let mut st = StationObs::new("WALK", "w.obs");
        st.approx = x0;
        st.t = (0..300).map(|k| 1.25e9 + f64::from(k)).collect();
        st.power_fail = vec![false; 300];
        for (sat, sp) in &sats {
            let rng: Vec<f64> = (0..300)
                .map(|k| norm3(sub3(*sp, at(f64::from(k)))))
                .collect();
            let mut sigs = SatSignals::new();
            for (order, code) in ["C1C", "L1C", "C2W", "L2W"].iter().enumerate() {
                let key = code_key(code.as_bytes(), 'G', 3, 3.04);
                let lam = CLIGHT
                    / if code.contains('1') {
                        1_575.42e6
                    } else {
                        1_227.60e6
                    };
                let values = rng
                    .iter()
                    .map(|v| if key.kind == 'L' { v / lam } else { *v })
                    .collect();
                let c: ObsCode = code.parse().unwrap_or_else(|_| unreachable!());
                sigs.insert(
                    c,
                    Signal {
                        key,
                        values,
                        lli: (key.kind == 'L').then(|| vec![0; 300]),
                        first_seen: u32::try_from(order).unwrap_or(0),
                    },
                );
            }
            st.sats.insert(*sat, sigs);
        }
        let orb = Orbits::new(vec![OrbitSource::External(Box::new(FixedSats(sats)))]);
        (st, orb)
    }

    #[test]
    fn a_static_receiver_does_not_wander() {
        let x0 = [4_395_000.0, 618_000.0, 4_566_000.0];
        let (st, orb) = walker(&|_| x0);
        let w = wander(&st, &orb, "G", &BTreeMap::new()).unwrap_or(Wander {
            n: 0,
            share_far: 1.0,
            max_m: 1e9,
        });
        assert!(w.n >= 100);
        assert!(w.max_m < 0.01, "{w:?}");
        assert!(!w.moved());
    }

    #[test]
    fn a_walking_receiver_has_moved() {
        // 300 s at 2 m/s along X: 600 m
        let (st, orb) = walker(&|k| [4_395_000.0 + 2.0 * k, 618_000.0, 4_566_000.0]);
        let w = wander(&st, &orb, "G", &BTreeMap::new()).unwrap_or(Wander {
            n: 0,
            share_far: 0.0,
            max_m: 0.0,
        });
        assert!(w.moved(), "{w:?}");
        assert!(w.max_m > 250.0);
    }

    /// A deterministic group: `nr` rows on `u` columns with the engine's layout (three rover columns, a
    /// diagonal per row, two more per row), small residuals plus a common offset (the SD clock).
    fn group(
        seed: u64,
        nr: usize,
        offset: f64,
        g0: usize,
    ) -> (Rows, Vec<usize>, Vec<f64>, Vec<f64>) {
        let mut st = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let mut rnd = || {
            st = st
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((st >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
        };
        let u = 3 + 3 * nr;
        let mut a = Rows::zeros(nr, u);
        for r in 0..nr {
            for c in 0..3 {
                a[(r, c)] = rnd();
            }
            a[(r, 3 + r)] = 1.0;
            a[(r, 3 + nr + 2 * r)] = -0.7 - rnd();
            a[(r, 3 + nr + 2 * r + 1)] = -0.3 + rnd();
        }
        let g: Vec<usize> = (0..3).chain((0..3 * nr).map(|k| g0 + k)).collect();
        let var: Vec<f64> = (0..nr).map(|_| 2e-6 * (1.0 + rnd().abs())).collect();
        let y: Vec<f64> = (0..nr).map(|_| offset + 0.004 * rnd()).collect();
        (a, g, var, y)
    }

    #[test]
    fn rover_block_matches_block_bit_for_bit() {
        for (seed, nr) in [(1, 2), (2, 5), (3, 9)] {
            let (a, g, var, y) = group(seed, nr, 30.0, 10);
            let full = block(&a, g.clone(), &var, &y, None, VtpvForm::Original);
            let rov = rover_block(&a, g.clone(), &var, &y, VtpvForm::Original);
            let u = g.len();
            for c1 in 0..ROVER {
                for c2 in 0..u {
                    assert_eq!(
                        full.m[c1 * u + c2].to_bits(),
                        rov.rows[c1 * u + c2].to_bits()
                    );
                }
            }
            for c1 in ROVER..u {
                for c2 in 0..ROVER {
                    assert_eq!(
                        full.m[c1 * u + c2].to_bits(),
                        rov.cols[(c1 - ROVER) * ROVER + c2].to_bits()
                    );
                }
            }
            for c1 in 0..u {
                assert_eq!(full.b[c1].to_bits(), rov.b[c1].to_bits());
            }
            assert_eq!(full.vtpv.to_bits(), rov.vtpv.to_bits());
            assert_eq!(full.nrow, rov.nrow);
        }
    }

    #[test]
    fn block_matches_the_per_element_formula() {
        let (a, g, var, y) = group(7, 6, 30.0, 4);
        let blk = block(&a, g.clone(), &var, &y, Some('p'), VtpvForm::Original);
        let u = g.len();
        let w: Vec<f64> = var.iter().map(|v| 1.0 / v).collect();
        let sw = np_sum(&w);
        let ca: Vec<f64> = (0..u)
            .map(|c| (0..a.nr).fold(0.0, |s, r| s + a[(r, c)] * w[r]))
            .collect();
        let wy: f64 = (0..a.nr).map(|r| w[r] * y[r]).sum();
        for c1 in 0..u {
            for c2 in 0..u {
                let s = (0..a.nr).fold(0.0, |s, r| s + a[(r, c1)] * (a[(r, c2)] * w[r]));
                assert_eq!(
                    blk.m[c1 * u + c2].to_bits(),
                    (s - ca[c1] * ca[c2] / sw).to_bits()
                );
            }
            let s = (0..a.nr).fold(0.0, |s, r| s + a[(r, c1)] * w[r] * y[r]);
            assert_eq!(blk.b[c1].to_bits(), (s - ca[c1] * (wy / sw)).to_bits());
        }
        let ywy: f64 = (0..a.nr).map(|r| y[r] * (w[r] * y[r])).sum();
        assert_eq!(blk.vtpv.to_bits(), (ywy - wy * wy / sw).max(0.0).to_bits());
        assert_eq!(blk.nrow, a.nr - 1);
    }

    #[test]
    fn assemble_rows_equals_the_sequential_scatter() {
        // enough groups for several rounds and a chunk tail; parameters shared between groups
        let p = 97;
        let groups: Vec<(usize, Vec<usize>)> = (0..3 * GROUP_ROUND + 17)
            .map(|e| (e, (0..(2 + e % 5)).collect()))
            .collect();
        let design = |e: usize, cj: &[usize]| {
            let nr = cj.len();
            let (a, _, var, y) = group(e as u64 + 1, nr, 2.0e5, 0);
            let g: Vec<usize> = (0..3)
                .chain((0..3 * nr).map(|k| 3 + (e * 7 + k * 11) % (p - 3)))
                .collect();
            (a, g, var, y)
        };
        let key = |e: usize| {
            if e.is_multiple_of(2) {
                Some('p')
            } else {
                Some('c')
            }
        };
        let mut n_seq = DMat::zeros(p);
        let mut b_seq = vec![0.0; p];
        let mut acc_seq = Acc::default();
        for (e, cj) in &groups {
            let (a, g, var, y) = design(*e, cj);
            let blk = block(&a, g, &var, &y, key(*e), VtpvForm::Original);
            scatter(&mut n_seq, &mut b_seq, &mut acc_seq, &blk);
        }
        let mut n_par = DMat::zeros(p);
        let mut b_par = vec![0.0; p];
        let mut acc_par = Acc::default();
        assemble_blocks(
            &groups,
            &mut n_par,
            &mut b_par,
            &mut acc_par,
            |e, cj, _sc, out| {
                let (a, g, var, y) = design(e, cj);
                out.push(block(&a, g, &var, &y, key(e), VtpvForm::Original));
            },
        );
        assert!(
            n_seq
                .data
                .iter()
                .zip(&n_par.data)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert!(
            b_seq
                .iter()
                .zip(&b_par)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert_eq!(acc_seq, acc_par);
        // the rover blocks through the same scatter, against their own sequential sum
        let mut n_seq = DMat::zeros(p);
        let mut b_seq = vec![0.0; p];
        let mut acc_seq = Acc::default();
        for (e, cj) in &groups {
            let (a, g, var, y) = design(*e, cj);
            let blk = rover_block(&a, g, &var, &y, VtpvForm::Original);
            for c1 in 0..blk.g.len() {
                let i = blk.g[c1];
                let (row, b) = (&mut n_seq.data[i * p..(i + 1) * p], &mut b_seq[i]);
                blk.add_row(c1, row, b);
            }
            blk.add_acc(&mut acc_seq);
        }
        let mut n_par = DMat::zeros(p);
        let mut b_par = vec![0.0; p];
        let mut acc_par = Acc::default();
        assemble_rows(
            &groups,
            &mut n_par,
            &mut b_par,
            &mut acc_par,
            |e, cj, _sc, out| {
                let (a, g, var, y) = design(e, cj);
                out.push(rover_block(&a, g, &var, &y, VtpvForm::Original));
            },
        );
        assert!(
            n_seq
                .data
                .iter()
                .zip(&n_par.data)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert!(
            b_seq
                .iter()
                .zip(&b_par)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
        assert_eq!(acc_seq, acc_par);
    }

    #[test]
    fn vtpv_forms_agree_without_a_clock_and_only_the_stable_one_holds_with_it() {
        let (a, g, var, y) = group(11, 8, 0.0, 3);
        let quiet = block(&a, g.clone(), &var, &y, None, VtpvForm::Original).vtpv;
        let stable = block(&a, g.clone(), &var, &y, None, VtpvForm::Stable).vtpv;
        assert!(quiet > 0.0 && ((quiet - stable) / stable).abs() < 1e-9);
        // the same residuals shifted by an SD receiver clock of 2.3e5 m: the quantity is unchanged
        let yc: Vec<f64> = y.iter().map(|v| v + 2.3e5).collect();
        let stable_c = block(&a, g.clone(), &var, &yc, None, VtpvForm::Stable).vtpv;
        let original_c = block(&a, g, &var, &yc, None, VtpvForm::Original).vtpv;
        assert!(
            ((stable_c - stable) / stable).abs() < 1e-6,
            "{stable_c} vs {stable}"
        );
        // the original expression is the last bits of ~1e16 terms here: further off than the stable one
        assert!(
            ((original_c - stable) / stable).abs() > 1e-6,
            "{original_c} vs {stable}"
        );
    }
}
