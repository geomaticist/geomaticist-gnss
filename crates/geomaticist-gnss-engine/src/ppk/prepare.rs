//! From raw observations to what the filter needs at one epoch.
//!
//! For every rover epoch: the satellite positions at the transmit time, the base's observation corrections
//! (observed minus modelled) brought to the rover's epoch by linear interpolation between the base epochs
//! around it, a code-only position and clock of the rover, a Doppler velocity, and the cycle-slip flags.
//!
//! Interpolating the base works because what is interpolated is the *correction*, not the observation: with
//! the geometry and the satellite clock taken out, what is left changes slowly (atmosphere, multipath) or is
//! the same for every satellite (the base receiver's clock, which the double difference removes whatever
//! its shape between the two epochs).
//!
//! Cycle slips. A new ambiguity is started (the `slip` flag of a [`PrepBand`]):
//!
//! 1. when a satellite appears, or returns after a break longer than the outage limit
//!    ([`Setup::outage_s`]);
//! 2. when a receiver's geometry-free phase jumps (two frequencies);
//! 3. when the epoch-to-epoch change of a single-difference phase does not fit the change every other
//!    satellite agrees on (three position components and one clock solved from all of them; limit
//!    `slip_td`). This time-difference test works on one frequency and for any motion.
//!
//! A loss-of-lock indicator of either receiver, or a short break in one satellite's phase, makes the phase
//! *suspect* rather than slipped. A suspect phase is tested against the motion of the unflagged phases
//! since it was last seen (the per-epoch motions are summed across the break): no jump clears it, a jump
//! of whole cycles (within [`REPAIR_MAX_FRAC`]) is taken out of this and all later phases of that
//! satellite and the ambiguity carries on, anything else starts a new ambiguity. Receivers keep the
//! indicator up for several epochs after one loss of lock, and lose and regain high satellites within a
//! fraction of a second under dynamic stress; restarting those ambiguities each time would cost the
//! satellites that matter most to the geometry.

use std::collections::BTreeMap;

use geomaticist_gnss_core::consts::{CLIGHT, OMGE};
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};
use geomaticist_gnss_products::{AntennaCal, receiver_correction};

use super::options::PpkOptions;
use super::plan::Plan;
use super::{PpkError, PpkInput, Progress, Stage};
use crate::arr::F2;
use crate::linalg::{DMat, add3, mat3_vec, mat3t_vec, norm3, sub3};
use crate::models::{azel, ecef_to_llh, enu_matrix, geometry, saastamoinen_zenith};
use crate::nputil::median;
use crate::orbits::Orbits;
use crate::solve::receiver_cal;
use crate::tropo::{TropoInputs, troposphere};

/// Most frequencies kept per satellite.
pub const MAX_BANDS: usize = 3;

/// One (satellite, frequency) of the plan.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KeyInfo {
    /// The satellite.
    pub sat: Sat,
    /// Index of the satellite in the plan.
    pub slot: usize,
    /// RINEX band number.
    pub band: u8,
    /// Group index.
    pub group: usize,
    /// Wavelength, m.
    pub lam: f64,
    /// Whether the ambiguity may be fixed.
    pub fixable: bool,
}

/// What stays the same over a session: the signal plan, the base antenna, the calibrations.
#[derive(Debug)]
pub struct Setup<'a> {
    /// The signal plan.
    pub plan: Plan,
    /// Per ambiguity key: satellite, band, group, wavelength.
    pub keys: Vec<KeyInfo>,
    /// The base antenna reference point, ECEF m.
    pub base_arp: [f64; 3],
    /// Base latitude, longitude (rad) and height (m).
    pub base_llh: (f64, f64, f64),
    /// The base's east / north / up rows.
    pub base_enu: [[f64; 3]; 3],
    /// Base antenna calibration.
    pub base_cal: Option<AntennaCal>,
    /// Rover antenna calibration.
    pub rover_cal: Option<AntennaCal>,
    /// Rover antenna reference point above its mark: height, east, north (m).
    pub rover_delta_hen: [f64; 3],
    /// The usual interval of the rover's epochs, s.
    pub interval: f64,
    /// The longest break in a satellite's data that its ambiguity may survive, s: `max_outage_s`, or two
    /// and a half intervals when that is longer (so that low-rate data are not all breaks).
    pub outage_s: f64,
    /// Troposphere grids.
    pub tropo: &'a TropoInputs<'a>,
}

/// Slant troposphere delays (m) at a site for elevations `el` (rad): the engine's troposphere model plus,
/// where the model gives no wet delay, the standard atmosphere's with relative humidity `humidity`.
#[must_use]
pub fn slant_tropo(
    model: &str,
    humidity: f64,
    inp: &TropoInputs,
    t: f64,
    llh: (f64, f64, f64),
    el: &[f64],
) -> Vec<f64> {
    let mut e = F2::zeros(1, el.len());
    for (j, v) in el.iter().enumerate() {
        e.set(0, j, *v);
    }
    let tr = troposphere(model, llh.0, llh.1, llh.2, &[t], &e, inp);
    let zwd = if tr.zwd0 > 0.0 {
        tr.zwd0
    } else {
        saastamoinen_zenith(llh.0, llh.2, humidity).1
    };
    (0..el.len())
        .map(|j| tr.zhd[0] * tr.mh.at(0, j) + zwd * tr.mw.at(0, j))
        .collect()
}

impl<'a> Setup<'a> {
    /// Builds the plan, places the base antenna and looks the antennas up.
    pub fn new(
        input: &mut PpkInput<'a>,
        opts: &PpkOptions,
        notes: &mut Vec<String>,
        warnings: &mut Vec<String>,
    ) -> Result<Self, PpkError> {
        let (base, rover) = (input.base, input.rover);
        if rover.t.is_empty() || base.t.is_empty() {
            return Err(PpkError::Input("no observations".into()));
        }
        let m = input.base_station.marker_xyz;
        if !(norm3(m) > 6.0e6 && norm3(m) < 7.0e6) {
            return Err(PpkError::Input(
                "the base coordinate is not a position on the Earth (ECEF metres expected)".into(),
            ));
        }
        let mut o = opts.clone();
        o.max_bands = o.max_bands.clamp(1, MAX_BANDS);
        let plan = Plan::build(rover, base, input.orbits, input.fcn, &o);
        if plan.sats.len() < 4 {
            return Err(PpkError::Input(format!(
                "the receivers share {} satellites with code and phase on a common frequency; at least 4 \
                 are needed",
                plan.sats.len()
            )));
        }
        let mut why: BTreeMap<&str, Vec<String>> = BTreeMap::new();
        for (s, w) in &plan.left_out {
            why.entry(w.as_str()).or_default().push(s.to_string());
        }
        for (w, sats) in why {
            notes.push(format!("left out ({w}): {}", sats.join(" ")));
        }
        let (lat, lon, _) = ecef_to_llh(m);
        let r = enu_matrix(lat, lon);
        let hen = input.base_station.delta_hen.unwrap_or(base.delta_hen);
        // ENU = (E, N, H)
        let base_arp = add3(m, mat3t_vec(&r, [hen[1], hen[2], hen[0]]));
        let base_llh = ecef_to_llh(base_arp);
        notes.push(format!(
            "base antenna reference point {:.4} m above the mark (east {:.4}, north {:.4})",
            hen[0], hen[1], hen[2]
        ));
        let mut cal_notes = Vec::new();
        let base_cal = receiver_cal(input.antex.as_deref_mut(), base, &mut cal_notes);
        let rover_cal = if rover.antenna.trim().is_empty() {
            None
        } else {
            receiver_cal(input.antex.as_deref_mut(), rover, &mut cal_notes)
        };
        notes.append(&mut cal_notes);
        if base_cal.is_none() {
            warnings.push(if base.antenna.trim().is_empty() {
                "the base file names no antenna type: no phase-centre correction at the base; heights \
                 can be off by the antenna's phase-centre offset (often 5 to 15 cm)"
                    .to_owned()
            } else if input.antex.is_none() {
                format!(
                    "no antenna calibration file given: the base antenna ({}) is used without its \
                     phase-centre correction; heights can be off by its offset (often 5 to 15 cm)",
                    base.antenna.trim()
                )
            } else {
                format!(
                    "the base antenna ({}) is not in the calibration file: no phase-centre correction; \
                     heights can be off by its offset (often 5 to 15 cm)",
                    base.antenna.trim()
                )
            });
        }
        if rover_cal.is_none() {
            notes.push(
                "rover antenna without a calibration: its positions are those of the antenna's phase \
                 centre"
                    .to_owned(),
            );
        }
        let interval = super::events::median_interval(&rover.t);
        let mut keys = vec![
            KeyInfo {
                sat: plan.sats[0].sat,
                slot: 0,
                band: 0,
                group: 0,
                lam: 0.0,
                fixable: false,
            };
            plan.n_keys
        ];
        for (slot, s) in plan.sats.iter().enumerate() {
            for b in &s.bands {
                keys[b.key] = KeyInfo {
                    sat: s.sat,
                    slot,
                    band: b.band,
                    group: b.group,
                    lam: b.lam,
                    fixable: plan.groups[b.group].fixable,
                };
            }
        }
        if plan.groups.iter().any(|g| !g.fixable) {
            notes.push(
                "GLONASS ambiguities are kept float (its satellites still strengthen the solution)"
                    .to_owned(),
            );
        }
        Ok(Setup {
            plan,
            keys,
            base_arp,
            base_llh,
            base_enu: enu_matrix(base_llh.0, base_llh.1),
            base_cal,
            rover_cal,
            rover_delta_hen: rover.delta_hen,
            interval,
            outage_s: opts.max_outage_s.max(2.5 * interval),
            tropo: input.tropo,
        })
    }
}

/// One frequency of one satellite at one epoch.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PrepBand {
    /// Ambiguity key.
    pub key: u16,
    /// A cycle slip (or a new arc) between the previous rover epoch and this one.
    pub slip: bool,
    /// The receiver flags a possible half-cycle ambiguity: not a candidate for fixing.
    pub half: bool,
    /// Rover phase (m) plus the satellite clock, less the base's phase correction: what is left is the
    /// rover's range, troposphere and antenna terms, the receivers' clock difference and the ambiguity.
    pub phase: f64,
    /// The same for the pseudorange (NaN when either code is missing).
    pub code: f64,
}

/// One satellite at one epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrepSat {
    /// Index of the satellite in the plan.
    pub slot: u16,
    /// Number of frequencies in `bands`.
    pub nb: u8,
    /// Satellite position at the transmit time of the rover's signal, ECEF m.
    pub xyz: [f64; 3],
    /// Elevation at the base, rad.
    pub el_base: f32,
    /// The frequencies.
    pub band: [PrepBand; MAX_BANDS],
}

impl PrepSat {
    /// The frequencies observed.
    #[must_use]
    pub fn bands(&self) -> &[PrepBand] {
        &self.band[..usize::from(self.nb)]
    }
}

/// Everything the filter needs at one epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedEpoch {
    /// Index of the rover epoch.
    pub index: usize,
    /// The epoch tag, GPS seconds.
    pub t_tag: f64,
    /// The receive time: the tag less the rover's clock offset.
    pub t: f64,
    /// Code-only position of the rover (ECEF m), when there is one.
    pub spp: Option<[f64; 3]>,
    /// Its covariance for unit range variance (ECEF), when there is a position.
    pub spp_q: [[f64; 3]; 3],
    /// Satellites in the code-only position.
    pub spp_n: u16,
    /// RMS of its residuals, m.
    pub spp_rms: f32,
    /// Seconds to the nearest base epoch used (NaN: no base data for this epoch).
    pub age: f32,
    /// Doppler velocity of the rover, ECEF m/s (NaN without).
    pub vel: [f32; 3],
    /// The satellites with single differences.
    pub sats: Vec<PrepSat>,
}

/// A satellite's positions and clocks at the transmit times of one receiver's epochs.
#[derive(Debug, Clone, Default)]
struct SatTrack {
    xyz: Vec<[f64; 3]>,
    /// Satellite clock, m.
    clk: Vec<f64>,
    ok: Vec<bool>,
    /// ECEF velocity, m/s, and clock drift, m/s (rover only, when Doppler is used).
    vel: Vec<[f64; 4]>,
}

fn sat_track(orb: &Orbits, sat: Sat, t: &[f64], pr: &[f64], with_vel: bool) -> SatTrack {
    let tr: Vec<f64> = t
        .iter()
        .zip(pr)
        .map(|(&x, &p)| x - if p.is_finite() { p } else { 0.07 * CLIGHT } / CLIGHT)
        .collect();
    let (_, clk1, o1, _) = orb.position(sat, &tr);
    let tr2: Vec<f64> = tr.iter().zip(&clk1).map(|(a, b)| a - b).collect();
    let (xyz, clk, o2, _) = orb.position(sat, &tr2);
    let ok: Vec<bool> = o1
        .iter()
        .zip(&o2)
        .zip(pr)
        .map(|((a, b), p)| *a && *b && p.is_finite())
        .collect();
    let mut vel = Vec::new();
    if with_vel {
        const H: f64 = 0.5;
        let ta: Vec<f64> = tr2.iter().map(|x| x - H).collect();
        let tb: Vec<f64> = tr2.iter().map(|x| x + H).collect();
        let (xa, ca, _, _) = orb.position(sat, &ta);
        let (xb, cb, _, _) = orb.position(sat, &tb);
        vel = (0..t.len())
            .map(|i| {
                [
                    (xb[i][0] - xa[i][0]) / (2.0 * H),
                    (xb[i][1] - xa[i][1]) / (2.0 * H),
                    (xb[i][2] - xa[i][2]) / (2.0 * H),
                    CLIGHT * (cb[i] - ca[i]) / (2.0 * H),
                ]
            })
            .collect();
    }
    SatTrack {
        xyz,
        clk: clk.iter().map(|c| c * CLIGHT).collect(),
        ok,
        vel,
    }
}

fn series(st: &StationObs, sat: Sat, code: geomaticist_gnss_core::ObsCode) -> Option<&[f64]> {
    st.sats
        .get(&sat)
        .and_then(|s| s.get(&code))
        .map(|s| s.values.as_slice())
}

fn lli_series(st: &StationObs, sat: Sat, code: geomaticist_gnss_core::ObsCode) -> Option<&[i8]> {
    st.sats
        .get(&sat)
        .and_then(|s| s.get(&code))
        .and_then(|s| s.lli.as_deref())
}

fn sane_code(p: f64) -> bool {
    p.is_finite() && p > 1.0e7 && p < 6.0e7
}

fn sane_phase(l: f64) -> bool {
    l.is_finite() && l.abs() > 1.0e3
}

/// A weighted least-squares code position: unknowns are the position and one clock per system present.
/// `rows`: (system index, satellite position, pseudorange + satellite clock in m). Returns the position,
/// the clocks (m, NaN for absent systems), the position cofactor, the number of satellites used and the
/// RMS of the residuals.
#[allow(clippy::type_complexity)]
pub fn code_position(
    rows: &[(usize, [f64; 3], f64)],
    x0: [f64; 3],
    elmask: f64,
    gross: f64,
) -> Option<([f64; 3], [f64; 5], [[f64; 3]; 3], usize, f64)> {
    let mut use_: Vec<bool> = vec![true; rows.len()];
    let mut x = x0;
    let mut clk = [0.0f64; 5];
    for _reject in 0..4 {
        let mut last: Option<(Vec<f64>, [[f64; 3]; 3], usize)> = None;
        let mut started = norm3(x) > 6.0e6;
        for it in 0..10 {
            let (lat, lon, h) = ecef_to_llh(x);
            let r = enu_matrix(lat, lon);
            let on_earth = norm3(x) > 6.0e6 && h.abs() < 50_000.0;
            let mut sys_on = [false; 5];
            let mut hr: Vec<([f64; 3], usize, f64, f64)> = Vec::with_capacity(rows.len());
            let mut idx: Vec<usize> = Vec::with_capacity(rows.len());
            for (k, (sys, sx, p)) in rows.iter().enumerate() {
                if !use_[k] {
                    continue;
                }
                let (rho, e) = geometry(*sx, x);
                let (mut w, mut trop) = (1.0, 0.0);
                if on_earth && started {
                    let (_, el, _) = azel(e, &r);
                    if el < elmask {
                        continue;
                    }
                    let s = el.max(0.05).sin();
                    w = s * s;
                    let (zh, zw) = saastamoinen_zenith(lat, h, 0.5);
                    trop = (zh + zw) / s;
                }
                sys_on[*sys] = true;
                hr.push((e, *sys, p - rho - trop - clk[*sys], w));
                idx.push(k);
            }
            let cols: Vec<usize> = (0..5).filter(|s| sys_on[*s]).collect();
            let nu = 3 + cols.len();
            if hr.len() < nu + usize::from(it > 0) {
                return None;
            }
            let mut n = DMat::zeros(nu);
            let mut b = vec![0.0; nu];
            for (e, sys, v, w) in &hr {
                let mut a = vec![0.0; nu];
                a[0] = -e[0];
                a[1] = -e[1];
                a[2] = -e[2];
                let c = cols.iter().position(|s| s == sys).unwrap_or(0);
                a[3 + c] = 1.0;
                for i in 0..nu {
                    b[i] += w * a[i] * v;
                    for j in 0..nu {
                        n.add(i, j, w * a[i] * a[j]);
                    }
                }
            }
            let ninv = n.inverse()?;
            let dx = ninv.matvec(&b);
            x = [x[0] + dx[0], x[1] + dx[1], x[2] + dx[2]];
            for (c, s) in cols.iter().enumerate() {
                clk[*s] += dx[3 + c];
            }
            let step = norm3([dx[0], dx[1], dx[2]]);
            if !started {
                started = true;
                continue;
            }
            if step < 1e-3 {
                // residuals after the step
                let res: Vec<f64> = hr
                    .iter()
                    .map(|(e, sys, v, _)| {
                        let c = cols.iter().position(|s| s == sys).unwrap_or(0);
                        v + e[0] * dx[0] + e[1] * dx[1] + e[2] * dx[2] - dx[3 + c]
                    })
                    .collect();
                let q = [
                    [ninv.at(0, 0), ninv.at(0, 1), ninv.at(0, 2)],
                    [ninv.at(1, 0), ninv.at(1, 1), ninv.at(1, 2)],
                    [ninv.at(2, 0), ninv.at(2, 1), ninv.at(2, 2)],
                ];
                let mut full = vec![f64::NAN; rows.len()];
                for (r_, k) in res.iter().zip(&idx) {
                    full[*k] = *r_;
                }
                last = Some((full, q, nu));
                break;
            }
        }
        let (res, q, nu) = last?;
        let used: Vec<usize> = (0..rows.len()).filter(|k| res[*k].is_finite()).collect();
        let n = used.len();
        let rms = (used.iter().map(|k| res[*k] * res[*k]).sum::<f64>() / n.max(1) as f64).sqrt();
        let worst = used
            .iter()
            .copied()
            .max_by(|a, b| res[*a].abs().total_cmp(&res[*b].abs()));
        match worst {
            Some(k) if res[k].abs() > gross && n > nu + 1 => use_[k] = false,
            _ => {
                let mut c = [f64::NAN; 5];
                for s in 0..5 {
                    if rows
                        .iter()
                        .zip(&res)
                        .any(|(r_, v)| r_.0 == s && v.is_finite())
                    {
                        c[s] = clk[s];
                    }
                }
                return Some((x, c, q, n, rms));
            }
        }
    }
    None
}

/// Index of a system among the clock unknowns of [`code_position`].
#[must_use]
pub fn sys_index(sys: GnssSystem) -> usize {
    match sys {
        GnssSystem::Glonass => 1,
        GnssSystem::Galileo => 2,
        GnssSystem::BeiDou => 3,
        GnssSystem::Qzss | GnssSystem::Gps | GnssSystem::Navic | GnssSystem::Sbas => 0,
    }
}

/// The base's corrections at one base epoch: per key, phase and code (m; NaN when missing).
#[derive(Debug, Clone, Default)]
struct BaseEpoch {
    phase: Vec<f64>,
    code: Vec<f64>,
    el: Vec<f32>,
}

struct KeyData<'s> {
    rl: &'s [f64],
    rlli: Option<&'s [i8]>,
    rc: &'s [f64],
    rd: Option<&'s [f64]>,
    bl: &'s [f64],
    bc: &'s [f64],
    /// Cumulative count of base loss-of-lock flags (index j + 1: up to and including base epoch j).
    b_lli_cum: Vec<u32>,
}

fn base_epoch(
    setup: &Setup,
    opts: &PpkOptions,
    t: f64,
    j: usize,
    tracks: &[SatTrack],
    kd: &[KeyData],
) -> BaseEpoch {
    let nk = setup.keys.len();
    let ns = setup.plan.sats.len();
    let mut out = BaseEpoch {
        phase: vec![f64::NAN; nk],
        code: vec![f64::NAN; nk],
        el: vec![f32::NAN; ns],
    };
    let mut az = vec![0.0; ns];
    let mut el = vec![0.0; ns];
    let mut rho = vec![f64::NAN; ns];
    for (s, tr) in tracks.iter().enumerate() {
        if !tr.ok[j] {
            continue;
        }
        let (r, e) = geometry(tr.xyz[j], setup.base_arp);
        let (a, l, _) = azel(e, &setup.base_enu);
        az[s] = a;
        el[s] = l;
        rho[s] = r;
        out.el[s] = l as f32;
    }
    let trop = slant_tropo(
        &opts.trop_model,
        opts.humidity,
        setup.tropo,
        t,
        setup.base_llh,
        &el,
    );
    for (s, sp) in setup.plan.sats.iter().enumerate() {
        if !rho[s].is_finite() || el[s] < 3f64.to_radians() {
            continue;
        }
        for b in &sp.bands {
            let k = &kd[b.key];
            let ant = receiver_correction(
                setup.base_cal.as_ref(),
                sp.sat.sys.letter() as u8,
                b.band,
                &az[s..=s],
                &el[s..=s],
                opts.pcv_azimuth,
            )[0];
            let model = rho[s] + trop[s] + ant - tracks[s].clk[j];
            if sane_phase(k.bl[j]) {
                out.phase[b.key] = k.bl[j] * b.lam - model;
            }
            if sane_code(k.bc[j]) {
                out.code[b.key] = k.bc[j] - model;
            }
        }
    }
    out
}

/// A phase is taken as whole (cleared, or repaired by whole cycles) when its jump against the motion of
/// the other phases is this close to a whole number of cycles.
pub const REPAIR_MAX_FRAC: f64 = 0.15;
/// Fewest agreeing phases behind the motion before another phase is cleared or repaired against it.
pub const REPAIR_MIN_ROWS: usize = 8;

/// The memory of the slip checks.
#[derive(Debug, Clone)]
struct Memory {
    /// Rover tag of the previous prepared epoch.
    t_prev: f64,
    /// The rover's code-only position at the previous epoch.
    x_prev: Option<[f64; 3]>,
    /// Per key: the tag at which it was last seen (NaN: never).
    t_last: Vec<f64>,
    /// Per key, as last seen: the single-difference phase less the modelled range.
    td_last: Vec<f64>,
    /// Per key, as last seen: the rover's phase in metres (for the geometry-free check).
    rphase_last: Vec<f64>,
    /// Per key, as last seen: the base's phase correction.
    bphase_last: Vec<f64>,
    /// Per key, as last seen: the half-cycle flag.
    half_last: Vec<bool>,
    /// Per key: the motion summed up to when it was last seen, and the chain that sum belongs to.
    cum_last: Vec<[f64; 4]>,
    chain_last: Vec<u32>,
    /// Per key: the latest base epoch used.
    base_j: Vec<usize>,
    /// Per key: whole cycles taken out of the rover phase by repairs so far.
    fix_cycles: Vec<f64>,
    /// The epoch-to-epoch motion (three position components and a clock, as the time-difference check
    /// solves them) summed since the chain began; a step without a solution starts a new chain.
    cum: [f64; 4],
    chain: u32,
}

/// What the preparation did beside flagging slips.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PrepStats {
    /// Cycle slips of whole cycles taken out of the phases (the ambiguity carried on).
    pub repaired: usize,
}

/// Prepares every rover epoch.
#[allow(clippy::too_many_lines)]
pub fn prepare_all(
    setup: &Setup,
    input: &PpkInput,
    opts: &PpkOptions,
    progress: &mut dyn FnMut(Progress) -> bool,
    warnings: &mut Vec<String>,
) -> Result<(Vec<PreparedEpoch>, PrepStats), PpkError> {
    let (rover, base, orb) = (input.rover, input.base, input.orbits);
    let (nr, nbe) = (rover.t.len(), base.t.len());
    let nk = setup.keys.len();
    let ns = setup.plan.sats.len();
    let empty_f: Vec<f64> = Vec::new();

    // the series of every key
    let mut kd: Vec<KeyData> = Vec::with_capacity(nk);
    for sp in &setup.plan.sats {
        for b in &sp.bands {
            let bl_lli = lli_series(base, sp.sat, b.base_l);
            let mut cum = vec![0u32; nbe + 1];
            for j in 0..nbe {
                let l = bl_lli.is_some_and(|v| v[j] & 1 != 0);
                cum[j + 1] = cum[j] + u32::from(l);
            }
            kd.push(KeyData {
                rl: series(rover, sp.sat, b.rover_l).unwrap_or(&empty_f),
                rlli: lli_series(rover, sp.sat, b.rover_l),
                rc: series(rover, sp.sat, b.rover_c).unwrap_or(&empty_f),
                rd: b.rover_d.and_then(|c| series(rover, sp.sat, c)),
                bl: series(base, sp.sat, b.base_l).unwrap_or(&empty_f),
                bc: series(base, sp.sat, b.base_c).unwrap_or(&empty_f),
                b_lli_cum: cum,
            });
        }
    }
    if kd
        .iter()
        .any(|k| k.rl.len() != nr || k.rc.len() != nr || k.bl.len() != nbe || k.bc.len() != nbe)
    {
        return Err(PpkError::Input(
            "observation series of unequal length".into(),
        ));
    }
    let with_vel = opts.velocity && kd.iter().any(|k| k.rd.is_some());
    // The sign of the Doppler. RINEX defines it as positive when the range decreases, so that it is minus
    // the rate of the phase; some receivers write the opposite. The file's own phases decide.
    let mut doppler_sign = 1.0;
    if with_vel {
        let (mut std_votes, mut inv_votes) = (0usize, 0usize);
        let step = (nr / 400).max(1);
        for k in kd.iter().filter(|k| k.rd.is_some()) {
            let Some(d) = k.rd else { continue };
            for i in (0..nr.saturating_sub(1)).step_by(step) {
                let dt = rover.t[i + 1] - rover.t[i];
                let (l0, l1, d0, d1) = (k.rl[i], k.rl[i + 1], d[i], d[i + 1]);
                if !(dt > 0.0 && dt < 2.5 && sane_phase(l0) && sane_phase(l1)) {
                    continue;
                }
                if !(d0.is_finite() && d1.is_finite()) || k.rlli.is_some_and(|v| v[i + 1] & 1 != 0)
                {
                    continue;
                }
                let (rate, mean) = ((l1 - l0) / dt, 0.5 * (d0 + d1));
                if rate.abs() < 50.0 {
                    continue;
                }
                if (rate + mean).abs() < (rate - mean).abs() {
                    std_votes += 1;
                } else {
                    inv_votes += 1;
                }
            }
        }
        if inv_votes > std_votes {
            doppler_sign = -1.0;
            warnings.push(
                "the rover file's Doppler has the opposite sign of the RINEX definition (positive when the \
                 range grows): its sign was turned for the velocities"
                    .to_owned(),
            );
        }
    }

    // satellite tracks at both receivers' transmit times (the first frequency with a code gives the time)
    let first_code =
        |kd: &[KeyData], sp: &super::plan::SatPlan, i: usize, rover_side: bool| -> f64 {
            for b in &sp.bands {
                let p = if rover_side {
                    kd[b.key].rc[i]
                } else {
                    kd[b.key].bc[i]
                };
                if sane_code(p) {
                    return p;
                }
            }
            f64::NAN
        };
    let mut rtracks: Vec<SatTrack> = Vec::with_capacity(ns);
    let mut btracks: Vec<SatTrack> = Vec::with_capacity(ns);
    for (s, sp) in setup.plan.sats.iter().enumerate() {
        let pr: Vec<f64> = (0..nr).map(|i| first_code(&kd, sp, i, true)).collect();
        rtracks.push(sat_track(orb, sp.sat, &rover.t, &pr, with_vel));
        let pb: Vec<f64> = (0..nbe).map(|j| first_code(&kd, sp, j, false)).collect();
        btracks.push(sat_track(orb, sp.sat, &base.t, &pb, false));
        if s % 8 == 0
            && !progress(Progress {
                stage: Stage::Prepare,
                done: s * nr / (2 * ns.max(1)),
                total: nr,
            })
        {
            return Err(PpkError::Cancelled);
        }
    }

    let mut out: Vec<PreparedEpoch> = Vec::with_capacity(nr);
    let mut mem = Memory {
        t_prev: f64::NEG_INFINITY,
        x_prev: None,
        t_last: vec![f64::NAN; nk],
        td_last: vec![f64::NAN; nk],
        rphase_last: vec![f64::NAN; nk],
        bphase_last: vec![f64::NAN; nk],
        half_last: vec![false; nk],
        cum_last: vec![[0.0; 4]; nk],
        chain_last: vec![0; nk],
        base_j: vec![usize::MAX; nk],
        fix_cycles: vec![0.0; nk],
        cum: [0.0; 4],
        chain: 1,
    };
    let mut n_repaired = 0usize;
    // base corrections of the two epochs around the rover epoch, cached
    let mut cache: [(usize, BaseEpoch); 2] = [
        (usize::MAX, BaseEpoch::default()),
        (usize::MAX, BaseEpoch::default()),
    ];
    let mut jb = 0usize;
    let mut x_last = setup.base_arp;
    let (mut n_no_base, mut n_td_all) = (0usize, 0usize);

    for i in 0..nr {
        let t = rover.t[i];
        // --- code-only position and clock of the rover
        let mut rows: Vec<(usize, [f64; 3], f64)> = Vec::with_capacity(ns);
        for (s, sp) in setup.plan.sats.iter().enumerate() {
            if !rtracks[s].ok[i] {
                continue;
            }
            let p = first_code(&kd, sp, i, true);
            rows.push((
                sys_index(sp.sat.sys),
                rtracks[s].xyz[i],
                p + rtracks[s].clk[i],
            ));
        }
        let spp = code_position(
            &rows,
            x_last,
            opts.elmask.to_radians() * 0.5,
            opts.code_gross,
        );
        let mut pe = PreparedEpoch {
            index: i,
            t_tag: t,
            t,
            spp: None,
            spp_q: [[0.0; 3]; 3],
            spp_n: 0,
            spp_rms: 0.0,
            age: f32::NAN,
            vel: [f32::NAN; 3],
            sats: Vec::new(),
        };
        if rows.is_empty() {
            // an epoch line with no observations at all (some tools write event marks that way): it is
            // not a break in tracking, so the slip memory is left alone
            out.push(pe);
            continue;
        }
        if let Some((x, clk, q, n, rms)) = spp {
            // a position far from the base (wrong by kilometres) is no position
            if norm3(sub3(x, setup.base_arp)) < 2.0e6 {
                pe.spp = Some(x);
                pe.spp_q = q;
                pe.spp_n = n as u16;
                pe.spp_rms = rms as f32;
                let c = clk.iter().copied().find(|c| c.is_finite()).unwrap_or(0.0);
                pe.t = t - c / CLIGHT;
                x_last = x;
            }
        }
        let x_model = pe.spp.unwrap_or(x_last);

        // --- Doppler velocity
        if with_vel && pe.spp.is_some() {
            pe.vel = doppler_velocity(setup, &kd, &rtracks, i, x_model, doppler_sign);
        }

        // --- the base epochs around this one
        while jb < nbe && base.t[jb] < t - opts.base_tol {
            jb += 1;
        }
        // (j0, j1, weight of j1): j1 alone when it is at the rover's time
        let pick: Option<(usize, usize, f64)> =
            if jb < nbe && (base.t[jb] - t).abs() <= opts.base_tol {
                Some((jb, jb, 1.0))
            } else if jb > 0 && jb < nbe && base.t[jb] - base.t[jb - 1] <= opts.base_max_gap {
                let w = (t - base.t[jb - 1]) / (base.t[jb] - base.t[jb - 1]);
                Some((jb - 1, jb, w))
            } else {
                None
            };
        let Some((j0, j1, w1)) = pick else {
            n_no_base += 1;
            // no base: nothing can be carried across this epoch
            mem.chain += 1;
            mem.cum = [0.0; 4];
            mem.t_last.fill(f64::NAN);
            mem.t_prev = t;
            mem.x_prev = pe.spp;
            out.push(pe);
            continue;
        };
        pe.age = (t - base.t[j0]).abs().min((t - base.t[j1]).abs()) as f32;
        for j in [j0, j1] {
            if cache[0].0 != j && cache[1].0 != j {
                // replace the entry this epoch does not need
                let slot = usize::from(cache[0].0 == j0 || cache[0].0 == j1);
                cache[slot] = (j, base_epoch(setup, opts, base.t[j], j, &btracks, &kd));
            }
        }
        let be = |j: usize| -> &BaseEpoch {
            if cache[0].0 == j {
                &cache[0].1
            } else {
                &cache[1].1
            }
        };
        let (b0, b1) = (be(j0), be(j1));

        // --- rover geometry at the code position, for the slip check
        let (rlat, rlon, rh) = ecef_to_llh(x_model);
        let renu = enu_matrix(rlat, rlon);
        let mut rel = vec![0.0; ns];
        let mut rrho = vec![f64::NAN; ns];
        let mut los = vec![[0.0; 3]; ns];
        for s in 0..ns {
            if rtracks[s].ok[i] {
                let (r, e) = geometry(rtracks[s].xyz[i], x_model);
                let (_, l, _) = azel(e, &renu);
                rel[s] = l;
                rrho[s] = r;
                los[s] = e;
            }
        }
        let rtrop = slant_tropo(
            &opts.trop_model,
            opts.humidity,
            setup.tropo,
            t,
            (rlat, rlon, rh),
            &rel,
        );

        let dt = t - mem.t_prev;
        let stepped = dt > 0.0 && dt <= setup.outage_s;
        let mask = opts.elmask.to_radians();
        // per key at this epoch: single-difference phase less the modelled range, rover phase (m), base
        // correction, half-cycle flag
        let mut td_now = vec![f64::NAN; nk];
        let mut rphase_now = vec![f64::NAN; nk];
        let mut bphase_now = vec![f64::NAN; nk];
        let mut half_now = vec![false; nk];
        // rows of the time-difference check, phases seen at the previous epoch and not flagged:
        // (key, slot, misfit, noise scale)
        let mut td_rows: Vec<(usize, usize, f64, f64)> = Vec::new();
        // phases to test against the motion of the others: flagged by a receiver, or back after a short
        // break: (key, slot)
        let mut td_cands: Vec<(usize, usize)> = Vec::new();
        for (s, sp) in setup.plan.sats.iter().enumerate() {
            if !rrho[s].is_finite() || !b0.el[s].is_finite() || !b1.el[s].is_finite() {
                continue;
            }
            let mut ps = PrepSat {
                slot: s as u16,
                nb: 0,
                xyz: rtracks[s].xyz[i],
                el_base: b1.el[s],
                band: [PrepBand::default(); MAX_BANDS],
            };
            for b in &sp.bands {
                let k = &kd[b.key];
                let l = k.rl[i];
                if !sane_phase(l) {
                    continue;
                }
                let (c0, c1) = (b0.phase[b.key], b1.phase[b.key]);
                if !(c0.is_finite() && c1.is_finite()) {
                    continue;
                }
                // a base slip between the two epochs breaks the interpolation
                if j0 != j1 && k.b_lli_cum[j1 + 1] > k.b_lli_cum[j0 + 1] {
                    continue;
                }
                let bph = (1.0 - w1) * c0 + w1 * c1;
                let clk = rtracks[s].clk[i];
                // the rover phase with the cycles of earlier repairs taken out
                let rph = (l - mem.fix_cycles[b.key]) * b.lam;
                let phase = rph + clk - bph;
                let (p, pc0, pc1) = (k.rc[i], b0.code[b.key], b1.code[b.key]);
                let code = if sane_code(p) && pc0.is_finite() && pc1.is_finite() {
                    p + clk - ((1.0 - w1) * pc0 + w1 * pc1)
                } else {
                    f64::NAN
                };
                let lli = k.rlli.map_or(0, |v| v[i]);
                let half = lli & 2 != 0;
                // A loss-of-lock flag (rover now, or base since the base epoch last used for this key)
                // makes the phase suspect, and so does a short break in it. A suspect phase is tested
                // below against the motion the other phases agree on; it starts a new arc only when it
                // fails that test or cannot be tested. (Receivers keep the flag up for several epochs
                // after one loss, on phases that are whole.)
                let jl = mem.base_j[b.key];
                let first_use = jl == usize::MAX || jl > j1;
                let flagged =
                    lli & 1 != 0 || (!first_use && k.b_lli_cum[j1 + 1] > k.b_lli_cum[jl + 1]);
                let since = t - mem.t_last[b.key];
                let seen_prev = stepped && mem.t_last[b.key] == mem.t_prev;
                let bridgeable = !seen_prev
                    && since > 0.0
                    && since <= setup.outage_s
                    && mem.chain_last[b.key] == mem.chain;
                let same_half = half == mem.half_last[b.key];
                let td = phase - rrho[s] - rtrop[s];
                let mut slip = true;
                if !first_use && same_half && mem.td_last[b.key].is_finite() {
                    if seen_prev && !flagged {
                        slip = false;
                        // below the mask the filter does not use the phase; low satellites are noisier
                        if rel[s] >= mask {
                            let scale = (0.5 / rel[s].sin().max(0.05)).max(1.0);
                            td_rows.push((b.key, s, td - mem.td_last[b.key], scale));
                        }
                    } else if (seen_prev || bridgeable) && rel[s] >= mask {
                        td_cands.push((b.key, s));
                    }
                }
                td_now[b.key] = td;
                rphase_now[b.key] = rph;
                bphase_now[b.key] = bph;
                half_now[b.key] = half;
                mem.base_j[b.key] = j1;
                ps.band[usize::from(ps.nb)] = PrepBand {
                    key: b.key as u16,
                    slip,
                    half,
                    phase,
                    code,
                };
                ps.nb += 1;
            }
            if ps.nb > 0 {
                pe.sats.push(ps);
            }
        }

        // --- the time-difference check: 3 position components + 1 clock from all unflagged phases
        let mut slipped: Vec<usize> = Vec::new();
        let mut cleared: Vec<(usize, f64)> = Vec::new(); // (key, whole cycles taken out)
        let mut motion: Option<[f64; 4]> = None;
        if stepped && pe.spp.is_some() && mem.x_prev.is_some() {
            let (bad, all, m) = td_check(&td_rows, &los, opts.slip_td);
            let n_clean = td_rows.len() - bad.len().min(td_rows.len());
            slipped.extend(bad);
            if all {
                n_td_all += 1;
            }
            motion = m;
            // the suspect phases against the motion of the others since they were last seen: no jump
            // clears them; a jump of whole cycles is taken out (repaired); anything else is a new arc
            if let Some(x) = m
                && n_clean >= REPAIR_MIN_ROWS
            {
                for &(key, s) in &td_cands {
                    let e = los[s];
                    let c = mem.cum_last[key];
                    let cum = [
                        mem.cum[0] + x[0] - c[0],
                        mem.cum[1] + x[1] - c[1],
                        mem.cum[2] + x[2] - c[2],
                        mem.cum[3] + x[3] - c[3],
                    ];
                    let pred = -e[0] * cum[0] - e[1] * cum[1] - e[2] * cum[2] + cum[3];
                    let cyc = (td_now[key] - mem.td_last[key] - pred) / setup.keys[key].lam;
                    if (cyc - cyc.round()).abs() <= REPAIR_MAX_FRAC {
                        cleared.push((key, cyc.round()));
                    }
                }
            }
        }
        // apply what was cleared or repaired
        if !cleared.is_empty() {
            for ps in &mut pe.sats {
                for q in 0..usize::from(ps.nb) {
                    let key = usize::from(ps.band[q].key);
                    if let Some(&(_, n)) = cleared.iter().find(|c| c.0 == key) {
                        let shift = n * setup.keys[key].lam;
                        ps.band[q].slip = false;
                        ps.band[q].phase -= shift;
                        td_now[key] -= shift;
                        rphase_now[key] -= shift;
                        mem.fix_cycles[key] += n;
                        if n != 0.0 {
                            n_repaired += 1;
                        }
                    }
                }
            }
        }

        // --- geometry-free jumps of either receiver (two frequencies), on what is taken as continuous
        for (s, sp) in setup.plan.sats.iter().enumerate() {
            if sp.bands.len() < 2 || rel[s] < mask {
                continue;
            }
            let (k1, k2) = (sp.bands[0].key, sp.bands[1].key);
            let since = t - mem.t_last[k1].min(mem.t_last[k2]);
            if !(since > 0.0 && since <= setup.outage_s) || mem.t_last[k1] != mem.t_last[k2] {
                continue;
            }
            // noisier phases low down, as in the time-difference check
            let lim = (opts.slip_gf + 0.003 * since) * (0.5 / rel[s].sin().max(0.05)).max(1.0);
            let gf = |v: &[f64]| v[k1] - v[k2];
            let dr = gf(&rphase_now) - gf(&mem.rphase_last);
            let db = gf(&bphase_now) - gf(&mem.bphase_last);
            if (dr.is_finite() && dr.abs() > lim) || (db.is_finite() && db.abs() > lim) {
                for b in &sp.bands {
                    slipped.push(b.key);
                }
            }
        }
        if !slipped.is_empty() {
            for ps in &mut pe.sats {
                for q in 0..usize::from(ps.nb) {
                    if slipped.contains(&usize::from(ps.band[q].key)) {
                        ps.band[q].slip = true;
                    }
                }
            }
        }

        // --- memory
        match motion {
            Some(x) => {
                for q in 0..4 {
                    mem.cum[q] += x[q];
                }
            }
            None => {
                // no motion for this step: nothing seen before it can be bridged across it
                mem.chain += 1;
                mem.cum = [0.0; 4];
            }
        }
        for key in 0..nk {
            if td_now[key].is_finite() {
                mem.t_last[key] = t;
                mem.td_last[key] = td_now[key];
                mem.rphase_last[key] = rphase_now[key];
                mem.bphase_last[key] = bphase_now[key];
                mem.half_last[key] = half_now[key];
                mem.cum_last[key] = mem.cum;
                mem.chain_last[key] = mem.chain;
            }
        }
        mem.t_prev = t;
        mem.x_prev = pe.spp;
        out.push(pe);
        if i % 1000 == 0
            && !progress(Progress {
                stage: Stage::Prepare,
                done: nr / 2 + i / 2,
                total: nr,
            })
        {
            return Err(PpkError::Cancelled);
        }
    }
    if n_no_base > 0 {
        warnings.push(format!(
            "{n_no_base} of {nr} rover epochs have no base data within {:.1} s on both sides: they get a \
             code-only position at best",
            opts.base_max_gap
        ));
    }
    if n_td_all > 0 {
        warnings.push(format!(
            "at {n_td_all} epochs the phases of most satellites jumped together (loss of lock of the whole \
             receiver): all ambiguities were restarted there"
        ));
    }
    Ok((
        out,
        PrepStats {
            repaired: n_repaired,
        },
    ))
}

/// The time-difference slip check. `rows`: (key, satellite slot, change of the single-difference phase
/// less the modelled range since the previous epoch, noise scale: 1 above 30° of elevation, growing as
/// 1 / sin below). Solves the change of position and clock all rows share and returns the keys that
/// misfit by more than `limit` times their scale, and whether the set as a whole was inconsistent (then
/// every key is returned), and the change of position and clock found (when the rows left agree on one).
fn td_check(
    rows: &[(usize, usize, f64, f64)],
    los: &[[f64; 3]],
    limit: f64,
) -> (Vec<usize>, bool, Option<[f64; 4]>) {
    let mut live: Vec<usize> = (0..rows.len()).collect();
    let mut bad: Vec<usize> = Vec::new();
    let min_rows = 6;
    loop {
        if live.len() < min_rows {
            // too few to tell: with slips already found, nothing left can be trusted
            if bad.is_empty() {
                return (bad, false, None);
            }
            return (rows.iter().map(|r| r.0).collect(), true, None);
        }
        let mut n = DMat::zeros(4);
        let mut b = [0.0; 4];
        for &r in &live {
            let e = los[rows[r].1];
            let a = [-e[0], -e[1], -e[2], 1.0];
            let w = 1.0 / (rows[r].3 * rows[r].3);
            for i in 0..4 {
                b[i] += w * a[i] * rows[r].2;
                for j in 0..4 {
                    n.add(i, j, w * a[i] * a[j]);
                }
            }
        }
        let Some(x) = n.solve(&b) else {
            return (bad, false, None);
        };
        // misfits in units of each row's noise scale
        let res: Vec<f64> = live
            .iter()
            .map(|&r| {
                let e = los[rows[r].1];
                (rows[r].2 - (-e[0] * x[0] - e[1] * x[1] - e[2] * x[2] + x[3])) / rows[r].3
            })
            .collect();
        let (worst, vmax) = res
            .iter()
            .enumerate()
            .map(|(i, v)| (i, v.abs()))
            .fold((0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
        if vmax <= limit {
            return (bad, false, Some([x[0], x[1], x[2], x[3]]));
        }
        // more than a third slipped at once: the whole receiver lost lock, or nothing can be told apart
        if bad.len() * 3 > rows.len() {
            return (rows.iter().map(|r| r.0).collect(), true, None);
        }
        bad.push(rows[live[worst]].0);
        live.remove(worst);
    }
}

/// The rover's velocity from the Doppler of the first frequency of each satellite (ECEF m/s).
fn doppler_velocity(
    setup: &Setup,
    kd: &[KeyData],
    tracks: &[SatTrack],
    i: usize,
    x: [f64; 3],
    sign: f64,
) -> [f32; 3] {
    let mut rows: Vec<([f64; 3], f64)> = Vec::new();
    for (s, sp) in setup.plan.sats.iter().enumerate() {
        if !tracks[s].ok[i] || tracks[s].vel.is_empty() {
            continue;
        }
        let Some(b) = sp.bands.iter().find(|b| kd[b.key].rd.is_some()) else {
            continue;
        };
        let d = kd[b.key].rd.map_or(f64::NAN, |v| v[i]);
        if !d.is_finite() {
            continue;
        }
        let sx = tracks[s].xyz[i];
        let v = tracks[s].vel[i];
        let (_, e) = geometry(sx, x);
        let sag = OMGE * (v[0] * x[1] - v[1] * x[0]) / CLIGHT;
        // range rate = e·(vs - vr) + sag + drift_r - drift_s  =>  -e·vr + drift_r = rate - e·vs - sag + drift_s
        let rate = -sign * b.lam * d;
        rows.push((
            e,
            rate - (e[0] * v[0] + e[1] * v[1] + e[2] * v[2]) - sag + v[3],
        ));
    }
    if rows.len() < 6 {
        return [f32::NAN; 3];
    }
    let mut keep = vec![true; rows.len()];
    for _ in 0..rows.len() / 2 + 1 {
        let mut n = DMat::zeros(4);
        let mut b = [0.0; 4];
        for ((e, y), k) in rows.iter().zip(&keep) {
            if !k {
                continue;
            }
            let a = [-e[0], -e[1], -e[2], 1.0];
            for p in 0..4 {
                b[p] += a[p] * y;
                for q in 0..4 {
                    n.add(p, q, a[p] * a[q]);
                }
            }
        }
        let Some(sol) = n.solve(&b) else {
            return [f32::NAN; 3];
        };
        let res: Vec<f64> = rows
            .iter()
            .map(|(e, y)| y - (-e[0] * sol[0] - e[1] * sol[1] - e[2] * sol[2] + sol[3]))
            .collect();
        let live: Vec<f64> = res
            .iter()
            .zip(&keep)
            .filter(|(_, k)| **k)
            .map(|(v, _)| v.abs())
            .collect();
        let lim = (4.0 * 1.4826 * median(&live)).max(0.2);
        let worst = res
            .iter()
            .zip(&keep)
            .enumerate()
            .filter(|(_, (_, k))| **k)
            .map(|(q, (v, _))| (q, v.abs()))
            .fold((0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
        if worst.1 <= lim {
            return [sol[0] as f32, sol[1] as f32, sol[2] as f32];
        }
        if live.len() <= 6 {
            // what is left does not agree: no velocity rather than a wrong one
            return [f32::NAN; 3];
        }
        keep[worst.0] = false;
    }
    [f32::NAN; 3]
}

/// East / north / up of an ECEF vector at a position.
#[must_use]
pub fn to_enu(at: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    let (lat, lon, _) = ecef_to_llh(at);
    mat3_vec(&enu_matrix(lat, lon), v)
}
