//! Kinematic post-processing (PPK): a moving receiver positioned epoch by epoch against a static base.
//!
//! The pieces, in the order the data flows:
//!
//! - [`options`]: [`PpkOptions`], every setting with its default.
//! - [`plan`]: which signals of each satellite both receivers share, and how they group for differencing.
//! - [`prepare`]: one rover epoch plus the base epochs around it become a [`prepare::PreparedEpoch`]
//!   (satellite positions at the transmit time, the base's corrections brought to the rover's epoch,
//!   cycle-slip flags, a code-only position, a Doppler velocity).
//! - [`filter`]: [`filter::PpkFilter`], the estimator. One call per epoch; float and fixed solutions.
//! - [`combine`]: the forward and backward passes merged into one track.
//! - [`events`]: positions at arbitrary times (camera exposures) by interpolation between epochs.
//! - [`pos`]: RTKLIB-compatible `.pos` text and a CSV for events.
//! - [`synth`]: simulated observations of a receiver on a known trajectory, for tests.
//!
//! [`solve`] runs all of it for two [`StationObs`]. The corrections are the static engine's own
//! ([`crate::models`], [`crate::tropo`], [`crate::orbits`], the antenna calibrations of
//! `geomaticist-gnss-products`); nothing is modelled twice.
//!
//! # Method
//!
//! Observations are single differences between the receivers (rover minus base) on each frequency, not
//! combined between frequencies. The states are the rover position, with no memory from one epoch to the next
//! (no assumption about how the receiver moves), and one single-difference ambiguity per satellite and
//! frequency, constant until a cycle slip. Each epoch the filter forms double differences against a reference
//! satellite per signal group, updates the float states, searches the integers (LAMBDA, the engine's
//! [`crate::ambiguity`]) and, when the checks pass, computes the fixed position. The float states are never
//! overwritten by a fix, so a wrong fix cannot feed back into later epochs.
//!
//! The model assumes the ionosphere is the same at both ends, which holds for short baselines (see
//! [`PpkOptions::max_baseline_km`]).

pub mod combine;
pub mod events;
pub mod filter;
pub mod options;
pub mod plan;
pub mod pos;
pub mod prepare;
pub mod synth;

use std::collections::BTreeMap;

use geomaticist_gnss_core::event::Event;
use geomaticist_gnss_core::{Sat, StationObs};
use geomaticist_gnss_products::Antex;
use serde::{Deserialize, Serialize};

pub use events::{EventInterp, EventResult, EventStatus, Occupation, OccupationResult};
pub use options::{ArMode, Direction, PpkOptions};

use crate::linalg::{norm3, sub3};
use crate::models::{ecef_to_llh, enu_matrix};
use crate::orbits::Orbits;
use crate::tropo::TropoInputs;

/// How an epoch was solved. The numbers are RTKLIB's quality flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum FixStatus {
    /// No position.
    #[default]
    None,
    /// Integer ambiguities fixed and validated (Q = 1).
    Fixed,
    /// Carrier-phase solution with float ambiguities (Q = 2).
    Float,
    /// Code-only position of the rover alone (Q = 5).
    Single,
    /// Precise point positioning: one receiver with precise orbits and clocks, float ambiguities (Q = 6).
    Ppp,
}

impl FixStatus {
    /// The RTKLIB quality flag: 1 fix, 2 float, 5 single, 0 none.
    #[must_use]
    pub fn q(self) -> u8 {
        match self {
            FixStatus::None => 0,
            FixStatus::Fixed => 1,
            FixStatus::Float => 2,
            FixStatus::Single => 5,
            FixStatus::Ppp => 6,
        }
    }

    fn rank(self) -> u8 {
        match self {
            FixStatus::Fixed => 0,
            FixStatus::Float => 1,
            FixStatus::Single => 3,
            FixStatus::Ppp => 2,
            FixStatus::None => 4,
        }
    }

    /// The worse of two statuses (`None` is worst, then `Single`, `Float`, `Fixed`).
    #[must_use]
    pub fn worst(self, other: FixStatus) -> FixStatus {
        if self.rank() >= other.rank() {
            self
        } else {
            other
        }
    }

    /// Whether `self` is a better solution than `other`.
    #[must_use]
    pub fn better_than(self, other: FixStatus) -> bool {
        self.rank() < other.rank()
    }

    /// `"fixed"`, `"float"`, `"single"` or `"none"`.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            FixStatus::None => "none",
            FixStatus::Fixed => "fixed",
            FixStatus::Float => "float",
            FixStatus::Single => "single",
            FixStatus::Ppp => "ppp",
        }
    }
}

/// The position of the rover antenna at one epoch.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EpochSolution {
    /// GPS seconds of the position: the receive time, that is the epoch tag less the receiver's clock offset.
    pub t: f64,
    /// The epoch tag as written in the observation file, GPS seconds.
    pub t_tag: f64,
    /// How the epoch was solved.
    pub status: FixStatus,
    /// ECEF metres, in the frame of the base coordinate (NaN when `status` is `None`).
    pub xyz: [f64; 3],
    /// Latitude and longitude in degrees, ellipsoidal height in metres.
    pub llh: [f64; 3],
    /// Realistic covariance in east, north, up, m².
    pub cov_enu: [[f64; 3]; 3],
    /// Satellites in the solution.
    pub n_sats: u16,
    /// Ratio of the ambiguity search (second-best over best), 0 when no search ran.
    pub ratio: f32,
    /// Seconds between the rover epoch and the nearest base epoch used.
    pub age: f32,
    /// Doppler velocity, east / north / up, m/s (NaN without Doppler).
    pub vel_enu: [f32; 3],
}

impl EpochSolution {
    /// An epoch without a position.
    #[must_use]
    pub fn none(t_tag: f64) -> Self {
        Self {
            t: t_tag,
            t_tag,
            status: FixStatus::None,
            xyz: [f64::NAN; 3],
            llh: [f64::NAN; 3],
            cov_enu: [[f64::NAN; 3]; 3],
            n_sats: 0,
            ratio: 0.0,
            age: 0.0,
            vel_enu: [f32::NAN; 3],
        }
    }

    /// Standard deviations east, north, up (m).
    #[must_use]
    pub fn sigma_enu(&self) -> [f64; 3] {
        [
            self.cov_enu[0][0].max(0.0).sqrt(),
            self.cov_enu[1][1].max(0.0).sqrt(),
            self.cov_enu[2][2].max(0.0).sqrt(),
        ]
    }

    /// Sets `xyz` and the matching `llh`.
    pub fn set_xyz(&mut self, xyz: [f64; 3]) {
        self.xyz = xyz;
        self.llh = llh_deg(xyz);
    }

    /// Whether the epoch has a position.
    #[must_use]
    pub fn has_position(&self) -> bool {
        self.status != FixStatus::None && self.xyz[0].is_finite()
    }
}

/// ECEF to latitude, longitude (degrees) and ellipsoidal height (m).
#[must_use]
pub fn llh_deg(xyz: [f64; 3]) -> [f64; 3] {
    if !xyz[0].is_finite() {
        return [f64::NAN; 3];
    }
    let (lat, lon, h) = ecef_to_llh(xyz);
    [lat.to_degrees(), lon.to_degrees(), h]
}

/// East / north / up of `xyz - origin` at `origin`.
#[must_use]
pub fn enu_at(origin: [f64; 3], xyz: [f64; 3]) -> [f64; 3] {
    let (lat, lon, _) = ecef_to_llh(origin);
    crate::linalg::mat3_vec(&enu_matrix(lat, lon), sub3(xyz, origin))
}

/// Where the base is.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BaseStation {
    /// ECEF metres of the mark. The result is in the frame and epoch of this coordinate.
    pub marker_xyz: [f64; 3],
    /// Antenna reference point above the mark: height, east, north (m). `None`: the base file's
    /// `ANTENNA: DELTA H/E/N`.
    pub delta_hen: Option<[f64; 3]>,
}

/// Everything [`solve`] needs.
#[derive(Debug)]
pub struct PpkInput<'a> {
    /// The static base.
    pub base: &'a StationObs,
    /// The moving receiver.
    pub rover: &'a StationObs,
    /// Orbits and clocks (broadcast, precise, or both).
    pub orbits: &'a Orbits,
    /// GLONASS frequency channels.
    pub fcn: &'a BTreeMap<Sat, i32>,
    /// The base coordinate.
    pub base_station: BaseStation,
    /// Receiver antenna calibrations.
    pub antex: Option<&'a mut Antex>,
    /// Troposphere grids for the `gpt3` / `vmf3` models (empty: the standard atmosphere is used).
    pub tropo: &'a TropoInputs<'a>,
    /// Times to position.
    pub events: &'a [Event],
    /// Stop-and-go occupations to average.
    pub occupations: &'a [Occupation],
}

/// What [`solve`] is doing, for a progress bar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Stage {
    /// Satellite positions, corrections, cycle slips.
    Prepare,
    /// The forward pass.
    Forward,
    /// The backward pass.
    Backward,
    /// Merging the passes.
    Combine,
    /// Interpolating the events.
    Events,
}

/// A progress report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Progress {
    /// The step.
    pub stage: Stage,
    /// Epochs (or events) done in this step.
    pub done: usize,
    /// Epochs (or events) in this step.
    pub total: usize,
}

/// Why [`solve`] returned no result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PpkError {
    /// The progress callback asked to stop.
    Cancelled,
    /// The data cannot be processed (the message says why).
    Input(String),
}

impl std::fmt::Display for PpkError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PpkError::Cancelled => f.write_str("cancelled"),
            PpkError::Input(s) => f.write_str(s),
        }
    }
}

impl std::error::Error for PpkError {}

/// Counts and checks over a whole run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct PpkSummary {
    /// Rover epochs.
    pub n_epochs: usize,
    /// Epochs fixed.
    pub n_fixed: usize,
    /// Epochs float (for precise point positioning: its epochs, which are float by nature).
    pub n_float: usize,
    /// Epochs with a code-only position.
    pub n_single: usize,
    /// Epochs without a position.
    pub n_none: usize,
    /// Fixed share of the epochs with a position, percent.
    pub fix_percent: f64,
    /// Shortest and longest distance between the receivers, m.
    pub baseline_m: [f64; 2],
    /// Epochs fixed in both passes.
    pub both_fixed: usize,
    /// RMS of the forward-minus-backward differences over those epochs: east, north, up (m).
    pub fb_rms_enu: [f64; 3],
    /// The largest such difference (3-D, m).
    pub fb_max: f64,
    /// Epochs whose passes disagreed and were demoted to float.
    pub fb_conflicts: usize,
    /// Fixed share of the forward pass, percent.
    pub fix_percent_forward: f64,
    /// Fixed share of the backward pass, percent.
    pub fix_percent_backward: f64,
    /// New ambiguities started (first appearances, returns after a break, cycle slips), all satellites
    /// and frequencies.
    pub slips: usize,
    /// Cycle slips of whole cycles that were taken out of the phase instead, the ambiguity carried on.
    pub slips_repaired: usize,
    /// Median interval of the rover epochs, s.
    pub interval_s: f64,
    /// The signals used per system, e.g. `"G L1C/L1C L2W/L2W"` (rover / base).
    pub signals: Vec<String>,
    /// Events with a position.
    pub events_ok: usize,
    /// Events interpolated across a gap.
    pub events_gap: usize,
    /// Events without a position.
    pub events_failed: usize,
}

/// The result of [`solve`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct PpkResult {
    /// The track asked for ([`PpkOptions::direction`]), one entry per rover epoch.
    pub epochs: Vec<EpochSolution>,
    /// The forward pass (kept when [`PpkOptions::keep_passes`]).
    pub forward: Vec<EpochSolution>,
    /// The backward pass (kept when [`PpkOptions::keep_passes`]).
    pub backward: Vec<EpochSolution>,
    /// One result per input event, in the input order.
    pub events: Vec<EventResult>,
    /// One result per input occupation.
    pub occupations: Vec<OccupationResult>,
    /// Counts and checks.
    pub summary: PpkSummary,
    /// Things the user must see.
    pub warnings: Vec<String>,
    /// What the solver did.
    pub notes: Vec<String>,
    /// The base antenna reference point used (ECEF m).
    pub base_arp: [f64; 3],
    /// The base mark (ECEF m), as given.
    pub base_marker: [f64; 3],
}

fn counts(track: &[EpochSolution]) -> (usize, usize, usize, usize) {
    let mut c = (0, 0, 0, 0);
    for e in track {
        match e.status {
            FixStatus::Fixed => c.0 += 1,
            FixStatus::Float | FixStatus::Ppp => c.1 += 1,
            FixStatus::Single => c.2 += 1,
            FixStatus::None => c.3 += 1,
        }
    }
    c
}

fn fix_percent(track: &[EpochSolution]) -> f64 {
    let (fx, fl, sg, _) = counts(track);
    let n = fx + fl + sg;
    if n == 0 {
        0.0
    } else {
        100.0 * fx as f64 / n as f64
    }
}

/// Kinematic processing of `input.rover` against `input.base`.
///
/// `progress` is called a few times per second; returning `false` stops the run with
/// [`PpkError::Cancelled`].
pub fn solve(
    input: &mut PpkInput<'_>,
    opts: &PpkOptions,
    progress: &mut dyn FnMut(Progress) -> bool,
) -> Result<PpkResult, PpkError> {
    let mut res = PpkResult::default();
    let setup = prepare::Setup::new(input, opts, &mut res.notes, &mut res.warnings)?;
    res.base_arp = setup.base_arp;
    res.base_marker = input.base_station.marker_xyz;
    res.summary.signals = setup.plan.describe();

    // 1. every epoch prepared once; both passes read the same list
    let (prepared, stats) = prepare::prepare_all(&setup, input, opts, progress, &mut res.warnings)?;
    res.summary.slips_repaired = stats.repaired;
    let n = prepared.len();
    res.summary.n_epochs = n;
    res.summary.slips = prepared
        .iter()
        .map(|p| {
            p.sats
                .iter()
                .map(|s| s.bands().iter().filter(|b| b.slip).count())
                .sum::<usize>()
        })
        .sum();
    res.summary.interval_s =
        events::median_interval(&prepared.iter().map(|p| p.t_tag).collect::<Vec<_>>());

    // 2. the passes
    let run = |forward: bool, progress: &mut dyn FnMut(Progress) -> bool| {
        let mut f = filter::PpkFilter::new(&setup, opts, !forward);
        let mut out = vec![EpochSolution::none(0.0); n];
        let stage = if forward {
            Stage::Forward
        } else {
            Stage::Backward
        };
        for k in 0..n {
            let i = if forward { k } else { n - 1 - k };
            out[i] = f.step(&prepared[i]);
            if k % 500 == 0
                && !progress(Progress {
                    stage,
                    done: k,
                    total: n,
                })
            {
                return Err(PpkError::Cancelled);
            }
        }
        Ok(out)
    };
    let (fwd, bwd) = match opts.direction {
        Direction::Forward => (run(true, progress)?, Vec::new()),
        Direction::Backward => (Vec::new(), run(false, progress)?),
        Direction::Combined => {
            // the two passes are independent: the backward one runs on another thread
            let (a, b) = std::thread::scope(|s| {
                let h = s.spawn(|| run(false, &mut |_| true));
                let a = run(true, progress);
                (a, h.join())
            });
            let b = b.map_err(|_| PpkError::Input("the backward pass failed".into()))?;
            (a?, b?)
        }
    };
    if !progress(Progress {
        stage: Stage::Combine,
        done: 0,
        total: n,
    }) {
        return Err(PpkError::Cancelled);
    }
    res.epochs = match opts.direction {
        Direction::Forward => fwd.clone(),
        Direction::Backward => bwd.clone(),
        Direction::Combined => combine::combine(&fwd, &bwd, opts, &mut res.summary),
    };
    res.summary.fix_percent_forward = fix_percent(&fwd);
    res.summary.fix_percent_backward = fix_percent(&bwd);
    let (fx, fl, sg, no) = counts(&res.epochs);
    res.summary.n_fixed = fx;
    res.summary.n_float = fl;
    res.summary.n_single = sg;
    res.summary.n_none = no;
    res.summary.fix_percent = fix_percent(&res.epochs);
    let mut lo = f64::INFINITY;
    let mut hi = 0.0f64;
    for e in res.epochs.iter().filter(|e| e.has_position()) {
        let d = norm3(sub3(e.xyz, setup.base_arp));
        lo = lo.min(d);
        hi = hi.max(d);
    }
    res.summary.baseline_m = if lo.is_finite() { [lo, hi] } else { [0.0, 0.0] };
    if hi > 1000.0 * opts.max_baseline_km {
        res.warnings.push(format!(
            "the receivers are up to {:.1} km apart: beyond about {:.0} km the ionosphere no longer cancels \
             between them, and this model loses accuracy and fixes",
            hi / 1000.0,
            opts.max_baseline_km
        ));
    }
    if n > 0 && fx + fl + sg == 0 {
        res.warnings.push("no epoch could be positioned".into());
    }

    // 3. events and occupations
    if !progress(Progress {
        stage: Stage::Events,
        done: 0,
        total: input.events.len(),
    }) {
        return Err(PpkError::Cancelled);
    }
    res.events = events::interpolate(
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
    res.occupations = events::average_occupations(&res.epochs, input.occupations);
    if opts.keep_passes {
        res.forward = fwd;
        res.backward = bwd;
    }
    Ok(res)
}
