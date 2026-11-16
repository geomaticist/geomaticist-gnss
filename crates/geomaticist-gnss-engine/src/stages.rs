//! Typed results of the seven stages of a solve, as structs. [`crate::dump`] writes them as
//! `NN_stage.npz / .json` when `Options { debug: true, debug_dir, .. }`.

use std::collections::BTreeMap;

use geomaticist_gnss_core::Sat;
use serde_json::{Map, Value};

use crate::arr::{B2, F2, I2, Vec3Arr};
use crate::linalg::DMat;
use crate::slips::SlipEvents;

/// Per-station preparation events.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct StationEvents {
    /// Codes off the receiver clock (dropped).
    pub code_ms: usize,
    /// Other phase-code jumps (new arcs).
    pub pc_jump: usize,
    /// Whole-millisecond phase steps put back.
    pub ms_repaired: usize,
}

/// Stage 1: the single-difference inputs per station on the common epochs.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)] // the names are the dump keys
pub struct ObservationsStage {
    pub t: Vec<f64>,
    pub cols: Vec<Sat>,
    pub ph_a: F2,
    pub ph_b: F2,
    pub cd_a: F2,
    pub cd_b: F2,
    pub gf_a: F2,
    pub gf_b: F2,
    pub mw_a: F2,
    pub mw_b: F2,
    pub ok_a: B2,
    pub ok_b: B2,
    pub valid: B2,
    pub lli: B2,
    pub hard: B2,
    pub half_a: B2,
    pub half_b: B2,
    pub ms_a: B2,
    pub ms_b: B2,
    pub el_a: F2,
    pub el_b: F2,
    pub az_a: F2,
    pub az_b: F2,
    pub rho_a: F2,
    pub rho_b: F2,
    pub sat_a: Vec3Arr,
    pub sat_b: Vec3Arr,
    pub clk_a: F2,
    pub clk_b: F2,
    pub tide_a: Vec<[f64; 3]>,
    pub tide_b: Vec<[f64; 3]>,
    pub dry_a: F2,
    pub dry_b: F2,
    pub map_a: F2,
    pub map_b: F2,
    pub ant_a: F2,
    pub ant_b: F2,
    pub lam1: Vec<f64>,
    pub lam2: Vec<f64>,
    pub aif: Vec<f64>,
    pub bif: Vec<f64>,
    // meta
    pub base: String,
    pub rover: String,
    pub base_xyz: [f64; 3],
    pub rover_xyz0: [f64; 3],
    pub pairs: BTreeMap<Sat, (u8, u8)>,
    /// `{sat: {band: [phase code, code code]}}` per station (null where none).
    pub signals_base: Value,
    pub signals_rover: Value,
    pub eph_age_min: BTreeMap<Sat, f64>,
    pub events_base: StationEvents,
    pub events_rover: StationEvents,
    pub notes: Vec<String>,
}

/// Stage 2: arcs.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct ArcsStage {
    pub arc: I2,
    pub starts: B2,
    pub valid: B2,
    pub corr: F2,
    pub corr_mw: F2,
    pub corr_gf: F2,
    pub half: B2,
    pub ms: B2,
    pub n_amb: usize,
    pub n_arc_raw: usize,
    pub wl_links: Vec<(usize, usize)>,
    pub slip_events: SlipEvents,
    /// Satellite -> (repairs, spikes).
    pub sat_slips: BTreeMap<Sat, (usize, usize)>,
    pub notes: Vec<String>,
}

/// The least-squares accumulators.
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize, serde::Deserialize)]
#[allow(missing_docs)]
pub struct Acc {
    pub vtpv: f64,
    pub nrow: usize,
    pub ss_ph: f64,
    pub n_ph: usize,
    pub ss_cd: f64,
    pub n_cd: usize,
}

/// Stage 3: the float solution.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct FloatStage {
    pub x: Vec<f64>,
    pub n: DMat,
    pub q: DMat,
    pub arc: I2,
    pub use_: B2,
    pub wf_ph: F2,
    pub wf_cd: F2,
    pub res_ph: F2,
    pub res_cd: F2,
    pub cf_ph: F2,
    pub var_ph: F2,
    pub var_cd: Option<F2>,
    pub s0: f64,
    pub vector: [f64; 3],
    pub cov: [[f64; 3]; 3],
    pub n_amb: usize,
    pub i_amb: usize,
    pub i_tb: usize,
    pub i_ta: usize,
    pub nn: usize,
    pub n_ta: usize,
    pub p: usize,
    pub acc: Acc,
    /// Column index -> code bias (m, rounded to 2 decimals).
    pub biased: BTreeMap<usize, f64>,
    pub notes: Vec<String>,
}

/// The wide lanes, one entry per ambiguity.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct WideLanes {
    pub nw: Vec<i64>,
    pub ok: Vec<bool>,
    pub span: Vec<f64>,
    pub cnt: Vec<f64>,
    pub jj: Vec<i64>,
    pub sys: Vec<String>,
    pub bias: BTreeMap<String, f64>,
    /// `[n_fix, n_cand]`.
    pub n: [usize; 2],
    pub groups: usize,
    pub mw: Vec<f64>,
    pub sig: Vec<f64>,
    pub cand: Vec<bool>,
}

/// Stage 4: the wide lanes.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct WideLaneStage {
    pub wl: Option<WideLanes>,
    pub keep_float: Vec<bool>,
    pub half_arcs: Option<Vec<bool>>,
    pub w_arc: Option<Vec<f64>>,
    pub pc_arc: Vec<f64>,
}

/// Stage 5: the iono-free narrow-lane fix.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct NarrowLaneStage {
    pub x: Option<Vec<f64>>,
    pub cov: Option<DMat>,
    /// The fix info dict with `s0` added.
    pub info: Map<String, Value>,
}

/// Stage 6: the validated answer.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct ValidationStage {
    pub x: Vec<f64>,
    pub cov: DMat,
    pub fix: Map<String, Value>,
    pub failed: String,
    pub s0: f64,
    pub vector: [f64; 3],
    pub notes: Vec<String>,
}

/// Stage 7: the realistic σ.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct SigmaStage {
    pub cov_formal: [[f64; 3]; 3],
    pub cov: [[f64; 3]; 3],
    pub sigma_q: Option<f64>,
    pub frag: Option<f64>,
    pub rej: Option<f64>,
    pub partial: Option<f64>,
    pub fdist: f64,
    pub s0: f64,
    pub hours: f64,
    pub length_m: f64,
    pub n_rejected: usize,
    pub n_dd: usize,
    pub notes: Vec<String>,
}

/// Every stage a solve went through (those after a failure are `None`).
#[derive(Debug, Clone, Default)]
#[allow(missing_docs)]
pub struct Stages {
    pub observations: Option<ObservationsStage>,
    pub arcs: Option<ArcsStage>,
    pub float: Option<FloatStage>,
    pub wide_lane: Option<WideLaneStage>,
    pub narrow_lane: Option<NarrowLaneStage>,
    pub validation: Option<ValidationStage>,
    pub sigma: Option<SigmaStage>,
}
