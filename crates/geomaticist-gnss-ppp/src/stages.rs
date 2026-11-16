//! Typed results of the seven stages of a PPP solve. [`crate::dump`] writes them as
//! `NN_stage.npz / .json` when debug output is on.

use std::collections::BTreeMap;

use geomaticist_gnss_core::Sat;
use geomaticist_gnss_engine::arr::{B2, F2, I2, Vec3Arr};
use geomaticist_gnss_engine::linalg::DMat;
use serde_json::{Map, Value};

use crate::options::Antenna;

/// Stage 1: the per-satellite inputs on the clock grid.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)] // the names are the dump keys
pub struct ObservationsStage {
    pub t: Vec<f64>,
    pub t_obs: Vec<f64>,
    pub cols: Vec<Sat>,
    pub phase: F2,
    pub code: F2,
    pub gf: F2,
    pub mw: F2,
    pub lli: B2,
    pub sxyz: Vec3Arr,
    pub sclk: F2,
    pub lamif: F2,
    pub ok: B2,
    pub ok_raw: B2,
    pub shadow: B2,
    pub pcvs: F2,
    pub az: F2,
    pub el: F2,
    pub windup: F2,
    pub disp: Vec<[f64; 3]>,
    // meta
    pub station: String,
    pub x0: [f64; 3],
    pub pairs: BTreeMap<Sat, (u8, u8)>,
    /// `"G L1C L2W"` per satellite (the phase signals the OSB was looked up for; `None` where unknown).
    pub signals: Vec<String>,
    pub a_if: Vec<f64>,
    pub b_if: Vec<f64>,
    pub lam1: Vec<f64>,
    pub lam2: Vec<f64>,
    /// Whether all four OSBs of the satellite's signals were found.
    pub osb: Vec<bool>,
    pub osb_hits: usize,
    pub osb_miss: usize,
    pub notes: Vec<String>,
}

/// Stage 2: arcs.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct ArcsStage {
    pub arc: I2,
    pub amb0: F2,
    pub ok: B2,
    pub phase_ok: B2,
    pub code_ok: B2,
    pub n_arcs: usize,
    pub antenna: Antenna,
    pub off_sys: Vec<char>,
    pub notes: Vec<String>,
}

/// Stage 3: the float solution.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct FloatStage {
    pub x: [f64; 3],
    pub n: DMat,
    pub ninv: DMat,
    pub sol: Vec<f64>,
    pub resp: F2,
    pub resc: F2,
    pub clk: Vec<f64>,
    pub arc: I2,
    pub use_p: B2,
    pub use_c: B2,
    pub rej_phase: B2,
    pub rej_code: B2,
    pub amb0: F2,
    pub sigma0: f64,
    pub cov: [[f64; 3]; 3],
    pub npar: usize,
    pub i_i: usize,
    pub i_g: usize,
    pub i_o: usize,
    pub i_z: usize,
    pub i_t: usize,
    pub i_a: usize,
    pub nz: usize,
    pub ng: usize,
    pub namb: usize,
    pub isb_sys: Vec<char>,
    pub ref_sys: char,
    /// Column indexes of the GLONASS satellites with their own code bias.
    pub glo_sats: Vec<usize>,
    pub off_sys: Vec<char>,
    pub n_rej: usize,
    pub n_split: usize,
    pub nobs: usize,
    pub notes: Vec<String>,
}

/// Stage 4: the wide lanes.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct WideLaneStage {
    pub ar: bool,
    /// `[fixed, tried]`.
    pub wl: Option<[usize; 2]>,
    pub why: Option<String>,
    pub ok: Option<bool>,
}

/// Stage 5: the narrow lanes.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct NarrowLaneStage {
    pub ar: Map<String, Value>,
}

/// Stage 6: the answer.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct ValidationStage {
    pub x: [f64; 3],
    pub cov: [[f64; 3]; 3],
    pub fixed: bool,
    pub sigma0: f64,
    pub n_rejected: usize,
    pub n_split: usize,
    pub notes: Vec<String>,
}

/// Stage 7: the realistic σ.
#[derive(Debug, Clone, PartialEq)]
#[allow(missing_docs)]
pub struct SigmaStage {
    pub cov_formal: [[f64; 3]; 3],
    pub cov_enu: [[f64; 3]; 3],
    pub cov_realistic_enu: [[f64; 3]; 3],
    pub sigma0: f64,
    pub hours: f64,
    pub rms_phase: f64,
    pub rms_code: f64,
    pub n_obs_phase: usize,
    pub n_obs_code: usize,
    pub n_amb: usize,
    pub n_rejected: usize,
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
