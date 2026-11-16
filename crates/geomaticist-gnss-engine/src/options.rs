//! The engine's options and result types, the realistic-σ calibration, and the antenna-name helpers.
//!
//! Field names and defaults are fixed, so options saved as JSON deserialise directly.

use std::collections::BTreeMap;

use geomaticist_gnss_core::Sat;
use serde::{Deserialize, Serialize};

use crate::linalg::{mat3_mul, mat3_t, rotate_cov3};

/// The realistic σ, fitted on the eight benchmark datasets, November 2026:
/// fixed solutions.
pub const SIGMA_FIXED: [(&str, f64); 18] = [
    ("k", 4.0),
    ("kU", 4.85),
    ("aH", 0.0009),
    ("aU", 0.0016),
    ("dH", 0.0018),
    ("dU", 0.0046),
    ("bH", 0.036),
    ("bU", 0.0),
    ("tau", 0.41),
    ("eU", 0.0),
    ("gf", 0.0),
    ("gr", 0.0),
    ("gs", 0.52),
    ("gd", 0.29),
    ("g0", 0.04),
    ("gx", 0.70),
    ("gu", 0.0),
    ("mp", 1.0),
];
/// Floats.
pub const SIGMA_FLOAT: [(&str, f64); 18] = [
    ("k", 4.6),
    ("kU", 3.96),
    ("aH", 0.0020),
    ("aU", 0.0),
    ("dH", 0.0010),
    ("dU", 0.0132),
    ("bH", 0.0),
    ("bU", 0.081),
    ("tau", 0.0),
    ("eU", 0.0),
    ("gf", 0.0),
    ("gr", 0.0),
    ("gs", 0.0),
    ("gd", 0.0),
    ("g0", 0.166),
    ("gx", 0.327),
    ("gu", 0.0),
    ("mp", 1.5),
];
/// The fixed-failure-rate ratio curve (n DDs, ratio), after RTKLIB-EX 2.5.1 (T. Takasu, BSD-2).
pub const RATIO_FFR: [(f64, f64); 14] = [
    (2.0, 9.79),
    (3.0, 6.87),
    (4.0, 5.33),
    (5.0, 4.41),
    (6.0, 3.79),
    (7.0, 3.35),
    (8.0, 3.02),
    (9.0, 2.77),
    (10.0, 2.57),
    (12.0, 2.27),
    (15.0, 1.98),
    (20.0, 1.70),
    (30.0, 1.40),
    (50.0, 1.20),
];

fn pairs_vec(p: &[(&str, f64)]) -> Vec<(String, f64)> {
    p.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
}

/// Slip detection options.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SlipOptions {
    /// Missing epochs that end an arc.
    pub max_gap: usize,
    /// m: the GF test's floor.
    pub slip_gf: f64,
    /// × the ionosphere's activity.
    pub gf_k: f64,
    /// Epochs in the GF line.
    pub gf_window: usize,
    /// × the MW scatter.
    pub mw_k: f64,
    /// Cycles: the MW test never tighter than this.
    pub mw_floor: f64,
    /// Cycles: the arc's MW σ never taken under this.
    pub mw_sigma_min: f64,
    /// ± epochs for the activity and the MW noise.
    pub activity_window: usize,
    /// Drop one-epoch spikes.
    pub spikes: bool,
    /// Repair certain integer steps.
    pub repair: bool,
    /// Epochs after a step used to measure it.
    pub repair_after: usize,
    /// Cycles: closeness to an integer.
    pub repair_frac: f64,
    /// Cycles: largest σ of each estimate.
    pub repair_sigma: f64,
    /// The zero test (B2).
    pub zero_test: bool,
    /// Epochs after the flag the jumps are measured on.
    pub zero_after: usize,
    /// GF jump within this × its σ of zero.
    pub zero_gf_k: f64,
    /// m: that σ never under this.
    pub zero_gf_floor: f64,
    /// m: the after-line's scatter over this → no decision.
    pub zero_fit_max: f64,
    /// Cycles: the MW jump under this → ΔNw = 0.
    pub zero_mw: f64,
    /// Epochs with MW after the flag it needs.
    pub zero_mw_n: usize,
    /// Cycles: the MW jump's σ over this → no decision.
    pub zero_mw_sigma: f64,
    /// Link the pieces for the wide lane.
    pub wl_link: bool,
    /// Cycles: the link's MW bound.
    pub wl_link_mw: f64,
}

impl Default for SlipOptions {
    fn default() -> Self {
        Self {
            max_gap: 5,
            slip_gf: 0.05,
            gf_k: 6.0,
            gf_window: 10,
            mw_k: 4.0,
            mw_floor: 1.0,
            mw_sigma_min: 0.15,
            activity_window: 10,
            spikes: true,
            repair: true,
            repair_after: 5,
            repair_frac: 0.25,
            repair_sigma: 0.15,
            zero_test: true,
            zero_after: 12,
            zero_gf_k: 3.0,
            zero_gf_floor: 0.003,
            zero_fit_max: 0.03,
            zero_mw: 0.4,
            zero_mw_n: 3,
            zero_mw_sigma: 0.25,
            wl_link: false,
            wl_link_mw: 0.35,
        }
    }
}

/// The solid tide model. `Degree2Step1` is the original model; `IersFull` the full IERS 2010 one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum TideModel {
    /// IERS 2010 step 1, degree 2 (`models::solid_tide_degree2`). The default.
    #[default]
    Degree2Step1,
    /// Full IERS 2010 (steps 1 and 2, degree 3): `tides::solid_tide`.
    IersFull,
}

/// How a group's vᵀPv is taken from its SD residuals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum VtpvForm {
    /// `yᵀWy − (1ᵀWy)² / 1ᵀW1`, floored at 0 (the original expression). With a large SD receiver clock in
    /// the residuals the two terms cancel to their last bits (σ0 off by ~1 % on free-running clocks).
    #[default]
    #[serde(alias = "Oracle")]
    Original,
    /// `Σ w (y − ȳ_w)²`: the same quantity without the cancellation.
    Stable,
}

/// Compatibility flags for the later improvements. Defaults reproduce the original results.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct Compat {
    /// The solid tide applied at both stations.
    pub tide: TideModel,
    /// The vᵀPv expression of the normal-equation groups.
    pub vtpv: VtpvForm,
}

/// The engine's options.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
#[allow(missing_docs)] // the field names are the established option names
pub struct Options {
    pub systems: String,
    pub elmask: f64,
    pub err_a: f64,
    pub err_b: f64,
    pub efact: BTreeMap<String, f64>,
    pub code_ratio: f64,
    pub trop: String,
    pub trop_step_h: f64,
    pub trop_rw: f64,
    pub trop_init: f64,
    pub trop_sigma: f64,
    pub trop_rel: f64,
    pub trop_rel_ppm: f64,
    pub trop_model: String,
    pub tides: bool,
    pub pcv_azimuth: bool,
    pub slip_mode: String,
    pub slips: SlipOptions,
    pub robust: String,
    pub igg_k0: f64,
    pub igg_k1: f64,
    pub step_k: f64,
    pub step_k_local: f64,
    pub step_min: f64,
    pub step_rounds: usize,
    pub slip_gf: f64,
    pub slip_pc: f64,
    pub phase_gross: f64,
    pub phase_gross_rate: f64,
    pub same_signals: bool,
    pub ms_repair: bool,
    pub ms_max_h: f64,
    pub ms_tol: f64,
    pub code_gross_tx: f64,
    pub max_gap: usize,
    pub min_arc: usize,
    pub code_gross: f64,
    pub outlier_k: f64,
    pub max_reject_rounds: usize,
    pub interval: f64,
    pub amb_prior: f64,
    pub amb_rw: f64,
    pub amb_seg_s: f64,
    pub rover_sigma: f64,
    pub fix: bool,
    pub fix_systems: String,
    pub glo_fix: String,
    pub nl_drop_group: bool,
    pub bds_split: bool,
    pub fix_geo: bool,
    pub wl_frac: f64,
    pub wl_sigma: f64,
    pub wl_corr_s: f64,
    pub nl_ratio: f64,
    pub nl_min: usize,
    pub nl_greedy: usize,
    pub nl_p0: f64,
    pub nl_floor: f64,
    pub nl_offset: f64,
    pub nl_sigma: f64,
    pub ar_el: f64,
    pub fix_min_s: f64,
    pub fix_qscale: f64,
    pub fix_s0_max: f64,
    pub nl_ratio_curve: Vec<(f64, f64)>,
    pub nl_small_scale: f64,
    pub iono: String,
    pub iono_km: f64,
    pub iono_sigma: f64,
    pub iono_ppm: f64,
    pub iono_step_s: f64,
    pub iono_max_nodes: usize,
    pub iono_data: f64,
    pub iono_rw: f64,
    pub iono_vce_rounds: usize,
    pub iono_vce_k: f64,
    pub iono_vce_red: f64,
    pub iono_refix_k: f64,
    pub iono_refix_rounds: usize,
    pub iono_final_k: f64,
    pub iono_alt: bool,
    pub iono_wl_float: bool,
    pub iono_wl_scale: f64,
    pub iono_z_max: f64,
    pub iono_sh_max: f64,
    pub iono_nl_min: usize,
    pub iono_share_min: f64,
    pub iono_agree: f64,
    pub iono_agree_k: f64,
    pub iono_conflict: String,
    pub arb_prior_k: f64,
    pub glo_check: String,
    pub cal_rule: String,
    pub cal_systems: String,
    pub fail_sigma: f64,
    pub fail_code: f64,
    pub code_sat_floor: f64,
    pub code_check_glo: usize,
    pub code_bias_m: f64,
    pub code_bias_k: f64,
    pub code_bias_n: usize,
    pub sigma_model: String,
    pub fix_share_min: f64,
    pub sigma_fixed: Vec<(String, f64)>,
    pub sigma_float: Vec<(String, f64)>,
    pub ms_split: bool,
    pub exclude: Vec<Sat>,
    pub epoch_tol: f64,
    pub half_var: f64,
    pub debug: bool,
    pub debug_dir: Option<String>,
    /// The compatibility flags (absent from older saved options).
    pub compat: Compat,
}

impl Default for Options {
    fn default() -> Self {
        let mut efact = BTreeMap::new();
        for (k, v) in [("G", 1.0), ("R", 1.5), ("E", 1.0), ("C", 1.0), ("J", 1.0)] {
            efact.insert(k.to_owned(), v);
        }
        Self {
            systems: "GRECJ".into(),
            elmask: 10.0,
            err_a: 0.003,
            err_b: 0.003,
            efact,
            code_ratio: 300.0,
            trop: "station".into(),
            trop_step_h: 2.0,
            trop_rw: 1e-4,
            trop_init: 0.15,
            trop_sigma: 0.3,
            trop_rel: 0.002,
            trop_rel_ppm: 0.25,
            trop_model: "niell".into(),
            tides: true,
            pcv_azimuth: false,
            slip_mode: "repair".into(),
            slips: SlipOptions::default(),
            robust: "igg3".into(),
            igg_k0: 2.0,
            igg_k1: 5.0,
            step_k: 8.0,
            step_k_local: 6.0,
            step_min: 0.03,
            step_rounds: 4,
            slip_gf: 0.05,
            slip_pc: 30.0,
            phase_gross: 100.0,
            phase_gross_rate: 1.0,
            same_signals: true,
            ms_repair: true,
            ms_max_h: 3.0,
            ms_tol: 15.0,
            code_gross_tx: 1000.0,
            max_gap: 5,
            min_arc: 4,
            code_gross: 30.0,
            outlier_k: 4.0,
            max_reject_rounds: 8,
            interval: 0.0,
            amb_prior: 100.0,
            amb_rw: 0.0,
            amb_seg_s: 300.0,
            rover_sigma: 0.0,
            fix: true,
            fix_systems: "GECJ".into(),
            glo_fix: "same".into(),
            nl_drop_group: true,
            bds_split: true,
            fix_geo: false,
            wl_frac: 0.25,
            wl_sigma: 0.15,
            wl_corr_s: 120.0,
            nl_ratio: 3.0,
            nl_min: 4,
            nl_greedy: 14,
            nl_p0: 0.999,
            nl_floor: 0.0,
            nl_offset: 0.0,
            nl_sigma: 0.35,
            ar_el: 15.0,
            fix_min_s: 600.0,
            fix_qscale: 1.0,
            fix_s0_max: 1.5,
            nl_ratio_curve: RATIO_FFR.to_vec(),
            nl_small_scale: 0.5,
            iono: "auto".into(),
            iono_km: 15.0,
            iono_sigma: 0.002,
            iono_ppm: 2.0,
            iono_step_s: 300.0,
            iono_max_nodes: 1500,
            iono_data: 0.0,
            iono_rw: 0.0,
            iono_vce_rounds: 2,
            iono_vce_k: 1.3,
            iono_vce_red: 4.0,
            iono_refix_k: 1.5,
            iono_refix_rounds: 2,
            iono_final_k: 1.0,
            iono_alt: false,
            iono_wl_float: true,
            iono_wl_scale: 1.5,
            iono_z_max: 0.0,
            iono_sh_max: 0.008,
            iono_nl_min: 5,
            iono_share_min: 0.85,
            iono_agree: 0.02,
            iono_agree_k: 4.0,
            iono_conflict: "float".into(),
            arb_prior_k: 3.0,
            glo_check: "float".into(),
            cal_rule: "own".into(),
            cal_systems: "EC".into(),
            fail_sigma: 0.02,
            fail_code: 1.0,
            code_sat_floor: 0.5,
            code_check_glo: 5,
            code_bias_m: 1.5,
            code_bias_k: 4.0,
            code_bias_n: 20,
            sigma_model: "floor".into(),
            fix_share_min: 0.15,
            sigma_fixed: pairs_vec(&SIGMA_FIXED),
            sigma_float: pairs_vec(&SIGMA_FLOAT),
            ms_split: true,
            exclude: Vec::new(),
            epoch_tol: 0.01,
            half_var: 0.01,
            debug: false,
            debug_dir: None,
            compat: Compat::default(),
        }
    }
}

impl Options {
    /// The σ parameter set as a map.
    #[must_use]
    pub fn sigma_params(&self, fixed: bool) -> BTreeMap<String, f64> {
        let src = if fixed {
            &self.sigma_fixed
        } else {
            &self.sigma_float
        };
        src.iter().cloned().collect()
    }

    /// `efact.get(s, 1.0)`.
    #[must_use]
    pub fn efact_of(&self, sys: char) -> f64 {
        self.efact.get(&sys.to_string()).copied().unwrap_or(1.0)
    }
}

/// The engine's result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StaticResult {
    /// Base ARP (held).
    pub base_xyz: [f64; 3],
    /// Rover ARP.
    pub rover_xyz: [f64; 3],
    /// Rover ARP − base ARP, ECEF (m).
    pub vector: [f64; 3],
    /// 3 × 3 ECEF (m²): realistic (`sigma_model`), what an adjustment weighs by.
    pub cov: [[f64; 3]; 3],
    /// 3 × 3 ECEF (m²): the least-squares covariance, scaled by σ0².
    pub cov_formal: [[f64; 3]; 3],
    /// σ0.
    pub sigma0: f64,
    /// DD phase residual RMS (m).
    pub rms_phase: f64,
    /// DD code residual RMS (m); NaN without code.
    pub rms_code: f64,
    /// Common epochs.
    pub n_epochs: usize,
    /// Phase DDs.
    pub n_dd: usize,
    /// Ambiguity parameters.
    pub n_amb: usize,
    /// Phase observations rejected.
    pub n_rejected: usize,
    /// System -> satellites used.
    pub sats: BTreeMap<String, Vec<Sat>>,
    /// Satellite -> (band1, band2).
    pub pairs: BTreeMap<Sat, (u8, u8)>,
    /// Troposphere: `rover_zwd`, `base_zwd` (rounded to 4 decimals).
    pub trop: BTreeMap<String, Vec<f64>>,
    /// Wall-clock seconds.
    pub seconds: f64,
    /// Code observations rejected.
    pub n_rejected_code: usize,
    /// Notes in order.
    pub notes: Vec<String>,
    /// Satellite -> {n, rms_mm, down, zero, arcs, repaired, spikes, eph_age_min}.
    pub sat_report: BTreeMap<Sat, serde_json::Value>,
    /// Slip event counts plus `split`.
    pub slip_events: BTreeMap<String, usize>,
    /// Whether vector / cov are the integer-fixed solution.
    pub fixed: bool,
    /// The fix record (status, ratio, wl / nl counts, dropped ...), as a JSON object.
    pub fix: serde_json::Map<String, serde_json::Value>,
    /// The float vector.
    pub vector_float: Option<[f64; 3]>,
    /// The float covariance (σ0²-scaled).
    pub cov_float: Option<[[f64; 3]; 3]>,
    /// Not to be trusted (`fix["failed"]` says why).
    pub failed: bool,
}

impl StaticResult {
    /// The result as JSON (the public fields but `seconds` and `debug`, plus `enu` and `cov_enu`; NaN as the string "NaN").
    #[must_use]
    pub fn to_json(&self) -> serde_json::Value {
        use crate::dump::{jfx, jsl};
        use serde_json::{Map, Value};
        let fl = |v: &[f64]| Value::Array(v.iter().map(|&x| jfx(x)).collect());
        let m3 = |m: &[[f64; 3]; 3]| Value::Array(m.iter().map(|r| fl(r)).collect());
        let mut m = Map::new();
        m.insert("base_xyz".into(), fl(&self.base_xyz));
        m.insert("rover_xyz".into(), fl(&self.rover_xyz));
        m.insert("vector".into(), fl(&self.vector));
        m.insert("cov".into(), m3(&self.cov));
        m.insert("cov_formal".into(), m3(&self.cov_formal));
        m.insert("sigma0".into(), jfx(self.sigma0));
        m.insert("rms_phase".into(), jfx(self.rms_phase));
        m.insert("rms_code".into(), jfx(self.rms_code));
        m.insert("n_epochs".into(), Value::from(self.n_epochs));
        m.insert("n_dd".into(), Value::from(self.n_dd));
        m.insert("n_amb".into(), Value::from(self.n_amb));
        m.insert("n_rejected".into(), Value::from(self.n_rejected));
        let sats: Map<String, Value> = self
            .sats
            .iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    jsl(&v.iter().map(ToString::to_string).collect::<Vec<_>>()),
                )
            })
            .collect();
        m.insert("sats".into(), Value::Object(sats));
        let pairs: Map<String, Value> = self
            .pairs
            .iter()
            .map(|(k, (a, b))| {
                (
                    k.to_string(),
                    Value::Array(vec![Value::from(*a), Value::from(*b)]),
                )
            })
            .collect();
        m.insert("pairs".into(), Value::Object(pairs));
        let trop: Map<String, Value> = self.trop.iter().map(|(k, v)| (k.clone(), fl(v))).collect();
        m.insert("trop".into(), Value::Object(trop));
        m.insert("n_rejected_code".into(), Value::from(self.n_rejected_code));
        m.insert("notes".into(), jsl(&self.notes));
        let rep: Map<String, Value> = self
            .sat_report
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect();
        m.insert("sat_report".into(), Value::Object(rep));
        let se: Map<String, Value> = self
            .slip_events
            .iter()
            .map(|(k, v)| (k.clone(), Value::from(*v)))
            .collect();
        m.insert("slip_events".into(), Value::Object(se));
        m.insert("fixed".into(), Value::Bool(self.fixed));
        m.insert("fix".into(), Value::Object(self.fix.clone()));
        m.insert(
            "vector_float".into(),
            self.vector_float.as_ref().map_or(Value::Null, |v| fl(v)),
        );
        m.insert(
            "cov_float".into(),
            self.cov_float.as_ref().map_or(Value::Null, m3),
        );
        m.insert("failed".into(), Value::Bool(self.failed));
        let (enu, cov_enu) = self.enu();
        m.insert("enu".into(), fl(&enu));
        m.insert("cov_enu".into(), m3(&cov_enu));
        Value::Object(m)
    }

    /// The vector and its covariance in the base's local frame (E, N, U).
    #[must_use]
    pub fn enu(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        let (lat, lon, _) = crate::models::ecef_to_llh(self.base_xyz);
        let r = crate::models::enu_matrix(lat, lon);
        (
            crate::linalg::mat3_vec(&r, self.vector),
            rotate_cov3(&r, &self.cov),
        )
    }
}

/// The covariance to weigh a vector by (`sigma_model "floor"`), ECEF.
#[must_use]
pub fn realistic_cov(
    cov_formal: &[[f64; 3]; 3],
    r: &[[f64; 3]; 3],
    length_m: f64,
    hours: f64,
    fixed: bool,
    opts: &Options,
    dh_m: f64,
    q: f64,
) -> [[f64; 3]; 3] {
    if opts.sigma_model == "formal" {
        return *cov_formal;
    }
    let p = opts.sigma_params(fixed);
    let g = |k: &str| p.get(k).copied().unwrap_or(0.0);
    let t = hours.max(0.1);
    let shrink = if g("tau") > 0.0 {
        1.0 / (1.0 + t / g("tau"))
    } else {
        1.0
    };
    let h2 = g("aH").powi(2) + g("dH").powi(2) / t + (g("bH") * 1e-6 * length_m).powi(2) * shrink;
    let u2 = g("aU").powi(2)
        + g("dU").powi(2) / t
        + (g("bU") * 1e-6 * length_m).powi(2) * shrink
        + (g("eU") * dh_m.abs() / 1000.0).powi(2) / t;
    let k = [[g("k"), 0.0, 0.0], [0.0, g("k"), 0.0], [0.0, 0.0, g("kU")]];
    let local = rotate_cov3(r, cov_formal);
    let mut c = mat3_mul(&mat3_mul(&k, &local), &k);
    c[0][0] += h2;
    c[1][1] += h2;
    c[2][2] += u2;
    let rt = mat3_t(r);
    let back = mat3_mul(&mat3_mul(&rt, &c), r);
    let mut out = back;
    for row in &mut out {
        for v in row.iter_mut() {
            *v *= q * q;
        }
    }
    out
}

/// `q`, the data-quality factor on the realistic σ.
#[must_use]
pub fn sigma_quality(
    p: &BTreeMap<String, f64>,
    fixed: bool,
    frag: f64,
    rej: f64,
    s0: f64,
    ratio: Option<f64>,
    fdist: f64,
    unfixed: bool,
    drop: f64,
    partial: bool,
) -> f64 {
    let g = |k: &str| p.get(k).copied().unwrap_or(0.0);
    let mut q = 1.0
        + g("gf") * (frag - 1.0).max(0.0)
        + g("gr") * rej / 0.05
        + g("g0") * (s0 / 0.2 - 1.0).max(0.0)
        + g("gx") * drop / 0.1;
    if fixed {
        q +=
            g("gs") * ratio.filter(|r| *r != 0.0).map_or(0.0, |r| 3.0 / r) + g("gd") * fdist / 10.0;
    } else if unfixed && !partial {
        q += g("gu");
    }
    q * if partial && g("mp") > 0.0 {
        g("mp")
    } else {
        1.0
    }
}

/// A BeiDou geostationary satellite (BDS-2 C01–C05, BDS-3 C59 on).
#[must_use]
pub fn bds_geo(sat: Sat) -> bool {
    sat.sys == geomaticist_gnss_core::GnssSystem::BeiDou && (sat.prn <= 5 || sat.prn >= 59)
}

/// The integer-fixing group of a satellite: its system letter, but `"C2"` for BeiDou-2 when `bds_split`.
#[must_use]
pub fn fix_group(sat: Sat, opts: &Options) -> String {
    if opts.bds_split && sat.sys == geomaticist_gnss_core::GnssSystem::BeiDou && sat.prn <= 18 {
        "C2".to_owned()
    } else {
        sat.sys.letter().to_string()
    }
}

/// A RINEX receiver type's maker.
#[must_use]
pub fn receiver_make(rec: &str) -> String {
    let w = rec.trim().to_uppercase();
    let Some(first) = w.split_whitespace().next() else {
        return String::new();
    };
    match first {
        "SEPTENTRIO" => "SEPT",
        "NOVATEL" => "NOV",
        "TOPCON" => "TPS",
        "ASHTECH" => "ASH",
        other => other,
    }
    .to_owned()
}

/// Whether GLONASS ambiguities are fixed on this pair (`glo_fix`), and why in words.
#[must_use]
pub fn glo_fixable(opts: &Options, rec_a: &str, rec_b: &str) -> (bool, String) {
    let (ma, mb) = (receiver_make(rec_a), receiver_make(rec_b));
    let or_q = |s: &str| {
        if s.is_empty() {
            "?".to_owned()
        } else {
            s.to_owned()
        }
    };
    let or_unknown = |s: &str| {
        if s.is_empty() {
            "unknown".to_owned()
        } else {
            s.to_owned()
        }
    };
    if opts.glo_fix == "off" {
        return (false, "GLONASS float (glo_fix off)".into());
    }
    if opts.glo_fix == "all" {
        return (
            true,
            format!(
                "GLONASS fixed with no inter-channel bias model ({} / {}, glo_fix all)",
                or_q(&ma),
                or_q(&mb)
            ),
        );
    }
    if !ma.is_empty() && ma == mb {
        return (
            true,
            format!("GLONASS fixed: {ma} receivers at both ends (inter-channel biases near zero)"),
        );
    }
    (
        false,
        format!(
            "GLONASS float: receivers of two makes ({} / {}), whose inter-channel biases would put its double \
             differences off their integers",
            or_unknown(&ma),
            or_unknown(&mb)
        ),
    )
}

/// `'TRMR10-2 NONE'` / `'TRMR10-2        NONE'` -> `'TRMR10-2        NONE'` (IGS 16 + 4), right-stripped.
#[must_use]
pub fn norm_ant(name: &str) -> String {
    let parts: Vec<&str> = name.split_whitespace().collect();
    let Some(ant) = parts.first() else {
        return String::new();
    };
    let rad = parts.get(1).copied().unwrap_or("NONE");
    format!("{ant:<16}{rad:>4}").trim_end().to_owned()
}

/// The ANTEX entry to model an antenna with -> (name, note): its own when calibrated, else the same antenna
/// with no radome (NONE), else `""`.
#[must_use]
pub fn antex_name(name: &str, types: &std::collections::BTreeSet<String>) -> (String, String) {
    let n = norm_ant(name);
    if n.is_empty() {
        return (String::new(), String::new());
    }
    if types.contains(&n) {
        return (n, String::new());
    }
    let chars: Vec<char> = n.chars().collect();
    let head: String = chars[..chars.len().min(16)].iter().collect();
    let bare = norm_ant(&format!("{} NONE", head.trim()));
    if types.contains(&bare) {
        let first = n.split_whitespace().next().unwrap_or("");
        let radome: String = if chars.len() > 16 {
            chars[16..].iter().collect::<String>().trim().to_owned()
        } else {
            String::new()
        };
        return (
            bare,
            format!("{first}: radome {radome} isn't calibrated - its NONE calibration used"),
        );
    }
    (String::new(), String::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_helpers() {
        let o = Options::default();
        assert_eq!(o.efact_of('R'), 1.5);
        assert_eq!(o.sigma_params(true)["kU"], 4.85);
        assert_eq!(receiver_make("SEPT POLARX5"), "SEPT");
        assert_eq!(receiver_make("Septentrio PolaRx5"), "SEPT");
        assert!(glo_fixable(&o, "SEPT POLARX5", "SEPT POLARX5").0);
        assert!(!glo_fixable(&o, "TRIMBLE R10-2", "SEPT POLARX5").0);
        assert_eq!(norm_ant("TRMR10-2 NONE"), "TRMR10-2        NONE");
        let mut types = std::collections::BTreeSet::new();
        types.insert("ASH701945G_M    NONE".to_owned());
        let (n, why) = antex_name("ASH701945G_M    JPLA", &types);
        assert_eq!(n, "ASH701945G_M    NONE");
        assert_eq!(
            why,
            "ASH701945G_M: radome JPLA isn't calibrated - its NONE calibration used"
        );
        assert!(bds_geo("C03".parse().unwrap()) && !bds_geo("C20".parse().unwrap()));
        assert_eq!(fix_group("C10".parse().unwrap(), &o), "C2");
        let v: serde_json::Value = serde_json::to_value(&o).unwrap();
        let back: Options = serde_json::from_value(v).unwrap();
        assert_eq!(back, o);
    }

    #[test]
    fn realistic_sigma_formula() {
        let o = Options::default();
        let r = crate::models::enu_matrix(1.1, -2.5);
        let cf = [[1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 4e-6]];
        let c = realistic_cov(&cf, &r, 20000.0, 1.0, true, &o, 0.0, 1.0);
        let local = rotate_cov3(&r, &c);
        let h2 = 0.0009f64.powi(2)
            + 0.0018f64.powi(2)
            + (0.036e-6 * 20000.0f64).powi(2) / (1.0 + 1.0 / 0.41);
        assert!((local[0][0] - (16.0 * 1e-6 + h2)).abs() < 1e-12);
        let p = o.sigma_params(true);
        let q = sigma_quality(&p, true, 1.0, 0.0, 0.18, Some(15.0), 4.0, false, 0.0, false);
        assert!((q - (1.0 + 0.52 * 3.0 / 15.0 + 0.29 * 0.4)).abs() < 1e-12);
    }
}
