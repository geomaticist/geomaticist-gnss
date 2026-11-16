//! Options, result and the realistic σ of the PPP.

use std::collections::BTreeMap;

use geomaticist_gnss_core::{ObsCode, Sat, StationObs};
use geomaticist_gnss_engine::dump::{jfx, jsl};
use geomaticist_gnss_engine::linalg::rotate_cov3;
use geomaticist_gnss_engine::models::{ecef_to_llh, enu_matrix};
use geomaticist_gnss_engine::options::SlipOptions;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

/// The order systems are reported and chosen as reference in.
pub const SYS_ORDER: &str = "GECJR";

/// System letter -> name.
#[must_use]
pub fn sys_name(s: char) -> &'static str {
    match s {
        'G' => "GPS",
        'R' => "GLONASS",
        'E' => "Galileo",
        'C' => "BeiDou",
        'J' => "QZSS",
        _ => "?",
    }
}

/// Position of a system letter in [`SYS_ORDER`] (unknown letters last).
#[must_use]
pub fn sys_rank(s: char) -> usize {
    SYS_ORDER.find(s).unwrap_or(SYS_ORDER.len())
}

/// The wide-lane σ rule of PPP-AR. `Ar1` is the default; `Legacy` the earlier rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum WlEstimator {
    /// `σ = max(std(MW) / sqrt(n / 10), 0.02)`: one independent value per ten epochs, whatever the rate.
    Legacy,
    /// `σ = max(std(MW) / sqrt(span / ar_wl_corr_s), 0.02)`: one independent value per correlation time.
    #[default]
    Ar1,
}

/// The PPP options. There is no `blq` field: the BLQ coefficients are handed to [`crate::solve::solve`] and applied when `otl` is on.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
#[allow(missing_docs)] // the field names are the established option names
pub struct PppOptions {
    pub systems: String,
    pub interval: f64,
    pub elmask: f64,
    pub err_phase: f64,
    pub code_ratio: f64,
    pub efact: BTreeMap<String, f64>,
    pub glo_code: bool,
    pub ztd_step_h: f64,
    pub ztd_rw: f64,
    pub ztd_prior: f64,
    pub ztd_sigma: f64,
    pub trop_model: String,
    pub trop_grad_apriori: bool,
    pub gradients: String,
    pub gradient_min_h: f64,
    pub gradient_sigma: f64,
    pub gradient_step_h: f64,
    pub gradient_rw: f64,
    pub tides: bool,
    pub pole_tide: bool,
    pub windup: bool,
    pub shapiro: bool,
    pub sat_pcv: bool,
    pub osb: bool,
    pub shadow: bool,
    pub shadow_after_s: f64,
    pub exclude: Vec<Sat>,
    pub min_arc_s: f64,
    pub slip_pc: f64,
    pub slips: SlipOptions,
    pub reject_k: f64,
    pub reject_rounds: usize,
    pub step_k: f64,
    pub step_k_local: f64,
    pub step_min: f64,
    pub code_gross: f64,
    pub max_iter: usize,
    pub uncal: String,
    pub uncal_sigma: f64,
    pub ar: bool,
    pub ar_systems: String,
    pub ar_min_arc_s: f64,
    pub ar_ratio: f64,
    pub ar_min_fix: usize,
    pub ar_wl_frac: f64,
    pub ar_wl_sig: f64,
    pub ar_wl_corr_s: f64,
    pub ar_nl_offset: f64,
    pub ar_nl_floor: f64,
    pub debug: bool,
    pub debug_dir: Option<String>,
    /// Rust only: apply the ocean tide loading of the BLQ coefficients handed to `solve`.
    pub otl: bool,
    /// Rust only: the wide-lane σ rule.
    pub wl_estimator: WlEstimator,
}

impl Default for PppOptions {
    fn default() -> Self {
        Self {
            systems: "GRE".into(),
            interval: 30.0,
            elmask: 10.0,
            err_phase: 0.003,
            code_ratio: 100.0,
            efact: [('G', 1.0), ('R', 1.5), ('E', 1.0), ('C', 1.5), ('J', 1.0)]
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
            glo_code: true,
            ztd_step_h: 1.0,
            ztd_rw: 1.5e-4,
            ztd_prior: 0.10,
            ztd_sigma: 0.5,
            trop_model: "niell".into(),
            trop_grad_apriori: false,
            gradients: "auto".into(),
            gradient_min_h: 3.0,
            gradient_sigma: 0.002,
            gradient_step_h: 6.0,
            gradient_rw: 0.0005,
            tides: true,
            pole_tide: true,
            windup: true,
            shapiro: true,
            sat_pcv: true,
            osb: true,
            shadow: true,
            shadow_after_s: 1800.0,
            exclude: Vec::new(),
            min_arc_s: 600.0,
            slip_pc: 30.0,
            slips: SlipOptions {
                max_gap: 4,
                slip_gf: 0.08,
                repair: false,
                spikes: true,
                ..SlipOptions::default()
            },
            reject_k: 4.0,
            reject_rounds: 8,
            step_k: 8.0,
            step_k_local: 6.0,
            step_min: 0.03,
            code_gross: 20.0,
            max_iter: 6,
            uncal: "drop".into(),
            uncal_sigma: 0.10,
            ar: false,
            ar_systems: "GEC".into(),
            ar_min_arc_s: 1200.0,
            ar_ratio: 3.0,
            ar_min_fix: 4,
            ar_wl_frac: 0.25,
            ar_wl_sig: 0.15,
            ar_wl_corr_s: 120.0,
            ar_nl_offset: 0.30,
            ar_nl_floor: 0.0,
            debug: false,
            debug_dir: None,
            otl: false,
            wl_estimator: WlEstimator::Ar1,
        }
    }
}

impl PppOptions {
    /// `efact.get(s, 1.0)`.
    #[must_use]
    pub fn efact_of(&self, sys: char) -> f64 {
        self.efact.get(&sys.to_string()).copied().unwrap_or(1.0)
    }
}

/// The realistic σ model of a float PPP (fitted on CORS station-days).
pub const SIGMA: [(&str, f64); 8] = [
    ("kH", 0.0),
    ("aH", 0.0042),
    ("dH", 0.0295),
    ("qH", 2.0),
    ("kU", 0.0),
    ("aU", 0.0083),
    ("dU", 0.0390),
    ("qU", 1.5),
];

/// The same for an ambiguity-fixed solution.
pub const SIGMA_FIXED: [(&str, f64); 8] = [
    ("kH", 0.0),
    ("aH", 0.0013),
    ("dH", 0.0200),
    ("qH", 1.0),
    ("kU", 0.0),
    ("aU", 0.0051),
    ("dU", 0.0385),
    ("qU", 1.0),
];

/// The ENU covariance to weigh a PPP position by (m²): the σ of `SIGMA` (`SIGMA_FIXED` when fixed) per
/// component, with the formal solution's correlations kept.
#[must_use]
pub fn realistic_cov(cov_enu: &[[f64; 3]; 3], hours: f64, fixed: bool) -> [[f64; 3]; 3] {
    let p: BTreeMap<&str, f64> = if fixed { SIGMA_FIXED } else { SIGMA }
        .into_iter()
        .collect();
    let g = |k: &str| p.get(k).copied().unwrap_or(0.0);
    let t = hours.max(0.25);
    let f: [f64; 3] = [0, 1, 2].map(|i| cov_enu[i][i].max(1e-30).sqrt());
    let h2 = |i: usize| {
        g("kH").powi(2) * f[i].powi(2) + g("aH").powi(2) + g("dH").powi(2) / t.powf(g("qH"))
    };
    let u2 = g("kU").powi(2) * f[2].powi(2) + g("aU").powi(2) + g("dU").powi(2) / t.powf(g("qU"));
    let s = [h2(0).sqrt(), h2(1).sqrt(), u2.sqrt()];
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = cov_enu[i][j] / (f[i] * f[j]) * (s[i] * s[j]);
        }
    }
    out
}

/// RINEX 3 observation name (`C1W`, `L2W` …) of a code as stored (RINEX 2 `P1`, `L2` mapped).
#[must_use]
pub fn obs_name(st: &StationObs, sat: Sat, code: Option<ObsCode>) -> Option<String> {
    let code = code?;
    let s = code.as_str();
    if s.len() == 3 {
        return Some(s.to_owned());
    }
    let k = st.sats.get(&sat)?.get(&code)?.key;
    let sys = sat.sys.letter();
    if k.kind == 'L' {
        // _R2_PHASE
        return match (sys, k.band) {
            ('G', 1) => Some("L1C"),
            ('G', 2) => Some("L2W"),
            ('R', 1) => Some("L1C"),
            ('R', 2) => Some("L2P"),
            ('E', 1) => Some("L1C"),
            ('E', 5) => Some("L5Q"),
            ('J', 1) => Some("L1C"),
            ('J', 2) => Some("L2L"),
            _ => None,
        }
        .map(str::to_owned);
    }
    // _R2_CODE with `attr or "C"`
    let attr = k.attr.unwrap_or('C');
    match (sys, attr, k.band) {
        ('G', 'P', 1) => Some("C1W"),
        ('G', 'C', 1) => Some("C1C"),
        ('G', 'P', 2) => Some("C2W"),
        ('G', 'C', 2) => Some("C2L"),
        ('R', 'P', 1) => Some("C1P"),
        ('R', 'C', 1) => Some("C1C"),
        ('R', 'P', 2) => Some("C2P"),
        ('R', 'C', 2) => Some("C2C"),
        ('E', 'C', 1) => Some("C1C"),
        ('E', 'C', 5) => Some("C5Q"),
        ('J', 'C', 1) => Some("C1C"),
        ('J', 'C', 2) => Some("C2L"),
        _ => None,
    }
    .map(str::to_owned)
}

/// Per system, the receiver antenna's status per signal, in the order the signals were met (the order
/// the notes list them in).
pub type Antenna = BTreeMap<char, Vec<(String, String)>>;

/// `antenna` as JSON: `{"E": {"E01": "same frequency (G01)", ...}}`.
#[must_use]
pub fn antenna_json(a: &Antenna) -> Value {
    let m: Map<String, Value> = a
        .iter()
        .map(|(s, v)| {
            let inner: Map<String, Value> = v
                .iter()
                .map(|(k, x)| (k.clone(), Value::String(x.clone())))
                .collect();
            (s.to_string(), Value::Object(inner))
        })
        .collect();
    Value::Object(m)
}

/// The PPP result.
#[derive(Debug, Clone, PartialEq)]
pub struct PppResult {
    /// ARP, ECEF (m), the products' frame at the data's epoch.
    pub xyz: [f64; 3],
    /// 3 × 3 ECEF (m²), scaled by σ0².
    pub cov_formal: [[f64; 3]; 3],
    /// σ0.
    pub sigma0: f64,
    /// First / last epoch on the grid (GPS s).
    pub t0: f64,
    /// Last epoch.
    pub t1: f64,
    /// Session length (h).
    pub hours: f64,
    /// Epochs on the grid.
    pub n_epochs: usize,
    /// Phase observations used.
    pub n_obs_phase: usize,
    /// Code observations used.
    pub n_obs_code: usize,
    /// Ambiguity parameters.
    pub n_amb: usize,
    /// Observations rejected (phase + code).
    pub n_rejected: usize,
    /// System -> satellites used.
    pub sats: BTreeMap<char, Vec<Sat>>,
    /// Phase residual RMS (m).
    pub rms_phase: f64,
    /// Code residual RMS (m).
    pub rms_code: f64,
    /// `(t, zenith total delay m, σ)` per node.
    pub ztd: Vec<(f64, f64, f64)>,
    /// System -> code bias against the reference (m).
    pub isb: BTreeMap<char, f64>,
    /// `(t, north, east)` per gradient node.
    pub gradients: Vec<(f64, f64, f64)>,
    /// Notes in order.
    pub notes: Vec<String>,
    /// Wall-clock seconds.
    pub seconds: f64,
    /// Products label (empty here).
    pub products: String,
    /// Whether code OSBs were applied to any satellite.
    pub osb_used: bool,
    /// System -> {signal: status}.
    pub antenna: Antenna,
    /// System -> (vertical antenna offset m, σ) when estimated.
    pub offsets: BTreeMap<char, (f64, f64)>,
    /// Ambiguities fixed (PPP-AR): `xyz` / `cov_formal` are the fixed solution.
    pub fixed: bool,
    /// The AR record (n_fixed, n_cand, ratio, wl, xyz_float, cov_float …), as a JSON object.
    pub ar: Map<String, Value>,
    /// The troposphere model used.
    pub troposphere: String,
}

impl PppResult {
    /// The ENU covariance at the solution.
    #[must_use]
    pub fn enu_cov(&self) -> [[f64; 3]; 3] {
        let (lat, lon, _) = ecef_to_llh(self.xyz);
        let r = enu_matrix(lat, lon);
        rotate_cov3(&r, &self.cov_formal)
    }

    /// The result as JSON (the public fields but `seconds` and `debug`, plus `enu_cov`; NaN as the string "NaN").
    #[must_use]
    pub fn to_json(&self) -> Value {
        let fl = |v: &[f64]| Value::Array(v.iter().map(|&x| jfx(x)).collect());
        let m3 = |m: &[[f64; 3]; 3]| Value::Array(m.iter().map(|r| fl(r)).collect());
        let mut m = Map::new();
        m.insert("xyz".into(), fl(&self.xyz));
        m.insert("cov_formal".into(), m3(&self.cov_formal));
        m.insert("sigma0".into(), jfx(self.sigma0));
        m.insert("t0".into(), jfx(self.t0));
        m.insert("t1".into(), jfx(self.t1));
        m.insert("hours".into(), jfx(self.hours));
        m.insert("n_epochs".into(), Value::from(self.n_epochs));
        m.insert("n_obs_phase".into(), Value::from(self.n_obs_phase));
        m.insert("n_obs_code".into(), Value::from(self.n_obs_code));
        m.insert("n_amb".into(), Value::from(self.n_amb));
        m.insert("n_rejected".into(), Value::from(self.n_rejected));
        let sats: Map<String, Value> = self
            .sats
            .iter()
            .map(|(k, v)| {
                (
                    k.to_string(),
                    jsl(&v.iter().map(ToString::to_string).collect::<Vec<_>>()),
                )
            })
            .collect();
        m.insert("sats".into(), Value::Object(sats));
        m.insert("rms_phase".into(), jfx(self.rms_phase));
        m.insert("rms_code".into(), jfx(self.rms_code));
        m.insert(
            "ztd".into(),
            Value::Array(self.ztd.iter().map(|(a, b, c)| fl(&[*a, *b, *c])).collect()),
        );
        let isb: Map<String, Value> = self
            .isb
            .iter()
            .map(|(k, v)| (k.to_string(), jfx(*v)))
            .collect();
        m.insert("isb".into(), Value::Object(isb));
        m.insert(
            "gradients".into(),
            Value::Array(
                self.gradients
                    .iter()
                    .map(|(a, b, c)| fl(&[*a, *b, *c]))
                    .collect(),
            ),
        );
        m.insert("notes".into(), jsl(&self.notes));
        m.insert("products".into(), Value::String(self.products.clone()));
        m.insert("osb_used".into(), Value::Bool(self.osb_used));
        m.insert("antenna".into(), antenna_json(&self.antenna));
        let off: Map<String, Value> = self
            .offsets
            .iter()
            .map(|(k, (v, s))| (k.to_string(), fl(&[*v, *s])))
            .collect();
        m.insert("offsets".into(), Value::Object(off));
        m.insert("fixed".into(), Value::Bool(self.fixed));
        m.insert("ar".into(), Value::Object(self.ar.clone()));
        m.insert(
            "troposphere".into(),
            Value::String(self.troposphere.clone()),
        );
        m.insert("enu_cov".into(), m3(&self.enu_cov()));
        Value::Object(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_and_sigma() {
        let o = PppOptions::default();
        assert_eq!(o.slips.max_gap, 4);
        assert!(!o.slips.repair && o.slips.zero_test);
        assert_eq!(o.wl_estimator, WlEstimator::Ar1);
        let v = serde_json::to_value(&o).unwrap();
        assert_eq!(v["ar_wl_sig"], 0.15);
        let back: PppOptions =
            serde_json::from_value(serde_json::json!({"ar": true, "blq": null})).unwrap();
        assert!(back.ar && back.interval == 30.0);
        // realistic σ at 4 h, float: H 8.5 mm, U 16 mm (P5_SIGMA)
        let c = realistic_cov(
            &[[1e-6, 0.0, 0.0], [0.0, 1e-6, 0.0], [0.0, 0.0, 4e-6]],
            4.0,
            false,
        );
        assert!((c[0][0].sqrt() - 0.0085).abs() < 0.0003);
        assert!((c[2][2].sqrt() - 0.016).abs() < 0.001);
        assert_eq!(c[0][1], 0.0);
    }
}
