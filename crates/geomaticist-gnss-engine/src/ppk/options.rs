//! The settings of kinematic processing. Every field has a default; a JSON object with some of the fields
//! overrides those and leaves the rest.

use geomaticist_gnss_core::Sat;
use serde::{Deserialize, Serialize};

use super::events::EventInterp;

/// Which track [`super::solve`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum Direction {
    /// The filter run forward in time only.
    Forward,
    /// The filter run backward in time only.
    Backward,
    /// Both runs, merged epoch by epoch (the smoothed solution).
    #[default]
    Combined,
}

/// Integer ambiguity resolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum ArMode {
    /// Float solutions only.
    Off,
    /// The integers are searched afresh at every epoch from the float states; a fix never changes the
    /// float states.
    #[default]
    Continuous,
}

/// Settings of [`super::solve`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PpkOptions {
    /// Systems to use, as RINEX letters (`G` GPS, `R` GLONASS, `E` Galileo, `C` BeiDou, `J` QZSS).
    pub systems: String,
    /// Satellites to leave out.
    pub exclude: Vec<Sat>,
    /// Elevation mask, degrees.
    pub elmask: f64,
    /// Most frequencies used per satellite.
    pub max_bands: usize,
    /// Which track to return.
    pub direction: Direction,
    /// Also return the forward and backward tracks.
    pub keep_passes: bool,

    /// Phase noise of one receiver at the zenith, m.
    pub err_phase: f64,
    /// Elevation-dependent phase noise, m: σ² = `err_phase`² + (`err_phase_el` / sin el)².
    pub err_phase_el: f64,
    /// Code noise as a multiple of the phase noise.
    pub code_ratio: f64,
    /// How long code errors last, s (receiver smoothing and multipath make the code of neighbouring
    /// epochs nearly the same error). At an interval shorter than this the code is down-weighted, so that
    /// the float solution does not claim more than the code can give.
    pub code_corr_s: f64,
    /// Extra noise factor of GLONASS code and phase (inter-channel biases).
    pub efact_glonass: f64,
    /// Troposphere model for both ends: `niell` (standard atmosphere), `gpt3`, `vmf3` (need grids).
    pub trop_model: String,
    /// Relative humidity of the standard atmosphere for the wet delay when the model gives none (0 to 1).
    pub humidity: f64,
    /// Use the azimuth-dependent part of the antenna calibrations.
    pub pcv_azimuth: bool,

    /// Integer ambiguity resolution.
    pub ar: ArMode,
    /// The least ratio (second-best over best candidate) for a fix.
    pub ar_ratio: f64,
    /// Elevation above which a satellite's ambiguities are candidates for fixing, degrees.
    pub ar_elmask: f64,
    /// Seconds an ambiguity must have been tracked without a slip before it is a candidate.
    pub ar_min_lock_s: f64,
    /// The least number of double-difference ambiguities in a fix.
    pub ar_min_dd: usize,
    /// The largest RMS of the double-difference phase residuals of a fixed solution, m.
    pub ar_max_res: f64,
    /// Seconds the same integers must keep coming out before an epoch is reported fixed.
    pub ar_confirm_s: f64,
    /// A fixed position may differ from the float one by at most this many float standard deviations
    /// (plus 5 cm).
    pub ar_max_shift_sigma: f64,
    /// Rounds of the partial search (the least certain ambiguities left out) when the full set fails.
    pub ar_partial_rounds: usize,
    /// Fix GLONASS ambiguities (needs receivers without inter-channel phase biases between them; off by
    /// default, GLONASS then helps as float).
    pub ar_glonass: bool,

    /// Standard deviation of a new ambiguity, cycles.
    pub amb_sigma0: f64,
    /// Random walk of the ambiguities, cycles per √s.
    pub amb_rw: f64,
    /// Standard deviation of the position before each epoch's update, m (no dynamics are assumed).
    pub pos_sigma0: f64,
    /// A satellite out of the data for longer than this (s) gets a new ambiguity. With data at a long
    /// interval, two and a half intervals are allowed instead.
    pub max_outage_s: f64,

    /// Cycle slip: jump of the geometry-free phase of one receiver between epochs, m.
    pub slip_gf: f64,
    /// Cycle slip: misfit of an epoch-to-epoch phase change against the motion all satellites agree on, m.
    pub slip_td: f64,
    /// Outliers: an innovation larger than this many standard deviations is rejected (code) or resets
    /// the ambiguity (phase).
    pub outlier_k: f64,
    /// A code residual larger than this (m) is rejected whatever its standard deviation.
    pub code_gross: f64,

    /// The base's observations are interpolated to the rover's epochs between base epochs at most this
    /// far apart, s.
    pub base_max_gap: f64,
    /// A base epoch this close to the rover epoch is used as it is, s.
    pub base_tol: f64,
    /// Beyond this distance between the receivers (km) a warning says the model is being stretched.
    pub max_baseline_km: f64,
    /// Fewest satellites for a phase solution.
    pub min_sats: usize,
    /// Compute Doppler velocities (needs Doppler in the rover file).
    pub velocity: bool,

    /// How a position at an event is interpolated between epochs.
    pub event_interp: EventInterp,
    /// How many missing epochs next to an event may be skipped (0: an event beside a gap gets no
    /// position).
    pub event_max_missing: usize,

    /// Realistic σ of a fixed epoch: the formal covariance is scaled by this factor squared ...
    pub sigma_fixed_scale: f64,
    /// ... and this horizontal floor (m) ...
    pub sigma_fixed_h: f64,
    /// ... and vertical floor (m) are added, with ...
    pub sigma_fixed_v: f64,
    /// ... this part per million of the baseline length.
    pub sigma_ppm: f64,
    /// Realistic σ of a float epoch: the formal covariance scaled by this factor squared.
    pub sigma_float_scale: f64,
    /// Realistic σ of a code-only epoch: the range error assumed, m.
    pub sigma_single: f64,
}

impl Default for PpkOptions {
    fn default() -> Self {
        Self {
            systems: "GRECJ".into(),
            exclude: Vec::new(),
            elmask: 12.0,
            max_bands: 3,
            direction: Direction::Combined,
            keep_passes: false,
            err_phase: 0.003,
            err_phase_el: 0.003,
            code_ratio: 100.0,
            code_corr_s: 5.0,
            efact_glonass: 1.5,
            trop_model: "niell".into(),
            humidity: 0.5,
            pcv_azimuth: true,
            ar: ArMode::Continuous,
            ar_ratio: 3.0,
            ar_elmask: 15.0,
            ar_min_lock_s: 2.0,
            ar_min_dd: 5,
            ar_max_res: 0.02,
            ar_confirm_s: 1.0,
            ar_max_shift_sigma: 5.0,
            ar_partial_rounds: 3,
            ar_glonass: false,
            amb_sigma0: 30.0,
            amb_rw: 1e-4,
            pos_sigma0: 30.0,
            max_outage_s: 1.0,
            slip_gf: 0.05,
            slip_td: 0.05,
            outlier_k: 5.0,
            code_gross: 30.0,
            base_max_gap: 2.5,
            base_tol: 5e-4,
            max_baseline_km: 20.0,
            min_sats: 5,
            velocity: true,
            event_interp: EventInterp::Cubic4,
            event_max_missing: 0,
            sigma_fixed_scale: 2.0,
            sigma_fixed_h: 0.005,
            sigma_fixed_v: 0.010,
            sigma_ppm: 1.0,
            sigma_float_scale: 3.0,
            sigma_single: 3.0,
        }
    }
}

impl PpkOptions {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_overrides() {
        let o = PpkOptions::default()
            .with_json(r#"{"elmask": 20, "direction": "forward", "event_interp": "linear2"}"#)
            .unwrap();
        assert_eq!(o.elmask, 20.0);
        assert_eq!(o.direction, Direction::Forward);
        assert_eq!(o.event_interp, EventInterp::Linear2);
        assert!(PpkOptions::default().with_json(r#"{"nope": 1}"#).is_err());
    }
}
