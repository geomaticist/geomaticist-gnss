//! Precise point positioning: one station, no base, the analysis centre's products.
//!
//! [`ppp_solve`] runs the static solution (one position for the session) or, with `kinematic=True`, one
//! position per epoch for a moving receiver, with camera events and the `.pos` texts as [`super::ppk`] gives
//! them. Products are file paths, read the way `gnss ppp` reads them.

use std::path::PathBuf;
use std::time::Instant;

use geomaticist_gnss_engine::dump::json_text;
use geomaticist_gnss_engine::models::ecef_to_llh;
use geomaticist_gnss_engine::orbits::{OrbitSource, Orbits};
use geomaticist_gnss_engine::ppk::pos::{PosHeader, epochs_pos, events_csv, events_pos};
use geomaticist_gnss_engine::ppk::{ArMode, PpkError, PpkOptions};
use geomaticist_gnss_engine::tropo::TropoInputs;
use geomaticist_gnss_ppp::kinematic::{KinematicInput, KinematicOptions};
use geomaticist_gnss_ppp::{PppOptions, realistic_cov};
use geomaticist_gnss_products::{Precise, read_blq, read_erp, read_osb};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use serde_json::{Value, json};

use super::ppk::{event_of, track_json};
use super::{TropoOwned, antex_of, bad, station_arc, tropo_owned, xyz_of};

fn paths(v: Option<Vec<String>>) -> Vec<PathBuf> {
    v.unwrap_or_default()
        .into_iter()
        .map(PathBuf::from)
        .collect()
}

fn runtime(e: impl std::fmt::Display) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// `PppOptions` with the fields of a JSON object replaced; an unknown field is an error.
fn ppp_options(json: &str) -> PyResult<PppOptions> {
    let over: serde_json::Map<String, Value> = serde_json::from_str(json).map_err(bad)?;
    let mut v = serde_json::to_value(PppOptions::default()).map_err(bad)?;
    for (k, val) in over {
        if k == "blq" {
            return Err(bad(
                "options: give the BLQ file with blq= (and set \"otl\": true)",
            ));
        }
        if v.get(&k).is_none() {
            return Err(bad(format!("options: unknown field {k:?}")));
        }
        v[k] = val;
    }
    serde_json::from_value(v).map_err(bad)
}

/// PPP of one station -> JSON.
///
/// - `station`: a station dict (`{"handle": read_obs(...), "name": .., "antenna": .., "receiver": ..}`,
///   as for `ppk_solve`). The antenna type selects the ANTEX model.
/// - `sp3`, `clk`: precise orbit and clock files (one per day the session spans). Required.
/// - `antex`: ANTEX file (e.g. `igs20.atx`) for satellite and receiver antennas.
/// - `osb`: Bias-SINEX files (code biases; phase biases too for PPP-AR, `"ar": true`).
/// - `erp`: IGS ERP files (pole tide). `nav`: broadcast files, for the GLONASS channels only.
/// - `blq`, `blq_station`: ocean loading coefficients (Onsala BLQ), applied with `"otl": true`.
/// - `xyz0`: start position, ECEF metres (default: a code solution).
/// - `options_json`: `PppOptions` fields (`systems` default `"GRE"`, `elmask`, `ar`, `otl`, `trop_model`, ...).
/// - `tropo`: `{"gpt3": path, "orography": path, "vmf3": [paths]}` for `trop_model` gpt3 / vmf3.
/// - `kinematic`: one position per epoch instead of one per session. Then also `kinematic_options_json`
///   (`KinematicOptions` fields except `ppp`), `events` (a list of `(gps_seconds, label, None or (north,
///   east, down))`, e.g. from `read_mrk` or `rinex_events`), and `outputs`, comma-separated: `epochs`, `pos`,
///   `events_pos`, `events_csv`.
///
/// Static returns `result` (the position `xyz`, `hours`, `fixed`, `sigma0`, `rms_phase`, `ztd`, `ar`, ...),
/// `llh_deg`, `cov_realistic_enu` (ENU, m², scaled for the session length) and `seconds`. Kinematic returns
/// `summary`, `warnings`, `notes`, `events`, `seconds` and what `outputs` asks for.
#[pyfunction]
#[pyo3(signature = (station, sp3, clk=None, antex=None, osb=None, erp=None, nav=None, blq=None, blq_station=None, xyz0=None, options_json="{}", tropo=None, kinematic=false, kinematic_options_json="{}", events=None, outputs="epochs"))]
#[allow(
    clippy::too_many_arguments,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::fn_params_excessive_bools
)]
pub fn ppp_solve(
    py: Python<'_>,
    station: &Bound<'_, PyDict>,
    sp3: Vec<String>,
    clk: Option<Vec<String>>,
    antex: Option<String>,
    osb: Option<Vec<String>>,
    erp: Option<Vec<String>>,
    nav: Option<Vec<String>>,
    blq: Option<String>,
    blq_station: Option<String>,
    xyz0: Option<Bound<'_, PyAny>>,
    options_json: &str,
    tropo: Option<Bound<'_, PyDict>>,
    kinematic: bool,
    kinematic_options_json: &str,
    events: Option<Bound<'_, PyList>>,
    outputs: &str,
) -> PyResult<String> {
    let t_all = Instant::now();
    if sp3.is_empty() {
        return Err(bad("sp3: at least one precise orbit file"));
    }
    let st = station_arc(station)?;
    let opts = ppp_options(options_json)?;
    let x0 = match xyz0 {
        Some(v) if !v.is_none() => Some(xyz_of(&v)?),
        _ => None,
    };
    let nav = geomaticist_gnss_rinex::nav::read_nav(&paths(nav));
    let fcn = nav.glo_fcn();
    let precise = Precise::open(&paths(Some(sp3)), &paths(clk)).map_err(runtime)?;
    let orb = Orbits::new(vec![OrbitSource::Precise(precise)]);
    let mut atx = antex_of(antex)?;
    let osb = read_osb(&paths(osb)).map_err(runtime)?;
    let erp_files = paths(erp);
    let erp = if erp_files.is_empty() {
        None
    } else {
        Some(read_erp(&erp_files).map_err(runtime)?)
    };
    let blq = match blq {
        None => None,
        Some(p) => {
            let text = std::fs::read(&p).map_err(|e| runtime(format!("{p}: {e}")))?;
            let map = read_blq(&text);
            let key = blq_station
                .unwrap_or_else(|| st.name.clone())
                .to_uppercase();
            let Some(c) = map.get(&key) else {
                return Err(runtime(format!(
                    "{p}: no station {key:?} in the BLQ file (it has {:?})",
                    map.keys().collect::<Vec<_>>()
                )));
            };
            Some(*c)
        }
    };
    let TropoOwned { gpt3, orog, vmf3 } = tropo_owned(tropo.as_ref())?;
    let mut ev = Vec::new();
    if let Some(list) = &events {
        for it in list.iter() {
            ev.push(event_of(&it)?);
        }
    }
    let ko = if kinematic {
        let mut k = KinematicOptions::default()
            .with_json(kinematic_options_json)
            .map_err(bad)?;
        k.ppp = opts.clone();
        Some(k)
    } else {
        None
    };
    let read_s = t_all.elapsed().as_secs_f64();

    let tin = |vmf3| TropoInputs {
        gpt3: gpt3.as_deref(),
        orography: orog.as_deref().map(Vec::as_slice),
        vmf3,
    };
    let Some(ko) = ko else {
        // static
        let t_solve = Instant::now();
        let (res, _stages) = py
            .detach(|| {
                geomaticist_gnss_ppp::solve(
                    &st,
                    &orb,
                    atx.as_mut(),
                    &opts,
                    &fcn,
                    &osb,
                    erp.as_ref(),
                    x0,
                    blq.as_ref(),
                    &tin(&vmf3),
                )
            })
            .map_err(runtime)?;
        let enu = res.enu_cov();
        let real = realistic_cov(&enu, res.hours, res.fixed);
        let (lat, lon, h) = ecef_to_llh(res.xyz);
        let v = json!({
            "result": res.to_json(),
            "llh_deg": [lat.to_degrees(), lon.to_degrees(), h],
            "cov_realistic_enu": real,
            "seconds": {"read": read_s, "solve": t_solve.elapsed().as_secs_f64()},
        });
        return Ok(json_text(&v));
    };

    // kinematic
    let t_solve = Instant::now();
    let res = py
        .detach(|| {
            let tr = tin(&vmf3);
            let mut input = KinematicInput {
                obs: &st,
                orbits: &orb,
                antex: atx.as_mut(),
                fcn: &fcn,
                osb: &osb,
                erp: erp.as_ref(),
                blq: blq.as_ref(),
                tropo: &tr,
                x0,
                events: &ev,
                occupations: &[],
            };
            geomaticist_gnss_ppp::kinematic::solve(&mut input, &ko, &mut |_| true)
        })
        .map_err(|e| match e {
            PpkError::Cancelled => runtime("cancelled"),
            PpkError::Input(s) => runtime(s),
        })?;
    let want = |k: &str| outputs.split(',').any(|w| w.trim() == k);
    let mut v = json!({
        "summary": res.summary,
        "warnings": res.warnings,
        "notes": res.notes,
        "events": res.events,
        "seconds": {"read": read_s, "solve": t_solve.elapsed().as_secs_f64()},
    });
    if want("epochs") {
        v["epochs"] = track_json(&res.epochs);
    }
    let mut h = PosHeader::new(&PpkOptions {
        systems: opts.systems.clone(),
        elmask: opts.elmask,
        direction: ko.direction,
        ar: ArMode::Off,
        trop_model: opts.trop_model.clone(),
        ..PpkOptions::default()
    });
    "PPP Kinematic".clone_into(&mut h.mode);
    if want("pos") {
        v["pos"] = json!(epochs_pos(&res, &h));
    }
    if want("events_pos") {
        let mut he = h.clone();
        he.time_decimals = 6;
        v["events_pos"] = json!(events_pos(&res, &he));
    }
    if want("events_csv") {
        v["events_csv"] = json!(events_csv(&res));
    }
    // NaN (an epoch or event without a position) is written as null
    serde_json::to_string(&v).map_err(bad)
}
