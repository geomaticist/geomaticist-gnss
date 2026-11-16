//! Kinematic processing and events for Python: [`ppk_solve`] and the event readers and writers.
//!
//! The functions take what [`super::static_solve`] takes (station dicts or `Station` handles from `read_obs`,
//! the list of orbit sources, the options as JSON) and return JSON text; the solve runs without the GIL.
//!
//! ```python
//! import json, geomaticist_gnss_py as g
//!
//! def station(h):
//!     return {"handle": h, "name": h.name, "antenna": h.antenna, "receiver": h.receiver,
//!             "approx": h.approx, "delta_hen": h.delta_hen}
//!
//! base = g.read_obs("base.obs")
//! rover = g.read_obs("rover.obs", kinds="CLSD")          # with Doppler
//! nav = g.read_nav(["base.nav"])
//! mrk = json.loads(g.read_mrk("flight_Timestamp.MRK", rover_first_epoch_gps_seconds))
//! out = json.loads(g.ppk_solve(
//!     station(base), [station(rover)], [("broadcast_handle", nav)], base_xyz,
//!     options_json='{"event_interp": "cubic4"}', antex="igs20.atx", fcn=nav_fcn,
//!     events=mrk["events"], outputs="epochs,pos,events_csv"))
//! print(out["summary"]["fix_percent"], len(out["events"]))
//! ```

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use geomaticist_gnss_core::StationObs;
use geomaticist_gnss_core::event::{Event, EventOffset, Occupation};
use geomaticist_gnss_core::time::fmt;
use geomaticist_gnss_engine::orbits::Orbits;
use geomaticist_gnss_engine::ppk::pos::{PosHeader, epochs_pos, events_csv, events_pos, track_pos};
use geomaticist_gnss_engine::ppk::{
    self as engine, BaseStation, EpochSolution, PpkError, PpkInput, PpkOptions, Stage,
};
use geomaticist_gnss_engine::tropo::TropoInputs;
use geomaticist_gnss_events::{mrk, rinex as evrinex};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use serde_json::{Value, json};

use super::{TropoOwned, antex_of, bad, fcn_of, orbit_sources, station_arc, tropo_owned, xyz_of};

/// One event from Python: `(t, label, offset)` with `offset` `None` or `(north, east, down)` metres, or
/// just a number (GPS seconds). Lists work as tuples, so the `events` of [`read_mrk`] / [`rinex_events`]
/// (JSON arrays) pass straight through.
pub(crate) fn event_of(o: &Bound<'_, PyAny>) -> PyResult<Event> {
    if let Ok(t) = o.extract::<f64>() {
        return Ok(Event::at(t));
    }
    let msg = || bad("an event is (gps_seconds, label, None or (north, east, down))");
    let parts: Vec<Bound<'_, PyAny>> = o.extract().map_err(|_| msg())?;
    let [t, label, off] = parts.as_slice() else {
        return Err(msg());
    };
    let t: f64 = t.extract().map_err(|_| msg())?;
    let label: String = label.extract().map_err(|_| msg())?;
    let off: Option<[f64; 3]> = if off.is_none() {
        None
    } else {
        Some(off.extract().map_err(|_| msg())?)
    };
    Ok(Event {
        t,
        label,
        offset: off.map(|[north, east, down]| EventOffset { north, east, down }),
    })
}

/// One occupation from Python: `(name, t0, t1)` (GPS seconds), as a tuple or a list.
fn occupation_of(o: &Bound<'_, PyAny>) -> PyResult<Occupation> {
    let msg = || bad("an occupation is (name, gps_t0, gps_t1)");
    let parts: Vec<Bound<'_, PyAny>> = o.extract().map_err(|_| msg())?;
    let [name, t0, t1] = parts.as_slice() else {
        return Err(msg());
    };
    Ok(Occupation {
        name: name.extract().map_err(|_| msg())?,
        t0: t0.extract().map_err(|_| msg())?,
        t1: t1.extract().map_err(|_| msg())?,
    })
}

fn events_json(ev: &[Event]) -> Value {
    Value::Array(
        ev.iter()
            .map(|e| json!([e.t, e.label, e.offset.map(|o| [o.north, o.east, o.down])]))
            .collect(),
    )
}

/// A track as columns (shorter than one object per epoch): `t`, `t_tag`, `q` (1 fix, 2 float, 5 single,
/// 0 none), `lat`, `lon`, `h`, `x`, `y`, `z`, `sd_e`, `sd_n`, `sd_u`, `ns`, `ratio`, `age`, `ve`, `vn`, `vu`.
pub(crate) fn track_json(tr: &[EpochSolution]) -> Value {
    let f = |g: &dyn Fn(&EpochSolution) -> f64| -> Value {
        Value::Array(tr.iter().map(|e| json!(g(e))).collect())
    };
    json!({
        "t": f(&|e| e.t),
        "t_tag": f(&|e| e.t_tag),
        "q": tr.iter().map(|e| e.status.q()).collect::<Vec<_>>(),
        "lat": f(&|e| e.llh[0]),
        "lon": f(&|e| e.llh[1]),
        "h": f(&|e| e.llh[2]),
        "x": f(&|e| e.xyz[0]),
        "y": f(&|e| e.xyz[1]),
        "z": f(&|e| e.xyz[2]),
        "sd_e": f(&|e| e.sigma_enu()[0]),
        "sd_n": f(&|e| e.sigma_enu()[1]),
        "sd_u": f(&|e| e.sigma_enu()[2]),
        "ns": tr.iter().map(|e| e.n_sats).collect::<Vec<_>>(),
        "ratio": f(&|e| f64::from(e.ratio)),
        "age": f(&|e| f64::from(e.age)),
        "ve": f(&|e| f64::from(e.vel_enu[0])),
        "vn": f(&|e| f64::from(e.vel_enu[1])),
        "vu": f(&|e| f64::from(e.vel_enu[2])),
    })
}

fn stage_name(s: Stage) -> &'static str {
    match s {
        Stage::Prepare => "prepare",
        Stage::Forward => "forward",
        Stage::Backward => "backward",
        Stage::Combine => "combine",
        Stage::Events => "events",
    }
}

/// Kinematic processing of a rover against a static base.
///
/// - `base`: a station dict (as for `static_solve`); `rover`: a station dict or a list of them (the files
///   of one receiver session, joined in time order).
/// - `sources`: orbit sources as for `static_solve`.
/// - `base_xyz`: ECEF metres of the base **mark**; `base_hen`: antenna height / east / north above it
///   (default: the base file's header).
/// - `options_json`: fields of `PpkOptions` to change (unknown fields are an error).
/// - `events`: a list of `(gps_seconds, label, None or (north, east, down))`.
/// - `occupations`: a list of `(name, t0, t1)` (GPS seconds); lists work too, so the output of
///   [`rinex_occupations`] passes straight in.
/// - `progress`: called as `progress(stage, done, total)`; returning `False` cancels (`RuntimeError`
///   "cancelled").
/// - `outputs`: which large parts the JSON holds, comma-separated: `epochs` (the track as columns),
///   `passes` (the forward and backward tracks), `pos`, `events_pos`, `events_csv` (the text files).
///
/// Returns JSON: `summary`, `warnings`, `notes`, `base_arp`, `base_marker`, `events`, `occupations`,
/// `seconds`, and what `outputs` asks for.
#[pyfunction]
#[pyo3(signature = (base, rover, sources, base_xyz, options_json="{}", antex=None, fcn=None, base_hen=None, events=None, occupations=None, tropo=None, progress=None, outputs="epochs"))]
#[allow(
    clippy::too_many_arguments,
    clippy::needless_pass_by_value,
    clippy::too_many_lines
)]
pub fn ppk_solve(
    py: Python<'_>,
    base: &Bound<'_, PyDict>,
    rover: &Bound<'_, PyAny>,
    sources: &Bound<'_, PyList>,
    base_xyz: &Bound<'_, PyAny>,
    options_json: &str,
    antex: Option<String>,
    fcn: Option<BTreeMap<String, i32>>,
    base_hen: Option<Bound<'_, PyAny>>,
    events: Option<Bound<'_, PyList>>,
    occupations: Option<Bound<'_, PyList>>,
    tropo: Option<Bound<'_, PyDict>>,
    progress: Option<Py<PyAny>>,
    outputs: &str,
) -> PyResult<String> {
    let t_conv = Instant::now();
    let base_st = station_arc(base)?;
    let rover_st: Arc<StationObs> = if let Ok(d) = rover.cast::<PyDict>() {
        station_arc(d)?
    } else {
        let list = rover
            .cast::<PyList>()
            .map_err(|_| bad("rover: a station dict or a list of station dicts"))?;
        let mut parts = Vec::new();
        for it in list.iter() {
            let d = it
                .cast::<PyDict>()
                .map_err(|_| bad("rover: not a station dict"))?;
            parts.push(station_arc(d)?.as_ref().clone());
        }
        Arc::new(StationObs::concat(&parts).ok_or_else(|| bad("rover: an empty list"))?)
    };
    let (srcs, externals) = orbit_sources(sources)?;
    let orb = Orbits::new(srcs);
    let marker_xyz = xyz_of(base_xyz)?;
    let delta_hen = match base_hen {
        Some(v) if !v.is_none() => Some(xyz_of(&v)?),
        _ => None,
    };
    let opts = PpkOptions::default().with_json(options_json).map_err(bad)?;
    let fcn = fcn_of(fcn)?;
    let mut atx = antex_of(antex)?;
    let TropoOwned { gpt3, orog, vmf3 } = tropo_owned(tropo.as_ref())?;
    let mut ev: Vec<Event> = Vec::new();
    if let Some(list) = &events {
        for it in list.iter() {
            ev.push(event_of(&it)?);
        }
    }
    let mut occ: Vec<Occupation> = Vec::new();
    if let Some(list) = &occupations {
        for it in list.iter() {
            occ.push(occupation_of(&it)?);
        }
    }
    let want = |k: &str| outputs.split(',').any(|w| w.trim() == k);
    let mut o = opts.clone();
    if want("passes") {
        o.keep_passes = true;
    }
    let convert_s = t_conv.elapsed().as_secs_f64();

    let out = py.detach(move || {
        let t_solve = Instant::now();
        let tin = TropoInputs {
            gpt3: gpt3.as_deref(),
            orography: orog.as_deref().map(Vec::as_slice),
            vmf3: &vmf3,
        };
        let mut input = PpkInput {
            base: &base_st,
            rover: &rover_st,
            orbits: &orb,
            fcn: &fcn,
            base_station: BaseStation {
                marker_xyz,
                delta_hen,
            },
            antex: atx.as_mut(),
            tropo: &tin,
            events: &ev,
            occupations: &occ,
        };
        let mut cb_error: Option<String> = None;
        let mut cb = |p: engine::Progress| -> bool {
            let Some(f) = &progress else { return true };
            Python::attach(|py| {
                match f.call1(py, (stage_name(p.stage), p.done, p.total)) {
                    // only an explicit False cancels
                    Ok(v) => v.extract::<bool>(py).unwrap_or(true) || v.is_none(py),
                    Err(e) => {
                        cb_error = Some(e.to_string());
                        false
                    }
                }
            })
        };
        let r = engine::solve(&mut input, &o, &mut cb);
        let ext_err = externals
            .iter()
            .find_map(|e| e.error.lock().ok().and_then(|g| g.clone()));
        (r, t_solve.elapsed().as_secs_f64(), ext_err, cb_error, o)
    });
    let (r, solve_s, ext_err, cb_error, o) = out;
    if let Some(e) = ext_err {
        return Err(PyRuntimeError::new_err(format!("orbit source: {e}")));
    }
    if let Some(e) = cb_error {
        return Err(PyRuntimeError::new_err(format!("progress callback: {e}")));
    }
    let res = r.map_err(|e| match e {
        PpkError::Cancelled => PyRuntimeError::new_err("cancelled"),
        PpkError::Input(s) => PyRuntimeError::new_err(s),
    })?;
    let mut v = json!({
        "summary": res.summary,
        "warnings": res.warnings,
        "notes": res.notes,
        "base_arp": res.base_arp,
        "base_marker": res.base_marker,
        "events": res.events,
        "occupations": res.occupations,
        "seconds": {"convert": convert_s, "solve": solve_s},
        "bridge": super::VERSION,
    });
    if want("epochs") {
        v["epochs"] = track_json(&res.epochs);
    }
    if want("passes") {
        v["forward"] = track_json(&res.forward);
        v["backward"] = track_json(&res.backward);
    }
    let h = PosHeader::new(&o);
    if want("pos") {
        v["pos"] = json!(epochs_pos(&res, &h));
        if want("passes") {
            v["forward_pos"] = json!(track_pos(&res, &res.forward, &h));
            v["backward_pos"] = json!(track_pos(&res, &res.backward, &h));
        }
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

/// Reads a DJI `Timestamp.MRK` file -> JSON `{"records": [...], "events": [[t, label, [north, east, down]],
/// ...], "quality": {flag: count}, "warnings": [...]}`. `near` (GPS seconds, the rover's first epoch) gives
/// lines without a GPS week their week.
#[pyfunction]
#[pyo3(signature = (path, near=None))]
pub fn read_mrk(path: &str, near: Option<f64>) -> PyResult<String> {
    let mut m = mrk::read(Path::new(path)).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    if let Some(t) = near {
        m.resolve_weeks(t);
    }
    let q: BTreeMap<String, usize> = m
        .quality_counts()
        .into_iter()
        .map(|(k, n)| (k.to_string(), n))
        .collect();
    serde_json::to_string(&json!({
        "records": m.records,
        "events": events_json(&m.events()),
        "quality": q,
        "warnings": m.warnings,
    }))
    .map_err(bad)
}

/// The event records (epoch flag 5) of a RINEX observation file -> JSON `{"events": [[t, label, null], ...],
/// "warnings": [...]}`.
#[pyfunction]
pub fn rinex_events(py: Python<'_>, path: &str) -> PyResult<String> {
    let r = py
        .detach(|| evrinex::read_events(Path::new(path)))
        .map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    serde_json::to_string(&json!({"events": events_json(&r.events), "warnings": r.warnings}))
        .map_err(bad)
}

/// The stop-and-go occupations (epoch flags 2 / 3) of a RINEX observation file -> JSON `[[name, t0, t1], ...]`.
#[pyfunction]
pub fn rinex_occupations(py: Python<'_>, path: &str) -> PyResult<String> {
    let r = py
        .detach(|| evrinex::read_occupations(Path::new(path)))
        .map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    let v: Vec<Value> = r.iter().map(|o| json!([o.name, o.t0, o.t1])).collect();
    serde_json::to_string(&v).map_err(bad)
}

/// Checks the epoch times of a RINEX observation file for known recorder glitches -> JSON `{"epochs",
/// "before_1990", "displaced", "not_later", "first", "last", "warnings"}`.
#[pyfunction]
pub fn check_times(py: Python<'_>, path: &str) -> PyResult<String> {
    let c = py
        .detach(|| {
            let data = geomaticist_gnss_rinex::open::read_plain(Path::new(path))
                .map_err(|e| e.to_string())?;
            evrinex::check_times(&data).map_err(|e| e.to_string())
        })
        .map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    serde_json::to_string(&json!({
        "epochs": c.epochs,
        "before_1990": c.before_1990,
        "displaced": c.displaced,
        "not_later": c.not_later,
        "first": c.span.map(|s| fmt(s.0, 7)),
        "last": c.span.map(|s| fmt(s.1, 7)),
        "warnings": c.warnings,
    }))
    .map_err(bad)
}

/// Writes a copy of a RINEX observation file with an event record (flag 5) per event -> JSON `{"path",
/// "inserted", "before", "after", "existing", "warnings"}`.
#[pyfunction]
pub fn insert_events(path: &str, events: &Bound<'_, PyList>, out_path: &str) -> PyResult<String> {
    let mut ev = Vec::new();
    for it in events.iter() {
        ev.push(event_of(&it)?);
    }
    let run = || -> Result<evrinex::Inserted, String> {
        let data =
            geomaticist_gnss_rinex::open::read_plain(Path::new(path)).map_err(|e| e.to_string())?;
        let ins = evrinex::insert_events(&data, &ev).map_err(|e| e.to_string())?;
        std::fs::write(out_path, &ins.data).map_err(|e| format!("{out_path}: {e}"))?;
        Ok(ins)
    };
    let ins = run().map_err(PyRuntimeError::new_err)?;
    serde_json::to_string(&json!({
        "path": out_path,
        "inserted": ins.inserted,
        "before": ins.before,
        "after": ins.after,
        "existing": ins.existing,
        "warnings": ins.warnings,
    }))
    .map_err(bad)
}

/// Writes a copy of a RINEX observation file marked as a moving receiver's, so office software that reads the
/// file's own words (Trimble Business Center among them) imports it as kinematic: RINEX 3 / 4 get the header's
/// `MARKER TYPE` set to `marker_type` (e.g. `AIRBORNE`), and a "start moving antenna" record (flag 2) is put
/// before the first epoch when none is there -> JSON `{"path", "flag_added", "blank_removed"}`.
#[pyfunction]
#[pyo3(signature = (path, out_path, marker_type="AIRBORNE"))]
pub fn mark_kinematic(
    py: Python<'_>,
    path: &str,
    out_path: &str,
    marker_type: &str,
) -> PyResult<String> {
    let k = py
        .detach(|| -> Result<evrinex::Kinematic, String> {
            let data = geomaticist_gnss_rinex::open::read_plain(Path::new(path))
                .map_err(|e| e.to_string())?;
            let k = evrinex::mark_kinematic(&data, marker_type).map_err(|e| e.to_string())?;
            std::fs::write(out_path, &k.data).map_err(|e| format!("{out_path}: {e}"))?;
            Ok(k)
        })
        .map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    serde_json::to_string(&json!({
        "path": out_path,
        "flag_added": k.flag_added,
        "blank_removed": k.blank_removed,
    }))
    .map_err(bad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geomaticist_gnss_engine::ppk::FixStatus;

    #[test]
    fn tracks_and_events_as_json() {
        let mut e = EpochSolution::none(100.0);
        e.status = FixStatus::Fixed;
        e.set_xyz([4_395_470.6, 617_743.1, 4_565_607.2]);
        e.cov_enu = [[1e-4, 0.0, 0.0], [0.0, 4e-4, 0.0], [0.0, 0.0, 9e-4]];
        let v = track_json(&[e, EpochSolution::none(101.0)]);
        assert_eq!(v["q"], json!([1, 0]));
        assert_eq!(v["t_tag"], json!([100.0, 101.0]));
        assert!((v["sd_n"][0].as_f64().unwrap() - 0.02).abs() < 1e-12);
        // NaN becomes null in the text
        assert!(serde_json::to_string(&v).unwrap().contains("null"));
        let ev = [Event {
            t: 5.0,
            label: "a".into(),
            offset: Some(EventOffset {
                north: 1.0,
                east: 2.0,
                down: 3.0,
            }),
        }];
        assert_eq!(events_json(&ev), json!([[5.0, "a", [1.0, 2.0, 3.0]]]));
        assert_eq!(stage_name(Stage::Backward), "backward");
    }
}
