//! `geomaticist-gnss-py`: Python bindings (pyo3) for the Geomaticist GNSS engine.
//!
//! Built by maturin as one abi3 wheel; the module is `geomaticist_gnss_py`. Every heavy call (reads,
//! decoding, solves) releases the GIL, so Python threads can run several at once. Results come back as
//! JSON text (`json.loads`).
//!
//! - Reading: [`read_obs`] and [`read_nav`] read RINEX with the engine's own reader into handles (`Station`,
//!   `Nav`) that keep the data in Rust; `read_sp3` reads precise orbits and clocks (`Precise`).
//! - Static baselines: [`static_solve`] takes two station dicts (`{"handle": Station, "name": ..,
//!   "antenna": .., "receiver": .., "approx": .., "delta_hen": ..}`; or the observations packed as
//!   little-endian `bytes`), a list of orbit sources (`("broadcast_handle", Nav)`, `("precise_handle",
//!   Precise)`, or any Python object with `has` / `position` methods), the options as JSON, and optional ANTEX
//!   and troposphere grid paths.
//! - Kinematic (PPK), events and stop-and-go: [`ppk`]. Precise point positioning, static and kinematic: [`ppp`].
//! - Raw receiver logs to RINEX: [`decode_raw`]. Simulated data for tests and examples: [`synth`].
//!
//! This is the only crate allowed to depend on pyo3 (`cargo xtask check-deps`, rule "py"), and no crate may
//! depend on it.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;

use geomaticist_gnss_core::obs::{Signal, StationObs};
use geomaticist_gnss_core::{GloEph, KeplerEph, NavData, ObsCode, Sat, SignalKey};
use geomaticist_gnss_engine::dump::{json_text, write_stages};
use geomaticist_gnss_engine::orbits::{Broadcast, ExternalOrbit, OrbitSource, Orbits};
use geomaticist_gnss_engine::tropo::{TropoInputs, Vmf3Grid};
use geomaticist_gnss_engine::{Options, solve};
use geomaticist_gnss_products::grids::{Gpt3Grid, load_orography, load_vmf3};
use geomaticist_gnss_products::{Antex, Precise};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyTuple};
use serde_json::{Value, json};

pub mod ppk;
pub mod ppp;
pub mod synth;

/// Version string shown in the processing log / report ("which engine ran").
pub const VERSION: &str = concat!("geomaticist-gnss-py ", env!("CARGO_PKG_VERSION"));

// ---- pure helpers (tested without Python) ---------------------------------------------------------------

/// Little-endian f64s from bytes (numpy `astype('<f8').tobytes()`).
pub fn f64s(b: &[u8]) -> Result<Vec<f64>, String> {
    let (chunks, rest) = b.as_chunks::<8>();
    if !rest.is_empty() {
        return Err(format!("f64 buffer of {} bytes", b.len()));
    }
    Ok(chunks.iter().map(|c| f64::from_le_bytes(*c)).collect())
}

/// Rows of `N` f64s.
pub fn rows<const N: usize>(b: &[u8]) -> Result<Vec<[f64; N]>, String> {
    let v = f64s(b)?;
    let (chunks, rest) = v.as_chunks::<N>();
    if !rest.is_empty() {
        return Err(format!("{} values are not rows of {N}", v.len()));
    }
    Ok(chunks.to_vec())
}

/// One Keplerian record from a row in `KEP_FIELDS` order (the inverse of `KeplerEph::to_array`).
#[must_use]
pub fn kepler(r: &[f64; 26]) -> KeplerEph {
    KeplerEph {
        toc: r[0],
        f0: r[1],
        f1: r[2],
        f2: r[3],
        iode: r[4],
        crs: r[5],
        deln: r[6],
        m0: r[7],
        cuc: r[8],
        e: r[9],
        cus: r[10],
        sqrta: r[11],
        toes: r[12],
        cic: r[13],
        omg0: r[14],
        cis: r[15],
        i0: r[16],
        crc: r[17],
        omg: r[18],
        omgd: r[19],
        idot: r[20],
        week: r[21],
        svh: r[22],
        tgd: r[23],
        toe: r[24],
        acc: r[25],
    }
}

/// One GLONASS record from a row: toe, taun, gamn, pos(3), vel(3), acc(3), svh, fcn.
#[must_use]
pub fn glonass(r: &[f64; 14]) -> GloEph {
    GloEph {
        toe: r[0],
        taun: r[1],
        gamn: r[2],
        pos: [r[3], r[4], r[5]],
        vel: [r[6], r[7], r[8]],
        acc: [r[9], r[10], r[11]],
        svh: r[12],
        fcn: r[13],
    }
}

/// A `(kind, band, attr)` code key -> `SignalKey`.
pub fn signal_key(kind: &str, band: i64, attr: &str) -> Result<SignalKey, String> {
    let mut k = kind.chars();
    let (Some(kind_c), None) = (k.next(), k.next()) else {
        return Err(format!("signal kind {kind:?}"));
    };
    let band = u8::try_from(band).map_err(|_| format!("band {band}"))?;
    let mut a = attr.chars();
    let attr = match (a.next(), a.next()) {
        (None, _) => None,
        (Some(c), None) => Some(c),
        _ => return Err(format!("attribute {attr:?}")),
    };
    Ok(SignalKey {
        kind: kind_c,
        band,
        attr,
    })
}

/// The engine's options from a JSON object of `Options` fields. Unknown fields are an error (a misspelt
/// field would otherwise be ignored silently).
pub fn options_from_json(text: &str) -> Result<Options, String> {
    let over: serde_json::Map<String, Value> =
        serde_json::from_str(text).map_err(|e| format!("options: {e}"))?;
    let mut ov = serde_json::to_value(Options::default()).map_err(|e| e.to_string())?;
    for (k, v) in over {
        if ov.get(&k).is_none() {
            return Err(format!("options: the Rust engine has no field {k:?}"));
        }
        ov[k] = v;
    }
    serde_json::from_value(ov).map_err(|e| format!("options: {e}"))
}

fn vmf3_epoch_of(path: &Path) -> Option<f64> {
    // VMF3_yyyymmdd.Hhh (as the CLI)
    let name = path.file_name()?.to_str()?;
    let rest = name.strip_prefix("VMF3_")?;
    let (date, hour) = rest.split_once(".H")?;
    if date.len() != 8 {
        return None;
    }
    let y: i64 = date[..4].parse().ok()?;
    let mo: i64 = date[4..6].parse().ok()?;
    let d: i64 = date[6..8].parse().ok()?;
    let h: i64 = hour.parse().ok()?;
    let dt = geomaticist_gnss_core::time::DateTime::new(y, mo, d, h, 0, 0)?;
    Some(geomaticist_gnss_core::time::gps_seconds(&dt))
}

// ---- caches of the big product files ----------------------------------------------------------------------

type Cache<T> = OnceLock<Mutex<HashMap<PathBuf, Arc<T>>>>;

fn cached<T>(
    cache: &'static Cache<T>,
    path: &Path,
    load: impl FnOnce(&Path) -> Result<T, String>,
) -> Result<Arc<T>, String> {
    let m = cache.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(v) = m.lock().map_err(|e| e.to_string())?.get(path) {
        return Ok(Arc::clone(v));
    }
    let v = Arc::new(load(path)?);
    m.lock()
        .map_err(|e| e.to_string())?
        .insert(path.to_path_buf(), Arc::clone(&v));
    Ok(v)
}

static ANTEX: Cache<Antex> = OnceLock::new();
static GPT3: Cache<Gpt3Grid> = OnceLock::new();
static OROG: Cache<Vec<f64>> = OnceLock::new();
static VMF3: Cache<Vec<[f64; 6]>> = OnceLock::new();

// ---- Python -> engine inputs ------------------------------------------------------------------------------

fn bad(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(format!("bridge: {e}"))
}

fn sat_of(s: &str) -> PyResult<Sat> {
    s.parse::<Sat>()
        .map_err(|_| bad(format!("satellite {s:?}")))
}

fn item<'py>(d: &Bound<'py, PyDict>, k: &str) -> PyResult<Bound<'py, PyAny>> {
    d.get_item(k)?
        .ok_or_else(|| bad(format!("missing key {k:?}")))
}

fn bytes_of(o: &Bound<'_, PyAny>) -> PyResult<Vec<u8>> {
    Ok(o.cast::<PyBytes>()
        .map_err(|_| bad("expected bytes"))?
        .as_bytes()
        .to_vec())
}

fn xyz_of(o: &Bound<'_, PyAny>) -> PyResult<[f64; 3]> {
    let v: Vec<f64> = o.extract()?;
    <[f64; 3]>::try_from(v.as_slice()).map_err(|_| bad("three numbers expected"))
}

/// The same numbers bit for bit (no tolerance: a changed approx / delta must reach the engine).
fn same_bits(a: &[f64; 3], b: &[f64; 3]) -> bool {
    a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// The station of a packed dict: a [`RustStation`] handle (`"handle"`, read by [`read_obs`]; its name, antenna
/// and receiver set from the dict, so the caller can override them) or the packed arrays (see the module docs).
fn station_arc(d: &Bound<'_, PyDict>) -> PyResult<Arc<StationObs>> {
    if let Some(h) = d.get_item("handle")? {
        let h = h
            .cast::<RustStation>()
            .map_err(|_| bad("handle: not a Station"))?;
        let mut g = h.get().st.lock().map_err(|e| bad(e.to_string()))?;
        let name: String = item(d, "name")?.extract()?;
        let antenna: String = item(d, "antenna")?.extract()?;
        let receiver: String = item(d, "receiver")?.extract()?;
        let approx = xyz_of(&item(d, "approx")?)?;
        let delta_hen = xyz_of(&item(d, "delta_hen")?)?;
        if g.name != name
            || g.antenna != antenna
            || g.receiver != receiver
            || !same_bits(&g.approx, &approx)
            || !same_bits(&g.delta_hen, &delta_hen)
        {
            // in place when the handle holds the only reference (the usual case)
            let s = Arc::make_mut(&mut g);
            s.name = name;
            s.antenna = antenna;
            s.receiver = receiver;
            s.approx = approx;
            s.delta_hen = delta_hen;
        }
        return Ok(Arc::clone(&g));
    }
    Ok(Arc::new(station(d)?))
}

/// A `StationObs` from the packed dict (see the module docs).
fn station(d: &Bound<'_, PyDict>) -> PyResult<StationObs> {
    let name: String = item(d, "name")?.extract()?;
    let path: String = item(d, "path")?.extract()?;
    let mut st = StationObs::new(name, path);
    st.antenna = item(d, "antenna")?.extract()?;
    st.receiver = item(d, "receiver")?.extract()?;
    st.approx = xyz_of(&item(d, "approx")?)?;
    st.delta_hen = xyz_of(&item(d, "delta_hen")?)?;
    st.t = f64s(&bytes_of(&item(d, "t")?)?).map_err(bad)?;
    st.power_fail = bytes_of(&item(d, "power_fail")?)?
        .iter()
        .map(|&b| b != 0)
        .collect();
    if st.power_fail.len() != st.t.len() {
        return Err(bad("power_fail and t differ in length"));
    }
    let ne = st.t.len();
    let mut seen: u32 = 0;
    for entry in item(d, "sats")?.cast::<PyList>()?.iter() {
        let (sat, sigs): (String, Bound<'_, PyList>) = entry.extract()?;
        let sat = sat_of(&sat)?;
        let mut out = BTreeMap::new();
        for s in sigs.iter() {
            let s = s.cast::<PyTuple>()?;
            let code: String = s.get_item(0)?.extract()?;
            let kind: String = s.get_item(1)?.extract()?;
            let band: i64 = s.get_item(2)?.extract()?;
            let attr: String = s.get_item(3)?.extract()?;
            let values = f64s(&bytes_of(&s.get_item(4)?)?).map_err(bad)?;
            let lli_o = s.get_item(5)?;
            let lli = if lli_o.is_none() {
                None
            } else {
                Some(
                    bytes_of(&lli_o)?
                        .iter()
                        .map(|&b| b.cast_signed())
                        .collect::<Vec<i8>>(),
                )
            };
            if values.len() != ne || lli.as_ref().is_some_and(|l| l.len() != ne) {
                return Err(bad(format!("{sat} {code}: array length differs from t")));
            }
            let oc: ObsCode = code.parse().map_err(bad)?;
            out.insert(
                oc,
                Signal {
                    key: signal_key(&kind, band, &attr).map_err(bad)?,
                    values,
                    lli,
                    first_seen: seen,
                },
            );
            seen += 1;
        }
        st.sats.insert(sat, out);
    }
    Ok(st)
}

/// A Python orbit object (anything with `has` / `position` / `center_of_mass`), called back with the GIL.
/// A Python exception in a call is kept and raised after the solve (the call answers "no orbit").
#[derive(Debug)]
struct PyOrbit {
    obj: Py<PyAny>,
    com: bool,
    error: Mutex<Option<String>>,
}

impl PyOrbit {
    fn fail(&self, e: &PyErr) {
        if let Ok(mut g) = self.error.lock()
            && g.is_none()
        {
            *g = Some(e.to_string());
        }
    }
}

/// Position, clock, ok per time (`ExternalOrbit::position`).
type Positions = (Vec<[f64; 3]>, Vec<f64>, Vec<bool>);

/// The engine's handle on a [`PyOrbit`]; the bindings keep another to read the error after the solve.
#[derive(Debug)]
struct SharedOrbit(Arc<PyOrbit>);

impl ExternalOrbit for SharedOrbit {
    fn has(&self, sat: Sat) -> bool {
        self.0.has(sat)
    }

    fn center_of_mass(&self) -> bool {
        self.0.com
    }

    fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        self.0.position(sat, t)
    }
}

impl PyOrbit {
    fn has(&self, sat: Sat) -> bool {
        Python::attach(|py| {
            self.obj
                .bind(py)
                .call_method1("has", (sat.to_string(),))
                .and_then(|r| r.extract::<bool>())
                .unwrap_or_else(|e| {
                    self.fail(&e);
                    false
                })
        })
    }

    fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        let n = t.len();
        let none = || (vec![[0.0; 3]; n], vec![0.0; n], vec![false; n]);
        Python::attach(|py| {
            let r = (|| -> PyResult<Positions> {
                let out = self
                    .obj
                    .bind(py)
                    .call_method1("position", (sat.to_string(), t.to_vec()))?;
                let tup = out.cast::<PyTuple>()?;
                let xyz: Vec<Vec<f64>> = tup.get_item(0)?.call_method0("tolist")?.extract()?;
                let clk: Vec<f64> = tup.get_item(1)?.call_method0("tolist")?.extract()?;
                let ok: Vec<bool> = tup.get_item(2)?.call_method0("tolist")?.extract()?;
                let xyz = xyz
                    .iter()
                    .map(|r| <[f64; 3]>::try_from(r.as_slice()).map_err(|_| bad("position rows")))
                    .collect::<PyResult<Vec<_>>>()?;
                if xyz.len() != n || clk.len() != n || ok.len() != n {
                    return Err(bad("position: lengths differ from t"));
                }
                Ok((xyz, clk, ok))
            })();
            r.unwrap_or_else(|e| {
                self.fail(&e);
                none()
            })
        })
    }
}

fn broadcast(src: &Bound<'_, PyTuple>) -> PyResult<OrbitSource> {
    let mut nav = NavData::default();
    for (k, v) in src.get_item(1)?.cast::<PyDict>()?.iter() {
        let sat = sat_of(&k.extract::<String>()?)?;
        let recs = rows::<26>(&bytes_of(&v)?).map_err(bad)?;
        nav.kep.insert(sat, recs.iter().map(kepler).collect());
    }
    for (k, v) in src.get_item(2)?.cast::<PyDict>()?.iter() {
        let sat = sat_of(&k.extract::<String>()?)?;
        let recs = rows::<14>(&bytes_of(&v)?).map_err(bad)?;
        nav.glo.insert(sat, recs.iter().map(glonass).collect());
    }
    Ok(OrbitSource::Broadcast(Broadcast::new(nav)))
}

fn precise(src: &Bound<'_, PyTuple>) -> PyResult<OrbitSource> {
    let mut p = Precise {
        t: f64s(&bytes_of(&src.get_item(1)?)?).map_err(bad)?,
        ..Precise::default()
    };
    for (k, v) in src.get_item(2)?.cast::<PyDict>()?.iter() {
        let sat = sat_of(&k.extract::<String>()?)?;
        p.pos.insert(sat, rows::<3>(&bytes_of(&v)?).map_err(bad)?);
    }
    for (k, v) in src.get_item(3)?.cast::<PyDict>()?.iter() {
        let sat = sat_of(&k.extract::<String>()?)?;
        p.clk_sp3.insert(sat, f64s(&bytes_of(&v)?).map_err(bad)?);
    }
    for (k, v) in src.get_item(4)?.cast::<PyDict>()?.iter() {
        let sat = sat_of(&k.extract::<String>()?)?;
        let tv = v.cast::<PyTuple>()?;
        let ct = f64s(&bytes_of(&tv.get_item(0)?)?).map_err(bad)?;
        let cv = f64s(&bytes_of(&tv.get_item(1)?)?).map_err(bad)?;
        p.clk.insert(sat, (ct, cv));
    }
    Ok(OrbitSource::Precise(p))
}

fn opt_path(o: Option<Bound<'_, PyAny>>) -> PyResult<Option<PathBuf>> {
    match o {
        Some(v) if !v.is_none() => Ok(Some(PathBuf::from(v.extract::<String>()?))),
        _ => Ok(None),
    }
}

// ---- RINEX reads and raw decoding in Rust -------------------------------------------------------------------

/// One station's observations read by Rust ([`read_obs`]), kept in Rust: `static_solve` takes the handle
/// directly (`{"handle": station, "name", "antenna", ...}`), so the arrays never become Python objects. Python
/// sees the header facts and the epochs (`t_bytes`) a caller needs around the solve.
#[pyclass(frozen, name = "Station", module = "geomaticist_gnss_py")]
struct RustStation {
    st: Mutex<Arc<StationObs>>,
    seconds: f64,
}

impl RustStation {
    fn with<T>(&self, f: impl FnOnce(&StationObs) -> T) -> PyResult<T> {
        let g = self.st.lock().map_err(|e| bad(e.to_string()))?;
        Ok(f(&g))
    }
}

#[pymethods]
impl RustStation {
    #[getter]
    fn name(&self) -> PyResult<String> {
        self.with(|s| s.name.clone())
    }

    #[getter]
    fn path(&self) -> PyResult<String> {
        self.with(|s| s.path.display().to_string())
    }

    #[getter]
    fn antenna(&self) -> PyResult<String> {
        self.with(|s| s.antenna.clone())
    }

    #[getter]
    fn receiver(&self) -> PyResult<String> {
        self.with(|s| s.receiver.clone())
    }

    #[getter]
    fn approx(&self) -> PyResult<Vec<f64>> {
        self.with(|s| s.approx.to_vec())
    }

    #[getter]
    fn delta_hen(&self) -> PyResult<Vec<f64>> {
        self.with(|s| s.delta_hen.to_vec())
    }

    /// Seconds the read took.
    #[getter]
    fn seconds(&self) -> f64 {
        self.seconds
    }

    /// Satellites with data.
    #[getter]
    fn n_sats(&self) -> PyResult<usize> {
        self.with(|s| s.sats.len())
    }

    /// The epochs (GPS seconds) as little-endian f64 bytes.
    fn t_bytes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let b = self.with(|s| {
            s.t.iter()
                .flat_map(|x| x.to_le_bytes())
                .collect::<Vec<u8>>()
        })?;
        Ok(PyBytes::new(py, &b))
    }

    /// Epoch flag 1 per epoch, one byte each (0 / 1).
    fn power_fail_bytes<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let b = self.with(|s| {
            s.power_fail
                .iter()
                .map(|&x| u8::from(x))
                .collect::<Vec<u8>>()
        })?;
        Ok(PyBytes::new(py, &b))
    }
}

/// Navigation data read by Rust ([`read_nav`]), kept in Rust: `static_solve` takes it as the source
/// `("broadcast_handle", nav)`.
#[pyclass(frozen, name = "Nav", module = "geomaticist_gnss_py")]
struct RustNav {
    nav: Arc<NavData>,
    seconds: f64,
}

#[pymethods]
impl RustNav {
    /// GLONASS slot -> frequency channel.
    fn glo_fcn(&self) -> BTreeMap<String, i32> {
        self.nav
            .glo_fcn()
            .iter()
            .map(|(s, v)| (s.to_string(), *v))
            .collect()
    }

    /// (file name, error) of files that could not be read.
    #[getter]
    fn bad(&self) -> Vec<(String, String)> {
        self.nav.bad.clone()
    }

    /// Records left out -> how many.
    #[getter]
    fn skipped(&self) -> BTreeMap<String, usize> {
        self.nav.skipped.clone()
    }

    /// Satellites with Keplerian / GLONASS records.
    #[getter]
    fn n_sats(&self) -> usize {
        self.nav.kep.len() + self.nav.glo.len()
    }

    /// Seconds the read took.
    #[getter]
    fn seconds(&self) -> f64 {
        self.seconds
    }
}

/// Precise orbits and clocks read by Rust ([`read_sp3`]), kept in Rust: `static_solve` and `ppk_solve` take
/// them as the source `("precise_handle", precise)`.
#[pyclass(frozen, name = "Precise", module = "geomaticist_gnss_py")]
struct RustPrecise {
    p: Arc<Precise>,
}

#[pymethods]
impl RustPrecise {
    /// The SP3 epochs (GPS seconds): first, last, count.
    #[getter]
    fn span(&self) -> (f64, f64, usize) {
        let t = &self.p.t;
        (
            t.first().copied().unwrap_or(f64::NAN),
            t.last().copied().unwrap_or(f64::NAN),
            t.len(),
        )
    }

    /// The satellites with positions, e.g. `["C01", "E01", "G01", ...]`.
    #[getter]
    fn sats(&self) -> Vec<String> {
        self.p.pos.keys().map(ToString::to_string).collect()
    }

    /// Interpolated positions of `sat` at GPS seconds `t` -> `[(x, y, z, clock_s, ok), ...]`: ECEF metres,
    /// the clock in seconds, `ok` false outside the orbit's span.
    #[allow(clippy::needless_pass_by_value)]
    fn position(&self, sat: &str, t: Vec<f64>) -> PyResult<Vec<SatPos>> {
        let (p, c, ok) = self.p.position(sat_of(sat)?, &t);
        Ok(p.iter()
            .zip(c)
            .zip(ok)
            .map(|((x, c), ok)| (x[0], x[1], x[2], c, ok))
            .collect())
    }
}

/// One interpolated SP3 position: x, y, z (ECEF metres), clock (seconds), inside the orbit's span.
type SatPos = (f64, f64, f64, f64, bool);

/// Reads SP3 orbit files and, optionally, RINEX clock files without the GIL -> a `Precise` handle (pass it
/// as `("precise_handle", precise)` in the orbit sources). Raises `RuntimeError` when unreadable.
#[pyfunction]
#[pyo3(signature = (sp3_paths, clk_paths=Vec::new()))]
#[allow(clippy::needless_pass_by_value)]
fn read_sp3(
    py: Python<'_>,
    sp3_paths: Vec<String>,
    clk_paths: Vec<String>,
) -> PyResult<RustPrecise> {
    let p = py
        .detach(|| Precise::open(&sp3_paths, &clk_paths))
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok(RustPrecise { p: Arc::new(p) })
}

/// The summary of a RINEX observation file (plain, gzip or Hatanaka) that `gnss rinex info --json` prints ->
/// JSON: the header, first / last epoch, interval, epochs and completeness, gaps, satellites per system,
/// receiver-flagged slips, event records and NGS OPUS checks.
#[pyfunction]
fn rinex_info(py: Python<'_>, path: &str) -> PyResult<String> {
    let v = py
        .detach(|| {
            geomaticist_gnss_rinex::reader::RinexObs::open(Path::new(path))
                .map(|obs| geomaticist_gnss_rinex::stats::summarize(&obs, None, None, None, None))
        })
        .map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    serde_json::to_string(&v).map_err(bad)
}

/// ECEF metres -> `(latitude, longitude, height)`: degrees and metres on WGS 84.
#[pyfunction]
fn ecef_to_llh(xyz: &Bound<'_, PyAny>) -> PyResult<(f64, f64, f64)> {
    let (lat, lon, h) = geomaticist_gnss_engine::models::ecef_to_llh(xyz_of(xyz)?);
    Ok((lat.to_degrees(), lon.to_degrees(), h))
}

/// Latitude, longitude (degrees) and height (metres) on WGS 84 -> ECEF `(x, y, z)` metres.
#[pyfunction]
fn llh_to_ecef(lat: f64, lon: f64, h: f64) -> (f64, f64, f64) {
    let x = geomaticist_gnss_engine::models::llh_to_ecef(lat.to_radians(), lon.to_radians(), h);
    (x[0], x[1], x[2])
}

/// Reads a RINEX 2/3/4 observation file (plain, gzip or Hatanaka) between GPS seconds `t0` and `t1` without the
/// GIL -> a `Station` handle. Raises `RuntimeError` when unreadable.
/// `kinds` selects the observation kinds kept: `"CLS"` (code, phase, signal strength) by default, `"CLSD"`
/// adds Doppler for [`ppk::ppk_solve`].
#[pyfunction]
#[pyo3(signature = (path, name="", t0=None, t1=None, systems="GRECJ", kinds="CLS"))]
fn read_obs(
    py: Python<'_>,
    path: &str,
    name: &str,
    t0: Option<f64>,
    t1: Option<f64>,
    systems: &str,
    kinds: &str,
) -> PyResult<RustStation> {
    let (r, seconds) = py.detach(|| {
        let t = Instant::now();
        // kinds "CLS" is the static engine's reader; kinematic processing asks for "CLSD" (with Doppler)
        let r = geomaticist_gnss_rinex::obs::read_obs_kinds(
            Path::new(path),
            name,
            t0,
            t1,
            systems,
            kinds,
        );
        (r, t.elapsed().as_secs_f64())
    });
    let st = r.map_err(|e| PyRuntimeError::new_err(format!("{path}: {e}")))?;
    Ok(RustStation {
        st: Mutex::new(Arc::new(st)),
        seconds,
    })
}

/// Reads navigation files without the GIL (unreadable files in `bad`, missing paths skipped) -> a `Nav`
/// handle.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
fn read_nav(py: Python<'_>, paths: Vec<String>) -> RustNav {
    let (nav, seconds) = py.detach(|| {
        let t = Instant::now();
        let nav = geomaticist_gnss_rinex::nav::read_nav(&paths);
        (nav, t.elapsed().as_secs_f64())
    });
    RustNav {
        nav: Arc::new(nav),
        seconds,
    }
}

/// The raw-log kind of a path for [`decode_raw`]: `"t04"`, `"packed"`, `"mdb"`, `"tps"`, `"greis"`, `"rtcm3"`,
/// `"ubx"`, `"novatel"`, `"sbf"` or `None`.
#[pyfunction]
fn raw_kind(path: &str) -> Option<&'static str> {
    use geomaticist_gnss_decoders::RawKind;
    geomaticist_gnss_decoders::raw_kind(Path::new(path)).map(|k| match k {
        RawKind::T04 => "t04",
        RawKind::Packed => "packed",
        RawKind::Mdb => "mdb",
        RawKind::Tps => "tps",
        RawKind::Greis => "greis",
        RawKind::Rtcm3 => "rtcm3",
        RawKind::Ubx => "ubx",
        RawKind::NovAtel => "novatel",
        RawKind::Sbf => "sbf",
    })
}

/// What [`decode_raw`] returns: obs path, nav paths, warnings, log lines.
type Decoded = (String, Vec<String>, Vec<String>, Vec<String>);

/// Decodes a raw log into `out_dir` without the GIL -> `(obs path, [nav paths], [warnings], [log lines])`.
/// Raises `RuntimeError` when it cannot be decoded.
#[pyfunction]
fn decode_raw(py: Python<'_>, path: &str, out_dir: &str) -> PyResult<Decoded> {
    let (r, lines) = py.detach(|| {
        let mut lines = Vec::new();
        let r = geomaticist_gnss_decoders::convert(Path::new(path), Path::new(out_dir), &mut |m| {
            lines.push(m.to_owned());
        });
        (r, lines)
    });
    let c = r.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok((
        c.obs.display().to_string(),
        c.nav.iter().map(|p| p.display().to_string()).collect(),
        c.warnings,
        lines,
    ))
}

// ---- the Python functions ---------------------------------------------------------------------------------

/// The bindings' version (`geomaticist-gnss-py X.Y.Z`), for logs and reports.
#[pyfunction]
fn version() -> &'static str {
    VERSION
}

/// One static baseline. Arguments: see the module docs. Returns JSON: `{"result": StaticResult::to_json(),
/// "seconds": {"convert", "solve"}, "stages_written": [...]}`. Raises `RuntimeError` for an engine error
/// and `ValueError` for an argument it cannot take.
#[pyfunction]
#[pyo3(signature = (base, rover, sources, base_xyz, options_json, antex=None, fcn=None, rover_xyz0=None, tropo=None, dump_dir=None))]
#[allow(
    clippy::too_many_arguments,
    clippy::needless_pass_by_value,
    clippy::too_many_lines
)]
fn static_solve(
    py: Python<'_>,
    base: &Bound<'_, PyDict>,
    rover: &Bound<'_, PyDict>,
    sources: &Bound<'_, PyList>,
    base_xyz: &Bound<'_, PyAny>,
    options_json: &str,
    antex: Option<String>,
    fcn: Option<BTreeMap<String, i32>>,
    rover_xyz0: Option<Bound<'_, PyAny>>,
    tropo: Option<Bound<'_, PyDict>>,
    dump_dir: Option<String>,
) -> PyResult<String> {
    let t_conv = Instant::now();
    let a = station_arc(base)?;
    let b = station_arc(rover)?;
    let (srcs, externals) = orbit_sources(sources)?;
    let orb = Orbits::new(srcs);
    let base_xyz = xyz_of(base_xyz)?;
    let rover0 = match rover_xyz0 {
        Some(v) if !v.is_none() => Some(xyz_of(&v)?),
        _ => None,
    };
    let opts = options_from_json(options_json).map_err(bad)?;
    let fcn = fcn_of(fcn)?;
    let mut atx = antex_of(antex)?;
    let TropoOwned { gpt3, orog, vmf3 } = tropo_owned(tropo.as_ref())?;
    let convert_s = t_conv.elapsed().as_secs_f64();

    // the solve itself without the GIL: other Python threads (e.g. parallel baselines) run meanwhile;
    // a Python orbit source re-attaches for its calls
    let dump = dump_dir.map(PathBuf::from);
    let out = py.detach(move || {
        let t_solve = Instant::now();
        let tin = TropoInputs {
            gpt3: gpt3.as_deref(),
            orography: orog.as_deref().map(Vec::as_slice),
            vmf3: &vmf3,
        };
        let r = solve(
            &a,
            &b,
            &orb,
            base_xyz,
            &opts,
            atx.as_mut(),
            &fcn,
            rover0,
            &tin,
        );
        let solve_s = t_solve.elapsed().as_secs_f64();
        let ext_err = externals
            .iter()
            .find_map(|e| e.error.lock().ok().and_then(|g| g.clone()));
        (r, solve_s, ext_err, dump)
    });
    let (r, solve_s, ext_err, dump) = out;
    if let Some(e) = ext_err {
        return Err(PyRuntimeError::new_err(format!("orbit source: {e}")));
    }
    let (res, stages) = r.map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mut written = Vec::new();
    if let Some(d) = dump {
        written = write_stages(&d, &stages)
            .map_err(|e| PyRuntimeError::new_err(format!("stage dumps: {e}")))?;
    }
    let v = json!({
        "result": res.to_json(),
        "seconds": {"convert": convert_s, "solve": solve_s, "engine": res.seconds},
        "stages_written": written.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
        "bridge": VERSION,
    });
    Ok(json_text(&v))
}

/// The orbit sources of a solve from the Python list of `(kind, ...)` tuples, and the Python-side sources
/// among them (their errors are collected after the solve).
fn orbit_sources(sources: &Bound<'_, PyList>) -> PyResult<(Vec<OrbitSource>, Vec<Arc<PyOrbit>>)> {
    let mut srcs = Vec::new();
    let mut externals: Vec<Arc<PyOrbit>> = Vec::new();
    for s in sources.iter() {
        let s = s.cast::<PyTuple>()?;
        let kind: String = s.get_item(0)?.extract()?;
        match kind.as_str() {
            "broadcast" => srcs.push(broadcast(s)?),
            "broadcast_handle" => {
                let h = s.get_item(1)?;
                let h = h
                    .cast::<RustNav>()
                    .map_err(|_| bad("broadcast_handle: not a Nav"))?;
                srcs.push(OrbitSource::Broadcast(Broadcast::new(
                    h.get().nav.as_ref().clone(),
                )));
            }
            "precise" => srcs.push(precise(s)?),
            "precise_handle" => {
                let h = s.get_item(1)?;
                let h = h
                    .cast::<RustPrecise>()
                    .map_err(|_| bad("precise_handle: not a Precise"))?;
                srcs.push(OrbitSource::Precise(h.get().p.as_ref().clone()));
            }
            "external" => {
                let obj = s.get_item(1)?;
                let com = obj
                    .getattr("center_of_mass")
                    .and_then(|c| c.extract::<bool>())
                    .unwrap_or(false);
                let e = Arc::new(PyOrbit {
                    obj: obj.unbind(),
                    com,
                    error: Mutex::new(None),
                });
                externals.push(Arc::clone(&e));
                srcs.push(OrbitSource::External(Box::new(SharedOrbit(e))));
            }
            k => return Err(bad(format!("orbit source kind {k:?}"))),
        }
    }
    Ok((srcs, externals))
}

/// GLONASS frequency channels from `{"R07": 5, ...}`.
fn fcn_of(fcn: Option<BTreeMap<String, i32>>) -> PyResult<BTreeMap<Sat, i32>> {
    fcn.unwrap_or_default()
        .iter()
        .map(|(k, v)| Ok((sat_of(k)?, *v)))
        .collect()
}

/// The ANTEX file of a path (parsed once per process), or none.
fn antex_of(antex: Option<String>) -> PyResult<Option<Antex>> {
    Ok(match antex {
        Some(p) if !p.is_empty() => Some(
            cached(&ANTEX, Path::new(&p), |p| {
                Antex::open(p).map_err(|e| e.to_string())
            })
            .map_err(bad)?
            .as_ref()
            .clone(),
        ),
        _ => None,
    })
}

/// The troposphere grids of a solve, owned.
struct TropoOwned {
    gpt3: Option<Arc<Gpt3Grid>>,
    orog: Option<Arc<Vec<f64>>>,
    vmf3: Vec<Vmf3Grid>,
}

/// The troposphere grids named in `{"gpt3": path, "orography": path, "vmf3": [paths]}`.
fn tropo_owned(tropo: Option<&Bound<'_, PyDict>>) -> PyResult<TropoOwned> {
    let (mut gpt3, mut orog, mut vmf3) = (None, None, Vec::new());
    if let Some(tr) = tropo {
        if let Some(p) = opt_path(tr.get_item("gpt3")?)? {
            gpt3 = Some(
                cached(&GPT3, &p, |p| Gpt3Grid::open(p).map_err(|e| e.to_string())).map_err(bad)?,
            );
        }
        if let Some(p) = opt_path(tr.get_item("orography")?)? {
            orog = Some(
                cached(&OROG, &p, |p| load_orography(p).map_err(|e| e.to_string())).map_err(bad)?,
            );
        }
        if let Some(l) = tr.get_item("vmf3")? {
            for p in l.extract::<Vec<String>>()? {
                let p = PathBuf::from(p);
                let tg = vmf3_epoch_of(&p)
                    .ok_or_else(|| bad(format!("{}: not a VMF3 file name", p.display())))?;
                let rows =
                    cached(&VMF3, &p, |p| load_vmf3(p).map_err(|e| e.to_string())).map_err(bad)?;
                vmf3.push(Vmf3Grid {
                    t: tg,
                    rows: rows.as_ref().clone(),
                });
            }
        }
    }
    Ok(TropoOwned { gpt3, orog, vmf3 })
}

/// The Python module `geomaticist_gnss_py`.
#[pymodule]
fn geomaticist_gnss_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(static_solve, m)?)?;
    m.add_function(wrap_pyfunction!(read_obs, m)?)?;
    m.add_function(wrap_pyfunction!(read_nav, m)?)?;
    m.add_function(wrap_pyfunction!(raw_kind, m)?)?;
    m.add_function(wrap_pyfunction!(decode_raw, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::ppk_solve, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::read_mrk, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::rinex_events, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::rinex_occupations, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::check_times, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::insert_events, m)?)?;
    m.add_function(wrap_pyfunction!(ppk::mark_kinematic, m)?)?;
    m.add_function(wrap_pyfunction!(ppp::ppp_solve, m)?)?;
    m.add_function(wrap_pyfunction!(read_sp3, m)?)?;
    m.add_function(wrap_pyfunction!(rinex_info, m)?)?;
    m.add_function(wrap_pyfunction!(ecef_to_llh, m)?)?;
    m.add_function(wrap_pyfunction!(llh_to_ecef, m)?)?;
    m.add_function(wrap_pyfunction!(synth::synth_session, m)?)?;
    m.add_class::<RustStation>()?;
    m.add_class::<RustNav>()?;
    m.add_class::<RustPrecise>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn buffers_and_records() {
        let b: Vec<u8> = [1.5f64, -2.0, 3.25]
            .iter()
            .flat_map(|x| x.to_le_bytes())
            .collect();
        assert_eq!(f64s(&b).unwrap(), vec![1.5, -2.0, 3.25]);
        assert!(f64s(&b[..7]).is_err());
        assert_eq!(rows::<3>(&b).unwrap(), vec![[1.5, -2.0, 3.25]]);
        assert!(rows::<2>(&b).is_err());
        let mut r = [0.0; 26];
        for (i, x) in (0u8..).zip(r.iter_mut()) {
            *x = f64::from(i);
        }
        assert_eq!(kepler(&r).to_array(), r);
        let g = glonass(&[
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
        ]);
        assert_eq!((g.toe, g.pos[2], g.acc[0], g.fcn), (1.0, 6.0, 10.0, 14.0));
    }

    #[test]
    fn keys_options_and_names() {
        let k = signal_key("L", 1, "C").unwrap();
        assert_eq!((k.kind, k.band, k.attr), ('L', 1, Some('C')));
        assert_eq!(signal_key("C", 2, "").unwrap().attr, None);
        assert!(signal_key("LC", 1, "").is_err());
        let o = options_from_json(r#"{"systems": "G", "elmask": 15.0}"#).unwrap();
        assert_eq!((o.systems.as_str(), o.elmask), ("G", 15.0));
        assert!(options_from_json(r#"{"no_such_field": 1}"#).is_err());
        assert!(vmf3_epoch_of(Path::new("VMF3_20240923.H18")).is_some());
        assert!(VERSION.starts_with("geomaticist-gnss-py "));
        assert!(same_bits(&[1.0, 0.0, 2.5], &[1.0, 0.0, 2.5]));
        assert!(!same_bits(&[0.0, 0.0, 0.0], &[-0.0, 0.0, 0.0]));
        assert_eq!(raw_kind("a/12341541.T04"), Some("t04"));
        assert_eq!(raw_kind("x.rnx"), None);
    }
}
