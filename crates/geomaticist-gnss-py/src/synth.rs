//! A small synthetic data set for the Python examples and their smoke test: simulated base and rover
//! observations written as RINEX 3, the simulated orbits as an SP3 file and, for a flight, a DJI-style
//! `Timestamp.MRK`. Nothing in it comes from a real receiver, so the examples can be shared freely.
//!
//! The simulation is the engine's own ([`geomaticist_gnss_engine::ppk::synth`]: circular orbits, a
//! troposphere, millimetre phase noise, integer ambiguities); the files are what any RINEX reader takes.

use std::fmt::Write as _;
use std::path::Path;

use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::obs::StationObs;
use geomaticist_gnss_core::time::{TICKS, from_ticks, gps_week_tow};
use geomaticist_gnss_engine::ppk::synth::{SynthOrbit, SynthReceiver, enu_point, simulate};
use geomaticist_gnss_rinex::station::{StationHeader, station_to_rinex3};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use serde_json::{Value, json};

use super::bad;

/// The session's first epoch: GPS week 2300, Thursday 00:00 (2024-02-08).
const T0: f64 = 2300.0 * 604_800.0 + 345_600.0;
/// The base mark, ECEF metres (about 46.0 N, 8.0 E, 500 m).
const BASE: [f64; 3] = [4_395_470.63, 617_743.11, 4_565_607.21];

/// Year, month, day, hour, minute, seconds of GPS seconds `t`.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
fn civil(t: f64) -> (i64, i64, i64, i64, i64, f64) {
    let c = from_ticks((t * TICKS as f64).round() as i64);
    let sec = c.sec_ticks as f64 / TICKS as f64;
    (c.year, c.month, c.day, c.hour, c.minute, sec)
}

/// The constellation as SP3 text: positions every 5 minutes from `t_from` to `t_to`, clocks in microseconds.
#[allow(clippy::many_single_char_names)] // y mo d h mi: the SP3 epoch fields
fn sp3(orbit: &SynthOrbit, t_from: f64, t_to: f64) -> String {
    let (y, mo, d, h, mi, s0) = civil(t_from);
    let mut s =
        format!("#cP{y:4} {mo:2} {d:2} {h:2} {mi:2} {s0:11.8}      25 ORBIT IGS20 FIT  SIM\n");
    let mut t = t_from;
    while t <= t_to + 1e-6 {
        let (y, mo, d, h, mi, sec) = civil(t);
        let _ = writeln!(s, "*  {y:4} {mo:2} {d:2} {h:2} {mi:2} {sec:11.8}");
        for sat in &orbit.sats {
            let p = sat.position(t - orbit.t0);
            let clk = sat.clk + sat.drift * (t - orbit.t0);
            let _ = writeln!(
                s,
                "P{}{:14.6}{:14.6}{:14.6}{:14.6}",
                sat.sat,
                p[0] / 1e3,
                p[1] / 1e3,
                p[2] / 1e3,
                clk * 1e6
            );
        }
        t += 300.0;
    }
    s.push_str("EOF\n");
    s
}

/// A RINEX 3 special record (epoch flag 2 or 3) at `t`, with `lines` of header records after it.
#[allow(clippy::many_single_char_names)] // y mo d h mi: the RINEX epoch fields
fn special_record(t: f64, flag: u8, lines: &[String]) -> String {
    let (y, mo, d, h, mi, sec) = civil(t);
    let mut s = format!(
        "> {y:4} {mo:02} {d:02} {h:02} {mi:02}{sec:11.7}  {flag}{:3}\n",
        lines.len()
    );
    for l in lines {
        let _ = writeln!(s, "{l}");
    }
    s
}

/// RINEX 3 text with special records put before the first epoch later than each record's time.
fn with_records(text: &str, epochs: &[f64], records: &[(f64, String)]) -> String {
    let mut out = String::with_capacity(text.len() + 1000);
    let mut k = 0; // epoch index
    let mut r = 0; // record index
    for line in text.split_inclusive('\n') {
        if line.starts_with("> ") {
            let t = epochs.get(k).copied().unwrap_or(f64::INFINITY);
            while r < records.len() && records[r].0 < t {
                out.push_str(&records[r].1);
                r += 1;
            }
            k += 1;
        }
        out.push_str(line);
    }
    for rec in &records[r..] {
        out.push_str(&rec.1);
    }
    out
}

fn header() -> StationHeader {
    StationHeader {
        program: "geomaticist synth".into(),
        comments: vec!["SYNTHETIC DATA: simulated, not from a receiver".into()],
        ..StationHeader::default()
    }
}

fn write(dir: &Path, name: &str, text: &str) -> Result<String, String> {
    let p = dir.join(name);
    std::fs::write(&p, text).map_err(|e| format!("{}: {e}", p.display()))?;
    Ok(p.display().to_string())
}

fn station_file(dir: &Path, name: &str, st: &StationObs) -> Result<String, String> {
    write(dir, name, &station_to_rinex3(st, &header()))
}

fn epochs(seconds: f64, step: f64) -> Vec<f64> {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = (seconds / step).round() as u32;
    (0..=n).map(|k| T0 + step * f64::from(k)).collect()
}

/// Three static stations, an hour at 15 s: `BASE` and the unknowns `PT01`, `PT02`.
fn static_set(dir: &Path, orbit: &SynthOrbit) -> Result<Value, String> {
    let ep = epochs(3600.0, 15.0);
    let pt01 = enu_point(BASE, [1234.567, -567.891, 12.345]);
    let pt02 = enu_point(BASE, [-845.210, 1502.334, -7.812]);
    let mut files = serde_json::Map::new();
    let mut truth = serde_json::Map::new();
    for (name, xyz, seed) in [("BASE", BASE, 11), ("PT01", pt01, 12), ("PT02", pt02, 13)] {
        let st = simulate(orbit, &SynthReceiver::new(name, seed), &|_| xyz, &ep);
        files.insert(
            name.into(),
            json!(station_file(dir, &format!("{name}.obs"), &st)?),
        );
        truth.insert(name.into(), json!(xyz));
    }
    Ok(json!({
        "stations": files,
        "sp3": write(dir, "orbit.sp3", &sp3(orbit, T0 - 3600.0, T0 + 7200.0))?,
        "base_xyz": BASE,
        "truth": truth,
    }))
}

/// A two-minute flight at 5 Hz with a photo every 10 s, logged in a `Timestamp.MRK`.
fn flight_set(dir: &Path, orbit: &SynthOrbit) -> Result<Value, String> {
    let path = |t: f64| {
        let s = t - T0;
        enu_point(
            BASE,
            [
                120.0 * (0.05 * s).sin() + 40.0,
                60.0 * (0.1 * s).sin(),
                90.0,
            ],
        )
    };
    let ep = epochs(120.0, 0.2);
    let base = simulate(orbit, &SynthReceiver::new("BASE", 3), &|_| BASE, &ep);
    let rover = simulate(orbit, &SynthReceiver::new("DRON", 4), &path, &ep);
    // the camera sits 45 mm east, 120 mm south and 193 mm below the antenna (the MRK's N / E / V in mm)
    let (week, tow) = gps_week_tow(T0);
    let mut mrk = String::new();
    let mut truth = Vec::new();
    for k in 1..=11 {
        let dt = 10.0 * f64::from(k) - 0.37;
        let t = T0 + dt;
        let _ = writeln!(
            mrk,
            "{k}\t{:.6}\t[{week}]\t  -120,N\t    45,E\t   193,V\t46.0,Lat\t8.0,Lon\t590.0,Ellh\t1.0, 1.0, 2.0\t50,Q",
            tow + dt
        );
        truth.push(json!({"t": t, "antenna_xyz": path(t)}));
    }
    Ok(json!({
        "base": station_file(dir, "base.obs", &base)?,
        "rover": station_file(dir, "rover.obs", &rover)?,
        "sp3": write(dir, "orbit.sp3", &sp3(orbit, T0 - 3600.0, T0 + 3600.0))?,
        "mrk": write(dir, "flight_Timestamp.MRK", &mrk)?,
        "base_xyz": BASE,
        "truth": truth,
    }))
}

/// Stop and go at 1 s: a minute on each of `PT01`, `PT02`, `PT03`, half a minute walking between them; the
/// rover file carries the flag 3 ("new site occupation") and flag 2 ("start moving") records.
fn stop_and_go_set(dir: &Path, orbit: &SynthOrbit) -> Result<Value, String> {
    let stops = [
        ("PT01", 0.0, 60.0, [25.0, 10.0, 1.8]),
        ("PT02", 90.0, 150.0, [55.0, -20.0, 2.3]),
        ("PT03", 180.0, 240.0, [30.0, -45.0, 1.1]),
    ];
    let enu = move |s: f64| -> [f64; 3] {
        for w in stops.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            if s <= a.2 {
                return a.3;
            }
            if s < b.1 {
                let f = (s - a.2) / (b.1 - a.2);
                return std::array::from_fn(|i| a.3[i] + f * (b.3[i] - a.3[i]));
            }
        }
        stops[stops.len() - 1].3
    };
    let path = move |t: f64| enu_point(BASE, enu(t - T0));
    let ep = epochs(240.0, 1.0);
    let base = simulate(orbit, &SynthReceiver::new("BASE", 5), &|_| BASE, &ep);
    let rover = simulate(orbit, &SynthReceiver::new("ROVR", 6), &path, &ep);
    let mut records = Vec::new();
    let mut truth = serde_json::Map::new();
    for (k, (name, t0, t1, e)) in stops.iter().enumerate() {
        let start = if k == 0 { T0 - 0.5 } else { T0 + t0 - 0.5 };
        records.push((
            start,
            special_record(start, 3, &[format!("{name:<60}MARKER NAME")]),
        ));
        if k + 1 < stops.len() {
            records.push((T0 + t1 + 0.5, special_record(T0 + t1 + 0.5, 2, &[])));
        }
        truth.insert((*name).into(), json!(enu_point(BASE, *e)));
    }
    let text = with_records(&station_to_rinex3(&rover, &header()), &rover.t, &records);
    Ok(json!({
        "base": station_file(dir, "base.obs", &base)?,
        "rover": write(dir, "rover_stopgo.obs", &text)?,
        "sp3": write(dir, "orbit.sp3", &sp3(orbit, T0 - 3600.0, T0 + 3600.0))?,
        "base_xyz": BASE,
        "truth": truth,
    }))
}

/// Writes a synthetic data set into `out_dir` (created if missing) -> JSON with the file paths, the base
/// mark's ECEF coordinates and the true positions to compare results with. `kind`:
/// - `"static"`: three static stations (`BASE`, `PT01`, `PT02`), one hour at 15 s, for baselines and loops;
/// - `"flight"`: a base and a drone rover at 5 Hz for two minutes with 11 photos in a `Timestamp.MRK`;
/// - `"stop_and_go"`: a base and a rover that stops on three points (flag 3 / flag 2 records).
///
/// Every set has an `orbit.sp3` with the simulated orbits; pass it to [`super::read_sp3`].
#[pyfunction]
#[pyo3(signature = (out_dir, kind="static"))]
pub fn synth_session(py: Python<'_>, out_dir: &str, kind: &str) -> PyResult<String> {
    let dir = Path::new(out_dir).to_path_buf();
    let kind = kind.to_owned();
    let v = py
        .detach(move || -> Result<Value, String> {
            std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
            let orbit = SynthOrbit::constellation(
                &[GnssSystem::Gps, GnssSystem::Galileo, GnssSystem::BeiDou],
                6,
                4,
                T0,
            );
            match kind.as_str() {
                "static" => static_set(&dir, &orbit),
                "flight" => flight_set(&dir, &orbit),
                "stop_and_go" => stop_and_go_set(&dir, &orbit),
                k => Err(format!(
                    "kind {k:?}: one of \"static\", \"flight\", \"stop_and_go\""
                )),
            }
        })
        .map_err(PyRuntimeError::new_err)?;
    serde_json::to_string(&v).map_err(bad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn special_records_go_before_later_epochs() {
        let text = "HDR\n> a\nx\n> b\ny\n";
        let r = special_record(T0 + 0.5, 2, &[]);
        let out = with_records(text, &[T0, T0 + 1.0], &[(T0 + 0.5, r.clone())]);
        assert_eq!(out, format!("HDR\n> a\nx\n{r}> b\ny\n"));
        assert!(r.starts_with("> 2024 02 08 00 00  0.5000000  2  0"), "{r}");
        let s = sp3(
            &SynthOrbit::constellation(&[GnssSystem::Gps], 1, 1, T0),
            T0,
            T0 + 300.0,
        );
        assert_eq!(s.lines().filter(|l| l.starts_with('*')).count(), 2);
    }
}
