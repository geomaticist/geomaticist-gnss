//! `gnss baseline`: one static baseline through the engine, from RINEX observation and navigation files
//! (and optional SP3 / CLK, ANTEX, troposphere grids) to the vector, its realistic covariance and the fix
//! record, as text or JSON; optionally the seven stage dumps.
//!
//! No downloads (rule E1: every product is a path given here).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Args;
use geomaticist_gnss_core::Sat;
use geomaticist_gnss_core::time::{parse_user_time, ticks_to_seconds};
use geomaticist_gnss_engine::dump::{json_text, write_stages};
use geomaticist_gnss_engine::options::Options;
use geomaticist_gnss_engine::orbits::{Broadcast, OrbitSource, Orbits};
use geomaticist_gnss_engine::tropo::{TropoInputs, Vmf3Grid};
use geomaticist_gnss_products::grids::{Gpt3Grid, load_orography, load_vmf3};
use geomaticist_gnss_products::{Antex, Precise};
use geomaticist_gnss_rinex::nav::read_nav;
use geomaticist_gnss_rinex::obs::read_obs;
use serde_json::{Value, json};

use crate::CliError;

/// `gnss baseline`.
#[derive(Debug, Args)]
pub struct BaselineArgs {
    /// Base station observation file (plain, .gz, Hatanaka).
    pub base: PathBuf,
    /// Rover observation file.
    pub rover: PathBuf,
    /// Base ARP, ECEF metres: `X,Y,Z` (held).
    #[arg(long = "base-xyz", value_name = "X,Y,Z")]
    pub base_xyz: String,
    /// Rover start position, ECEF metres (default: the file's APPROX POSITION, else a code solution).
    #[arg(long = "rover-xyz0", value_name = "X,Y,Z")]
    pub rover_xyz0: Option<String>,
    /// Broadcast navigation files (RINEX 2 / 3 / 4), repeatable; merged.
    #[arg(long, value_name = "FILE")]
    pub nav: Vec<PathBuf>,
    /// SP3 precise orbit files (used before the broadcast ephemerides), repeatable.
    #[arg(long, value_name = "FILE")]
    pub sp3: Vec<PathBuf>,
    /// RINEX clock files for the SP3 orbits, repeatable.
    #[arg(long, value_name = "FILE")]
    pub clk: Vec<PathBuf>,
    /// ANTEX file (igs20.atx).
    #[arg(long, value_name = "FILE")]
    pub antex: Option<PathBuf>,
    /// First epoch (inclusive), e.g. `2024-09-23 18:00:00`.
    #[arg(long)]
    pub start: Option<String>,
    /// Last epoch (inclusive).
    #[arg(long)]
    pub end: Option<String>,
    /// Systems to read (default GRECJ).
    #[arg(long, default_value = "GRECJ")]
    pub systems: String,
    /// Solver options as a JSON object (`Options` field names); missing fields keep their defaults.
    #[arg(long, value_name = "FILE")]
    pub options: Option<PathBuf>,
    /// GPT3 1° grid (`gpt3_1.grd`) for `trop_model` gpt3 / vmf3.
    #[arg(long, value_name = "FILE")]
    pub gpt3: Option<PathBuf>,
    /// VMF3 orography (`orography_ell_1x1`).
    #[arg(long, value_name = "FILE")]
    pub orography: Option<PathBuf>,
    /// VMF3 6-hourly grids (`VMF3_yyyymmdd.Hhh`), repeatable.
    #[arg(long, value_name = "FILE")]
    pub vmf3: Vec<PathBuf>,
    /// Write the seven stage dumps (`NN_stage.npz / .json`) into this folder.
    #[arg(long = "dump-dir", value_name = "DIR")]
    pub dump_dir: Option<PathBuf>,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

pub(crate) fn xyz(s: &str, what: &str) -> Result<[f64; 3], CliError> {
    let v: Vec<f64> = s
        .split(',')
        .map(|x| x.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .map_err(|e| CliError(format!("{what}: {e}")))?;
    if v.len() != 3 {
        return Err(CliError(format!("{what}: three numbers X,Y,Z expected")));
    }
    Ok([v[0], v[1], v[2]])
}

pub(crate) fn vmf3_epoch_of(path: &Path) -> Option<f64> {
    // VMF3_yyyymmdd.Hhh
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

/// Runs the baseline and returns the text or JSON to print.
#[allow(clippy::too_many_lines)]
pub fn baseline(a: &BaselineArgs) -> Result<String, CliError> {
    let t_all = Instant::now();
    let e = |x: &dyn std::fmt::Display| CliError(x.to_string());
    let base_xyz = xyz(&a.base_xyz, "--base-xyz")?;
    let rover0 = a
        .rover_xyz0
        .as_deref()
        .map(|s| xyz(s, "--rover-xyz0"))
        .transpose()?;
    let t = |s: &Option<String>| -> Result<Option<f64>, CliError> {
        s.as_deref()
            .map(|x| parse_user_time(x).map(ticks_to_seconds))
            .transpose()
            .map_err(CliError)
    };
    let (t0, t1) = (t(&a.start)?, t(&a.end)?);
    let mut opts = Options::default();
    if let Some(p) = &a.options {
        let text = std::fs::read_to_string(p).map_err(|x| e(&x))?;
        let over: serde_json::Map<String, Value> =
            serde_json::from_str(&text).map_err(|x| e(&x))?;
        let mut ov = serde_json::to_value(&opts).map_err(|x| e(&x))?;
        for (k, v) in over {
            if ov.get(&k).is_none() {
                return Err(CliError(format!("--options: unknown field {k:?}")));
            }
            ov[k] = v;
        }
        opts = serde_json::from_value(ov).map_err(|x| e(&x))?;
    }
    if a.dump_dir.is_some() {
        opts.debug = true;
    }
    let t_read = Instant::now();
    let stem = |p: &Path| {
        p.file_stem().map_or_default(|s| {
            s.to_string_lossy()
                .chars()
                .take(4)
                .collect::<String>()
                .to_uppercase()
        })
    };
    let base = read_obs(&a.base, &stem(&a.base), t0, t1, &a.systems).map_err(|x| e(&x))?;
    let rover = read_obs(&a.rover, &stem(&a.rover), t0, t1, &a.systems).map_err(|x| e(&x))?;
    if base.t.is_empty() || rover.t.is_empty() {
        return Err(CliError(
            "no observations in the window (check --start / --end)".into(),
        ));
    }
    let nav = read_nav(&a.nav);
    for (name, why) in &nav.bad {
        eprintln!("note: navigation file {name} not readable ({why}): left out");
    }
    let fcn: BTreeMap<Sat, i32> = nav.glo_fcn();
    let mut srcs = Vec::new();
    if !a.sp3.is_empty() {
        srcs.push(OrbitSource::Precise(
            Precise::open(&a.sp3, &a.clk).map_err(|x| e(&x))?,
        ));
    }
    srcs.push(OrbitSource::Broadcast(Broadcast::new(nav)));
    let orb = Orbits::new(srcs);
    let mut atx = a
        .antex
        .as_deref()
        .map(Antex::open)
        .transpose()
        .map_err(|x| e(&x))?;
    let gpt3 = a
        .gpt3
        .as_deref()
        .map(Gpt3Grid::open)
        .transpose()
        .map_err(|x| e(&x))?;
    let orog = a
        .orography
        .as_deref()
        .map(load_orography)
        .transpose()
        .map_err(|x| e(&x))?;
    let mut vmf3 = Vec::new();
    for p in &a.vmf3 {
        let Some(tg) = vmf3_epoch_of(p) else {
            return Err(CliError(format!(
                "{}: not a VMF3_yyyymmdd.Hhh file name",
                p.display()
            )));
        };
        vmf3.push(Vmf3Grid {
            t: tg,
            rows: load_vmf3(p).map_err(|x| e(&x))?,
        });
    }
    let tropo = TropoInputs {
        gpt3: gpt3.as_ref(),
        orography: orog.as_deref(),
        vmf3: &vmf3,
    };
    let read_s = t_read.elapsed().as_secs_f64();
    let t_solve = Instant::now();
    let (res, stages) = geomaticist_gnss_engine::solve(
        &base,
        &rover,
        &orb,
        base_xyz,
        &opts,
        atx.as_mut(),
        &fcn,
        rover0,
        &tropo,
    )
    .map_err(|x| e(&x))?;
    let solve_s = t_solve.elapsed().as_secs_f64();
    let mut written = Vec::new();
    if let Some(d) = &a.dump_dir {
        written = write_stages(d, &stages).map_err(|x| e(&x))?;
    }
    let (enu, cov_enu) = res.enu();
    if a.json {
        let v = json!({
            "result": res.to_json(),
            "seconds": {"read": read_s, "solve": solve_s, "total": t_all.elapsed().as_secs_f64()},
            "stages_written": written.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
        });
        return Ok(json_text(&v).trim_end().to_owned());
    }
    let status = if res.failed {
        "Failed"
    } else if res.fixed {
        "Fixed"
    } else {
        "Float"
    };
    let mut s = String::new();
    let _ = writeln!(
        s,
        "{} -> {}: {} ({} epochs, {} DD, {} ambiguities, σ0 {:.2})",
        base.name, rover.name, status, res.n_epochs, res.n_dd, res.n_amb, res.sigma0
    );
    let _ = writeln!(
        s,
        "vector XYZ  {:14.4} {:14.4} {:14.4} m",
        res.vector[0], res.vector[1], res.vector[2]
    );
    let _ = writeln!(
        s,
        "vector ENU  {:14.4} {:14.4} {:14.4} m",
        enu[0], enu[1], enu[2]
    );
    let _ = writeln!(
        s,
        "σ ENU       {:14.4} {:14.4} {:14.4} m (realistic)",
        cov_enu[0][0].max(0.0).sqrt(),
        cov_enu[1][1].max(0.0).sqrt(),
        cov_enu[2][2].max(0.0).sqrt()
    );
    let _ = writeln!(
        s,
        "length {:.3} m, rms phase {:.1} mm, rejected {} phase / {} code",
        geomaticist_gnss_engine::linalg::norm3(res.vector),
        1000.0 * res.rms_phase,
        res.n_rejected,
        res.n_rejected_code
    );
    let fx = &res.fix;
    let _ = writeln!(
        s,
        "fix: {} ({}), ratio {}, wide lanes {}, narrow lanes {}{}",
        fx.get("status").and_then(Value::as_str).unwrap_or(""),
        fx.get("mode").and_then(Value::as_str).unwrap_or(""),
        fx.get("ratio").map_or("-".to_owned(), ToString::to_string),
        fx.get("wl").map_or("-".to_owned(), ToString::to_string),
        fx.get("nl").map_or("-".to_owned(), ToString::to_string),
        fx.get("why")
            .and_then(Value::as_str)
            .filter(|w| !w.is_empty())
            .map_or(String::new(), |w| format!("; {w}")),
    );
    for n in &res.notes {
        let _ = writeln!(s, "note: {n}");
    }
    let _ = writeln!(s, "time: read {read_s:.2} s, solve {solve_s:.2} s");
    for p in &written {
        let _ = writeln!(s, "wrote {}", p.display());
    }
    Ok(s.trim_end().to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_xyz_and_vmf3_names() {
        assert_eq!(xyz("1, 2,3", "x").unwrap(), [1.0, 2.0, 3.0]);
        assert!(xyz("1,2", "x").is_err());
        let t = vmf3_epoch_of(Path::new("C:/grids/VMF3_20240923.H18")).unwrap();
        assert!((t - 1_411_149_600.0).abs() < 1e-9);
        assert!(vmf3_epoch_of(Path::new("gpt3_1.grd")).is_none());
    }
}
