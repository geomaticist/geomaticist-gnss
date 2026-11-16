//! `gnss ppp`: one station's static PPP from its RINEX observation file and the analysis centre's
//! products (SP3 / CLK, Bias-SINEX OSB, ERP, ANTEX; a broadcast navigation file for the GLONASS
//! frequency numbers; an optional BLQ for the ocean loading), to the ARP, its realistic covariance and the
//! AR record, as text or JSON; optionally the seven stage dumps.
//!
//! No downloads (rule E1: every product is a path given here).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use geomaticist_gnss_core::Sat;
use geomaticist_gnss_core::time::{parse_user_time, ticks_to_seconds};
use geomaticist_gnss_engine::dump::json_text;
use geomaticist_gnss_engine::models::ecef_to_llh;
use geomaticist_gnss_engine::orbits::{OrbitSource, Orbits};
use geomaticist_gnss_engine::tropo::{TropoInputs, Vmf3Grid};
use geomaticist_gnss_ppp::dump::write_stages;
use geomaticist_gnss_ppp::{PppOptions, realistic_cov};
use geomaticist_gnss_products::grids::{Gpt3Grid, load_orography, load_vmf3};
use geomaticist_gnss_products::{Antex, Precise, read_blq, read_erp, read_osb};
use geomaticist_gnss_rinex::nav::read_nav;
use geomaticist_gnss_rinex::obs::read_obs;
use serde_json::{Value, json};

use crate::CliError;
use crate::baseline::{vmf3_epoch_of, xyz};

/// `gnss ppp`.
#[derive(Debug, Args)]
#[allow(clippy::struct_excessive_bools)] // command-line switches
pub struct PppArgs {
    /// Observation file (plain, .gz, Hatanaka).
    pub obs: PathBuf,
    /// SP3 precise orbit files (one per day the session spans), repeatable.
    #[arg(long, value_name = "FILE", required = true)]
    pub sp3: Vec<PathBuf>,
    /// RINEX clock files for the SP3 orbits, repeatable.
    #[arg(long, value_name = "FILE")]
    pub clk: Vec<PathBuf>,
    /// ANTEX file (igs20.atx): satellite and receiver antennas.
    #[arg(long, value_name = "FILE")]
    pub antex: Option<PathBuf>,
    /// Bias-SINEX OSB files (code biases; phase biases for PPP-AR), repeatable.
    #[arg(long, value_name = "FILE")]
    pub osb: Vec<PathBuf>,
    /// IGS ERP files (pole tide), repeatable.
    #[arg(long, value_name = "FILE")]
    pub erp: Vec<PathBuf>,
    /// Broadcast navigation files, for the GLONASS frequency numbers only, repeatable.
    #[arg(long, value_name = "FILE")]
    pub nav: Vec<PathBuf>,
    /// BLQ ocean loading coefficients (Onsala format); applied when the options set `otl`.
    #[arg(long, value_name = "FILE")]
    pub blq: Option<PathBuf>,
    /// Station name in the BLQ file (default: the station name).
    #[arg(long = "blq-station", value_name = "NAME")]
    pub blq_station: Option<String>,
    /// Station name (default: the file name's first four characters, upper case).
    #[arg(long, value_name = "NAME")]
    pub name: Option<String>,
    /// Start position, ECEF metres `X,Y,Z` (default: a code solution).
    #[arg(long, value_name = "X,Y,Z")]
    pub xyz0: Option<String>,
    /// First epoch (inclusive), e.g. `2024-09-23 16:00:00`.
    #[arg(long)]
    pub start: Option<String>,
    /// Last epoch (inclusive).
    #[arg(long)]
    pub end: Option<String>,
    /// Systems to use (default GRE).
    #[arg(long, default_value = "GRE")]
    pub systems: String,
    /// PPP options as a JSON object (`PppOptions` field names); missing fields keep their defaults.
    #[arg(long, value_name = "FILE")]
    pub options: Option<PathBuf>,
    /// Resolve the ambiguities (PPP-AR; needs OSB phase biases). Same as `"ar": true` in --options.
    #[arg(long)]
    pub ar: bool,
    /// Apply the ocean loading of --blq. Same as `"otl": true` in --options.
    #[arg(long)]
    pub otl: bool,
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
    /// Kinematic: one position per epoch for a moving receiver (forward, backward and merged), instead of
    /// one for the session. Events and the `.pos` outputs work as for `gnss ppk`.
    #[arg(long)]
    pub kinematic: bool,
    /// Kinematic: DJI `Timestamp.MRK` file(s) with the events, repeatable.
    #[arg(long, value_name = "FILE")]
    pub mrk: Vec<PathBuf>,
    /// Kinematic: take the events from the event records (epoch flag 5) of the observation file.
    #[arg(long = "rinex-events")]
    pub rinex_events: bool,
    /// Kinematic: a text file of events (`WEEK SECONDS [LABEL]` or `GPS_SECONDS [LABEL]` per line).
    #[arg(long = "events", value_name = "FILE")]
    pub events: Option<PathBuf>,
    /// Kinematic: event interpolation, `linear2`, `cubic4` (default) or `hermite2`.
    #[arg(long, value_name = "METHOD")]
    pub interp: Option<String>,
    /// Kinematic: missing epochs next to an event that may be skipped (default 0).
    #[arg(long = "max-missing", value_name = "N")]
    pub max_missing: Option<usize>,
    /// Kinematic: `forward`, `backward` or `combined` (default).
    #[arg(long)]
    pub direction: Option<String>,
    /// Kinematic: options of the kinematic filter as a JSON object (`KinematicOptions` field names
    /// except `ppp`, which --options sets).
    #[arg(long = "kinematic-options", value_name = "FILE")]
    pub kinematic_options: Option<PathBuf>,
    /// Kinematic: write the track as a `.pos` file (events beside it as `NAME_events.pos` / `.csv`).
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Kinematic: also write the forward and backward tracks.
    #[arg(long = "keep-passes")]
    pub keep_passes: bool,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// Runs the PPP and returns the text or JSON to print.
#[allow(clippy::too_many_lines)]
pub fn ppp(a: &PppArgs) -> Result<String, CliError> {
    let t_all = Instant::now();
    let e = |x: &dyn std::fmt::Display| CliError(x.to_string());
    let x0 = a.xyz0.as_deref().map(|s| xyz(s, "--xyz0")).transpose()?;
    let t = |s: &Option<String>| -> Result<Option<f64>, CliError> {
        s.as_deref()
            .map(|x| parse_user_time(x).map(ticks_to_seconds))
            .transpose()
            .map_err(CliError)
    };
    let (t0, t1) = (t(&a.start)?, t(&a.end)?);
    let mut opts = PppOptions {
        systems: a.systems.clone(),
        ..PppOptions::default()
    };
    if let Some(p) = &a.options {
        let text = std::fs::read_to_string(p).map_err(|x| e(&x))?;
        let over: serde_json::Map<String, Value> =
            serde_json::from_str(&text).map_err(|x| e(&x))?;
        let mut ov = serde_json::to_value(&opts).map_err(|x| e(&x))?;
        for (k, v) in over {
            if k == "blq" {
                // the coefficients are handed in with --blq
                if !v.is_null() {
                    return Err(CliError(
                        "--options: hand the BLQ coefficients in with --blq (and set otl)".into(),
                    ));
                }
                continue;
            }
            if ov.get(&k).is_none() {
                return Err(CliError(format!("--options: unknown field {k:?}")));
            }
            ov[k] = v;
        }
        opts = serde_json::from_value(ov).map_err(|x| e(&x))?;
    }
    if a.ar {
        opts.ar = true;
    }
    if a.otl {
        opts.otl = true;
    }
    if a.dump_dir.is_some() {
        opts.debug = true;
    }
    let t_read = Instant::now();
    let name = a.name.clone().unwrap_or_else(|| {
        a.obs
            .file_stem()
            .map(|s| {
                s.to_string_lossy()
                    .chars()
                    .take(4)
                    .collect::<String>()
                    .to_uppercase()
            })
            .unwrap_or_default()
    });
    let st = read_obs(&a.obs, &name, t0, t1, &opts.systems).map_err(|x| e(&x))?;
    if st.t.is_empty() {
        return Err(CliError(
            "no observations in the window (check --start / --end)".into(),
        ));
    }
    let nav = read_nav(&a.nav);
    for (n, why) in &nav.bad {
        eprintln!("note: navigation file {n} not readable ({why}): left out");
    }
    let fcn: BTreeMap<Sat, i32> = nav.glo_fcn();
    let precise = Precise::open(&a.sp3, &a.clk).map_err(|x| e(&x))?;
    let orb = Orbits::new(vec![OrbitSource::Precise(precise)]);
    let mut atx = a
        .antex
        .as_deref()
        .map(Antex::open)
        .transpose()
        .map_err(|x| e(&x))?;
    let osb = read_osb(&a.osb).map_err(|x| e(&x))?;
    let erp = if a.erp.is_empty() {
        None
    } else {
        Some(read_erp(&a.erp).map_err(|x| e(&x))?)
    };
    let blq = match &a.blq {
        None => None,
        Some(p) => {
            let text = std::fs::read(p).map_err(|x| e(&x))?;
            let map = read_blq(&text);
            let key = a
                .blq_station
                .clone()
                .unwrap_or_else(|| name.clone())
                .to_uppercase();
            let Some(c) = map.get(&key) else {
                return Err(CliError(format!(
                    "{}: no station {key:?} in the BLQ file (it has {:?})",
                    p.display(),
                    map.keys().collect::<Vec<_>>()
                )));
            };
            Some(*c)
        }
    };
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
    if a.kinematic {
        return crate::ppk::ppp_kinematic(
            a,
            &st,
            &orb,
            atx.as_mut(),
            &opts,
            &fcn,
            &osb,
            erp.as_ref(),
            blq.as_ref(),
            &tropo,
            x0,
        );
    }
    let t_solve = Instant::now();
    let (res, stages) = geomaticist_gnss_ppp::solve(
        &st,
        &orb,
        atx.as_mut(),
        &opts,
        &fcn,
        &osb,
        erp.as_ref(),
        x0,
        blq.as_ref(),
        &tropo,
    )
    .map_err(|x| e(&x))?;
    let solve_s = t_solve.elapsed().as_secs_f64();
    let mut written = Vec::new();
    if let Some(d) = &a.dump_dir {
        written = write_stages(d, &stages).map_err(|x| e(&x))?;
    }
    let enu = res.enu_cov();
    let real = realistic_cov(&enu, res.hours, res.fixed);
    let (lat, lon, h) = ecef_to_llh(res.xyz);
    if a.json {
        let v = json!({
            "result": res.to_json(),
            "llh_deg": [lat.to_degrees(), lon.to_degrees(), h],
            "cov_realistic_enu": real,
            "seconds": {"read": read_s, "solve": solve_s, "total": t_all.elapsed().as_secs_f64()},
            "stages_written": written.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
        });
        return Ok(json_text(&v).trim_end().to_owned());
    }
    let mut s = String::new();
    let _ = writeln!(
        s,
        "{}: {} ({:.2} h, {} epochs, {} satellites, {} arcs, σ0 {:.2})",
        st.name,
        if res.fixed { "Fixed" } else { "Float" },
        res.hours,
        res.n_epochs,
        res.sats.values().map(Vec::len).sum::<usize>(),
        res.n_amb,
        res.sigma0
    );
    let _ = writeln!(
        s,
        "ARP XYZ     {:14.4} {:14.4} {:14.4} m ({})",
        res.xyz[0],
        res.xyz[1],
        res.xyz[2],
        if res.products.is_empty() {
            "products' frame"
        } else {
            res.products.as_str()
        }
    );
    let _ = writeln!(
        s,
        "lat lon h   {:14.9} {:14.9} {:10.4}",
        lat.to_degrees(),
        lon.to_degrees(),
        h
    );
    let _ = writeln!(
        s,
        "σ ENU       {:14.4} {:14.4} {:14.4} m (realistic; formal {:.4} {:.4} {:.4})",
        real[0][0].max(0.0).sqrt(),
        real[1][1].max(0.0).sqrt(),
        real[2][2].max(0.0).sqrt(),
        enu[0][0].max(0.0).sqrt(),
        enu[1][1].max(0.0).sqrt(),
        enu[2][2].max(0.0).sqrt()
    );
    let _ = writeln!(
        s,
        "rms phase {:.1} mm, code {:.2} m, rejected {}, troposphere {}",
        1000.0 * res.rms_phase,
        res.rms_code,
        res.n_rejected,
        res.troposphere
    );
    if let Some((_, z, sz)) = res.ztd.first() {
        let end = res.ztd.last().map_or((0.0, 0.0, 0.0), |x| *x);
        let _ = writeln!(
            s,
            "ztd {:.4} ± {:.4} m at the start, {:.4} ± {:.4} m at the end",
            z, sz, end.1, end.2
        );
    }
    if opts.ar {
        let _ = writeln!(
            s,
            "ar: ok {}, wide lanes {}, narrow lanes {} of {}, ratio {}{}",
            res.ar.get("ok").map_or("-".to_owned(), ToString::to_string),
            res.ar.get("wl").map_or("-".to_owned(), ToString::to_string),
            res.ar
                .get("n_fixed")
                .map_or("-".to_owned(), ToString::to_string),
            res.ar
                .get("n_cand")
                .map_or("-".to_owned(), ToString::to_string),
            res.ar
                .get("ratio")
                .map_or("-".to_owned(), ToString::to_string),
            res.ar
                .get("why")
                .and_then(Value::as_str)
                .filter(|w| !w.is_empty())
                .map_or(String::new(), |w| format!("; {w}")),
        );
    }
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
    use crate::{Cli, Command};
    use clap::Parser;

    #[test]
    fn parses_a_ppp_command_line() {
        let c = Cli::try_parse_from([
            "gnss", "ppp", "fair.24o", "--sp3", "a.sp3", "--clk", "a.clk", "--osb", "a.bia",
            "--ar", "--otl", "--blq", "x.blq", "--json",
        ])
        .unwrap();
        let Command::Ppp(a) = c.command else {
            panic!("ppp expected");
        };
        assert!(a.ar && a.otl && a.json);
        assert_eq!(a.sp3.len(), 1);
        assert!(Cli::try_parse_from(["gnss", "ppp", "fair.24o"]).is_err());
    }
}
