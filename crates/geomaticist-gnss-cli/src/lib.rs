//! `geomaticist-gnss-cli`: the `gnss` command line.
//!
//! `gnss rinex info | trim | convert [--json]` for batch work, tests, CI and users; `gnss baseline`
//! runs the static engine ([`baseline`]), `gnss ppp` the PPP ([`ppp`]). `rinex info` prints
//! a file summary; `trim` and `convert` call the RINEX writer. Takes RINEX only (rule R2: no
//! dependency on `geomaticist-gnss-decoders`); `gnss convert <raw> -o out` runs the decoders' own
//! `gnss-decode` program as a subprocess.

#![allow(
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]

use std::fmt::Write as _;
use std::path::PathBuf;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};

pub mod baseline;
pub mod ppk;
pub mod ppp;
use geomaticist_gnss_core::time::{TICKS, fmt, parse_user_time};
use geomaticist_gnss_rinex::reader::RinexObs;
use geomaticist_gnss_rinex::stats::summarize;
use geomaticist_gnss_rinex::writer::{WriteOptions, WriteResult, write_obs};
use serde_json::{Value, json};

/// Crate name, as Cargo knows it.
pub const CRATE: &str = env!("CARGO_PKG_NAME");

/// The `gnss` command line.
#[derive(Debug, Parser)]
#[command(name = "gnss", version, about = "Geomaticist GNSS tools")]
pub struct Cli {
    /// What to do.
    #[command(subcommand)]
    pub command: Command,
}

/// Top-level commands.
#[derive(Debug, Subcommand)]
#[allow(clippy::large_enum_variant)] // parsed once; the baseline's many paths are the point
pub enum Command {
    /// RINEX observation files: describe, trim, convert.
    #[command(subcommand)]
    Rinex(RinexCommand),
    /// A receiver raw log (.T04, .T01/.T00/.T02/.R00, Leica .m00, Topcon .tps) or a logged RTCM 3 stream (.rtcm3) to RINEX, through the `gnss-decode`
    /// program installed beside `gnss`.
    Convert(ConvertArgs),
    /// One static baseline through the engine (RINEX in, vector and covariance out).
    Baseline(baseline::BaselineArgs),
    /// One station's static PPP (RINEX and precise products in, ARP and covariance out).
    Ppp(ppp::PppArgs),
    /// Kinematic processing: a moving receiver against a static base, with positions at events.
    Ppk(ppk::PpkArgs),
    /// Events: list them, check a file's epoch times, insert event records into a RINEX file.
    #[command(subcommand)]
    Events(ppk::EventsCommand),
}

/// `gnss convert`.
#[derive(Debug, Args)]
pub struct ConvertArgs {
    /// Raw log.
    pub file: PathBuf,
    /// Observation file to write (navigation files go beside it), or an existing folder.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Trimble runpkr00 to use for a .T01 / .T00 / .T02 / .R00 the built-in reader doesn't know (default:
    /// $RUNPKR00, beside gnss, PATH).
    #[arg(long)]
    pub runpkr00: Option<PathBuf>,
    /// RTCM 3 logs: a date near the data, `YYYY-MM-DD` (the stream carries no week number; default: the
    /// file's modification time).
    #[arg(long)]
    pub date: Option<String>,
    /// RTCM 3 logs: a navigation file with the GLONASS frequency channels (for MSM4 / MSM6 streams).
    #[arg(long)]
    pub nav: Vec<PathBuf>,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// `gnss rinex ...`.
#[derive(Debug, Subcommand)]
pub enum RinexCommand {
    /// Describe a file (or a selection): epochs, satellites, gaps, OPUS checks.
    Info(InfoArgs),
    /// Cut a time window (and optionally decimate), keeping the version and every kind.
    Trim(WriteArgs),
    /// Write another version (2.11, 3.02-3.05, 4.00-4.02) with chosen kinds and systems.
    Convert(WriteArgs),
}

/// Selection shared by every subcommand.
#[derive(Debug, Args)]
pub struct Selection {
    /// First epoch (inclusive), e.g. `2024-09-23 18:00:00`.
    #[arg(long)]
    pub start: Option<String>,
    /// Last epoch (inclusive).
    #[arg(long)]
    pub end: Option<String>,
    /// New interval in seconds (a multiple of the data's).
    #[arg(long)]
    pub interval: Option<f64>,
}

/// `gnss rinex info`.
#[derive(Debug, Args)]
pub struct InfoArgs {
    /// Observation file (plain, .gz, Hatanaka).
    pub file: PathBuf,
    /// Selection.
    #[command(flatten)]
    pub sel: Selection,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// `gnss rinex trim` / `convert`.
#[derive(Debug, Args)]
pub struct WriteArgs {
    /// Observation file (plain, .gz, Hatanaka).
    pub file: PathBuf,
    /// Output file.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Selection.
    #[command(flatten)]
    pub sel: Selection,
    /// Target version (convert; default 3.04, trim: the source's).
    #[arg(long = "rinex-version")]
    pub rinex_version: Option<String>,
    /// Observation kinds: C code, L phase, D Doppler, S SNR (convert default CL, trim CLDS).
    #[arg(long)]
    pub kinds: Option<String>,
    /// Systems to keep (default all: GRECJIS).
    #[arg(long)]
    pub systems: Option<String>,
    /// Line ending: lf or crlf (default: the source's).
    #[arg(long)]
    pub eol: Option<String>,
    /// Do not add COMMENT lines describing the conversion.
    #[arg(long)]
    pub no_comments: bool,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// A command's failure, printed as `error: ...`.
#[derive(Debug)]
pub struct CliError(pub String);

impl std::fmt::Display for CliError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CliError {}

fn err(e: impl std::fmt::Display) -> CliError {
    CliError(e.to_string())
}

struct Window {
    start: Option<i64>,
    end: Option<i64>,
    interval: Option<i64>,
}

fn window(sel: &Selection) -> Result<Window, CliError> {
    let t = |s: &Option<String>| s.as_deref().map(parse_user_time).transpose().map_err(err);
    let interval = match sel.interval {
        Some(s) if !(s.is_finite() && s > 0.0) => {
            return Err(CliError(format!("interval {s} s must be positive")));
        }
        Some(s) => Some(geomaticist_gnss_core::text::round_half_even(s * TICKS as f64) as i64),
        None => None,
    };
    Ok(Window {
        start: t(&sel.start)?,
        end: t(&sel.end)?,
        interval,
    })
}

fn open(file: &std::path::Path) -> Result<(RinexObs, f64), CliError> {
    let t = Instant::now();
    let obs = RinexObs::open(file).map_err(err)?;
    Ok((obs, t.elapsed().as_secs_f64()))
}

fn result_json(r: &WriteResult, read_s: f64, write_s: f64) -> Value {
    let counts: serde_json::Map<String, Value> = r
        .counts
        .iter()
        .map(|(s, v)| (s.to_string(), json!(v)))
        .collect();
    json!({
        "path": r.path.display().to_string(),
        "version": r.plan.target,
        "systems": r.plan.systems,
        "epochs": r.epochs,
        "first": r.first.map(|t| fmt(t, 7)),
        "last": r.last.map(|t| fmt(t, 7)),
        "first_ticks": r.first,
        "last_ticks": r.last,
        "satellites": r.satellites,
        "counts": counts,
        "notes": r.notes,
        "read_s": read_s,
        "write_s": write_s,
    })
}

fn write(a: &WriteArgs, trim: bool) -> Result<String, CliError> {
    let w = window(&a.sel)?;
    let (obs, read_s) = open(&a.file)?;
    let version = a.rinex_version.clone().unwrap_or_else(|| {
        if trim {
            obs.header.version_str.clone()
        } else {
            "3.04".into()
        }
    });
    let eol = match a.eol.as_deref() {
        None => None,
        Some("lf") => Some("\n".to_owned()),
        Some("crlf") => Some("\r\n".to_owned()),
        Some(x) => return Err(CliError(format!("--eol must be lf or crlf, not {x}"))),
    };
    let opts = WriteOptions {
        version,
        start: w.start,
        end: w.end,
        interval: w.interval,
        kinds: a
            .kinds
            .clone()
            .unwrap_or_else(|| if trim { "CLDS" } else { "CL" }.into()),
        systems: a
            .systems
            .clone()
            .unwrap_or_else(|| geomaticist_gnss_core::system::SYSTEM_ORDER.into()),
        eol,
        add_comments: !a.no_comments,
        ..WriteOptions::default()
    };
    let t = Instant::now();
    let r = write_obs(&obs, &a.output, &opts, None).map_err(err)?;
    let write_s = t.elapsed().as_secs_f64();
    if a.json {
        return serde_json::to_string_pretty(&result_json(&r, read_s, write_s)).map_err(err);
    }
    let mut s = format!(
        "wrote {} (RINEX {}, {}): {} epochs, {} satellites, {} - {}",
        r.path.display(),
        r.plan.target,
        r.plan.systems,
        r.epochs,
        r.satellites,
        r.first.map_or_else(String::new, |t| fmt(t, 1)),
        r.last.map_or_else(String::new, |t| fmt(t, 1)),
    );
    for n in &r.notes {
        s.push_str("\n  note: ");
        s.push_str(n);
    }
    Ok(s)
}

fn info(a: &InfoArgs) -> Result<String, CliError> {
    let w = window(&a.sel)?;
    let (obs, read_s) = open(&a.file)?;
    let mut v = summarize(&obs, w.start, w.end, w.interval, None);
    v["read_s"] = json!(read_s);
    if a.json {
        return serde_json::to_string_pretty(&v).map_err(err);
    }
    let g = |k: &str| match &v[k] {
        Value::String(s) => s.clone(),
        Value::Null => "-".into(),
        x => x.to_string(),
    };
    let mut s = format!(
        "{}\n  RINEX {}  {} - {}  ({}), interval {} s\n  {} epochs ({}% complete), {} satellites, \
         {} cycle slips, {} events",
        g("file"),
        g("rinex_version"),
        g("first_epoch"),
        g("last_epoch"),
        g("duration"),
        g("interval_actual_s"),
        g("epochs"),
        g("completeness_pct"),
        g("satellites"),
        g("cycle_slips"),
        g("events"),
    );
    if let Some(sys) = v["systems"].as_object() {
        // RINEX order (the JSON object's keys are sorted)
        let order = geomaticist_gnss_core::system::SYSTEM_ORDER;
        let mut keys: Vec<&String> = sys.keys().collect();
        keys.sort_by_key(|k| order.find(k.as_str()).unwrap_or(99));
        for k in keys {
            let d = &sys[k];
            let _ = write!(
                s,
                "\n  {k} {:<8} {:>3} sats  {:>7} sat-epochs",
                d["name"].as_str().unwrap_or(""),
                d["count"].as_u64().unwrap_or(0),
                d["epochs"].as_u64().unwrap_or(0)
            );
        }
    }
    if let Some(checks) = v["opus"].as_array() {
        for c in checks {
            let _ = write!(
                s,
                "\n  OPUS {:<4} {}",
                c["level"].as_str().unwrap_or(""),
                c["text"].as_str().unwrap_or("")
            );
        }
    }
    Ok(s)
}

/// The `gnss-decode` program: beside this executable, else on `PATH`.
fn decode_exe() -> Option<PathBuf> {
    let name = if cfg!(windows) {
        "gnss-decode.exe"
    } else {
        "gnss-decode"
    };
    let beside = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|d| d.join(name)))
        .filter(|p| p.is_file());
    beside.or_else(|| {
        std::env::split_paths(&std::env::var_os("PATH")?)
            .map(|d| d.join(name))
            .find(|p| p.is_file())
    })
}

fn convert_raw(a: &ConvertArgs) -> Result<String, CliError> {
    let exe = decode_exe().ok_or_else(|| {
        CliError(
            "gnss-decode (the raw-log decoders) isn't installed beside gnss or on PATH; for RINEX \
             input use `gnss rinex convert`"
                .into(),
        )
    })?;
    let mut cmd = std::process::Command::new(exe);
    cmd.arg(&a.file).arg("-o").arg(&a.output);
    if let Some(r) = &a.runpkr00 {
        cmd.arg("--runpkr00").arg(r);
    }
    if let Some(d) = &a.date {
        cmd.arg("--date").arg(d);
    }
    for n in &a.nav {
        cmd.arg("--nav").arg(n);
    }
    if a.json {
        cmd.arg("--json");
    }
    let out = cmd
        .output()
        .map_err(|e| CliError(format!("gnss-decode didn't start ({e})")))?;
    let text = String::from_utf8_lossy(&out.stdout).trim_end().to_owned();
    if out.status.success() {
        Ok(text)
    } else {
        let e = String::from_utf8_lossy(&out.stderr);
        Err(CliError(e.trim().trim_start_matches("error: ").to_owned()))
    }
}

/// Runs a parsed command line and returns what to print.
pub fn run(cli: &Cli) -> Result<String, CliError> {
    match &cli.command {
        Command::Rinex(RinexCommand::Info(a)) => info(a),
        Command::Rinex(RinexCommand::Trim(a)) => write(a, true),
        Command::Rinex(RinexCommand::Convert(a)) => write(a, false),
        Command::Convert(a) => convert_raw(a),
        Command::Baseline(a) => baseline::baseline(a),
        Command::Ppp(a) => ppp::ppp(a),
        Command::Ppk(a) => ppk::ppk(a),
        Command::Events(c) => ppk::events(c),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn parses_command_lines() {
        Cli::command().debug_assert();
        let c = Cli::try_parse_from([
            "gnss",
            "rinex",
            "convert",
            "a.rnx",
            "-o",
            "b.rnx",
            "--rinex-version",
            "2.11",
            "--interval",
            "30",
            "--json",
        ])
        .unwrap();
        let Command::Rinex(RinexCommand::Convert(a)) = c.command else {
            panic!("convert expected");
        };
        assert_eq!(a.rinex_version.as_deref(), Some("2.11"));
        assert!(a.json);
        assert_eq!(window(&a.sel).unwrap().interval, Some(300_000_000));
        assert!(Cli::try_parse_from(["gnss", "rinex", "info"]).is_err());
    }

    #[test]
    fn round_trip_on_a_small_file() {
        let dir = std::env::temp_dir().join(format!("gnss-cli-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("t.24o");
        std::fs::write(
            &src,
            "     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
     2    L1    C1                                          # / TYPES OF OBSERV
                                                            END OF HEADER
 24  9 23  0  0  0.0000000  0  1G01
 100000000.123   20000000.000
 24  9 23  0  0 30.0000000  0  1G01
 100000001.123   20000001.000
",
        )
        .unwrap();
        let out = dir.join("t.rnx");
        let c = Cli::try_parse_from([
            "gnss",
            "rinex",
            "convert",
            src.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--json",
        ])
        .unwrap();
        let v: Value = serde_json::from_str(&run(&c).unwrap()).unwrap();
        assert_eq!(v["epochs"], 2);
        assert_eq!(v["version"], "3.04");
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("G    2 L1C C1C"));
        let i = Cli::try_parse_from(["gnss", "rinex", "info", src.to_str().unwrap(), "--json"])
            .unwrap();
        let v: Value = serde_json::from_str(&run(&i).unwrap()).unwrap();
        assert_eq!(v["epochs"], 2);
        assert_eq!(v["interval_actual_s"], 30.0);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
