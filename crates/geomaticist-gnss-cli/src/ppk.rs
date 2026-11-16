//! `gnss ppk`: a moving receiver against a static base, epoch by epoch, with positions at events;
//! `gnss events`: list, check and insert events.
//!
//! Inputs are RINEX (decode raw receiver logs with `gnss convert` first). Outputs are a `.pos` file in
//! RTKLIB's layout for the epochs, and for the events the same layout (`*_events.pos`) and a CSV
//! (`*_events.csv`). See `docs/PPK.md`.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::{Args, Subcommand};
use geomaticist_gnss_core::event::{Event, Occupation};
use geomaticist_gnss_core::time::{fmt, parse_user_time, ticks_to_seconds};
use geomaticist_gnss_core::{Sat, StationObs};
use geomaticist_gnss_engine::models::llh_to_ecef;
use geomaticist_gnss_engine::orbits::{Broadcast, OrbitSource, Orbits};
use geomaticist_gnss_engine::ppk::pos::{PosHeader, epochs_pos, events_csv, events_pos, gpst_text};
use geomaticist_gnss_engine::ppk::{
    self, BaseStation, Direction, EventInterp, EventStatus, PpkInput, PpkOptions, PpkResult,
};
use geomaticist_gnss_engine::tropo::{TropoInputs, Vmf3Grid};
use geomaticist_gnss_events::{mrk, rinex as evrinex};
use geomaticist_gnss_products::grids::{Gpt3Grid, load_orography, load_vmf3};
use geomaticist_gnss_products::{Antex, Precise};
use geomaticist_gnss_rinex::nav::read_nav;
use geomaticist_gnss_rinex::obs::read_obs_bytes_kinds;
use geomaticist_gnss_rinex::open::read_plain;
use serde_json::json;

use crate::CliError;
use crate::baseline::{vmf3_epoch_of, xyz};

/// `gnss ppk`.
#[derive(Debug, Args)]
#[allow(clippy::struct_excessive_bools)] // command-line switches
pub struct PpkArgs {
    /// Base station observation file (RINEX; plain, .gz, Hatanaka).
    pub base: PathBuf,
    /// Rover observation file(s). Several files of one receiver session (a file per flight) are joined.
    #[arg(required = true)]
    pub rover: Vec<PathBuf>,
    /// Base mark, ECEF metres: `X,Y,Z`.
    #[arg(long = "base-xyz", value_name = "X,Y,Z", conflicts_with = "base_llh")]
    pub base_xyz: Option<String>,
    /// Base mark as latitude, longitude (degrees) and ellipsoidal height (m): `LAT,LON,H`.
    #[arg(
        long = "base-llh",
        value_name = "LAT,LON,H",
        allow_hyphen_values = true
    )]
    pub base_llh: Option<String>,
    /// Base antenna reference point above the mark: `H,E,N` metres (default: the base file's header).
    #[arg(long = "base-hen", value_name = "H,E,N", allow_hyphen_values = true)]
    pub base_hen: Option<String>,
    /// Broadcast navigation files (RINEX 2 / 3 / 4), repeatable; merged.
    #[arg(long, value_name = "FILE")]
    pub nav: Vec<PathBuf>,
    /// SP3 precise orbit files (used before the broadcast ephemerides), repeatable.
    #[arg(long, value_name = "FILE")]
    pub sp3: Vec<PathBuf>,
    /// RINEX clock files for the SP3 orbits, repeatable.
    #[arg(long, value_name = "FILE")]
    pub clk: Vec<PathBuf>,
    /// ANTEX file (igs20.atx) for the receiver antennas.
    #[arg(long, value_name = "FILE")]
    pub antex: Option<PathBuf>,
    /// DJI `Timestamp.MRK` file(s) with the events (photo times and antenna-to-camera offsets), repeatable.
    #[arg(long, value_name = "FILE")]
    pub mrk: Vec<PathBuf>,
    /// Take the events from the event records (epoch flag 5) of the rover file(s).
    #[arg(long = "rinex-events")]
    pub rinex_events: bool,
    /// A text file of events: one per line, `WEEK SECONDS [LABEL]` or `GPS_SECONDS [LABEL]`.
    #[arg(long = "events", value_name = "FILE")]
    pub events: Option<PathBuf>,
    /// Do not apply the antenna-to-camera offsets of the MRK file(s): positions stay at the antenna.
    #[arg(long = "no-offsets")]
    pub no_offsets: bool,
    /// Average the stop-and-go occupations marked in the rover file(s) (epoch flags 2 / 3).
    #[arg(long)]
    pub occupations: bool,
    /// Event interpolation: `linear2` (2-point), `cubic4` (4-point, the default) or `hermite2` (2-point with
    /// Doppler velocities, for low-rate data).
    #[arg(long, value_name = "METHOD")]
    pub interp: Option<String>,
    /// Missing epochs next to an event that may be skipped (default 0).
    #[arg(long = "max-missing", value_name = "N")]
    pub max_missing: Option<usize>,
    /// `forward`, `backward` or `combined` (default).
    #[arg(long)]
    pub direction: Option<String>,
    /// Elevation mask, degrees.
    #[arg(long)]
    pub elmask: Option<f64>,
    /// Systems to use (default GRECJ).
    #[arg(long)]
    pub systems: Option<String>,
    /// Float solutions only (no integer ambiguity fixing).
    #[arg(long = "no-fix")]
    pub no_fix: bool,
    /// Options as a JSON object (`PpkOptions` field names); missing fields keep their defaults.
    #[arg(long, value_name = "FILE")]
    pub options: Option<PathBuf>,
    /// First epoch (inclusive), e.g. `2024-09-23 18:00:00`.
    #[arg(long)]
    pub start: Option<String>,
    /// Last epoch (inclusive).
    #[arg(long)]
    pub end: Option<String>,
    /// Write the epochs as a `.pos` file (the events go beside it as `NAME_events.pos` / `.csv`).
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Decimals of the seconds in the `.pos` time column (default 3; events 6).
    #[arg(long = "time-decimals", value_name = "N")]
    pub time_decimals: Option<usize>,
    /// Add velocity columns (north, east, up, m/s) to the epoch `.pos`.
    #[arg(long = "pos-velocity")]
    pub pos_velocity: bool,
    /// Also write the forward and backward tracks (`NAME_forward.pos`, `NAME_backward.pos`).
    #[arg(long = "keep-passes")]
    pub keep_passes: bool,
    /// GPT3 1° grid for `trop_model` gpt3 / vmf3.
    #[arg(long, value_name = "FILE")]
    pub gpt3: Option<PathBuf>,
    /// VMF3 orography.
    #[arg(long, value_name = "FILE")]
    pub orography: Option<PathBuf>,
    /// VMF3 6-hourly grids, repeatable.
    #[arg(long, value_name = "FILE")]
    pub vmf3: Vec<PathBuf>,
    /// Print JSON (summary, warnings, events; the epochs with `--json-epochs`).
    #[arg(long)]
    pub json: bool,
    /// Include every epoch in the JSON.
    #[arg(long = "json-epochs")]
    pub json_epochs: bool,
}

fn e(x: impl std::fmt::Display) -> CliError {
    CliError(x.to_string())
}

fn stem4(p: &Path) -> String {
    p.file_stem()
        .map(|s| {
            s.to_string_lossy()
                .chars()
                .take(4)
                .collect::<String>()
                .to_uppercase()
        })
        .unwrap_or_default()
}

fn beside(out: &Path, suffix: &str, ext: &str) -> PathBuf {
    let stem = out
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    out.with_file_name(format!("{stem}{suffix}.{ext}"))
}

/// Events from a plain text file: `WEEK SECONDS [LABEL]` or `GPS_SECONDS [LABEL]` per line; `#` starts a
/// comment; commas, tabs and blanks separate.
pub fn read_event_list(text: &str) -> Result<Vec<Event>, String> {
    let mut out = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let w: Vec<&str> = line
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .collect();
        let num = |s: &str| s.parse::<f64>().ok().filter(|v| v.is_finite());
        let bad = || {
            format!(
                "line {}: expected WEEK SECONDS [LABEL] or GPS_SECONDS [LABEL]",
                n + 1
            )
        };
        let a = num(w[0]).ok_or_else(bad)?;
        let (t, rest) = match w.get(1).and_then(|s| num(s)) {
            // a week is a small whole number, followed by seconds of the week
            Some(b) if a < 10_000.0 && a.fract() == 0.0 && (0.0..604_800.0).contains(&b) => {
                (a * 604_800.0 + b, &w[2..])
            }
            _ if a > 1.0e8 => (a, &w[1..]),
            _ => return Err(bad()),
        };
        out.push(Event {
            t,
            label: if rest.is_empty() {
                (out.len() + 1).to_string()
            } else {
                rest.join(" ")
            },
            offset: None,
        });
    }
    Ok(out)
}

fn three(s: &str, what: &str) -> Result<[f64; 3], CliError> {
    xyz(s, what)
}

struct Loaded {
    base: StationObs,
    rover: StationObs,
    events: Vec<Event>,
    occupations: Vec<Occupation>,
    warnings: Vec<String>,
    inputs: Vec<String>,
}

fn load(a: &PpkArgs, systems: &str) -> Result<Loaded, CliError> {
    let t = |s: &Option<String>| -> Result<Option<f64>, CliError> {
        s.as_deref()
            .map(|x| parse_user_time(x).map(ticks_to_seconds))
            .transpose()
            .map_err(CliError)
    };
    let (t0, t1) = (t(&a.start)?, t(&a.end)?);
    let mut warnings = Vec::new();
    let mut inputs = Vec::new();
    let mut parts = Vec::new();
    let mut events: Vec<Event> = Vec::new();
    let mut occupations = Vec::new();
    for p in &a.rover {
        let data = read_plain(p).map_err(e)?;
        let name = p
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        let chk = evrinex::check_times(&data).map_err(e)?;
        warnings.extend(chk.warnings.iter().map(|w| format!("{name}: {w}")));
        parts.push(read_obs_bytes_kinds(p, &data, &stem4(p), t0, t1, systems, "CLSD").map_err(e)?);
        if a.rinex_events {
            let ev = evrinex::events_from_bytes(&data).map_err(e)?;
            warnings.extend(ev.warnings.iter().map(|w| format!("{name}: {w}")));
            events.extend(ev.events);
        }
        if a.occupations {
            occupations.extend(evrinex::occupations_from_bytes(&data).map_err(e)?);
        }
        inputs.push(p.display().to_string());
    }
    let rover = StationObs::concat(&parts).ok_or_else(|| CliError("no rover file".into()))?;
    let base_data = read_plain(&a.base).map_err(e)?;
    // the base is only needed around the rover's span
    let pad = 120.0;
    let (b0, b1) = (
        rover.t.first().map(|x| x - pad),
        rover.t.last().map(|x| x + pad),
    );
    let base = read_obs_bytes_kinds(&a.base, &base_data, &stem4(&a.base), b0, b1, systems, "CLS")
        .map_err(e)?;
    inputs.push(a.base.display().to_string());
    if rover.t.is_empty() {
        return Err(CliError(
            "no rover observations in the window (check --start / --end)".into(),
        ));
    }
    if base.t.is_empty() {
        return Err(CliError(
            "the base file has no observations during the rover's session".into(),
        ));
    }
    let many = a.mrk.len() > 1;
    for p in &a.mrk {
        let mut m = mrk::read(p).map_err(e)?;
        m.resolve_weeks(rover.t[0]);
        let name = p
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        warnings.extend(m.warnings.iter().map(|w| format!("{name}: {w}")));
        let tag = p
            .file_stem()
            .map(|s| s.to_string_lossy().replace("_Timestamp", ""))
            .unwrap_or_default();
        for mut ev in m.events() {
            if many {
                ev.label = format!("{tag}:{}", ev.label);
            }
            if a.no_offsets {
                ev.offset = None;
            }
            events.push(ev);
        }
        inputs.push(p.display().to_string());
    }
    if let Some(p) = &a.events {
        let text = std::fs::read_to_string(p).map_err(|x| e(format!("{}: {x}", p.display())))?;
        events.extend(read_event_list(&text).map_err(|x| e(format!("{}: {x}", p.display())))?);
        inputs.push(p.display().to_string());
    }
    Ok(Loaded {
        base,
        rover,
        events,
        occupations,
        warnings,
        inputs,
    })
}

fn options(a: &PpkArgs) -> Result<PpkOptions, CliError> {
    let mut o = PpkOptions::default();
    if let Some(p) = &a.options {
        let text = std::fs::read_to_string(p).map_err(|x| e(format!("{}: {x}", p.display())))?;
        o = o
            .with_json(&text)
            .map_err(|x| e(format!("--options: {x}")))?;
    }
    if let Some(s) = &a.systems {
        o.systems.clone_from(s);
    }
    if let Some(v) = a.elmask {
        o.elmask = v;
    }
    if let Some(d) = a.direction.as_deref() {
        o.direction = match d {
            "forward" => Direction::Forward,
            "backward" => Direction::Backward,
            "combined" => Direction::Combined,
            x => {
                return Err(CliError(format!(
                    "--direction: forward, backward or combined, not {x}"
                )));
            }
        };
    }
    if let Some(m) = a.interp.as_deref() {
        o.event_interp = match m {
            "linear2" | "linear" | "2" => EventInterp::Linear2,
            "cubic4" | "cubic" | "4" => EventInterp::Cubic4,
            "hermite2" | "hermite" => EventInterp::Hermite2,
            x => {
                return Err(CliError(format!(
                    "--interp: linear2, cubic4 or hermite2, not {x}"
                )));
            }
        };
    }
    if let Some(n) = a.max_missing {
        o.event_max_missing = n;
    }
    if a.no_fix {
        o.ar = ppk::ArMode::Off;
    }
    if a.keep_passes {
        o.keep_passes = true;
    }
    Ok(o)
}

fn summary_text(res: &PpkResult, l: &Loaded, secs: f64) -> String {
    let s = &res.summary;
    let mut out = String::new();
    let _ = writeln!(
        out,
        "{} -> {}: {} epochs at {:.2} s, fixed {:.1}% ({} fixed, {} float, {} single, {} none)",
        l.base.name,
        l.rover.name,
        s.n_epochs,
        s.interval_s,
        s.fix_percent,
        s.n_fixed,
        s.n_float,
        s.n_single,
        s.n_none
    );
    if let (Some(a), Some(b)) = (res.epochs.first(), res.epochs.last()) {
        let _ = writeln!(
            out,
            "span {} - {} GPST",
            gpst_text(a.t_tag, 1),
            gpst_text(b.t_tag, 1)
        );
    }
    let _ = writeln!(out, "signals (rover/base): {}", s.signals.join("; "));
    let _ = writeln!(
        out,
        "baseline {:.1} - {:.1} m; new ambiguities {} (first appearances, breaks, slips); slips repaired {}",
        s.baseline_m[0], s.baseline_m[1], s.slips, s.slips_repaired
    );
    if s.fix_percent_forward > 0.0 || s.fix_percent_backward > 0.0 {
        let _ = writeln!(
            out,
            "forward fixed {:.1}%, backward {:.1}%; both fixed at {} epochs: difference rms E {:.1} N {:.1} U {:.1} mm, \
             max {:.1} mm; conflicts {}",
            s.fix_percent_forward,
            s.fix_percent_backward,
            s.both_fixed,
            1e3 * s.fb_rms_enu[0],
            1e3 * s.fb_rms_enu[1],
            1e3 * s.fb_rms_enu[2],
            1e3 * s.fb_max,
            s.fb_conflicts
        );
    }
    if !res.events.is_empty() {
        let fixed = res
            .events
            .iter()
            .filter(|ev| ev.status != EventStatus::NoSolution && ev.fix == ppk::FixStatus::Fixed)
            .count();
        let _ = writeln!(
            out,
            "events: {} ({} with a position, {} of them fixed; {} across a gap; {} without)",
            res.events.len(),
            s.events_ok + s.events_gap,
            fixed,
            s.events_gap,
            s.events_failed
        );
    }
    for oc in &res.occupations {
        let _ = writeln!(
            out,
            "occupation {}: {} ({} of {} epochs) {:.9} {:.9} {:.4}  σ E {:.4} N {:.4} U {:.4} m",
            oc.name,
            oc.fix.name(),
            oc.n,
            oc.n_all,
            oc.llh[0],
            oc.llh[1],
            oc.llh[2],
            oc.sigma_enu[0],
            oc.sigma_enu[1],
            oc.sigma_enu[2]
        );
    }
    for w in &res.warnings {
        let _ = writeln!(out, "warning: {w}");
    }
    for n in &res.notes {
        let _ = writeln!(out, "note: {n}");
    }
    let _ = write!(out, "time: {secs:.2} s");
    out
}

/// Runs `gnss ppk` and returns the text or JSON to print.
#[allow(clippy::too_many_lines)]
pub fn ppk(a: &PpkArgs) -> Result<String, CliError> {
    let t_all = Instant::now();
    let opts = options(a)?;
    let l = load(a, &opts.systems)?;
    let nav = read_nav(&a.nav);
    for (name, why) in &nav.bad {
        eprintln!("note: navigation file {name} not readable ({why}): left out");
    }
    if a.nav.is_empty() && a.sp3.is_empty() {
        return Err(CliError(
            "no orbits: give --nav (broadcast) and / or --sp3".into(),
        ));
    }
    let fcn: BTreeMap<Sat, i32> = nav.glo_fcn();
    let mut srcs = Vec::new();
    if !a.sp3.is_empty() {
        srcs.push(OrbitSource::Precise(
            Precise::open(&a.sp3, &a.clk).map_err(e)?,
        ));
    }
    srcs.push(OrbitSource::Broadcast(Broadcast::new(nav)));
    let orb = Orbits::new(srcs);
    let mut atx = a.antex.as_deref().map(Antex::open).transpose().map_err(e)?;
    let gpt3 = a
        .gpt3
        .as_deref()
        .map(Gpt3Grid::open)
        .transpose()
        .map_err(e)?;
    let orog = a
        .orography
        .as_deref()
        .map(load_orography)
        .transpose()
        .map_err(e)?;
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
            rows: load_vmf3(p).map_err(e)?,
        });
    }
    let tropo = TropoInputs {
        gpt3: gpt3.as_ref(),
        orography: orog.as_deref(),
        vmf3: &vmf3,
    };
    let mut extra = l.warnings.clone();
    let marker_xyz = if let Some(s) = &a.base_xyz {
        three(s, "--base-xyz")?
    } else if let Some(s) = &a.base_llh {
        let v = three(s, "--base-llh")?;
        if v[0].abs() > 90.0 || v[1].abs() > 360.0 {
            return Err(CliError(
                "--base-llh: LAT,LON,H in degrees and metres expected".into(),
            ));
        }
        llh_to_ecef(v[0].to_radians(), v[1].to_radians(), v[2])
    } else {
        extra.push(
            "no base coordinate given: the approximate position in the base file's header is used. The \
             positions are then only as good as that number (often metres off), though their shape is right"
                .into(),
        );
        l.base.approx
    };
    let delta_hen = a
        .base_hen
        .as_deref()
        .map(|s| three(s, "--base-hen"))
        .transpose()?;
    let mut input = PpkInput {
        base: &l.base,
        rover: &l.rover,
        orbits: &orb,
        fcn: &fcn,
        base_station: BaseStation {
            marker_xyz,
            delta_hen,
        },
        antex: atx.as_mut(),
        tropo: &tropo,
        events: &l.events,
        occupations: &l.occupations,
    };
    let mut res = ppk::solve(&mut input, &opts, &mut |_| true).map_err(e)?;
    extra.append(&mut res.warnings);
    res.warnings = extra;

    let mut written: Vec<String> = Vec::new();
    if let Some(out) = &a.output {
        let mut h = PosHeader::new(&opts);
        h.inputs.clone_from(&l.inputs);
        h.inputs
            .extend(a.nav.iter().map(|p| p.display().to_string()));
        h.time_decimals = a.time_decimals.unwrap_or(3);
        h.velocity = a.pos_velocity;
        let mut put = |p: &Path, text: String| -> Result<(), CliError> {
            std::fs::write(p, text).map_err(|x| e(format!("{}: {x}", p.display())))?;
            written.push(p.display().to_string());
            Ok(())
        };
        put(out, epochs_pos(&res, &h))?;
        if opts.keep_passes {
            put(
                &beside(out, "_forward", "pos"),
                ppk::pos::track_pos(&res, &res.forward, &h),
            )?;
            put(
                &beside(out, "_backward", "pos"),
                ppk::pos::track_pos(&res, &res.backward, &h),
            )?;
        }
        if !res.events.is_empty() {
            let mut he = h.clone();
            he.time_decimals = a.time_decimals.unwrap_or(6);
            he.velocity = false;
            put(&beside(out, "_events", "pos"), events_pos(&res, &he))?;
            put(&beside(out, "_events", "csv"), events_csv(&res))?;
        }
    }
    let elapsed = t_all.elapsed().as_secs_f64();
    if a.json {
        let mut v = json!({
            "summary": res.summary,
            "warnings": res.warnings,
            "notes": res.notes,
            "base_marker": res.base_marker,
            "base_arp": res.base_arp,
            "events": res.events,
            "occupations": res.occupations,
            "written": written,
            "seconds": elapsed,
        });
        if a.json_epochs {
            v["epochs"] = serde_json::to_value(&res.epochs).map_err(e)?;
        }
        // NaN is not JSON: serde_json writes null for it
        return serde_json::to_string_pretty(&v).map_err(e);
    }
    let mut s = summary_text(&res, &l, elapsed);
    for w in &written {
        let _ = write!(s, "\nwrote {w}");
    }
    Ok(s)
}

/// `gnss ppp --kinematic`: the kinematic precise point solution of an observation file already read, with
/// the events and outputs of `gnss ppk`.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn ppp_kinematic(
    a: &crate::ppp::PppArgs,
    st: &StationObs,
    orb: &Orbits,
    antex: Option<&mut Antex>,
    popts: &geomaticist_gnss_ppp::PppOptions,
    fcn: &BTreeMap<Sat, i32>,
    osb: &geomaticist_gnss_ppp::solve::Osb,
    erp: Option<&geomaticist_gnss_products::bias::Erp>,
    blq: Option<&[[f64; 11]; 6]>,
    tropo: &TropoInputs,
    x0: Option<[f64; 3]>,
) -> Result<String, CliError> {
    use geomaticist_gnss_ppp::kinematic::{KinematicInput, KinematicOptions, solve};
    let t_all = Instant::now();
    let mut ko = KinematicOptions::default();
    if let Some(p) = &a.kinematic_options {
        let text = std::fs::read_to_string(p).map_err(|x| e(format!("{}: {x}", p.display())))?;
        ko = ko
            .with_json(&text)
            .map_err(|x| e(format!("--kinematic-options: {x}")))?;
    }
    ko.ppp = popts.clone();
    if let Some(d) = a.direction.as_deref() {
        ko.direction = match d {
            "forward" => Direction::Forward,
            "backward" => Direction::Backward,
            "combined" => Direction::Combined,
            x => {
                return Err(CliError(format!(
                    "--direction: forward, backward or combined, not {x}"
                )));
            }
        };
    }
    if let Some(m) = a.interp.as_deref() {
        ko.event_interp = match m {
            "linear2" | "linear" => EventInterp::Linear2,
            "cubic4" | "cubic" => EventInterp::Cubic4,
            "hermite2" | "hermite" => EventInterp::Hermite2,
            x => {
                return Err(CliError(format!(
                    "--interp: linear2, cubic4 or hermite2, not {x}"
                )));
            }
        };
    }
    if let Some(n) = a.max_missing {
        ko.event_max_missing = n;
    }
    ko.keep_passes = a.keep_passes;
    let mut warnings = Vec::new();
    let mut events: Vec<Event> = Vec::new();
    if a.rinex_events {
        let ev = evrinex::read_events(&a.obs).map_err(e)?;
        warnings.extend(ev.warnings);
        events.extend(ev.events);
    }
    let many = a.mrk.len() > 1;
    for p in &a.mrk {
        let mut m = mrk::read(p).map_err(e)?;
        m.resolve_weeks(st.t[0]);
        warnings.append(&mut m.warnings);
        let tag = p
            .file_stem()
            .map(|s| s.to_string_lossy().replace("_Timestamp", ""))
            .unwrap_or_default();
        for mut ev in m.events() {
            if many {
                ev.label = format!("{tag}:{}", ev.label);
            }
            events.push(ev);
        }
    }
    if let Some(p) = &a.events {
        let text = std::fs::read_to_string(p).map_err(|x| e(format!("{}: {x}", p.display())))?;
        events.extend(read_event_list(&text).map_err(|x| e(format!("{}: {x}", p.display())))?);
    }
    let mut input = KinematicInput {
        obs: st,
        orbits: orb,
        antex,
        fcn,
        osb,
        erp,
        blq,
        tropo,
        x0,
        events: &events,
        occupations: &[],
    };
    let mut res = solve(&mut input, &ko, &mut |_| true).map_err(e)?;
    warnings.append(&mut res.warnings);
    res.warnings = warnings;
    let mut written: Vec<String> = Vec::new();
    if let Some(out) = &a.output {
        let mut h = PosHeader::new(&PpkOptions {
            systems: popts.systems.clone(),
            elmask: popts.elmask,
            direction: ko.direction,
            ar: ppk::ArMode::Off,
            trop_model: popts.trop_model.clone(),
            ..PpkOptions::default()
        });
        h.inputs = vec![a.obs.display().to_string()];
        "PPP Kinematic".clone_into(&mut h.mode);
        let mut put = |p: &Path, text: String| -> Result<(), CliError> {
            std::fs::write(p, text).map_err(|x| e(format!("{}: {x}", p.display())))?;
            written.push(p.display().to_string());
            Ok(())
        };
        put(out, epochs_pos(&res, &h))?;
        if ko.keep_passes {
            put(
                &beside(out, "_forward", "pos"),
                ppk::pos::track_pos(&res, &res.forward, &h),
            )?;
            put(
                &beside(out, "_backward", "pos"),
                ppk::pos::track_pos(&res, &res.backward, &h),
            )?;
        }
        if !res.events.is_empty() {
            let mut he = h.clone();
            he.time_decimals = 6;
            put(&beside(out, "_events", "pos"), events_pos(&res, &he))?;
            put(&beside(out, "_events", "csv"), events_csv(&res))?;
        }
    }
    let elapsed = t_all.elapsed().as_secs_f64();
    if a.json {
        return serde_json::to_string_pretty(&json!({
            "summary": res.summary,
            "warnings": res.warnings,
            "notes": res.notes,
            "events": res.events,
            "written": written,
            "seconds": elapsed,
        }))
        .map_err(e);
    }
    let s = &res.summary;
    let mut out = format!(
        "{}: kinematic PPP, {} epochs at {:.2} s, {} with a position, {} arcs",
        st.name, s.n_epochs, s.interval_s, s.n_float, s.slips
    );
    let mut sig: Vec<f64> = res
        .epochs
        .iter()
        .filter(|x| x.has_position())
        .map(|x| {
            let g = x.sigma_enu();
            (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt()
        })
        .collect();
    if !sig.is_empty() {
        sig.sort_by(f64::total_cmp);
        let _ = write!(
            out,
            "\n3-D sigma: median {:.3} m, largest {:.3} m",
            sig[sig.len() / 2],
            sig[sig.len() - 1]
        );
    }
    if !res.events.is_empty() {
        let _ = write!(
            out,
            "\nevents: {} ({} with a position, {} across a gap, {} without)",
            res.events.len(),
            s.events_ok + s.events_gap,
            s.events_gap,
            s.events_failed
        );
    }
    for w in &res.warnings {
        let _ = write!(out, "\nwarning: {w}");
    }
    for n in &res.notes {
        let _ = write!(out, "\nnote: {n}");
    }
    let _ = write!(out, "\ntime: {elapsed:.2} s");
    for w in &written {
        let _ = write!(out, "\nwrote {w}");
    }
    Ok(out)
}

/// `gnss events ...`.
#[derive(Debug, Subcommand)]
pub enum EventsCommand {
    /// List the events of a DJI Timestamp.MRK file or of a RINEX observation file (flag 5).
    List(EventsListArgs),
    /// Check the epoch times of an observation file for known recorder glitches.
    Check(EventsCheckArgs),
    /// Write a copy of an observation file with an event record (flag 5) for each event.
    Insert(EventsInsertArgs),
}

/// `gnss events list`.
#[derive(Debug, Args)]
pub struct EventsListArgs {
    /// A `Timestamp.MRK` file or a RINEX observation file.
    pub file: PathBuf,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// `gnss events check`.
#[derive(Debug, Args)]
pub struct EventsCheckArgs {
    /// RINEX observation file.
    pub file: PathBuf,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

/// `gnss events insert`.
#[derive(Debug, Args)]
pub struct EventsInsertArgs {
    /// RINEX observation file (the rover's).
    pub file: PathBuf,
    /// DJI `Timestamp.MRK` file(s) with the events, repeatable.
    #[arg(long, value_name = "FILE")]
    pub mrk: Vec<PathBuf>,
    /// A text file of events (`WEEK SECONDS [LABEL]` or `GPS_SECONDS [LABEL]` per line).
    #[arg(long = "events", value_name = "FILE")]
    pub events: Option<PathBuf>,
    /// The file to write.
    #[arg(short, long)]
    pub output: PathBuf,
    /// Print JSON.
    #[arg(long)]
    pub json: bool,
}

fn is_mrk(p: &Path) -> bool {
    p.extension().is_some_and(|x| x.eq_ignore_ascii_case("mrk"))
}

fn events_list(a: &EventsListArgs) -> Result<String, CliError> {
    if is_mrk(&a.file) {
        let m = mrk::read(&a.file).map_err(e)?;
        if a.json {
            return serde_json::to_string_pretty(
                &json!({"records": m.records, "warnings": m.warnings}),
            )
            .map_err(e);
        }
        let mut s = format!("{}: {} photos", a.file.display(), m.records.len());
        for (q, n) in m.quality_counts() {
            let name = match q {
                50 => "fixed",
                34 => "float",
                16 => "single point",
                _ => "other",
            };
            let _ = write!(s, "\n  aircraft position quality {q} ({name}): {n}");
        }
        let ev = m.events();
        if let (Some(f), Some(l)) = (ev.first(), ev.last()) {
            let _ = write!(s, "\n  {} - {} GPST", gpst_text(f.t, 3), gpst_text(l.t, 3));
        }
        for w in &m.warnings {
            let _ = write!(s, "\n  warning: {w}");
        }
        return Ok(s);
    }
    let r = evrinex::read_events(&a.file).map_err(e)?;
    if a.json {
        return serde_json::to_string_pretty(&json!({"events": r.events, "warnings": r.warnings}))
            .map_err(e);
    }
    let mut s = format!(
        "{}: {} event records (flag 5)",
        a.file.display(),
        r.events.len()
    );
    for ev in &r.events {
        let _ = write!(s, "\n  {:>5}  {}", ev.label, gpst_text(ev.t, 7));
    }
    for w in &r.warnings {
        let _ = write!(s, "\n  warning: {w}");
    }
    Ok(s)
}

fn events_check(a: &EventsCheckArgs) -> Result<String, CliError> {
    let data = read_plain(&a.file).map_err(e)?;
    let c = evrinex::check_times(&data).map_err(e)?;
    if a.json {
        return serde_json::to_string_pretty(&json!({
            "epochs": c.epochs,
            "before_1990": c.before_1990,
            "displaced": c.displaced,
            "not_later": c.not_later,
            "first": c.span.map(|s| fmt(s.0, 7)),
            "last": c.span.map(|s| fmt(s.1, 7)),
            "warnings": c.warnings,
        }))
        .map_err(e);
    }
    let mut s = format!("{}: {} epochs", a.file.display(), c.epochs);
    if let Some((f, l)) = c.span {
        let _ = write!(s, ", {} - {}", fmt(f, 1), fmt(l, 1));
    }
    if c.warnings.is_empty() {
        s.push_str("\n  epoch times look right");
    }
    for w in &c.warnings {
        let _ = write!(s, "\n  warning: {w}");
    }
    Ok(s)
}

fn events_insert(a: &EventsInsertArgs) -> Result<String, CliError> {
    let data = read_plain(&a.file).map_err(e)?;
    let chk = evrinex::check_times(&data).map_err(e)?;
    let mut warnings = chk.warnings.clone();
    let first = chk
        .span
        .map(|s| ticks_to_seconds(s.0))
        .ok_or_else(|| CliError("the file has no observation epochs".into()))?;
    let mut events = Vec::new();
    for p in &a.mrk {
        let mut m = mrk::read(p).map_err(e)?;
        m.resolve_weeks(first);
        warnings.append(&mut m.warnings);
        events.extend(m.events());
    }
    if let Some(p) = &a.events {
        let text = std::fs::read_to_string(p).map_err(|x| e(format!("{}: {x}", p.display())))?;
        events.extend(read_event_list(&text).map_err(|x| e(format!("{}: {x}", p.display())))?);
    }
    if events.is_empty() {
        return Err(CliError("no events: give --mrk or --events".into()));
    }
    let mut ins = evrinex::insert_events(&data, &events).map_err(e)?;
    warnings.append(&mut ins.warnings);
    std::fs::write(&a.output, &ins.data).map_err(|x| e(format!("{}: {x}", a.output.display())))?;
    if a.json {
        return serde_json::to_string_pretty(&json!({
            "path": a.output.display().to_string(),
            "inserted": ins.inserted,
            "before": ins.before,
            "after": ins.after,
            "existing": ins.existing,
            "warnings": warnings,
        }))
        .map_err(e);
    }
    let mut s = format!(
        "wrote {}: {} of {} events inserted",
        a.output.display(),
        ins.inserted,
        events.len()
    );
    for w in &warnings {
        let _ = write!(s, "\n  warning: {w}");
    }
    Ok(s)
}

/// Runs `gnss events ...`.
pub fn events(c: &EventsCommand) -> Result<String, CliError> {
    match c {
        EventsCommand::List(a) => events_list(a),
        EventsCommand::Check(a) => events_check(a),
        EventsCommand::Insert(a) => events_insert(a),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn event_lists() {
        let ev = read_event_list(
            "# photos\n2300 345600.25 IMG_0001\n2300,345602.5\n1391385600.75 a b\n",
        )
        .unwrap();
        assert_eq!(ev.len(), 3);
        assert!((ev[0].t - (2300.0 * 604_800.0 + 345_600.25)).abs() < 1e-6);
        assert_eq!(ev[0].label, "IMG_0001");
        assert_eq!(ev[1].label, "2");
        assert_eq!(ev[2].t, 1_391_385_600.75);
        assert_eq!(ev[2].label, "a b");
        assert!(read_event_list("hello\n").is_err());
        assert!(read_event_list("12.5\n").is_err());
    }

    #[test]
    fn output_names() {
        assert_eq!(
            beside(Path::new("out/flight.pos"), "_events", "csv"),
            Path::new("out/flight_events.csv")
        );
        assert!(is_mrk(Path::new("x_Timestamp.MRK")));
    }
}
