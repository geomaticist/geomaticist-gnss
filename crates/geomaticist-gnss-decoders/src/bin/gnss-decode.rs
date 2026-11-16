//! `gnss-decode <raw> -o <out> [--json] [--runpkr00 <exe>]`: one raw log (.T04, .T01/.T00/.T02/.R00, Leica .m00, Topcon .tps)
//! to RINEX. `gnss convert` runs this program for raw logs, so the `gnss` CLI itself never links the decoders
//! (rule R2).
//!
//! A packed Trimble log is read by the built-in reader; `--runpkr00` names the user's installed copy of
//! Trimble's runpkr00 for a file that reader doesn't know (default: `$RUNPKR00`, beside this program, `PATH`).
//!
//! `-o` names the observation file (the navigation files go beside it under the decoder's names), or an existing
//! folder (every file under the decoder's names, the usual RINEX short names).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use geomaticist_gnss_decoders::{Converted, DecodeError, RawKind, raw_kind};
use serde_json::json;

const USAGE: &str = "usage: gnss-decode <raw log> -o <out.rnx | folder> [--json] [--runpkr00 <exe>]      [--date YYYY-MM-DD] [--nav <GLONASS navigation file>]   (the last two for RTCM 3 logs)";

struct Args {
    raw: PathBuf,
    out: PathBuf,
    json: bool,
    runpkr00: Option<PathBuf>,
    /// RTCM 3 logs: a date near the data (the stream has no week number).
    date: Option<String>,
    /// RTCM 3 logs: a navigation file with the GLONASS frequency channels.
    nav: Vec<PathBuf>,
}

fn parse(mut it: impl Iterator<Item = String>) -> Result<Args, String> {
    let (mut raw, mut out, mut json, mut runpkr00) = (None, None, false, None);
    let (mut date, mut nav) = (None, Vec::new());
    while let Some(a) = it.next() {
        match a.as_str() {
            "-o" | "--output" => out = Some(PathBuf::from(it.next().ok_or("-o needs a path")?)),
            "--json" => json = true,
            "--runpkr00" => {
                runpkr00 = Some(PathBuf::from(it.next().ok_or("--runpkr00 needs a path")?));
            }
            "--date" => date = Some(it.next().ok_or("--date needs YYYY-MM-DD")?),
            "--nav" => nav.push(PathBuf::from(it.next().ok_or("--nav needs a path")?)),
            "-h" | "--help" => return Err(USAGE.into()),
            s if s.starts_with('-') => return Err(format!("unknown option {s}\n{USAGE}")),
            _ if raw.is_none() => raw = Some(PathBuf::from(a)),
            _ => return Err(format!("one raw log at a time\n{USAGE}")),
        }
    }
    Ok(Args {
        raw: raw.ok_or(USAGE)?,
        out: out.ok_or(USAGE)?,
        json,
        runpkr00,
        date,
        nav,
    })
}

fn decode(a: &Args, dir: &Path, log: &mut dyn FnMut(&str)) -> Result<Converted, DecodeError> {
    let ext = a
        .raw
        .extension()
        .map_or_default(|e| e.to_string_lossy().into_owned());
    if geomaticist_gnss_decoders::rtcm3::is_rtcm_ext(&ext) {
        let near = match &a.date {
            None => None,
            Some(d) => Some(
                geomaticist_gnss_core::time::parse_user_time(&format!("{d} 12:00:00"))
                    .map(geomaticist_gnss_core::time::ticks_to_seconds)
                    .map_err(|e| DecodeError::Format(format!("--date: {e}")))?,
            ),
        };
        let fcn = geomaticist_gnss_rinex::nav::read_nav(&a.nav).glo_fcn();
        return geomaticist_gnss_decoders::rtcm3::convert(&a.raw, dir, log, near, &fcn);
    }
    match raw_kind(&a.raw) {
        Some(RawKind::Packed) => {
            geomaticist_gnss_decoders::t01::convert(&a.raw, dir, log, a.runpkr00.as_deref())
        }
        _ => geomaticist_gnss_decoders::convert(&a.raw, dir, log),
    }
}

fn run(a: &Args) -> Result<String, String> {
    let to_dir = a.out.is_dir();
    let dir = if to_dir {
        a.out.clone()
    } else {
        a.out
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
    };
    let mut lines = Vec::new();
    let t = Instant::now();
    let mut res = decode(a, &dir, &mut |m| lines.push(m.to_owned())).map_err(|e| e.to_string())?;
    if !to_dir && res.obs != a.out {
        if a.out.exists() {
            std::fs::remove_file(&a.out).map_err(|e| format!("{}: {e}", a.out.display()))?;
        }
        std::fs::rename(&res.obs, &a.out).map_err(|e| format!("{}: {e}", a.out.display()))?;
        res.obs.clone_from(&a.out);
    }
    let secs = t.elapsed().as_secs_f64();
    if a.json {
        let v = json!({
            "obs": res.obs.display().to_string(),
            "nav": res.nav.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
            "warnings": res.warnings,
            "info": res.info,
            "log": lines,
            "seconds": secs,
        });
        return serde_json::to_string_pretty(&v).map_err(|e| e.to_string());
    }
    let mut s = lines.join("\n");
    for w in &res.warnings {
        let _ = write!(s, "\nwarning: {w}");
    }
    let _ = write!(s, "\nwrote {}", res.obs.display());
    for n in &res.nav {
        let _ = write!(s, "\nwrote {}", n.display());
    }
    let _ = write!(s, "\n{secs:.2} s");
    Ok(s.trim_start().to_owned())
}

fn main() -> ExitCode {
    let args = match parse(std::env::args().skip(1)) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(s) => {
            println!("{s}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn args() {
        let a = parse(
            ["x.T04", "-o", "y.26o", "--json"]
                .iter()
                .map(|s| (*s).to_owned()),
        )
        .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(a.raw, PathBuf::from("x.T04"));
        assert!(a.json);
        assert!(parse(["x.T04"].iter().map(|s| (*s).to_owned())).is_err());
    }
}
