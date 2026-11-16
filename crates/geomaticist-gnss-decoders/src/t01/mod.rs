//! Packed Trimble raw logs (.T01 ...) to RINEX 2.11: the log is unpacked into
//! Trimble's record chain, by the built-in reader ([`native`]) or, for a file it doesn't know, by Trimble's
//! runpkr00 when that program is on the machine ([`runpkr`]); we read the records ([`records`]) and write
//! observation and navigation files ([`rinex2`], [`nav2`]). A clean-room decoder written from our own
//! receiver logs.

pub mod native;
pub mod nav2;
pub mod records;
pub mod rinex2;
pub mod runpkr;

use std::path::Path;

use serde_json::json;

use crate::{Converted, DecodeError};

/// PGM of the files this decoder writes.
pub const PROGRAM: &str = "geomaticist-gnss T01";
/// RUN BY.
pub const RUN_BY: &str = "geomaticist-gnss";

/// (IGS name, dome, model name) for an antenna code, from the T04 decoder's table (same Trimble codes).
#[must_use]
pub fn antenna_names(code: Option<i64>) -> (&'static str, &'static str, &'static str) {
    code.and_then(crate::antenna_lib::antenna)
        .map_or(("", "", ""), |a| (a.0, a.1, a.3))
}

fn py_opt_float(v: Option<f64>) -> String {
    v.map_or_else(|| "?".to_owned(), |h| format!("{h:?}"))
}

/// Whether a packed log can be turned into RINEX on this machine: the built-in reader knows it, or runpkr00 is
/// there to try.
#[must_use]
pub fn can_convert(raw: &Path) -> bool {
    use std::io::Read as _;
    let mut head = Vec::new();
    let known = std::fs::File::open(raw)
        .and_then(|f| f.take(16 * 1024).read_to_end(&mut head))
        .is_ok()
        && native::recognises(&head);
    known || runpkr::find_runpkr00(None).is_some()
}

fn names(raw: &Path) -> (String, String) {
    (
        raw.file_name()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned()),
        raw.file_stem()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned()),
    )
}

/// Unpacks `raw` and writes RINEX 2.11 into `outdir`. Names follow the usual
/// RINEX 2 short names: 12341230.T01 -> 12341230.15o, .15n (GPS), .15g (GLONASS) for a log of 2015.
///
/// The built-in reader is tried first. runpkr00 (`exe`, or the one found on the machine) is used instead when
/// the built-in reader doesn't know the file, or lost epochs in it and the program is there to try; when it
/// isn't, what the built-in reader got is written, with a warning.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    exe: Option<&Path>,
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let (name, _) = names(raw);
    let found = || runpkr::find_runpkr00(exe);
    let now = rinex2::pgm_date(&crate::now_utc());
    match native::unpack(&data) {
        Ok(u) if u.complete() => write_native(raw, &u, outdir, log, &now),
        Ok(u) => {
            if let Some(exe) = found() {
                log(&format!(
                    "{name}: the built-in reader could not read all of it ({}); using runpkr00",
                    u.notes.join("; ")
                ));
                if let Ok(c) = convert_runpkr00(raw, outdir, log, Some(&exe)) {
                    return Ok(c);
                }
            }
            write_native(raw, &u, outdir, log, &now)
        }
        Err(e) => match found() {
            Some(exe) => convert_runpkr00(raw, outdir, log, Some(&exe)),
            None => Err(match e {
                DecodeError::Format(m) => DecodeError::Runpkr(format!(
                    "{name}: {m}, and Trimble's runpkr00 isn't on this machine to try instead"
                )),
                other => other,
            }),
        },
    }
}

/// [`convert`] with the built-in reader only.
pub fn convert_native(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    let now = rinex2::pgm_date(&crate::now_utc());
    convert_native_with(raw, outdir, log, &now)
}

/// [`convert_native`] with the PGM date text given (`%d-%b-%y %H:%M UTC`).
pub fn convert_native_with(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let (name, _) = names(raw);
    let u = native::unpack(&data).map_err(|e| match e {
        DecodeError::Format(m) => DecodeError::Format(format!("{name}: {m}")),
        other => other,
    })?;
    write_native(raw, &u, outdir, log, date)
}

fn write_native(
    raw: &Path,
    u: &native::Unpacked,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
) -> Result<Converted, DecodeError> {
    let (name, stem) = names(raw);
    let comment = format!("Unpacked from {name} by {RUN_BY}");
    let blob = u.records.as_slice();
    let mut c = write(&[blob], blob, &name, &stem, outdir, log, date, &comment, "")?;
    if !u.complete() {
        let mut w = format!(
            "{name}: part of the log could not be read ({} epochs read, {} lost",
            u.epochs, u.lost_epochs
        );
        if let Some(n) = u.notes.first() {
            w.push_str("; ");
            w.push_str(n);
        }
        w.push_str(").");
        c.warnings.push(w);
    }
    Ok(c)
}

/// [`convert`] with Trimble's runpkr00 only (`exe`, or the one found on the machine).
pub fn convert_runpkr00(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    exe: Option<&Path>,
) -> Result<Converted, DecodeError> {
    let (dat, tgd) = runpkr::unpack(raw, exe)?;
    let (name, stem) = names(raw);
    let now = rinex2::pgm_date(&crate::now_utc());
    convert_unpacked(&[&dat, &tgd], &tgd, &name, &stem, outdir, log, &now)
}

/// The writer half of [`convert_runpkr00`]: `header_blobs` (the .DAT and / or .TGD) and the .TGD already
/// unpacked, with the PGM date text given (`%d-%b-%y %H:%M UTC`).
pub fn convert_unpacked(
    header_blobs: &[&[u8]],
    tgd: &[u8],
    raw_name: &str,
    stem: &str,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
) -> Result<Converted, DecodeError> {
    let comment = format!("Unpacked from {raw_name} by Trimble runpkr00");
    write(
        header_blobs,
        tgd,
        raw_name,
        stem,
        outdir,
        log,
        date,
        &comment,
        " (unpacked by runpkr00)",
    )
}

/// Unpacked records to RINEX: `comment` goes in the observation header, `how` ends the log line.
#[allow(clippy::too_many_arguments)]
fn write(
    header_blobs: &[&[u8]],
    tgd: &[u8],
    raw_name: &str,
    stem: &str,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
    comment: &str,
    how: &str,
) -> Result<Converted, DecodeError> {
    std::fs::create_dir_all(outdir).map_err(|e| DecodeError::io(outdir, &e))?;
    let fmt_err = |e: String| DecodeError::Format(format!("{raw_name}: {e}"));
    let hdr = records::read_header(header_blobs).map_err(fmt_err)?;
    let (ant, dome, model) = antenna_names(hdr.antenna_code);
    let mut warnings = Vec::new();
    if ant.is_empty() {
        warnings.push(format!(
            "{raw_name}: the antenna isn't known (code {}); set the antenna type before processing.",
            hdr.antenna_code.map_or_else(|| "None".to_owned(), |c| c.to_string())
        ));
    }
    let first = records::epochs(tgd)
        .map_err(fmt_err)?
        .into_iter()
        .next()
        .ok_or_else(|| DecodeError::Format(format!("{raw_name} has no measurement epochs")))?;
    let t0 = rinex2::epoch_time(&first);
    let yy = t0.year.rem_euclid(100);
    let obs = outdir.join(format!("{stem}.{yy:02}o"));
    let opt = rinex2::ObsOptions {
        program: PROGRAM,
        run_by: RUN_BY,
        date,
        comment,
        antenna: (ant, dome),
        line_end: "\r\n",
    };
    let (_head, text) = rinex2::write_obs(tgd, &hdr, &opt).map_err(fmt_err)?;
    crate::write_atomic(&obs, &text)?;
    let leap = rinex2::leap_seconds(&t0);
    let mut nav = Vec::new();
    if let Some(t) = nav2::gps_nav(tgd, PROGRAM, RUN_BY, date, "\r\n").map_err(fmt_err)? {
        let p = outdir.join(format!("{stem}.{yy:02}n"));
        crate::write_atomic(&p, &t)?;
        nav.push(p);
    }
    if let Some(t) = nav2::glonass_nav(tgd, leap, PROGRAM, RUN_BY, date, "\r\n").map_err(fmt_err)? {
        let p = outdir.join(format!("{stem}.{yy:02}g"));
        crate::write_atomic(&p, &t)?;
        nav.push(p);
    }
    log(&format!(
        "{raw_name}: {} #{} fw {}, antenna {}, ARP height {} m, marker '{}'{how}",
        hdr.model,
        hdr.serial,
        hdr.firmware,
        if ant.is_empty() { "?" } else { ant },
        py_opt_float(hdr.arp_height),
        hdr.marker
    ));
    let info = json!({
        "receiver": hdr.model,
        "serial": hdr.serial,
        "firmware": hdr.firmware,
        "antenna": ant,
        "antenna_model": model,
        "antenna_code": hdr.antenna_code,
        "arp_height": hdr.arp_height,
        "marker": hdr.marker,
    });
    Ok(Converted {
        obs,
        nav,
        warnings,
        info,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Dir(std::path::PathBuf);

    impl Dir {
        fn new(tag: &str) -> Self {
            let d = std::env::temp_dir().join(format!("t01_mod_{tag}_{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&d);
            std::fs::create_dir_all(&d).unwrap();
            Dir(d)
        }
    }

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn a_synthetic_log_becomes_rinex_without_runpkr00() {
        let dir = Dir::new("ok");
        let raw = dir.0.join("12341230.T01");
        let log = native::synth::log(5);
        std::fs::write(&raw, native::synth::file(&log.stream)).unwrap();
        assert!(can_convert(&raw));
        let mut lines = Vec::new();
        let c = convert_native_with(
            &raw,
            &dir.0.join("out"),
            &mut |m| lines.push(m.to_owned()),
            "01-Jan-26 00:00 UTC",
        )
        .unwrap();
        assert!(c.warnings.is_empty(), "{:?}", c.warnings);
        assert_eq!(c.info["serial"], "1234567890");
        assert_eq!(c.info["firmware"], "4.61");
        assert!(lines[0].starts_with("12341230.T01: R8 Model 2 #1234567890 fw 4.61"));
        assert!(!lines[0].contains("runpkr00"));
        let obs = String::from_utf8(std::fs::read(&c.obs).unwrap()).unwrap();
        assert!(obs.contains("Unpacked from 12341230.T01 by geomaticist-gnss"));
        assert!(obs.contains("TEST"));
        // five epochs of two satellites
        assert_eq!(obs.matches("  2G06R09").count(), 5);
        // the GPS and the GLONASS navigation file
        assert_eq!(c.nav.len(), 2);
        // the dispatcher takes the same path: the file is whole, so runpkr00 is not asked
        let mut lines = Vec::new();
        let c2 = convert(
            &raw,
            &dir.0.join("out2"),
            &mut |m| lines.push(m.to_owned()),
            None,
        )
        .unwrap();
        assert_eq!(lines.len(), 1);
        let obs2 = String::from_utf8(std::fs::read(&c2.obs).unwrap()).unwrap();
        let body = |s: &str| {
            s.lines()
                .filter(|l| !l.contains("PGM / RUN BY"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(body(&obs), body(&obs2));
    }

    #[test]
    fn a_damaged_log_is_written_with_a_warning() {
        let dir = Dir::new("damaged");
        let raw = dir.0.join("12341231.T01");
        let log = native::synth::log(6);
        // damage in the last packet of differences: the last epoch is lost
        let mut stream = log.stream.clone();
        let n = stream.len();
        stream[n - 6] ^= 0x04;
        std::fs::write(&raw, native::synth::file(&stream)).unwrap();
        let c = convert_native_with(&raw, &dir.0.join("out"), &mut |_| {}, "01-Jan-26 00:00 UTC")
            .unwrap();
        assert_eq!(c.warnings.len(), 1);
        assert!(
            c.warnings[0].contains("5 epochs read, 1 lost"),
            "{}",
            c.warnings[0]
        );
        // a file that is no packed log at all
        let other = dir.0.join("notes.T01");
        std::fs::write(&other, vec![0x20; 5000]).unwrap();
        assert!(convert_native(&other, &dir.0.join("out3"), &mut |_| {}).is_err());
        assert!(convert(&other, &dir.0.join("out4"), &mut |_| {}, None).is_err());
    }
}
