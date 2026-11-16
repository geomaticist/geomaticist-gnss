//! `geomaticist-gnss-decoders`: receiver raw-format decoders (Trimble T04 and packed T01, Leica MDB/DBX,
//! Topcon GREIS) to RINEX.
//!
//! Modules:
//! - Trimble .T04 -> [`t04`]
//! - Trimble .T01 / .T00 / .T02 / .R00 (the unpacked records to RINEX) -> [`t01`]. A built-in reader unpacks
//!   the log ([`t01::native`]); Trimble's runpkr00 (located on the machine, never bundled) is only the
//!   fallback for a file that reader doesn't know
//! - Leica .m00 -> [`mdb`]
//! - raw-log header readers -> [`makers`]; Leica Captivate job tables -> [`dbx`]
//! - Topcon / Javad GREIS logs (.tps) -> [`topcon`]: written here from our own logs, with the
//!   public GREIS framing
//! - Standard GREIS messages (.jps) -> [`topcon::greis`], RTCM 3 -> [`rtcm3`], u-blox UBX -> [`ubx`], NovAtel
//!   OEM4/OEMV/OEM6 -> [`novatel`], and the recognise-only list of other receiver logs -> [`sniff`]: written here
//!   from the makers' published message references. Where a reference leaves a rule open (scale offsets,
//!   loss-of-lock tests, signal tables), the choice follows RTKLIB 2.4.3 (T. Takasu, BSD-2-Clause), whose
//!   logic is restated in our own code with credit in each module; no RTKLIB code is copied
//! - Septentrio SBF -> [`sbf`]
//! - antenna codes, measure-up methods and the antennas table as data -> [`antenna_lib`], [`antennas`]
//!
//! The decoders live in the engine repository as this crate, but rule R2 holds: no
//! processing crate depends on `geomaticist-gnss-decoders`, so the RINEX-only path stays clean.
//!
//! Clean-room: decoders written from sample receiver logs.
//!
//! Every writer's RINEX text is deterministic, byte for byte. The `PGM / RUN BY / DATE` date is the
//! only wall-clock field; the `*_with` entry points take it as a parameter so tests can fix it.

// Numeric code: casts between integer fields, counts and f64 are deliberate and bounded (the record layouts
// fix the widths); exact float comparisons are intended; the field-by-field decoding keeps short names and
// long functions so each can be read against its record layout.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::cast_lossless,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::struct_excessive_bools,
    // docs cite record fields verbatim; constants are transcribed as published
    clippy::doc_markdown,
    clippy::unreadable_literal,
    clippy::type_complexity,
    // the writers build text line by line; clarity over the last allocation
    clippy::format_push_string,
    clippy::format_collect,
    clippy::assigning_clones,
    // signal and code tables are listed entry by entry; small helpers keep an Option shape
    clippy::match_same_arms,
    clippy::unnecessary_wraps
)]

pub mod antenna_lib;
pub mod antennas;
pub mod dbx;
pub mod decimal;
pub(crate) mod epochs;
pub mod makers;
pub mod mdb;
pub mod novatel;
pub mod rtcm3;
pub mod sbf;
pub mod sniff;
pub mod t01;
pub mod t04;
pub mod topcon;
pub mod ubx;
pub mod util;

use std::path::{Path, PathBuf};

use geomaticist_gnss_core::time::DateTime;

/// What a decoder wrote.
#[derive(Debug, Clone)]
pub struct Converted {
    /// The observation file.
    pub obs: PathBuf,
    /// The navigation files.
    pub nav: Vec<PathBuf>,
    /// Warnings for the user.
    pub warnings: Vec<String>,
    /// Header facts (receiver, antenna, height ...).
    pub info: serde_json::Value,
}

/// Decoder errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecodeError {
    /// The file isn't one the decoder can read, or decoding failed.
    Format(String),
    /// A file could not be read or written.
    Io(String),
    /// Trimble's runpkr00 was needed (a packed log the built-in reader doesn't know) and is not available, or
    /// failed.
    Runpkr(String),
    /// Not a raw log any decoder here handles.
    Unsupported(String),
}

impl DecodeError {
    fn io(path: &Path, e: &std::io::Error) -> Self {
        DecodeError::Io(format!("{}: {e}", path.display()))
    }
}

impl std::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecodeError::Format(s)
            | DecodeError::Io(s)
            | DecodeError::Runpkr(s)
            | DecodeError::Unsupported(s) => f.write_str(s),
        }
    }
}

impl std::error::Error for DecodeError {}

/// The current UTC time.
#[must_use]
pub fn now_utc() -> DateTime {
    let unix_us = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| i64::try_from(d.as_micros()).unwrap_or(0));
    // 1970-01-01 is 315 964 800 s before 1980-01-06
    DateTime::from_micros_since_gps0(unix_us - 315_964_800_000_000)
}

/// `f"{now:%Y%m%d %H%M%S} UTC"` (the RINEX 3 `PGM / RUN BY / DATE` date).
#[must_use]
pub fn now_utc_text() -> String {
    let t = now_utc();
    format!(
        "{:04}{:02}{:02} {:02}{:02}{:02} UTC",
        t.year, t.month, t.day, t.hour, t.minute, t.second
    )
}

/// Writes `bytes` to `path` through a `.part` file renamed into place.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), DecodeError> {
    let mut part = path.as_os_str().to_owned();
    part.push(".part");
    let part = PathBuf::from(part);
    std::fs::write(&part, bytes).map_err(|e| DecodeError::io(&part, &e))?;
    if path.exists() {
        std::fs::remove_file(path).map_err(|e| DecodeError::io(path, &e))?;
    }
    std::fs::rename(&part, path).map_err(|e| {
        let _ = std::fs::remove_file(&part);
        DecodeError::io(path, &e)
    })
}

/// The kind of raw log a file name says it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RawKind {
    /// Trimble .T04 (built-in decoder).
    T04,
    /// Trimble packed log (.T01 / .T00 / .T02 / .R00): the built-in reader (runpkr00 for a file it doesn't know),
    /// then our RINEX 2.11 writer.
    Packed,
    /// Leica receiver log (.m00 ...).
    Mdb,
    /// Topcon / Javad GREIS log (.tps) with `rD` measurements.
    Tps,
    /// JAVAD / Topcon GREIS log with the standard measurement messages (.jps, or a .tps without `rD`).
    Greis,
    /// Logged RTCM 3 stream (.rtcm3 ...).
    Rtcm3,
    /// u-blox UBX log (RXM-RAW / RXM-RAWX).
    Ubx,
    /// NovAtel OEM4 / OEMV / OEM6 binary log (RANGEB / RANGECMPB).
    NovAtel,
    /// Septentrio SBF log (MeasEpoch, navigation and setup blocks).
    Sbf,
}

/// The raw-log kind of a path, by its extension (`None`: not a raw log decoded here).
#[must_use]
pub fn raw_kind(path: &Path) -> Option<RawKind> {
    let ext = path.extension()?.to_string_lossy().to_ascii_lowercase();
    match ext.as_str() {
        "t04" => Some(RawKind::T04),
        e if t01::runpkr::RAW_EXTS.contains(&e) => Some(RawKind::Packed),
        e if mdb::is_mdb_ext(e) => Some(RawKind::Mdb),
        e if topcon::EXTS.contains(&e) => Some(RawKind::Tps),
        "jps" => Some(RawKind::Greis),
        e if rtcm3::is_rtcm_ext(e) => Some(RawKind::Rtcm3),
        "ubx" => Some(RawKind::Ubx),
        "sbf" => Some(RawKind::Sbf),
        _ => None,
    }
}

/// Extensions receiver logs often carry whatever made them (`.gps`, `.log`, `.bin`, ...): the content decides.
pub const SNIFF_EXTS: [&str; 6] = ["gps", "log", "bin", "raw", "sbas", ""];

/// What a file is, by its content (the first 256 KiB), with its extension deciding first for the formats that
/// have their own (`.T04`, `.m00`, ...). A `.tps` without `rD` measurements is read as standard GREIS; an RTCM 2
/// stream named `.rtcm2` is named as such.
#[must_use]
pub fn identify(path: &Path) -> sniff::Sniffed {
    use std::io::Read as _;
    let ext = path
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    if ext == "rtcm2" || ext == "rtc2" {
        return sniff::Sniffed::Unsupported("RTCM 2");
    }
    let by_ext = raw_kind(path);
    if matches!(by_ext, Some(RawKind::T04 | RawKind::Packed | RawKind::Mdb)) {
        return sniff::Sniffed::Supported(by_ext.unwrap_or(RawKind::T04));
    }
    let mut head = Vec::new();
    if let Ok(f) = std::fs::File::open(path) {
        let _ = f.take(256 * 1024).read_to_end(&mut head);
    }
    if head.is_empty() {
        return sniff::Sniffed::Unknown;
    }
    let seen = sniff::sniff(&head);
    match (by_ext, seen) {
        // the content is the better judge of which GREIS reader, and of a mislabelled log
        (_, sniff::Sniffed::Supported(_) | sniff::Sniffed::Unsupported(_)) => seen,
        (Some(k), sniff::Sniffed::Unknown) => sniff::Sniffed::Supported(k),
        (None, sniff::Sniffed::Unknown) => sniff::Sniffed::Unknown,
    }
}

/// Decodes any supported raw log into `outdir`: `.T04` to RINEX 3.05,
/// packed Trimble logs to RINEX 2.11 (the built-in reader; runpkr00, when found, for a file it doesn't know),
/// Leica `.m00` to RINEX 3.04, Topcon `.tps` to RINEX 3.05.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    let near = near_from_name(raw);
    match identify(raw) {
        sniff::Sniffed::Supported(RawKind::T04) => t04::convert(raw, outdir, log),
        sniff::Sniffed::Supported(RawKind::Packed) => t01::convert(raw, outdir, log, None),
        sniff::Sniffed::Supported(RawKind::Mdb) => mdb::convert(raw, outdir, log),
        sniff::Sniffed::Supported(RawKind::Tps) => topcon::convert(raw, outdir, log),
        sniff::Sniffed::Supported(RawKind::Greis) => topcon::greis::convert(raw, outdir, near, log),
        sniff::Sniffed::Supported(RawKind::Rtcm3) => {
            rtcm3::convert(raw, outdir, log, near, &std::collections::BTreeMap::new())
        }
        sniff::Sniffed::Supported(RawKind::Ubx) => ubx::convert(raw, outdir, log),
        sniff::Sniffed::Supported(RawKind::NovAtel) => novatel::convert(raw, outdir, log),
        sniff::Sniffed::Supported(RawKind::Sbf) => sbf::convert(raw, outdir, log),
        sniff::Sniffed::Unsupported(what) => Err(DecodeError::Unsupported(format!(
            "{}: recognised as {what}: not supported",
            raw.display()
        ))),
        sniff::Sniffed::Unknown => Err(DecodeError::Unsupported(format!(
            "{}: not a raw log this decoder reads (Trimble .T04/.T01/.T00/.T02/.R00, Leica .m00, Topcon / \
             JAVAD .tps/.jps, RTCM 3, u-blox UBX, NovAtel OEM4/OEMV/OEM6, \n             Septentrio SBF)",
            raw.display()
        ))),
    }
}

/// An approximate GPS time (s) of a log from a `YYYYMMDD` date in its file name (`ubx_20080526.ubx`), else
/// `None` (the decoders then take the file's modification time).
#[must_use]
pub fn near_from_name(raw: &Path) -> Option<f64> {
    let stem = raw.file_stem()?.to_string_lossy().into_owned();
    let b = stem.as_bytes();
    let mut i = 0;
    while i + 8 <= b.len() {
        let run = b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        if run == 8 {
            let num = |r: std::ops::Range<usize>| stem[r].parse::<i64>().ok();
            let (y, m, d) = (num(i..i + 4)?, num(i + 4..i + 6)?, num(i + 6..i + 8)?);
            if (1990..=2100).contains(&y)
                && let Some(dt) = geomaticist_gnss_core::time::DateTime::new(y, m, d, 12, 0, 0)
            {
                return Some(geomaticist_gnss_core::time::utc_to_gps_seconds(&dt));
            }
        }
        i += run.max(1);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinds() {
        assert_eq!(raw_kind(Path::new("a/12341541.T04")), Some(RawKind::T04));
        assert_eq!(raw_kind(Path::new("x.t01")), Some(RawKind::Packed));
        assert_eq!(
            raw_kind(Path::new("0808_1002_160139.m00")),
            Some(RawKind::Mdb)
        );
        assert_eq!(raw_kind(Path::new("x.TPS")), Some(RawKind::Tps));
        assert_eq!(raw_kind(Path::new("x.rnx")), None);
        assert_eq!(raw_kind(Path::new("x.jps")), Some(RawKind::Greis));
        assert_eq!(raw_kind(Path::new("x.ubx")), Some(RawKind::Ubx));
        assert_eq!(raw_kind(Path::new("x.SBF")), Some(RawKind::Sbf));
        let t = near_from_name(Path::new("ubx_20080526.ubx")).unwrap();
        let dt = geomaticist_gnss_core::time::DateTime::new(2008, 5, 26, 12, 0, 14).unwrap();
        assert_eq!(t, geomaticist_gnss_core::time::gps_seconds(&dt));
        assert_eq!(near_from_name(Path::new("oemv_200911218.gps")), None);
        assert_eq!(now_utc_text().len(), 19);
    }
}
