//! Trimble .T04 files, decoded with no Trimble software.
//!
//! A clean-room decoder written from sample receiver logs. It writes a full-rate RINEX 3.05 observation file
//! (code, phase, Doppler, SNR) and RINEX 3.05 navigation files.

pub mod avro;
pub mod container;
pub mod header;
pub mod measurements;
pub mod nav;
pub mod nav_writer;
pub mod obs_writer;
pub mod rec35;
pub mod signals;

use std::fmt;
use std::path::Path;

use serde_json::json;

use crate::{Converted, DecodeError};

/// PGM of the files this decoder writes.
pub const PROGRAM: &str = "geomaticist-gnss T04";
/// RUN BY.
pub const RUN_BY: &str = "geomaticist-gnss";

/// A RINEX satellite id as the decoder spells it (`f"{sys}{prn:02d}"`), ordered as the string.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SatId([u8; 8]);

impl SatId {
    /// `sys` + PRN with at least two digits.
    #[must_use]
    pub fn new(sys: char, prn: i64) -> Self {
        let s = format!("{sys}{prn:02}");
        let mut b = [0u8; 8];
        for (d, c) in b.iter_mut().zip(s.bytes()) {
            *d = c;
        }
        Self(b)
    }

    /// The system letter.
    #[must_use]
    pub fn sys(&self) -> char {
        char::from(self.0[0])
    }

    /// The text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        let n = self.0.iter().position(|&c| c == 0).unwrap_or(8);
        std::str::from_utf8(&self.0[..n]).unwrap_or("?")
    }
}

impl fmt::Display for SatId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl fmt::Debug for SatId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Errors of the T04 decoder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum T04Error {
    /// The data could not be decoded.
    Decode(String),
    /// A file could not be read or written.
    Io(String, String),
}

impl fmt::Display for T04Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            T04Error::Decode(s) => f.write_str(s),
            T04Error::Io(p, e) => write!(f, "{p}: {e}"),
        }
    }
}

impl std::error::Error for T04Error {}

impl From<container::ContainerError> for T04Error {
    fn from(e: container::ContainerError) -> Self {
        T04Error::Decode(e.0)
    }
}

impl From<rec35::Rec35Error> for T04Error {
    fn from(e: rec35::Rec35Error) -> Self {
        T04Error::Decode(e.0)
    }
}

/// True for a `.t04` name.
#[must_use]
pub fn is_t04(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("t04"))
}

fn stem_of(path: &Path) -> String {
    path.file_stem()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned())
}

fn name_of(path: &Path) -> String {
    path.file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned())
}

/// Decodes one .T04 to a RINEX 3.05 observation file plus navigation files in `outdir`.
/// Names follow the usual RINEX short names: `12341541.T04` -> `12341541.26o`, `.26n`,
/// `.26g`, `.26l`, `.26c`, `.26q`, `.26h`.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    convert_with(raw, outdir, log, &crate::now_utc_text())
}

/// [`convert`] with the `PGM / RUN BY / DATE` date text given.
pub fn convert_with(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
) -> Result<Converted, DecodeError> {
    let name = name_of(raw);
    std::fs::create_dir_all(outdir).map_err(|e| DecodeError::io(outdir, &e))?;
    let data = std::fs::read(raw)
        .map_err(|e| DecodeError::Format(format!("{name} isn't a readable T04 file ({e})")))?;
    let inner = container::inner_stream(&data);
    drop(data);
    for d in &inner.damage {
        log(&format!(
            "Damaged data in {name} at byte {} ({}): skipped {} bytes, as Trimble's converter does",
            group3(d.file_offset),
            d.reason,
            group3(d.bytes_skipped)
        ));
    }
    let hi =
        header::read_header(&inner).map_err(|e| DecodeError::Format(format!("{name}: {e}")))?;
    let info = json!({
        "receiver": hi.model, "serial": hi.serial, "firmware": hi.firmware, "antenna": hi.antenna_type,
        "radome": hi.antenna_dome, "antenna_code": hi.antenna_code, "arp_height": hi.antenna_height_arp,
        "marker": hi.marker_name,
    });
    let mut warnings = Vec::new();
    if !inner.damage.is_empty() {
        let n = inner.damage.len();
        warnings.push(format!(
            "{name} has {n} damaged block{} (written that way by the receiver). Decoding skipped over it, as \
             Trimble's converter does: expect a gap of a few minutes in the data there.",
            if n > 1 { "s" } else { "" }
        ));
    }
    if hi.antenna_type.starts_with("UNKNOWN") {
        warnings.push(format!(
            "Antenna code {} isn't in the antenna table, so the header says {}. Set the antenna type on the \
             Header page (e.g. TRMR10-2).",
            hi.antenna_code.map_or_else(|| "None".to_owned(), |c| c.to_string()),
            hi.antenna_type
        ));
    }
    log(&format!(
        "Built-in T04 decoder: {} #{} fw {}, antenna {} {}, ARP height {:.4} m, marker '{}'",
        hi.model,
        hi.serial,
        hi.firmware,
        hi.antenna_type,
        hi.antenna_dome,
        hi.antenna_height_arp,
        hi.marker_name
    ));
    let opt = obs_writer::ObsOptions {
        program: PROGRAM.into(),
        run_by: RUN_BY.into(),
        date: date.to_owned(),
        ..obs_writer::ObsOptions::default()
    };
    let mut text = Vec::new();
    let (_head, first) = obs_writer::write_obs(&inner, &hi, &name, &opt, &mut text).map_err(|e| {
        DecodeError::Format(format!(
            "the built-in decoder couldn't read {name} ({e}). If Trimble's converter is installed, turn on \
             \"Use Trimble's converter for .T04\" in Settings, and please send the file so the decoder can learn it."
        ))
    })?;
    let stem = stem_of(raw);
    let obs = outdir.join(format!("{stem}.{:02}o", first.year.rem_euclid(100)));
    crate::write_atomic(&obs, &text)?;
    drop(text);
    let nav = nav_writer::write_nav(&inner, &name, outdir, &stem, PROGRAM, RUN_BY, date)
        .map_err(|e| DecodeError::Format(format!("{name}: {e}")))?;
    let navs: Vec<_> = nav.into_iter().map(|(_, p)| p).collect();
    log(&format!(
        "Built-in T04 decoder wrote {}{}",
        name_of(&obs),
        if navs.is_empty() {
            String::new()
        } else {
            format!(" and {} navigation file(s)", navs.len())
        }
    ));
    Ok(Converted {
        obs,
        nav: navs,
        warnings,
        info,
    })
}

/// `f"{n:,}"`.
#[must_use]
pub fn group3(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sat_ids_sort_as_strings() {
        let mut v = [
            SatId::new('R', 3),
            SatId::new('C', 21),
            SatId::new('G', 5),
            SatId::new('G', 100),
        ];
        v.sort();
        let s: Vec<String> = v.iter().map(ToString::to_string).collect();
        assert_eq!(s, vec!["C21", "G05", "G100", "R03"]);
        assert_eq!(group3(1_234_567), "1,234,567");
        assert_eq!(group3(12), "12");
    }
}
