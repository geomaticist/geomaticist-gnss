//! Receiver makers behind one interface: which raw logs each maker owns, whether they can be decoded here,
//! and what to tell the user when they can't.
//!
//! The dispatch: `maker_of`, `is_raw`, `can_decode`, `undecodable_note`, `to_rinex`. The Leica DBX table
//! reader is [`crate::dbx`]; job importers (Captivate, Trimble JobXML) need a job / CRS model and are not
//! part of this crate.

use std::path::Path;

use crate::{Converted, DecodeError};

/// A receiver maker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Maker {
    /// Trimble (and Spectra Geospatial).
    Trimble,
    /// Leica Geosystems.
    Leica,
    /// Topcon (and Javad / Sokkia GREIS receivers).
    Topcon,
    /// NovAtel OEM4 / OEMV / OEM6 binary logs.
    NovAtel,
    /// u-blox UBX logs.
    Ublox,
    /// Septentrio SBF logs.
    Septentrio,
    /// Logged RTCM 3 streams (any receiver or caster).
    Rtcm,
}

/// Every maker, in a fixed order.
pub const MAKERS: [Maker; 7] = [
    Maker::Trimble,
    Maker::Leica,
    Maker::Topcon,
    Maker::NovAtel,
    Maker::Ublox,
    Maker::Septentrio,
    Maker::Rtcm,
];

fn ext_lower(path: &Path) -> String {
    path.extension()
        .map_or_else(String::new, |e| e.to_string_lossy().to_lowercase())
}

fn ext_upper(path: &Path) -> String {
    path.extension()
        .map_or_else(String::new, |e| e.to_string_lossy().to_uppercase())
}

fn file_name(path: &Path) -> String {
    path.file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned())
}

impl Maker {
    /// The maker as the user knows it.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Maker::Trimble => "Trimble",
            Maker::Leica => "Leica",
            Maker::Topcon => "Topcon",
            Maker::NovAtel => "NovAtel",
            Maker::Ublox => "u-blox",
            Maker::Septentrio => "Septentrio",
            Maker::Rtcm => "RTCM",
        }
    }

    /// The raw-log suffixes it owns, lower case, without the dot.
    #[must_use]
    pub fn raw_exts(self) -> &'static [&'static str] {
        match self {
            Maker::Trimble => &["t04", "t02", "t01", "t00", "dat", "t0c"],
            Maker::Leica => &["m00"],
            Maker::Topcon => &["tps", "jps"],
            Maker::NovAtel => &[],
            Maker::Ublox => &["ubx"],
            Maker::Septentrio => &["sbf"],
            Maker::Rtcm => &["rtcm3", "rtcm", "rtc", "rt3"],
        }
    }

    /// The file is one of its raw logs (by suffix).
    #[must_use]
    pub fn is_raw(self, path: &Path) -> bool {
        self.raw_exts().contains(&ext_lower(path).as_str())
    }

    /// Its raw log can be turned into RINEX here: .T04 always, packed Trimble logs the built-in reader knows
    /// (any of them when runpkr00 is found), a Leica .m00 that holds observation records, a Topcon .tps with
    /// `rD` measurements.
    #[must_use]
    pub fn can_decode(self, path: &Path) -> bool {
        let e = ext_lower(path);
        match self {
            Maker::Trimble => {
                e == "t04"
                    || (crate::t01::runpkr::RAW_EXTS.contains(&e.as_str())
                        && crate::t01::can_convert(path))
            }
            Maker::Leica => {
                is_mdb(path) && std::fs::read(path).is_ok_and(|d| crate::mdb::has_epochs(&d))
            }
            Maker::Topcon | Maker::NovAtel | Maker::Ublox | Maker::Septentrio | Maker::Rtcm => {
                matches!(crate::identify(path), crate::sniff::Sniffed::Supported(_))
            }
        }
    }

    /// What the user should do with a raw log it can't decode.
    #[must_use]
    pub fn undecodable_note(self, path: &Path) -> String {
        let (name, up) = (file_name(path), ext_upper(path));
        match self {
            Maker::Trimble if crate::t01::runpkr::RAW_EXTS.contains(&ext_lower(path).as_str()) => {
                format!(
                    "{name}: this {up} isn't one the built-in reader knows, and Trimble's runpkr00, which could \
                 unpack it, wasn't found (put it beside this program or on PATH, or set RUNPKR00 to its path). \
                 Otherwise convert it to RINEX with Trimble's Convert to RINEX or Business Center and import \
                 the RINEX file."
                )
            }
            Maker::Trimble => format!(
                "{name}: only .T04 raw files can be decoded. Convert {up} files to RINEX (e.g. with Trimble's \
                 Convert to RINEX or Business Center) and import the RINEX file."
            ),
            Maker::Leica => format!(
                "{name} holds no GNSS observations (a Captivate job's own .m00 has only its points). Import the \
                 base receiver's .m00 log, or RINEX converted with Leica Infinity / Leica Geo Office."
            ),
            Maker::Topcon => format!(
                "{name}: this GREIS log has no measurement messages the built-in reader knows (the compact \
                 rD / rP records, or the standard SI with pseudorange and phase messages). Convert it to RINEX \
                 with Topcon's MAGNET Tools or Javad's software and import the RINEX file."
            ),
            Maker::NovAtel | Maker::Ublox | Maker::Septentrio | Maker::Rtcm => {
                match crate::identify(path) {
                    crate::sniff::Sniffed::Unsupported(what) => format!(
                        "{name}: recognised as {what}: not supported. Convert it to RINEX with the receiver \
                     maker's software and import the RINEX file."
                    ),
                    _ => format!(
                        "{name}: this {} log has no raw measurement messages the built-in reader knows. Log the \
                     raw observations (RANGEB, RXM-RAWX, MeasEpoch, MSM) or convert it to RINEX and import that.",
                        self.name()
                    ),
                }
            }
        }
    }

    /// Decodes its raw log into `folder`.
    pub fn to_rinex(
        self,
        path: &Path,
        folder: &Path,
        log: &mut dyn FnMut(&str),
    ) -> Result<Converted, DecodeError> {
        match self {
            Maker::Trimble if ext_lower(path) == "t04" => crate::t04::convert(path, folder, log),
            Maker::Trimble => crate::t01::convert(path, folder, log, None),
            Maker::Leica => crate::mdb::convert(path, folder, log),
            Maker::Topcon | Maker::NovAtel | Maker::Ublox | Maker::Septentrio | Maker::Rtcm => {
                crate::convert(path, folder, log)
            }
        }
    }
}

/// A Leica MDB log: the framing sync then the header record.
#[must_use]
pub fn is_mdb(path: &Path) -> bool {
    use std::io::Read as _;
    let mut head = [0u8; 8];
    let Ok(mut f) = std::fs::File::open(path) else {
        return false;
    };
    let mut n = 0;
    while n < head.len() {
        match f.read(&mut head[n..]) {
            Ok(0) | Err(_) => break,
            Ok(k) => n += k,
        }
    }
    n > 4 && head[..2] == [0x9c, 0xae] && head[4] == crate::mdb::T_HEADER
}

/// The maker whose raw log this is, or `None` (RINEX, or not a raw log). A file whose extension says nothing
/// (`.gps`, `.log`, `.bin`, ...) is judged by its content.
#[must_use]
pub fn maker_of(path: &Path) -> Option<Maker> {
    if let Some(m) = MAKERS.into_iter().find(|m| m.is_raw(path)) {
        return Some(m);
    }
    if !crate::SNIFF_EXTS.contains(&ext_lower(path).as_str()) {
        return None;
    }
    match crate::identify(path) {
        crate::sniff::Sniffed::Supported(k) => Some(match k {
            crate::RawKind::T04 | crate::RawKind::Packed => Maker::Trimble,
            crate::RawKind::Mdb => Maker::Leica,
            crate::RawKind::Tps | crate::RawKind::Greis => Maker::Topcon,
            crate::RawKind::NovAtel => Maker::NovAtel,
            crate::RawKind::Ubx => Maker::Ublox,
            crate::RawKind::Sbf => Maker::Septentrio,
            crate::RawKind::Rtcm3 => Maker::Rtcm,
        }),
        _ => None,
    }
}

/// True for any maker's raw log.
#[must_use]
pub fn is_raw(path: &Path) -> bool {
    maker_of(path).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatch() {
        assert_eq!(maker_of(Path::new("a.T04")), Some(Maker::Trimble));
        assert_eq!(maker_of(Path::new("a.m00")), Some(Maker::Leica));
        assert_eq!(maker_of(Path::new("a.TPS")), Some(Maker::Topcon));
        assert_eq!(maker_of(Path::new("a.25o")), None);
        assert_eq!(maker_of(Path::new("a.jps")), Some(Maker::Topcon));
        assert_eq!(maker_of(Path::new("a.ubx")), Some(Maker::Ublox));
        assert_eq!(maker_of(Path::new("a.SBF")), Some(Maker::Septentrio));
        assert_eq!(maker_of(Path::new("a.rtcm3")), Some(Maker::Rtcm));
        // a .gps that isn't there is nobody's
        assert_eq!(maker_of(Path::new("no such file.gps")), None);
        assert!(!Maker::Topcon.can_decode(Path::new("no such file.tps")));
        assert!(Maker::Trimble.can_decode(Path::new("x.t04")));
        assert!(
            Maker::Trimble
                .undecodable_note(Path::new("x.dat"))
                .contains("only .T04")
        );
    }
}
