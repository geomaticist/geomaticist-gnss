//! Observation codes, signal keys, attribute priorities, signal classes and frequencies.
//!
//! Some rules the results depend on: a blank attribute scores priority 0 in `best_code`; BeiDou
//! band renaming for RINEX < 3.025.

use std::fmt;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::consts::{base_freq, glo_df};
use crate::system::{GnssSystem, Sat};

/// An observation code as written in the file: RINEX 2 two characters (`"L1"`, `"P2"`), RINEX 3/4
/// three (`"L1C"`, `"C2W"`). `Ord` is byte-wise string order.
/// The default is the empty code.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct ObsCode {
    bytes: [u8; 3],
    len: u8,
}

impl ObsCode {
    /// From up to three ASCII bytes (`None` for empty or longer input).
    #[must_use]
    pub fn new(code: &[u8]) -> Option<Self> {
        if code.is_empty() || code.len() > 3 {
            return None;
        }
        let mut bytes = [0u8; 3];
        bytes[..code.len()].copy_from_slice(code);
        Some(Self {
            bytes,
            len: u8::try_from(code.len()).ok()?,
        })
    }

    /// The code's bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[..usize::from(self.len)]
    }

    /// The code as text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).unwrap_or("?")
    }

    /// First character: `C`, `L`, `D`, `S` (or RINEX 2 `P`).
    #[must_use]
    pub fn kind_char(&self) -> char {
        char::from(self.bytes[0])
    }

    /// Number of characters (2 or 3).
    #[must_use]
    #[allow(clippy::len_without_is_empty)] // never empty by construction
    pub fn len(&self) -> usize {
        usize::from(self.len)
    }
}

impl PartialOrd for ObsCode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ObsCode {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_bytes().cmp(other.as_bytes())
    }
}

impl fmt::Debug for ObsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.as_str())
    }
}

impl fmt::Display for ObsCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for ObsCode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ObsCode::new(s.as_bytes()).ok_or_else(|| format!("not an observation code: {s:?}"))
    }
}

impl Serialize for ObsCode {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for ObsCode {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        s.parse().map_err(serde::de::Error::custom)
    }
}

/// What a code measures, by band and tracking attribute, from
/// [`code_key`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct SignalKey {
    /// `'C'` code, `'L'` phase, `'S'` signal strength, `'D'` Doppler (RINEX 2 `'P'` becomes `'C'`).
    pub kind: char,
    /// RINEX 3 band number (0 when the code has no band digit).
    pub band: u8,
    /// Tracking attribute; `None` when blank (RINEX 2 phases, Doppler, SNR).
    pub attr: Option<char>,
}

/// The kinds the engine keeps: code, phase, signal strength (Doppler not used yet).
pub const KINDS: &str = "CLS";

/// `'P2'` / `'L1C'` -> (kind, band, attribute). RINEX 2 `P` codes are kind `'C'` attribute `'P'`;
/// RINEX 2 `C` codes attribute `'C'`. BeiDou band 1 is band 2 (B1I) in RINEX 2 and below 3.025
/// (IGSMAIL-7287).
#[must_use]
pub fn code_key(code: &[u8], sys: char, major: i64, version: f64) -> SignalKey {
    let mut kind = code.first().map_or(' ', |&b| char::from(b));
    let mut attr = if major >= 3 {
        code.get(2).map(|&b| char::from(b))
    } else {
        None
    };
    if kind == 'P' {
        kind = 'C';
        attr = Some('P');
    } else if major < 3 && kind == 'C' {
        attr = Some('C');
    }
    let Some(band) = code
        .get(1)
        .and_then(|&b| char::from(b).to_digit(10))
        .and_then(|d| u8::try_from(d).ok())
    else {
        return SignalKey {
            kind,
            band: 0,
            attr,
        };
    };
    let band = if sys == 'C' && band == 1 && (major < 3 || version < 3.025) {
        2
    } else {
        band
    };
    SignalKey { kind, band, attr }
}

/// Attribute priority within a band: precise / pilot first (RINEX 3), `P` before `C` (RINEX 2).
/// `""` where the table has no entry.
#[must_use]
pub fn attr_priority(sys: char, band: u8) -> &'static str {
    match (sys, band) {
        ('G', 1) => "PWYCMLSXN",
        ('G', 2) => "WPYXLSCDMN",
        ('G' | 'J', 5) | ('R', 3) | ('E', 5 | 7 | 8) => "QXI",
        ('J', 1) => "CSLXZ",
        ('J', 2) => "XLS",
        ('R', 1 | 2) => "PC",
        ('E', 1 | 6) => "CXBAZ",
        ('C', 1) => "PXDSLZ",
        ('C', 2 | 6) => "IXQ",
        ('C', 5 | 8) => "PXD",
        ('C', 7) => "IXQPDZ",
        _ => "",
    }
}

/// The priority score of an attribute in `best_code`: its index in the band's priority string,
/// `len + 1` when absent. A blank attribute scores 0, as if found at the
/// front of every priority string.
#[must_use]
pub fn attr_score(pri: &str, attr: Option<char>) -> usize {
    match attr {
        None => 0,
        Some(a) => pri.find(a).unwrap_or(pri.len() + 1),
    }
}

/// Signal classes: the same tracking at two receivers.
#[must_use]
pub fn class_priority(sys: char, kind: char, band: u8) -> Option<&'static str> {
    Some(match (sys, kind, band) {
        ('G' | 'R', 'C', 1) | ('R', 'C', 2) => "PC",
        ('G', 'C' | 'L', 2) => "PL",
        ('G' | 'J' | 'R', 'L', 1) | ('J', 'C', 1) => "C",
        ('J', 'C' | 'L', 2) => "L",
        ('R', 'L', 2) => "P",
        _ => return None,
    })
}

/// What a receiver tracked, comparable between RINEX 2 and 3: GPS / GLONASS `"P"` (P(Y), W), `"C"`
/// (C/A), `"L"` (L2C); RINEX 2 phases carry no attribute: L1 is C/A, L2 is P(Y). Other systems: one
/// class (`""`).
#[must_use]
pub fn signal_class(sys: char, band: u8, attr: Option<char>) -> String {
    if !"GJR".contains(sys) {
        return String::new();
    }
    match attr {
        Some('P' | 'W' | 'Y') => "P".into(),
        None => if band == 1 { "C" } else { "P" }.into(),
        Some('C') => if "GJ".contains(sys) && band == 2 {
            "L"
        } else {
            "C"
        }
        .into(),
        Some('L' | 'S' | 'X') => "L".into(),
        Some(a) => a.to_string(),
    }
}

/// Error from [`freq`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FreqError {
    /// GLONASS FDMA band without a frequency channel.
    GlonassChannelUnknown(Sat),
    /// No such band for the system.
    NoBand(Sat, u8),
}

impl fmt::Display for FreqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FreqError::GlonassChannelUnknown(s) => {
                write!(f, "{s}: GLONASS frequency channel unknown")
            }
            FreqError::NoBand(s, b) => write!(f, "{s}: no band {b}"),
        }
    }
}

impl std::error::Error for FreqError {}

/// Carrier frequency (Hz) of `band` on satellite `sat`; GLONASS bands 1 and 2 need the channel
/// (`fcn`).
pub fn freq(sat: Sat, band: u8, fcn: Option<i32>) -> Result<f64, FreqError> {
    if sat.sys == GnssSystem::Glonass
        && let Some(df) = glo_df(band)
    {
        let k = fcn.ok_or(FreqError::GlonassChannelUnknown(sat))?;
        let f0 = base_freq(GnssSystem::Glonass, band).ok_or(FreqError::NoBand(sat, band))?;
        return Ok(f0 + f64::from(k) * df);
    }
    base_freq(sat.sys, band).ok_or(FreqError::NoBand(sat, band))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn k(code: &str, sys: char, major: i64, v: f64) -> SignalKey {
        code_key(code.as_bytes(), sys, major, v)
    }

    #[test]
    fn code_keys_as_python() {
        assert_eq!(
            k("P2", 'G', 2, 2.11),
            SignalKey {
                kind: 'C',
                band: 2,
                attr: Some('P')
            }
        );
        assert_eq!(
            k("C1", 'G', 2, 2.11),
            SignalKey {
                kind: 'C',
                band: 1,
                attr: Some('C')
            }
        );
        assert_eq!(
            k("L1", 'G', 2, 2.11),
            SignalKey {
                kind: 'L',
                band: 1,
                attr: None
            }
        );
        assert_eq!(k("L1C", 'G', 3, 3.04).attr, Some('C'));
        assert_eq!(k("C1I", 'C', 3, 3.02).band, 2); // B1I renamed
        assert_eq!(k("C1P", 'C', 3, 3.04).band, 1); // B1C stays
        assert_eq!(k("C1", 'C', 2, 2.12).band, 2);
        assert_eq!(k("CX", 'G', 3, 3.04).band, 0);
    }

    #[test]
    fn priorities_and_classes() {
        assert_eq!(attr_score("PWYCMLSXN", Some('C')), 3);
        assert_eq!(attr_score("PWYCMLSXN", None), 0);
        assert_eq!(attr_score("QXI", Some('Z')), 4);
        assert_eq!(attr_score("", Some('X')), 1);
        assert_eq!(signal_class('G', 2, Some('C')), "L");
        assert_eq!(signal_class('R', 2, Some('C')), "C");
        assert_eq!(signal_class('G', 2, None), "P");
        assert_eq!(signal_class('G', 1, Some('W')), "P");
        assert_eq!(signal_class('E', 1, Some('C')), "");
        assert_eq!(signal_class('G', 5, Some('Q')), "Q");
        assert_eq!(class_priority('G', 'L', 2), Some("PL"));
        assert_eq!(class_priority('E', 'L', 1), None);
    }

    #[test]
    fn frequencies() {
        let r: Sat = "R07".parse().unwrap();
        assert_eq!(freq(r, 1, Some(-2)).unwrap(), 1602.0e6 - 2.0 * 0.5625e6);
        assert!(freq(r, 1, None).is_err());
        assert_eq!(freq(r, 3, None).unwrap(), 1202.025e6);
        let g: Sat = "G01".parse().unwrap();
        assert!(freq(g, 7, None).is_err());
    }
}
