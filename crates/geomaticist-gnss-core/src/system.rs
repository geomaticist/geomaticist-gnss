//! Satellite systems and satellite ids.
//!
//! A satellite is a typed [`Sat`] whose `Ord` is the order of its id strings (`"G05"`, `"R12"`),
//! so every map keyed by `Sat` iterates in sorted id order.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// A GNSS constellation, by its RINEX letter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum GnssSystem {
    /// `G` GPS.
    Gps,
    /// `R` GLONASS.
    Glonass,
    /// `E` Galileo.
    Galileo,
    /// `C` BeiDou.
    BeiDou,
    /// `J` QZSS.
    Qzss,
    /// `I` NavIC / IRNSS.
    Navic,
    /// `S` SBAS.
    Sbas,
}

/// System order in RINEX headers and listings.
pub const SYSTEM_ORDER: &str = "GRECJIS";

impl GnssSystem {
    /// Every system in [`SYSTEM_ORDER`].
    pub const ALL: [GnssSystem; 7] = [
        GnssSystem::Gps,
        GnssSystem::Glonass,
        GnssSystem::Galileo,
        GnssSystem::BeiDou,
        GnssSystem::Qzss,
        GnssSystem::Navic,
        GnssSystem::Sbas,
    ];

    /// The RINEX letter.
    #[must_use]
    pub const fn letter(self) -> char {
        match self {
            GnssSystem::Gps => 'G',
            GnssSystem::Glonass => 'R',
            GnssSystem::Galileo => 'E',
            GnssSystem::BeiDou => 'C',
            GnssSystem::Qzss => 'J',
            GnssSystem::Navic => 'I',
            GnssSystem::Sbas => 'S',
        }
    }

    /// From a RINEX letter (upper case only).
    #[must_use]
    pub const fn from_letter(c: char) -> Option<Self> {
        Some(match c {
            'G' => GnssSystem::Gps,
            'R' => GnssSystem::Glonass,
            'E' => GnssSystem::Galileo,
            'C' => GnssSystem::BeiDou,
            'J' => GnssSystem::Qzss,
            'I' => GnssSystem::Navic,
            'S' => GnssSystem::Sbas,
            _ => return None,
        })
    }

    /// Display name (`header.SYSTEM_NAMES`).
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            GnssSystem::Gps => "GPS",
            GnssSystem::Glonass => "GLONASS",
            GnssSystem::Galileo => "Galileo",
            GnssSystem::BeiDou => "BeiDou",
            GnssSystem::Qzss => "QZSS",
            GnssSystem::Navic => "NavIC",
            GnssSystem::Sbas => "SBAS",
        }
    }

    /// Index in [`SYSTEM_ORDER`] (`sat_sort_key`'s first element).
    #[must_use]
    pub const fn order(self) -> usize {
        match self {
            GnssSystem::Gps => 0,
            GnssSystem::Glonass => 1,
            GnssSystem::Galileo => 2,
            GnssSystem::BeiDou => 3,
            GnssSystem::Qzss => 4,
            GnssSystem::Navic => 5,
            GnssSystem::Sbas => 6,
        }
    }
}

impl fmt::Display for GnssSystem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.letter())
    }
}

/// A satellite: system and PRN / slot (0..=99, RINEX's two digits).
///
/// `Ord` is the order of the id strings (`"C01" < "E02" < "G01" < "R03"`): by letter, then PRN.
/// Use [`Sat::rinex_order`] for RINEX order (system order `GRECJIS`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Sat {
    /// Constellation.
    pub sys: GnssSystem,
    /// PRN or slot number.
    pub prn: u8,
}

impl Sat {
    /// A satellite (PRN above 99 is clamped by RINEX's two digits; callers pass 0..=99).
    #[must_use]
    pub const fn new(sys: GnssSystem, prn: u8) -> Self {
        Self { sys, prn }
    }

    /// RINEX sort key: (system order, PRN).
    #[must_use]
    pub const fn rinex_order(self) -> (usize, u8) {
        (self.sys.order(), self.prn)
    }

    /// The three-character id as bytes (`b"G05"`).
    #[must_use]
    pub fn id_bytes(self) -> [u8; 3] {
        [
            self.sys.letter() as u8,
            b'0' + self.prn / 10 % 10,
            b'0' + self.prn % 10,
        ]
    }
}

impl PartialOrd for Sat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Sat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.sys.letter(), self.prn).cmp(&(other.sys.letter(), other.prn))
    }
}

impl fmt::Display for Sat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{:02}", self.sys.letter(), self.prn)
    }
}

/// Error parsing a satellite id.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseSatError(pub String);

impl fmt::Display for ParseSatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "not a satellite id: {:?}", self.0)
    }
}

impl std::error::Error for ParseSatError {}

impl FromStr for Sat {
    type Err = ParseSatError;

    /// `"G05"`, `"R12"`, `"E1"`: a system letter and a number 0..=99.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let err = || ParseSatError(s.to_owned());
        let mut ch = s.chars();
        let sys = ch
            .next()
            .and_then(GnssSystem::from_letter)
            .ok_or_else(err)?;
        let prn: u8 = ch.as_str().trim().parse().map_err(|_| err())?;
        if prn > 99 {
            return Err(err());
        }
        Ok(Sat::new(sys, prn))
    }
}

impl Serialize for Sat {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_str(self)
    }
}

impl<'de> Deserialize<'de> for Sat {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        s.parse().map_err(serde::de::Error::custom)
    }
}

impl Serialize for GnssSystem {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_str(self)
    }
}

impl<'de> Deserialize<'de> for GnssSystem {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let s = String::deserialize(d)?;
        let mut c = s.chars();
        match (c.next().and_then(GnssSystem::from_letter), c.next()) {
            (Some(g), None) => Ok(g),
            _ => Err(serde::de::Error::custom(format!(
                "not a system letter: {s:?}"
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_is_python_string_order() {
        let mut v: Vec<Sat> = ["R03", "G10", "E02", "C01", "G02", "J01", "S20", "I05"]
            .iter()
            .map(|s| s.parse().unwrap())
            .collect();
        v.sort();
        let s: Vec<String> = v.iter().map(ToString::to_string).collect();
        let mut py = s.clone();
        py.sort();
        assert_eq!(s, py);
        assert_eq!(s[0], "C01");
        let g: Sat = "G5".parse().unwrap();
        assert_eq!(g.to_string(), "G05");
        assert_eq!(&g.id_bytes(), b"G05");
        assert!("X01".parse::<Sat>().is_err());
        assert!(
            Sat::new(GnssSystem::Gps, 1).rinex_order()
                < Sat::new(GnssSystem::Glonass, 1).rinex_order()
        );
    }

    #[test]
    fn serde_as_string() {
        let s: Sat = "E24".parse().unwrap();
        let j = serde_json::to_string(&s).unwrap();
        assert_eq!(j, "\"E24\"");
        assert_eq!(serde_json::from_str::<Sat>(&j).unwrap(), s);
    }
}
