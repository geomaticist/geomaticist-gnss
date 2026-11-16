//! Broadcast navigation data as the engine consumes it: Keplerian records (GPS, Galileo, BeiDou,
//! QZSS) and GLONASS state-vector records, per satellite in file order.
//!
//! Record layouts: Keplerian fields in [`KEP_FIELDS`] order; GLONASS `toe, taun, gamn, pos, vel,
//! acc, svh, fcn`. Each record is a struct and a satellite's records a `Vec` in file order.
//! The reader is `geomaticist_gnss_rinex::nav::read_nav`; orbit models are in the engine.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::system::Sat;

/// Field names of a Keplerian record, in record order.
pub const KEP_FIELDS: [&str; 26] = [
    "toc", "f0", "f1", "f2", "iode", "crs", "deln", "m0", "cuc", "e", "cus", "sqrta", "toes",
    "cic", "omg0", "cis", "i0", "crc", "omg", "omgd", "idot", "week", "svh", "tgd", "toe", "acc",
];

/// One Keplerian ephemeris (GPS / QZSS LNAV, Galileo I/NAV and F/NAV, BeiDou D1/D2), as the reader
/// stores it: times in GPS seconds (BeiDou's toc and toe moved from BDT by +14 s), the rest
/// as broadcast (s, m, rad, rad/s), Galileo health zeroed when its DVS/HS bits are clear.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[allow(missing_docs)] // field names are the RINEX / ICD symbols, listed in KEP_FIELDS
pub struct KeplerEph {
    pub toc: f64,
    pub f0: f64,
    pub f1: f64,
    pub f2: f64,
    pub iode: f64,
    pub crs: f64,
    pub deln: f64,
    pub m0: f64,
    pub cuc: f64,
    pub e: f64,
    pub cus: f64,
    pub sqrta: f64,
    pub toes: f64,
    pub cic: f64,
    pub omg0: f64,
    pub cis: f64,
    pub i0: f64,
    pub crc: f64,
    pub omg: f64,
    pub omgd: f64,
    pub idot: f64,
    pub week: f64,
    pub svh: f64,
    pub tgd: f64,
    pub toe: f64,
    pub acc: f64,
}

impl KeplerEph {
    /// The record in [`KEP_FIELDS`] order.
    #[must_use]
    pub fn to_array(&self) -> [f64; 26] {
        [
            self.toc, self.f0, self.f1, self.f2, self.iode, self.crs, self.deln, self.m0, self.cuc,
            self.e, self.cus, self.sqrta, self.toes, self.cic, self.omg0, self.cis, self.i0,
            self.crc, self.omg, self.omgd, self.idot, self.week, self.svh, self.tgd, self.toe,
            self.acc,
        ]
    }

    /// A field by its name in [`KEP_FIELDS`].
    #[must_use]
    pub fn field(&self, name: &str) -> Option<f64> {
        let i = KEP_FIELDS.iter().position(|f| *f == name)?;
        Some(self.to_array()[i])
    }
}

/// One GLONASS ephemeris: tb as GPS seconds (UTC epoch + leap seconds), -`TauN`, `GammaN`, position,
/// velocity, lunisolar acceleration (m, m/s, m/s²), health, frequency channel.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GloEph {
    /// tb in GPS seconds.
    pub toe: f64,
    /// `TauN` (s), sign flipped from RINEX's -`TauN`.
    pub taun: f64,
    /// `GammaN`.
    pub gamn: f64,
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Acceleration (m/s²).
    pub acc: [f64; 3],
    /// Health.
    pub svh: f64,
    /// Frequency channel number.
    pub fcn: f64,
}

/// Navigation data from one or more files.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct NavData {
    /// Keplerian records per satellite, in file order.
    pub kep: BTreeMap<Sat, Vec<KeplerEph>>,
    /// GLONASS records per satellite, in file order.
    pub glo: BTreeMap<Sat, Vec<GloEph>>,
    /// RINEX 4 records left out (`"G CNAV"`, `"STO"`, `"unreadable"` ...) -> how many.
    pub skipped: BTreeMap<String, usize>,
    /// (file name, error) of files that could not be read.
    pub bad: Vec<(String, String)>,
}

impl NavData {
    /// GLONASS slot -> frequency channel: the most common value broadcast (ties: the smallest).
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn glo_fcn(&self) -> BTreeMap<Sat, i32> {
        let mut out = BTreeMap::new();
        for (s, recs) in &self.glo {
            let mut cnt: BTreeMap<i64, usize> = BTreeMap::new();
            for r in recs {
                *cnt.entry(r.fcn as i64).or_default() += 1; // astype(int): toward zero
            }
            let mut best: Option<(i64, usize)> = None;
            for (v, n) in cnt {
                if best.is_none_or(|(_, bn)| n > bn) {
                    best = Some((v, n));
                }
            }
            if let Some((v, _)) = best {
                out.insert(*s, i32::try_from(v).unwrap_or(0));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn glo(fcn: f64) -> GloEph {
        GloEph {
            toe: 0.0,
            taun: 0.0,
            gamn: 0.0,
            pos: [0.0; 3],
            vel: [0.0; 3],
            acc: [0.0; 3],
            svh: 0.0,
            fcn,
        }
    }

    #[test]
    fn fcn_mode_with_ties_to_smallest() {
        let mut n = NavData::default();
        let r: Sat = "R01".parse().unwrap();
        n.glo
            .insert(r, vec![glo(1.0), glo(-4.0), glo(1.0), glo(-4.0), glo(2.0)]);
        assert_eq!(n.glo_fcn()[&r], -4);
        let k = KeplerEph {
            toc: 1.0,
            f0: 0.0,
            f1: 0.0,
            f2: 0.0,
            iode: 0.0,
            crs: 0.0,
            deln: 0.0,
            m0: 0.0,
            cuc: 0.0,
            e: 0.01,
            cus: 0.0,
            sqrta: 5153.0,
            toes: 0.0,
            cic: 0.0,
            omg0: 0.0,
            cis: 0.0,
            i0: 0.0,
            crc: 0.0,
            omg: 0.0,
            omgd: 0.0,
            idot: 0.0,
            week: 2333.0,
            svh: 0.0,
            tgd: 0.0,
            toe: 2.0,
            acc: 2.0,
        };
        assert_eq!(k.field("sqrta"), Some(5153.0));
        assert_eq!(k.to_array()[24], 2.0);
    }
}
