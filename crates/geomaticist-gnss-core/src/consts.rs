//! Physical constants and carrier frequencies.
//!
//! Speed of light, Earth rotation, the WGS 84 ellipsoid, GPS and BeiDou week constants, carrier
//! frequencies per system and band, GLONASS channel spacing and the dual-frequency pairs.

use crate::system::GnssSystem;

/// Speed of light (m/s).
pub const CLIGHT: f64 = 299_792_458.0;
/// Earth rotation rate, WGS 84 / IS-GPS (rad/s).
pub const OMGE: f64 = 7.292_115_146_7e-5;
/// WGS 84 semi-major axis (m).
pub const RE_WGS84: f64 = 6_378_137.0;
/// WGS 84 flattening.
pub const FE_WGS84: f64 = 1.0 / 298.257_223_563;
/// Seconds in a GPS week.
pub const WEEK: f64 = 604_800.0;
/// BDT = GPST - 14 s.
pub const BDT_GPST: f64 = 14.0;
/// BDT week 0 = GPS week 1356.
pub const BDS_WEEK0: i64 = 1356;

/// GLONASS FDMA channel spacing (Hz) of G1 (band 1) and G2 (band 2).
#[must_use]
pub fn glo_df(band: u8) -> Option<f64> {
    match band {
        1 => Some(0.5625e6),
        2 => Some(0.4375e6),
        _ => None,
    }
}

/// Carrier frequency (Hz) of a RINEX band number for a system, without the GLONASS channel
/// offset; `None` where the table has no entry.
///
/// Systems without a row of their own (NavIC, SBAS) use the GPS row.
#[must_use]
pub fn base_freq(sys: GnssSystem, band: u8) -> Option<f64> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Navic, Qzss, Sbas};
    let f = match (sys, band) {
        (Gps | Navic | Sbas | Qzss | Galileo | BeiDou, 1) => 1575.42e6,
        (Gps | Navic | Sbas | Qzss, 2) => 1227.60e6,
        (Gps | Navic | Sbas | Qzss | Galileo | BeiDou, 5) => 1176.45e6,
        (Qzss | Galileo, 6) => 1278.75e6,
        (Galileo | BeiDou, 7) => 1207.14e6,
        (Galileo | BeiDou, 8) => 1191.795e6,
        (BeiDou, 2) => 1561.098e6,
        (BeiDou, 6) => 1268.52e6,
        (Glonass, 1) => 1602.0e6,
        (Glonass, 2) => 1246.0e6,
        (Glonass, 3) => 1202.025e6,
        _ => return None,
    };
    Some(f)
}

/// Dual-frequency pairs, most wanted first (RINEX 3 band numbers; BeiDou B1I is band 2, B3I band 6).
/// Empty for systems the engine does not pair (NavIC, SBAS).
#[must_use]
pub fn pairs(sys: GnssSystem) -> &'static [(u8, u8)] {
    match sys {
        GnssSystem::Gps | GnssSystem::Qzss => &[(1, 2), (1, 5)],
        GnssSystem::Glonass => &[(1, 2)],
        GnssSystem::Galileo => &[(1, 5), (1, 8), (1, 7)],
        GnssSystem::BeiDou => &[(2, 6), (2, 7), (1, 5), (2, 5)],
        GnssSystem::Navic | GnssSystem::Sbas => &[],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_matches_python() {
        assert_eq!(base_freq(GnssSystem::Gps, 2), Some(1227.60e6));
        assert_eq!(base_freq(GnssSystem::BeiDou, 6), Some(1268.52e6));
        assert_eq!(base_freq(GnssSystem::Qzss, 6), Some(1278.75e6));
        assert_eq!(base_freq(GnssSystem::Gps, 6), None);
        assert_eq!(base_freq(GnssSystem::Sbas, 1), Some(1575.42e6)); // GPS row fallback
        assert_eq!(base_freq(GnssSystem::Glonass, 3), Some(1202.025e6));
        assert_eq!(pairs(GnssSystem::BeiDou)[0], (2, 6));
        assert!((OMGE - 7.292_115_146_7e-5).abs() < 1e-20);
    }
}
