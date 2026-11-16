//! Signal identification table for rec35 sub19 measurement blocks.
//!
//! `(rt_sattype, blockType, trackType)` -> RINEX 3 system letter, obs code suffix, carrier frequency. GLONASS
//! frequencies depend on the FDMA channel k.

/// Speed of light, m/s.
pub const C_LIGHT: f64 = 299_792_458.0;

/// `rt_sattype` -> (RINEX system letter, PRN offset added to svId).
#[must_use]
pub fn sattype(st: i64) -> (char, i64) {
    match st {
        0 => ('G', 0),
        1 => ('S', -100), // SBAS: svId 131 -> S31
        2 => ('R', 0),
        3 => ('E', 0),
        4 => ('J', -192), // QZSS: svId 195 -> J03
        10 => ('C', 0),
        _ => ('?', 0),
    }
}

/// A carrier: a fixed frequency or a GLONASS FDMA band.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Freq {
    /// Fixed frequency, Hz.
    Hz(f64),
    /// GLONASS G1 (1602 MHz + k * 0.5625 MHz).
    G1,
    /// GLONASS G2 (1246 MHz + k * 0.4375 MHz).
    G2,
}

/// `(sattype, blockType, trackType)` -> (RINEX code suffix, carrier).
#[must_use]
pub fn signal(st: i64, block: i64, track: i64) -> Option<(&'static str, Freq)> {
    Some(match (st, block, track) {
        (0 | 1 | 4, 0, 0) => ("1C", Freq::Hz(1575.42e6)),
        (0, 1, 2) => ("2W", Freq::Hz(1227.60e6)),
        (0 | 4, 1, 5) => ("2X", Freq::Hz(1227.60e6)), // QZSS: +0.25 cycle phase shift (header)
        (0 | 4 | 10, 2, 8) => ("5X", Freq::Hz(1176.45e6)),
        (2, 0, 0) => ("1C", Freq::G1),
        (2, 1, 0) => ("2C", Freq::G2),
        (3, 0, 23) => ("1X", Freq::Hz(1575.42e6)),
        (3, 4, 14) => ("8X", Freq::Hz(1191.795e6)),
        (10, 0, 20) => ("1X", Freq::Hz(1575.42e6)),
        (10, 3, 28) => ("7I", Freq::Hz(1207.14e6)),
        (10, 6, 26) => ("2I", Freq::Hz(1561.098e6)),
        (10, 7, 29) => ("6I", Freq::Hz(1268.52e6)),
        _ => return None,
    })
}

/// Carrier frequency, Hz.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn frequency(f: Freq, fdma_k: i64) -> f64 {
    match f {
        Freq::Hz(v) => v,
        Freq::G1 => 1602.0e6 + fdma_k as f64 * 0.5625e6,
        Freq::G2 => 1246.0e6 + fdma_k as f64 * 0.4375e6,
    }
}

/// Wavelength, m.
#[must_use]
pub fn wavelength(f: Freq, fdma_k: i64) -> f64 {
    C_LIGHT / frequency(f, fdma_k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table() {
        assert_eq!(signal(0, 1, 2).map(|s| s.0), Some("2W"));
        assert!(signal(0, 9, 9).is_none());
        assert_eq!(frequency(Freq::G1, -7), 1602.0e6 - 7.0 * 0.5625e6);
        assert_eq!(sattype(4), ('J', -192));
    }
}
