//! Small helpers shared by the decoders: header lines, GPS time to civil time, big- and
//! little-endian field readers, latin-1 text.

use geomaticist_gnss_core::time::{DateTime, GPS0_ORDINAL, ordinal};

/// Microseconds in a GPS week.
pub const WEEK_US: i64 = 604_800 * 1_000_000;

/// `f"{content:<60.60}{label:<20}"`: a RINEX header line.
#[must_use]
pub fn hline(content: &str, label: &str) -> String {
    format!("{content:<60.60}{label:<20}")
}

/// `GPS0 + timedelta(weeks=week, milliseconds=ms)`.
#[must_use]
pub fn gps_week_ms(week: i64, ms: i64) -> DateTime {
    DateTime::from_micros_since_gps0(week * WEEK_US + ms * 1000)
}

/// `GPS0 + timedelta(weeks=week, seconds=sec)` with whole seconds.
#[must_use]
pub fn gps_week_s(week: i64, sec: i64) -> DateTime {
    DateTime::from_micros_since_gps0(week * WEEK_US + sec * 1_000_000)
}

/// `GPS0 + timedelta(seconds=float)` (also with whole weeks added): `timedelta` keeps the whole seconds and
/// rounds the fraction to microseconds, half to even.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn gps_week_fsec(week: i64, sec: f64) -> DateTime {
    let whole = sec.trunc();
    let us = geomaticist_gnss_core::text::round_half_even((sec - whole) * 1e6) as i64;
    DateTime::from_micros_since_gps0(week * WEEK_US + whole as i64 * 1_000_000 + us)
}

/// Microseconds since GPS0 of midnight at the start of a civil date (`datetime(y, m, d)`); `None` when
/// the date is invalid.
#[must_use]
pub fn date_us(year: i64, month: i64, day: i64) -> Option<i64> {
    Some((ordinal(year, month, day)? - GPS0_ORDINAL) * 86_400_000_000)
}

/// `(t - GPS0).days`, floor.
#[must_use]
pub fn days_since_gps0(t: &DateTime) -> i64 {
    t.micros_since_gps0().div_euclid(86_400_000_000)
}

/// `t.second + t.microsecond / 1e6` (the float the writers format with `.7f`).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn sec_f(t: &DateTime) -> f64 {
    t.second as f64 + t.micro as f64 / 1e6
}

/// `f"{t:%Y %m %d %H %M %S}"`.
#[must_use]
pub fn ymdhms(t: &DateTime) -> String {
    format!(
        "{:04} {:02} {:02} {:02} {:02} {:02}",
        t.year, t.month, t.day, t.hour, t.minute, t.second
    )
}

/// Bytes of a string whose characters are all latin-1 (what the decoders read with `decode("latin-1")`);
/// any other character becomes `?`.
#[must_use]
pub fn latin1_bytes(s: &str) -> Vec<u8> {
    s.chars()
        .map(|c| u8::try_from(u32::from(c)).unwrap_or(b'?'))
        .collect()
}

/// Latin-1 decoding.
#[must_use]
pub fn latin1(b: &[u8]) -> String {
    b.iter().map(|&c| char::from(c)).collect()
}

/// A latin-1 string with surrounding whitespace removed (whitespace incl. `\x1c`-`\x1f`, `\x85`, `\xa0`).
#[must_use]
pub fn py_strip(s: &str) -> &str {
    s.trim_matches(py_space)
}

/// Whether one latin-1 character is whitespace.
#[must_use]
pub fn py_space(c: char) -> bool {
    matches!(
        c,
        ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c' | '\x1c'..='\x1f' | '\u{85}' | '\u{a0}'
    )
}

/// A decimal integer string parsed (surrounding whitespace, sign, digits, `_` between digits).
#[must_use]
pub fn py_int(s: &str) -> Option<i64> {
    geomaticist_gnss_core::time::py_int(&latin1_bytes(s))
}

/// A float from a text field: whitespace around, then a decimal or `inf` / `nan`.
#[must_use]
pub fn py_float(s: &str) -> Option<f64> {
    geomaticist_gnss_core::text::py_float(&latin1_bytes(s))
}

/// Field readers over a byte slice; out-of-range reads give 0 (callers check record lengths first).
pub mod be {
    /// Big-endian u16.
    #[must_use]
    pub fn u16(b: &[u8], o: usize) -> i64 {
        b.get(o..o + 2)
            .map_or(0, |x| i64::from(u16::from_be_bytes([x[0], x[1]])))
    }
    /// Big-endian u32.
    #[must_use]
    pub fn u32(b: &[u8], o: usize) -> i64 {
        b.get(o..o + 4).map_or(0, |x| {
            i64::from(u32::from_be_bytes([x[0], x[1], x[2], x[3]]))
        })
    }
    /// Big-endian f64.
    #[must_use]
    pub fn f64(b: &[u8], o: usize) -> f64 {
        b.get(o..o + 8).map_or(0.0, |x| {
            let mut a = [0u8; 8];
            a.copy_from_slice(x);
            f64::from_be_bytes(a)
        })
    }
}

/// Little-endian field readers.
pub mod le {
    /// u16.
    #[must_use]
    pub fn u16(b: &[u8], o: usize) -> i64 {
        b.get(o..o + 2)
            .map_or(0, |x| i64::from(u16::from_le_bytes([x[0], x[1]])))
    }
    /// i16.
    #[must_use]
    pub fn i16(b: &[u8], o: usize) -> i64 {
        b.get(o..o + 2)
            .map_or(0, |x| i64::from(i16::from_le_bytes([x[0], x[1]])))
    }
    /// u32.
    #[must_use]
    pub fn u32(b: &[u8], o: usize) -> i64 {
        b.get(o..o + 4).map_or(0, |x| {
            i64::from(u32::from_le_bytes([x[0], x[1], x[2], x[3]]))
        })
    }
    /// f32 as f64 (exact).
    #[must_use]
    pub fn f32(b: &[u8], o: usize) -> f64 {
        b.get(o..o + 4).map_or(0.0, |x| {
            f64::from(f32::from_le_bytes([x[0], x[1], x[2], x[3]]))
        })
    }
    /// f64.
    #[must_use]
    pub fn f64(b: &[u8], o: usize) -> f64 {
        b.get(o..o + 8).map_or(0.0, |x| {
            let mut a = [0u8; 8];
            a.copy_from_slice(x);
            f64::from_le_bytes(a)
        })
    }
}

/// WGS84 geodetic (radians) to ECEF, as the decoders compute it (`math.sin(la) ** 2` = `sin * sin`).
#[must_use]
pub fn llh_rad_to_xyz(la: f64, lo: f64, h: f64) -> (f64, f64, f64) {
    let a = 6_378_137.0;
    let f = 1.0 / 298.257_223_563;
    let e2 = f * (2.0 - f);
    let s = la.sin();
    let n = a / (1.0 - e2 * s.powi(2)).sqrt();
    (
        (n + h) * la.cos() * lo.cos(),
        (n + h) * la.cos() * lo.sin(),
        (n * (1.0 - e2) + h) * la.sin(),
    )
}

/// The same from degrees (`math.radians`).
#[must_use]
pub fn llh_deg_to_xyz(lat: f64, lon: f64, h: f64) -> (f64, f64, f64) {
    llh_rad_to_xyz(lat.to_radians(), lon.to_radians(), h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_line() {
        let l = hline("abc", "COMMENT");
        assert_eq!(l.len(), 80);
        assert!(l.starts_with("abc "));
        assert_eq!(&l[60..67], "COMMENT");
        assert_eq!(hline(&"x".repeat(70), "L").len(), 80);
    }

    #[test]
    fn gps_times() {
        let t = gps_week_ms(1, 86_400_000 + 1500);
        assert_eq!((t.year, t.month, t.day), (1980, 1, 14));
        assert_eq!((t.second, t.micro), (1, 500_000));
        assert!((sec_f(&t) - 1.5).abs() < 1e-12);
    }
}
