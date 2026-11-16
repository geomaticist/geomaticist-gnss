//! Broadcast ephemerides in a Leica MDB log -> RINEX 3.04 navigation files.
//!
//! GPS, record 0x73 (129 bytes, little-endian; offsets from the type byte): 2 u8 PRN, 5 u16 week, 7 u8 URA
//! index, 19 u16 IODE (also written as IODC), 21 f32 TGD, 25 f64 M0, 33 f32 dn, 37 f64 e, 45 f64 sqrtA,
//! 53 f64 Omega0, 61 f64 i0, 69 f64 omega, 77 f32 OmegaDot, 81 f32 IDOT, 85 f32 Cuc, 89 f32 Cus, 93 f32 Crc,
//! 97 f32 Crs, 101 f32 Cic, 105 f32 Cis, 109 u32 toc, 113 u32 toe, 117 f32 af0, 121 f32 af1, 125 f32 af2.
//! GLONASS, record 0x8d (147 bytes): 2 u8 slot; 9 i8 frequency channel; 11 u16 week + 13 u32 ms: toe (GPS
//! time); 19 u8 health; 29 / 37 / 45 f64 X Y Z (m), 53 / 61 / 69 velocity, 77 / 85 / 93 acceleration, 101 f64
//! tauN (RINEX writes -tauN), 117 f64 gammaN, 125 u32 tk.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use crate::decimal::Dec;
use crate::util::{le, ymdhms};
use crate::{DecodeError, util};

const T_GPS: u8 = 0x73;
const T_GLO: u8 = 0x8d;
const URA: [f64; 16] = [
    2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0,
    6144.0,
];
/// GPS - UTC from 2017 on; the logs hold no UTC record we decode.
pub const LEAP: i64 = 18;

/// One ephemeris.
#[derive(Debug, Clone)]
pub struct Eph {
    /// "G05" / "R04".
    pub sat: String,
    /// Week.
    pub week: i64,
    /// GPS: toc (s of week); GLONASS: toe (ms of week).
    pub t: i64,
    /// GPS toe (dedupe key).
    pub toe: Option<i64>,
    /// GPS IODE (dedupe key).
    pub iode: Option<i64>,
    /// Values in RINEX order.
    pub vals: Vec<f64>,
}

#[allow(clippy::cast_precision_loss)]
fn gps_eph(b: &[u8]) -> Option<Eph> {
    if b.len() < 129 {
        return None;
    }
    let f = |o| le::f32(b, o);
    let d = |o| le::f64(b, o);
    let week = le::u16(b, 5);
    let iode = le::u16(b, 19);
    let (toc, toe) = (le::u32(b, 109), le::u32(b, 113));
    let ura = URA.get(usize::from(b[7])).copied().unwrap_or(6144.0);
    Some(Eph {
        sat: format!("G{:02}", b[2]),
        week,
        t: toc,
        toe: Some(toe),
        iode: Some(iode),
        vals: vec![
            f(117),
            f(121),
            f(125),
            iode as f64,
            f(97),
            f(33),
            d(25),
            f(85),
            d(37),
            f(89),
            d(45),
            toe as f64,
            f(101),
            d(53),
            f(105),
            d(61),
            f(93),
            d(69),
            f(77),
            f(81),
            1.0,
            week as f64,
            0.0,
            ura,
            0.0,
            f(21),
            iode as f64,
            toc as f64,
        ],
    })
}

#[allow(clippy::cast_precision_loss)]
fn glo_eph(b: &[u8]) -> Option<Eph> {
    if b.len() < 129 {
        return None;
    }
    let d = |o| le::f64(b, o);
    let (week, ms) = (le::u16(b, 11), le::u32(b, 13));
    let k = i64::from(b[9] as i8);
    let tk = le::u32(b, 125);
    Some(Eph {
        sat: format!("R{:02}", b[2]),
        week,
        t: ms,
        toe: None,
        iode: None,
        vals: vec![
            -d(101),
            d(117),
            tk as f64,
            d(29) / 1000.0,
            d(53) / 1000.0,
            d(77) / 1000.0,
            f64::from(b[19]),
            d(37) / 1000.0,
            d(61) / 1000.0,
            d(85) / 1000.0,
            k as f64,
            d(45) / 1000.0,
            d(69) / 1000.0,
            d(93) / 1000.0,
            0.0,
        ],
    })
}

/// `(gps, glonass)`: the decoded ephemerides, one per (satellite, toe), in log order; records whose week is
/// more than a week from `week_near` (a stale or empty record) are left out.
#[must_use]
pub fn ephemerides(data: &[u8], week_near: Option<i64>) -> (Vec<Eph>, Vec<Eph>) {
    let (mut gps, mut glo) = (Vec::new(), Vec::new());
    let mut seen: HashSet<(String, i64, i64, Option<i64>)> = HashSet::new();
    for (_p, typ, body) in super::records(data) {
        let e = match typ {
            T_GPS => gps_eph(body),
            T_GLO => glo_eph(body),
            _ => None,
        };
        let Some(e) = e else { continue };
        if let Some(w) = week_near.filter(|w| *w != 0)
            && (e.week - w).abs() > 1
        {
            continue;
        }
        let key = (e.sat.clone(), e.week, e.toe.unwrap_or(e.t), e.iode);
        if !seen.insert(key) {
            continue;
        }
        if e.sat.starts_with('G') {
            gps.push(e);
        } else {
            glo.push(e);
        }
    }
    (gps, glo)
}

/// D19.12 of the exact double, half away from zero (ties too, e.g. GLONASS km on a 2^-11 grid).
#[must_use]
pub fn num(v: f64) -> String {
    if v == 0.0 {
        let s = if v.is_sign_negative() { '-' } else { ' ' };
        return format!("{s}0.000000000000D+00");
    }
    let x = Dec::exact(v);
    let mut e = x.adjusted();
    let mut m = x.scaleb(-e).round_at(-12);
    if m.adjusted() >= 1 {
        m = m.scaleb(-1).round_at(-12);
        e += 1;
    }
    let body = format!(
        "{}{}D{}{:02}",
        if m.neg { '-' } else { ' ' },
        m.abs_fixed(),
        if e < 0 { '-' } else { '+' },
        e.abs()
    );
    format!("{body:>19}")
}

fn lines(
    sat: &str,
    t: &geomaticist_gnss_core::time::DateTime,
    vals: &[f64],
    out: &mut Vec<String>,
) {
    let mut first = format!("{sat} {}", ymdhms(t));
    for v in vals.iter().take(3) {
        first.push_str(&num(*v));
    }
    out.push(first);
    for chunk in vals.get(3..).unwrap_or(&[]).chunks(4) {
        let mut l = String::from("    ");
        for v in chunk {
            l.push_str(&num(*v));
        }
        out.push(l);
    }
}

/// Writes `<stem>.<yy>n` (GPS) and `<stem>.<yy>g` (GLONASS) for the ephemerides the log holds.
#[allow(clippy::too_many_arguments)]
pub fn write_nav(
    data: &[u8],
    out_dir: &Path,
    stem: &str,
    year: i64,
    program: &str,
    run_by: &str,
    date: &str,
    week_near: Option<i64>,
) -> Result<Vec<PathBuf>, DecodeError> {
    let (gps, glo) = ephemerides(data, week_near);
    let mut out = Vec::new();
    for (letter, name, ext, mut ephs) in [('G', "GPS", 'n', gps), ('R', "GLONASS", 'g', glo)] {
        if ephs.is_empty() {
            continue;
        }
        ephs.sort_by(|a, b| (&a.sat, a.week, a.t).cmp(&(&b.sat, b.week, b.t)));
        let mut l = vec![
            format!(
                "{:>9}{:11}{:<20}{:<20}RINEX VERSION / TYPE",
                "3.04",
                "",
                "N: GNSS NAV DATA",
                format!("{letter}: {name}")
            ),
            format!("{program:<20.20}{run_by:<20.20}{date} PGM / RUN BY / DATE"),
            format!("{LEAP:6}{:54}LEAP SECONDS", ""),
            format!("{:60}END OF HEADER", ""),
        ];
        for e in &ephs {
            let t = if letter == 'G' {
                util::gps_week_s(e.week, e.t)
            } else {
                geomaticist_gnss_core::time::DateTime::from_micros_since_gps0(
                    e.week * util::WEEK_US + e.t * 1000 - LEAP * 1_000_000,
                )
            };
            lines(&e.sat, &t, &e.vals, &mut l);
        }
        // CRLF line ends, as text files written on Windows have them
        let mut text = l.join("\r\n");
        text.push_str("\r\n");
        let p = out_dir.join(format!("{stem}.{:02}{ext}", year.rem_euclid(100)));
        crate::write_atomic(&p, &util::latin1_bytes(&text))?;
        out.push(p);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn number_format() {
        assert_eq!(num(0.0), " 0.000000000000D+00");
        assert_eq!(num(-0.0), "-0.000000000000D+00");
        assert_eq!(num(1.0), " 1.000000000000D+00");
        assert_eq!(num(-2.5e-5), "-2.500000000000D-05");
        assert_eq!(num(18.0), " 1.800000000000D+01");
    }
}
