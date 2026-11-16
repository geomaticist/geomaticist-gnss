//! RINEX 2.11 navigation files (GPS .yyn, GLONASS .yyg) from runpkr00's .TGD.
//! Records: 0x15 GPS ephemeris, 0x0F
//! ionosphere / UTC, 0x1C GLONASS ephemeris; little-endian, angles in semicircles.

use std::collections::{HashMap, HashSet};

use geomaticist_gnss_core::time::DateTime;

use super::records as rec;
use crate::decimal::Dec;
use crate::util::{gps_week_s, hline, le};

/// Semicircles -> radians (IS-GPS-200's pi).
#[allow(clippy::approx_constant)] // the interface specification's value, not f64::consts::PI
const GPS_PI: f64 = 3.141_592_653_589_8;
const URA_M: [f64; 16] = [
    2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0,
    6144.0,
];

/// Fortran D(digits+7).digits, `' 0.123456789012D+01'`: the exact binary value rounded
/// half up.
#[must_use]
pub fn dfmt(x: f64, digits: i32) -> String {
    let nd = usize::try_from(digits).unwrap_or(0);
    if x == 0.0 {
        return format!(" 0.{}D+00", "0".repeat(nd));
    }
    let d = Dec::exact(x.abs());
    let mut e = d.adjusted() + 1;
    let mut m = d.scaleb(-e).round_at(-digits);
    if m.adjusted() >= 0 {
        // rounding carried into the next power
        e += 1;
        m = d.scaleb(-e).round_at(-digits);
    }
    format!(
        "{}{}D{}{:02}",
        if x < 0.0 { '-' } else { ' ' },
        m.abs_fixed(),
        if e < 0 { '-' } else { '+' },
        e.abs()
    )
}

fn d12(x: f64) -> String {
    dfmt(x, 12)
}

/// One GPS ephemeris (0x15).
#[derive(Debug, Clone)]
struct GpsEph {
    prn: u8,
    week: i64,
    iodc: i64,
    iode: u8,
    ttm: i64,
    toc: i64,
    toe: i64,
    accuracy: f64,
    l2codes: u8,
    health: u8,
    /// tgd af2 af1 af0 crs dn m0 cuc e cus sqrta cic omega0 cis i0 crc omega omegadot idot
    d: [f64; 19],
}

// offset of each double and whether it is in semicircles, in the order of `GpsEph::d`
const GPS_DOUBLES: [(usize, bool); 19] = [
    (24, false),
    (32, false),
    (40, false),
    (48, false),
    (56, false),
    (64, true),
    (72, true),
    (80, true),
    (88, false),
    (96, true),
    (104, false),
    (112, true),
    (120, true),
    (128, true),
    (136, true),
    (144, false),
    (152, true),
    (160, true),
    (168, true),
];

fn gps_ephemerides(tgd: &[u8]) -> Result<Vec<GpsEph>, String> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for (i, t, _s, ln) in rec::records(tgd)? {
        if t != 0x15 {
            continue;
        }
        let r = &tgd[i..i + ln];
        if r.len() < 181 {
            return Err("ephemeris record too short".into());
        }
        let mut d = [0.0; 19];
        for (k, (o, sc)) in GPS_DOUBLES.iter().enumerate() {
            d[k] = le::f64(r, *o) * if *sc { GPS_PI } else { 1.0 };
        }
        let mut e = GpsEph {
            prn: r[4],
            week: le::u16(r, 6),
            iodc: le::u16(r, 8),
            iode: r[11],
            ttm: le::u32(r, 12),
            toc: le::u32(r, 16),
            toe: le::u32(r, 20),
            accuracy: URA_M[usize::from((r[177] >> 3).min(15))], // inferred, not confirmed
            l2codes: r[178],
            health: r[180],
            d,
        };
        if e.toc < e.ttm - 302_400 {
            // toc / toe early in the next week than the one sent in
            e.week += 1;
        }
        if seen.insert((e.prn, e.toe, e.iode)) {
            out.push(e);
        }
    }
    Ok(out)
}

/// The file's last ionosphere / UTC record: (alpha, beta, a1, a0, tot, leap).
#[allow(clippy::type_complexity)]
fn ion_utc(tgd: &[u8]) -> Result<Option<([f64; 4], [f64; 4], f64, f64, f64, f64)>, String> {
    let mut out = None;
    for (i, t, _s, _ln) in rec::records(tgd)? {
        if t == 0x0F {
            let v: Vec<f64> = (0..12).map(|k| le::f64(tgd, i + 4 + 8 * k)).collect();
            out = Some((
                [v[0], v[1], v[2], v[3]],
                [v[4], v[5], v[6], v[7]],
                v[8],
                v[9],
                v[10],
                v[11],
            ));
        }
    }
    Ok(out)
}

#[derive(Debug, Clone)]
struct GloEph {
    slot: u8,
    tb: i64, // micros since GPS0 (UTC clock reading)
    tk: i64,
    freq: i64,
    health: i64,
    /// x vx ax y vy ay z vz az tau_c tau_n gamma_n
    d: [f64; 12],
    order: usize,
}

const DAY_US: i64 = 86_400_000_000;
const HOUR_US: i64 = 3_600_000_000;

/// The last record of each (slot, tb) is kept.
fn glonass_ephemerides(tgd: &[u8], leap: i64) -> Result<Vec<GloEph>, String> {
    let mut keep: HashMap<(u8, i64), GloEph> = HashMap::new();
    for (i, t, _s, ln) in rec::records(tgd)? {
        if t != 0x1C || tgd.get(i + 4) != Some(&2) {
            continue;
        }
        let r = &tgd[i..i + ln];
        if r.len() < 27 + 13 * 8 {
            return Err("GLONASS ephemeris record too short".into());
        }
        let week = le::u16(r, 6);
        let tk_tow = le::u32(r, 8);
        let tb_i = i64::from(r[14]);
        let mut d = [0.0; 12];
        for (k, v) in d.iter_mut().enumerate() {
            *v = le::f64(r, 27 + 8 * k);
        }
        // tb is the quarter hour of the Moscow day; its date is the UTC day of the frame time tk
        let tk_utc = gps_week_s(week, tk_tow - leap).micros_since_gps0();
        let day0 = tk_utc.div_euclid(DAY_US) * DAY_US;
        let mut tb = day0 + (tb_i * 900 - 3 * 3600) * 1_000_000;
        if tb - tk_utc > 12 * HOUR_US {
            tb -= DAY_US;
        } else if tk_utc - tb > 12 * HOUR_US {
            tb += DAY_US;
        }
        let e = GloEph {
            slot: r[5],
            tb,
            tk: (tk_tow - leap).rem_euclid(86_400),
            freq: i64::from(r[23] as i8),
            health: i64::from(r[24] != 0),
            d,
            order: i,
        };
        keep.insert((e.slot, e.tb), e);
    }
    Ok(keep.into_values().collect())
}

fn line_start(prn: u8, t: &DateTime) -> String {
    #[allow(clippy::cast_precision_loss)]
    let sec = t.second as f64;
    format!(
        "{prn:2} {:02} {:02} {:02} {:02} {:02}{sec:5.1}",
        t.year.rem_euclid(100),
        t.month,
        t.day,
        t.hour,
        t.minute
    )
}

fn finish(lines: &[String], line_end: &str) -> Vec<u8> {
    let mut s = lines.join(line_end);
    s.push_str(line_end);
    crate::util::latin1_bytes(&s)
}

/// The GPS navigation file text, or `None` without ephemerides.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub fn gps_nav(
    tgd: &[u8],
    program: &str,
    run_by: &str,
    date: &str,
    line_end: &str,
) -> Result<Option<Vec<u8>>, String> {
    let mut eph = gps_ephemerides(tgd)?;
    if eph.is_empty() {
        return Ok(None);
    }
    eph.sort_by_key(|e| (e.week, e.toc, e.prn));
    let iu = ion_utc(tgd)?;
    let mut l = vec![
        hline(
            &format!(
                "{:>9}{:11}{:<20}{:<20}",
                "2.11", "", "NAVIGATION DATA", "GPS(GPS)"
            ),
            "RINEX VERSION / TYPE",
        ),
        hline(
            &format!("{program:<20.20}{run_by:<20.20}{date}"),
            "PGM / RUN BY / DATE",
        ),
    ];
    if let Some((alpha, beta, a1, a0, tot, leap)) = iu {
        let d4 = |x: f64| format!("{:>12}", dfmt(x, 4).trim());
        l.push(hline(
            &format!("  {}", alpha.iter().map(|x| d4(*x)).collect::<String>()),
            "ION ALPHA",
        ));
        l.push(hline(
            &format!("  {}", beta.iter().map(|x| d4(*x)).collect::<String>()),
            "ION BETA",
        ));
        l.push(hline(
            &format!("   {}{}{:9}{:9}", d12(a0), d12(a1), tot as i64, eph[0].week),
            "DELTA-UTC: A0,A1,T,W",
        ));
        l.push(hline(&format!("{:6}", leap as i64), "LEAP SECONDS"));
    }
    l.push(hline("", "END OF HEADER"));
    for e in &eph {
        let t = gps_week_s(e.week, e.toc);
        let d = &e.d;
        l.push(format!(
            "{}{}{}{}",
            line_start(e.prn, &t),
            d12(d[3]),
            d12(d[2]),
            d12(d[1])
        ));
        let rows: [[f64; 4]; 7] = [
            [f64::from(e.iode), d[4], d[5], d[6]],
            [d[7], d[8], d[9], d[10]],
            [e.toe as f64, d[11], d[12], d[13]],
            [d[14], d[15], d[16], d[17]],
            [d[18], f64::from(e.l2codes), e.week as f64, 0.0],
            [e.accuracy, f64::from(e.health), d[0], e.iodc as f64],
            [e.ttm as f64, 4.0, 0.0, 0.0],
        ];
        for row in rows {
            l.push(format!(
                "   {}",
                row.iter().map(|x| d12(*x)).collect::<String>()
            ));
        }
    }
    Ok(Some(finish(&l, line_end)))
}

/// The GLONASS navigation file text, or `None` without ephemerides.
#[allow(clippy::cast_precision_loss)]
pub fn glonass_nav(
    tgd: &[u8],
    leap: i64,
    program: &str,
    run_by: &str,
    date: &str,
    line_end: &str,
) -> Result<Option<Vec<u8>>, String> {
    let mut eph = glonass_ephemerides(tgd, leap)?;
    if eph.is_empty() {
        return Ok(None);
    }
    eph.sort_by_key(|e| (e.tb, e.slot));
    // the header's CORR TO SYSTEM TIME: the file's last record
    let last = eph
        .iter()
        .max_by_key(|e| e.order)
        .map_or((0, 0.0), |e| (e.tb, e.d[9]));
    let t0 = DateTime::from_micros_since_gps0(last.0);
    let mut l = vec![
        hline(
            &format!("{:>9}{:11}{:<20}{:<20}", "2.11", "", "GLONASS NAV DATA", ""),
            "RINEX VERSION / TYPE",
        ),
        hline(
            &format!("{program:<20.20}{run_by:<20.20}{date}"),
            "PGM / RUN BY / DATE",
        ),
        hline(
            &format!("{:6}{:6}{:6}   {}", t0.year, t0.month, t0.day, d12(-last.1)),
            "CORR TO SYSTEM TIME",
        ),
        hline(&format!("{leap:6}"), "LEAP SECONDS"),
        hline("", "END OF HEADER"),
    ];
    for e in &eph {
        let t = DateTime::from_micros_since_gps0(e.tb);
        let d = &e.d;
        l.push(format!(
            "{}{}{}{}",
            line_start(e.slot, &t),
            d12(-d[10]),
            d12(d[11]),
            d12(e.tk as f64)
        ));
        let rows = [
            [d[0], d[1], d[2], e.health as f64],
            [d[3], d[4], d[5], e.freq as f64],
            [d[6], d[7], d[8], 0.0],
        ];
        for row in rows {
            let txt: String = row
                .iter()
                .enumerate()
                .map(|(k, v)| d12(if k < 3 { v / 1000.0 } else { *v }))
                .collect();
            l.push(format!("   {txt}"));
        }
    }
    Ok(Some(finish(&l, line_end)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d_format() {
        assert_eq!(dfmt(0.0, 12), " 0.000000000000D+00");
        assert_eq!(dfmt(1.0, 12), " 0.100000000000D+01");
        assert_eq!(dfmt(-0.5, 12), "-0.500000000000D+00");
        assert_eq!(dfmt(1.862645149231e-8, 4).trim(), "0.1863D-07");
        assert_eq!(dfmt(0.99999999999999, 4), " 0.1000D+01");
    }
}
