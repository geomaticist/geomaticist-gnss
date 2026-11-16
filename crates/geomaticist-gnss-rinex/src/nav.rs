//! Navigation reader: RINEX 2 (GPS `N`, GLONASS `G`), RINEX 3 and RINEX 4 broadcast ephemerides.
//!
//! Includes the B1 fix (record-type dispatch for RINEX 4: `> EPH sys TYPE` heads, STO / EOP / ION
//! and the modern message types skipped and counted, a record of the wrong line count skipped).
//! Records go to [`geomaticist_gnss_core::NavData`] in file order, files in the order given.

use std::path::Path;

use geomaticist_gnss_core::consts::{BDS_WEEK0, BDT_GPST, WEEK};
use geomaticist_gnss_core::nav::{GloEph, KeplerEph, NavData};
use geomaticist_gnss_core::text::{py_float, py_splitlines, py_strip, sl};
use geomaticist_gnss_core::time::{DateTime, gps_seconds, py_int, utc_to_gps_seconds};
use geomaticist_gnss_core::{GnssSystem, Sat};

/// Message types whose fields the models read, per system (RINEX 4.01 Tables A5-A30).
#[must_use]
pub fn legacy_types(sys: u8) -> &'static [&'static str] {
    match sys {
        b'G' | b'J' | b'I' => &["LNAV"],
        b'E' => &["INAV", "FNAV"],
        b'C' => &["D1", "D2"],
        b'R' => &["FDMA"],
        b'S' => &["SBAS"],
        _ => &[],
    }
}

/// Line counts of RINEX 4 ephemeris records, head included.
#[must_use]
pub fn rnx4_lines(mtype: &str) -> Option<usize> {
    Some(match mtype {
        "LNAV" | "INAV" | "FNAV" | "D1" | "D2" => 8,
        "FDMA" => 5,
        "SBAS" => 4,
        "CNAV" | "CNV3" => 9,
        "CNV1" | "CNV2" => 10,
        _ => return None,
    })
}

#[derive(Debug)]
struct Unreadable;

/// One number field: strip, D/d -> E, empty -> 0.0, else parse as a float.
fn num(s: &[u8]) -> Result<f64, Unreadable> {
    let t = py_strip(s);
    if t.is_empty() {
        return Ok(0.0);
    }
    let v: Vec<u8> = t
        .iter()
        .map(|&b| if b == b'D' || b == b'd' { b'E' } else { b })
        .collect();
    py_float(&v).ok_or(Unreadable)
}

/// `n` fields of 19 characters from `start`.
fn nums(line: &[u8], start: usize, n: usize, out: &mut Vec<f64>) -> Result<(), Unreadable> {
    for i in 0..n {
        out.push(num(sl(line, start + 19 * i, start + 19 * (i + 1)))?);
    }
    Ok(())
}

fn count(nav: &mut NavData, key: impl Into<String>) {
    *nav.skipped.entry(key.into()).or_default() += 1;
}

/// Reads navigation files. A file that cannot be read is left out and named in `bad`; a record that
/// cannot be read is counted in `skipped["unreadable"]`. Missing paths are skipped silently.
pub fn read_nav<P: AsRef<Path>>(paths: &[P]) -> NavData {
    let mut nav = NavData::default();
    for p in paths {
        let p = p.as_ref();
        if !p.exists() {
            continue;
        }
        let name = p
            .file_name()
            .map_or_default(|n| n.to_string_lossy().into_owned());
        match crate::open::read_text(p) {
            Ok(text) => {
                let lines = py_splitlines(&text);
                if !lines.is_empty()
                    && let Err(e) = read_nav_lines(&lines, &mut nav)
                {
                    nav.bad.push((name, e));
                }
            }
            Err(e) => nav.bad.push((name, format!("OSError: {e}"))),
        }
    }
    nav
}

/// Reads navigation text already in memory (one file).
pub fn read_nav_text(text: &[u8], nav: &mut NavData) -> Result<(), String> {
    let lines = py_splitlines(text);
    if lines.is_empty() {
        return Ok(());
    }
    read_nav_lines(&lines, nav)
}

fn read_nav_lines(lines: &[&[u8]], nav: &mut NavData) -> Result<(), String> {
    let l0 = lines[0];
    let vtxt = sl(l0, 0, 9);
    let ver = if vtxt.is_empty() {
        0.0
    } else {
        py_float(vtxt).ok_or_else(|| {
            format!(
                "ValueError: could not convert string to float: {:?}",
                String::from_utf8_lossy(vtxt)
            )
        })?
    };
    let ftype = sl(l0, 20, 21).first().map_or(0, u8::to_ascii_uppercase);
    let fsys = sl(l0, 40, 41).first().map_or(0, u8::to_ascii_uppercase);
    let mut i = 0;
    while i < lines.len() && !contains(sl(lines[i], 60, usize::MAX), b"END OF HEADER") {
        i += 1;
    }
    i += 1;
    if ver >= 3.0 {
        while i < lines.len() {
            let ln = lines[i];
            if ln.starts_with(b">") {
                // RINEX 4: "> EPH G18 LNAV", "> STO G LNAV", "> EOP ...", "> ION ..."; the record is
                // the head's following lines up to the next head (the first carries the satellite,
                // the others start with four blanks)
                let head: Vec<&[u8]> = ln[1..]
                    .split(|b| geomaticist_gnss_core::text::py_isspace(*b))
                    .filter(|w| !w.is_empty())
                    .collect();
                i += 1;
                let mut body: Vec<&[u8]> = Vec::new();
                while i < lines.len() && !lines[i].starts_with(b">") {
                    if !body.is_empty() && !lines[i].starts_with(b"    ") {
                        break;
                    }
                    body.push(lines[i]);
                    i += 1;
                }
                let rtype = head
                    .first()
                    .map_or("?".to_owned(), |w| String::from_utf8_lossy(w).into_owned());
                let mtype = head
                    .get(2)
                    .map_or(String::new(), |w| String::from_utf8_lossy(w).into_owned());
                if rtype != "EPH" {
                    count(nav, rtype); // STO / EOP / ION: no orbit in them
                    continue;
                }
                let Some(&sysc) = body.first().and_then(|b| b.first()) else {
                    count(nav, "EPH ?");
                    continue;
                };
                if !b"GRECJIS".contains(&sysc) {
                    count(nav, "EPH ?");
                    continue;
                }
                let sc = char::from(sysc);
                if !legacy_types(sysc).contains(&mtype.as_str()) {
                    let m = if mtype.is_empty() {
                        "?"
                    } else {
                        mtype.as_str()
                    };
                    count(nav, format!("{sc} {m}"));
                    continue;
                }
                if Some(body.len()) != rnx4_lines(&mtype) {
                    count(nav, format!("{sc} {mtype} ({} lines)", body.len()));
                    continue;
                }
                record3(&body, nav);
            } else if geomaticist_gnss_core::text::is_blank(ln) || !b"GRECJIS".contains(&ln[0]) {
                i += 1;
            } else {
                let sysc = ln[0];
                let nlines = match sysc {
                    b'R' => {
                        if ver < 3.05 {
                            4
                        } else {
                            5
                        }
                    }
                    b'S' => 4,
                    _ => 8,
                };
                if i + nlines > lines.len() {
                    break;
                }
                record3(&lines[i..i + nlines], nav);
                i += nlines;
            }
        }
    } else {
        let sysc = if ftype == b'G' || fsys == b'R' {
            b'R'
        } else {
            b'G'
        };
        let nlines = if sysc == b'R' { 4 } else { 8 };
        while i + nlines <= lines.len() {
            let ln = lines[i];
            if geomaticist_gnss_core::text::is_blank(ln) {
                i += 1;
                continue;
            }
            let rec = &lines[i..i + nlines];
            i += nlines;
            let parsed = (|| -> Result<(i64, DateTime, Vec<f64>), Unreadable> {
                let prn = py_int(sl(ln, 0, 2)).ok_or(Unreadable)?;
                let yy = py_int(sl(ln, 3, 5)).ok_or(Unreadable)?;
                let y = if yy < 80 { 2000 + yy } else { 1900 + yy };
                let sec = py_float(sl(ln, 17, 22)).ok_or(Unreadable)?;
                let when = DateTime::new(
                    y,
                    py_int(sl(ln, 6, 8)).ok_or(Unreadable)?,
                    py_int(sl(ln, 9, 11)).ok_or(Unreadable)?,
                    py_int(sl(ln, 12, 14)).ok_or(Unreadable)?,
                    py_int(sl(ln, 15, 17)).ok_or(Unreadable)?,
                    sec.trunc() as i64,
                )
                .ok_or(Unreadable)?;
                let mut vals = Vec::with_capacity(31);
                nums(ln, 22, 3, &mut vals)?;
                for r in &rec[1..] {
                    nums(r, 3, 4, &mut vals)?;
                }
                Ok((prn, when, vals))
            })();
            match parsed {
                Ok((prn, when, vals)) => {
                    let sys = GnssSystem::from_letter(char::from(sysc)).unwrap_or(GnssSystem::Gps);
                    match u8::try_from(prn).ok().filter(|p| *p <= 99) {
                        Some(p) => store(nav, Sat::new(sys, p), &when, &vals),
                        None => count(nav, "unreadable"),
                    }
                }
                Err(Unreadable) => count(nav, "unreadable"),
            }
        }
    }
    Ok(())
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

/// One RINEX 3 / 4 ephemeris record (its lines).
fn record3(rec: &[&[u8]], nav: &mut NavData) {
    let ln = rec[0];
    let parsed = (|| -> Result<(Sat, DateTime, Vec<f64>), Unreadable> {
        let sys = GnssSystem::from_letter(char::from(ln[0])).ok_or(Unreadable)?;
        let prn = py_int(sl(ln, 1, 3)).ok_or(Unreadable)?;
        let prn = u8::try_from(prn)
            .ok()
            .filter(|p| *p <= 99)
            .ok_or(Unreadable)?;
        let field = |a, b| py_int(sl(ln, a, b)).ok_or(Unreadable);
        let (y, mo, d, h, mi) = (
            field(4, 8)?,
            field(9, 11)?,
            field(12, 14)?,
            field(15, 17)?,
            field(18, 20)?,
        );
        let s = py_float(sl(ln, 21, 23)).ok_or(Unreadable)?.trunc() as i64;
        let mut vals = Vec::with_capacity(31);
        nums(ln, 23, 3, &mut vals)?;
        for r in &rec[1..] {
            nums(r, 4, 4, &mut vals)?;
        }
        let when = DateTime::new(y, mo, d, h, mi, s).ok_or(Unreadable)?;
        Ok((Sat::new(sys, prn), when, vals))
    })();
    match parsed {
        Ok((sat, when, vals)) => store(nav, sat, &when, &vals),
        Err(Unreadable) => count(nav, "unreadable"),
    }
}

/// One record's values into the GLONASS or Keplerian tables.
fn store(nav: &mut NavData, sat: Sat, when: &DateTime, v: &[f64]) {
    if sat.sys == GnssSystem::Glonass {
        // v: -taun, +gamn, tk | x, vx, ax, health | y, vy, ay, fcn | z, vz, az, age  (km, km/s, km/s²)
        if v.len() < 15 {
            return;
        }
        let toe = utc_to_gps_seconds(when);
        nav.glo.entry(sat).or_default().push(GloEph {
            toe,
            taun: -v[0],
            gamn: v[1],
            pos: [v[3] * 1e3, v[7] * 1e3, v[11] * 1e3],
            vel: [v[4] * 1e3, v[8] * 1e3, v[12] * 1e3],
            acc: [v[5] * 1e3, v[9] * 1e3, v[13] * 1e3],
            svh: v[6],
            fcn: v[10],
        });
        return;
    }
    let kepler = matches!(
        sat.sys,
        GnssSystem::Gps | GnssSystem::Qzss | GnssSystem::Galileo | GnssSystem::BeiDou
    );
    if !kepler || v.len() < 28 {
        return;
    }
    let mut toc = gps_seconds(when);
    if sat.sys == GnssSystem::BeiDou {
        toc += BDT_GPST;
    }
    let week = v[21];
    let toes = v[11];
    let toe = if sat.sys == GnssSystem::BeiDou {
        (week + BDS_WEEK0 as f64) * WEEK + toes + BDT_GPST
    } else {
        // a week number written mod 1024 by old software: put toe next to toc
        let mut toe = week * WEEK + toes;
        while toe - toc > WEEK / 2.0 {
            toe -= WEEK;
        }
        while toc - toe > WEEK / 2.0 {
            toe += WEEK;
        }
        toe
    };
    let mut svh = v[24];
    if sat.sys == GnssSystem::Galileo && (svh as i64).trailing_zeros() >= 9 {
        svh = 0.0; // Galileo health bits (DVS / HS); 0 = usable
    }
    nav.kep.entry(sat).or_default().push(KeplerEph {
        toc,
        f0: v[0],
        f1: v[1],
        f2: v[2],
        iode: v[3],
        crs: v[4],
        deln: v[5],
        m0: v[6],
        cuc: v[7],
        e: v[8],
        cus: v[9],
        sqrta: v[10],
        toes,
        cic: v[12],
        omg0: v[13],
        cis: v[14],
        i0: v[15],
        crc: v[16],
        omg: v[17],
        omgd: v[18],
        idot: v[19],
        week,
        svh,
        tgd: v[25],
        toe,
        acc: v[23],
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    const RNX4: &str =
        "     4.01           N: GNSS NAV DATA    M: MIXED            RINEX VERSION / TYPE
                                                            END OF HEADER
> EPH G05 LNAV
G05 2024 09 23 00 00 00 1.234567890123E-04 1.000000000000E-12 0.000000000000E+00
     1.000000000000E+01 2.000000000000E+01 3.000000000000E-09 1.000000000000E+00
     1.000000000000E-06 5.000000000000E-03 2.000000000000E-06 5.153700000000E+03
     8.640000000000E+04 1.000000000000E-07 1.000000000000E+00 2.000000000000E-07
     9.000000000000E-01 2.500000000000E+02 1.000000000000E+00 -8.00000000000E-09
     1.000000000000E-10 1.000000000000E+00 2.333000000000E+03 0.000000000000E+00
     2.000000000000E+00 0.000000000000E+00-1.000000000000E-08 1.000000000000E+01
     8.000000000000E+04 4.000000000000E+00
> EPH G06 CNAV
G06 2024 09 23 00 00 00 1.0E-04 1.0E-12 0.0E+00
     1 2 3 4
     1 2 3 4
     1 2 3 4
     1 2 3 4
     1 2 3 4
     1 2 3 4
     1 2 3 4
     1 2 3 4
> STO G GPUT
    2024 09 23 00 00 00 GPUT
     1 2 3 4
> EPH R07 FDMA
R07 2024 09 23 00 15 00-1.000000000000E-05 0.000000000000E+00 5.400000000000E+04
     1.000000000000E+04 1.000000000000E+00 0.000000000000E+00 0.000000000000E+00
     2.000000000000E+04 2.000000000000E+00 0.000000000000E+00-4.000000000000E+00
     3.000000000000E+03 3.000000000000E+00 0.000000000000E+00 0.000000000000E+00
     0.000000000000E+00 0.000000000000E+00 0.000000000000E+00 0.000000000000E+00
";

    #[test]
    fn rinex4_dispatch() {
        let mut nav = NavData::default();
        read_nav_text(RNX4.as_bytes(), &mut nav).unwrap();
        let g: Sat = "G05".parse().unwrap();
        let k = nav.kep[&g][0];
        assert_eq!(k.f0, 1.234_567_890_123e-4);
        assert_eq!(k.sqrta, 5153.7);
        assert_eq!(k.week, 2333.0);
        assert_eq!(k.toe, 2333.0 * WEEK + 86_400.0);
        assert_eq!(k.toc, 1_411_084_800.0);
        assert_eq!(k.acc, 2.0);
        assert_eq!(k.tgd, -1e-8);
        assert_eq!(nav.skipped["G CNAV"], 1);
        assert_eq!(nav.skipped["STO"], 1);
        let r: Sat = "R07".parse().unwrap();
        let e = nav.glo[&r][0];
        assert_eq!(e.taun, 1e-5);
        assert_eq!(e.pos, [1e7, 2e7, 3e6]);
        assert_eq!(e.fcn, -4.0);
        assert_eq!(e.toe, 1_411_084_800.0 + 900.0 + 18.0);
        assert_eq!(nav.glo_fcn()[&r], -4);
    }

    #[test]
    fn unreadable_file_is_bad() {
        let dir = std::env::temp_dir().join("geomaticist_gnss_nav_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("bad.24n");
        std::fs::write(&p, b"xx.x     \n").unwrap();
        let nav = read_nav(&[&p, &dir.join("missing.24n")]);
        assert_eq!(nav.bad.len(), 1);
        assert!(nav.kep.is_empty());
    }
}
