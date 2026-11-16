//! RINEX 3.05 navigation writer, one file per system.
//!
//! Header sources (the LAST instance of each record in the file is used):
//! - 0x0f: GPS iono alpha/beta, UTC A0/A1/tot, dtLS/dtLSF, 8-bit WNt, WN_LSF, DN
//! - 0x20 sub 0x05: QZSS iono alpha/beta, UTC A0/A1/tot, WNt
//! - 0x20 sub 0x04: Galileo ai0..ai2, GST-UTC A0/A1/tot/WNt, GST-GPS A0G/A1G/t0G/WN0G
//! - 0x20 sub 0x06: BDS iono alpha/beta, BDT-UTC A0/A1 (+ a 14.0 written as T); fw 5.41: alpha/beta only
//! - 0x20 sub 0x02: GLONASS tau_c (RINEX GLUT A0 = -tau_c; T/W copied from GPS UTC)
//!
//! A missing record still gets a BDUT / GLUT line, with zeros.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use super::T04Error;
use super::container::{Inner, Records};
use super::nav::{self, NavRecord, fmt_d19, fmt_e};
use crate::util::{be, days_since_gps0, hline, ymdhms};

const ORDER: &str = "GREJCS";

fn sysname(s: char) -> &'static str {
    match s {
        'G' => "G: GPS",
        'R' => "R: GLONASS",
        'E' => "E: GALILEO",
        'C' => "C: BDS",
        'J' => "J: QZSS",
        _ => "S: SBAS Payload",
    }
}

fn ext(s: char) -> char {
    match s {
        'G' => 'n',
        'R' => 'g',
        'E' => 'l',
        'C' => 'c',
        'J' => 'q',
        _ => 'h',
    }
}

/// GPS week / day number of each leap-second event, by the new GPS-UTC value (public IERS history).
fn leap_event(dtlsf: i64) -> Option<(i64, i64)> {
    match dtlsf {
        18 => Some((1929, 7)),
        17 => Some((1851, 3)),
        16 => Some((1694, 7)),
        15 => Some((1512, 7)),
        14 => Some((1356, 7)),
        _ => None,
    }
}

/// Full GPS week from an 8-bit week, nearest to `ref_week`.
#[must_use]
pub fn expand_week(w8: i64, ref_week: i64) -> i64 {
    let base = ref_week - ref_week.rem_euclid(256) + w8;
    let mut best = base - 256;
    for w in [base, base + 256] {
        if (w - ref_week).abs() < (best - ref_week).abs() {
            best = w;
        }
    }
    best
}

/// The navigation files' text, by system letter (in "GREJCS" order), with their file names.
pub struct NavText {
    /// `(system, file name, text)`.
    pub files: Vec<(char, String, String)>,
}

fn header_sources(inner: &Inner) -> Result<HashMap<i32, Vec<u8>>, T04Error> {
    let mut last: HashMap<i32, Vec<u8>> = HashMap::new();
    for r in Records::new(inner) {
        let r = r?;
        let b = r.raw;
        if r.rtype == 0x0F && b.len() >= 101 {
            last.insert(-1, b.to_vec()); // "gps"
        } else if r.rtype == 0x20 && b.len() > 2 && [2, 4, 5, 6].contains(&b[2]) {
            last.insert(i32::from(b[2]), b.to_vec());
        }
    }
    Ok(last)
}

/// Builds the navigation files (without writing them).
#[allow(clippy::too_many_lines)]
pub fn nav_text(
    inner: &Inner,
    raw_name: &str,
    stem: &str,
    program: &str,
    run_by: &str,
    date: &str,
    line_end: &str,
) -> Result<NavText, T04Error> {
    let all = nav::iter_nav(inner)?;
    let mut order_seen: Vec<char> = Vec::new();
    let mut recs: HashMap<char, Vec<&NavRecord>> = HashMap::new();
    for r in &all {
        if !order_seen.contains(&r.sys) {
            order_seen.push(r.sys);
        }
        recs.entry(r.sys).or_default().push(r);
    }
    let src = header_sources(inner)?;
    let ref_week = order_seen
        .first()
        .and_then(|s| recs.get(s))
        .and_then(|v| v.first())
        .map(|r| days_since_gps0(&r.epoch).div_euclid(7));
    let gps = src.get(&-1);
    let mut gput: Option<(f64, f64, i64, i64)> = None;
    let mut leap: Option<String> = None;
    if let Some(g) = gps {
        let wnt = expand_week(i64::from(g[98]), ref_week.unwrap_or(0));
        let (dtls, dtlsf) = (be::u16(g, 90), be::u16(g, 92));
        let (wn_lsf, dn) =
            leap_event(dtlsf).unwrap_or((expand_week(i64::from(g[99]), wnt), i64::from(g[100])));
        #[allow(clippy::cast_possible_truncation)]
        let tot = be::f64(g, 82) as i64;
        gput = Some((be::f64(g, 74), be::f64(g, 66), tot, wnt));
        leap = Some(format!("{dtls:6}{dtlsf:6}{wn_lsf:6}{dn:6}"));
    }
    let mut files = Vec::new();
    for s in ORDER.chars() {
        let Some(list) = recs.get(&s) else { continue };
        let yy = list[0].epoch.year.rem_euclid(100);
        let mut l = vec![
            hline(
                &format!(
                    "{:>9}{:11}{:<20}{:<20}",
                    "3.05",
                    "",
                    "N: GNSS NAV DATA",
                    sysname(s)
                ),
                "RINEX VERSION / TYPE",
            ),
            hline(
                &format!("{program:<20.20}{run_by:<20.20}{date}"),
                "PGM / RUN BY / DATE",
            ),
            hline(
                &format!("Decoded from {raw_name} by an open clean-room decoder"),
                "COMMENT",
            ),
        ];
        let iono = |l: &mut Vec<String>, tag: &str, vals: &[Option<f64>]| {
            let body: String = vals
                .iter()
                .map(|v| v.map_or_else(|| " ".repeat(12), |v| fmt_e(v, 4, 12)))
                .collect();
            l.push(hline(&format!("{tag:<4} {body}"), "IONOSPHERIC CORR"));
        };
        let tcorr = |l: &mut Vec<String>, tag: &str, a0: f64, a1: f64, t: i64, w: i64| {
            l.push(hline(
                &format!(
                    "{tag:<4} {}{} {t:6} {w:4}",
                    fmt_e(a0, 10, 17),
                    fmt_e(a1, 9, 16)
                ),
                "TIME SYSTEM CORR",
            ));
        };
        let f4 = |b: &[u8], o: usize| -> Vec<Option<f64>> {
            (0..4).map(|k| Some(be::f64(b, o + 8 * k))).collect()
        };
        match s {
            'G' => {
                if let (Some(g), Some(gp)) = (gps, gput) {
                    iono(&mut l, "GPSA", &f4(g, 2));
                    iono(&mut l, "GPSB", &f4(g, 34));
                    tcorr(&mut l, "GPUT", gp.0, gp.1, gp.2, gp.3);
                    l.push(hline(leap.as_deref().unwrap_or(""), "LEAP SECONDS"));
                }
            }
            'R' => {
                if let Some(gp) = gput {
                    // no record (fw 5.41): tau_c 0
                    let a0 = src.get(&2).map_or(0.0, |g2| -be::f64(g2, 10));
                    tcorr(&mut l, "GLUT", a0, 0.0, gp.2, gp.3);
                }
            }
            'E' => {
                if let Some(b) = src.get(&4) {
                    iono(
                        &mut l,
                        "GAL",
                        &[
                            Some(be::f64(b, 12)),
                            Some(be::f64(b, 20)),
                            Some(be::f64(b, 28)),
                            None,
                        ],
                    );
                    tcorr(
                        &mut l,
                        "GAUT",
                        be::f64(b, 37),
                        be::f64(b, 45),
                        be::u32(b, 55),
                        be::u16(b, 59),
                    );
                    tcorr(
                        &mut l,
                        "GAGP",
                        be::f64(b, 67),
                        be::f64(b, 75),
                        be::u32(b, 83),
                        be::u16(b, 87),
                    );
                    if let Some(lp) = &leap {
                        l.push(hline(lp, "LEAP SECONDS"));
                    }
                }
            }
            'J' => {
                if let Some(b) = src.get(&5) {
                    iono(&mut l, "QZSA", &f4(b, 3));
                    iono(&mut l, "QZSB", &f4(b, 35));
                    #[allow(clippy::cast_possible_truncation)]
                    let t = be::f64(b, 83) as i64;
                    tcorr(
                        &mut l,
                        "QZUT",
                        be::f64(b, 75),
                        be::f64(b, 67),
                        t,
                        be::u16(b, 95),
                    );
                    if let Some(lp) = &leap {
                        l.push(hline(lp, "LEAP SECONDS"));
                    }
                }
            }
            'C' => {
                let bweek = ref_week.unwrap_or(0) - 1356;
                match src.get(&6) {
                    None => tcorr(&mut l, "BDUT", 0.0, 0.0, 14, bweek), // defaults
                    Some(b) => {
                        iono(&mut l, "BDSA", &f4(b, 14));
                        iono(&mut l, "BDSB", &f4(b, 46));
                        if b.len() >= 104 {
                            #[allow(clippy::cast_possible_truncation)]
                            let t = be::f64(b, 96) as i64;
                            tcorr(&mut l, "BDUT", be::f64(b, 80), be::f64(b, 88), t, bweek);
                        } else {
                            // fw 5.41's 78 bytes: iono only; zeros, GPS week
                            tcorr(&mut l, "BDUT", 0.0, 0.0, 0, ref_week.unwrap_or(0));
                        }
                    }
                }
                if let Some(g) = gps {
                    l.push(hline(
                        &format!("{:6}{:18}BDS", be::u16(g, 90) - 14, ""),
                        "LEAP SECONDS",
                    ));
                }
            }
            _ => {
                // 'S'
                if let (Some(g), Some(gp)) = (gps, gput) {
                    // positive A0 layout assumed (not met yet)
                    l.push(hline(
                        &format!(
                            "SBUT   {:>19}{}",
                            fmt_e(gp.0, 12, 19).trim(),
                            fmt_e(gp.1, 12, 19)
                        ),
                        "CORR TO SYSTEM TIME",
                    ));
                    l.push(hline(&format!("{:6}", be::u16(g, 90)), "LEAP SECONDS"));
                }
            }
        }
        l.push(hline("", "END OF HEADER"));
        for r in list {
            let mut v = r.vals.clone();
            if "GEJ".contains(s) {
                v.extend([Some(0.0), Some(0.0)]);
            } else if s == 'C' {
                v.extend([None, None]);
            }
            let mut first = format!("{} {}", r.sat, ymdhms(&r.epoch));
            for x in v.iter().take(3) {
                first.push_str(&fmt_d19(*x));
            }
            l.push(first);
            for chunk in v.get(3..).unwrap_or(&[]).chunks(4) {
                let mut line = String::from("    ");
                for x in chunk {
                    line.push_str(&fmt_d19(*x));
                }
                l.push(line);
            }
        }
        let mut text = l.join(line_end);
        text.push_str(line_end);
        files.push((s, format!("{stem}.{yy:02}{}", ext(s)), text));
    }
    Ok(NavText { files })
}

/// Writes the navigation files into `out_dir` -> `(system, path)` in "GREJCS" order.
pub fn write_nav(
    inner: &Inner,
    raw_name: &str,
    out_dir: &Path,
    stem: &str,
    program: &str,
    run_by: &str,
    date: &str,
) -> Result<Vec<(char, PathBuf)>, T04Error> {
    let t = nav_text(inner, raw_name, stem, program, run_by, date, "\r\n")?;
    let mut out = Vec::new();
    for (s, name, text) in t.files {
        let p = out_dir.join(name);
        std::fs::write(&p, crate::util::latin1_bytes(&text))
            .map_err(|e| T04Error::Io(p.display().to_string(), e.to_string()))?;
        out.push((s, p));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eight_bit_week() {
        assert_eq!(expand_week(0x2a, 2386), 2346);
        assert_eq!(expand_week(80, 2386), 2384);
        assert_eq!(expand_week(255, 2304), 2303);
    }
}
