//! A RINEX 3 observation file from a [`StationObs`]: the way out for observations that did not come from
//! a RINEX file (a decoded correction stream, simulated data).
//!
//! The file is plain RINEX 3.04: one `SYS / # / OBS TYPES` record per system with
//! the codes the station holds, epochs in time order, values in F14.3 with the loss-of-lock indicator on
//! phases, blank fields where a value is missing. [`crate::obs::read_obs_kinds`] reads it back to the same
//! values (to the three decimals of the format).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use geomaticist_gnss_core::time::{TICKS, header_time, v3_epoch_time};
use geomaticist_gnss_core::{ObsCode, Sat, StationObs};

/// Header facts a [`StationObs`] does not carry.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StationHeader {
    /// `PGM / RUN BY / DATE`: program name.
    pub program: String,
    /// `PGM / RUN BY / DATE`: who ran it.
    pub run_by: String,
    /// `PGM / RUN BY / DATE`: date text (`yyyymmdd hhmmss UTC`).
    pub date: String,
    /// `COMMENT` lines.
    pub comments: Vec<String>,
    /// GLONASS frequency channels for `GLONASS SLOT / FRQ #`.
    pub glonass_fcn: BTreeMap<Sat, i32>,
}

fn line(out: &mut String, content: &str, label: &str) {
    let _ = writeln!(out, "{content:<60.60}{label}");
}

fn ticks(t: f64) -> i64 {
    (t * TICKS as f64).round() as i64
}

/// The RINEX 3.04 text of `st` (LF line ends).
#[must_use]
pub fn station_to_rinex3(st: &StationObs, h: &StationHeader) -> String {
    // codes per system, in the order C L D S per band and attribute
    let mut types: BTreeMap<char, Vec<ObsCode>> = BTreeMap::new();
    for (sat, sigs) in &st.sats {
        let v = types.entry(sat.sys.letter()).or_default();
        for code in sigs.keys() {
            if !v.contains(code) {
                v.push(*code);
            }
        }
    }
    let kind_rank = |c: &ObsCode| "CLDS".find(c.kind_char()).unwrap_or(9);
    for v in types.values_mut() {
        v.sort_by(|a, b| {
            let (ta, tb) = (&a.as_str()[1..], &b.as_str()[1..]);
            ta.cmp(tb).then(kind_rank(a).cmp(&kind_rank(b)))
        });
    }
    let mut out = String::new();
    line(
        &mut out,
        "     3.04           OBSERVATION DATA    M",
        "RINEX VERSION / TYPE",
    );
    line(
        &mut out,
        &format!("{:<20.20}{:<20.20}{:<20.20}", h.program, h.run_by, h.date),
        "PGM / RUN BY / DATE",
    );
    for c in &h.comments {
        line(&mut out, c, "COMMENT");
    }
    line(&mut out, &st.name, "MARKER NAME");
    line(&mut out, "", "OBSERVER / AGENCY");
    line(
        &mut out,
        &format!("{:<20}{:<20.20}", "", st.receiver),
        "REC # / TYPE / VERS",
    );
    line(
        &mut out,
        &format!("{:<20}{:<20.20}", "", st.antenna),
        "ANT # / TYPE",
    );
    line(
        &mut out,
        &format!(
            "{:14.4}{:14.4}{:14.4}",
            st.approx[0], st.approx[1], st.approx[2]
        ),
        "APPROX POSITION XYZ",
    );
    line(
        &mut out,
        &format!(
            "{:14.4}{:14.4}{:14.4}",
            st.delta_hen[0], st.delta_hen[1], st.delta_hen[2]
        ),
        "ANTENNA: DELTA H/E/N",
    );
    let order = geomaticist_gnss_core::system::SYSTEM_ORDER;
    let mut systems: Vec<char> = types.keys().copied().collect();
    systems.sort_by_key(|s| order.find(*s).unwrap_or(99));
    for s in &systems {
        let v = &types[s];
        for (k, chunk) in v.chunks(13).enumerate() {
            let mut c = if k == 0 {
                format!("{s}  {:3}", v.len())
            } else {
                "      ".to_owned()
            };
            for code in chunk {
                let _ = write!(c, " {}", code.as_str());
            }
            line(&mut out, &c, "SYS / # / OBS TYPES");
        }
    }
    if let (Some(a), Some(b)) = (st.t.first(), st.t.last()) {
        line(
            &mut out,
            &header_time(ticks(*a), "GPS"),
            "TIME OF FIRST OBS",
        );
        line(&mut out, &header_time(ticks(*b), "GPS"), "TIME OF LAST OBS");
    }
    if !h.glonass_fcn.is_empty() {
        let all: Vec<(&Sat, &i32)> = h.glonass_fcn.iter().collect();
        for (k, chunk) in all.chunks(8).enumerate() {
            let mut c = if k == 0 {
                format!("{:3} ", all.len())
            } else {
                "    ".to_owned()
            };
            for (sat, f) in chunk {
                let _ = write!(c, "{sat} {f:2} ");
            }
            line(&mut out, &c, "GLONASS SLOT / FRQ #");
        }
    }
    line(&mut out, "", "END OF HEADER");
    let sats: Vec<Sat> = {
        let mut v: Vec<Sat> = st.sats.keys().copied().collect();
        v.sort_by_key(|s| s.rinex_order());
        v
    };
    for (i, &t) in st.t.iter().enumerate() {
        let present: Vec<Sat> = sats
            .iter()
            .copied()
            .filter(|s| st.sats[s].values().any(|g| g.values[i].is_finite()))
            .collect();
        let flag = u8::from(st.power_fail.get(i).copied().unwrap_or(false));
        let _ = writeln!(
            out,
            "> {}  {flag}{:3}",
            v3_epoch_time(ticks(t)),
            present.len()
        );
        for sat in present {
            let sigs = &st.sats[&sat];
            let mut rec = sat.to_string();
            for code in &types[&sat.sys.letter()] {
                match sigs.get(code).filter(|g| g.values[i].is_finite()) {
                    Some(g) => {
                        let lli = g.lli.as_ref().map_or(0, |l| l[i]);
                        let _ = write!(rec, "{:14.3}", g.values[i]);
                        if lli > 0 {
                            let _ = write!(rec, "{lli} ");
                        } else {
                            rec.push_str("  ");
                        }
                    }
                    None => rec.push_str("                "),
                }
            }
            out.push_str(rec.trim_end());
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::obs::read_obs_bytes_kinds;
    use geomaticist_gnss_core::obs::Signal;
    use geomaticist_gnss_core::signal::code_key;
    use std::path::Path;

    #[test]
    fn writes_what_the_reader_reads_back() {
        let mut st = StationObs::new("TEST", "t.obs");
        st.antenna = "TRM59800.00     NONE".into();
        st.approx = [4_395_470.6, 617_743.1, 4_565_607.2];
        st.delta_hen = [1.5, 0.0, 0.0];
        st.t = vec![1_391_385_600.0, 1_391_385_600.2, 1_391_385_600.4];
        st.power_fail = vec![false; 3];
        let mut order = 0;
        for (sat, codes) in [
            ("G05", vec!["C1C", "L1C", "D1C", "S1C", "C2W", "L2W"]),
            ("R07", vec!["C1C", "L1C"]),
        ] {
            let sat: Sat = sat.parse().unwrap();
            for code in codes {
                let key = code_key(code.as_bytes(), sat.sys.letter(), 3, 3.04);
                let base = if key.kind == 'L' { 1.1e8 } else { 2.1e7 };
                let mut values = vec![base + 0.125, base + 1.5, f64::NAN];
                if sat.sys.letter() == 'R' {
                    values[0] = f64::NAN;
                }
                st.sats.entry(sat).or_default().insert(
                    code.parse().unwrap(),
                    Signal {
                        key,
                        values,
                        lli: (key.kind == 'L').then(|| vec![0, 1, 0]),
                        first_seen: order,
                    },
                );
                order += 1;
            }
        }
        let mut h = StationHeader {
            program: "test".into(),
            ..StationHeader::default()
        };
        h.glonass_fcn.insert("R07".parse().unwrap(), -4);
        let text = station_to_rinex3(&st, &h);
        assert!(text.contains("G    6 C1C L1C D1C S1C C2W L2W"));
        assert!(text.contains("  1 R07 -4 "));
        assert!(
            text.lines()
                .all(|l| l.len() <= 80 || !l.contains("OBS TYPES"))
        );
        let back = read_obs_bytes_kinds(
            Path::new("t.obs"),
            text.as_bytes(),
            "",
            None,
            None,
            "GR",
            "CLSD",
        )
        .unwrap();
        // the last epoch has no data: it is written with no satellites and read as an epoch
        assert_eq!(back.t.len(), 3);
        assert_eq!(back.antenna, st.antenna);
        assert_eq!(back.delta_hen, st.delta_hen);
        for (sat, sigs) in &st.sats {
            for (code, g) in sigs {
                let b = &back.sats[sat][code];
                for i in 0..3 {
                    assert!(
                        (g.values[i].is_nan() && b.values[i].is_nan())
                            || (g.values[i] - b.values[i]).abs() < 6e-4,
                        "{sat} {code} {i}"
                    );
                }
                assert_eq!(b.lli, g.lli);
            }
        }
    }
}
