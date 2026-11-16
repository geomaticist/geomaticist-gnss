//! RINEX 2.11 observations from the unpacked .TGD records, written by these rules.
//!
//! - Types: C1 L1 L2 P2, plus C2 when any satellite has L2C, C5 L5 when any has L5; alphabetical.
//! - L2 is the L2C phase when a satellite has both L2C and L2 P(Y), else the P(Y) phase.
//! - LLI: 1 when the "tracked since the last epoch" flag (0x04) is clear, the slip flag (0x08) is set, or the
//!   phase was not written the epoch before; GPS L2 P(Y) phase and P2 always carry 4 (anti-spoofing).
//! - A satellite's first line is padded to one field short of full; other lines end at their last value.
//! - SSI: fixed thresholds per signal.
//! - PRN / # OF OBS counts L2 once per L2-band phase.

use std::collections::{HashMap, HashSet};

use geomaticist_gnss_core::time::DateTime;

use super::records::{self as rec, Epoch, Header, Sat, Signal};
use crate::decimal::Dec;
use crate::util::{gps_week_fsec, hline, llh_deg_to_xyz};

const SYS_ORDER: &str = "GR";

/// GPS - UTC by date (`LEAP_TABLE`).
#[must_use]
pub fn leap_seconds(t: &DateTime) -> i64 {
    let key = (t.year, t.month, t.day);
    for (start, n) in [
        ((2017, 1, 1), 18),
        ((2015, 7, 1), 17),
        ((2012, 7, 1), 16),
        ((2009, 1, 1), 15),
        ((2006, 1, 1), 14),
    ] {
        if key >= start {
            return n;
        }
    }
    13
}

/// Lowest SNR (0.1 dB-Hz) for each SSI by (system id, band, code); unseen steps are the gap's midpoint.
fn ssi_from(sysid: u8, band: u8, code: u8) -> &'static [(i64, i64)] {
    match (sysid, band, code) {
        (0, 0, 0) => &[
            (2, 339),
            (3, 357),
            (4, 374),
            (5, 392),
            (6, 425),
            (7, 476),
            (8, 525),
            (9, 575),
        ],
        (2, 0, 0) => &[
            (2, 339),
            (3, 357),
            (4, 374),
            (5, 392),
            (6, 419),
            (7, 457),
            (8, 494),
            (9, 532),
        ],
        (0, 1, 2) => &[(5, 0), (6, 188), (7, 263), (8, 338), (9, 413)],
        (2, 1, 1) => &[(6, 0), (7, 262), (8, 338), (9, 413)],
        (0, 1, 5) => &[(8, 0), (9, 413)],
        (0, 2, 8) => &[(8, 0), (9, 444)],
        _ => &[],
    }
}

/// SSI of a signal.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn ssi(sysid: u8, sig: &Signal) -> i64 {
    let x = (sig.snr * 10.0).round_ties_even() as i64;
    let mut out = 1;
    for &(k, lo) in ssi_from(sysid, sig.band, sig.code) {
        if x >= lo {
            out = k;
        }
    }
    out
}

/// F14.3 + LLI + SSI: the exact binary value, ties away from zero.
#[must_use]
pub fn fmt(v: Option<(f64, i64, Option<i64>)>) -> String {
    let Some((v, lli, s)) = v else {
        return " ".repeat(16);
    };
    let d = Dec::exact(v).round_at(-3);
    let sign = if v < 0.0 && !d.is_zero() { "-" } else { "" };
    let txt = format!("{:>14}", format!("{sign}{}", d.abs_fixed()));
    if txt.len() > 14 {
        return " ".repeat(16);
    }
    let l = if lli == 0 {
        " ".to_owned()
    } else {
        lli.to_string()
    };
    let q = s.map_or_else(|| " ".to_owned(), |s| s.to_string());
    format!("{txt}{l}{q}")
}

type ObsMap = HashMap<&'static str, (f64, i64, Option<i64>)>;

/// `{type: (value, lli, ssi)}` for one satellite in one epoch, and the signals whose phase was written.
fn sat_obs(sat: &Sat, written_before: Option<&HashSet<(u8, u8)>>) -> (ObsMap, HashSet<(u8, u8)>) {
    let sysid = sat.sysid;
    // dict semantics: first position, last value
    let mut keys: Vec<(u8, u8)> = Vec::new();
    let mut by: HashMap<(u8, u8), &Signal> = HashMap::new();
    for s in &sat.signals {
        if rec::signal(sysid, s.band, s.code).is_none() {
            continue;
        }
        let k = (s.band, s.code);
        if by.insert(k, s).is_none() {
            keys.push(k);
        }
    }
    let mut out = ObsMap::new();
    let mut phases = HashSet::new();
    let l2c = by.get(&(1, 5)).copied();
    for key in keys {
        let s = by[&key];
        let Some((c_t, l_t)) = rec::signal(sysid, key.0, key.1) else {
            continue;
        };
        let q = ssi(sysid, s);
        let as_flag = if sat.sys == 'G' && key == (1, 2) {
            4
        } else {
            0
        };
        out.insert(c_t, (s.range, as_flag, Some(q)));
        if s.flags & 0x20 == 0 {
            continue;
        }
        phases.insert(key);
        if l_t == "L2" && key != (1, 5) && l2c.is_some_and(|c| c.flags & 0x20 != 0) {
            continue; // L2 = the L2C phase when there is one
        }
        let new = (s.flags & 0x0C) != 0x04 || !written_before.is_some_and(|w| w.contains(&key));
        out.insert(l_t, (s.phase, as_flag | i64::from(new), Some(q)));
    }
    (out, phases)
}

/// `GPS0 + timedelta(weeks=week, seconds=tow)`.
#[must_use]
pub fn epoch_time(ep: &Epoch) -> DateTime {
    gps_week_fsec(ep.week, ep.tow)
}

fn sat_order(sid: &str) -> (usize, String) {
    let c = sid.chars().next().unwrap_or('?');
    (SYS_ORDER.find(c).unwrap_or(9), sid.to_owned())
}

/// `%d-%b-%y %H:%M UTC`.
#[must_use]
pub fn pgm_date(t: &DateTime) -> String {
    const MON: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let m = usize::try_from(t.month - 1).unwrap_or(0).min(11);
    format!(
        "{:02}-{}-{:02} {:02}:{:02} UTC",
        t.day,
        MON[m],
        t.year.rem_euclid(100),
        t.hour,
        t.minute
    )
}

/// Options of [`write_obs`].
#[derive(Debug, Clone)]
pub struct ObsOptions<'a> {
    /// PGM.
    pub program: &'a str,
    /// RUN BY.
    pub run_by: &'a str,
    /// The PGM date text (`%d-%b-%y %H:%M UTC`).
    pub date: &'a str,
    /// COMMENT line, when not empty.
    pub comment: &'a str,
    /// (IGS name, dome).
    pub antenna: (&'a str, &'a str),
    /// Line end.
    pub line_end: &'a str,
}

fn sec_of(t: &DateTime) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let v = t.second as f64 + t.micro as f64 / 1e6;
    v
}

/// The RINEX 2.11 observation file text of a .TGD; returns (header text, whole file bytes).
pub fn write_obs(
    tgd: &[u8],
    hdr: &Header,
    opt: &ObsOptions<'_>,
) -> Result<(String, Vec<u8>), String> {
    let eps = rec::epochs(tgd)?;
    let mut have: HashSet<&'static str> = HashSet::new();
    let mut counts: HashMap<String, HashMap<&'static str, i64>> = HashMap::new();
    let mut prev: HashMap<String, HashSet<(u8, u8)>> = HashMap::new();
    let mut bodies: Vec<(&Epoch, Vec<(String, ObsMap)>)> = Vec::new();
    for ep in &eps {
        let mut lines = Vec::new();
        let mut cur = HashMap::new();
        let mut sats: Vec<&Sat> = ep.sats.iter().collect();
        sats.sort_by_cached_key(|s| sat_order(&s.id()));
        for sat in sats {
            let id = sat.id();
            let (o, ph) = sat_obs(sat, prev.get(&id));
            let n_l2 = ph.iter().filter(|k| k.0 == 1).count();
            cur.insert(id.clone(), ph);
            if o.is_empty() {
                continue;
            }
            let c = counts.entry(id.clone()).or_default();
            for t in o.keys() {
                have.insert(t);
                *c.entry(t).or_insert(0) += 1;
            }
            let extra = i64::try_from(n_l2).unwrap_or(0) - i64::from(o.contains_key("L2"));
            *c.entry("L2").or_insert(0) += extra.max(0);
            lines.push((id, o));
        }
        prev = cur;
        bodies.push((ep, lines));
    }
    if bodies.is_empty() {
        return Err("no measurement epochs in the file".into());
    }
    for t in ["C1", "L1", "L2", "P2"] {
        have.insert(t);
    }
    let mut types: Vec<&'static str> = have.into_iter().collect();
    types.sort_unstable();

    let first = epoch_time(bodies[0].0);
    let last = epoch_time(bodies[bodies.len() - 1].0);
    let mut l = vec![
        hline(
            &format!(
                "{:>9}{:11}{:<20}{:<20}",
                "2.11", "", "OBSERVATION DATA", "Mixed(MIXED)"
            ),
            "RINEX VERSION / TYPE",
        ),
        hline(
            &format!("{:<20.20}{:<20.20}{}", opt.program, opt.run_by, opt.date),
            "PGM / RUN BY / DATE",
        ),
    ];
    if !opt.comment.is_empty() {
        l.push(hline(opt.comment, "COMMENT"));
    }
    l.push(hline(&hdr.marker, "MARKER NAME"));
    l.push(hline("", "MARKER NUMBER"));
    l.push(hline("", "OBSERVER / AGENCY"));
    l.push(hline(
        &format!(
            "{:<20.20}{:<20.20}{:<20.20}",
            hdr.serial, hdr.model, hdr.firmware
        ),
        "REC # / TYPE / VERS",
    ));
    l.push(hline(
        &format!("{:<20}{:<16.16}{:<4.4}", "", opt.antenna.0, opt.antenna.1),
        "ANT # / TYPE",
    ));
    let (x, y, z) = hdr
        .approx_llh
        .map_or((0.0, 0.0, 0.0), |(a, b, c)| llh_deg_to_xyz(a, b, c));
    l.push(hline(
        &format!("{x:14.4}{y:14.4}{z:14.4}"),
        "APPROX POSITION XYZ",
    ));
    let h = hdr.arp_height.filter(|h| *h != 0.0).unwrap_or(0.0);
    l.push(hline(
        &format!("{h:14.4}{:14.4}{:14.4}", 0.0, 0.0),
        "ANTENNA: DELTA H/E/N",
    ));
    l.push(hline(
        &format!("{:6}{:6}{:6}", 1, 1, 0),
        "WAVELENGTH FACT L1/2",
    ));
    for (k, chunk) in types.chunks(9).enumerate() {
        let mut s = if k == 0 {
            format!("{:6}", types.len())
        } else {
            " ".repeat(6)
        };
        for t in chunk {
            s.push_str(&format!("{t:>6}"));
        }
        l.push(hline(&s, "# / TYPES OF OBSERV"));
    }
    if let Some(iv) = hdr.interval_s.filter(|v| *v != 0.0) {
        l.push(hline(&format!("{iv:10.3}"), "INTERVAL"));
    }
    for (t, lab) in [(&first, "TIME OF FIRST OBS"), (&last, "TIME OF LAST OBS")] {
        l.push(hline(
            &format!(
                "{:6}{:6}{:6}{:6}{:6}{:13.7}     GPS",
                t.year,
                t.month,
                t.day,
                t.hour,
                t.minute,
                sec_of(t)
            ),
            lab,
        ));
    }
    l.push(hline(&format!("{:6}", 0), "RCV CLOCK OFFS APPL"));
    l.push(hline(
        &format!("{:6}", leap_seconds(&first)),
        "LEAP SECONDS",
    ));
    let mut sats: Vec<&String> = counts.keys().collect();
    sats.sort_by_cached_key(|s| sat_order(s));
    l.push(hline(&format!("{:6}", sats.len()), "# OF SATELLITES"));
    for s in sats {
        let vals: Vec<i64> = types
            .iter()
            .map(|t| counts[s].get(t).copied().unwrap_or(0))
            .collect();
        for (k, chunk) in vals.chunks(9).enumerate() {
            let mut line = if k == 0 {
                format!("   {s}")
            } else {
                " ".repeat(6)
            };
            for v in chunk {
                line.push_str(&format!("{v:6}"));
            }
            l.push(hline(&line, "PRN / # OF OBS"));
        }
    }
    l.push(hline("", "END OF HEADER"));
    let le = opt.line_end;
    let mut head = l.join(le);
    head.push_str(le);

    let mut f = String::with_capacity(head.len() + bodies.len() * 400);
    f.push_str(&head);
    let pad_first = 16 * (types.len().min(5) - 1);
    for (ep, lines) in &bodies {
        let t = epoch_time(ep);
        let ids: String = lines.iter().map(|(s, _)| s.as_str()).collect();
        f.push_str(&format!(
            " {:2} {:2} {:2} {:2} {:2}{:11.7}  0{:3}{}{le}",
            t.year.rem_euclid(100),
            t.month,
            t.day,
            t.hour,
            t.minute,
            sec_of(&t),
            lines.len(),
            &ids[..ids.len().min(36)]
        ));
        let mut k = 36;
        while k < ids.len() {
            f.push_str(&" ".repeat(32));
            f.push_str(&ids[k..ids.len().min(k + 36)]);
            f.push_str(le);
            k += 36;
        }
        for (_sid, o) in lines {
            let flds: Vec<String> = types.iter().map(|t| fmt(o.get(t).copied())).collect();
            for (k, chunk) in flds.chunks(5).enumerate() {
                let joined = chunk.concat();
                let mut txt = joined.trim_end_matches(crate::util::py_space).to_owned();
                if k == 0 {
                    // a short first line is padded to one field short
                    while txt.len() < pad_first {
                        txt.push(' ');
                    }
                } else if txt.is_empty() {
                    // ... and writes an empty one as one blank field
                    txt = " ".repeat(16);
                }
                f.push_str(&txt);
                f.push_str(le);
            }
        }
    }
    Ok((head, crate::util::latin1_bytes(&f)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_format() {
        assert_eq!(fmt(None), " ".repeat(16));
        assert_eq!(fmt(Some((2.0625, 1, Some(7)))), "         2.06317");
        assert_eq!(fmt(Some((-0.0004, 0, None))), "         0.000  ");
        assert_eq!(fmt(Some((-2.5, 4, Some(5)))), "        -2.50045");
        assert_eq!(fmt(Some((1e12, 0, None))), " ".repeat(16));
    }

    #[test]
    fn dates() {
        let t = DateTime::from_micros_since_gps0(0);
        assert_eq!(pgm_date(&t), "06-Jan-80 00:00 UTC");
        assert_eq!(leap_seconds(&t), 13);
    }
}
