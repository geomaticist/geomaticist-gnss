//! Events and occupations in RINEX observation files, a check of the epoch times, and a writer that
//! inserts event records.
//!
//! RINEX marks things that are not observations with an *epoch flag* on the epoch line: 2 "start moving
//! antenna", 3 "new site occupation" (followed by header records such as `MARKER NAME`), 5 "external
//! event" (the epoch is the time of the event). This module reads flags 5 as [`Event`]s, flags 2 and 3 as
//! stop-and-go [`Occupation`]s, and writes flag-5 records for a list of events into a copy of a file.
//!
//! The reader is forgiving about how an event line was written. A line such as
//! `> 2024 09 23 18 04 05.123456 5  0`, written by scripts with six decimals and so with the flag one
//! column to the left of where the format puts it, is read as the event it is meant to be. The writer
//! writes the standard form: `> 2024 09 23 18 04  5.1234560  5  0`.

use std::path::Path;

use geomaticist_gnss_core::event::{Event, Occupation};
use geomaticist_gnss_core::text::{latin1, py_strip, sl};
use geomaticist_gnss_core::time::{
    TICKS, fmt, parse_seconds, py_int, ticks_to_seconds, to_ticks, v2_epoch_time,
};
use geomaticist_gnss_rinex::header::RinexHeader;
use geomaticist_gnss_rinex::open::read_plain;
use geomaticist_gnss_rinex::reader::read_header;
use geomaticist_gnss_rinex::writer::epoch_line_v3;

use crate::EventError;

fn err(e: impl std::fmt::Display) -> EventError {
    EventError(e.to_string())
}

/// One epoch line of a file with what follows it.
#[derive(Debug, Clone, PartialEq)]
struct Mark {
    /// Epoch time in ticks (0.1 µs since 1980-01-06), `None` when the line has no time.
    tick: Option<i64>,
    flag: u8,
    /// Byte offset of the epoch line in the file.
    at: usize,
    /// The records after an event line (flags 2 to 5), latin-1.
    records: Vec<String>,
}

/// The lines of `data` with their byte offsets (without the line ending).
fn lines(data: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    let mut pos = 0usize;
    std::iter::from_fn(move || {
        if pos >= data.len() {
            return None;
        }
        let start = pos;
        let end = data[pos..]
            .iter()
            .position(|&b| b == b'\n')
            .map_or(data.len(), |p| pos + p);
        pos = end + 1;
        let mut line = &data[start..end];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        Some((start, line))
    })
}

fn ticks_of(y: i64, mo: i64, d: i64, h: i64, mi: i64, sec: &[u8]) -> Option<i64> {
    to_ticks(y, mo, d, h, mi, parse_seconds(sec)?)
}

/// A RINEX 3 / 4 epoch line by its words, so that a flag written a column off is still found.
fn v3_epoch(line: &[u8]) -> Option<(Option<i64>, u8, usize)> {
    let text = latin1(line);
    let w: Vec<&str> = text[1..].split_whitespace().collect();
    match w.len() {
        // an event without a time: flag and count only
        2 => Some((None, w[0].parse().ok()?, w[1].parse().ok()?)),
        n if n >= 8 => {
            let i = |k: usize| w[k].parse::<i64>().ok();
            let tick = ticks_of(i(0)?, i(1)?, i(2)?, i(3)?, i(4)?, w[5].as_bytes());
            Some((tick, w[6].parse().ok()?, w[7].parse().ok()?))
        }
        _ => None,
    }
}

fn year2(yy: i64) -> i64 {
    if yy < 80 { 2000 + yy } else { 1900 + yy }
}

/// The epoch lines of a file, in file order.
fn marks(data: &[u8], hdr: &RinexHeader) -> Vec<Mark> {
    let mut out: Vec<Mark> = Vec::new();
    let mut it = lines(data).skip_while(|(_, l)| py_strip(sl(l, 60, 80)) != b"END OF HEADER");
    it.next();
    if hdr.major() >= 3 {
        for (at, line) in it {
            if line.first() == Some(&b'>') {
                if let Some((tick, flag, _)) = v3_epoch(line) {
                    out.push(Mark {
                        tick,
                        flag,
                        at,
                        records: Vec::new(),
                    });
                }
            } else if let Some(m) = out.last_mut()
                && (2..=5).contains(&m.flag)
            {
                m.records.push(latin1(line));
            }
        }
        return out;
    }
    // RINEX 2: the structure has to be followed line by line
    let per_sat = hdr.v2_types.len().div_ceil(5).max(1);
    while let Some((at, line)) = it.next() {
        if py_strip(line).is_empty() {
            continue;
        }
        let flag = py_int(sl(line, 28, 29)).unwrap_or(0);
        let n = py_int(sl(line, 29, 32)).unwrap_or(0).max(0) as usize;
        let flag = u8::try_from(flag).unwrap_or(u8::MAX);
        let tick = (|| {
            ticks_of(
                year2(py_int(sl(line, 1, 3))?),
                py_int(sl(line, 4, 6))?,
                py_int(sl(line, 7, 9))?,
                py_int(sl(line, 10, 12))?,
                py_int(sl(line, 13, 15))?,
                sl(line, 15, 26),
            )
        })();
        let mut m = Mark {
            tick,
            flag,
            at,
            records: Vec::new(),
        };
        if matches!(flag, 0 | 1 | 6) {
            let skip = n.saturating_sub(1) / 12 + n * per_sat;
            for _ in 0..skip {
                if it.next().is_none() {
                    break;
                }
            }
        } else {
            for _ in 0..n {
                match it.next() {
                    Some((_, l)) => m.records.push(latin1(l)),
                    None => break,
                }
            }
        }
        out.push(m);
    }
    out
}

fn open(data: &[u8]) -> Result<RinexHeader, EventError> {
    let (hdr, _) = read_header(data).map_err(err)?;
    Ok(hdr)
}

/// The events of a file with what the reader noticed.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RinexEvents {
    /// The flag-5 records with a time, in file order; labels are their running numbers from 1.
    pub events: Vec<Event>,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// The external events (epoch flag 5) of RINEX observation bytes.
pub fn events_from_bytes(data: &[u8]) -> Result<RinexEvents, EventError> {
    let hdr = open(data)?;
    let mut out = RinexEvents::default();
    let mut no_time = 0usize;
    for m in marks(data, &hdr).iter().filter(|m| m.flag == 5) {
        match m.tick {
            Some(t) => {
                let k = out.events.len() + 1;
                out.events.push(Event {
                    t: ticks_to_seconds(t),
                    label: k.to_string(),
                    offset: None,
                });
            }
            None => no_time += 1,
        }
    }
    if no_time > 0 {
        out.warnings.push(format!(
            "{no_time} event records have no time and were left out"
        ));
    }
    let ts = hdr.time_system();
    if ts != "GPS" && !out.events.is_empty() {
        out.warnings.push(format!(
            "the file's time system is {ts}, not GPS: the event times are taken as written"
        ));
    }
    Ok(out)
}

/// The external events (epoch flag 5) of a RINEX observation file (plain, gzip or Hatanaka).
pub fn read_events(path: &Path) -> Result<RinexEvents, EventError> {
    events_from_bytes(&read_plain(path).map_err(err)?)
}

fn marker_name(records: &[String]) -> Option<String> {
    records
        .iter()
        .find(|r| r.get(60..).is_some_and(|l| l.trim() == "MARKER NAME"))
        .map(|r| r[..60.min(r.len())].trim().to_owned())
        .filter(|s| !s.is_empty())
}

/// The stop-and-go occupations of RINEX observation bytes: from each "new site occupation" record (flag
/// 3) to the last epoch before the antenna moves again (flag 2) or before the next occupation. A file
/// that starts static and then moves gives a first occupation named by the header's `MARKER NAME`.
pub fn occupations_from_bytes(data: &[u8]) -> Result<Vec<Occupation>, EventError> {
    let hdr = open(data)?;
    let ms = marks(data, &hdr);
    if !ms.iter().any(|m| m.flag == 2 || m.flag == 3) {
        return Ok(Vec::new());
    }
    let mut out: Vec<Occupation> = Vec::new();
    // an occupation being collected: (name, first epoch, last epoch)
    let mut cur: Option<(String, Option<i64>, Option<i64>)> = None;
    let first_is_move = ms
        .iter()
        .find(|m| m.flag == 2 || m.flag == 3)
        .is_some_and(|m| m.flag == 2);
    if first_is_move {
        cur = Some((hdr.marker_name(), None, None));
    }
    let mut n_unnamed = 0usize;
    let close = |cur: &mut Option<(String, Option<i64>, Option<i64>)>,
                 out: &mut Vec<Occupation>| {
        if let Some((name, Some(a), Some(b))) = cur.take() {
            out.push(Occupation {
                name,
                t0: ticks_to_seconds(a),
                t1: ticks_to_seconds(b),
            });
        }
    };
    for m in &ms {
        match m.flag {
            0 | 1 => {
                if let (Some(c), Some(t)) = (cur.as_mut(), m.tick) {
                    c.1.get_or_insert(t);
                    c.2 = Some(t);
                }
            }
            2 => close(&mut cur, &mut out),
            3 => {
                close(&mut cur, &mut out);
                let name = marker_name(&m.records).unwrap_or_else(|| {
                    n_unnamed += 1;
                    format!("STOP-{n_unnamed}")
                });
                cur = Some((name, None, None));
            }
            _ => {}
        }
    }
    close(&mut cur, &mut out);
    Ok(out)
}

/// The stop-and-go occupations of a RINEX observation file.
pub fn read_occupations(path: &Path) -> Result<Vec<Occupation>, EventError> {
    occupations_from_bytes(&read_plain(path).map_err(err)?)
}

/// What [`check_times`] found.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TimeCheck {
    /// Observation epochs in the file.
    pub epochs: usize,
    /// Epochs dated before 1990 (a recorder that had no date yet writes 1980).
    pub before_1990: usize,
    /// Epochs displaced by whole days from their neighbours (a recorder's date glitch).
    pub displaced: usize,
    /// Epochs not later than the epoch before them (repeated or out of order).
    pub not_later: usize,
    /// First and last epoch in ticks, glitched epochs left out.
    pub span: Option<(i64, i64)>,
    /// One sentence per finding.
    pub warnings: Vec<String>,
}

/// Looks at the epoch times of an observation file for what recorders are known to get wrong:
///
/// - a first epoch (or more) dated 1980, written before the receiver knew the date;
/// - a run of epochs whose date is some days off while the time of day continues (seen shortly before
///   a GPS week change);
/// - repeated epochs;
/// - a `TIME OF FIRST OBS` in the header that is not the first epoch.
///
/// Nothing is changed: the findings are returned as sentences for the user.
pub fn check_times(data: &[u8]) -> Result<TimeCheck, EventError> {
    let hdr = open(data)?;
    let ticks: Vec<i64> = marks(data, &hdr)
        .iter()
        .filter(|m| m.flag <= 1)
        .filter_map(|m| m.tick)
        .collect();
    let mut c = TimeCheck {
        epochs: ticks.len(),
        ..TimeCheck::default()
    };
    let y1990 = to_ticks(1990, 1, 1, 0, 0, 0).unwrap_or(0);
    let day = 86_400 * TICKS;
    // the good timeline: epochs after 1990 that follow each other
    let mut good: Vec<i64> = Vec::with_capacity(ticks.len());
    let mut run: Vec<i64> = Vec::new();
    let mut days_off: Vec<i64> = Vec::new();
    for &t in &ticks {
        if t < y1990 {
            c.before_1990 += 1;
            continue;
        }
        let Some(&last) = good.last() else {
            good.push(t);
            continue;
        };
        let d = t - last;
        // off by a whole number of days (within 10 minutes), forward or back
        let k = (d as f64 / day as f64).round() as i64;
        if k != 0 && (d - k * day).abs() < 600 * TICKS {
            run.push(t);
            days_off.push(k);
            // a long run is not a glitch: the file really continues days later
            if run.len() > 2000 {
                good.append(&mut run);
                days_off.clear();
            }
            continue;
        }
        if !run.is_empty() {
            c.displaced += run.len();
            run.clear();
        }
        if d <= 0 {
            c.not_later += 1;
            continue;
        }
        good.push(t);
    }
    // a run that reaches the end of the file is the file's real time, not a glitch
    if !run.is_empty() {
        good.append(&mut run);
        days_off.clear();
    }
    c.span = good.first().copied().zip(good.last().copied());
    if c.before_1990 > 0 {
        c.warnings.push(format!(
            "{} epoch{} dated before 1990 (the recorder had no date yet): not usable as written",
            c.before_1990,
            if c.before_1990 == 1 { " is" } else { "s are" }
        ));
    }
    if c.displaced > 0 {
        days_off.sort_unstable();
        days_off.dedup();
        let days: Vec<String> = days_off.iter().map(|d| format!("{d:+}")).collect();
        c.warnings.push(format!(
            "{} epochs carry a date {} day(s) off their neighbours while the time of day continues (a \
             recorder glitch, seen near GPS week changes): correct the dates before processing",
            c.displaced,
            days.join(", ")
        ));
    }
    if c.not_later > 0 {
        c.warnings.push(format!(
            "{} epoch{} not later than the epoch before (repeated or out of order)",
            c.not_later,
            if c.not_later == 1 { " is" } else { "s are" }
        ));
    }
    if let (Some((first, _)), Some(h)) = (c.span, hdr.get("TIME OF FIRST OBS")) {
        let head = (|| {
            ticks_of(
                py_int(sl(h, 0, 6))?,
                py_int(sl(h, 6, 12))?,
                py_int(sl(h, 12, 18))?,
                py_int(sl(h, 18, 24))?,
                py_int(sl(h, 24, 30))?,
                sl(h, 30, 43),
            )
        })();
        if let Some(head) = head
            && (head - first).abs() > TICKS
        {
            c.warnings.push(format!(
                "the header's TIME OF FIRST OBS ({}) is not the first epoch ({})",
                fmt(head, 1),
                fmt(first, 1)
            ));
        }
    }
    Ok(c)
}

/// The result of [`insert_events`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Inserted {
    /// The new file's bytes.
    pub data: Vec<u8>,
    /// Events written.
    pub inserted: usize,
    /// Events before the first epoch (not written).
    pub before: usize,
    /// Events after the last epoch (not written).
    pub after: usize,
    /// Event records the file already had (kept).
    pub existing: usize,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// A copy of RINEX observation bytes with a flag-5 event record (`> date time  5  0`, no records
/// following) inserted for each of `events`, each before the first observation epoch that is later than
/// it. Events outside the file's epochs are not written (they could not be positioned anyway) and are
/// counted in the result. Everything else is copied byte for byte.
pub fn insert_events(data: &[u8], events: &[Event]) -> Result<Inserted, EventError> {
    let hdr = open(data)?;
    let v3 = hdr.major() >= 3;
    let ms = marks(data, &hdr);
    let obs: Vec<&Mark> = ms
        .iter()
        .filter(|m| m.flag <= 1 && m.tick.is_some())
        .collect();
    let (Some(first), Some(last)) = (
        obs.first().and_then(|m| m.tick),
        obs.last().and_then(|m| m.tick),
    ) else {
        return Err(EventError("the file has no observation epochs".into()));
    };
    let mut out = Inserted {
        existing: ms.iter().filter(|m| m.flag == 5).count(),
        ..Inserted::default()
    };
    let ts = hdr.time_system();
    if ts != "GPS" {
        return Err(EventError(format!(
            "the file's time system is {ts}; events are GPS times and are only written into GPS-time files"
        )));
    }
    if out.existing > 0 {
        out.warnings.push(format!(
            "the file already has {} event records: they are kept, the new ones are added",
            out.existing
        ));
    }
    let eol: &[u8] = if data.windows(2).any(|w| w == b"\r\n") {
        b"\r\n"
    } else {
        b"\n"
    };
    let mut ev: Vec<i64> = events
        .iter()
        .filter(|e| e.t.is_finite())
        .map(|e| (e.t * TICKS as f64).round() as i64)
        .collect();
    ev.sort_unstable();
    out.before = ev.iter().filter(|t| **t <= first).count();
    out.after = ev.iter().filter(|t| **t > last).count();
    let ev: Vec<i64> = ev
        .into_iter()
        .filter(|t| *t > first && *t <= last)
        .collect();
    let line_of = |t: i64| -> Vec<u8> {
        if v3 {
            epoch_line_v3(Some(t), 5, 0, b"")
        } else {
            format!("{}  5  0", v2_epoch_time(t)).into_bytes()
        }
    };
    out.data = Vec::with_capacity(data.len() + ev.len() * 40);
    let mut copied = 0usize;
    let mut k = 0usize;
    for m in &obs {
        let t = m.tick.unwrap_or(first);
        if k < ev.len() && ev[k] <= t {
            out.data.extend_from_slice(&data[copied..m.at]);
            copied = m.at;
            while k < ev.len() && ev[k] <= t {
                out.data.extend(line_of(ev[k]));
                out.data.extend_from_slice(eol);
                k += 1;
                out.inserted += 1;
            }
        }
    }
    out.data.extend_from_slice(&data[copied..]);
    if out.before + out.after > 0 {
        out.warnings.push(format!(
            "{} events lie outside the file's epochs ({} before the first, {} after the last) and were \
             not written",
            out.before + out.after,
            out.before,
            out.after
        ));
    }
    Ok(out)
}

/// The result of [`mark_kinematic`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Kinematic {
    /// The new file's bytes.
    pub data: Vec<u8>,
    /// A "start moving antenna" record (flag 2) was put before the first epoch (false: one was there).
    pub flag_added: bool,
    /// Blank or whitespace-only lines dropped from the data section (RINEX 3 / 4 only).
    pub blank_removed: usize,
}

/// A copy of RINEX observation bytes marked as a moving receiver's, so that office software that reads
/// the file's own words (Trimble Business Center among them) takes it as kinematic without asking:
/// - RINEX 3 / 4: the header's `MARKER TYPE` is set to `marker_type` (e.g. `AIRBORNE`), the record added
///   before `END OF HEADER` when there is none;
/// - a "start moving antenna" record (epoch flag 2, no date, no records) goes just before the first
///   observation epoch, unless the file already has one before it;
/// - RINEX 3 / 4: blank and whitespace-only lines in the data section are dropped (DJI writes a long
///   line of spaces after `END OF HEADER`). In RINEX 2 a blank line can be a satellite's empty
///   observations, so there they are kept.
///
/// Only for a receiver that moves the whole time (a drone, a vehicle); never for a base.
pub fn mark_kinematic(data: &[u8], marker_type: &str) -> Result<Kinematic, EventError> {
    let hdr = open(data)?;
    let v3 = hdr.major() >= 3;
    let ms = marks(data, &hdr);
    let first = ms
        .iter()
        .position(|m| m.flag <= 1 && m.tick.is_some())
        .ok_or_else(|| EventError("the file has no observation epochs".into()))?;
    let need_flag = !ms[..first].iter().any(|m| m.flag == 2);
    let first_at = ms[first].at;
    let eol: &[u8] = if data.windows(2).any(|w| w == b"\r\n") {
        b"\r\n"
    } else {
        b"\n"
    };
    let mt_line = format!("{:<60}{:<20}", marker_type.trim(), "MARKER TYPE");
    let flag_line = if v3 {
        epoch_line_v3(None, 2, 0, b"")
    } else {
        format!("{}2  0", " ".repeat(28)).into_bytes()
    };
    let mut out = Kinematic {
        data: Vec::with_capacity(data.len() + 200),
        ..Kinematic::default()
    };
    let mut in_header = true;
    let mut have_mt = false;
    for (at, line) in lines(data) {
        if in_header {
            let label = py_strip(sl(line, 60, 80));
            if v3 && label == b"MARKER TYPE" {
                out.data.extend_from_slice(mt_line.as_bytes());
                out.data.extend_from_slice(eol);
                have_mt = true;
                continue;
            }
            if label == b"END OF HEADER" {
                if v3 && !have_mt {
                    out.data.extend_from_slice(mt_line.as_bytes());
                    out.data.extend_from_slice(eol);
                }
                in_header = false;
            }
            out.data.extend_from_slice(line);
            out.data.extend_from_slice(eol);
            continue;
        }
        if at == first_at && need_flag {
            out.data.extend_from_slice(&flag_line);
            out.data.extend_from_slice(eol);
            out.flag_added = true;
        }
        if v3 && py_strip(line).is_empty() {
            out.blank_removed += 1;
            continue;
        }
        out.data.extend_from_slice(line);
        out.data.extend_from_slice(eol);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geomaticist_gnss_core::time::gps_week_tow;

    const V3: &str =
        "     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
BASE                                                        MARKER NAME
G    2 C1C L1C                                              SYS / # / OBS TYPES
  2024     9    23    18     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
> 2024 09 23 18 00  0.0000000  0  2
G01  20000000.000   100000000.123
G02  21000000.000   110000000.123
> 2024 09 23 18 00  1.0000000  0  1
G01  20000001.000   100000001.123
> 2024 09 23 18 00  1.2500000  5  0
> 2024 09 23 18 00 01.750000 5  0
> 2024 09 23 18 00  2.0000000  0  1
G01  20000002.000   100000002.123
> 2024 09 23 18 00  3.0000000  0  1
G01  20000003.000   100000003.123
";

    fn t(sec: f64) -> f64 {
        // 2024-09-23 18:00:00 GPS
        ticks_to_seconds(to_ticks(2024, 9, 23, 18, 0, 0).unwrap()) + sec
    }

    #[test]
    fn reads_events_in_both_spellings() {
        let r = events_from_bytes(V3.as_bytes()).unwrap();
        assert_eq!(r.warnings.len(), 0);
        assert_eq!(r.events.len(), 2);
        assert!((r.events[0].t - t(1.25)).abs() < 1e-7);
        assert!((r.events[1].t - t(1.75)).abs() < 1e-7);
        assert_eq!(r.events[1].label, "2");
    }

    #[test]
    fn inserts_events_in_order_and_reads_them_back() {
        // the file without its events
        let plain = without(V3, " 5  0");
        assert_eq!(events_from_bytes(plain.as_bytes()).unwrap().events.len(), 0);
        let ev = [
            Event::at(t(2.5)),
            Event::at(t(0.4)),
            Event::at(t(0.123_456)),
            Event::at(t(-5.0)),
            Event::at(t(3.5)),
            Event::at(t(3.0)),
        ];
        let ins = insert_events(plain.as_bytes(), &ev).unwrap();
        assert_eq!(
            (ins.inserted, ins.before, ins.after, ins.existing),
            (4, 1, 1, 0)
        );
        assert_eq!(ins.warnings.len(), 1);
        let text = String::from_utf8(ins.data.clone()).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        // standard columns: flag in column 32, count in 33-35
        let e0 = lines.iter().position(|l| l.contains("  5  0")).unwrap();
        assert_eq!(lines[e0], "> 2024 09 23 18 00  0.1234560  5  0");
        assert_eq!(lines[e0 + 1], "> 2024 09 23 18 00  0.4000000  5  0");
        assert_eq!(lines[e0 + 2], "> 2024 09 23 18 00  1.0000000  0  1");
        // an event at an epoch's own time goes before that epoch
        let e3 = lines
            .iter()
            .position(|l| l.starts_with("> 2024 09 23 18 00  3.0000000  5"))
            .unwrap();
        assert_eq!(lines[e3 + 1], "> 2024 09 23 18 00  3.0000000  0  1");
        let back = events_from_bytes(&ins.data).unwrap();
        let times: Vec<f64> = back.events.iter().map(|e| e.t - t(0.0)).collect();
        assert_eq!(times.len(), 4);
        for (a, b) in times.iter().zip([0.123_456, 0.4, 2.5, 3.0]) {
            assert!((a - b).abs() < 2e-7, "{times:?}");
        }
        // the observation records are untouched, and the standard reader still reads the file
        assert_eq!(without(&text, "  5  0"), plain);
        let st = geomaticist_gnss_rinex::obs::read_obs_bytes(
            Path::new("x.obs"),
            &ins.data,
            "",
            None,
            None,
            "G",
        )
        .unwrap();
        assert_eq!(st.t.len(), 4);
    }

    const V2: &str =
        "     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
PT01                                                        MARKER NAME
     6    L1    L2    C1    P2    S1    S2                  # / TYPES OF OBSERV
                                                            END OF HEADER
 24  9 23 18  0  0.0000000  0  2G01G02
 100000000.123    78000000.456    20000000.000    20000004.000        45.000
        40.000
 110000000.123    88000000.456    21000000.000    21000004.000        44.000
        39.000
 24  9 23 18  0  1.0000000  0  1G01
 100000001.123    78000001.456    20000001.000    20000005.000        45.000
        40.000
                            2  0
 24  9 23 18  0  2.0000000  0  1G01
 100000002.123    78000002.456    20000002.000    20000006.000        45.000
        40.000
 24  9 23 18  0  2.5000000  3  2
PT02                                                        MARKER NAME
        1.8000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
 24  9 23 18  0  3.0000000  0  1G01
 100000003.123    78000003.456    20000003.000    20000007.000        45.000
        40.000
 24  9 23 18  0  4.0000000  0  1G01
 100000004.123    78000004.456    20000004.000    20000008.000        45.000
        40.000
";

    #[test]
    fn occupations_from_flags_2_and_3() {
        let oc = occupations_from_bytes(V2.as_bytes()).unwrap();
        assert_eq!(oc.len(), 2);
        assert_eq!(oc[0].name, "PT01");
        assert!((oc[0].t0 - t(0.0)).abs() < 1e-7 && (oc[0].t1 - t(1.0)).abs() < 1e-7);
        assert_eq!(oc[1].name, "PT02");
        assert!((oc[1].t0 - t(3.0)).abs() < 1e-7 && (oc[1].t1 - t(4.0)).abs() < 1e-7);
        assert_eq!(occupations_from_bytes(V3.as_bytes()).unwrap().len(), 0);
        // RINEX 2 events are written in the RINEX 2 layout
        let ins = insert_events(V2.as_bytes(), &[Event::at(t(0.5))]).unwrap();
        let text = String::from_utf8(ins.data).unwrap();
        assert!(
            text.contains(
                "\n 24  9 23 18  0  0.5000000  5  0\n 24  9 23 18  0  1.0000000  0  1G01\n"
            )
        );
    }

    /// The lines of `text` that do not contain `pat`.
    fn without(text: &str, pat: &str) -> String {
        let mut s = String::new();
        for l in text.lines().filter(|l| !l.contains(pat)) {
            s.push_str(l);
            s.push('\n');
        }
        s
    }

    fn file(epochs: &[&str]) -> String {
        let mut s = String::from(
            "     3.05           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    1 C1C                                                  SYS / # / OBS TYPES
  2024     9    28    23    40    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
",
        );
        for e in epochs {
            s.push_str("> ");
            s.push_str(e);
            s.push_str("  0  1\nG01  20000000.000\n");
        }
        s
    }

    #[test]
    fn finds_the_known_time_glitches() {
        // a 1980 first epoch, then three epochs seven days ahead just before the week change, one repeated
        let f = file(&[
            "1980 01 06 00 00  0.0000000",
            "2024 09 28 23 59 57.0000000",
            "2024 09 28 23 59 58.0000000",
            "2024 10 05 23 59 59.0000000",
            "2024 10 06 00 00  0.0000000",
            "2024 10 06 00 00  1.0000000",
            "2024 09 29 00 00  2.0000000",
            "2024 09 29 00 00  2.0000000",
            "2024 09 29 00 00  3.0000000",
        ]);
        let c = check_times(f.as_bytes()).unwrap();
        assert_eq!(
            (c.epochs, c.before_1990, c.displaced, c.not_later),
            (9, 1, 3, 1)
        );
        assert_eq!(c.warnings.len(), 4, "{:?}", c.warnings);
        assert!(c.warnings[1].contains("+7"));
        assert!(c.warnings[3].contains("TIME OF FIRST OBS"));
        let (a, b) = c.span.unwrap();
        assert_eq!(
            gps_week_tow(ticks_to_seconds(b)).0 - gps_week_tow(ticks_to_seconds(a)).0,
            1
        );
        // a clean file: nothing to say
        let ok = file(&["2024 09 28 23 40  0.0000000", "2024 09 28 23 40  1.0000000"]);
        let c = check_times(ok.as_bytes()).unwrap();
        assert!(c.warnings.is_empty(), "{:?}", c.warnings);
        // a file that really continues days later is not a glitch
        let two = file(&[
            "2024 09 28 23 40  0.0000000",
            "2024 09 30 23 40  0.0000000",
            "2024 09 30 23 40  1.0000000",
        ]);
        assert_eq!(check_times(two.as_bytes()).unwrap().displaced, 0);
    }

    #[test]
    fn marks_a_drone_file_kinematic_for_office_software() {
        // DJI-like: no MARKER TYPE, a long line of spaces after END OF HEADER, CRLF; events as the drone
        // export writes them
        let plain = without(V3, " 5  0");
        let with_ev =
            insert_events(plain.as_bytes(), &[Event::at(t(0.5)), Event::at(t(2.5))]).unwrap();
        let dji = String::from_utf8(with_ev.data)
            .unwrap()
            .replace(
                "                                                            END OF HEADER\n",
                &format!(
                    "                                                            END OF HEADER\n{}\n",
                    " ".repeat(300)
                ),
            )
            .replace('\n', "\r\n");
        let k = mark_kinematic(dji.as_bytes(), "AIRBORNE").unwrap();
        assert!(k.flag_added);
        assert_eq!(k.blank_removed, 1);
        let text = String::from_utf8(k.data.clone()).unwrap();
        let lines: Vec<&str> = text.split("\r\n").collect();
        let eoh = lines
            .iter()
            .position(|l| l.ends_with("END OF HEADER"))
            .unwrap();
        assert_eq!(
            lines[eoh - 1],
            format!("{:<60}{:<20}", "AIRBORNE", "MARKER TYPE")
        );
        // flag 2 in column 32, count 0, blank date: right before the first epoch
        assert_eq!(lines[eoh + 1], format!(">{}2  0", " ".repeat(30)));
        assert_eq!(lines[eoh + 2], "> 2024 09 23 18 00  0.0000000  0  2");
        assert!(!text.contains("\r\n\r\n"), "no blank lines left");
        assert!(text.ends_with("\r\n") && !text.contains("\n\n"));
        // the events and the observations read as before
        assert_eq!(events_from_bytes(&k.data).unwrap().events.len(), 2);
        let st = geomaticist_gnss_rinex::obs::read_obs_bytes(
            Path::new("x.obs"),
            &k.data,
            "",
            None,
            None,
            "G",
        )
        .unwrap();
        assert_eq!(st.t.len(), 4);
        // twice is the same as once; an existing MARKER TYPE is replaced, not doubled
        let again = mark_kinematic(&k.data, "AIRBORNE").unwrap();
        assert!(!again.flag_added);
        assert_eq!(again.data, k.data);
        let geo = V3.replace(
            "BASE                                                        MARKER NAME\n",
            "BASE                                                        MARKER NAME\n\
             GEODETIC                                                    MARKER TYPE\n",
        );
        let g =
            String::from_utf8(mark_kinematic(geo.as_bytes(), "AIRBORNE").unwrap().data).unwrap();
        assert_eq!(g.matches("MARKER TYPE").count(), 1);
        assert!(g.contains("AIRBORNE      "));
    }
}
