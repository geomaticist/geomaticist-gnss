//! Static and kinematic segments of a Leica MDB log, and the point records a Captivate job's .m00 holds.
//!
//! A rover that logs raw data while it measures RTK points writes a moving (kinematic) log with stops on the
//! points. The log says so itself, in its event records (type 0x83, 21 bytes from the type byte):
//!
//! ```text
//! [0x83][u16 GPS week][u32 ms of week][u32 link][u16 0][u8 code][u16 week][u32 ms][u8 3]
//! ```
//!
//! | code | meaning | link | time (bytes 14..20) |
//! |---|---|---|---|
//! | 0x0c | kinematic segment starts | the record's own file offset | its first epoch |
//! | 0x0d | static segment (an occupation) starts | its own offset | its first epoch |
//! | 0x0e | kinematic segment ends | offset of the 0x0c that started it | its last epoch |
//! | 0x0f | static segment ends | offset of the 0x0d that started it | its last epoch |
//!
//! A base log has one 0x0d at its first epoch and one 0x0f at its last. A rover job log starts with 0x0c, then
//! for each point a pair (0x0e, 0x0d) at the occupation start and a pair (0x0f, 0x0c) at its end, the links
//! pointing back as above. Each occupation has a point record (0x71) without a name at its start and one with
//! the point's name and its X01 record number at its end (`point_record`).
//!
//! Type 0x82 (8 bytes) is logging on (byte 7 = 0) / off (1): a rover log can hold several logging sessions.

use super::{T_POINT, records, t_of, utf16};
use crate::util::le;

/// The event record type (segment starts and ends).
pub const T_EVENT: u8 = 0x83;
/// Logging on / off.
pub const T_LOGGING: u8 = 0x82;
/// Event code: a kinematic segment starts.
pub const EV_KIN_START: u8 = 0x0c;
/// Event code: a static segment (occupation) starts.
pub const EV_STATIC_START: u8 = 0x0d;
/// Event code: a kinematic segment ends.
pub const EV_KIN_END: u8 = 0x0e;
/// Event code: a static segment ends.
pub const EV_STATIC_END: u8 = 0x0f;

/// GPS seconds of a timed record (week and ms of week after the type byte); `None` for an implausible week.
#[must_use]
pub fn record_time(body: &[u8]) -> Option<f64> {
    if body.len() < 7 {
        return None;
    }
    let week = le::u16(body, 1);
    if !(1000..5000).contains(&week) {
        return None;
    }
    Some((week * 604_800) as f64 + le::u32(body, 3) as f64 / 1000.0)
}

/// One event record (0x83).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Event {
    /// File offset of the record.
    pub offset: usize,
    /// Its code (`EV_*`).
    pub code: u8,
    /// The segment's first epoch (a start) or last epoch (an end), GPS seconds.
    pub time: f64,
    /// For a start, its own offset; for an end, the offset of the start record it closes.
    pub link: usize,
}

/// The event record whose body (from the type byte) is `body`, at `offset`.
#[must_use]
pub fn event(offset: usize, body: &[u8]) -> Option<Event> {
    if body.first() != Some(&T_EVENT) || body.len() < 20 {
        return None;
    }
    let code = body[13];
    if !(EV_KIN_START..=EV_STATIC_END).contains(&code) {
        return None;
    }
    let time = record_time(&body[13..]).or_else(|| record_time(body))?;
    Some(Event {
        offset,
        code,
        time,
        link: le::u32(body, 7) as usize,
    })
}

/// A point record (0x71): what the receiver stored for a point.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointRecord {
    /// File offset of the record (X01 field 0x20 points here).
    pub offset: usize,
    /// GPS seconds.
    pub time: Option<f64>,
    /// The point's name (empty on the record a job log writes at an occupation's start).
    pub name: String,
    /// The X01 record number (0 when none).
    pub record: u32,
    /// Antenna height as entered, m (to the height hook on a tripod; the pole height on a pole).
    pub antenna_height: f64,
    /// Latitude, longitude (degrees), ellipsoid height (m), when the record holds them.
    pub llh: Option<(f64, f64, f64)>,
    /// Covariance, m² (north, east, up order: nn ne nu ee eu uu), when the record holds it.
    pub cov: Option<[f64; 6]>,
}

/// The point record at `offset` of the log (`None` if there's none there).
#[must_use]
pub fn point_record(data: &[u8], offset: usize) -> Option<PointRecord> {
    if data.get(offset..offset + 2)? != b"\x9c\xae" {
        return None;
    }
    let n = usize::from(u16::from_le_bytes([
        *data.get(offset + 2)?,
        *data.get(offset + 3)?,
    ]));
    let body = data.get(offset + 4..offset + n)?;
    point_from_body(offset, body)
}

fn point_from_body(offset: usize, body: &[u8]) -> Option<PointRecord> {
    if body.first() != Some(&T_POINT) || body.len() < 67 {
        return None;
    }
    let mut p = PointRecord {
        offset,
        time: record_time(body),
        name: utf16(&body[25..57]).trim().to_owned(),
        record: u32::try_from(le::u32(body, 59)).unwrap_or(0),
        antenna_height: le::f32(body, 63),
        ..PointRecord::default()
    };
    if body.len() >= 103 {
        let (lat, lon, h) = (le::f64(body, 79), le::f64(body, 87), le::f64(body, 95));
        if lat.abs() <= std::f64::consts::FRAC_PI_2 && lon.abs() <= 2.0 * std::f64::consts::PI {
            p.llh = Some((lat.to_degrees(), lon.to_degrees(), h));
        }
    }
    if body.len() >= 132 {
        let mut c = [0.0; 6];
        for (k, v) in c.iter_mut().enumerate() {
            *v = le::f32(body, 108 + 4 * k);
        }
        p.cov = Some(c);
    }
    Some(p)
}

/// Every point record of the log, in file order.
#[must_use]
pub fn point_records(data: &[u8]) -> Vec<PointRecord> {
    records(data)
        .into_iter()
        .filter(|r| r.1 == T_POINT)
        .filter_map(|(pos, _, body)| point_from_body(pos, body))
        .collect()
}

/// Static (the antenna on a point) or kinematic (moving).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentKind {
    /// An occupation.
    Static,
    /// The antenna moves.
    Kinematic,
}

/// A stretch of the log between a start and an end event.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    /// Static or kinematic.
    pub kind: SegmentKind,
    /// First epoch, GPS seconds.
    pub start: f64,
    /// Last epoch, GPS seconds.
    pub end: f64,
    /// File offset of its start record.
    pub start_offset: usize,
    /// The named point record stored for it (a static segment of a job log), if any.
    pub point: Option<PointRecord>,
    /// The first point record inside it (a job log writes one without a name at the occupation's start).
    pub first_record: Option<PointRecord>,
}

impl Segment {
    /// Seconds from its first to its last epoch.
    #[must_use]
    pub fn seconds(&self) -> f64 {
        self.end - self.start
    }

    /// The point name the receiver stored for it ("" when none).
    #[must_use]
    pub fn point_name(&self) -> &str {
        self.point.as_ref().map_or("", |p| p.name.as_str())
    }

    /// The antenna height as entered for it (the named record's, else the first record's).
    #[must_use]
    pub fn antenna_height(&self) -> Option<f64> {
        self.point
            .as_ref()
            .or(self.first_record.as_ref())
            .map(|p| p.antenna_height)
    }

    /// A name for a stop with no point name: `STOP-<day of year>-<hhmmss>` (UTC) of its first epoch.
    #[must_use]
    pub fn generated_name(&self) -> String {
        const CUM: [i64; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
        let t = t_of(self.start - leap_at(self.start));
        let leap = (t.year % 4 == 0 && t.year % 100 != 0) || t.year % 400 == 0;
        let m = usize::try_from(t.month - 1).unwrap_or(0).min(11);
        let doy = CUM[m] + t.day + i64::from(leap && t.month > 2);
        format!("STOP-{doy:03}-{:02}{:02}{:02}", t.hour, t.minute, t.second)
    }

    /// The stored point name, or [`Segment::generated_name`].
    #[must_use]
    pub fn name(&self) -> String {
        let n = self.point_name();
        if n.is_empty() {
            self.generated_name()
        } else {
            n.to_owned()
        }
    }
}

/// GPS − UTC at a GPS time (the leap seconds since 2006 are enough for these logs).
fn leap_at(gps: f64) -> f64 {
    // 2017-01-01, 2015-07-01, 2012-07-01, 2009-01-01 as GPS seconds
    if gps >= 1_167_264_018.0 {
        18.0
    } else if gps >= 1_119_744_017.0 {
        17.0
    } else if gps >= 1_025_136_016.0 {
        16.0
    } else if gps >= 914_803_215.0 {
        15.0
    } else {
        14.0
    }
}

/// The segments of a log and what the events say about it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LogMotion {
    /// Segments in time order.
    pub segments: Vec<Segment>,
    /// Event records read.
    pub events: usize,
    /// Logging sessions (0x82 on records).
    pub logging_sessions: usize,
    /// First and last epoch of the log.
    pub span: Option<(f64, f64)>,
}

impl LogMotion {
    /// The log has at least one kinematic segment: the antenna moved (a rover's job log).
    #[must_use]
    pub fn is_kinematic(&self) -> bool {
        self.segments
            .iter()
            .any(|s| s.kind == SegmentKind::Kinematic)
    }

    /// The static segments (occupations).
    pub fn stops(&self) -> impl Iterator<Item = &Segment> {
        self.segments
            .iter()
            .filter(|s| s.kind == SegmentKind::Static)
    }

    /// Seconds spent moving.
    #[must_use]
    pub fn moving_seconds(&self) -> f64 {
        self.segments
            .iter()
            .filter(|s| s.kind == SegmentKind::Kinematic)
            .map(Segment::seconds)
            .sum()
    }
}

/// Reads the segments of a log from its event records (see the module docs). A segment still open at a
/// logging-off record or at the end of the log ends at the last epoch before it. A log without event records
/// gives no segments (older firmware: treat it as static).
#[must_use]
pub fn motion(data: &[u8]) -> LogMotion {
    let mut out = LogMotion::default();
    let mut open: Option<Segment> = None;
    let mut last_epoch: Option<f64> = None;
    let mut first_epoch: Option<f64> = None;
    let close = |seg: Option<Segment>, end: Option<f64>, out: &mut LogMotion| {
        if let Some(mut s) = seg {
            if let Some(e) = end {
                s.end = e.max(s.start);
            }
            out.segments.push(s);
        }
    };
    for (pos, typ, body) in records(data) {
        match typ {
            super::epochs::T_EPOCH => {
                if let Some(t) = record_time(body) {
                    first_epoch.get_or_insert(t);
                    last_epoch = Some(t);
                }
            }
            T_LOGGING => {
                if body.get(7) == Some(&0) {
                    out.logging_sessions += 1;
                } else {
                    close(open.take(), last_epoch, &mut out);
                }
            }
            T_POINT => {
                let Some(p) = point_from_body(pos, body) else {
                    continue;
                };
                let named = !p.name.is_empty() || p.record != 0;
                match open.as_mut() {
                    Some(s) if s.kind == SegmentKind::Static => {
                        if s.first_record.is_none() {
                            s.first_record = Some(p.clone());
                        }
                        if named {
                            s.point = Some(p);
                        }
                    }
                    other => {
                        // a job log stores the named record just after the occupation's end events
                        if named
                            && let Some(s) = out.segments.last_mut()
                            && s.kind == SegmentKind::Static
                            && s.point.is_none()
                        {
                            s.point = Some(p);
                        } else if let Some(s) = other
                            && s.first_record.is_none()
                        {
                            s.first_record = Some(p);
                        }
                    }
                }
            }
            T_EVENT => {
                let Some(ev) = event(pos, body) else {
                    continue;
                };
                out.events += 1;
                match ev.code {
                    EV_KIN_START | EV_STATIC_START => {
                        close(open.take(), last_epoch, &mut out);
                        open = Some(Segment {
                            kind: if ev.code == EV_KIN_START {
                                SegmentKind::Kinematic
                            } else {
                                SegmentKind::Static
                            },
                            start: ev.time,
                            end: ev.time,
                            start_offset: pos,
                            point: None,
                            first_record: None,
                        });
                    }
                    _ => {
                        if let Some(s) = open.as_ref()
                            && (s.start_offset == ev.link
                                || (s.kind == SegmentKind::Kinematic) == (ev.code == EV_KIN_END))
                        {
                            close(open.take(), Some(ev.time), &mut out);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    close(open.take(), last_epoch, &mut out);
    if let (Some(a), Some(b)) = (first_epoch, last_epoch) {
        out.span = Some((a, b));
    }
    out
}

/// True when the log has a kinematic segment (reads only the event records; a cheap whole-file scan).
#[must_use]
pub fn is_kinematic(data: &[u8]) -> bool {
    records(data)
        .into_iter()
        .filter(|r| r.1 == T_EVENT)
        .any(|(pos, _, body)| event(pos, body).is_some_and(|e| e.code == EV_KIN_START))
}

#[doc(hidden)]
pub mod synth {
    //! Synthetic MDB records for tests (no real data); public so the app's tests can build them too.

    /// A framed record: `[9C AE][u16 n][body][u16 check]`.
    #[must_use]
    pub fn frame(body: &[u8]) -> Vec<u8> {
        let n = u16::try_from(body.len() + 4).unwrap_or(u16::MAX);
        let mut v = vec![0x9c, 0xae];
        v.extend_from_slice(&n.to_le_bytes());
        v.extend_from_slice(body);
        v.extend_from_slice(&[0, 0]);
        v
    }

    /// Week and ms of week of GPS seconds.
    #[must_use]
    pub fn wk(t: f64) -> (u16, u32) {
        let week = (t / 604_800.0).floor();
        (
            week as u16,
            ((t - week * 604_800.0) * 1000.0).round() as u32,
        )
    }

    /// An event record body.
    #[must_use]
    pub fn event(t_rec: f64, link: u32, code: u8, t_seg: f64) -> Vec<u8> {
        let (w, ms) = wk(t_rec);
        let (w2, ms2) = wk(t_seg);
        let mut b = vec![super::T_EVENT];
        b.extend_from_slice(&w.to_le_bytes());
        b.extend_from_slice(&ms.to_le_bytes());
        b.extend_from_slice(&link.to_le_bytes());
        b.extend_from_slice(&[0, 0, code]);
        b.extend_from_slice(&w2.to_le_bytes());
        b.extend_from_slice(&ms2.to_le_bytes());
        b.push(3);
        b
    }

    /// A point record body (job-log layout, 238 bytes).
    #[must_use]
    pub fn point(
        t: f64,
        name: &str,
        record: u32,
        height: f32,
        llh_deg: (f64, f64, f64),
    ) -> Vec<u8> {
        let (w, ms) = wk(t);
        let mut b = vec![0u8; 238];
        b[0] = super::T_POINT;
        b[1..3].copy_from_slice(&w.to_le_bytes());
        b[3..7].copy_from_slice(&ms.to_le_bytes());
        for (k, u) in name.encode_utf16().take(15).enumerate() {
            b[25 + 2 * k..27 + 2 * k].copy_from_slice(&u.to_le_bytes());
        }
        b[59..63].copy_from_slice(&record.to_le_bytes());
        b[63..67].copy_from_slice(&height.to_le_bytes());
        b[79..87].copy_from_slice(&llh_deg.0.to_radians().to_le_bytes());
        b[87..95].copy_from_slice(&llh_deg.1.to_radians().to_le_bytes());
        b[95..103].copy_from_slice(&llh_deg.2.to_le_bytes());
        for k in 0..6 {
            let v: f32 = if k == 0 || k == 3 {
                1e-4
            } else if k == 5 {
                4e-4
            } else {
                0.0
            };
            b[108 + 4 * k..112 + 4 * k].copy_from_slice(&v.to_le_bytes());
        }
        b
    }

    /// A minimal observation record body (no satellites) at `t`.
    #[must_use]
    pub fn epoch(t: f64) -> Vec<u8> {
        let (w, ms) = wk(t);
        let mut b = vec![super::super::epochs::T_EPOCH];
        b.extend_from_slice(&w.to_le_bytes());
        b.extend_from_slice(&ms.to_le_bytes());
        b.extend_from_slice(&[3, 0, 0xb4, 0x18, 0x78, 0x66, 0]);
        b
    }

    /// A job log: kinematic from `t0`, then one occupation per `(start, seconds, name)`, moving between them;
    /// epochs every second; ends `tail` seconds after the last stop. → the log's bytes.
    #[must_use]
    pub fn job_log(t0: f64, stops: &[(f64, f64, &str)], tail: f64) -> Vec<u8> {
        let mut d = frame(&[super::super::T_HEADER, 0, 0, 0, 0, 0, 0]);
        let end = stops.last().map_or(t0, |s| s.0 + s.1) + tail;
        let off = |d: &Vec<u8>| u32::try_from(d.len()).unwrap_or(0);
        let mut kin_start = off(&d);
        d.extend(frame(&event(t0, kin_start, super::EV_KIN_START, t0)));
        let mut stat_start = 0u32;
        let mut t = t0;
        let mut k = 0usize;
        while t <= end {
            if let Some(&(s0, secs, name)) = stops.get(k) {
                if (t - s0).abs() < 1e-6 {
                    d.extend(frame(&event(t, kin_start, super::EV_KIN_END, t - 1.0)));
                    stat_start = off(&d);
                    d.extend(frame(&event(t, stat_start, super::EV_STATIC_START, t)));
                    d.extend(frame(&point(t, "", 0, 1.8, (46.0, 8.0, 100.0))));
                } else if (t - (s0 + secs)).abs() < 1e-6 {
                    d.extend(frame(&event(t, stat_start, super::EV_STATIC_END, t - 1.0)));
                    kin_start = off(&d);
                    d.extend(frame(&event(t, kin_start, super::EV_KIN_START, t)));
                    let rec = u32::try_from(2 * k + 6).unwrap_or(0);
                    let lat = 46.0 + 0.001 * k as f64;
                    d.extend(frame(&point(s0, name, rec, 1.8, (lat, 8.0, 100.0))));
                    k += 1;
                }
            }
            d.extend(frame(&epoch(t)));
            t += 1.0;
        }
        d
    }
}

#[cfg(test)]
mod tests {
    use super::synth::{epoch, frame, job_log, point};
    use super::*;

    const T0: f64 = 1_411_084_800.0;

    #[test]
    fn job_log_segments() {
        let d = job_log(
            T0,
            &[(T0 + 100.0, 4.0, "2001"), (T0 + 300.0, 700.0, "")],
            50.0,
        );
        let m = motion(&d);
        assert!(m.is_kinematic());
        assert!(is_kinematic(&d));
        let kinds: Vec<SegmentKind> = m.segments.iter().map(|s| s.kind).collect();
        assert_eq!(
            kinds,
            [
                SegmentKind::Kinematic,
                SegmentKind::Static,
                SegmentKind::Kinematic,
                SegmentKind::Static,
                SegmentKind::Kinematic
            ]
        );
        let stops: Vec<&Segment> = m.stops().collect();
        assert_eq!(stops.len(), 2);
        assert_eq!(stops[0].start, T0 + 100.0);
        assert_eq!(stops[0].end, T0 + 103.0);
        assert_eq!(stops[0].point_name(), "2001");
        assert_eq!(stops[0].point.as_ref().map(|p| p.record), Some(6));
        assert!((stops[0].antenna_height().unwrap_or(0.0) - 1.8).abs() < 1e-6);
        assert_eq!(stops[1].seconds(), 699.0);
        // the second stop has no point name: a generated one
        assert_eq!(stops[1].point_name(), "");
        assert!(stops[1].name().starts_with("STOP-"));
        assert_eq!(m.span, Some((T0, T0 + 1050.0)));
        assert!(m.moving_seconds() > 300.0);
        assert_eq!(point_records(&d).len(), 4);
    }

    #[test]
    fn base_log_is_static() {
        // a base log: one static segment over the whole log, no kinematic start
        let mut d = frame(&[crate::mdb::T_HEADER, 0, 0, 0, 0, 0, 0]);
        let at = u32::try_from(d.len()).unwrap_or(0);
        d.extend(frame(&super::synth::event(T0, at, EV_STATIC_START, T0)));
        d.extend(frame(&point(T0, "103", 0, 1.25, (46.0, 8.0, 100.0))));
        for k in 0..100 {
            d.extend(frame(&epoch(T0 + f64::from(k))));
        }
        d.extend(frame(&super::synth::event(
            T0 + 99.0,
            at,
            EV_STATIC_END,
            T0 + 99.0,
        )));
        let m = motion(&d);
        assert!(!m.is_kinematic());
        assert!(!is_kinematic(&d));
        assert_eq!(m.segments.len(), 1);
        assert_eq!(m.segments[0].seconds(), 99.0);
        assert_eq!(m.segments[0].point_name(), "103");
    }

    #[test]
    fn generated_names() {
        // 2024-09-23 00:07:11 UTC = GPS 1411085249 (18 s leap)
        let s = Segment {
            kind: SegmentKind::Static,
            start: 1_411_085_249.0,
            end: 1_411_085_253.0,
            start_offset: 0,
            point: None,
            first_record: None,
        };
        assert_eq!(s.generated_name(), "STOP-267-000711");
    }

    #[test]
    fn point_record_at_offset() {
        let mut d = frame(&[crate::mdb::T_HEADER, 0, 0, 0, 0, 0, 0]);
        let off = d.len();
        d.extend(frame(&point(T0, "K7", 12, 2.0, (46.5, 8.25, 55.0))));
        let p = point_record(&d, off).unwrap_or_default();
        assert_eq!(p.name, "K7");
        assert_eq!(p.record, 12);
        assert_eq!(p.time, Some(T0));
        let (la, lo, h) = p.llh.unwrap_or_default();
        assert!((la - 46.5).abs() < 1e-12 && (lo - 8.25).abs() < 1e-12 && h == 55.0);
        assert!(point_record(&d, 0).is_none());
        assert!(point_record(&d, off + 1).is_none());
    }
}
