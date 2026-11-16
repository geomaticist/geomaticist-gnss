//! DJI `Timestamp.MRK` files: one line per photo.
//!
//! A line holds tab-separated fields (the tabs are shown as blanks here; the GPS week is absent in some
//! firmware):
//!
//! ```text
//! 1   345600.250000   [2300]     -120,N       45,E      193,V   45.00000001,Lat   7.00000002,Lon   512.345,Ellh   0.012000, 0.011000, 0.025000   50,Q
//! ```
//!
//! | field | meaning |
//! |---|---|
//! | `1` | sequence number of the photo in the flight |
//! | `345600.250000` | exposure time, GPS seconds of the week |
//! | `[2300]` | GPS week |
//! | `-120,N` `45,E` `193,V` | offset from the antenna phase centre to the camera's sensor centre at the exposure, millimetres, in the local level frame: north, east, and **down** (`V` is positive when the camera is below the antenna, which it is on these aircraft) |
//! | `...,Lat` `...,Lon` `...,Ellh` | the aircraft's own position at the exposure, degrees and metres |
//! | three numbers | its standard deviations: north, east, vertical, metres |
//! | `50,Q` | its quality: 50 fixed RTK, 34 float RTK, 16 single point (no corrections), 0 none |
//!
//! The offset changes from photo to photo because it is the fixed lever arm between antenna and camera
//! rotated by the aircraft's attitude and the gimbal angles at that moment. Applying it is simple vector
//! addition: camera = antenna + (N, E) horizontally, and the camera's height = antenna height − V.
//!
//! Times: the seconds of the week restart at 0 when a new GPS week begins (Sunday 00:00:00 GPS time). A
//! flight across that moment has lines with a small time after lines with a large one. [`read`] handles
//! it (the week is advanced, with a warning, when the file's own week number did not change).

use std::path::Path;

use geomaticist_gnss_core::consts::WEEK;
use geomaticist_gnss_core::event::{Event, EventOffset};
use serde::{Deserialize, Serialize};

use crate::EventError;

/// One line of a `Timestamp.MRK` file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MrkRecord {
    /// Sequence number (as written).
    pub seq: String,
    /// GPS week (after [`read`]: resolved, also across a week change).
    pub week: Option<i64>,
    /// GPS seconds of the week.
    pub tow: f64,
    /// Antenna phase centre to camera, mm: north, east, down. `None` when the line has no offset.
    pub offset_mm: Option<[f64; 3]>,
    /// The aircraft's own latitude, degrees.
    pub lat: f64,
    /// Its longitude, degrees.
    pub lon: f64,
    /// Its ellipsoidal height, m.
    pub ellh: f64,
    /// Its standard deviations: north, east, vertical (m).
    pub sigma: [f64; 3],
    /// Its quality flag: 50 fixed, 34 float, 16 single point; -1 when absent.
    pub q: i32,
}

impl MrkRecord {
    /// GPS seconds since 1980-01-06 (`None` while the week is unknown).
    #[must_use]
    pub fn gps_seconds(&self) -> Option<f64> {
        self.week.map(|w| w as f64 * WEEK + self.tow)
    }

    /// The offset as the solver takes it.
    #[must_use]
    pub fn offset(&self) -> Option<EventOffset> {
        self.offset_mm.map(|o| EventOffset {
            north: o[0] / 1000.0,
            east: o[1] / 1000.0,
            down: o[2] / 1000.0,
        })
    }

    /// What the quality flag means.
    #[must_use]
    pub fn quality(&self) -> &'static str {
        match self.q {
            50 => "fixed",
            34 => "float",
            16 => "single",
            0 => "none",
            _ => "unknown",
        }
    }
}

/// A `Timestamp.MRK` file.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Mrk {
    /// The photos, in file order.
    pub records: Vec<MrkRecord>,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

impl Mrk {
    /// The events of the file (records whose week is unknown are left out; [`Mrk::resolve_weeks`] first).
    /// Labels are the sequence numbers; the offsets are the antenna-to-camera vectors.
    #[must_use]
    pub fn events(&self) -> Vec<Event> {
        self.records
            .iter()
            .filter_map(|r| {
                Some(Event {
                    t: r.gps_seconds()?,
                    label: r.seq.clone(),
                    offset: r.offset(),
                })
            })
            .collect()
    }

    /// Gives every record without a week the week that puts it nearest to `near` (GPS seconds, the first
    /// epoch of the rover file for example), and re-checks the order.
    pub fn resolve_weeks(&mut self, near: f64) {
        for r in &mut self.records {
            if r.week.is_none() {
                r.week = Some(((near - r.tow) / WEEK).round() as i64);
            }
        }
        self.order_weeks();
    }

    /// Advances the week where the seconds of the week restart without the week number following.
    fn order_weeks(&mut self) {
        let mut add = 0i64;
        let mut fixed = 0usize;
        let mut prev: Option<f64> = None;
        for r in &mut self.records {
            let Some(w) = r.week else { continue };
            let mut t = (w + add) as f64 * WEEK + r.tow;
            if let Some(p) = prev
                && t < p - 0.5 * WEEK
            {
                add += 1;
                t += WEEK;
            }
            if add > 0 {
                r.week = Some(w + add);
                fixed += 1;
            }
            prev = Some(t);
        }
        if fixed > 0 {
            self.warnings.push(format!(
                "the GPS week changed during the flight but the file's week number did not follow: \
                 {fixed} photos moved to the next week"
            ));
        }
    }

    /// How many records have each quality flag, e.g. `[(16, 695)]`.
    #[must_use]
    pub fn quality_counts(&self) -> Vec<(i32, usize)> {
        let mut m = std::collections::BTreeMap::new();
        for r in &self.records {
            *m.entry(r.q).or_insert(0usize) += 1;
        }
        m.into_iter().collect()
    }
}

fn number(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok().filter(|v| v.is_finite())
}

fn parse_line(line: &str) -> Option<MrkRecord> {
    let fields: Vec<String> = if line.contains('\t') {
        line.split('\t').map(|s| s.trim().to_owned()).collect()
    } else {
        // some tools rewrite the tabs as blanks: a field that ends with a comma (inside the triple of
        // standard deviations) continues in the next word
        let mut out: Vec<String> = Vec::new();
        for word in line.split_whitespace() {
            match out.last_mut() {
                Some(last) if last.ends_with(',') => {
                    last.push(' ');
                    last.push_str(word);
                }
                _ => out.push(word.to_owned()),
            }
        }
        out
    };
    if fields.len() < 2 {
        return None;
    }
    let tow = number(&fields[1])?;
    if !(0.0..WEEK).contains(&tow) {
        return None;
    }
    let mut r = MrkRecord {
        seq: fields[0].clone(),
        week: None,
        tow,
        offset_mm: None,
        lat: f64::NAN,
        lon: f64::NAN,
        ellh: f64::NAN,
        sigma: [f64::NAN; 3],
        q: -1,
    };
    let mut off = [f64::NAN; 3];
    for f in &fields[2..] {
        if let Some(w) = f.strip_prefix('[').and_then(|x| x.strip_suffix(']')) {
            r.week = w
                .trim()
                .parse::<i64>()
                .ok()
                .filter(|w| (1..10_000).contains(w));
            continue;
        }
        let parts: Vec<&str> = f.split(',').map(str::trim).collect();
        match parts.as_slice() {
            [v, "N"] => off[0] = number(v).unwrap_or(f64::NAN),
            [v, "E"] => off[1] = number(v).unwrap_or(f64::NAN),
            [v, "V"] => off[2] = number(v).unwrap_or(f64::NAN),
            [v, "Lat"] => r.lat = number(v).unwrap_or(f64::NAN),
            [v, "Lon"] => r.lon = number(v).unwrap_or(f64::NAN),
            [v, "Ellh"] => r.ellh = number(v).unwrap_or(f64::NAN),
            [v, "Q"] => r.q = number(v).map_or(-1, |q| q as i32),
            [a, b, c] => {
                if let (Some(a), Some(b), Some(c)) = (number(a), number(b), number(c)) {
                    r.sigma = [a, b, c];
                }
            }
            _ => {}
        }
    }
    if off.iter().all(|v| v.is_finite()) {
        r.offset_mm = Some(off);
    }
    Some(r)
}

/// Parses the text of a `Timestamp.MRK` file.
#[must_use]
pub fn parse(text: &str) -> Mrk {
    let mut mrk = Mrk::default();
    let mut bad = 0usize;
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        match parse_line(line) {
            Some(r) => mrk.records.push(r),
            None => bad += 1,
        }
    }
    if bad > 0 {
        mrk.warnings
            .push(format!("{bad} lines could not be read and were left out"));
    }
    let no_week = mrk.records.iter().filter(|r| r.week.is_none()).count();
    if no_week > 0 {
        mrk.warnings.push(format!(
            "{no_week} of {} lines have no GPS week: it is taken from the observation file's time",
            mrk.records.len()
        ));
    }
    let no_off = mrk.records.iter().filter(|r| r.offset_mm.is_none()).count();
    if no_off > 0 && !mrk.records.is_empty() {
        mrk.warnings.push(format!(
            "{no_off} of {} lines have no antenna-to-camera offset: those positions stay at the antenna",
            mrk.records.len()
        ));
    }
    // a lever arm of several metres is not a lever arm
    let wild = mrk
        .records
        .iter()
        .filter(|r| {
            r.offset_mm
                .is_some_and(|o| o.iter().any(|v| v.abs() > 5000.0))
        })
        .count();
    if wild > 0 {
        mrk.warnings.push(format!(
            "{wild} lines have an antenna-to-camera offset above 5 m: check the file"
        ));
    }
    mrk.order_weeks();
    // times must increase
    let t: Vec<f64> = mrk
        .records
        .iter()
        .filter_map(MrkRecord::gps_seconds)
        .collect();
    let back = t.windows(2).filter(|w| w[1] <= w[0]).count();
    if back > 0 {
        mrk.warnings.push(format!(
            "{back} photos are not later than the photo before them"
        ));
    }
    mrk
}

/// Reads a `Timestamp.MRK` file.
pub fn read(path: &Path) -> Result<Mrk, EventError> {
    let bytes = std::fs::read(path).map_err(|e| EventError(format!("{}: {e}", path.display())))?;
    let text = String::from_utf8_lossy(&bytes);
    let mrk = parse(&text);
    if mrk.records.is_empty() {
        return Err(EventError(format!(
            "{}: no photo lines found (is this a Timestamp.MRK file?)",
            path.display()
        )));
    }
    Ok(mrk)
}

#[cfg(test)]
mod tests {
    use super::*;

    // synthetic lines in the layout of the format
    const TEXT: &str = "1\t345600.250000\t[2300]\t  -120,N\t    45,E\t   193,V\t45.00000001,Lat\t7.00000002,Lon\t512.345,Ellh\t0.012000, 0.011000, 0.025000\t50,Q\n\
2\t345602.750000\t[2300]\t   118,N\t   -44,E\t   195,V\t45.00010001,Lat\t7.00000002,Lon\t512.300,Ellh\t1.300000, 0.900000, 2.400000\t16,Q\n";

    #[test]
    fn reads_the_fields() {
        let m = parse(TEXT);
        assert!(m.warnings.is_empty(), "{:?}", m.warnings);
        assert_eq!(m.records.len(), 2);
        let r = &m.records[0];
        assert_eq!(r.seq, "1");
        assert_eq!(r.week, Some(2300));
        assert_eq!(r.offset_mm, Some([-120.0, 45.0, 193.0]));
        assert_eq!(
            (r.lat, r.lon, r.ellh),
            (45.000_000_01, 7.000_000_02, 512.345)
        );
        assert_eq!(r.sigma, [0.012, 0.011, 0.025]);
        assert_eq!((r.q, r.quality()), (50, "fixed"));
        assert_eq!(m.records[1].quality(), "single");
        assert_eq!(m.quality_counts(), vec![(16, 1), (50, 1)]);
        let ev = m.events();
        assert!((ev[0].t - (2300.0 * 604_800.0 + 345_600.25)).abs() < 1e-6);
        assert_eq!(ev[1].label, "2");
        let o = ev[0].offset.unwrap();
        // north, east as written; V is down
        assert_eq!((o.north, o.east, o.down), (-0.120, 0.045, 0.193));
        assert_eq!(o.enu(), [0.045, -0.120, -0.193]);
    }

    #[test]
    fn week_changes_and_missing_weeks() {
        // the week number follows the rollover: nothing to warn about
        let ok = "1\t604799.500000\t[2300]\t1,N\t2,E\t3,V\n2\t1.500000\t[2301]\t1,N\t2,E\t3,V\n";
        let m = parse(ok);
        assert!(m.warnings.is_empty(), "{:?}", m.warnings);
        let ev = m.events();
        assert!((ev[1].t - ev[0].t - 2.0).abs() < 1e-6);
        // the week number does not follow: advanced, with a warning
        let stuck = "1\t604799.500000\t[2300]\t1,N\t2,E\t3,V\n2\t1.500000\t[2300]\t1,N\t2,E\t3,V\n";
        let m = parse(stuck);
        assert_eq!(m.records[1].week, Some(2301));
        assert!(m.warnings.iter().any(|w| w.contains("GPS week changed")));
        // no week at all: resolved from a nearby time
        let none = "1\t604799.500000\t-5,N\t2,E\t3,V\n2\t1.500000\t1,N\t2,E\t3,V\n";
        let mut m = parse(none);
        assert_eq!(m.events().len(), 0);
        m.resolve_weeks(2300.0 * 604_800.0 + 604_000.0);
        let ev = m.events();
        assert_eq!(ev.len(), 2);
        assert!((ev[1].t - ev[0].t - 2.0).abs() < 1e-6, "{:?}", m.records);
        assert_eq!(m.records[0].week, Some(2300));
    }

    #[test]
    fn rejects_junk_and_reports() {
        let m = parse(
            "hello world\n\n3\tnot-a-time\t[2300]\n4\t100.0\t[2300]\t9000,N\t0,E\t0,V\n5\t101.0\t[2300]\n",
        );
        assert_eq!(m.records.len(), 2);
        assert!(
            m.warnings
                .iter()
                .any(|w| w.contains("2 lines could not be read"))
        );
        assert!(m.warnings.iter().any(|w| w.contains("above 5 m")));
        assert!(
            m.warnings
                .iter()
                .any(|w| w.contains("no antenna-to-camera offset"))
        );
        assert!(m.records[1].offset().is_none());
    }

    #[test]
    fn spaces_instead_of_tabs() {
        let line = "7   345600.250000   [2300]     -120,N       45,E      193,V   45.1,Lat   7.2,Lon   512.345,Ellh   0.012000, 0.011000, 0.025000   50,Q";
        let m = parse(line);
        let r = &m.records[0];
        assert_eq!(r.offset_mm, Some([-120.0, 45.0, 193.0]));
        assert_eq!(r.sigma, [0.012, 0.011, 0.025]);
        assert_eq!(r.q, 50);
    }
}
