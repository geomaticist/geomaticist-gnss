//! GNSS time: exact integer RINEX ticks, civil dates, GPS seconds, the leap-second table, BDT.
//!
//! Two time scales:
//! - RINEX epochs as an integer number of ticks of 1e-7 s since 1980-01-06 00:00:00 in the file's
//!   own time system ([`TICKS`], [`to_ticks`], [`from_ticks`], parsing and formatting);
//! - engine times as `f64` GPS seconds since 1980-01-06 GPST ([`gps_seconds`], [`leap_seconds`],
//!   [`utc_to_gps_seconds`], [`gps_datetime`], [`leap_at`]).
//!
//! Why not `hifitime`: the engine works in naive date-times (proleptic Gregorian, microseconds,
//! no leap seconds inside the arithmetic) plus integer ticks, and converts with a small leap table
//! that defaults to 13 s before 1999. Doing that exactly (including rounding a float number of
//! seconds to microseconds half to even, and `leap_at` comparing a GPS date against UTC dates) is
//! simpler and safer with ~300 lines of our own than through a TAI-based library.

use std::fmt;

use serde::{Deserialize, Serialize};

/// Ticks per second (RINEX writes epoch seconds with 7 decimals).
pub const TICKS: i64 = 10_000_000;
/// Ticks per day.
pub const DAY: i64 = 86_400 * TICKS;

/// Days from 0001-01-01 (day 1) in the proleptic Gregorian calendar. `None` for an invalid
/// date.
#[must_use]
pub fn ordinal(year: i64, month: i64, day: i64) -> Option<i64> {
    if !(1..=9999).contains(&year) || !(1..=12).contains(&month) || day < 1 {
        return None;
    }
    if day > days_in_month(year, month) {
        return None;
    }
    let y = year - 1;
    let mut days = y * 365 + y / 4 - y / 100 + y / 400;
    for m in 1..month {
        days += days_in_month(year, m);
    }
    Some(days + day)
}

fn is_leap_year(y: i64) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}

fn days_in_month(y: i64, m: i64) -> i64 {
    match m {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if is_leap_year(y) => 29,
        _ => 28,
    }
}

/// Ordinal day (see [`ordinal`]) -> (year, month, day).
#[must_use]
pub fn from_ordinal(n: i64) -> (i64, i64, i64) {
    // Split into 400-, 100-, 4- and 1-year cycles.
    let n0 = n - 1;
    let (n400, n) = (n0.div_euclid(146_097), n0.rem_euclid(146_097));
    let year = n400 * 400 + 1;
    let (n100, n) = (n / 36_524, n % 36_524);
    let (n4, n) = (n / 1461, n % 1461);
    let (n1, n) = (n / 365, n % 365);
    let mut year = year + n100 * 100 + n4 * 4 + n1;
    if n1 == 4 || n100 == 4 {
        year -= 1;
        return (year, 12, 31);
    }
    let mut month = 1;
    let mut rem = n;
    loop {
        let dim = days_in_month(year, month);
        if rem < dim {
            break;
        }
        rem -= dim;
        month += 1;
        if month > 12 {
            month = 1;
            year += 1;
        }
    }
    (year, month, rem + 1)
}

/// Ordinal of 1980-01-06, the GPS epoch.
pub const GPS0_ORDINAL: i64 = 722_820;

/// Ticks of a civil time (the file's own time system); `None` for an
/// invalid date.
#[must_use]
pub fn to_ticks(y: i64, mo: i64, d: i64, h: i64, mi: i64, sec_ticks: i64) -> Option<i64> {
    let days = ordinal(y, mo, d)? - GPS0_ORDINAL;
    Some(days * DAY + (h * 3600 + mi * 60) * TICKS + sec_ticks)
}

/// Civil time with the seconds in ticks (floor division for negative tick counts).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TickTime {
    /// Year.
    pub year: i64,
    /// Month 1..=12.
    pub month: i64,
    /// Day 1..=31.
    pub day: i64,
    /// Hour.
    pub hour: i64,
    /// Minute.
    pub minute: i64,
    /// Seconds of the minute in ticks.
    pub sec_ticks: i64,
}

/// Civil time of a tick count.
#[must_use]
pub fn from_ticks(t: i64) -> TickTime {
    let days = t.div_euclid(DAY);
    let rem = t.rem_euclid(DAY);
    let (year, month, day) = from_ordinal(GPS0_ORDINAL + days);
    let hour = rem / (3600 * TICKS);
    let rem = rem % (3600 * TICKS);
    TickTime {
        year,
        month,
        day,
        hour,
        minute: rem / (60 * TICKS),
        sec_ticks: rem % (60 * TICKS),
    }
}

/// An integer from a text field: surrounding whitespace allowed, optional sign, digits
/// (underscores between digits accepted). `None` when the field is not an integer.
#[must_use]
pub fn py_int(s: &[u8]) -> Option<i64> {
    let s = crate::text::py_strip(s);
    let (neg, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || digits[0] == b'_' || digits[digits.len() - 1] == b'_' {
        return None;
    }
    let mut v: i64 = 0;
    let mut prev_us = false;
    for &b in digits {
        if b == b'_' {
            if prev_us {
                return None;
            }
            prev_us = true;
            continue;
        }
        prev_us = false;
        if !b.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add(i64::from(b - b'0'))?;
    }
    Some(if neg { -v } else { v })
}

/// `'5.0000000'` -> 50000000 ticks, exactly (more than 7 decimals: the 8th rounds half up).
/// `None` when the text is not a number of seconds.
#[must_use]
pub fn parse_seconds(text: &[u8]) -> Option<i64> {
    let s = crate::text::py_strip(text);
    if s.is_empty() {
        return Some(0);
    }
    let neg = s[0] == b'-';
    let mut start = 0;
    while start < s.len() && (s[start] == b'+' || s[start] == b'-') {
        start += 1;
    }
    let s = &s[start..];
    let (whole, frac) = match s.iter().position(|&b| b == b'.') {
        Some(p) => (&s[..p], &s[p + 1..]),
        None => (s, &b""[..]),
    };
    let whole_i = if whole.is_empty() { 0 } else { py_int(whole)? };
    let mut f = frac.to_vec();
    f.extend_from_slice(b"0000000");
    let mut ticks = whole_i * TICKS + py_int(&f[..7])?;
    if f.len() > 14 && f[7] >= b'5' {
        ticks += 1;
    }
    Some(if neg { -ticks } else { ticks })
}

/// `"{whole}.{frac:07d}"` right-justified.
#[must_use]
pub fn sec_str(sec_ticks: i64, width: usize) -> String {
    let whole = sec_ticks.div_euclid(TICKS);
    let frac = sec_ticks.rem_euclid(TICKS);
    format!("{:>width$}", format!("{whole}.{frac:07}"))
}

/// `'2024 09 23 16 32  5.0000000'`.
#[must_use]
pub fn v3_epoch_time(t: i64) -> String {
    let c = from_ticks(t);
    format!(
        "{:4} {:02} {:02} {:02} {:02}{}",
        c.year,
        c.month,
        c.day,
        c.hour,
        c.minute,
        sec_str(c.sec_ticks, 11)
    )
}

/// `' 26  6  3 16 32  5.0000000'`.
#[must_use]
pub fn v2_epoch_time(t: i64) -> String {
    let c = from_ticks(t);
    format!(
        " {:02} {:2} {:2} {:2} {:2}{}",
        c.year.rem_euclid(100),
        c.month,
        c.day,
        c.hour,
        c.minute,
        sec_str(c.sec_ticks, 11)
    )
}

/// TIME OF FIRST / LAST OBS content: 5I6, F13.7, 5X, A3.
#[must_use]
pub fn header_time(t: i64, system: &str) -> String {
    let c = from_ticks(t);
    format!(
        "{:6}{:6}{:6}{:6}{:6}{}     {:<3}",
        c.year,
        c.month,
        c.day,
        c.hour,
        c.minute,
        sec_str(c.sec_ticks, 13),
        system
    )
}

/// `'2024-09-23 16:32:05.0'` for display.
#[must_use]
pub fn fmt(t: i64, decimals: usize) -> String {
    let c = from_ticks(t);
    let whole = c.sec_ticks.div_euclid(TICKS);
    let frac = c.sec_ticks.rem_euclid(TICKS);
    let mut out = format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
        c.year, c.month, c.day, c.hour, c.minute, whole
    );
    if decimals > 0 {
        let f = format!("{frac:07}");
        out.push('.');
        out.push_str(&f[..decimals.min(7)]);
    }
    out
}

/// `'1 h 05 min 00 s'`.
#[must_use]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub fn fmt_duration(ticks: i64) -> String {
    let secs = ticks as f64 / TICKS as f64;
    let r = crate::text::round_half_even(secs) as i64;
    let (h, rem) = (r.div_euclid(3600), r.rem_euclid(3600));
    let (m, s) = (rem.div_euclid(60), rem.rem_euclid(60));
    if h != 0 {
        format!("{h} h {m:02} min {s:02} s")
    } else if m != 0 {
        format!("{m} min {s:02} s")
    } else {
        format!("{} s", crate::text::format_g(secs))
    }
}

/// `'30 s'`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn fmt_interval(ticks: i64) -> String {
    format!("{} s", crate::text::format_g(ticks as f64 / TICKS as f64))
}

/// Day of year of a tick time.
#[must_use]
pub fn doy(t: i64) -> i64 {
    let c = from_ticks(t);
    ordinal(c.year, c.month, c.day).unwrap_or(0) - ordinal(c.year, 1, 1).unwrap_or(0) + 1
}

/// `'2024-09-23 16:40:00'`, `'2024-09-23T16:40:00.2'`, `'2024 09 23 16 40 00'` -> ticks.
pub fn parse_user_time(text: &str) -> Result<i64, String> {
    let s = text.trim().replace('T', " ").replace('/', "-");
    let s = s.replace(['-', ':'], " ");
    let parts: Vec<&str> = s.split_whitespace().collect();
    let bad = || format!("not a date and time: {text:?} (use YYYY-MM-DD hh:mm:ss)");
    if parts.len() < 5 {
        return Err(bad());
    }
    let mut v = [0i64; 5];
    for (i, p) in parts[..5].iter().enumerate() {
        v[i] = py_int(p.as_bytes()).ok_or_else(bad)?;
    }
    let sec = if parts.len() > 5 {
        parse_seconds(parts[5].as_bytes()).ok_or_else(bad)?
    } else {
        0
    };
    to_ticks(v[0], v[1], v[2], v[3], v[4], sec).ok_or_else(bad)
}

// ---- naive date-time with microseconds ---------------------------------------------------------

/// A naive civil date-time with microseconds (proleptic Gregorian, no time zone, no leap
/// seconds).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct DateTime {
    /// Year 1..=9999.
    pub year: i64,
    /// Month.
    pub month: i64,
    /// Day.
    pub day: i64,
    /// Hour.
    pub hour: i64,
    /// Minute.
    pub minute: i64,
    /// Second 0..=59.
    pub second: i64,
    /// Microsecond.
    pub micro: i64,
}

impl DateTime {
    /// A validated date-time; `None` for an invalid date or time.
    #[must_use]
    pub fn new(
        year: i64,
        month: i64,
        day: i64,
        hour: i64,
        minute: i64,
        second: i64,
    ) -> Option<Self> {
        Self::with_micro(year, month, day, hour, minute, second, 0)
    }

    /// With microseconds.
    #[must_use]
    pub fn with_micro(
        year: i64,
        month: i64,
        day: i64,
        hour: i64,
        minute: i64,
        second: i64,
        micro: i64,
    ) -> Option<Self> {
        ordinal(year, month, day)?;
        if !(0..24).contains(&hour)
            || !(0..60).contains(&minute)
            || !(0..60).contains(&second)
            || !(0..1_000_000).contains(&micro)
        {
            return None;
        }
        Some(Self {
            year,
            month,
            day,
            hour,
            minute,
            second,
            micro,
        })
    }

    /// Microseconds since 1980-01-06 00:00:00.
    #[must_use]
    pub fn micros_since_gps0(&self) -> i64 {
        let days = ordinal(self.year, self.month, self.day).unwrap_or(GPS0_ORDINAL) - GPS0_ORDINAL;
        (days * 86_400 + self.hour * 3600 + self.minute * 60 + self.second) * 1_000_000 + self.micro
    }

    /// From microseconds since 1980-01-06.
    #[must_use]
    pub fn from_micros_since_gps0(us: i64) -> Self {
        let days = us.div_euclid(86_400_000_000);
        let rem = us.rem_euclid(86_400_000_000);
        let (year, month, day) = from_ordinal(GPS0_ORDINAL + days);
        let secs = rem / 1_000_000;
        Self {
            year,
            month,
            day,
            hour: secs / 3600,
            minute: secs % 3600 / 60,
            second: secs % 60,
            micro: rem % 1_000_000,
        }
    }

    /// Day of the year (`timetuple().tm_yday`).
    #[must_use]
    pub fn yday(&self) -> i64 {
        ordinal(self.year, self.month, self.day).unwrap_or(0)
            - ordinal(self.year, 1, 1).unwrap_or(0)
            + 1
    }
}

impl fmt::Display for DateTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            self.year, self.month, self.day, self.hour, self.minute, self.second
        )?;
        if self.micro != 0 {
            write!(f, ".{:06}", self.micro)?;
        }
        Ok(())
    }
}

/// Tick count of a date-time (naive = the file's time system).
#[must_use]
pub fn ticks_from_datetime(dt: &DateTime) -> i64 {
    dt.micros_since_gps0() * 10
}

/// Date-time of a tick count (microsecond resolution, floor).
#[must_use]
pub fn datetime_from_ticks(t: i64) -> DateTime {
    DateTime::from_micros_since_gps0(t.div_euclid(10))
}

/// GPS seconds of a GPS-time date-time: microseconds since GPS0, divided by 1e6.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn gps_seconds(dt: &DateTime) -> f64 {
    // One correctly rounded division; i64 -> f64 is exact below 2**53 µs.
    dt.micros_since_gps0() as f64 / 1e6
}

/// Date-time of GPS seconds `t`: GPS0 plus whole seconds plus the fraction rounded to
/// microseconds, half to even.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn gps_datetime(t: f64) -> DateTime {
    let whole = t.trunc();
    let frac = t - whole;
    let us = crate::text::round_half_even(frac * 1e6) as i64;
    DateTime::from_micros_since_gps0(whole as i64 * 1_000_000 + us)
}

/// GPS - UTC leap seconds table: (first UTC day in force, seconds), newest
/// first; before the oldest entry the leap is 13 s.
#[derive(Debug, Clone, PartialEq)]
pub struct LeapTable {
    entries: Vec<(DateTime, i32)>,
    before: i32,
}

impl Default for LeapTable {
    fn default() -> Self {
        let e = |y, m, n| {
            (
                DateTime {
                    year: y,
                    month: m,
                    day: 1,
                    hour: 0,
                    minute: 0,
                    second: 0,
                    micro: 0,
                },
                n,
            )
        };
        Self {
            entries: vec![
                e(2017, 1, 18),
                e(2015, 7, 17),
                e(2012, 7, 16),
                e(2009, 1, 15),
                e(2006, 1, 14),
                e(1999, 1, 13),
            ],
            before: 13,
        }
    }
}

impl LeapTable {
    /// Adds a leap second (first UTC day in force, new GPS - UTC); keeps newest first.
    pub fn push(&mut self, from: DateTime, leap: i32) {
        self.entries.push((from, leap));
        self.entries.sort_by_key(|e| std::cmp::Reverse(e.0));
    }

    /// GPS - UTC (s) at a UTC date-time.
    #[must_use]
    pub fn at(&self, utc: &DateTime) -> i32 {
        self.entries
            .iter()
            .find(|(when, _)| utc >= when)
            .map_or(self.before, |&(_, n)| n)
    }
}

/// GPS - UTC (s) at a UTC date-time with the built-in table.
#[must_use]
pub fn leap_seconds(utc: &DateTime) -> i32 {
    LeapTable::default().at(utc)
}

/// GPS seconds of a UTC date-time.
#[must_use]
pub fn utc_to_gps_seconds(utc: &DateTime) -> f64 {
    gps_seconds(utc) + f64::from(leap_seconds(utc))
}

/// GPS - UTC (s) at GPS seconds `t`. The table is looked up with the
/// GPS date-time (not UTC).
#[must_use]
pub fn leap_at(t_gps: f64) -> i32 {
    leap_seconds(&gps_datetime(t_gps))
}

/// GPS seconds of a BDT time in seconds (BDT = GPST - 14 s).
#[must_use]
pub fn bdt_to_gps(t_bdt: f64) -> f64 {
    t_bdt + crate::consts::BDT_GPST
}

/// (GPS week, seconds of week) of GPS seconds.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn gps_week_tow(t: f64) -> (i64, f64) {
    let w = (t / crate::consts::WEEK).floor();
    (w as i64, t - w * crate::consts::WEEK)
}

/// Engine time (`f64` GPS seconds) of a tick count: `ticks / TICKS` as numpy's
/// `int64 -> float64 / 10_000_000` (`obs.read_obs`).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn ticks_to_seconds(t: i64) -> f64 {
    t as f64 / TICKS as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordinals_match_python() {
        assert_eq!(ordinal(1980, 1, 6), Some(GPS0_ORDINAL));
        assert_eq!(ordinal(1, 1, 1), Some(1));
        assert_eq!(ordinal(2024, 2, 30), None);
        for n in [1, 365, 366, 730_000, GPS0_ORDINAL, 739_000, 3_652_059] {
            let (y, m, d) = from_ordinal(n);
            assert_eq!(ordinal(y, m, d), Some(n), "{n}");
        }
    }

    #[test]
    fn ticks_round_trip() {
        let t = to_ticks(2024, 9, 23, 16, 32, 50_000_000).unwrap();
        assert_eq!(v3_epoch_time(t), "2024 09 23 16 32  5.0000000");
        assert_eq!(v2_epoch_time(t), " 24  9 23 16 32  5.0000000");
        assert_eq!(
            header_time(t, "GPS"),
            "  2024     9    23    16    32    5.0000000     GPS"
        );
        assert_eq!(fmt(t, 1), "2024-09-23 16:32:05.0");
        assert_eq!(fmt(t, 0), "2024-09-23 16:32:05");
        assert_eq!(from_ticks(-1).year, 1980);
        assert_eq!(doy(t), 267);
    }

    #[test]
    fn seconds_parse_like_python() {
        assert_eq!(parse_seconds(b" 5.0000000"), Some(50_000_000));
        assert_eq!(parse_seconds(b"30"), Some(300_000_000));
        assert_eq!(parse_seconds(b""), Some(0));
        assert_eq!(parse_seconds(b"-0.5"), Some(-5_000_000));
        assert_eq!(parse_seconds(b"1.00000005"), Some(10_000_001)); // 8 decimals: the 8th rounds
        assert_eq!(parse_seconds(b"1.00000004"), Some(10_000_000));
        assert_eq!(parse_seconds(b"1.000000050000000"), Some(10_000_001));
        assert_eq!(parse_seconds(b"1.x"), None);
        assert_eq!(py_int(b"  12 "), Some(12));
        assert_eq!(py_int(b"+7"), Some(7));
        assert_eq!(py_int(b""), None);
        assert_eq!(py_int(b"1 2"), None);
        assert_eq!(py_int(b"1_000"), Some(1000));
    }

    #[test]
    fn durations() {
        assert_eq!(fmt_duration(3600 * TICKS + 300 * TICKS), "1 h 05 min 00 s");
        assert_eq!(fmt_duration(90 * TICKS), "1 min 30 s");
        assert_eq!(fmt_duration(TICKS / 2), "0.5 s");
        assert_eq!(fmt_interval(30 * TICKS), "30 s");
        assert_eq!(
            parse_user_time("2024-09-23T16:40:00.2"),
            to_ticks(2024, 9, 23, 16, 40, 2_000_000).ok_or(String::new())
        );
    }

    #[test]
    fn gps_time_and_leaps() {
        let dt = DateTime::new(2024, 9, 23, 0, 0, 0).unwrap();
        assert_eq!(gps_seconds(&dt), 1_411_084_800.0);
        assert_eq!(leap_seconds(&dt), 18);
        assert_eq!(utc_to_gps_seconds(&dt), 1_411_084_818.0);
        assert_eq!(
            leap_seconds(&DateTime::new(1998, 12, 31, 23, 59, 59).unwrap()),
            13
        );
        assert_eq!(
            leap_seconds(&DateTime::new(2016, 12, 31, 0, 0, 0).unwrap()),
            17
        );
        let back = gps_datetime(1_411_084_800.25);
        assert_eq!(
            (back.year, back.month, back.day, back.micro),
            (2024, 9, 23, 250_000)
        );
        assert_eq!(leap_at(1_411_084_800.0), 18);
        let mut t = LeapTable::default();
        t.push(DateTime::new(2030, 1, 1, 0, 0, 0).unwrap(), 19);
        assert_eq!(t.at(&DateTime::new(2031, 1, 1, 0, 0, 0).unwrap()), 19);
        assert_eq!(gps_week_tow(1_411_084_800.0), (2333, 86_400.0));
        assert_eq!(ticks_to_seconds(14_110_848_000_000_000), 1_411_084_800.0);
    }
}
