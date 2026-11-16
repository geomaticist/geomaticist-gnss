//! Navigation data of a Topcon `.tps` log: broadcast ephemerides (`GE` GPS, `QE` QZSS, `NE` GLONASS, `LE`
//! Galileo, `DE` BeiDou), GPS UTC parameters (`UO`) and ionosphere (`IO`), into RINEX 3.05 records.
//!
//! `GE` / `QE`, `NE`, `UO` and `IO` follow the public GREIS layouts (GPS ephemeris with angles in
//! semi-circles; GLONASS ephemeris with positions in km). `LE` and `DE` are laid out differently from GREIS's
//! Galileo and BeiDou messages.

use geomaticist_gnss_core::time::DateTime;

use crate::t04::nav::{NavRecord, sisa_m};
use crate::util::{WEEK_US, le};

/// GPS's value of pi for semi-circles (IS-GPS-200).
#[allow(clippy::approx_constant)] // the interface specification's value, not f64::consts::PI
const GPS_PI: f64 = 3.141_592_653_589_8;
/// GPS / QZSS / BeiDou URA index to metres.
const URA_M: [f64; 16] = [
    2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0,
    6144.0, 6144.0,
];
/// BeiDou week 0 in GPS weeks.
const BDS_WEEK0: i64 = 1356;

fn ura(i: i64) -> f64 {
    usize::try_from(i)
        .ok()
        .and_then(|i| URA_M.get(i).copied())
        .unwrap_or(6144.0)
}

fn i8_at(b: &[u8], o: usize) -> i64 {
    b.get(o).map_or(0, |&x| i64::from(x as i8))
}

fn u8_at(b: &[u8], o: usize) -> i64 {
    b.get(o).map_or(0, |&x| i64::from(x))
}

fn i32_at(b: &[u8], o: usize) -> i64 {
    le::u32(b, o) as u32 as i32 as i64
}

/// The full week of a week number known modulo `modulus`, nearest to `ref_week`.
#[must_use]
pub fn expand_week(w: i64, modulus: i64, ref_week: i64) -> i64 {
    let base = ref_week - ref_week.rem_euclid(modulus) + w.rem_euclid(modulus);
    [base - modulus, base, base + modulus]
        .into_iter()
        .min_by_key(|x| (x - ref_week).abs())
        .unwrap_or(base)
}

/// The week of a time `t` (s of week) that belongs with a reference time of week `tow` in week `week`.
fn week_of(week: i64, tow: i64, t: i64) -> i64 {
    if t - tow > 302_400 {
        week - 1
    } else if t - tow < -302_400 {
        week + 1
    } else {
        week
    }
}

fn gps_time(week: i64, sec: i64) -> DateTime {
    DateTime::from_micros_since_gps0(week * WEEK_US + sec * 1_000_000)
}

/// An ephemeris with the key it is told apart by when it is logged again.
#[derive(Debug, Clone)]
pub struct Eph {
    /// The RINEX record.
    pub rec: NavRecord,
    /// Issue of data (GLONASS: tb), with the Galileo message type.
    pub iod: i64,
}

impl Eph {
    /// Satellite, epoch and issue: two records with the same key are the same ephemeris logged twice.
    #[must_use]
    pub fn key(&self) -> (String, i64, i64) {
        (
            self.rec.sat.clone(),
            self.rec.epoch.micros_since_gps0(),
            self.iod,
        )
    }
}

/// Where and when a navigation message was logged: the GPS week of the log (full) and the GPS time of the
/// last epoch before it, s since 1980-01-06, and GPS - UTC.
#[derive(Debug, Clone, Copy)]
pub struct NavContext {
    /// Full GPS week of the log.
    pub week: i64,
    /// Reception time, GPS seconds since 1980-01-06.
    pub t_gps: f64,
    /// GPS - UTC, s.
    pub leap: i64,
}

/// `GE` (GPS) or `QE` (QZSS) body, without its checksum (122 bytes).
#[must_use]
pub fn gps_like(b: &[u8], qzss: bool, ctx: &NavContext) -> Option<Eph> {
    if b.len() < 122 {
        return None;
    }
    let sv = u8_at(b, 0);
    let (sys, prn) = if qzss { ('J', sv - 192) } else { ('G', sv) };
    if !(1..=99).contains(&prn) {
        return None;
    }
    let tow = le::u32(b, 1);
    let flags = u8_at(b, 5);
    let iodc = le::i16(b, 6);
    let toc = i32_at(b, 8);
    let ura_i = i8_at(b, 12);
    let health = u8_at(b, 13);
    let wn = le::i16(b, 14);
    let f4 = |o| le::f32(b, o);
    let f8 = |o| le::f64(b, o);
    let toe = i32_at(b, 32);
    let iode = le::i16(b, 36);
    let week_tx = expand_week(wn, 1024, ctx.week);
    let week = week_of(week_tx, tow, toe);
    let week_c = week_of(week_tx, tow, toc);
    let fit = if qzss {
        (flags & 1) as f64
    } else if flags & 1 == 0 {
        4.0
    } else {
        6.0
    };
    let vals = vec![
        Some(f4(28)),
        Some(f4(24)),
        Some(f4(20)),
        Some(iode as f64),
        Some(f4(102)),
        Some(f4(86) * GPS_PI),
        Some(f8(54) * GPS_PI),
        Some(f4(106)),
        Some(f8(46)),
        Some(f4(110)),
        Some(f8(38)),
        Some(toe as f64),
        Some(f4(114)),
        Some(f8(62) * GPS_PI),
        Some(f4(118)),
        Some(f8(70) * GPS_PI),
        Some(f4(98)),
        Some(f8(78) * GPS_PI),
        Some(f4(90) * GPS_PI),
        Some(f4(94) * GPS_PI),
        Some(((flags >> 2) & 3) as f64),
        Some(week as f64),
        Some(((flags >> 1) & 1) as f64),
        Some(ura(ura_i)),
        Some(health as f64),
        Some(f4(16)),
        Some(iodc as f64),
        Some((tow + (week_tx - week) * 604_800) as f64),
        Some(fit),
    ];
    Some(Eph {
        rec: NavRecord {
            sys,
            sat: format!("{sys}{prn:02}"),
            epoch: gps_time(week_c, toc),
            vals,
        },
        iod: iode,
    })
}

/// `NE` (GLONASS) body, without its checksum (94 bytes).
#[must_use]
pub fn glonass(b: &[u8], ctx: &NavContext) -> Option<Eph> {
    if b.len() < 94 {
        return None;
    }
    let slot = u8_at(b, 0);
    if !(1..=99).contains(&slot) {
        return None;
    }
    let fcn = i8_at(b, 1);
    let tk = i32_at(b, 4);
    let tb = i32_at(b, 8);
    let health = i64::from(u8_at(b, 12) & 0x07 != 0);
    let age = u8_at(b, 13);
    let f8 = |o| le::f64(b, o);
    let f4 = |o| le::f32(b, o);
    // tb and tk count Moscow time (UTC + 3 h) from the start of the day
    let rx_utc_us = ((ctx.t_gps - ctx.leap as f64) * 1e6).round() as i64;
    let day_us = 86_400_000_000_i64;
    let tod_us = (tb - 10_800).rem_euclid(86_400) * 1_000_000;
    let day0 = rx_utc_us.div_euclid(day_us) * day_us;
    let epoch_us = [day0 - day_us, day0, day0 + day_us]
        .into_iter()
        .map(|d| d + tod_us)
        .min_by_key(|t| (t - rx_utc_us).abs())
        .unwrap_or(day0 + tod_us);
    // message frame time: seconds of the UTC week
    let tk_tod = (tk - 10_800).rem_euclid(86_400);
    let frame_us = {
        let d = epoch_us.div_euclid(day_us) * day_us;
        let mut t = d + tk_tod * 1_000_000;
        if t - epoch_us > day_us / 2 {
            t -= day_us;
        } else if epoch_us - t > day_us / 2 {
            t += day_us;
        }
        t
    };
    let frame_tow = frame_us.rem_euclid(WEEK_US) as f64 / 1e6;
    let vals = vec![
        Some(-f4(71)),
        Some(f4(75)),
        Some(frame_tow),
        Some(f8(15)),
        Some(f4(39)),
        Some(f4(51)),
        Some(health as f64),
        Some(f8(23)),
        Some(f4(43)),
        Some(f4(55)),
        Some(fcn as f64),
        Some(f8(31)),
        Some(f4(47)),
        Some(f4(59)),
        Some(age as f64),
        None,
        Some(0.999_999_999_999e9),
        Some(15.0),
        None,
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'R',
            sat: format!("R{slot:02}"),
            epoch: DateTime::from_micros_since_gps0(epoch_us),
            vals,
        },
        iod: tb,
    })
}

/// `LE` (Galileo) body, without its checksum (129 bytes).
#[must_use]
pub fn galileo(b: &[u8], ctx: &NavContext) -> Option<Eph> {
    if b.len() < 129 {
        return None;
    }
    let prn = u8_at(b, 0);
    if !(1..=99).contains(&prn) {
        return None;
    }
    let tow = le::u32(b, 1);
    let f4 = |o| le::f32(b, o);
    let f8 = |o| le::f64(b, o);
    let toc = i32_at(b, 17);
    let iod = le::i16(b, 37);
    let toe = i32_at(b, 123);
    let sisa = sisa_m(b[127]);
    // 1: I/NAV (E1-B / E5b, clock for E5b/E1), 4: F/NAV (E5a, clock for E5a/E1)
    let source = match b[128] {
        1 => 517.0,
        4 => 258.0,
        _ => return None,
    };
    let week_rx = (ctx.t_gps / 604_800.0).floor() as i64;
    let tow_rx = (ctx.t_gps - week_rx as f64 * 604_800.0) as i64;
    let week_tx = week_of(week_rx, tow_rx, tow);
    let week = week_of(week_tx, tow, toe);
    let week_c = week_of(week_tx, tow, toc);
    let vals = vec![
        Some(f8(21)),
        Some(f4(29)),
        Some(f4(33)),
        Some(iod as f64),
        Some(f4(111)),
        Some(f4(47) * GPS_PI),
        Some(f8(39) * GPS_PI),
        Some(f4(99)),
        Some(f8(51)),
        Some(f4(103)),
        Some(f8(59)),
        Some(toe as f64),
        Some(f4(115)),
        Some(f8(67) * GPS_PI),
        Some(f4(119)),
        Some(f8(75) * GPS_PI),
        Some(f4(107)),
        Some(f8(83) * GPS_PI),
        Some(f4(91) * GPS_PI),
        Some(f4(95) * GPS_PI),
        Some(source),
        Some(week as f64),
        None,
        Some(sisa),
        Some(0.0),
        Some(f4(9)),
        Some(f4(13)),
        Some((tow + (week_tx - week) * 604_800) as f64),
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'E',
            sat: format!("E{prn:02}"),
            epoch: gps_time(week_c, toc),
            vals,
        },
        iod: iod * 16 + i64::from(b[128]),
    })
}

/// `DE` (BeiDou) body, without its checksum (155 bytes). Times are BDT.
#[must_use]
pub fn beidou(b: &[u8]) -> Option<Eph> {
    if b.len() < 155 {
        return None;
    }
    let prn = u8_at(b, 0);
    if !(1..=99).contains(&prn) {
        return None;
    }
    let tow = le::u32(b, 1);
    let week = le::u16(b, 5);
    let f4 = |o| le::f32(b, o);
    let f8 = |o| le::f64(b, o);
    let toc = i32_at(b, 14);
    let toe = i32_at(b, 31);
    let aode = u8_at(b, 30);
    let week_e = week_of(week, tow, toe);
    let week_c = week_of(week, tow, toc);
    let vals = vec![
        Some(f4(18)),
        Some(f4(22)),
        Some(f4(26)),
        Some(aode as f64),
        Some(f4(107)),
        Some(f4(59) * GPS_PI),
        Some(f8(63) * GPS_PI),
        Some(f4(95)),
        Some(f8(43)),
        Some(f4(99)),
        Some(f8(35)),
        Some(toe as f64),
        Some(f4(111)),
        Some(f8(71) * GPS_PI),
        Some(f4(115)),
        Some(f8(83) * GPS_PI),
        Some(f4(103)),
        Some(f8(51) * GPS_PI),
        Some(f4(79) * GPS_PI),
        Some(f4(91) * GPS_PI),
        None,
        Some(week_e as f64),
        None,
        Some(ura(u8_at(b, 7))),
        Some(u8_at(b, 8) as f64),
        Some(f4(9) * 1e-9),
        Some(f4(119) * 1e-9),
        Some((tow + (week - week_e) * 604_800) as f64),
        Some(u8_at(b, 13) as f64),
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'C',
            sat: format!("C{prn:02}"),
            epoch: gps_time(week_c + BDS_WEEK0, toc),
            vals,
        },
        iod: aode,
    })
}

/// The BeiDou ionosphere (Klobuchar alpha, beta) a `DE` body carries.
#[must_use]
pub fn beidou_iono(b: &[u8]) -> Option<([f64; 4], [f64; 4])> {
    if b.len() < 155 {
        return None;
    }
    let mut a = [0.0; 4];
    let mut g = [0.0; 4];
    for k in 0..4 {
        a[k] = le::f32(b, 123 + 4 * k);
        g[k] = le::f32(b, 139 + 4 * k);
    }
    Some((a, g))
}

/// GPS UTC parameters (`UO`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpsUtc {
    /// A0, s.
    pub a0: f64,
    /// A1, s/s.
    pub a1: f64,
    /// Reference time, s of week.
    pub tot: i64,
    /// Reference week (full).
    pub wnt: i64,
    /// GPS - UTC now, s.
    pub dtls: i64,
    /// Day of the leap second event.
    pub dn: i64,
    /// Week of the leap second event (full).
    pub wnlsf: i64,
    /// GPS - UTC after the event, s.
    pub dtlsf: i64,
}

/// `UO` body, without its checksum (23 bytes).
#[must_use]
pub fn gps_utc(b: &[u8], ref_week: i64) -> Option<GpsUtc> {
    if b.len() < 23 {
        return None;
    }
    Some(GpsUtc {
        a0: le::f64(b, 0),
        a1: le::f32(b, 8),
        tot: le::u32(b, 12),
        wnt: expand_week(le::u16(b, 16), 1024, ref_week),
        dtls: i8_at(b, 18),
        dn: u8_at(b, 19),
        wnlsf: expand_week(le::u16(b, 20), 1024, ref_week),
        dtlsf: i8_at(b, 22),
    })
}

/// `IO` body, without its checksum (38 bytes): GPS Klobuchar alpha and beta.
#[must_use]
pub fn gps_iono(b: &[u8]) -> Option<([f64; 4], [f64; 4])> {
    if b.len() < 38 {
        return None;
    }
    let mut a = [0.0; 4];
    let mut g = [0.0; 4];
    for k in 0..4 {
        a[k] = le::f32(b, 6 + 4 * k);
        g[k] = le::f32(b, 22 + 4 * k);
    }
    Some((a, g))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weeks() {
        assert_eq!(expand_week(362, 1024, 2410), 2410);
        assert_eq!(expand_week(1023, 1024, 2048), 2047);
        assert_eq!(expand_week(0, 1024, 2047), 2048);
        assert_eq!(week_of(100, 600_000, 1000), 101);
        assert_eq!(week_of(100, 1000, 600_000), 99);
        assert_eq!(week_of(100, 1000, 7200), 100);
    }
}
