//! Broadcast navigation records -> RINEX 3 navigation values.
//!
//! T04 inner record types (all big-endian): 0x15 GPS LNAV ephemeris (sub = PRN), 187 bytes; 0x1c sub 0x01
//! SBAS, 100 bytes; sub 0x02 GLONASS, 137 bytes; sub 0x04 Galileo, 186 bytes (byte 4: 0 = I/NAV E1-B,
//! 1 = I/NAV E5b, 2 = F/NAV); sub 0x05 QZSS, 179 bytes; sub 0x0b BeiDou, 179 bytes. One RINEX record is
//! written per T04 record, in file order. Angles and the four harmonic corrections Cuc/Cus/Cic/Cis are
//! stored in semicircles; RINEX gets x*pi.

use geomaticist_gnss_core::time::DateTime;

use super::T04Error;
use super::container::{Inner, Records};
use crate::decimal::Dec;
use crate::util::{be, date_us, gps_week_s};

/// pi in full double precision (not the GPS-ICD 3.1415926535898).
const PI: f64 = std::f64::consts::PI;
const URA_M: [f64; 16] = [
    2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0,
    6144.0,
];

/// One ephemeris: system letter, RINEX satellite, epoch, values in RINEX 3.05 order (`None` = blank field).
#[derive(Debug, Clone)]
pub struct NavRecord {
    /// System letter.
    pub sys: char,
    /// RINEX satellite ("G05").
    pub sat: String,
    /// Epoch (toc; GLONASS: UTC; BDS: BDT).
    pub epoch: DateTime,
    /// Values.
    pub vals: Vec<Option<f64>>,
}

/// A record logged before the week rollover can already carry the next week's toe (00:00 Sunday) while its
/// week number is still the week it was sent in. Written as: the toe's week, and
/// for GPS the transmission time counted from that week (negative); for QZSS the transmission time stays in the
/// week it was sent (assumed: only one such QZSS record so far).
#[must_use]
pub fn roll_week(week: i64, ttm: i64, toe: i64, shift_ttm: bool) -> (i64, i64) {
    if toe - ttm < -302_400 {
        (week + 1, if shift_ttm { ttm - 604_800 } else { ttm })
    } else {
        (week, ttm)
    }
}

/// The 19 doubles common to GPS/QZSS/BDS (base = offset of TGD), angles in radians.
struct Kepler {
    tgd: f64,
    af2: f64,
    af1: f64,
    af0: f64,
    crs: f64,
    dn: f64,
    m0: f64,
    cuc: f64,
    e: f64,
    cus: f64,
    sqrta: f64,
    cic: f64,
    om0: f64,
    cis: f64,
    i0: f64,
    crc: f64,
    w: f64,
    omd: f64,
    idot: f64,
}

fn kepler(b: &[u8], base: usize) -> Kepler {
    let g = |k: usize| be::f64(b, base + 8 * k);
    Kepler {
        tgd: g(0),
        af2: g(1),
        af1: g(2),
        af0: g(3),
        crs: g(4),
        dn: g(5) * PI,
        m0: g(6) * PI,
        cuc: g(7) * PI,
        e: g(8),
        cus: g(9) * PI,
        sqrta: g(10),
        cic: g(11) * PI,
        om0: g(12) * PI,
        cis: g(13) * PI,
        i0: g(14) * PI,
        crc: g(15),
        w: g(16) * PI,
        omd: g(17) * PI,
        idot: g(18) * PI,
    }
}

fn ura_meters(idx: u8) -> f64 {
    URA_M.get(usize::from(idx)).copied().unwrap_or(6144.0)
}

#[allow(clippy::cast_precision_loss)]
fn f(v: i64) -> Option<f64> {
    Some(v as f64)
}

fn some(v: &[f64]) -> Vec<Option<f64>> {
    v.iter().copied().map(Some).collect()
}

/// GPS LNAV (0x15).
#[allow(clippy::cast_precision_loss)]
fn decode_gps(b: &[u8]) -> NavRecord {
    let k = kepler(b, 21);
    let (week, iodc, iode, ttm, toc, toe) = (
        be::u16(b, 3),
        be::u16(b, 5),
        be::u16(b, 7),
        be::u32(b, 9),
        be::u32(b, 13),
        be::u32(b, 17),
    );
    let (week, ttm) = roll_week(week, ttm, toe, true);
    let flags = b[176]; // bit0 fit flag, bits1-2 L2 codes, bit4 L2P flag
    let (l2codes, l2p, fitflag) = ((flags >> 1) & 3, (flags >> 4) & 1, flags & 1);
    let ura = ura_meters(b[175] >> 3); // URA index in the top 5 bits
    let mut vals = some(&[k.af0, k.af1, k.af2]);
    vals.extend([f(iode), Some(k.crs), Some(k.dn), Some(k.m0)]);
    vals.extend(some(&[k.cuc, k.e, k.cus, k.sqrta]));
    vals.extend([f(toe), Some(k.cic), Some(k.om0), Some(k.cis)]);
    vals.extend(some(&[k.i0, k.crc, k.w, k.omd]));
    vals.extend([
        Some(k.idot),
        Some(f64::from(l2codes)),
        f(week),
        Some(f64::from(l2p)),
    ]);
    vals.extend([Some(ura), Some(f64::from(b[173])), Some(k.tgd), f(iodc)]);
    // health byte 173; fit flag 1 -> 6 h (assumed, not met yet)
    vals.extend([f(ttm), Some(if fitflag == 0 { 4.0 } else { 6.0 })]);
    NavRecord {
        sys: 'G',
        sat: format!("G{:02}", b[2]),
        epoch: gps_week_s(week, toc),
        vals,
    }
}

/// QZSS (0x1c sub 5).
fn decode_qzss(b: &[u8]) -> NavRecord {
    let k = kepler(b, 23);
    let (week, iodc, iode, ttm, toc, toe) = (
        be::u16(b, 5),
        be::u16(b, 7),
        be::u16(b, 9),
        be::u32(b, 11),
        be::u32(b, 15),
        be::u32(b, 19),
    );
    let (week, ttm) = roll_week(week, ttm, toe, false);
    // byte 178: bit 0 fit flag, bits 1-2 codes on L2 as in GPS, but bit 4 is health and the L2P flag is byte 176
    // (fw 6.x has both 1; fw 5.41 has 176 = 1 and bit 4 = 0, written as L2P 1, health 0)
    let mut vals = some(&[k.af0, k.af1, k.af2]);
    vals.extend([f(iode), Some(k.crs), Some(k.dn), Some(k.m0)]);
    vals.extend(some(&[k.cuc, k.e, k.cus, k.sqrta]));
    vals.extend([f(toe), Some(k.cic), Some(k.om0), Some(k.cis)]);
    vals.extend(some(&[k.i0, k.crc, k.w, k.omd]));
    vals.extend([
        Some(k.idot),
        Some(f64::from((b[178] >> 1) & 3)),
        f(week),
        Some(f64::from(b[176] & 1)),
    ]);
    vals.extend([
        Some(ura_meters(b[177] >> 3)),
        Some(f64::from((b[178] >> 4) & 1)),
        Some(k.tgd),
        f(iodc),
    ]);
    vals.extend([f(ttm), Some(f64::from(b[178] & 1))]); // QZSS writes the fit flag
    NavRecord {
        sys: 'J',
        sat: format!("J{:02}", i64::from(b[3]) - 192),
        epoch: gps_week_s(week, toc),
        vals,
    }
}

/// BeiDou (0x1c sub 0x0b).
fn decode_bds(b: &[u8]) -> NavRecord {
    let k = kepler(b, 23);
    let (week, aodc, aode) = (be::u16(b, 5), be::u16(b, 7), be::u16(b, 9));
    let (ttm, toc, toe) = (
        be::u32(b, 11) - 14,
        be::u32(b, 15) - 14,
        be::u32(b, 19) - 14,
    ); // GPS s -> BDT s
    let t10 = ((i64::from(b[175]) << 8) | i64::from(b[176])) & 0x3FF; // TGD2, 10-bit signed, 0.1 ns
    #[allow(clippy::cast_precision_loss)]
    let tgd2 = (if t10 & 0x200 != 0 { t10 - 1024 } else { t10 }) as f64 * 1e-10;
    let mut vals = some(&[k.af0, k.af1, k.af2]);
    vals.extend([f(aode), Some(k.crs), Some(k.dn), Some(k.m0)]);
    vals.extend(some(&[k.cuc, k.e, k.cus, k.sqrta]));
    vals.extend([f(toe), Some(k.cic), Some(k.om0), Some(k.cis)]);
    vals.extend(some(&[k.i0, k.crc, k.w, k.omd]));
    vals.extend([Some(k.idot), None, f(week - 1356), None]);
    // URA 2.0 assumed (the only value so far; low nibble of b[178] = 3)
    vals.extend([
        Some(2.0),
        Some(f64::from((b[178] >> 4) & 1)),
        Some(k.tgd),
        Some(tgd2),
    ]);
    vals.extend([f(ttm), f(aodc)]);
    NavRecord {
        sys: 'C',
        sat: format!("C{:02}", b[3]),
        // gps_time(week, toc + 14) - 14 s: the BDT epoch, printed as such
        epoch: gps_week_s(week, toc),
        vals,
    }
}

/// Galileo signal-in-space accuracy, m.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn sisa_m(n: u8) -> f64 {
    let n = i64::from(n);
    if n < 50 {
        n as f64 * 0.01
    } else if n < 75 {
        0.5 + (n - 50) as f64 * 0.02
    } else if n < 100 {
        1.0 + (n - 75) as f64 * 0.04
    } else if n < 126 {
        2.0 + (n - 100) as f64 * 0.16
    } else {
        -1.0
    }
}

/// Galileo (0x1c sub 4).
fn decode_gal(b: &[u8]) -> NavRecord {
    let kind = b[4];
    let (week, ttm, iod, toe) = (be::u16(b, 5), be::u32(b, 7), be::u16(b, 11), be::u32(b, 13));
    let g = |o: usize| be::f64(b, o);
    let toc = be::u32(b, 140);
    let (bgd_a, bgd_b) = if kind == 2 {
        (g(168), 0.0)
    } else {
        (g(177), g(168))
    };
    let source = match kind {
        0 => 513.0,
        1 => 516.0,
        2 => 258.0,
        _ => 0.0,
    };
    let mut vals = some(&[g(144), g(152), g(160)]);
    vals.extend([f(iod), Some(g(17)), Some(g(25) * PI), Some(g(33) * PI)]);
    vals.extend(some(&[g(41) * PI, g(49), g(57) * PI, g(65)]));
    vals.extend([f(toe), Some(g(73) * PI), Some(g(81) * PI), Some(g(89) * PI)]);
    vals.extend(some(&[g(97) * PI, g(105), g(113) * PI, g(121) * PI]));
    vals.extend(some(&[g(129) * PI, source]));
    vals.extend([f(week), Some(0.0)]);
    vals.extend([
        Some(geomaticist_gnss_core::text::round_to(sisa_m(b[137]), 2)),
        Some(0.0),
        Some(bgd_a),
        Some(bgd_b),
    ]);
    vals.extend([f(ttm), Some(0.0)]); // health 0 assumed (the only value so far)
    NavRecord {
        sys: 'E',
        sat: format!("E{:02}", b[3]),
        epoch: gps_week_s(week, toc),
        vals,
    }
}

/// GLONASS (0x1c sub 2).
fn decode_glo(b: &[u8]) -> Result<NavRecord, T04Error> {
    let leap = i64::from(b[13]);
    let nt = be::u16(b, 10); // day within the 4-year interval
    let tb_idx = i64::from(b[12]);
    let week = be::u16(b, 4);
    let tk = be::u32(b, 6) - leap; // RINEX: message frame time, UTC seconds of week
    let r = gps_week_s(week, 0);
    let mut y4 = r.year - r.year.rem_euclid(4); // start of the 4-year interval (leap year)
    let ref_us = r.micros_since_gps0();
    let bad = || T04Error::Decode("GLONASS ephemeris date out of range".into());
    if date_us(y4, 1, 1).ok_or_else(bad)? + (nt - 1) * 86_400_000_000 > ref_us + 8 * 86_400_000_000
    {
        y4 -= 4;
    }
    let msk =
        date_us(y4, 1, 1).ok_or_else(bad)? + (nt - 1) * 86_400_000_000 + 15 * tb_idx * 60_000_000;
    let epoch = DateTime::from_micros_since_gps0(msk - 3 * 3_600_000_000); // UTC
    let g = |o: usize| be::f64(b, o);
    let k = i64::from(b[21] as i8);
    // byte 22 = 3 on the only unhealthy records so far; written as health 1 and flags 3 (else 0, 7)
    let (health, flags) = if b[22] != 0 { (1.0, 3.0) } else { (0.0, 7.0) };
    let mut vals = some(&[-g(105), g(113)]);
    vals.push(f(tk));
    vals.extend(some(&[
        g(25) / 1000.0,
        g(33) / 1000.0,
        g(41) / 1000.0,
        health,
    ]));
    vals.extend(some(&[g(49) / 1000.0, g(57) / 1000.0, g(65) / 1000.0]));
    vals.push(f(k));
    vals.extend(some(&[g(73) / 1000.0, g(81) / 1000.0, g(89) / 1000.0, 0.0]));
    vals.extend(some(&[
        f64::from(b[14] | 0x80),
        g(129),
        f64::from(b[24]),
        flags,
    ])); // age 0 assumed
    Ok(NavRecord {
        sys: 'R',
        sat: format!("R{:02}", b[3]),
        epoch,
        vals,
    })
}

/// SBAS (0x1c sub 1).
fn decode_sbas(b: &[u8]) -> NavRecord {
    let (week, t) = (be::u16(b, 4), be::u32(b, 6));
    let g = |o: usize| be::f64(b, o);
    let mut vals = some(&[g(84), g(92)]);
    vals.push(f(t));
    vals.extend(some(&[g(12) / 1000.0, g(20) / 1000.0, g(28) / 1000.0, 0.0]));
    vals.extend(some(&[g(36) / 1000.0, g(44) / 1000.0, g(52) / 1000.0, 0.0]));
    // health, URA, IODN 0 assumed (the only values so far)
    vals.extend(some(&[g(60) / 1000.0, g(68) / 1000.0, g(76) / 1000.0, 0.0]));
    NavRecord {
        sys: 'S',
        sat: format!("S{:02}", i64::from(b[3]) - 100),
        epoch: gps_week_s(week, t),
        vals,
    }
}

/// Every ephemeris record, in file order.
pub fn iter_nav(inner: &Inner) -> Result<Vec<NavRecord>, T04Error> {
    let mut out = Vec::new();
    for r in Records::new(inner) {
        let r = r?;
        let b = r.raw;
        if r.rtype == 0x15 && b.len() == 187 {
            out.push(decode_gps(b));
        } else if r.rtype == 0x1C && b.len() > 2 {
            match (b[2], b.len()) {
                (0x01, 100) => out.push(decode_sbas(b)),
                (0x02, 137) => out.push(decode_glo(b)?),
                (0x04, 186) => out.push(decode_gal(b)),
                (0x05, 179) => out.push(decode_qzss(b)),
                (0x0B, 179) => out.push(decode_bds(b)),
                _ => {}
            }
        }
    }
    Ok(out)
}

/// RINEX D19.12: the double is first rounded to 15 significant digits (the
/// .NET default), then to 13 (ties away from zero); no negative zero. `None` -> 19 blanks.
#[must_use]
pub fn fmt_d19(v: Option<f64>) -> String {
    let Some(v) = v else {
        return " ".repeat(19);
    };
    if v == 0.0 {
        return " 0.000000000000e+00".into();
    }
    let d = Dec::exact(v).round_sig(15);
    let mut exp = d.adjusted();
    let mut q = d.scaleb(-exp).round_at(-12);
    if q.adjusted() >= 1 {
        exp += 1;
        q = d.scaleb(-exp).round_at(-12);
    }
    let sign = if q.neg { '-' } else { ' ' };
    format!("{sign}{}e{exp:+03}", q.abs_fixed())
}

/// D<width>.<mant> with the same rounding (15 significant digits, then ties away), right-aligned
///.
#[must_use]
pub fn fmt_e(v: f64, mant: i32, width: usize) -> String {
    if v == 0.0 {
        let body = format!("0.{}e+00", "0".repeat(usize::try_from(mant).unwrap_or(0)));
        return format!("{body:>width$}");
    }
    let d = Dec::exact(v).round_sig(15);
    let mut exp = d.adjusted();
    let mut q = d.scaleb(-exp).round_at(-mant);
    if q.adjusted() >= 1 {
        exp += 1;
        q = d.scaleb(-exp).round_at(-mant);
    }
    let body = format!(
        "{}{}e{exp:+03}",
        if q.neg { "-" } else { "" },
        q.abs_fixed()
    );
    format!("{body:>width$}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d19_format() {
        assert_eq!(fmt_d19(Some(0.0)), " 0.000000000000e+00");
        assert_eq!(fmt_d19(Some(-1.5)), "-1.500000000000e+00");
        assert_eq!(fmt_d19(Some(4.0)), " 4.000000000000e+00");
        assert_eq!(fmt_d19(Some(9.999_999_999_999_9)), " 1.000000000000e+01");
        assert_eq!(fmt_d19(None), " ".repeat(19));
        assert_eq!(fmt_e(1.2345e-8, 4, 12), "  1.2345e-08");
        assert_eq!(fmt_e(0.0, 4, 12), "  0.0000e+00");
    }

    #[test]
    fn week_rollover() {
        assert_eq!(roll_week(2300, 604_000, 0, true), (2301, -800));
        assert_eq!(roll_week(2300, 604_000, 0, false), (2301, 604_000));
        assert_eq!(roll_week(2300, 1000, 7200, true), (2300, 1000));
    }
}
