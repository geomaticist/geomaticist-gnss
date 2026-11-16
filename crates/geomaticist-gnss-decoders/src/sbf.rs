//! Septentrio SBF logs: MeasEpoch into the engine's observations; GPSNav, GLONav, GALNav, BDSNav and QZSNav
//! into RINEX 3 navigation; ReceiverSetup and the PVT blocks into the header facts.
//!
//! A clean-room decoder written from the published SBF reference.
//!
//! Every block carries the GPS week, so no outside date is needed. Only the main antenna's measurements
//! are kept.
//!
//! Known limits:
//! - Tested on synthetic blocks only: no openly licensed SBF sample log was found.
//! - GPS L1P(Y) carrier phase and Doppler are written (`L1W`, `D1W`).
//! - A phase whose half-cycle flag is set is written with LLI bit 1.
//! - GPS, QZSS and BeiDou SV accuracy is the upper bound of the URA index range (2.4 m for index 0),
//!   as in the other decoders here.
//! - MeasExtra is not read, so there is no multipath or code-variance information.
//! - Left out, and counted in a warning: QZSS L6 (no RINEX 3 code here) and L1C/B, and any signal
//!   number not in `sig_tail`. Left out silently: L-band (MSS) satellites, SVIDs outside the table, and
//!   the second and third antennas.
//! - SBAS and NavIC observations are kept, but their navigation blocks (GEONav, NavICRaw and the like)
//!   are not decoded; nor are GPS CNAV, Galileo FNAV/INAV raw pages or BeiDou CNAV.
//! - GLONav carries no t_k, so the GLONASS record's message frame time is the start of the 30-second
//!   frame in which the block was received (UTC seconds of the week).
//! - GLOTime is not read, so the navigation header has no GLUT/GLGP time-system corrections.
//! - The receiver serial and firmware go only into the info JSON and are not written to the RINEX
//!   `REC # / TYPE / VERS` line: the shared observation writer takes the receiver type alone.
//! - The observation header's program name is the one `write_obs` sets.
//! - Galileo E6 is written as 6B or 6C for the whole epoch from the block's E6B flag, not per signal.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use geomaticist_gnss_core::consts::{CLIGHT, WEEK};
use geomaticist_gnss_core::obs::norm_antenna;
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::time::{DateTime, gps_datetime, leap_at};
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};

use crate::epochs::{Epochs, stem_and_name, write_obs};
use crate::t04::nav::{NavRecord, sisa_m};
use crate::topcon::nav::{Eph, GpsUtc, expand_week};
use crate::topcon::rinex::{NavHeader, nav_text};
use crate::util::WEEK_US;
use crate::{Converted, DecodeError};

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
/// The do-not-use value of SBF floats.
const F_DNU: f64 = -2e10;

// block numbers
const MEAS_EPOCH: u16 = 4027;
const GPS_NAV: u16 = 5891;
const GPS_ION: u16 = 5893;
const GPS_UTC: u16 = 5894;
const GLO_NAV: u16 = 4004;
const GAL_NAV: u16 = 4002;
const BDS_NAV: u16 = 4081;
const BDS_ION: u16 = 4120;
const QZS_NAV: u16 = 4095;
const PVT_CART: u16 = 4006;
const PVT_GEOD: u16 = 4007;
const RX_TIME: u16 = 5914;
const RX_SETUP: u16 = 5902;

/// CRC-16-CCITT (polynomial 0x1021, seed 0) as SBF uses it.
fn crc16(d: &[u8]) -> u16 {
    let mut crc = 0u16;
    for &b in d {
        crc ^= u16::from(b) << 8;
        for _ in 0..8 {
            crc = if crc & 0x8000 == 0 {
                crc << 1
            } else {
                (crc << 1) ^ 0x1021
            };
        }
    }
    crc
}

/// One SBF block: number, revision and the whole block from its sync bytes.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Block<'a> {
    pub num: u16,
    pub rev: u8,
    pub b: &'a [u8],
}

/// The SBF blocks of a stream with a good CRC, and the number of bad ones.
pub(crate) fn blocks(data: &[u8]) -> (Vec<Block<'_>>, usize) {
    let mut out = Vec::new();
    let mut bad = 0;
    let mut i = 0;
    while i + 8 <= data.len() {
        if data[i] != b'$' || data[i + 1] != b'@' {
            i += 1;
            continue;
        }
        let len = usize::from(u2(data, i + 6));
        let end = i + len;
        if len < 14 || len % 4 != 0 || end > data.len() {
            i += 1;
            continue;
        }
        if crc16(&data[i + 4..end]) == u2(data, i + 2) {
            let id = u2(data, i + 4);
            out.push(Block {
                num: id & 0x1FFF,
                rev: (id >> 13) as u8,
                b: &data[i..end],
            });
            i = end;
        } else {
            bad += 1;
            i += 1;
        }
    }
    (out, bad)
}

fn u2(d: &[u8], o: usize) -> u16 {
    u16::from_le_bytes([d[o], d[o + 1]])
}
fn u4(d: &[u8], o: usize) -> u32 {
    u32::from_le_bytes([d[o], d[o + 1], d[o + 2], d[o + 3]])
}
fn i4(d: &[u8], o: usize) -> i32 {
    u4(d, o) as i32
}
fn f4(d: &[u8], o: usize) -> f64 {
    f64::from(f32::from_bits(u4(d, o)))
}
fn f8(d: &[u8], o: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[o..o + 8]);
    f64::from_le_bytes(b)
}
/// An f4 that may be do-not-use.
fn f4_opt(d: &[u8], o: usize) -> Option<f64> {
    let v = f4(d, o);
    (v != F_DNU as f32 as f64).then_some(v)
}
/// Zero-padded text.
fn text(d: &[u8], o: usize, n: usize) -> String {
    let s = d.get(o..o + n).unwrap_or_default();
    let end = s.iter().position(|&c| c == 0).unwrap_or(s.len());
    s[..end]
        .iter()
        .map(|&c| char::from(c))
        .collect::<String>()
        .trim()
        .to_owned()
}

/// The block's time: full GPS week and time of week (ms), or `None` when do-not-use.
fn block_time(b: &[u8]) -> Option<(i64, i64)> {
    let tow = u4(b, 8);
    let wn = u2(b, 12);
    (tow != u32::MAX && wn != u16::MAX).then(|| (i64::from(wn), i64::from(tow)))
}

fn ura(i: u8) -> f64 {
    URA_M.get(usize::from(i)).copied().unwrap_or(6144.0)
}

fn gps_time(week: i64, sec: i64) -> DateTime {
    DateTime::from_micros_since_gps0(week * WEEK_US + sec * 1_000_000)
}

/// The satellite of an SBF SVID; `None` for L-band, unknown GLONASS slots and unused numbers.
fn svid_sat(svid: u8) -> Option<Sat> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Navic, Qzss, Sbas};
    let (sys, prn) = match svid {
        1..=37 => (Gps, svid),
        38..=61 => (Glonass, svid - 37),
        63..=68 => (Glonass, svid - 38),
        71..=106 => (Galileo, svid - 70),
        120..=140 => (Sbas, svid - 100),
        141..=180 => (BeiDou, svid - 140),
        181..=187 => (Qzss, svid - 180),
        191..=197 => (Navic, svid - 190),
        198..=215 => (Sbas, svid - 57),
        216..=222 => (Navic, svid - 208),
        223..=245 => (BeiDou, svid - 182),
        _ => return None,
    };
    Some(Sat::new(sys, prn))
}

/// The RINEX 3 code tail of an SBF signal number on a satellite of `sys`; `None` when unknown here.
fn sig_tail(sig: u8, sys: GnssSystem, e6b: bool) -> Option<&'static str> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Navic, Qzss, Sbas};
    Some(match (sig, sys) {
        (0, Gps) | (6, Qzss) | (8, Glonass) | (17, Galileo) | (24, Sbas) => "1C",
        (1, Gps) => "1W",
        (2, Gps) => "2W",
        (3, Gps) | (7, Qzss) => "2L",
        (4, Gps) | (26, Qzss) => "5Q",
        (5, Gps) | (32, Qzss) => "1L",
        (33, Qzss) => "1Z",
        (39, Qzss) | (14, BeiDou) => "5P",
        (9, Glonass) => "1P",
        (10, Glonass) => "2P",
        (11, Glonass) => "2C",
        (12, Glonass) => "3Q",
        (13, BeiDou) => "1P",
        (28, BeiDou) => "2I",
        (29, BeiDou) => "7I",
        (30, BeiDou) => "6I",
        (34, BeiDou) => "7D",
        (15, Navic) => "5A",
        (19, Galileo) if e6b => "6B",
        (19, Galileo) => "6C",
        (20, Galileo) => "5Q",
        (21, Galileo) => "7Q",
        (22, Galileo) => "8Q",
        (25, Sbas) => "5I",
        _ => return None,
    })
}

/// C/N0 (dB-Hz) of a raw `CN0` field for signal number `sig`.
fn cn0(raw: u8, sig: u8) -> f64 {
    if raw == 255 {
        f64::NAN
    } else if sig == 1 || sig == 2 {
        f64::from(raw) * 0.25
    } else {
        f64::from(raw) * 0.25 + 10.0
    }
}

/// The signal number from a sub-block's type byte and `ObsInfo`, and the GLONASS frequency number when the
/// `ObsInfo` carries one.
fn signal_of(typ: u8, info: u8) -> (u8, Option<i32>) {
    let lo = typ & 31;
    if lo == 31 {
        ((info >> 3) + 32, None)
    } else if (8..=11).contains(&lo) {
        (lo, Some(i32::from(info >> 3) - 8))
    } else {
        (lo, None)
    }
}

/// A sign-extended field of `bits` bits.
fn signed(v: u8, bits: u32) -> i64 {
    let v = i64::from(v);
    let half = 1i64 << (bits - 1);
    if v >= half { v - (1 << bits) } else { v }
}

/// Loss-of-lock state per satellite and signal: last epoch with a phase (s) and its lock time (s).
type LockMap = BTreeMap<(Sat, &'static str), (f64, u16)>;

/// One measurement of a sub-block.
struct Meas {
    sat: Sat,
    tail: &'static str,
    code: f64,
    phase: f64,
    doppler: f64,
    snr: f64,
    lock: Option<u16>,
    half: bool,
}

/// Counts of what [`decode`] left out.
#[derive(Debug, Default, Clone, Copy)]
struct Dropped {
    unknown_signals: usize,
    scrambled: usize,
}

/// The measurements of one MeasEpoch block; GLONASS frequency numbers go into `fcn`.
#[allow(clippy::similar_names)]
fn meas_epoch(b: &[u8], fcn: &mut BTreeMap<Sat, i32>, dropped: &mut Dropped) -> Vec<Meas> {
    let mut out = Vec::new();
    if b.len() < 20 {
        return out;
    }
    let n1 = usize::from(b[14]);
    let sb1 = usize::from(b[15]);
    let sb2 = usize::from(b[16]);
    let flags = b[17];
    if flags & 0x80 != 0 {
        dropped.scrambled += 1;
        return out;
    }
    let e6b = flags & 0x40 != 0;
    if sb1 < 20 {
        return out;
    }
    let mut o = 20;
    for _ in 0..n1 {
        if o + sb1 > b.len() {
            break;
        }
        let s = &b[o..o + sb1];
        let n2 = usize::from(s[19]);
        o += sb1;
        if n2 > 0 && (sb2 < 12 || o + n2 * sb2 > b.len()) {
            break;
        }
        let subs = &b[o..o + n2 * sb2];
        o += n2 * sb2;
        // the main antenna only
        if s[1] >> 5 != 0 {
            continue;
        }
        let Some(sat) = svid_sat(s[2]) else {
            continue;
        };
        let (sig1, k1) = signal_of(s[1], s[18]);
        if let Some(k) = k1 {
            fcn.insert(sat, k);
        }
        let k = fcn.get(&sat).copied();
        let code_raw = u64::from(s[3] & 15) << 32 | u64::from(u4(s, 4));
        let pr1 = if code_raw == 0 {
            f64::NAN
        } else {
            code_raw as f64 * 0.001
        };
        let d_raw = i4(s, 8);
        let d1 = if d_raw == i32::MIN {
            f64::NAN
        } else {
            f64::from(d_raw) * 1e-4
        };
        let f1 =
            sig_tail(sig1, sat.sys, e6b).and_then(|t| freq(sat, t.as_bytes()[0] - b'0', k).ok());
        if let Some(tail) = sig_tail(sig1, sat.sys, e6b) {
            let carr_msb = s[14] as i8;
            let carr_lsb = u2(s, 12);
            let phase = match f1 {
                Some(f) if pr1.is_finite() && !(carr_msb == -128 && carr_lsb == 0) => {
                    pr1 * f / CLIGHT
                        + (i64::from(carr_msb) * 65536 + i64::from(carr_lsb)) as f64 * 0.001
                }
                _ => f64::NAN,
            };
            let lock = u2(s, 16);
            out.push(Meas {
                sat,
                tail,
                code: pr1,
                phase,
                doppler: d1,
                snr: cn0(s[15], sig1),
                lock: (lock != u16::MAX).then_some(lock),
                half: s[18] & 4 != 0,
            });
        } else {
            dropped.unknown_signals += 1;
        }
        for q in subs.chunks_exact(sb2) {
            if q[0] >> 5 != 0 {
                continue;
            }
            let (sig2, k2) = signal_of(q[0], q[5]);
            if let Some(k) = k2 {
                fcn.insert(sat, k);
            }
            let Some(tail) = sig_tail(sig2, sat.sys, e6b) else {
                dropped.unknown_signals += 1;
                continue;
            };
            let f2 = freq(sat, tail.as_bytes()[0] - b'0', fcn.get(&sat).copied()).ok();
            let code_msb = signed(q[3] & 7, 3);
            let code_lsb = i64::from(u2(q, 6));
            let pr2 = if pr1.is_finite() && !(code_msb == -4 && code_lsb == 0) {
                pr1 + (code_msb * 65536 + code_lsb) as f64 * 0.001
            } else {
                f64::NAN
            };
            let carr_msb = q[4] as i8;
            let carr_lsb = u2(q, 8);
            let phase = match f2 {
                Some(f) if pr2.is_finite() && !(carr_msb == -128 && carr_lsb == 0) => {
                    pr2 * f / CLIGHT
                        + (i64::from(carr_msb) * 65536 + i64::from(carr_lsb)) as f64 * 0.001
                }
                _ => f64::NAN,
            };
            let dop_msb = signed(q[3] >> 3, 5);
            let dop_lsb = i64::from(u2(q, 10));
            let doppler = match (f1, f2) {
                (Some(fa), Some(fb)) if d1.is_finite() && !(dop_msb == -16 && dop_lsb == 0) => {
                    d1 * fb / fa + (dop_msb * 65536 + dop_lsb) as f64 * 1e-4
                }
                _ => f64::NAN,
            };
            out.push(Meas {
                sat,
                tail,
                code: pr2,
                phase,
                doppler,
                snr: cn0(q[2], sig2),
                lock: (q[1] != 255).then_some(u16::from(q[1])),
                half: q[5] & 4 != 0,
            });
        }
    }
    out
}

/// GPSNav or QZSNav.
fn gps_nav(b: &[u8], qzss: bool) -> Option<Eph> {
    if b.len() < 140 {
        return None;
    }
    let (wnc, tow_ms) = block_time(b)?;
    let sat = svid_sat(b[14])?;
    let want = if qzss {
        GnssSystem::Qzss
    } else {
        GnssSystem::Gps
    };
    if sat.sys != want {
        return None;
    }
    let toc = i64::from(u4(b, 32));
    let toe = i64::from(u4(b, 88));
    let week = expand_week(i64::from(u2(b, 138)), 1024, wnc);
    let week_c = expand_week(i64::from(u2(b, 136)), 1024, wnc);
    let ttm = tow_ms / 1000 + (wnc - week) * 604_800;
    let iode = i64::from(b[24]);
    let fit = if qzss {
        f64::from(b[26])
    } else if b[26] == 0 {
        4.0
    } else {
        6.0
    };
    let vals = vec![
        Some(f4(b, 44)),
        Some(f4(b, 40)),
        Some(f4(b, 36)),
        Some(iode as f64),
        Some(f4(b, 48)),
        Some(f4(b, 52) * GPS_PI),
        Some(f8(b, 56) * GPS_PI),
        Some(f4(b, 64)),
        Some(f8(b, 68)),
        Some(f4(b, 76)),
        Some(f8(b, 80)),
        Some(toe as f64),
        Some(f4(b, 92)),
        Some(f8(b, 96) * GPS_PI),
        Some(f4(b, 104)),
        Some(f8(b, 108) * GPS_PI),
        Some(f4(b, 116)),
        Some(f8(b, 120) * GPS_PI),
        Some(f4(b, 128) * GPS_PI),
        Some(f4(b, 132) * GPS_PI),
        Some(f64::from(b[18])),
        Some(week as f64),
        Some(f64::from(b[21])),
        Some(ura(b[19])),
        Some(f64::from(b[20])),
        Some(f4(b, 28)),
        Some(f64::from(u2(b, 22))),
        Some(ttm as f64),
        Some(fit),
    ];
    let sys = sat.sys.letter();
    Some(Eph {
        rec: NavRecord {
            sys,
            sat: format!("{sys}{:02}", sat.prn),
            epoch: gps_time(week_c, toc),
            vals,
        },
        iod: iode,
    })
}

/// GLONav; `leap` is GPS - UTC.
fn glo_nav(b: &[u8], leap: i64, fcn: &mut BTreeMap<Sat, i32>) -> Option<Eph> {
    if b.len() < 96 {
        return None;
    }
    let (wnc, tow_ms) = block_time(b)?;
    let sat = svid_sat(b[14])?;
    if sat.sys != GnssSystem::Glonass {
        return None;
    }
    let k = i32::from(b[15]) - 8;
    fcn.entry(sat).or_insert(k);
    let week = expand_week(i64::from(u2(b, 80)), 1024, wnc);
    // t_oe is the GPS time of tb; RINEX wants UTC
    let epoch_us = (week * 604_800 + i64::from(u4(b, 76)) - leap) * 1_000_000;
    // GLONASS frames start on 30 s boundaries of the UTC week
    let frame_tow = ((tow_ms / 1000 - leap).rem_euclid(604_800) / 30 * 30) as f64;
    let health = i64::from(b[85] & 4 != 0 || b[90] != 0);
    let vals = vec![
        Some(-f4(b, 68)),
        Some(f4(b, 64)),
        Some(frame_tow),
        Some(f8(b, 16)),
        Some(f4(b, 40)),
        Some(f4(b, 52)),
        Some(health as f64),
        Some(f8(b, 24)),
        Some(f4(b, 44)),
        Some(f4(b, 56)),
        Some(f64::from(k)),
        Some(f8(b, 32)),
        Some(f4(b, 48)),
        Some(f4(b, 60)),
        Some(f64::from(b[84])),
        None,
        Some(0.999_999_999_999e9),
        Some(15.0),
        None,
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'R',
            sat: format!("R{:02}", sat.prn),
            epoch: DateTime::from_micros_since_gps0(epoch_us),
            vals,
        },
        iod: i64::from(u2(b, 86)),
    })
}

/// The RINEX Galileo health bits of SBF's `Health_OSSOL`.
fn gal_health(h: u16) -> u16 {
    let mut r = 0;
    if h & 1 != 0 {
        r |= ((h >> 1) & 1) | (((h >> 2) & 3) << 1);
    }
    if h & 0x100 != 0 {
        r |= (((h >> 9) & 1) << 3) | (((h >> 10) & 3) << 4);
    }
    if h & 0x10 != 0 {
        r |= (((h >> 5) & 1) << 6) | (((h >> 6) & 3) << 7);
    }
    r
}

/// GALNav.
fn gal_nav(b: &[u8]) -> Option<Eph> {
    if b.len() < 148 {
        return None;
    }
    let (wnc, tow_ms) = block_time(b)?;
    let sat = svid_sat(b[14])?;
    if sat.sys != GnssSystem::Galileo {
        return None;
    }
    // 2: I/NAV (clock for E5b/E1), 16: F/NAV (clock for E5a/E1)
    let (source, sisa) = match b[15] {
        2 => (517.0, b[134]),
        16 => (258.0, b[133]),
        _ => return None,
    };
    let toe = i64::from(u4(b, 100));
    let toc = i64::from(u4(b, 104));
    let week = expand_week(i64::from(u2(b, 124)), 4096, wnc);
    let week_c = expand_week(i64::from(u2(b, 126)), 4096, wnc);
    let ttm = tow_ms / 1000 + (wnc - week) * 604_800;
    let iod = i64::from(u2(b, 128));
    let vals = vec![
        Some(f8(b, 116)),
        Some(f4(b, 112)),
        Some(f4(b, 108)),
        Some(iod as f64),
        Some(f4(b, 88)),
        Some(f4(b, 72) * GPS_PI),
        Some(f8(b, 24) * GPS_PI),
        Some(f4(b, 76)),
        Some(f8(b, 32)),
        Some(f4(b, 80)),
        Some(f8(b, 16)),
        Some(toe as f64),
        Some(f4(b, 92)),
        Some(f8(b, 56) * GPS_PI),
        Some(f4(b, 96)),
        Some(f8(b, 40) * GPS_PI),
        Some(f4(b, 84)),
        Some(f8(b, 48) * GPS_PI),
        Some(f4(b, 64) * GPS_PI),
        Some(f4(b, 68) * GPS_PI),
        Some(source),
        Some(week as f64),
        None,
        Some(sisa_m(sisa)),
        Some(f64::from(gal_health(u2(b, 130)))),
        f4_opt(b, 136),
        f4_opt(b, 140),
        Some(ttm as f64),
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'E',
            sat: format!("E{:02}", sat.prn),
            epoch: gps_time(week_c, toc),
            vals,
        },
        iod: iod * 32 + i64::from(b[15]),
    })
}

/// BDSNav. Times are BDT.
fn bds_nav(b: &[u8]) -> Option<Eph> {
    if b.len() < 140 {
        return None;
    }
    let (wnc, tow_ms) = block_time(b)?;
    let sat = svid_sat(b[14])?;
    if sat.sys != GnssSystem::BeiDou {
        return None;
    }
    let toc = i64::from(u4(b, 32));
    let toe = i64::from(u4(b, 88));
    // reception time in BDT (GPST - 14 s), s since BDT week 0
    let t_rx = (wnc - BDS_WEEK0) * 604_800 + tow_ms / 1000 - 14;
    let ref_week = t_rx.div_euclid(604_800);
    let week = expand_week(i64::from(u2(b, 138)), 8192, ref_week);
    let week_c = expand_week(i64::from(u2(b, 136)), 8192, ref_week);
    let aode = i64::from(b[21]);
    let vals = vec![
        Some(f4(b, 44)),
        Some(f4(b, 40)),
        Some(f4(b, 36)),
        Some(aode as f64),
        Some(f4(b, 48)),
        Some(f4(b, 52) * GPS_PI),
        Some(f8(b, 56) * GPS_PI),
        Some(f4(b, 64)),
        Some(f8(b, 68)),
        Some(f4(b, 76)),
        Some(f8(b, 80)),
        Some(toe as f64),
        Some(f4(b, 92)),
        Some(f8(b, 96) * GPS_PI),
        Some(f4(b, 104)),
        Some(f8(b, 108) * GPS_PI),
        Some(f4(b, 116)),
        Some(f8(b, 120) * GPS_PI),
        Some(f4(b, 128) * GPS_PI),
        Some(f4(b, 132) * GPS_PI),
        None,
        Some(week as f64),
        None,
        Some(ura(b[18])),
        Some(f64::from(b[19])),
        f4_opt(b, 24),
        f4_opt(b, 28),
        Some((t_rx - week * 604_800) as f64),
        Some(f64::from(b[20])),
    ];
    Some(Eph {
        rec: NavRecord {
            sys: 'C',
            sat: format!("C{:02}", sat.prn),
            epoch: gps_time(week_c + BDS_WEEK0, toc),
            vals,
        },
        iod: aode,
    })
}

/// Klobuchar alpha and beta of GPSIon / BDSIon.
fn iono(b: &[u8]) -> Option<([f64; 4], [f64; 4])> {
    if b.len() < 48 {
        return None;
    }
    let mut a = [0.0; 4];
    let mut g = [0.0; 4];
    for k in 0..4 {
        a[k] = f4(b, 16 + 4 * k);
        g[k] = f4(b, 32 + 4 * k);
    }
    Some((a, g))
}

/// GPSUtc.
fn gps_utc(b: &[u8]) -> Option<GpsUtc> {
    if b.len() < 37 {
        return None;
    }
    let (wnc, _) = block_time(b)?;
    Some(GpsUtc {
        a0: f8(b, 20),
        a1: f4(b, 16),
        tot: i64::from(u4(b, 28)),
        wnt: expand_week(i64::from(b[32]), 256, wnc),
        dtls: i64::from(b[33] as i8),
        dn: i64::from(b[35]),
        wnlsf: expand_week(i64::from(b[34]), 256, wnc),
        dtlsf: i64::from(b[36] as i8),
    })
}

/// WGS84 latitude, longitude (rad) and ellipsoidal height (m) to ECEF.
fn geodetic_to_ecef(lat: f64, lon: f64, h: f64) -> [f64; 3] {
    let a = 6_378_137.0;
    let f = 1.0 / 298.257_223_563;
    let e2 = f * (2.0 - f);
    let n = a / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    [
        (n + h) * lat.cos() * lon.cos(),
        (n + h) * lat.cos() * lon.sin(),
        (n * (1.0 - e2) + h) * lat.sin(),
    ]
}

/// The receiver and antenna facts of ReceiverSetup.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Setup {
    /// Marker name.
    pub marker: String,
    /// Marker number.
    pub marker_number: String,
    /// Observer.
    pub observer: String,
    /// Agency.
    pub agency: String,
    /// Receiver serial number.
    pub rx_serial: String,
    /// Receiver GNSS engine name.
    pub rx_name: String,
    /// Receiver firmware version.
    pub rx_version: String,
    /// Product name (revision 3 and later).
    pub product: String,
    /// Antenna serial number.
    pub ant_serial: String,
    /// Antenna type as logged.
    pub ant_type: String,
    /// Antenna delta H, E, N (m).
    pub delta_hen: [f64; 3],
}

fn setup(b: &[u8], rev: u8) -> Option<Setup> {
    if b.len() < 268 {
        return None;
    }
    Some(Setup {
        marker: text(b, 16, 60),
        marker_number: text(b, 76, 20),
        observer: text(b, 96, 20),
        agency: text(b, 116, 40),
        rx_serial: text(b, 156, 20),
        rx_name: text(b, 176, 20),
        rx_version: text(b, 196, 20),
        product: if rev >= 3 && b.len() >= 368 {
            text(b, 328, 40)
        } else {
            String::new()
        },
        ant_serial: text(b, 216, 20),
        ant_type: text(b, 236, 20),
        delta_hen: [f4(b, 256), f4(b, 260), f4(b, 264)],
    })
}

/// What [`decode`] found.
#[derive(Debug)]
pub struct Sbf {
    /// The observations, with the header facts the log carries.
    pub obs: StationObs,
    /// GLONASS frequency numbers.
    pub glonass_fcn: BTreeMap<Sat, i32>,
    /// Broadcast ephemerides, each once.
    pub nav: Vec<Eph>,
    /// Ionosphere, UTC and leap seconds for the navigation header.
    pub nav_header: NavHeader,
    /// ReceiverSetup, when logged.
    pub setup: Option<Setup>,
    /// Blocks seen, by block number.
    pub blocks: BTreeMap<u16, usize>,
    /// Blocks that failed their CRC.
    pub bad_crc: usize,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// GPS - UTC from ReceiverTime, else GPSUtc.
fn leap_of(blocks: &[Block<'_>]) -> Option<i64> {
    let from_time = blocks.iter().find_map(|k| {
        (k.num == RX_TIME && k.b.len() >= 22 && k.b[20] as i8 != -128)
            .then(|| i64::from(k.b[20] as i8))
    });
    from_time.or_else(|| {
        blocks
            .iter()
            .find_map(|k| (k.num == GPS_UTC).then(|| gps_utc(k.b)).flatten())
            .map(|u| u.dtls)
    })
}

/// Decodes an SBF log.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn decode(data: &[u8], name: &str) -> Sbf {
    let (blocks, bad_crc) = blocks(data);
    let leap = leap_of(&blocks);
    let mut counts: BTreeMap<u16, usize> = BTreeMap::new();
    let mut epochs = Epochs::default();
    let mut fcn: BTreeMap<Sat, i32> = BTreeMap::new();
    let mut lock: LockMap = BTreeMap::new();
    let mut dropped = Dropped::default();
    let mut nav: Vec<Eph> = Vec::new();
    let mut seen: BTreeSet<(String, i64, i64)> = BTreeSet::new();
    let mut nh = NavHeader {
        program: geomaticist_gnss_rinex::writer::PROGRAM.to_owned(),
        date: crate::now_utc_text(),
        leap,
        ..NavHeader::default()
    };
    let mut set: Option<Setup> = None;
    let mut cart: Vec<[f64; 3]> = Vec::new();
    let mut geod: Vec<[f64; 3]> = Vec::new();
    for k in &blocks {
        *counts.entry(k.num).or_insert(0) += 1;
        let b = k.b;
        let eph = match k.num {
            MEAS_EPOCH => {
                let Some((wn, tow)) = block_time(b) else {
                    continue;
                };
                let t = wn as f64 * WEEK + tow as f64 * 0.001;
                for m in meas_epoch(b, &mut fcn, &mut dropped) {
                    let mut lli = 0i8;
                    if m.phase.is_finite() {
                        let prev = m.lock.and_then(|l| lock.insert((m.sat, m.tail), (t, l)));
                        let slip = match (m.lock, prev) {
                            (Some(l), None) => l == 0,
                            (Some(l), Some((tp, lp))) => l < lp || f64::from(l) + 1.0 < t - tp,
                            (None, _) => false,
                        };
                        lli = i8::from(slip) | (i8::from(m.half) << 1);
                    }
                    let c = epochs.cell(t, m.sat, m.tail);
                    c.code = m.code;
                    c.phase = m.phase;
                    c.doppler = m.doppler;
                    c.snr = m.snr;
                    c.lli = lli;
                }
                None
            }
            GPS_NAV => gps_nav(b, false),
            QZS_NAV => gps_nav(b, true),
            GLO_NAV => {
                let l = leap.or_else(|| {
                    block_time(b)
                        .map(|(w, t)| i64::from(leap_at(w as f64 * WEEK + t as f64 * 0.001)))
                });
                l.and_then(|l| glo_nav(b, l, &mut fcn))
            }
            GAL_NAV => gal_nav(b),
            BDS_NAV => bds_nav(b),
            GPS_ION => {
                nh.gps_iono = iono(b).or(nh.gps_iono);
                None
            }
            BDS_ION => {
                nh.bds_iono = iono(b).or(nh.bds_iono);
                None
            }
            GPS_UTC => {
                nh.gps_utc = gps_utc(b).or(nh.gps_utc);
                None
            }
            RX_SETUP => {
                set = setup(b, k.rev).or(set);
                None
            }
            PVT_CART if b.len() >= 40 && b[14] & 15 != 0 => {
                let p = [f8(b, 16), f8(b, 24), f8(b, 32)];
                if p.iter().all(|v| *v != F_DNU) {
                    cart.push(p);
                }
                None
            }
            PVT_GEOD if b.len() >= 40 && b[14] & 15 != 0 => {
                let (lat, lon, h) = (f8(b, 16), f8(b, 24), f8(b, 32));
                if lat != F_DNU && lon != F_DNU && h != F_DNU {
                    geod.push(geodetic_to_ecef(lat, lon, h));
                }
                None
            }
            _ => None,
        };
        if let Some(e) = eph
            && seen.insert(e.key())
        {
            nav.push(e);
        }
    }
    let mut warnings = Vec::new();
    if bad_crc > 0 {
        warnings.push(format!(
            "{bad_crc} SBF blocks failed their CRC and were left out"
        ));
    }
    if dropped.unknown_signals > 0 {
        warnings.push(format!(
            "{} SBF measurements of signals not known here were left out",
            dropped.unknown_signals
        ));
    }
    if dropped.scrambled > 0 {
        warnings.push(format!(
            "{} MeasEpoch blocks were scrambled (no measurement permission on the receiver) and were \
             left out",
            dropped.scrambled
        ));
    }
    let marker = set.as_ref().map_or_default(|s| s.marker.clone());
    let mut obs = epochs.into_obs(if marker.is_empty() { name } else { &marker });
    if let Some(s) = &set {
        obs.antenna = norm_antenna(&s.ant_type);
        let model = if s.product.is_empty() {
            &s.rx_name
        } else {
            &s.product
        };
        obs.receiver = format!("SEPT {}", model.to_uppercase()).trim().to_owned();
        obs.delta_hen = s.delta_hen;
    }
    let fixes = if cart.is_empty() { &geod } else { &cart };
    if !fixes.is_empty() {
        let n = fixes.len() as f64;
        for (i, a) in obs.approx.iter_mut().enumerate() {
            *a = fixes.iter().map(|p| p[i]).sum::<f64>() / n;
        }
    }
    fcn.retain(|s, _| obs.sats.contains_key(s));
    Sbf {
        obs,
        glonass_fcn: fcn,
        nav,
        nav_header: nh,
        setup: set,
        blocks: counts,
        bad_crc,
        warnings,
    }
}

/// Whether a stream reads as SBF: two good blocks in its head.
#[must_use]
pub fn recognises(head: &[u8]) -> bool {
    blocks(head).0.len() >= 2
}

/// Decodes an SBF log into `outdir` as RINEX 3.04 `{stem}.obs`, and RINEX 3.05 `{stem}.{yy}p` when it
/// carries navigation blocks.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let (stem, name) = stem_and_name(raw);
    let r = decode(&data, &name);
    if r.obs.t.is_empty() {
        return Err(DecodeError::Format(format!(
            "{}: no Septentrio MeasEpoch blocks found; log them on the receiver to get data that can \
             be processed",
            raw.display()
        )));
    }
    let list: Vec<String> = r.blocks.iter().map(|(n, c)| format!("{n} x{c}")).collect();
    log(&format!(
        "Septentrio: {} epochs, {} satellites, {} navigation records; blocks {}",
        r.obs.t.len(),
        r.obs.sats.len(),
        r.nav.len(),
        list.join(", ")
    ));
    let s = r.setup.clone().unwrap_or_default();
    let info = serde_json::json!({
        "receiver": r.obs.receiver,
        "receiver_serial": s.rx_serial,
        "firmware": s.rx_version,
        "antenna": r.obs.antenna,
        "antenna_serial": s.ant_serial,
        "marker": s.marker,
        "marker_number": s.marker_number,
        "observer": s.observer,
        "agency": s.agency,
        "delta_hen": r.obs.delta_hen,
        "epochs": r.obs.t.len(),
        "nav_records": r.nav.len(),
        "bad_crc": r.bad_crc,
    });
    let mut out = write_obs(
        raw,
        outdir,
        &r.obs,
        r.glonass_fcn,
        "Septentrio SBF",
        r.warnings,
        info,
    )?;
    if !r.nav.is_empty() {
        let yy = gps_datetime(r.obs.t[0]).year.rem_euclid(100);
        let p = outdir.join(format!("{stem}.{yy:02}p"));
        crate::write_atomic(&p, nav_text(&r.nav, &r.nav_header).as_bytes())?;
        out.nav.push(p);
    }
    Ok(out)
}

/// An SBF block of number `num`, revision `rev` and body `body` (the bytes after `WNc`), with its TOW (ms),
/// week, padding and CRC (for tests and for writing logs).
#[must_use]
pub fn block(num: u16, rev: u8, tow_ms: u32, wnc: u16, body: &[u8]) -> Vec<u8> {
    let mut out = vec![b'$', b'@', 0, 0];
    out.extend((num | (u16::from(rev) << 13)).to_le_bytes());
    out.extend([0, 0]);
    out.extend(tow_ms.to_le_bytes());
    out.extend(wnc.to_le_bytes());
    out.extend_from_slice(body);
    while out.len() % 4 != 0 {
        out.push(0);
    }
    let len = out.len() as u16;
    out[6..8].copy_from_slice(&len.to_le_bytes());
    let crc = crc16(&out[4..]);
    out[2..4].copy_from_slice(&crc.to_le_bytes());
    out
}

#[cfg(test)]
mod tests;
