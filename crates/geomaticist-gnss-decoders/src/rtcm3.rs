//! RTCM 3 logged streams: Multiple Signal Messages (MSM4, MSM5, MSM6, MSM7) of GPS, GLONASS, Galileo, QZSS
//! and BeiDou, and the legacy observables 1002 / 1004 (GPS) and 1010 / 1012 (GLONASS), into the engine's observations, with the station messages 1005 / 1006 (reference point),
//! 1007 / 1008 / 1033 (antenna and receiver descriptors).
//!
//! Written from the message layouts of the RTCM 10403 standard: the frame
//! (preamble `0xD3`, 10-bit length, CRC-24Q), the MSM header with its satellite, signal and cell masks, the
//! satellite and signal data blocks, and the signal tables that map an MSM signal number to a RINEX code.
//! The rule that turns the lock-time indicator into a loss-of-lock flag (a lock time that went down, or
//! that is zero twice in a row) is the one RTKLIB uses (T. Takasu, BSD-2-Clause, `rtcm3.c`); no RTKLIB code
//! is copied.
//!
//! What a stream does not say and the caller has to: **the week.** MSM times are milliseconds of the week
//! (of the day for GLONASS), so [`decode`] takes an approximate time and puts every epoch in the week
//! (day) nearest to it. GLONASS phases need the satellite's frequency channel: MSM5 and MSM7 carry it,
//! MSM4 and MSM6 do not, and then it must be given (from a navigation file) or the GLONASS phases are left
//! out with a warning.
//!
//! [`encode`] writes the same messages, for tests and for turning observations back into a stream.

use std::collections::BTreeMap;

use geomaticist_gnss_core::consts::{CLIGHT, WEEK};
use geomaticist_gnss_core::obs::norm_antenna;
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::time::leap_at;
use geomaticist_gnss_core::{GnssSystem, ObsCode, Sat, StationObs};

use crate::epochs::{Cell, Epochs};

/// Range of one millisecond of light time, m.
const RANGE_MS: f64 = CLIGHT * 0.001;

/// CRC-24Q of `data` (polynomial 0x1864CFB), as RTCM 3 frames carry it.
#[must_use]
pub fn crc24q(data: &[u8]) -> u32 {
    let mut crc: u32 = 0;
    for &b in data {
        crc ^= u32::from(b) << 16;
        for _ in 0..8 {
            crc <<= 1;
            if crc & 0x0100_0000 != 0 {
                crc ^= 0x0186_4CFB;
            }
        }
    }
    crc & 0x00FF_FFFF
}

/// A reader of big-endian bit fields.
struct Bits<'a> {
    d: &'a [u8],
    pos: usize,
}

impl Bits<'_> {
    fn left(&self) -> usize {
        (self.d.len() * 8).saturating_sub(self.pos)
    }

    fn u(&mut self, n: usize) -> u64 {
        let mut v = 0u64;
        for _ in 0..n {
            let byte = self.d.get(self.pos / 8).copied().unwrap_or(0);
            v = (v << 1) | u64::from((byte >> (7 - self.pos % 8)) & 1);
            self.pos += 1;
        }
        v
    }

    fn i(&mut self, n: usize) -> i64 {
        let v = self.u(n);
        if n < 64 && v & (1 << (n - 1)) != 0 {
            v as i64 - (1i64 << n)
        } else {
            v as i64
        }
    }

    fn text(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| char::from(self.u(8) as u8))
            .collect::<String>()
            .trim()
            .to_owned()
    }
}

/// A writer of big-endian bit fields.
#[derive(Default)]
struct BitsOut {
    d: Vec<u8>,
    n: usize,
}

impl BitsOut {
    fn u(&mut self, bits: usize, v: u64) {
        for k in (0..bits).rev() {
            if self.n.is_multiple_of(8) {
                self.d.push(0);
            }
            if (v >> k) & 1 != 0 {
                let last = self.d.len() - 1;
                self.d[last] |= 1 << (7 - self.n % 8);
            }
            self.n += 1;
        }
    }

    fn i(&mut self, bits: usize, v: i64) {
        self.u(bits, (v as u64) & (u64::MAX >> (64 - bits)));
    }

    fn text(&mut self, s: &str) {
        let b: Vec<u8> = s.bytes().take(31).collect();
        self.u(8, b.len() as u64);
        for c in b {
            self.u(8, u64::from(c));
        }
    }

    /// The payload wrapped in a frame: preamble, length, CRC.
    fn frame(self) -> Vec<u8> {
        let len = self.d.len();
        let mut out = vec![0xD3, ((len >> 8) & 0x03) as u8, (len & 0xFF) as u8];
        out.extend(self.d);
        let crc = crc24q(&out);
        out.extend([(crc >> 16) as u8, (crc >> 8) as u8, crc as u8]);
        out
    }
}

/// The MSM message number base of a system (`1070` GPS ... `1120` BeiDou), by its place in the standard.
fn msm_system(msg: u16) -> Option<GnssSystem> {
    match msg / 10 {
        107 => Some(GnssSystem::Gps),
        108 => Some(GnssSystem::Glonass),
        109 => Some(GnssSystem::Galileo),
        111 => Some(GnssSystem::Qzss),
        112 => Some(GnssSystem::BeiDou),
        _ => None,
    }
}

fn msm_base(sys: GnssSystem) -> Option<u16> {
    match sys {
        GnssSystem::Gps => Some(1070),
        GnssSystem::Glonass => Some(1080),
        GnssSystem::Galileo => Some(1090),
        GnssSystem::Qzss => Some(1110),
        GnssSystem::BeiDou => Some(1120),
        GnssSystem::Navic | GnssSystem::Sbas => None,
    }
}

/// The RINEX band and attribute (`"1C"`, `"2W"` ...) of MSM signal number `id` (1 to 32) of a system;
/// `None` for numbers the standard reserves.
#[must_use]
pub fn msm_signal(sys: GnssSystem, id: usize) -> Option<&'static str> {
    const GPS: [&str; 32] = [
        "", "1C", "1P", "1W", "", "", "", "2C", "2P", "2W", "", "", "", "", "2S", "2L", "2X", "",
        "", "", "", "5I", "5Q", "5X", "", "", "", "", "", "1S", "1L", "1X",
    ];
    const GLO: [&str; 32] = [
        "", "1C", "1P", "", "", "", "", "2C", "2P", "", "", "", "", "", "", "", "", "", "", "", "",
        "", "", "", "", "", "", "", "", "", "", "",
    ];
    const GAL: [&str; 32] = [
        "", "1C", "1A", "1B", "1X", "1Z", "", "6C", "6A", "6B", "6X", "6Z", "", "7I", "7Q", "7X",
        "", "8I", "8Q", "8X", "", "5I", "5Q", "5X", "", "", "", "", "", "", "", "",
    ];
    const QZS: [&str; 32] = [
        "", "1C", "", "", "", "", "", "", "6S", "6L", "6X", "", "", "", "2S", "2L", "2X", "", "",
        "", "", "5I", "5Q", "5X", "", "", "", "", "", "1S", "1L", "1X",
    ];
    const BDS: [&str; 32] = [
        "", "2I", "2Q", "2X", "", "", "", "6I", "6Q", "6X", "", "", "", "7I", "7Q", "7X", "", "",
        "", "", "", "5D", "5P", "5X", "7D", "", "", "", "", "1D", "1P", "1X",
    ];
    let t = match sys {
        GnssSystem::Gps => &GPS,
        GnssSystem::Glonass => &GLO,
        GnssSystem::Galileo => &GAL,
        GnssSystem::Qzss => &QZS,
        GnssSystem::BeiDou => &BDS,
        GnssSystem::Navic | GnssSystem::Sbas => return None,
    };
    t.get(id.wrapping_sub(1)).copied().filter(|s| !s.is_empty())
}

fn msm_signal_id(sys: GnssSystem, tail: &str) -> Option<usize> {
    (1..=32).find(|&id| msm_signal(sys, id) == Some(tail))
}

/// The shortest lock time (ms) an MSM lock-time indicator stands for: the 4-bit one of MSM4 / MSM5, or the
/// 10-bit one of MSM6 / MSM7.
#[must_use]
pub fn lock_time_ms(indicator: u32, extended: bool) -> u32 {
    if !extended {
        return if indicator == 0 {
            0
        } else {
            1 << (indicator + 4)
        };
    }
    match indicator {
        0..=63 => indicator,
        64..=95 => 2 * indicator - 64,
        96..=127 => 4 * indicator - 256,
        128..=159 => 8 * indicator - 768,
        160..=191 => 16 * indicator - 2048,
        192..=223 => 32 * indicator - 5120,
        224..=255 => 64 * indicator - 12288,
        256..=287 => 128 * indicator - 28672,
        288..=319 => 256 * indicator - 65536,
        320..=351 => 512 * indicator - 147_456,
        352..=383 => 1024 * indicator - 327_680,
        384..=415 => 2048 * indicator - 720_896,
        416..=447 => 4096 * indicator - 1_572_864,
        448..=479 => 8192 * indicator - 3_407_872,
        480..=511 => 16384 * indicator - 7_340_032,
        512..=543 => 32768 * indicator - 15_728_640,
        544..=575 => 65536 * indicator - 33_554_432,
        576..=607 => 131_072 * indicator - 71_303_168,
        608..=639 => 262_144 * indicator - 150_994_944,
        640..=671 => 524_288 * indicator - 318_767_104,
        672..=703 => 1_048_576 * indicator - 671_088_640,
        _ => 67_108_864,
    }
}

/// The station facts a stream carries.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StationInfo {
    /// Reference station id.
    pub id: u16,
    /// Antenna reference point, ECEF m (1005 / 1006).
    pub arp: Option<[f64; 3]>,
    /// Antenna height above the mark, m (1006).
    pub height: Option<f64>,
    /// Antenna descriptor (1007 / 1008 / 1033).
    pub antenna: String,
    /// Antenna serial number (1008 / 1033).
    pub antenna_serial: String,
    /// Receiver type (1033).
    pub receiver: String,
    /// Receiver firmware (1033).
    pub firmware: String,
    /// Receiver serial number (1033).
    pub receiver_serial: String,
}

/// What [`decode`] found.
#[derive(Debug, Clone, PartialEq)]
pub struct Rtcm3 {
    /// The observations (`approx`, `delta_hen`, `antenna`, `receiver` from the station messages).
    pub obs: StationObs,
    /// The station messages.
    pub station: StationInfo,
    /// GLONASS frequency channels learnt from MSM5 / MSM7.
    pub glonass_fcn: BTreeMap<Sat, i32>,
    /// Messages decoded, by number.
    pub messages: BTreeMap<u16, usize>,
    /// Frames whose CRC failed.
    pub crc_errors: usize,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// The phase (cycles) nearest to `prev` among `cp` and `cp ± k·1500` cycles: the legacy messages carry
/// the phase minus the range, which wraps (RTKLIB `adjcp`).
fn adjust_cycles(prev: Option<f64>, cp: f64) -> f64 {
    let Some(p) = prev.filter(|p| p.is_finite()) else {
        return cp;
    };
    if cp < p - 750.0 {
        cp + 1500.0
    } else if cp > p + 750.0 {
        cp - 1500.0
    } else {
        cp
    }
}

/// State the legacy (1002 / 1004 / 1010 / 1012) messages carry from epoch to epoch.
#[derive(Default)]
struct Legacy {
    /// Last phase - range (cycles) per satellite and signal, for [`adjust_cycles`].
    cp: BTreeMap<(Sat, &'static str), f64>,
    /// Last lock-time indicator per satellite and signal.
    lock: BTreeMap<(Sat, &'static str), u64>,
}

impl Legacy {
    /// Puts one legacy signal into the epoch table. `base` is the L1 pseudorange: both phase-range fields of
    /// 1004 / 1012 are counted from it (DF018 / DF048), not from the signal's own range.
    #[allow(clippy::too_many_arguments)]
    fn put(
        &mut self,
        epochs: &mut Epochs,
        t: f64,
        sat: Sat,
        tail: &'static str,
        fr: Option<f64>,
        pr: f64,
        base: f64,
        ppr: Option<f64>,
        lock: u64,
        cnr: Option<f64>,
    ) {
        let c = epochs.cell(t, sat, tail);
        c.code = pr;
        if let (Some(fr), Some(ppr)) = (fr, ppr) {
            let lam = CLIGHT / fr;
            let cp = adjust_cycles(self.cp.get(&(sat, tail)).copied(), ppr / lam);
            self.cp.insert((sat, tail), cp);
            c.phase = base / lam + cp;
        }
        if let Some(s) = cnr.filter(|s| *s > 0.0) {
            c.snr = s;
        }
        let prev = self.lock.insert((sat, tail), lock);
        let slip = prev.is_some_and(|p| (lock == 0 && p == 0) || lock < p);
        c.lli = i8::from(slip);
    }
}

/// The frames of a byte stream: payloads with a valid CRC, and the number of failed ones.
fn frames(data: &[u8]) -> (Vec<&[u8]>, usize) {
    let mut out = Vec::new();
    let mut bad = 0;
    let mut i = 0;
    while i + 6 <= data.len() {
        if data[i] != 0xD3 || data[i + 1] & 0xFC != 0 {
            i += 1;
            continue;
        }
        let len = (usize::from(data[i + 1] & 3) << 8) | usize::from(data[i + 2]);
        let end = i + 3 + len + 3;
        if end > data.len() {
            i += 1;
            continue;
        }
        let crc = (u32::from(data[end - 3]) << 16)
            | (u32::from(data[end - 2]) << 8)
            | u32::from(data[end - 1]);
        if crc24q(&data[i..end - 3]) == crc {
            out.push(&data[i + 3..end - 3]);
            i = end;
        } else {
            bad += 1;
            i += 1;
        }
    }
    (out, bad)
}

/// GPS seconds of an MSM epoch field, in the week (day) nearest to `near`.
fn msm_time(sys: GnssSystem, field: u64, near: f64) -> f64 {
    let nearest = |t: f64, period: f64| t + ((near - t) / period).round() * period;
    match sys {
        GnssSystem::Glonass => {
            // day of the week (0 Sunday, 7 unknown) and milliseconds of the day, Moscow time (UTC + 3 h)
            let tod = (field & 0x07FF_FFFF) as f64 * 0.001;
            let dow = (field >> 27) & 7;
            let utc = if dow < 7 {
                nearest(dow as f64 * 86_400.0 + tod - 10_800.0, WEEK)
            } else {
                nearest(tod - 10_800.0, 86_400.0)
            };
            utc + f64::from(leap_at(utc))
        }
        GnssSystem::BeiDou => nearest(field as f64 * 0.001 + 14.0, WEEK),
        _ => nearest(field as f64 * 0.001, WEEK),
    }
}

/// Decodes a logged RTCM 3 stream.
///
/// - `near`: GPS seconds within three days of the data (within half a day for a stream with GLONASS
///   only); it settles the week.
/// - `fcn`: GLONASS frequency channels, for streams whose messages do not carry them.
/// - `name`: the station name to give the observations.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn decode(data: &[u8], near: f64, fcn: &BTreeMap<Sat, i32>, name: &str) -> Rtcm3 {
    let (frames, crc_errors) = frames(data);
    let mut station = StationInfo::default();
    let mut messages: BTreeMap<u16, usize> = BTreeMap::new();
    let mut warnings = Vec::new();
    let mut glo_fcn: BTreeMap<Sat, i32> = fcn.clone();
    let mut epochs = Epochs::default();
    let mut lock: BTreeMap<(Sat, &'static str), u32> = BTreeMap::new();
    let mut legacy = Legacy::default();
    let (mut no_fcn, mut unknown_sig) = (0usize, 0usize);
    for f in frames {
        let mut b = Bits { d: f, pos: 0 };
        if b.left() < 12 {
            continue;
        }
        let msg = b.u(12) as u16;
        *messages.entry(msg).or_insert(0) += 1;
        match msg {
            1005 | 1006 => {
                if b.left() < 140 {
                    continue;
                }
                station.id = b.u(12) as u16;
                b.u(6 + 4);
                let x = b.i(38) as f64 / 1e4;
                b.u(2);
                let y = b.i(38) as f64 / 1e4;
                b.u(2);
                let z = b.i(38) as f64 / 1e4;
                station.arp = Some([x, y, z]);
                if msg == 1006 && b.left() >= 16 {
                    station.height = Some(b.u(16) as f64 / 1e4);
                }
            }
            1007 | 1008 | 1033 => {
                station.id = b.u(12) as u16;
                let n = b.u(8) as usize;
                station.antenna = b.text(n);
                b.u(8);
                if msg != 1007 {
                    let m = b.u(8) as usize;
                    station.antenna_serial = b.text(m);
                }
                if msg == 1033 {
                    let n = b.u(8) as usize;
                    station.receiver = b.text(n);
                    let n = b.u(8) as usize;
                    station.firmware = b.text(n);
                    let n = b.u(8) as usize;
                    station.receiver_serial = b.text(n);
                }
            }
            1002 | 1004 => {
                // legacy GPS L1 (/ L2) observables with the ambiguity, layout of RTCM 10403 3.5.1
                if b.left() < 52 {
                    continue;
                }
                station.id = b.u(12) as u16;
                let t = msm_time(GnssSystem::Gps, b.u(30), near);
                b.u(1);
                let nsat = b.u(5) as usize;
                b.u(1 + 3);
                let size = if msg == 1002 { 74 } else { 125 };
                if b.left() < nsat * size {
                    continue;
                }
                for _ in 0..nsat {
                    let prn = b.u(6) as u8;
                    let code1 = b.u(1);
                    let pr1 = b.u(24) as f64;
                    let ppr1 = b.i(20);
                    let lock1 = b.u(7);
                    let amb = b.u(8) as f64;
                    let cnr1 = b.u(8) as f64 * 0.25;
                    let l2 = (msg == 1004).then(|| (b.u(2), b.i(14), b.i(20), b.u(7), b.u(8)));
                    if prn == 0 || prn >= 40 {
                        continue; // SBAS PRNs 40+ are left out
                    }
                    let sat = Sat::new(GnssSystem::Gps, prn);
                    let pr = pr1 * 0.02 + amb * 299_792.458;
                    let tail1 = if code1 == 1 { "1P" } else { "1C" };
                    let ppr = (ppr1 != -(1 << 19)).then_some(ppr1 as f64 * 0.0005);
                    let f1 = freq(sat, 1, None).ok();
                    legacy.put(
                        &mut epochs,
                        t,
                        sat,
                        tail1,
                        f1,
                        pr,
                        pr,
                        ppr,
                        lock1,
                        Some(cnr1),
                    );
                    if let Some((code2, pr21, ppr2, lock2, cnr2)) = l2 {
                        if pr21 == -(1 << 13) {
                            continue;
                        }
                        let tail2 = ["2X", "2P", "2W", "2W"][code2 as usize];
                        let ppr = (ppr2 != -(1 << 19)).then_some(ppr2 as f64 * 0.0005);
                        let f2 = freq(sat, 2, None).ok();
                        let p2 = pr + pr21 as f64 * 0.02;
                        let cnr = Some(cnr2 as f64 * 0.25);
                        legacy.put(&mut epochs, t, sat, tail2, f2, p2, pr, ppr, lock2, cnr);
                    }
                }
            }
            1010 | 1012 => {
                // legacy GLONASS L1 (/ L2) observables with the ambiguity
                if b.left() < 49 {
                    continue;
                }
                station.id = b.u(12) as u16;
                let tod = b.u(27);
                let t = msm_time(GnssSystem::Glonass, tod | (7 << 27), near);
                b.u(1);
                let nsat = b.u(5) as usize;
                b.u(1 + 3);
                let size = if msg == 1010 { 79 } else { 130 };
                if b.left() < nsat * size {
                    continue;
                }
                for _ in 0..nsat {
                    let prn = b.u(6) as u8;
                    let code1 = b.u(1);
                    let fcn = b.u(5) as i32 - 7;
                    let pr1 = b.u(25) as f64;
                    let ppr1 = b.i(20);
                    let lock1 = b.u(7);
                    let amb = b.u(7) as f64;
                    let cnr1 = b.u(8) as f64 * 0.25;
                    let l2 = (msg == 1012).then(|| (b.u(2), b.i(14), b.i(20), b.u(7), b.u(8)));
                    if prn == 0 {
                        continue;
                    }
                    let sat = Sat::new(GnssSystem::Glonass, prn);
                    glo_fcn.insert(sat, fcn);
                    let pr = pr1 * 0.02 + amb * 599_584.916;
                    let tail1 = if code1 == 1 { "1P" } else { "1C" };
                    let ppr = (ppr1 != -(1 << 19)).then_some(ppr1 as f64 * 0.0005);
                    let f1 = freq(sat, 1, Some(fcn)).ok();
                    legacy.put(
                        &mut epochs,
                        t,
                        sat,
                        tail1,
                        f1,
                        pr,
                        pr,
                        ppr,
                        lock1,
                        Some(cnr1),
                    );
                    if let Some((code2, pr21, ppr2, lock2, cnr2)) = l2 {
                        if pr21 == -(1 << 13) {
                            continue;
                        }
                        let tail2 = if code2 == 0 { "2C" } else { "2P" };
                        let ppr = (ppr2 != -(1 << 19)).then_some(ppr2 as f64 * 0.0005);
                        let f2 = freq(sat, 2, Some(fcn)).ok();
                        let p2 = pr + pr21 as f64 * 0.02;
                        let cnr = Some(cnr2 as f64 * 0.25);
                        legacy.put(&mut epochs, t, sat, tail2, f2, p2, pr, ppr, lock2, cnr);
                    }
                }
            }
            _ => {
                let (Some(sys), kind) = (msm_system(msg), msg % 10) else {
                    continue;
                };
                if !(4..=7).contains(&kind) || b.left() < 157 {
                    continue;
                }
                station.id = b.u(12) as u16;
                let t = msm_time(sys, b.u(30), near);
                b.u(1 + 3 + 7 + 2 + 2 + 1 + 3);
                let sats: Vec<u8> = (1..=64u8).filter(|_| b.u(1) == 1).collect();
                let sigs: Vec<usize> = (1..=32usize).filter(|_| b.u(1) == 1).collect();
                let (ns, ng) = (sats.len(), sigs.len());
                if ns * ng > 64 || b.left() < ns * ng {
                    continue;
                }
                let cells: Vec<bool> = (0..ns * ng).map(|_| b.u(1) == 1).collect();
                let nc = cells.iter().filter(|c| **c).count();
                let rate = kind == 5 || kind == 7;
                let ext = kind >= 6;
                let need = ns * (18 + if rate { 18 } else { 0 })
                    + nc * (if ext {
                        20 + 24 + 10 + 1 + 10
                    } else {
                        15 + 22 + 4 + 1 + 6
                    } + if rate { 15 } else { 0 });
                if b.left() < need {
                    continue;
                }
                // satellite data, field by field
                let int_ms: Vec<u64> = (0..ns).map(|_| b.u(8)).collect();
                let info: Vec<u64> = if rate {
                    (0..ns).map(|_| b.u(4)).collect()
                } else {
                    vec![15; ns]
                };
                let mod_ms: Vec<u64> = (0..ns).map(|_| b.u(10)).collect();
                let rough_rate: Vec<i64> = if rate {
                    (0..ns).map(|_| b.i(14)).collect()
                } else {
                    vec![-8192; ns]
                };
                // signal data, field by field
                let (pr_bits, ph_bits, lock_bits, snr_bits) = if ext {
                    (20, 24, 10, 10)
                } else {
                    (15, 22, 4, 6)
                };
                let fine_pr: Vec<i64> = (0..nc).map(|_| b.i(pr_bits)).collect();
                let fine_ph: Vec<i64> = (0..nc).map(|_| b.i(ph_bits)).collect();
                let lock_ind: Vec<u32> = (0..nc).map(|_| b.u(lock_bits) as u32).collect();
                let half: Vec<bool> = (0..nc).map(|_| b.u(1) == 1).collect();
                let cnr: Vec<u64> = (0..nc).map(|_| b.u(snr_bits)).collect();
                let fine_rate: Vec<i64> = if rate {
                    (0..nc).map(|_| b.i(15)).collect()
                } else {
                    vec![-16384; nc]
                };
                let (pr_scale, ph_scale) = if ext {
                    (2f64.powi(-29), 2f64.powi(-31))
                } else {
                    (2f64.powi(-24), 2f64.powi(-29))
                };
                let key = (t * 1000.0).round() as i64;
                let mut c = 0usize;
                for (si, &prn) in sats.iter().enumerate() {
                    let sat = Sat::new(sys, prn);
                    if sys == GnssSystem::Glonass && info[si] != 15 {
                        glo_fcn.insert(sat, info[si] as i32 - 7);
                    }
                    for (gi, &sig) in sigs.iter().enumerate() {
                        if !cells[si * ng + gi] {
                            continue;
                        }
                        let q = c;
                        c += 1;
                        let Some(tail) = msm_signal(sys, sig) else {
                            unknown_sig += 1;
                            continue;
                        };
                        if int_ms[si] == 255 {
                            continue;
                        }
                        let rough = int_ms[si] as f64 + mod_ms[si] as f64 / 1024.0;
                        let band = tail.as_bytes()[0] - b'0';
                        let band = if sys == GnssSystem::BeiDou && band == 2 {
                            2
                        } else {
                            band
                        };
                        let fr = freq(sat, band, glo_fcn.get(&sat).copied()).ok();
                        let pr_bad = fine_pr[q] == -(1 << (pr_bits - 1));
                        let ph_bad = fine_ph[q] == -(1 << (ph_bits - 1));
                        let code = if pr_bad {
                            f64::NAN
                        } else {
                            RANGE_MS * (rough + fine_pr[q] as f64 * pr_scale)
                        };
                        let mut phase = f64::NAN;
                        let mut doppler = f64::NAN;
                        match fr {
                            Some(fr) => {
                                let lam = CLIGHT / fr;
                                if !ph_bad {
                                    phase = RANGE_MS * (rough + fine_ph[q] as f64 * ph_scale) / lam;
                                }
                                if rate && rough_rate[si] != -8192 && fine_rate[q] != -16384 {
                                    doppler =
                                        -(rough_rate[si] as f64 + fine_rate[q] as f64 * 1e-4) / lam;
                                }
                            }
                            None => no_fcn += 1,
                        }
                        let snr = if cnr[q] == 0 {
                            f64::NAN
                        } else if ext {
                            cnr[q] as f64 / 16.0
                        } else {
                            cnr[q] as f64
                        };
                        // loss of lock: the lock time went down, or is zero twice in a row
                        let lt = lock_time_ms(lock_ind[q], ext);
                        let prev = lock.insert((sat, tail), lt);
                        let slip = prev.is_some_and(|p| lt < p || (lt == 0 && p == 0));
                        let lli = i8::from(slip) | (i8::from(half[q]) << 1);
                        *epochs.cell(key as f64 * 0.001, sat, tail) = Cell {
                            code,
                            phase,
                            doppler,
                            snr,
                            lli,
                        };
                    }
                }
            }
        }
    }
    if crc_errors > 0 {
        warnings.push(format!(
            "{crc_errors} frames failed their checksum and were left out"
        ));
    }
    if no_fcn > 0 {
        warnings.push(format!(
            "{no_fcn} GLONASS phases left out: the stream does not say the satellites' frequency channels \
             (give a GLONASS navigation file)"
        ));
    }
    if unknown_sig > 0 {
        warnings.push(format!(
            "{unknown_sig} observations of signals the standard does not name were left out"
        ));
    }
    let no_amb: usize = [1001u16, 1003, 1009, 1011]
        .iter()
        .filter_map(|m| messages.get(m))
        .sum();
    if no_amb > 0 {
        warnings.push(format!(
            "{no_amb} legacy observation messages without the range ambiguity (1001 / 1003 / 1009 / 1011) \
             were left out"
        ));
    }
    let mut st = epochs.into_obs(name);
    st.approx = station.arp.unwrap_or([0.0; 3]);
    st.delta_hen = [station.height.unwrap_or(0.0), 0.0, 0.0];
    st.antenna = norm_antenna(&station.antenna);
    st.receiver.clone_from(&station.receiver);
    let glonass_fcn: BTreeMap<Sat, i32> = glo_fcn
        .into_iter()
        .filter(|(s, _)| st.sats.contains_key(s))
        .collect();
    Rtcm3 {
        obs: st,
        station,
        glonass_fcn,
        messages,
        crc_errors,
        warnings,
    }
}

/// Which MSM to write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Msm {
    /// Compact pseudorange and phase, 1 dB-Hz signal strength.
    Msm4,
    /// MSM4 plus Doppler and GLONASS channels.
    Msm5,
    /// Full resolution, with Doppler and GLONASS channels.
    Msm7,
}

/// Station messages of a stream: 1006 (reference point and height) and 1033 (descriptors).
#[must_use]
pub fn encode_station(info: &StationInfo) -> Vec<u8> {
    let mut out = Vec::new();
    if let Some(p) = info.arp {
        let mut b = BitsOut::default();
        b.u(12, 1006);
        b.u(12, u64::from(info.id));
        b.u(6, 0);
        b.u(4, 0b1110);
        b.i(38, (p[0] * 1e4).round() as i64);
        b.u(2, 0);
        b.i(38, (p[1] * 1e4).round() as i64);
        b.u(2, 0);
        b.i(38, (p[2] * 1e4).round() as i64);
        b.u(16, (info.height.unwrap_or(0.0) * 1e4).round() as u64);
        out.extend(b.frame());
    }
    if !info.antenna.is_empty() || !info.receiver.is_empty() {
        let mut b = BitsOut::default();
        b.u(12, 1033);
        b.u(12, u64::from(info.id));
        b.text(&info.antenna);
        b.u(8, 0);
        b.text(&info.antenna_serial);
        b.text(&info.receiver);
        b.text(&info.firmware);
        b.text(&info.receiver_serial);
        out.extend(b.frame());
    }
    out
}

/// The observations of `st` as MSM frames, one message per system and epoch. Phases are written as the
/// standard wants them, as a range near the pseudorange: a whole number of cycles is taken out of each
/// satellite's phase once (and again, with the lock time restarted, if the phase drifts out of the
/// message's reach), so a decoded phase differs from the original by a constant number of cycles per arc.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn encode(st: &StationObs, msm: Msm, station_id: u16, fcn: &BTreeMap<Sat, i32>) -> Vec<u8> {
    let mut out = Vec::new();
    let kind: u16 = match msm {
        Msm::Msm4 => 4,
        Msm::Msm5 => 5,
        Msm::Msm7 => 7,
    };
    let (rate, ext) = (kind != 4, kind == 7);
    // per (satellite, signal): cycles taken out of the phase, and the epoch its lock began
    let mut offset: BTreeMap<(Sat, &'static str), f64> = BTreeMap::new();
    let mut lock0: BTreeMap<(Sat, &'static str), f64> = BTreeMap::new();
    let systems = [
        GnssSystem::Gps,
        GnssSystem::Glonass,
        GnssSystem::Galileo,
        GnssSystem::Qzss,
        GnssSystem::BeiDou,
    ];
    for (i, &t) in st.t.iter().enumerate() {
        let present: Vec<GnssSystem> = systems
            .iter()
            .copied()
            .filter(|s| {
                st.sats
                    .iter()
                    .any(|(sat, g)| sat.sys == *s && g.values().any(|x| x.values[i].is_finite()))
            })
            .collect();
        for (pi, &sys) in present.iter().enumerate() {
            let Some(base) = msm_base(sys) else { continue };
            // the cells of this system at this epoch
            let mut sat_list: Vec<Sat> = Vec::new();
            let mut sig_list: Vec<usize> = Vec::new();
            let mut vals: BTreeMap<(u8, usize), (f64, f64, f64, f64, i8, &'static str)> =
                BTreeMap::new();
            for (sat, sigs) in st
                .sats
                .iter()
                .filter(|(s, _)| s.sys == sys && (1..=64).contains(&s.prn))
            {
                for (code, g) in sigs.iter().filter(|(c, _)| c.kind_char() == 'C') {
                    let tail_owned = code.as_str()[1..].to_owned();
                    let Some(id) = msm_signal_id(sys, &tail_owned) else {
                        continue;
                    };
                    let Some(tail) = msm_signal(sys, id) else {
                        continue;
                    };
                    let get = |k: char| {
                        format!("{k}{tail}")
                            .parse::<ObsCode>()
                            .ok()
                            .and_then(|c| sigs.get(&c))
                    };
                    let p = g.values[i];
                    if !p.is_finite() {
                        continue;
                    }
                    let l = get('L').map_or(f64::NAN, |x| x.values[i]);
                    let lli = get('L').and_then(|x| x.lli.as_ref()).map_or(0, |x| x[i]);
                    let d = get('D').map_or(f64::NAN, |x| x.values[i]);
                    let s = get('S').map_or(f64::NAN, |x| x.values[i]);
                    if !sat_list.contains(sat) {
                        sat_list.push(*sat);
                    }
                    if !sig_list.contains(&id) {
                        sig_list.push(id);
                    }
                    vals.insert((sat.prn, id), (p, l, d, s, lli, tail));
                }
            }
            sig_list.sort_unstable();
            if sat_list.is_empty() || sat_list.len() * sig_list.len() > 64 {
                continue;
            }
            let mut b = BitsOut::default();
            b.u(12, u64::from(base + kind));
            b.u(12, u64::from(station_id));
            let field = match sys {
                GnssSystem::Glonass => {
                    let utc = t - f64::from(leap_at(t)) + 10_800.0;
                    let day = (utc / 86_400.0).floor();
                    let dow = (day as i64).rem_euclid(7) as u64;
                    (dow << 27) | ((utc - day * 86_400.0) * 1000.0).round() as u64
                }
                GnssSystem::BeiDou => ((t - 14.0).rem_euclid(WEEK) * 1000.0).round() as u64,
                _ => (t.rem_euclid(WEEK) * 1000.0).round() as u64,
            };
            b.u(30, field);
            b.u(1, u64::from(pi + 1 < present.len()));
            b.u(3 + 7 + 2 + 2 + 1 + 3, 0);
            for prn in 1..=64u8 {
                b.u(1, u64::from(sat_list.iter().any(|s| s.prn == prn)));
            }
            for id in 1..=32usize {
                b.u(1, u64::from(sig_list.contains(&id)));
            }
            for sat in &sat_list {
                for id in &sig_list {
                    b.u(1, u64::from(vals.contains_key(&(sat.prn, *id))));
                }
            }
            // rough range per satellite: from its first pseudorange
            let rough: Vec<(u64, u64, f64)> = sat_list
                .iter()
                .map(|sat| {
                    let p = sig_list
                        .iter()
                        .find_map(|id| vals.get(&(sat.prn, *id)).map(|v| v.0))
                        .unwrap_or(0.0);
                    let r = (p / RANGE_MS * 1024.0).round() / 1024.0;
                    let int = r.floor();
                    (int as u64, ((r - int) * 1024.0).round() as u64, r)
                })
                .collect();
            let chan = |sat: &Sat| fcn.get(sat).copied();
            let lam_of = |sat: &Sat, tail: &str| -> Option<f64> {
                freq(*sat, tail.as_bytes()[0] - b'0', chan(sat))
                    .ok()
                    .map(|f| CLIGHT / f)
            };
            let rough_rate: Vec<i64> = sat_list
                .iter()
                .map(|sat| {
                    sig_list
                        .iter()
                        .find_map(|id| {
                            let v = vals.get(&(sat.prn, *id))?;
                            let lam = lam_of(sat, v.5)?;
                            v.2.is_finite().then(|| (-v.2 * lam).round() as i64)
                        })
                        .unwrap_or(-8192)
                })
                .collect();
            for r in &rough {
                b.u(8, r.0);
            }
            if rate {
                for sat in &sat_list {
                    b.u(4, chan(sat).map_or(15, |f| (f + 7) as u64));
                }
            }
            for r in &rough {
                b.u(10, r.1);
            }
            if rate {
                for r in &rough_rate {
                    b.i(14, *r);
                }
            }
            let (pr_bits, ph_bits, lock_bits, snr_bits) = if ext {
                (20, 24, 10, 10)
            } else {
                (15, 22, 4, 6)
            };
            let (pr_scale, ph_scale) = if ext {
                (2f64.powi(29), 2f64.powi(31))
            } else {
                (2f64.powi(24), 2f64.powi(29))
            };
            let mut cells: Vec<(i64, i64, u64, u64, u64, i64)> = Vec::new();
            for (si, sat) in sat_list.iter().enumerate() {
                for id in &sig_list {
                    let Some(&(p, l, d, s, lli, tail)) = vals.get(&(sat.prn, *id)) else {
                        continue;
                    };
                    let fine_pr = ((p / RANGE_MS - rough[si].2) * pr_scale).round() as i64;
                    let mut fine_ph = -(1i64 << (ph_bits - 1));
                    let mut lock_ind = 0u64;
                    let mut fine_rate = -16384i64;
                    if let Some(lam) = lam_of(sat, tail) {
                        if l.is_finite() {
                            let limit = (1i64 << (ph_bits - 1)) - 1;
                            let k = (*sat, tail);
                            let need = |off: f64| {
                                (((l - off) * lam / RANGE_MS - rough[si].2) * ph_scale).round()
                                    as i64
                            };
                            let mut off = offset.get(&k).copied();
                            if off.is_none_or(|o| need(o).abs() > limit) || lli & 1 != 0 {
                                off = Some(((l * lam - p) / lam).round());
                                lock0.insert(k, t);
                            }
                            let off = off.unwrap_or(0.0);
                            offset.insert(k, off);
                            fine_ph = need(off);
                            let locked_ms = ((t - lock0.get(&k).copied().unwrap_or(t)) * 1000.0)
                                .max(0.0) as u32;
                            let n_ind: u32 = if ext { 704 } else { 16 };
                            lock_ind = u64::from(
                                (0..n_ind)
                                    .rev()
                                    .find(|q| lock_time_ms(*q, ext) <= locked_ms)
                                    .unwrap_or(0),
                            );
                        }
                        if d.is_finite() && rough_rate[si] != -8192 {
                            fine_rate = ((-d * lam - rough_rate[si] as f64) * 1e4).round() as i64;
                        }
                    }
                    let snr = if !s.is_finite() {
                        0
                    } else if ext {
                        (s * 16.0).round() as u64
                    } else {
                        s.round() as u64
                    };
                    cells.push((
                        fine_pr,
                        fine_ph,
                        lock_ind,
                        u64::from(lli & 2 != 0),
                        snr,
                        fine_rate,
                    ));
                }
            }
            for c in &cells {
                b.i(pr_bits, c.0);
            }
            for c in &cells {
                b.i(ph_bits, c.1);
            }
            for c in &cells {
                b.u(lock_bits, c.2);
            }
            for c in &cells {
                b.u(1, c.3);
            }
            for c in &cells {
                b.u(snr_bits, c.4);
            }
            if rate {
                for c in &cells {
                    b.i(15, c.5);
                }
            }
            out.extend(b.frame());
        }
    }
    out
}

/// Whether the start of a file is an RTCM 3 stream: at least two frames pass their CRC, more than fail.
#[must_use]
pub fn recognises(head: &[u8]) -> bool {
    let (good, bad) = frames(head);
    good.len() >= 2 && good.len() > bad
}

/// Whether a file name's extension says it is a logged RTCM 3 stream (`.rtcm3`, `.rtcm`, `.rtc`, `.rt3`).
#[must_use]
pub fn is_rtcm_ext(ext: &str) -> bool {
    matches!(
        ext.to_ascii_lowercase().as_str(),
        "rtcm3" | "rtcm" | "rtc" | "rt3"
    )
}

/// Decodes a logged RTCM 3 stream into `outdir` as a RINEX 3.04 observation file named after the log
/// (`name.rtcm3` -> `name.obs`). `near` is an approximate GPS time of the data (it settles the week); when
/// `None`, the file's modification time is taken, which is right for a log that was not copied in a way
/// that changed it. `fcn` gives GLONASS frequency channels for streams that do not carry them.
///
/// Navigation messages are not decoded: use a broadcast navigation file for the day.
pub fn convert(
    raw: &std::path::Path,
    outdir: &std::path::Path,
    log: &mut dyn FnMut(&str),
    near: Option<f64>,
    fcn: &BTreeMap<Sat, i32>,
) -> Result<crate::Converted, crate::DecodeError> {
    use crate::DecodeError;
    let data =
        std::fs::read(raw).map_err(|e| DecodeError::Io(format!("{}: {e}", raw.display())))?;
    let mut warnings = Vec::new();
    let near = if let Some(t) = near {
        t
    } else {
        let unix = std::fs::metadata(raw)
            .and_then(|m| m.modified())
            .ok()
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map_or(0.0, |d| d.as_secs_f64());
        warnings.push(
            "no date given: the week of the data was taken from the file's modification time"
                .to_owned(),
        );
        // 1970-01-01 is 315 964 800 s before 1980-01-06; GPS runs ahead of UTC by the leap seconds
        unix - 315_964_800.0 + 18.0
    };
    let stem = raw
        .file_stem()
        .map_or_default(|s| s.to_string_lossy().into_owned());
    let name: String = stem.chars().take(4).collect::<String>().to_uppercase();
    let mut r = decode(&data, near, fcn, &name);
    if r.obs.t.is_empty() {
        return Err(DecodeError::Format(format!(
            "{}: no observation messages found (MSM4 to MSM7, or legacy 1002 / 1004 / 1010 / 1012) ({} frames with a valid checksum)",
            raw.display(),
            r.messages.values().sum::<usize>()
        )));
    }
    warnings.append(&mut r.warnings);
    if r.station.arp.is_none() {
        warnings.push("the stream has no station position message (1005 / 1006)".to_owned());
    }
    let list: Vec<String> = r
        .messages
        .iter()
        .map(|(m, n)| format!("{m} x{n}"))
        .collect();
    log(&format!(
        "RTCM 3: {} epochs, {} satellites; messages {}",
        r.obs.t.len(),
        r.obs.sats.len(),
        list.join(", ")
    ));
    let header = geomaticist_gnss_rinex::station::StationHeader {
        program: "geomaticist-gnss".to_owned(),
        run_by: String::new(),
        date: crate::now_utc_text(),
        comments: vec![format!(
            "Decoded from the RTCM 3 log {}",
            raw.file_name()
                .map_or_default(|s| s.to_string_lossy().into_owned())
        )],
        glonass_fcn: r.glonass_fcn.clone(),
    };
    let text = geomaticist_gnss_rinex::station::station_to_rinex3(&r.obs, &header);
    let obs = outdir.join(format!("{stem}.obs"));
    crate::write_atomic(&obs, text.as_bytes())?;
    let info = serde_json::json!({
        "marker": name,
        "antenna": r.station.antenna,
        "receiver": r.station.receiver,
        "firmware": r.station.firmware,
        "serial": r.station.receiver_serial,
        "arp_height": r.station.height,
        "arp_xyz": r.station.arp,
        "station_id": r.station.id,
        "epochs": r.obs.t.len(),
    });
    Ok(crate::Converted {
        obs,
        nav: Vec::new(),
        warnings,
        info,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use geomaticist_gnss_core::obs::Signal;
    use geomaticist_gnss_core::signal::code_key;

    const T0: f64 = 2300.0 * 604_800.0 + 345_600.0;

    #[test]
    fn a_log_becomes_a_rinex_file() {
        let (st, fcn) = station();
        let dir = std::env::temp_dir().join(format!("rtcm3-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("base.rtcm3");
        let info = StationInfo {
            id: 7,
            arp: Some([4_395_470.631_2, 617_743.110_5, 4_565_607.209_9]),
            height: Some(1.5),
            antenna: "TRM59800.00     NONE".into(),
            ..StationInfo::default()
        };
        let mut bytes = encode_station(&info);
        bytes.extend(encode(&st, Msm::Msm7, 7, &fcn));
        std::fs::write(&raw, bytes).unwrap();
        let mut lines = Vec::new();
        let c = convert(
            &raw,
            &dir,
            &mut |m| lines.push(m.to_owned()),
            Some(T0),
            &BTreeMap::new(),
        )
        .unwrap();
        assert!(lines[0].contains("3 epochs, 5 satellites"), "{lines:?}");
        assert_eq!(c.info["arp_height"], 1.5);
        let text = std::fs::read(&c.obs).unwrap();
        let back = geomaticist_gnss_rinex::obs::read_obs_bytes_kinds(
            &c.obs, &text, "", None, None, "GREC", "CLSD",
        )
        .unwrap();
        assert_eq!(back.t, st.t);
        assert_eq!(back.antenna, "TRM59800.00     NONE");
        assert_eq!(back.delta_hen, [1.5, 0.0, 0.0]);
        let g5: Sat = "G05".parse().unwrap();
        let c1: ObsCode = "C1C".parse().unwrap();
        assert!((back.sats[&g5][&c1].values[0] - st.sats[&g5][&c1].values[0]).abs() < 1e-3);
        assert!(String::from_utf8_lossy(&text).contains("GLONASS SLOT / FRQ #"));
        assert!(is_rtcm_ext("RTCM3") && !is_rtcm_ext("t04"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A small station: GPS dual frequency, Galileo, BeiDou and GLONASS, three epochs, one slip.
    fn station() -> (StationObs, BTreeMap<Sat, i32>) {
        let mut st = StationObs::new("TEST", "");
        st.t = vec![T0, T0 + 1.0, T0 + 2.0];
        st.power_fail = vec![false; 3];
        let mut fcn = BTreeMap::new();
        fcn.insert("R07".parse::<Sat>().unwrap(), -4);
        let mut order = 0;
        for (k, (sat, tails)) in [
            ("G05", vec!["1C", "2W"]),
            ("G17", vec!["1C", "2W", "5Q"]),
            ("E11", vec!["1C", "7Q"]),
            ("C23", vec!["2I", "6I"]),
            ("R07", vec!["1C", "2P"]),
        ]
        .into_iter()
        .enumerate()
        {
            let sat: Sat = sat.parse().unwrap();
            for tail in tails {
                let band = tail.as_bytes()[0] - b'0';
                let lam = CLIGHT / freq(sat, band, fcn.get(&sat).copied()).unwrap();
                let range: Vec<f64> = (0..3)
                    .map(|i| 2.1e7 + 1.3e6 * k as f64 + 612.3 * f64::from(i))
                    .collect();
                for kind in ['C', 'L', 'D', 'S'] {
                    let text = format!("{kind}{tail}");
                    let values: Vec<f64> = range
                        .iter()
                        .enumerate()
                        .map(|(i, r)| match kind {
                            'C' => r + 0.37,
                            'L' => {
                                r / lam
                                    + 12_345.25
                                    + if sat.prn == 17 && i == 2 { 40.0 } else { 0.0 }
                            }
                            'D' => -612.3 / lam,
                            _ => 43.25,
                        })
                        .collect();
                    st.sats.entry(sat).or_default().insert(
                        text.parse().unwrap(),
                        Signal {
                            key: code_key(text.as_bytes(), sat.sys.letter(), 3, 3.04),
                            values,
                            lli: (kind == 'L').then(|| vec![0, 0, i8::from(sat.prn == 17)]),
                            first_seen: order,
                        },
                    );
                    order += 1;
                }
            }
        }
        (st, fcn)
    }

    #[test]
    fn crc_and_tables() {
        // the CRC of a frame followed by its CRC is zero
        let mut f = vec![0xD3, 0x00, 0x02, 0x12, 0x34];
        let c = crc24q(&f);
        f.extend([(c >> 16) as u8, (c >> 8) as u8, c as u8]);
        assert_eq!(crc24q(&f), 0);
        assert_eq!(msm_signal(GnssSystem::Gps, 2), Some("1C"));
        assert_eq!(msm_signal(GnssSystem::Gps, 10), Some("2W"));
        assert_eq!(msm_signal(GnssSystem::Galileo, 15), Some("7Q"));
        assert_eq!(msm_signal(GnssSystem::BeiDou, 2), Some("2I"));
        assert_eq!(msm_signal(GnssSystem::Glonass, 9), Some("2P"));
        assert_eq!(msm_signal(GnssSystem::Gps, 1), None);
        assert_eq!(lock_time_ms(0, false), 0);
        assert_eq!(lock_time_ms(1, false), 32);
        assert_eq!(lock_time_ms(15, false), 524_288);
        assert_eq!(lock_time_ms(63, true), 63);
        assert_eq!(lock_time_ms(64, true), 64);
        assert_eq!(lock_time_ms(96, true), 128);
        assert_eq!(lock_time_ms(703, true), 66_060_288);
    }

    #[test]
    fn round_trips_every_msm() {
        let (st, fcn) = station();
        for (msm, pr_tol, ph_tol_m) in [
            (Msm::Msm4, 0.02, 0.0006),
            (Msm::Msm5, 0.02, 0.0006),
            (Msm::Msm7, 0.0006, 0.0002),
        ] {
            let info = StationInfo {
                id: 42,
                arp: Some([4_395_470.631_2, 617_743.110_5, 4_565_607.209_9]),
                height: Some(1.5432),
                antenna: "TRM59800.00     NONE".into(),
                antenna_serial: "A1".into(),
                receiver: "RECEIVER X".into(),
                firmware: "1.2".into(),
                receiver_serial: "0001".into(),
            };
            let mut bytes = encode_station(&info);
            bytes.extend(encode(&st, msm, 42, &fcn));
            // junk between frames and a damaged frame are skipped
            bytes.extend([0xD3, 0x00, 0x01, 0xFF, 0x00, 0x00, 0x00, 0x11, 0x22]);
            // without channels MSM4 cannot give GLONASS phases; MSM5 and MSM7 carry them
            let given = if msm == Msm::Msm4 {
                fcn.clone()
            } else {
                BTreeMap::new()
            };
            let r = decode(&bytes, T0 + 86_400.0, &given, "BACK");
            assert_eq!(r.station, info);
            assert_eq!(
                r.obs.approx,
                [4_395_470.631_2, 617_743.110_5, 4_565_607.209_9]
            );
            assert_eq!(r.obs.delta_hen, [1.5432, 0.0, 0.0]);
            assert_eq!(r.obs.antenna, "TRM59800.00     NONE");
            assert_eq!(r.obs.t, st.t, "{msm:?}");
            assert_eq!(r.glonass_fcn.get(&"R07".parse().unwrap()), Some(&-4));
            assert_eq!(r.obs.sats.len(), st.sats.len());
            for (sat, sigs) in &st.sats {
                for (code, g) in sigs {
                    if code.kind_char() == 'D' && msm == Msm::Msm4 {
                        assert!(!r.obs.sats[sat].contains_key(code));
                        continue;
                    }
                    let back = &r.obs.sats[sat][code];
                    let band = code.as_str().as_bytes()[1] - b'0';
                    let lam = CLIGHT / freq(*sat, band, fcn.get(sat).copied()).unwrap();
                    match code.kind_char() {
                        'C' => {
                            for i in 0..3 {
                                assert!(
                                    (back.values[i] - g.values[i]).abs() < pr_tol,
                                    "{msm:?} {sat} {code}"
                                );
                            }
                        }
                        'L' => {
                            // equal up to whole cycles, constant while the lock holds
                            let d: Vec<f64> =
                                (0..3).map(|i| back.values[i] - g.values[i]).collect();
                            for v in &d {
                                assert!(
                                    (v - v.round()).abs() * lam < ph_tol_m,
                                    "{msm:?} {sat} {code} {v}"
                                );
                            }
                            assert!((d[1] - d[0]).abs() < 1e-3, "{msm:?} {sat} {code} {d:?}");
                            let lli = back.lli.as_ref().unwrap();
                            assert_eq!(lli[1], 0);
                            assert_eq!(lli[2], i8::from(sat.prn == 17), "{msm:?} {sat} {code}");
                        }
                        'D' => assert!((back.values[0] - g.values[0]).abs() < 0.01),
                        _ => {
                            let tol = if msm == Msm::Msm7 { 0.04 } else { 0.51 };
                            assert!((back.values[0] - g.values[0]).abs() < tol);
                        }
                    }
                }
            }
            assert_eq!(r.crc_errors, 1);
            assert_eq!(r.messages[&1006], 1);
        }
    }

    #[test]
    fn glonass_phase_needs_a_channel_and_bad_frames_are_counted() {
        let (st, fcn) = station();
        let mut bytes = encode(&st, Msm::Msm4, 1, &fcn);
        let r = decode(&bytes, T0, &BTreeMap::new(), "X");
        let r07: Sat = "R07".parse().unwrap();
        assert!(r.obs.sats[&r07].contains_key(&"C1C".parse::<ObsCode>().unwrap()));
        assert!(!r.obs.sats[&r07].contains_key(&"L1C".parse::<ObsCode>().unwrap()));
        assert!(r.warnings.iter().any(|w| w.contains("frequency channels")));
        // one flipped bit in a payload: that frame is dropped and counted
        bytes[10] ^= 0x10;
        let r = decode(&bytes, T0, &fcn, "X");
        assert!(r.crc_errors >= 1);
        assert!(r.warnings.iter().any(|w| w.contains("checksum")));
        // the week comes from the hint
        let later = decode(
            &encode(&st, Msm::Msm7, 1, &fcn),
            T0 + 604_800.0 + 5.0,
            &fcn,
            "X",
        );
        assert!((later.obs.t[0] - (T0 + 604_800.0)).abs() < 1e-6);
    }

    #[test]
    fn legacy_1004_decodes() {
        let mut b = BitsOut {
            d: Vec::new(),
            n: 0,
        };
        b.u(12, 1004);
        b.u(12, 7); // station
        b.u(30, 345_600_000); // GPS time of week, ms
        b.u(1, 0);
        b.u(5, 1); // one satellite
        b.u(4, 0);
        b.u(6, 5); // G05
        b.u(1, 0); // C/A
        b.u(24, 726_397); // pr modulo one light-millisecond, 0.02 m
        b.i(20, 100); // phase - range, 0.0005 m
        b.u(7, 40);
        b.u(8, 70); // whole light-milliseconds
        b.u(8, 180); // C/N0 45 dB-Hz
        b.u(2, 3); // L2 P(Y) cross-correlated
        b.i(14, 50); // L2 - L1 range, 0.02 m
        b.i(20, -200);
        b.u(7, 40);
        b.u(8, 160);
        let r = decode(&b.frame(), T0, &BTreeMap::new(), "X");
        assert_eq!(r.messages[&1004], 1);
        assert!((r.obs.t[0] - T0).abs() < 1e-6);
        let g05: Sat = "G05".parse().unwrap();
        let get = |c: &str| r.obs.sats[&g05][&c.parse::<ObsCode>().unwrap()].values[0];
        let pr = 726_397.0 * 0.02 + 70.0 * 299_792.458;
        assert!((get("C1C") - pr).abs() < 1e-6);
        assert!((get("C2W") - (pr + 1.0)).abs() < 1e-6);
        let l1 = CLIGHT / freq(g05, 1, None).unwrap();
        assert!((get("L1C") - (pr + 0.05) / l1).abs() < 1e-3);
        let l2 = CLIGHT / freq(g05, 2, None).unwrap();
        assert!((get("L2W") - (pr - 0.1) / l2).abs() < 1e-3);
        assert!((get("S1C") - 45.0).abs() < 1e-9);
    }
}
