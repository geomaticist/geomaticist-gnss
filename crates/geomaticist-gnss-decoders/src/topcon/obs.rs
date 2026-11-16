//! Layer 2 of a Topcon `.tps` log: the compact measurement messages `rD` (all measurements of an epoch) and
//! `rP` (the signal list of the same epoch).
//!
//! `rD` body (after the framing, without its CRC):
//!
//! - a 22-byte header: epoch counter (u16), GPS week modulo 1024 (u16), GPS time of week in ms (u32), then
//!   fields of which only the receiver clock offset is known: bits 91..=120 of the little-endian header, a
//!   signed 30-bit count of 1/1024 ns;
//! - per satellite: the reference pseudorange (u32, 0.02 m), the satellite index (USI, u8), the number of
//!   signals minus one (u8), then 14 bytes per signal.
//!
//! A signal block, as a little-endian 112-bit number:
//!
//! | bits | meaning |
//! |---|---|
//! | 0..=15 | first signal: bits 10..=15 the GLONASS slot (0 elsewhere), bits 0..=9 not known; other signals: signed pseudorange minus the reference, 0.02 m |
//! | 16..=21 | carrier-to-noise density, whole dB-Hz |
//! | 22..=47 | signed carrier phase (in metres) minus the reference pseudorange, 1/128000 m |
//! | 48..=57 | lock time, 0.2 s, stops at 1023 |
//! | 60..=63 | carrier-to-noise density, sixteenths of dB-Hz |
//! | 72 | set from about the third epoch of a lock on (phase settled, assumed) |
//! | 73 | not known: the same over a whole track of most satellites, set on some and clear on others |
//! | 87..=110 | signed Doppler, -1 mHz (the receiver's sign is the opposite of RINEX's) |
//!
//! Bits 58..=59 were always 3, bits 64..=71 always 0x53 and bit 79 always set; bits 74..=78, 80..=86 and
//! 111 vary and are not known.
//!
//! The signal of each block is named by the `rP` message of the same epoch (a code byte per signal). A
//! firmware that logs no `rP` leaves the signals to be inferred: the receiver keeps a fixed order per system,
//! and the Doppler of each block, which scales with its carrier frequency, tells which signals of that order are
//! present ([`Inference`]).

use std::collections::HashMap;

use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::consts::{CLIGHT, base_freq, glo_df};

/// Pseudorange unit, m.
pub const RANGE_UNIT: f64 = 0.02;
/// Phase unit, m.
pub const PHASE_UNIT: f64 = 1.0 / 128_000.0;
/// Phase units per pseudorange unit.
pub const PHASE_PER_RANGE: i64 = 2560;
/// Receiver clock unit, s.
pub const CLOCK_UNIT: f64 = 1e-9 / 1024.0;
/// Length of the `rD` header.
pub const HEADER_LEN: usize = 22;
/// Length of a signal block.
pub const BLOCK_LEN: usize = 14;

/// The satellite a USI (receiver satellite index) names: system, PRN (0 for GLONASS: the slot is in the
/// first signal block) and the GLONASS frequency channel. `None` for an index the reader doesn't know.
#[must_use]
pub fn usi_sat(usi: u8) -> Option<(GnssSystem, u8, Option<i8>)> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Qzss, Sbas};
    Some(match usi {
        1..=37 => (Gps, usi, None),
        // GLONASS by frequency channel: 45 + channel, channels -7..=+6
        38..=51 => (Glonass, 0, Some(i8::try_from(i16::from(usi) - 45).ok()?)),
        71..=106 => (Galileo, usi - 70, None),
        // BeiDou PRN 38 and up (assumed beyond C38)
        107..=119 => (BeiDou, usi - 69, None),
        120..=158 => (Sbas, usi - 100, None),
        193..=202 => (Qzss, usi - 192, None),
        203..=239 => (BeiDou, usi - 202, None),
        _ => return None,
    })
}

/// The USI of a satellite (inverse of [`usi_sat`]; GLONASS by its channel).
#[must_use]
pub fn sat_usi(sys: GnssSystem, prn: u8, fcn: Option<i8>) -> Option<u8> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Qzss, Sbas};
    let u = match sys {
        Gps if (1..=37).contains(&prn) => i16::from(prn),
        Glonass => 45 + i16::from(fcn?),
        Galileo if (1..=36).contains(&prn) => 70 + i16::from(prn),
        BeiDou if (1..=37).contains(&prn) => 202 + i16::from(prn),
        BeiDou if (38..=50).contains(&prn) => 69 + i16::from(prn),
        Sbas if (20..=58).contains(&prn) => 100 + i16::from(prn),
        Qzss if (1..=10).contains(&prn) => 192 + i16::from(prn),
        _ => return None,
    };
    u8::try_from(u).ok()
}

/// The RINEX band and attribute of an `rP` signal code byte.
#[must_use]
pub fn rp_signal(sys: GnssSystem, code: u8) -> Option<[u8; 2]> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps};
    Some(match (sys, code) {
        (Gps, 8) => *b"1C",
        (Gps, 25) => *b"1W",
        (Gps, 74) => *b"2W",
        (Gps, 131) => *b"2X",
        (Gps, 188) => *b"5X",
        (Glonass, 8) => *b"1C",
        (Glonass, 17) => *b"1P",
        (Glonass, 66) => *b"2P",
        (Glonass, 59) => *b"2C",
        (Glonass, 116) => *b"3X",
        (Galileo, 32) => *b"1X",
        (Galileo, 188) => *b"5X",
        (Galileo, 125) => *b"7X",
        (Galileo, 158) => *b"8X",
        (BeiDou, 8) => *b"2I",
        (BeiDou, 107) => *b"7I",
        (BeiDou, 194) => *b"7D",
        (BeiDou, 181) => *b"5P",
        _ => return None,
    })
}

/// The `rP` code byte of a signal (inverse of [`rp_signal`]).
#[must_use]
pub fn rp_code(sys: GnssSystem, sig: [u8; 2]) -> Option<u8> {
    (0..=255u8).find(|&c| rp_signal(sys, c) == Some(sig))
}

/// The receiver's signal order of a system, used when there is no `rP` (BeiDou: PRN 1 to 18 are BDS-2).
#[must_use]
pub fn signal_order(sys: GnssSystem, prn: u8) -> &'static [[u8; 2]] {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps};
    const G: [[u8; 2]; 5] = [*b"1C", *b"1W", *b"2W", *b"2X", *b"5X"];
    const R: [[u8; 2]; 5] = [*b"1C", *b"1P", *b"2P", *b"2C", *b"3X"];
    const E: [[u8; 2]; 4] = [*b"1X", *b"5X", *b"7X", *b"8X"];
    const C2: [[u8; 2]; 2] = [*b"2I", *b"7I"];
    const C3: [[u8; 2]; 3] = [*b"2I", *b"7D", *b"5P"];
    match sys {
        Gps => &G,
        Glonass => &R,
        Galileo => &E,
        BeiDou if prn <= 18 => &C2,
        BeiDou => &C3,
        _ => &[],
    }
}

/// Carrier frequency of a RINEX band, Hz (GLONASS FDMA bands need the channel).
#[must_use]
pub fn carrier(sys: GnssSystem, band: u8, fcn: Option<i8>) -> Option<f64> {
    let f = base_freq(sys, band)?;
    if sys == GnssSystem::Glonass && (band == 1 || band == 2) {
        Some(f + f64::from(fcn?) * glo_df(band)?)
    } else {
        Some(f)
    }
}

/// One signal of one satellite at one epoch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Signal {
    /// RINEX band and attribute (`*b"1C"`).
    pub sig: [u8; 2],
    /// Pseudorange, units of [`RANGE_UNIT`].
    pub range: i64,
    /// Carrier phase in metres, units of [`PHASE_UNIT`].
    pub phase: i64,
    /// Doppler in RINEX's sign, mHz.
    pub doppler: i32,
    /// Carrier-to-noise density, 1/16 dB-Hz.
    pub snr16: u16,
    /// Lock time, 0.2 s (stops at 1023).
    pub lock: u16,
    /// Bit 72 of the block.
    pub settled: bool,
    /// Bit 73 of the block (meaning not known; kept under its first guess, a half-cycle flag, which the data
    /// doesn't bear out: it holds over whole tracks).
    pub half_cycle_known: bool,
}

impl Signal {
    /// Pseudorange, m.
    #[must_use]
    pub fn range_m(&self) -> f64 {
        self.range as f64 * RANGE_UNIT
    }

    /// Carrier phase, cycles of the carrier `f` (Hz).
    #[must_use]
    pub fn phase_cycles(&self, f: f64) -> f64 {
        self.phase as f64 * PHASE_UNIT * f / CLIGHT
    }

    /// Doppler, Hz.
    #[must_use]
    pub fn doppler_hz(&self) -> f64 {
        f64::from(self.doppler) / 1000.0
    }

    /// Carrier-to-noise density, dB-Hz.
    #[must_use]
    pub fn snr(&self) -> f64 {
        f64::from(self.snr16) / 16.0
    }
}

/// One satellite at one epoch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SatObs {
    /// System.
    pub sys: GnssSystem,
    /// PRN (GLONASS: slot; 0 if the receiver didn't know it).
    pub prn: u8,
    /// The receiver's satellite index.
    pub usi: u8,
    /// GLONASS frequency channel.
    pub fcn: Option<i8>,
    /// Bits 0..=9 of the first block (not known).
    pub info: u16,
    /// Signals, in the receiver's order.
    pub signals: Vec<Signal>,
}

/// One measurement epoch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Epoch {
    /// The receiver's epoch counter.
    pub counter: u16,
    /// GPS week modulo 1024.
    pub week10: u16,
    /// GPS time of week, ms.
    pub tow_ms: u32,
    /// Receiver clock offset, units of [`CLOCK_UNIT`].
    pub clock: i32,
    /// The rest of the 22-byte header (bytes 8..=21, with the clock bits zeroed), kept for the encoder.
    pub header_rest: [u8; 14],
    /// Satellites, in the receiver's order.
    pub sats: Vec<SatObs>,
}

impl Epoch {
    /// Receiver clock offset, s.
    #[must_use]
    pub fn clock_s(&self) -> f64 {
        f64::from(self.clock) * CLOCK_UNIT
    }
}

fn sext(v: u128, bits: u32) -> i64 {
    let v = v as i64;
    let s = 64 - bits;
    (v << s) >> s
}

fn field(v: u128, lo: u32, bits: u32) -> u128 {
    (v >> lo) & ((1u128 << bits) - 1)
}

/// The raw fields of one signal block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Block {
    /// Bits 0..=15 (signed for the second and later signals).
    pub head: u16,
    /// Bits 16..=21.
    pub snr_int: u8,
    /// Bits 22..=47, signed.
    pub phase: i32,
    /// Bits 48..=57.
    pub lock: u16,
    /// Bits 60..=63.
    pub snr_frac: u8,
    /// Bits 72 and 73.
    pub settled: bool,
    /// Bit 73.
    pub half_cycle_known: bool,
    /// Bits 87..=110, signed (receiver sign).
    pub doppler: i32,
    /// The bits not decoded (58..=59, 64..=71, 74..=86, 111), in place.
    pub rest: u128,
}

const KNOWN_MASK: u128 =
    ((1u128 << 58) - 1) | (0xf << 60) | (0b11 << 72) | (((1u128 << 24) - 1) << 87);

impl Block {
    /// Decodes 14 bytes.
    #[must_use]
    pub fn read(b: &[u8]) -> Self {
        let mut a = [0u8; 16];
        a[..BLOCK_LEN].copy_from_slice(&b[..BLOCK_LEN]);
        let v = u128::from_le_bytes(a);
        Self {
            head: field(v, 0, 16) as u16,
            snr_int: field(v, 16, 6) as u8,
            phase: sext(field(v, 22, 26), 26) as i32,
            lock: field(v, 48, 10) as u16,
            snr_frac: field(v, 60, 4) as u8,
            settled: field(v, 72, 1) == 1,
            half_cycle_known: field(v, 73, 1) == 1,
            doppler: sext(field(v, 87, 24), 24) as i32,
            rest: v & !KNOWN_MASK,
        }
    }

    /// The 14 bytes (inverse of [`Block::read`]).
    #[must_use]
    pub fn bytes(&self) -> [u8; BLOCK_LEN] {
        let mut v = self.rest & !KNOWN_MASK;
        v |= u128::from(self.head);
        v |= u128::from(self.snr_int & 0x3f) << 16;
        v |= (u128::from(self.phase as u32) & ((1 << 26) - 1)) << 22;
        v |= u128::from(self.lock & 0x3ff) << 48;
        v |= u128::from(self.snr_frac & 0xf) << 60;
        v |= u128::from(self.settled) << 72;
        v |= u128::from(self.half_cycle_known) << 73;
        v |= (u128::from(self.doppler as u32) & ((1 << 24) - 1)) << 87;
        let a = v.to_le_bytes();
        let mut out = [0u8; BLOCK_LEN];
        out.copy_from_slice(&a[..BLOCK_LEN]);
        out
    }
}

/// The raw layout of an `rD` body: header and, per satellite, (reference range, USI, blocks).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawEpoch {
    /// The 22 header bytes.
    pub header: [u8; HEADER_LEN],
    /// Per satellite: reference pseudorange (0.02 m), USI, signal blocks.
    pub sats: Vec<(u32, u8, Vec<Block>)>,
}

/// Splits an `rD` body (without its CRC) into its parts.
pub fn split_rd(body: &[u8]) -> Result<RawEpoch, String> {
    if body.len() < HEADER_LEN {
        return Err(format!("rD body of {} bytes", body.len()));
    }
    let mut header = [0u8; HEADER_LEN];
    header.copy_from_slice(&body[..HEADER_LEN]);
    let mut sats = Vec::new();
    let mut j = HEADER_LEN;
    while j < body.len() {
        if j + 6 > body.len() {
            return Err("rD: a satellite record is cut".into());
        }
        let range = u32::from_le_bytes([body[j], body[j + 1], body[j + 2], body[j + 3]]);
        let usi = body[j + 4];
        let n = usize::from(body[j + 5]) + 1;
        let end = j + 6 + n * BLOCK_LEN;
        if end > body.len() {
            return Err("rD: signal blocks run past the end".into());
        }
        let blocks = (0..n)
            .map(|k| Block::read(&body[j + 6 + k * BLOCK_LEN..]))
            .collect();
        sats.push((range, usi, blocks));
        j = end;
    }
    Ok(RawEpoch { header, sats })
}

/// The signal code bytes of an `rP` body (without its CRC), per satellite; `None` if its layout doesn't hold.
#[must_use]
pub fn split_rp(body: &[u8]) -> Option<Vec<Vec<u8>>> {
    let n = usize::from(*body.first()?);
    let mut j = 1;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let k = usize::from(*body.get(j)?) + 1;
        let codes = (0..k)
            .map(|q| body.get(j + 1 + 4 * q + 3).copied())
            .collect::<Option<Vec<u8>>>()?;
        out.push(codes);
        j += 1 + 4 * k;
    }
    (j == body.len()).then_some(out)
}

/// Signal inference for logs without `rP`: per satellite, the signal list last chosen.
#[derive(Debug, Clone, Default)]
pub struct Inference {
    last: HashMap<(GnssSystem, u8), Vec<[u8; 2]>>,
}

/// Doppler (Hz) above which the frequency ratios decide.
const RATIO_DOPPLER: f64 = 30.0;

impl Inference {
    /// The signals of `n` blocks with Dopplers `dop` (Hz, block 0 first) of a satellite: the subset of the
    /// receiver's order (first entry always included) whose carrier frequencies best explain the Dopplers. When
    /// the Doppler is too small to tell, the list chosen last for this satellite (same count) is kept.
    pub fn signals(
        &mut self,
        sys: GnssSystem,
        prn: u8,
        fcn: Option<i8>,
        dop: &[f64],
    ) -> Option<Vec<[u8; 2]>> {
        let order = signal_order(sys, prn);
        let n = dop.len();
        if n == 0 || n > order.len() {
            return None;
        }
        let key = (sys, prn);
        if dop[0].abs() <= RATIO_DOPPLER
            && let Some(v) = self.last.get(&key).filter(|v| v.len() == n)
        {
            return Some(v.clone());
        }
        let freq = |s: [u8; 2]| carrier(sys, s[0] - b'0', fcn);
        let f0 = freq(order[0])?;
        let rest: Vec<f64> = order[1..].iter().map(|&s| freq(s)).collect::<Option<_>>()?;
        let mut best: Option<(f64, Vec<usize>)> = None;
        let mut idx = Vec::with_capacity(n - 1);
        search(&rest, n - 1, 0, &mut idx, &mut |idx: &[usize]| {
            let score: f64 = idx
                .iter()
                .enumerate()
                .map(|(k, &i)| (dop[k + 1] - dop[0] * rest[i] / f0).abs())
                .sum();
            if best.as_ref().is_none_or(|b| score < b.0 - 1e-9) {
                best = Some((score, idx.to_vec()));
            }
        });
        let (_, idx) = best?;
        let mut v = vec![order[0]];
        v.extend(idx.iter().map(|&i| order[1 + i]));
        self.last.insert(key, v.clone());
        Some(v)
    }
}

/// Calls `f` with every increasing index list of length `m` over `0..rest.len()`, in lexicographic order.
fn search(rest: &[f64], m: usize, from: usize, idx: &mut Vec<usize>, f: &mut dyn FnMut(&[usize])) {
    if idx.len() == m {
        f(idx);
        return;
    }
    for i in from..rest.len() {
        idx.push(i);
        search(rest, m, i + 1, idx, f);
        idx.pop();
    }
}

/// What decoding an epoch left out.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Dropped {
    /// Satellites with a USI the reader doesn't know.
    pub unknown_sats: usize,
    /// GLONASS satellites whose slot the receiver didn't know.
    pub glonass_no_slot: usize,
    /// Signals whose `rP` code the reader doesn't know, or that couldn't be inferred.
    pub unknown_signals: usize,
    /// Epochs whose `rP` didn't fit their `rD` (the signals were inferred instead).
    pub rp_mismatch: usize,
}

/// Decodes an `rD` body with the `rP` body of the same epoch, if any.
pub fn decode_epoch(
    rd: &[u8],
    rp: Option<&[u8]>,
    inf: &mut Inference,
    dropped: &mut Dropped,
) -> Result<Epoch, String> {
    let raw = split_rd(rd)?;
    let h = &raw.header;
    let hv = {
        let mut a = [0u8; 16];
        a[..HEADER_LEN.min(16)].copy_from_slice(&h[..16]);
        u128::from_le_bytes(a)
    };
    // the clock field runs from bit 91 to bit 120: bytes 11..=15
    let clock = sext(field(hv, 91, 30), 30) as i32;
    let mut header_rest = [0u8; 14];
    header_rest.copy_from_slice(&h[8..22]);
    // zero the clock bits (header bits 91..=120 = rest bits 27..=56)
    {
        let mut a = [0u8; 16];
        a[..14].copy_from_slice(&header_rest);
        let v = u128::from_le_bytes(a) & !(((1u128 << 30) - 1) << 27);
        header_rest.copy_from_slice(&v.to_le_bytes()[..14]);
    }
    let codes = rp.and_then(split_rp).filter(|c| {
        c.len() == raw.sats.len()
            && c.iter()
                .zip(&raw.sats)
                .all(|(a, (_, _, b))| a.len() == b.len())
    });
    if rp.is_some() && codes.is_none() {
        dropped.rp_mismatch += 1;
    }
    let mut sats = Vec::with_capacity(raw.sats.len());
    for (i, (range, usi, blocks)) in raw.sats.iter().enumerate() {
        let Some((sys, mut prn, fcn)) = usi_sat(*usi) else {
            dropped.unknown_sats += 1;
            continue;
        };
        let first = blocks[0];
        if sys == GnssSystem::Glonass {
            prn = (first.head >> 10) as u8;
            if prn == 0 {
                dropped.glonass_no_slot += 1;
                continue;
            }
        }
        let names: Vec<Option<[u8; 2]>> = if let Some(c) = &codes {
            c[i].iter().map(|&b| rp_signal(sys, b)).collect()
        } else {
            let dop: Vec<f64> = blocks
                .iter()
                .map(|b| -f64::from(b.doppler) / 1000.0)
                .collect();
            match inf.signals(sys, prn, fcn, &dop) {
                Some(v) => v.into_iter().map(Some).collect(),
                None => vec![None; blocks.len()],
            }
        };
        let base = i64::from(*range);
        let mut signals = Vec::with_capacity(blocks.len());
        for (k, (b, name)) in blocks.iter().zip(names).enumerate() {
            let Some(sig) = name else {
                dropped.unknown_signals += 1;
                continue;
            };
            let dc = if k == 0 { 0 } else { i64::from(b.head as i16) };
            signals.push(Signal {
                sig,
                range: base + dc,
                phase: base * PHASE_PER_RANGE + i64::from(b.phase),
                doppler: -b.doppler,
                snr16: u16::from(b.snr_int) * 16 + u16::from(b.snr_frac),
                lock: b.lock,
                settled: b.settled,
                half_cycle_known: b.half_cycle_known,
            });
        }
        if signals.is_empty() {
            continue;
        }
        sats.push(SatObs {
            sys,
            prn,
            usi: *usi,
            fcn,
            info: first.head & 0x3ff,
            signals,
        });
    }
    Ok(Epoch {
        counter: u16::from_le_bytes([h[0], h[1]]),
        week10: u16::from_le_bytes([h[2], h[3]]),
        tow_ms: u32::from_le_bytes([h[4], h[5], h[6], h[7]]),
        clock,
        header_rest,
        sats,
    })
}

/// The `rD` body of an epoch (inverse of [`decode_epoch`] for what it decodes; for tests and tools). The
/// first signal of each satellite sets the reference pseudorange; every signal's range must lie within
/// ±655 m of it and its phase within ±262 m of it.
pub fn encode_rd(e: &Epoch) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    out.extend_from_slice(&e.counter.to_le_bytes());
    out.extend_from_slice(&e.week10.to_le_bytes());
    out.extend_from_slice(&e.tow_ms.to_le_bytes());
    let mut a = [0u8; 16];
    a[..14].copy_from_slice(&e.header_rest);
    let mut v = u128::from_le_bytes(a) & !(((1u128 << 30) - 1) << 27);
    v |= (u128::from(e.clock as u32) & ((1 << 30) - 1)) << 27;
    out.extend_from_slice(&v.to_le_bytes()[..14]);
    for s in &e.sats {
        let first = s.signals.first().ok_or("a satellite without signals")?;
        let base = first.range;
        out.extend_from_slice(
            &u32::try_from(base)
                .map_err(|_| "pseudorange out of range")?
                .to_le_bytes(),
        );
        let usi = if s.usi != 0 {
            s.usi
        } else {
            sat_usi(s.sys, s.prn, s.fcn).ok_or("no USI for the satellite")?
        };
        out.push(usi);
        out.push(u8::try_from(s.signals.len() - 1).map_err(|_| "too many signals")?);
        for (k, g) in s.signals.iter().enumerate() {
            let head = if k == 0 {
                (u16::from(if s.sys == GnssSystem::Glonass {
                    s.prn
                } else {
                    0
                }) << 10)
                    | (s.info & 0x3ff)
            } else {
                i16::try_from(g.range - base)
                    .map_err(|_| "pseudorange too far from the reference")? as u16
            };
            let ph = g.phase - base * PHASE_PER_RANGE;
            if !(-(1 << 25)..(1 << 25)).contains(&ph) {
                return Err("phase too far from the reference pseudorange".into());
            }
            let b = Block {
                head,
                snr_int: u8::try_from(g.snr16 / 16).map_err(|_| "snr")?,
                phase: ph as i32,
                lock: g.lock,
                snr_frac: (g.snr16 % 16) as u8,
                settled: g.settled,
                half_cycle_known: g.half_cycle_known,
                doppler: -g.doppler,
                rest: (0b11u128 << 58) | (0x53u128 << 64) | (1u128 << 79),
            };
            out.extend_from_slice(&b.bytes());
        }
    }
    Ok(out)
}

/// The `rP` body of an epoch (signal codes; the other bytes as the receiver writes them: coarse SNR in half
/// dB-Hz and 30000).
pub fn encode_rp(e: &Epoch) -> Result<Vec<u8>, String> {
    let mut out = vec![u8::try_from(e.sats.len()).map_err(|_| "too many satellites")?];
    for s in &e.sats {
        out.push(u8::try_from(s.signals.len() - 1).map_err(|_| "too many signals")?);
        for g in &s.signals {
            out.push(u8::try_from(g.snr16 / 8).unwrap_or(255));
            out.extend_from_slice(&30_000u16.to_le_bytes());
            out.push(rp_code(s.sys, g.sig).ok_or("a signal without an rP code")?);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usi_table_round_trips() {
        for usi in 0..=255u8 {
            if let Some((sys, prn, fcn)) = usi_sat(usi) {
                let p = if sys == GnssSystem::Glonass { 1 } else { prn };
                assert_eq!(sat_usi(sys, p, fcn), Some(usi), "usi {usi}");
            }
        }
        assert_eq!(usi_sat(43), Some((GnssSystem::Glonass, 0, Some(-2))));
        assert_eq!(usi_sat(208), Some((GnssSystem::BeiDou, 6, None)));
        assert_eq!(usi_sat(107), Some((GnssSystem::BeiDou, 38, None)));
        assert_eq!(usi_sat(106), Some((GnssSystem::Galileo, 36, None)));
        assert_eq!(usi_sat(70), None);
    }

    #[test]
    fn block_round_trips() {
        let b = Block {
            head: 0xfff9,
            snr_int: 45,
            phase: -1_234_567,
            lock: 1023,
            snr_frac: 9,
            settled: true,
            half_cycle_known: false,
            doppler: -5_347_123,
            rest: (0b11u128 << 58) | (0x53u128 << 64) | (1u128 << 111),
        };
        assert_eq!(Block::read(&b.bytes()), b);
        // a block put together bit by bit (not through the encoder): Doppler -10.695 Hz (receiver +10695 mHz)
        let v: u128 = 3
            | (45 << 16)
            | (0b11 << 58)
            | (1 << 60)
            | (0x53 << 64)
            | (1 << 79)
            | (0x29c7_u128 << 87); // 10695
        let mut raw = [0u8; BLOCK_LEN];
        raw.copy_from_slice(&v.to_le_bytes()[..BLOCK_LEN]);
        let r = Block::read(&raw);
        assert_eq!(r.snr_int, 45);
        assert_eq!(r.snr_frac, 1);
        assert_eq!(r.doppler, 10695);
        assert_eq!(r.head, 3);
        assert_eq!(r.bytes(), raw);
    }

    #[test]
    fn rp_codes_round_trip() {
        for sys in [
            GnssSystem::Gps,
            GnssSystem::Glonass,
            GnssSystem::Galileo,
            GnssSystem::BeiDou,
        ] {
            for prn in [5u8, 30] {
                for &s in signal_order(sys, prn) {
                    let c = rp_code(sys, s).expect("a code");
                    assert_eq!(rp_signal(sys, c), Some(s));
                }
            }
        }
    }

    #[test]
    fn inference_by_doppler() {
        let mut inf = Inference::default();
        let g = GnssSystem::Gps;
        let d0 = 1234.567;
        let r2 = 1227.60 / 1575.42;
        let r5 = 1176.45 / 1575.42;
        // four signals: L5 rather than L2C as the fourth
        let v = inf
            .signals(g, 9, None, &[d0, d0, d0 * r2, d0 * r5])
            .expect("signals");
        assert_eq!(v, [*b"1C", *b"1W", *b"2W", *b"5X"]);
        let v = inf
            .signals(g, 10, None, &[d0, d0, d0 * r2, d0 * r2])
            .expect("signals");
        assert_eq!(v, [*b"1C", *b"1W", *b"2W", *b"2X"]);
        // too small a Doppler: the last list of the satellite is kept
        let v = inf
            .signals(g, 9, None, &[0.5, 0.5, 0.39, 0.37])
            .expect("signals");
        assert_eq!(v, [*b"1C", *b"1W", *b"2W", *b"5X"]);
        // one signal, and all five
        assert_eq!(inf.signals(g, 3, None, &[100.0]).expect("one"), [*b"1C"]);
        assert_eq!(
            inf.signals(g, 4, None, &[d0, d0, d0 * r2, d0 * r2, d0 * r5])
                .expect("five")
                .len(),
            5
        );
        // Galileo E5a / E5 / E5b
        let e = GnssSystem::Galileo;
        let f = |b: u8| carrier(e, b, None).expect("f") / carrier(e, 1, None).expect("f");
        let v = inf
            .signals(e, 11, None, &[-2000.0, -2000.0 * f(5), -2000.0 * f(8)])
            .expect("signals");
        assert_eq!(v, [*b"1X", *b"5X", *b"8X"]);
        // BeiDou: BDS-2 and BDS-3 lists
        let c = GnssSystem::BeiDou;
        let r7 = 1207.14 / 1561.098;
        assert_eq!(
            inf.signals(c, 9, None, &[500.0, 500.0 * r7]).expect("bds2"),
            [*b"2I", *b"7I"]
        );
        assert_eq!(
            inf.signals(c, 30, None, &[500.0, 500.0 * r7])
                .expect("bds3"),
            [*b"2I", *b"7D"]
        );
        assert!(inf.signals(GnssSystem::Qzss, 1, None, &[1.0]).is_none());
    }
}
