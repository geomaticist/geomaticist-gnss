//! Layer 3 of a packed log: the measurement epochs.
//!
//! A full epoch (packet 0x26) lists every satellite and signal. An epoch packet 0x25 holds the same content as
//! differences from the epoch before it, slot by slot (the n-th satellite of the list against the n-th of the
//! list before, whichever satellite that was): ranges, phases and the clock as second differences in
//! sign-and-magnitude bytes of a size named per field, or as the value itself when the size is the full one;
//! small signed bit fields for the counters, their most negative value meaning "the value follows".
//!
//! Its two check bytes are the check of the epoch once rebuilt, so every epoch decoded from differences is
//! verified against what the receiver had.

use std::collections::HashMap;

use super::packet::{record_check, sum16};

const WEEK_MS: u32 = 604_800_000;

/// One signal of a satellite. The range of the first signal is in 1/128 m; the others hold their range as an
/// offset from it in 1/256 m. Phases are in 1/32768 cycle.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Sig {
    /// Length of the signal block in a full epoch (0x11 for the first signal, 0x0F for the others).
    pub len: u8,
    /// Band (0 L1, 1 L2, 2 L5).
    pub band: u8,
    /// Code.
    pub code: u8,
    /// Signal-to-noise ratio, 0.1 dB-Hz.
    pub snr: u16,
    /// Range (first signal) or range offset (the others).
    pub rng: i64,
    /// Phase.
    pub ph: i64,
    /// Slip counter.
    pub slip: u8,
    /// Flags.
    pub fl: u8,
    /// The last difference of `rng`.
    pub d_rng: i64,
    /// The last difference of `ph`.
    pub d_ph: i64,
}

/// One satellite (one slot of the list).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Sat {
    /// PRN / slot number.
    pub prn: u8,
    /// System (0 GPS, 2 GLONASS).
    pub sys: u8,
    /// Receiver channel (GPS) or frequency channel (GLONASS, signed).
    pub ch: u8,
    /// Elevation, degrees.
    pub elev: u8,
    /// Azimuth, 2 degrees.
    pub az: u8,
    /// The last byte of the satellite block (2 in every sample).
    pub b7: u8,
    /// Signals.
    pub sigs: Vec<Sig>,
}

/// The measurements of one epoch, and what the next packet of differences builds on.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Meas {
    /// GPS week.
    pub week: u16,
    /// Time of week, ms.
    pub tow: u32,
    /// Receiver clock offset, 2^-19 ms.
    pub clk: i32,
    /// The last difference of `clk`.
    pub d_clk: i32,
    /// Header byte (0x12 in every sample).
    pub b12: u8,
    /// Clock drift, 2^-22 (unit not known).
    pub drift: i32,
    /// Header byte (0x60 in every sample).
    pub b60: u8,
    /// Step of the time of week, ms.
    pub dt: u16,
    /// Satellites.
    pub sats: Vec<Sat>,
}

struct Cur<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Cur<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        let s = self
            .b
            .get(self.i..self.i + n)
            .ok_or_else(|| "the epoch packet is too short".to_owned())?;
        self.i += n;
        Ok(s)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }

    /// Big-endian unsigned.
    fn uint(&mut self, n: usize) -> Result<i64, String> {
        Ok(self
            .take(n)?
            .iter()
            .fold(0i64, |v, &x| v << 8 | i64::from(x)))
    }

    /// Big-endian two's complement.
    fn int(&mut self, n: usize) -> Result<i64, String> {
        Ok(wrap(self.uint(n)?, 8 * n as u32))
    }

    /// Sign and magnitude: the first bit is the sign.
    fn sign_mag(&mut self, n: usize) -> Result<i64, String> {
        if n == 0 {
            return Ok(0);
        }
        let v = self.uint(n)?;
        let sign = 1i64 << (8 * n - 1);
        Ok(if v & sign != 0 { -(v & (sign - 1)) } else { v })
    }

    fn done(&self) -> Result<(), String> {
        if self.i == self.b.len() {
            Ok(())
        } else {
            Err(format!(
                "{} bytes of the epoch packet are left over",
                self.b.len() - self.i
            ))
        }
    }
}

/// `v` as a two's complement number of `bits` bits.
fn wrap(v: i64, bits: u32) -> i64 {
    let s = 64 - bits;
    (v << s) >> s
}

/// A signed bit field of `bits` bits whose most negative value means "the value follows".
fn small(field: u8, bits: u32) -> Option<i64> {
    let v = wrap(i64::from(field), bits);
    (v != -(1 << (bits - 1))).then_some(v)
}

impl Meas {
    /// A full epoch (the payload of a packet 0x26).
    pub(crate) fn full(payload: &[u8], dt: u16) -> Result<Self, String> {
        let mut c = Cur { b: payload, i: 0 };
        let hl = c.u8()?;
        if hl != 0x10 {
            return Err(format!("an epoch header of {hl} bytes"));
        }
        let mut m = Meas {
            week: c.uint(2)? as u16,
            tow: c.uint(4)? as u32,
            clk: c.int(3)? as i32,
            dt,
            ..Meas::default()
        };
        let nsat = c.u8()?;
        m.b12 = c.u8()?;
        m.drift = c.int(3)? as i32;
        m.b60 = c.u8()?;
        for _ in 0..nsat {
            let sl = c.u8()?;
            if sl != 8 {
                return Err(format!("a satellite block of {sl} bytes"));
            }
            let s = c.take(7)?;
            let mut sat = Sat {
                prn: s[0],
                sys: s[1],
                ch: s[2],
                elev: s[4],
                az: s[5],
                b7: s[6],
                sigs: Vec::with_capacity(usize::from(s[3])),
            };
            for j in 0..s[3] {
                let len = c.u8()?;
                if len != if j == 0 { 0x11 } else { 0x0f } {
                    return Err(format!("a signal block of {len} bytes"));
                }
                let (band, code) = (c.u8()?, c.u8()?);
                let snr = c.uint(2)? as u16;
                let rng = if j == 0 { c.uint(4)? } else { c.int(2)? };
                sat.sigs.push(Sig {
                    len,
                    band,
                    code,
                    snr,
                    rng,
                    ph: c.int(6)?,
                    slip: c.u8()?,
                    fl: c.u8()?,
                    d_rng: 0,
                    d_ph: 0,
                });
            }
            m.sats.push(sat);
        }
        c.done()?;
        Ok(m)
    }

    /// The epoch after this one, from a packet 0x25 (its payload and its two check bytes). On an error `self`
    /// is no longer an epoch to build on.
    pub(crate) fn apply_delta(&mut self, payload: &[u8], check: u16) -> Result<(), String> {
        let mut c = Cur { b: payload, i: 0 };
        let f = c.u8()?;
        if f & 0x20 != 0 {
            return Err(format!("epoch flags {f:02x}"));
        }
        if f & 0x80 != 0 {
            let b = c.take(8)?;
            if b[0] != 0x10 {
                return Err(format!(
                    "an epoch header block that starts with {:02x}",
                    b[0]
                ));
            }
            self.week = u16::from(b[1]) << 8 | u16::from(b[2]);
            self.b12 = b[3];
            self.drift = wrap(
                i64::from(b[4]) << 16 | i64::from(b[5]) << 8 | i64::from(b[6]),
                24,
            ) as i32;
            self.b60 = b[7];
        }
        if f & 0x40 != 0 {
            self.dt = c.uint(2)? as u16;
        }
        let size = usize::from(f >> 3 & 3);
        let d = if size == 3 {
            c.int(3)? - i64::from(self.clk)
        } else {
            i64::from(self.d_clk) + c.sign_mag(size)?
        };
        self.d_clk = d as i32;
        self.clk = wrap(i64::from(self.clk) + d, 24) as i32;
        let nsat = match small(f & 7, 3) {
            Some(d) => usize::try_from(self.sats.len() as i64 + d)
                .map_err(|_| "fewer than no satellites".to_owned())?,
            None => usize::from(c.u8()?),
        };
        self.tow += u32::from(self.dt);
        if self.tow >= WEEK_MS {
            self.tow -= WEEK_MS;
            self.week = self.week.wrapping_add(1);
        }
        self.sats.truncate(nsat);
        for k in 0..nsat {
            let sp = c.u8()?;
            if k == self.sats.len() {
                if sp & 0x81 != 0x81 {
                    return Err(format!("a new satellite slot described by {sp:02x}"));
                }
                self.sats.push(Sat::default());
            }
            let sat = &mut self.sats[k];
            if sp & 0x80 != 0 {
                let sl = c.u8()?;
                if sl != 8 {
                    return Err(format!("a satellite block of {sl} bytes"));
                }
                sat.b7 = c.u8()?;
            }
            if sp & 1 != 0 {
                let s = c.take(4)?;
                (sat.prn, sat.sys, sat.ch) = (s[0], s[1], s[2]);
                sat.sigs.resize_with(usize::from(s[3]), Sig::default);
            }
            sat.elev = match small(sp >> 4 & 7, 3) {
                Some(d) => sat.elev.wrapping_add(d as u8),
                None => c.u8()?,
            };
            sat.az = match small(sp >> 1 & 7, 3) {
                Some(d) => sat.az.wrapping_add(d as u8),
                None => c.u8()?,
            };
            for j in 0..sat.sigs.len() {
                let (p, sizes) = (c.u8()?, c.u8()?);
                let x = &mut sat.sigs[j];
                if p & 0x80 != 0 {
                    let s = c.take(4)?;
                    (x.len, x.band, x.code, x.fl) = (s[0], s[1], s[2], s[3]);
                }
                x.slip = match small(p & 3, 2) {
                    Some(d) => x.slip.wrapping_add(d as u8),
                    None => c.u8()?,
                };
                x.snr = match small(p >> 2 & 31, 5) {
                    Some(d) => x.snr.wrapping_add(d as u16),
                    None => c.uint(2)? as u16,
                };
                let (rn, pn) = (usize::from(sizes >> 4), usize::from(sizes & 15));
                // the range: its own (first signal), or its offset from the first signal's
                let (full, most, bits) = if j == 0 { (4, 3, 32) } else { (2, 1, 16) };
                let d = if rn == full {
                    let v = if j == 0 { c.uint(4)? } else { c.int(2)? };
                    v - x.rng
                } else if rn <= most {
                    x.d_rng + c.sign_mag(rn)?
                } else {
                    return Err(format!("a range of {rn} bytes (signal {j})"));
                };
                x.d_rng = d;
                x.rng = if j == 0 {
                    (x.rng + d) & 0xffff_ffff
                } else {
                    wrap(x.rng + d, bits)
                };
                let d = match pn {
                    6 => c.int(6)? - x.ph,
                    0..=5 => x.d_ph + c.sign_mag(pn)?,
                    _ => return Err(format!("a phase of {pn} bytes")),
                };
                x.d_ph = d;
                x.ph = wrap(x.ph + d, 48);
            }
        }
        c.done()?;
        if self.check() != check {
            return Err("the rebuilt epoch doesn't pass the packet's check".to_owned());
        }
        Ok(())
    }

    /// The payload of the full epoch packet (0x26) that holds this epoch.
    pub(crate) fn full_payload(&self) -> Vec<u8> {
        let mut b = Vec::with_capacity(16 + self.sats.len() * 40);
        b.push(0x10);
        b.extend_from_slice(&self.week.to_be_bytes());
        b.extend_from_slice(&self.tow.to_be_bytes());
        b.extend_from_slice(&self.clk.to_be_bytes()[1..]);
        b.extend_from_slice(&[self.sats.len() as u8, self.b12]);
        b.extend_from_slice(&self.drift.to_be_bytes()[1..]);
        b.push(self.b60);
        for s in &self.sats {
            b.extend_from_slice(&[
                8,
                s.prn,
                s.sys,
                s.ch,
                s.sigs.len() as u8,
                s.elev,
                s.az,
                s.b7,
            ]);
            for (j, x) in s.sigs.iter().enumerate() {
                b.extend_from_slice(&[x.len, x.band, x.code]);
                b.extend_from_slice(&x.snr.to_be_bytes());
                if j == 0 {
                    b.extend_from_slice(&(x.rng as u32).to_be_bytes());
                } else {
                    b.extend_from_slice(&(x.rng as i16).to_be_bytes());
                }
                b.extend_from_slice(&x.ph.to_be_bytes()[2..]);
                b.extend_from_slice(&[x.slip, x.fl]);
            }
        }
        b
    }

    /// The check bytes a packet of differences carries for this epoch.
    pub(crate) fn check(&self) -> u16 {
        let body = self.full_payload();
        record_check(sum16(&body), body.len() + 6, 0x1B)
    }
}

/// Carrier frequency of a signal, Hz (0 when not known).
fn frequency(sys: u8, ch: u8, band: u8) -> f64 {
    let k = f64::from(ch as i8);
    match (sys, band) {
        (2, 0) => 1602.0e6 + k * 562_500.0,
        (2, 1) => 1246.0e6 + k * 437_500.0,
        (_, 0) => 1575.42e6,
        (_, 1) => 1227.6e6,
        (_, 2) => 1176.45e6,
        _ => 0.0,
    }
}

#[derive(Clone, Copy)]
struct Track {
    ph: i64,
    fl: u8,
    doppler: Option<f32>,
    doppler_before: Option<f32>,
}

/// The Doppler values the unpacked epoch records carry. They are not in the file: they are derived from the
/// epoch before, and flag 0x04 is set on a signal that has one (the RINEX writer reads that flag).
#[derive(Default)]
pub(crate) struct Doppler {
    prev: HashMap<(u8, u8, u8, u8), Track>,
    tow: Option<u32>,
    clk: i32,
}

impl Doppler {
    /// The Doppler of every signal of `m`, in its order.
    pub(crate) fn next(&mut self, m: &Meas) -> Vec<Option<f32>> {
        // the time of week alone is differenced (so the values are off across a week rollover)
        let dt = self
            .tow
            .map_or(0.0, |t| (f64::from(m.tow) - f64::from(t)) / 1000.0);
        // a whole-millisecond clock jump moves every phase by its frequency times the jump
        let jump = (f64::from(m.clk - self.clk) * 2f64.powi(-19)).round() * 1e-3;
        let mut cur = HashMap::with_capacity(self.prev.len());
        let mut out = Vec::new();
        for s in &m.sats {
            for x in &s.sigs {
                let key = (s.sys, s.prn, x.band, x.code);
                let p = self.prev.get(&key);
                let doppler = p.and_then(|p| {
                    if x.fl & 0x08 != 0 {
                        // a slip: the two Doppler values before, continued
                        Some(2.0 * p.doppler? - p.doppler_before?)
                    } else if dt - jump == 0.0 {
                        None
                    } else {
                        // half a cycle less when the phase was not valid before (no flag 0x20) and the
                        // signal now carries flag 0x40
                        let half = if p.fl & 0x20 == 0 && x.fl & 0x40 != 0 {
                            -0.5
                        } else {
                            0.0
                        };
                        let cycles = (x.ph - p.ph) as f64 / 32768.0 + half;
                        Some(
                            ((cycles + frequency(s.sys, s.ch, x.band) * jump) / (dt - jump)) as f32,
                        )
                    }
                });
                cur.insert(
                    key,
                    Track {
                        ph: x.ph,
                        fl: x.fl,
                        doppler,
                        doppler_before: p.and_then(|p| p.doppler),
                    },
                );
                out.push(doppler);
            }
        }
        self.prev = cur;
        self.tow = Some(m.tow);
        self.clk = m.clk;
        out
    }
}

/// The unpacked record of an epoch (type 0x1B, the .TGD layout): header, then per
/// satellite a 10-byte block and one 27-byte block per signal, little-endian.
pub(crate) fn record(m: &Meas, doppler: &[Option<f32>]) -> Vec<u8> {
    let mut b = vec![0x74, 0, 0, 0x1b, 0x16];
    b.extend_from_slice(&m.week.to_le_bytes());
    b.extend_from_slice(&f64::from(m.tow).to_le_bytes());
    b.extend_from_slice(&((f64::from(m.clk) * 2f64.powi(-19)) as f32).to_le_bytes());
    b.extend_from_slice(&[m.sats.len() as u8, m.b12]);
    b.extend_from_slice(&((f64::from(m.drift) * 2f64.powi(-22)) as f32).to_le_bytes());
    b.push(m.b60);
    let mut dop = doppler.iter();
    for s in &m.sats {
        b.extend_from_slice(&[0x0a, s.prn, s.sys, s.ch, s.sigs.len() as u8, s.elev, 0]);
        b.extend_from_slice(&(i16::from(s.az) * 2).to_le_bytes());
        b.push(s.b7);
        let first = s.sigs.first().map_or(0, |x| x.rng);
        for (j, x) in s.sigs.iter().enumerate() {
            let range = if j == 0 {
                x.rng as f64 / 128.0
            } else {
                (first * 2 + x.rng) as f64 / 256.0
            };
            let d = dop.next().copied().flatten();
            b.extend_from_slice(&[0x1b, x.band, x.code]);
            b.extend_from_slice(&x.snr.to_le_bytes());
            b.extend_from_slice(&range.to_le_bytes());
            b.extend_from_slice(&(x.ph as f64 / 32768.0).to_le_bytes());
            b.extend_from_slice(&d.unwrap_or(0.0).to_le_bytes());
            b.extend_from_slice(&[x.slip, x.fl | if d.is_some() { 0x04 } else { 0 }]);
        }
    }
    let len = b.len();
    b[1] = (len & 0xff) as u8;
    b[2] = (len >> 8) as u8;
    b
}

/// The packet of differences (payload and check bytes) that takes a receiver from `prev` to `next`, written
/// the way the receiver does: the smallest size that holds each difference. `next` gets the differences its
/// successor builds on.
#[cfg(test)]
pub(crate) fn encode_delta(prev: &Meas, next: &mut Meas) -> (Vec<u8>, u16) {
    fn sign_mag(v: i64, n: usize) -> Vec<u8> {
        let m = v.unsigned_abs() | if v < 0 { 1u64 << (8 * n - 1) } else { 0 };
        m.to_be_bytes()[8 - n..].to_vec()
    }
    /// Bytes a second difference takes, `None` past `most`.
    fn size(v: i64, most: usize) -> Option<usize> {
        (0..=most).find(|&n| v == 0 || (n > 0 && v.unsigned_abs() < 1u64 << (8 * n - 1)))
    }
    fn field(d: i64, bits: u32) -> (u8, bool) {
        let lim = 1i64 << (bits - 1);
        if d > -lim && d < lim {
            ((d & ((1 << bits) - 1)) as u8, false)
        } else {
            (lim as u8, true)
        }
    }
    let mut b = vec![0u8];
    let mut f = 0u8;
    if (prev.week, prev.b12, prev.drift, prev.b60) != (next.week, next.b12, next.drift, next.b60)
        && !(prev.tow + u32::from(next.dt) >= WEEK_MS
            && next.week == prev.week.wrapping_add(1)
            && (prev.b12, prev.drift, prev.b60) == (next.b12, next.drift, next.b60))
    {
        f |= 0x80;
        b.push(0x10);
        b.extend_from_slice(&next.week.to_be_bytes());
        b.push(next.b12);
        b.extend_from_slice(&next.drift.to_be_bytes()[1..]);
        b.push(next.b60);
    }
    if next.dt != prev.dt {
        f |= 0x40;
        b.extend_from_slice(&next.dt.to_be_bytes());
    }
    let d = i64::from(next.clk) - i64::from(prev.clk);
    next.d_clk = d as i32;
    if let Some(n) = size(d - i64::from(prev.d_clk), 2) {
        f |= (n as u8) << 3;
        b.extend(sign_mag(d - i64::from(prev.d_clk), n));
    } else {
        f |= 3 << 3;
        b.extend_from_slice(&next.clk.to_be_bytes()[1..]);
    }
    let (bits, escaped) = field(next.sats.len() as i64 - prev.sats.len() as i64, 3);
    f |= bits;
    if escaped {
        b.push(next.sats.len() as u8);
    }
    b[0] = f;
    let empty = Sat::default();
    for (k, sat) in next.sats.iter_mut().enumerate() {
        let old = prev.sats.get(k).unwrap_or(&empty);
        let fresh = k >= prev.sats.len();
        let same = !fresh
            && (old.prn, old.sys, old.ch, old.sigs.len())
                == (sat.prn, sat.sys, sat.ch, sat.sigs.len());
        let (e, e_esc) = field(i64::from(sat.elev) - i64::from(old.elev), 3);
        let (a, a_esc) = field(i64::from(sat.az) - i64::from(old.az), 3);
        b.push(if fresh { 0x80 } else { 0 } | e << 4 | a << 1 | u8::from(!same));
        if fresh {
            b.extend_from_slice(&[8, sat.b7]);
        }
        if !same {
            b.extend_from_slice(&[sat.prn, sat.sys, sat.ch, sat.sigs.len() as u8]);
        }
        if e_esc {
            b.push(sat.elev);
        }
        if a_esc {
            b.push(sat.az);
        }
        let none = Sig::default();
        for (j, x) in sat.sigs.iter_mut().enumerate() {
            let o = old.sigs.get(j).unwrap_or(&none);
            let extra = (o.len, o.band, o.code, o.fl) != (x.len, x.band, x.code, x.fl);
            let (slip, slip_esc) = field(i64::from(x.slip.wrapping_sub(o.slip) as i8), 2);
            let (snr, snr_esc) = field(i64::from(x.snr) - i64::from(o.snr), 5);
            let (full, most) = if j == 0 { (4, 3) } else { (2, 1) };
            x.d_rng = x.rng - o.rng;
            x.d_ph = x.ph - o.ph;
            let rn = size(x.d_rng - o.d_rng, most);
            let pn = size(x.d_ph - o.d_ph, 5);
            b.push(u8::from(extra) << 7 | snr << 2 | slip);
            b.push((rn.unwrap_or(full) as u8) << 4 | pn.unwrap_or(6) as u8);
            if extra {
                b.extend_from_slice(&[x.len, x.band, x.code, x.fl]);
            }
            if slip_esc {
                b.push(x.slip);
            }
            if snr_esc {
                b.extend_from_slice(&x.snr.to_be_bytes());
            }
            match rn {
                Some(n) => b.extend(sign_mag(x.d_rng - o.d_rng, n)),
                None if j == 0 => b.extend_from_slice(&(x.rng as u32).to_be_bytes()),
                None => b.extend_from_slice(&(x.rng as i16).to_be_bytes()),
            }
            match pn {
                Some(n) => b.extend(sign_mag(x.d_ph - o.d_ph, n)),
                None => b.extend_from_slice(&x.ph.to_be_bytes()[2..]),
            }
        }
    }
    (b, next.check())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sig(j: usize, band: u8, code: u8, snr: u16, rng: i64, ph: i64) -> Sig {
        Sig {
            len: if j == 0 { 0x11 } else { 0x0f },
            band,
            code,
            snr,
            rng,
            ph,
            slip: 7,
            fl: 0x2b,
            d_rng: 0,
            d_ph: 0,
        }
    }

    fn sat(prn: u8, sys: u8, ch: u8, sigs: Vec<Sig>) -> Sat {
        Sat {
            prn,
            sys,
            ch,
            elev: 40,
            az: 100,
            b7: 2,
            sigs,
        }
    }

    fn first() -> Meas {
        Meas {
            week: 2000,
            tow: 100_000_000,
            clk: 20,
            b12: 0x12,
            drift: 7,
            b60: 0x60,
            sats: vec![
                sat(
                    6,
                    0,
                    0,
                    vec![
                        sig(0, 0, 0, 450, 2_560_000_128, -3_276_800_000_000),
                        sig(1, 1, 2, 350, 384, -2_555_904_000_000),
                    ],
                ),
                sat(
                    9,
                    2,
                    0xfe,
                    vec![
                        sig(0, 0, 0, 440, 2_688_000_000, -3_604_480_000_000),
                        sig(1, 1, 1, 360, -2, -2_785_280_000_000),
                    ],
                ),
            ],
            ..Meas::default()
        }
    }

    /// Decodes what `encode_delta` writes and checks it gives `next` back.
    fn step(state: &mut Meas, next: &mut Meas) -> Vec<u8> {
        let (payload, check) = encode_delta(state, next);
        state.apply_delta(&payload, check).unwrap();
        assert_eq!(state, next);
        payload
    }

    #[test]
    fn a_full_epoch_round_trips() {
        let m = first();
        let payload = m.full_payload();
        assert_eq!(payload.len(), 16 + 2 * (8 + 17 + 15));
        assert_eq!(Meas::full(&payload, 0).unwrap(), m);
        // a cut or padded packet is refused
        assert!(Meas::full(&payload[..payload.len() - 1], 0).is_err());
        let mut long = payload.clone();
        long.push(0);
        assert!(Meas::full(&long, 0).is_err());
    }

    #[test]
    fn the_unpacked_record() {
        let m = first();
        let r = record(&m, &[None, None, Some(-602.5), None]);
        assert_eq!(r.len(), 4 + 22 + 2 * (10 + 2 * 27));
        assert_eq!(&r[..4], [0x74, 154, 0, 0x1b]);
        assert_eq!(u16::from_le_bytes([r[5], r[6]]), 2000);
        let f = |i: usize| f64::from_le_bytes(r[i..i + 8].try_into().unwrap());
        assert_eq!(f(7), 100_000_000.0);
        // satellite 1: block at 26, signals at 36 and 63
        assert_eq!(&r[26..36], [0x0a, 6, 0, 0, 2, 40, 0, 200, 0, 2]);
        assert_eq!(f(36 + 5), 20_000_001.0);
        assert_eq!(f(36 + 13), -100_000_000.0);
        assert_eq!(f(63 + 5), 20_000_002.5); // (2 * 2560000128 + 384) / 256
        assert_eq!(r[36 + 26], 0x2b);
        // satellite 2: the Doppler and its flag on the first signal
        let s2 = 26 + 10 + 54 + 10;
        assert_eq!(r[s2 + 21..s2 + 25], (-602.5f32).to_le_bytes());
        assert_eq!(r[s2 + 26], 0x2f);
        assert_eq!(f(s2 + 27 + 5), 20_999_999.992_187_5); // a negative offset
    }

    #[test]
    fn differences_follow_a_moving_constellation() {
        let mut state = first();
        // 1: first differences after a full epoch, a flags change, the time step
        let mut n = state.clone();
        n.tow += 15_000;
        n.dt = 15_000;
        n.clk = -980;
        for s in &mut n.sats {
            for x in &mut s.sigs {
                x.fl = 0x23;
                x.snr -= 9;
                x.ph += 250_000_000;
            }
            s.sigs[0].rng += 200_000;
            s.sigs[1].rng += 168;
        }
        let p = step(&mut state, &mut n);
        // the time step, then the clock difference in two bytes, sign and magnitude (-1000)
        assert_eq!(p[..5], [0x40 | 2 << 3, 0x3a, 0x98, 0x83, 0xe8]);
        // 2: second differences in the small sizes
        let mut n = state.clone();
        n.tow += 15_000;
        n.clk += -1000 - 5;
        for s in &mut n.sats {
            s.sigs[0].rng += 200_000 + 3000;
            s.sigs[0].ph += 250_000_000 - 4_000_000;
            s.sigs[1].ph += 250_000_000 + 100;
            s.sigs[1].snr += 3;
        }
        n.sats[1].elev += 1;
        n.sats[1].az -= 1;
        let p = step(&mut state, &mut n);
        assert_eq!(p[..2], [1 << 3, 0x85]); // a one-byte clock difference, sign and magnitude
        assert_eq!(p[2], 0); // the first satellite did not move
        assert_eq!(p[4], 0x23); // two range bytes, three phase bytes
        // 3: a satellite joins in front (every slot is described again against what it held), a slip, a
        //    signal-to-noise jump, a third signal
        let mut n = state.clone();
        n.tow += 15_000;
        let mut new = sat(
            17,
            0,
            5,
            vec![
                sig(0, 0, 0, 300, 3_000_000_000, -4_000_000_000_000),
                sig(1, 1, 2, 250, 900, -3_100_000_000_000),
                sig(2, 1, 5, 380, 940, -3_100_000_000_512),
            ],
        );
        new.elev = 5;
        new.az = 170;
        n.sats.insert(0, new);
        n.sats[1].sigs[0].slip += 1;
        n.sats[1].sigs[0].fl = 0x2b;
        n.sats[2].sigs[1].snr += 40;
        n.sats[2].sigs[1].slip = n.sats[2].sigs[1].slip.wrapping_sub(1);
        n.clk += 1 << 19; // a millisecond jump: the clock itself, in three bytes
        let p = step(&mut state, &mut n);
        assert_eq!(p[0], 3 << 3 | 1); // and one satellite more
        // 4: the list shrinks from the front, far more than the small field holds; a week rollover
        let mut n = state.clone();
        n.tow = WEEK_MS - 15_000;
        state.tow = WEEK_MS - 30_000;
        for _ in 0..5 {
            let mut s = n.sats[2].clone();
            s.prn += 1 + n.sats.len() as u8;
            n.sats.push(s);
        }
        step(&mut state, &mut n);
        let mut n = state.clone();
        n.tow = 0;
        n.week += 1;
        n.sats.drain(..6);
        let p = step(&mut state, &mut n);
        assert_eq!(p[0] & 7, 4); // the count itself follows
        assert_eq!((state.week, state.tow, state.sats.len()), (2001, 0, 2));
        // 5: nothing changes at all
        let mut n = state.clone();
        n.tow += 15_000;
        for s in &mut n.sats {
            for x in &mut s.sigs {
                x.rng += x.d_rng;
                x.ph += x.d_ph;
            }
        }
        n.clk += n.d_clk;
        let p = step(&mut state, &mut n);
        assert_eq!(
            p.len(),
            1 + n.sats.iter().map(|s| 1 + 2 * s.sigs.len()).sum::<usize>()
        );
    }

    #[test]
    fn damage_is_caught_by_the_check() {
        let mut state = first();
        let mut n = state.clone();
        n.tow += 15_000;
        n.dt = 15_000;
        n.sats[0].sigs[0].rng += 1000;
        let (payload, check) = encode_delta(&state, &mut n);
        let mut ok = state.clone();
        ok.apply_delta(&payload, check).unwrap();
        // any one bit flipped: refused, by the layout or by the check
        for bit in 0..payload.len() * 8 {
            let mut bad = payload.clone();
            bad[bit / 8] ^= 0x80 >> (bit % 8);
            assert!(state.clone().apply_delta(&bad, check).is_err(), "bit {bit}");
        }
        assert!(state.apply_delta(&payload, check ^ 1).is_err());
    }

    #[test]
    fn doppler_rules() {
        let mut d = Doppler::default();
        let mut m = first();
        assert_eq!(d.next(&m), [None; 4]);
        // plain: the phase step over the time step
        m.tow += 15_000;
        for s in &mut m.sats {
            for x in &mut s.sigs {
                x.fl = 0x23;
                x.ph -= 15 * 1000 * 32768;
            }
        }
        assert_eq!(d.next(&m), [Some(-1000.0); 4]);
        // a clock jump of one millisecond: the phases move by a millisecond of carrier, the step lasts a
        // millisecond less
        m.tow += 15_000;
        m.clk += 1 << 19;
        for s in &mut m.sats {
            for x in &mut s.sigs {
                let f = frequency(s.sys, s.ch, x.band);
                x.ph += ((-1000.0 * 14.999 - f * 1e-3) * 32768.0).round() as i64;
            }
        }
        let v = d.next(&m);
        assert!(
            v.iter().all(|x| (x.unwrap() + 1000.0).abs() < 1e-3),
            "{v:?}"
        );
        assert_eq!(frequency(2, 0xfe, 0), 1_600_875_000.0);
        assert_eq!(frequency(2, 0xfe, 1), 1_245_125_000.0);
        // a slip: the two values before, continued; with one value before, none
        m.tow += 15_000;
        m.sats[0].sigs[0].fl = 0x2b;
        m.sats[0].sigs[0].ph = 12345;
        let w = d.next(&m);
        assert_eq!(w[0], Some(2.0 * v[0].unwrap() - -1000.0));
        m.sats.remove(1);
        m.tow += 15_000;
        d.next(&m);
        m.sats.push(first().sats.remove(1));
        m.tow += 15_000;
        assert_eq!(d.next(&m)[2..], [None, None]); // back after an epoch away: no epoch before
        m.tow += 15_000;
        m.sats[1].sigs[0].fl = 0x2b;
        assert_eq!(d.next(&m)[2], None);
        // the same time twice: no value
        m.sats[0].sigs[0].fl = 0x23;
        m.sats[1].sigs[0].fl = 0x23;
        assert!(d.next(&m).iter().all(Option::is_none));
    }
}
