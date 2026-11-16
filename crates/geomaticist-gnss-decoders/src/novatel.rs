//! NovAtel OEM4 / OEMV / OEM6 / OEM7 binary logs: the range logs RANGECMPB (id 140) and RANGEB (id 43)
//! into the engine's observations.
//!
//! The frame (sync `AA 44 12`, a header whose length is its fourth byte, the message, a CRC-32) and the
//! two layouts follow NovAtel's public firmware reference; the use of the channel tracking status (which
//! bits mean code lock, phase lock, parity known and half-cycle added, the system and signal-type fields),
//! the rebuilding of the compressed carrier phase from its rolled-over value, and the lock-time rule for
//! loss of lock follow RTKLIB's `novatel.c` (T. Takasu, BSD-2-Clause). No RTKLIB code is copied.
//!
//! Both logs carry the GPS week, so no outside date is needed. Older OEM3 logs (sync `AA 44 11`) are
//! recognised but not read. Navigation logs are not decoded.

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::consts::{CLIGHT, WEEK};
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};

use crate::epochs::{Epochs, stem_and_name, write_obs};
use crate::{Converted, DecodeError};

/// RANGECMPB phase rolls over at this many cycles.
const MAXVAL: f64 = 8_388_608.0;

/// NovAtel's CRC-32 (reflected polynomial 0xEDB88320, starting at 0, no final xor).
#[must_use]
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0u32;
    for &b in data {
        crc ^= u32::from(b);
        for _ in 0..8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xEDB8_8320
            } else {
                crc >> 1
            };
        }
    }
    crc
}

/// A binary OEM4-family message: id, GPS week, time status and milliseconds of the week, body.
pub(crate) struct Msg<'a> {
    pub id: u16,
    pub week: u16,
    pub status: u8,
    pub ms: u32,
    pub body: &'a [u8],
}

/// The OEM4-family binary messages of a stream with a good CRC, and the number of bad ones.
pub(crate) fn messages(data: &[u8]) -> (Vec<Msg<'_>>, usize) {
    let mut out = Vec::new();
    let mut bad = 0;
    let mut i = 0;
    while i + 28 <= data.len() {
        if data[i..i + 3] != [0xAA, 0x44, 0x12] {
            i += 1;
            continue;
        }
        let h = &data[i..];
        let hlen = usize::from(h[3]);
        let mlen = usize::from(u16::from_le_bytes([h[8], h[9]]));
        let end = hlen + mlen + 4;
        if hlen < 28 || end > h.len() {
            i += 1;
            continue;
        }
        let crc = u32::from_le_bytes([h[end - 4], h[end - 3], h[end - 2], h[end - 1]]);
        if crc32(&h[..end - 4]) != crc {
            bad += 1;
            i += 1;
            continue;
        }
        if (h[6] >> 4).trailing_zeros() >= 2 {
            out.push(Msg {
                id: u16::from_le_bytes([h[4], h[5]]),
                status: h[13],
                week: u16::from_le_bytes([h[14], h[15]]),
                ms: u32::from_le_bytes([h[16], h[17], h[18], h[19]]),
                body: &h[hlen..hlen + mlen],
            });
        }
        i += end;
    }
    (out, bad)
}

fn u2(d: &[u8], o: usize) -> u16 {
    u16::from_le_bytes([d[o], d[o + 1]])
}
fn u4(d: &[u8], o: usize) -> u32 {
    u32::from_le_bytes([d[o], d[o + 1], d[o + 2], d[o + 3]])
}
fn r4(d: &[u8], o: usize) -> f64 {
    f64::from(f32::from_bits(u4(d, o)))
}
fn r8(d: &[u8], o: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[o..o + 8]);
    f64::from_le_bytes(b)
}

/// What a channel tracking status says.
struct Track {
    sys: GnssSystem,
    tail: &'static str,
    plock: bool,
    parity: bool,
    clock: bool,
}

fn track(stat: u32) -> Option<Track> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Qzss, Sbas};
    let sys = match (stat >> 16) & 7 {
        0 => Gps,
        1 => Glonass,
        2 => Sbas,
        3 => Galileo,
        4 => BeiDou,
        5 => Qzss,
        _ => return None,
    };
    let sig = (stat >> 21) & 0x1F;
    let tail = match (sys, sig) {
        (Gps | Glonass | Qzss | Sbas, 0) | (Galileo, 2) => "1C",
        (Gps | Glonass, 5) => "2P",
        (Gps, 9) => "2W",
        (Gps | Qzss, 14) | (Galileo, 12) => "5Q",
        (Gps | Qzss, 16) => "1L",
        (Gps | Qzss, 17) => "2S",
        (Glonass, 1) => "2C",
        (Glonass, 6) => "3Q",
        (Galileo, 6) => "6B",
        (Galileo, 7) => "6C",
        (Galileo, 17) => "7Q",
        (Galileo, 20) => "8Q",
        (Qzss, 27) => "6L",
        (BeiDou, 0 | 4) => "2I",
        (BeiDou, 1 | 5) => "7I",
        (BeiDou, 2 | 6) => "6I",
        (BeiDou, 7) => "1P",
        (BeiDou, 9) => "5P",
        (BeiDou, 11) => "7D",
        (Sbas, 6) => "5I",
        _ => return None,
    };
    Some(Track {
        sys,
        tail,
        plock: stat & (1 << 10) != 0,
        parity: stat & (1 << 11) != 0,
        clock: stat & (1 << 12) != 0,
    })
}

/// The satellite of a NovAtel PRN field in a system.
fn sat_of(sys: GnssSystem, prn: u32) -> Option<Sat> {
    let prn = match sys {
        GnssSystem::Glonass => prn.checked_sub(37)?,
        GnssSystem::Sbas => prn.checked_sub(100)?,
        GnssSystem::Qzss if prn > 192 => prn - 192,
        _ => prn,
    };
    (1..100).contains(&prn).then(|| Sat::new(sys, prn as u8))
}

/// What [`decode`] found.
#[derive(Debug)]
pub struct NovAtel {
    /// The observations.
    pub obs: StationObs,
    /// GLONASS frequency channels from RANGEB.
    pub glonass_fcn: BTreeMap<Sat, i32>,
    /// Messages seen, by id.
    pub messages: BTreeMap<u16, usize>,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// Decodes a NovAtel binary log.
#[must_use]
pub fn decode(data: &[u8], name: &str) -> NovAtel {
    let (msgs, bad) = messages(data);
    let mut seen: BTreeMap<u16, usize> = BTreeMap::new();
    let mut epochs = Epochs::default();
    let mut fcn: BTreeMap<Sat, i32> = BTreeMap::new();
    // last epoch time and lock time per satellite and signal
    let mut last: BTreeMap<(Sat, &'static str), (f64, f64)> = BTreeMap::new();
    let mut unknown = 0usize;
    for m in msgs {
        *seen.entry(m.id).or_insert(0) += 1;
        if m.week == 0 || m.status == 20 || !matches!(m.id, 43 | 140) || m.body.len() < 4 {
            continue;
        }
        let t = f64::from(m.week) * WEEK + f64::from(m.ms) * 0.001;
        let n = u4(m.body, 0) as usize;
        let size = if m.id == 140 { 24 } else { 44 };
        if m.body.len() < 4 + n * size {
            continue;
        }
        for k in 0..n {
            let p = &m.body[4 + k * size..4 + (k + 1) * size];
            let stat = if m.id == 140 { u4(p, 0) } else { u4(p, 40) };
            let Some(tr) = track(stat) else {
                unknown += 1;
                continue;
            };
            let prn = if m.id == 140 {
                u32::from(p[17])
            } else {
                u32::from(u2(p, 0))
            };
            let Some(sat) = sat_of(tr.sys, prn) else {
                continue;
            };
            if tr.sys == GnssSystem::Glonass && !tr.parity {
                continue;
            }
            let band = tr.tail.as_bytes()[0] - b'0';
            let (psr, adr, dop, snr, lockt);
            if m.id == 140 {
                dop = f64::from(((u4(p, 4) & 0x0FFF_FFFF) << 4) as i32 >> 4) / 256.0;
                psr = f64::from(u4(p, 7) >> 4) / 128.0 + f64::from(p[11]) * 2_097_152.0;
                let raw = f64::from(u4(p, 12) as i32) / 256.0;
                let fr = freq(sat, band, Some(0)).unwrap_or(1.5e9);
                let rolls = (psr * fr / CLIGHT + raw) / MAXVAL;
                adr = -raw + MAXVAL * (rolls + if rolls <= 0.0 { -0.5 } else { 0.5 }).floor();
                lockt = f64::from(u4(p, 18) & 0x1F_FFFF) / 32.0;
                snr = f64::from((u2(p, 20) & 0x3FF) >> 5) + 20.0;
            } else {
                if tr.sys == GnssSystem::Glonass {
                    fcn.insert(sat, i32::from(u2(p, 2)) - 7);
                }
                psr = r8(p, 4);
                adr = -r8(p, 16);
                dop = r4(p, 28);
                snr = r4(p, 32);
                lockt = r4(p, 36);
            }
            let prev = last.insert((sat, tr.tail), (t, lockt));
            let slip =
                prev.is_some_and(|(t0, l0)| lockt < 65_535.968 && lockt - l0 + 0.05 <= t - t0);
            let c = epochs.cell(t, sat, tr.tail);
            if tr.clock {
                c.code = psr;
            }
            if tr.plock {
                c.phase = adr;
                c.doppler = dop;
                c.lli = i8::from(slip) | (i8::from(!tr.parity) << 1);
            }
            c.snr = snr;
        }
    }
    let mut warnings = Vec::new();
    if bad > 0 {
        warnings.push(format!(
            "{bad} NovAtel messages failed their checksum and were left out"
        ));
    }
    if unknown > 0 {
        warnings.push(format!(
            "{unknown} range observations of signals not known here were left out"
        ));
    }
    let obs = epochs.into_obs(name);
    fcn.retain(|s, _| obs.sats.contains_key(s));
    NovAtel {
        obs,
        glonass_fcn: fcn,
        messages: seen,
        warnings,
    }
}

/// Whether a stream reads as an OEM4-family binary log: a message with a good CRC in its head.
#[must_use]
pub fn recognises(head: &[u8]) -> bool {
    !messages(head).0.is_empty()
}

/// Whether a stream holds OEM3 binary messages (sync `AA 44 11`), which are not read here.
#[must_use]
pub fn is_oem3(head: &[u8]) -> bool {
    head.windows(3).filter(|w| *w == [0xAA, 0x44, 0x11]).count() >= 2
}

/// Decodes a NovAtel binary log into `outdir` as RINEX 3.04 `{stem}.obs`.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let (_, name) = stem_and_name(raw);
    let r = decode(&data, &name);
    if r.obs.t.is_empty() {
        return Err(DecodeError::Format(format!(
            "{}: no NovAtel range logs (RANGEB / RANGECMPB) found; log one of them in the receiver to \
             record data that can be processed",
            raw.display()
        )));
    }
    let list: Vec<String> = r
        .messages
        .iter()
        .map(|(m, n)| format!("{m} x{n}"))
        .collect();
    log(&format!(
        "NovAtel: {} epochs, {} satellites; messages {}",
        r.obs.t.len(),
        r.obs.sats.len(),
        list.join(", ")
    ));
    let info = serde_json::json!({ "receiver": "NovAtel", "epochs": r.obs.t.len() });
    write_obs(
        raw,
        outdir,
        &r.obs,
        r.glonass_fcn,
        "NovAtel",
        r.warnings,
        info,
    )
}

/// A binary message around `body` (for tests and for writing logs).
#[must_use]
pub fn frame(id: u16, week: u16, ms: u32, body: &[u8]) -> Vec<u8> {
    let mut h = vec![0u8; 28];
    h[..4].copy_from_slice(&[0xAA, 0x44, 0x12, 28]);
    h[4..6].copy_from_slice(&id.to_le_bytes());
    h[8..10].copy_from_slice(&(body.len() as u16).to_le_bytes());
    h[13] = 180; // fine steering
    h[14..16].copy_from_slice(&week.to_le_bytes());
    h[16..20].copy_from_slice(&ms.to_le_bytes());
    h.extend_from_slice(body);
    let crc = crc32(&h);
    h.extend(crc.to_le_bytes());
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tracking status: phase, code and parity lock, system, signal type.
    fn stat(sys: u32, sig: u32) -> u32 {
        (1 << 10) | (1 << 11) | (1 << 12) | (sys << 16) | (sig << 21)
    }

    fn rangeb(obs: &[(u16, u16, f64, f64, f32, u32)]) -> Vec<u8> {
        let mut b = (obs.len() as u32).to_le_bytes().to_vec();
        for &(prn, gfrq, psr, adr, lockt, st) in obs {
            let mut q = vec![0u8; 44];
            q[0..2].copy_from_slice(&prn.to_le_bytes());
            q[2..4].copy_from_slice(&gfrq.to_le_bytes());
            q[4..12].copy_from_slice(&psr.to_le_bytes());
            q[16..24].copy_from_slice(&adr.to_le_bytes());
            q[28..32].copy_from_slice(&(-800.0f32).to_le_bytes());
            q[32..36].copy_from_slice(&44.0f32.to_le_bytes());
            q[36..40].copy_from_slice(&lockt.to_le_bytes());
            q[40..44].copy_from_slice(&st.to_le_bytes());
            b.extend(q);
        }
        b
    }

    #[test]
    fn rangeb_reads_and_flags_a_slip() {
        let a = rangeb(&[
            (7, 0, 2.2e7, -1.15e8, 100.0, stat(0, 0)),
            (7, 0, 2.2e7 + 2.0, -9.0e7, 100.0, stat(0, 9)),
            (40, 9, 2.0e7, -1.07e8, 50.0, stat(1, 0)),
        ]);
        let b = rangeb(&[(7, 0, 2.2e7, -1.15e8, 0.5, stat(0, 0))]);
        let d = [
            b"Com1>".to_vec(),
            frame(43, 2100, 1000, &a),
            frame(43, 2100, 2000, &b),
        ]
        .concat();
        assert!(recognises(&d));
        let r = decode(&d, "X");
        assert_eq!(r.obs.t.len(), 2);
        let g7 = &r.obs.sats[&Sat::new(GnssSystem::Gps, 7)];
        let l1 = g7.iter().find(|(c, _)| c.to_string() == "L1C").unwrap().1;
        assert_eq!(l1.values[0], 1.15e8);
        assert_eq!(l1.lli.as_ref().unwrap()[1] & 1, 1);
        assert!(g7.keys().any(|c| c.to_string() == "C2W"));
        assert_eq!(r.glonass_fcn[&Sat::new(GnssSystem::Glonass, 3)], 2);
    }

    #[test]
    fn rangecmp_phase_rebuilt() {
        let psr: f64 = 21_000_000.0;
        let f1 = 1575.42e6;
        let cycles = psr * f1 / CLIGHT; // a phase that agrees with the range
        let raw = -(cycles - MAXVAL * (cycles / MAXVAL).round());
        let mut q = vec![0u8; 24];
        q[0..4].copy_from_slice(&stat(0, 0).to_le_bytes());
        let dop = (-1000.0f64 * 256.0) as i32 as u32 & 0x0FFF_FFFF;
        let psr_q = (psr * 128.0).round() as u64;
        let w = u32::from_le_bytes([q[4], q[5], q[6], q[7]]) | dop;
        q[4..8].copy_from_slice(&w.to_le_bytes());
        let lo = ((psr_q & 0x0FFF_FFFF) as u32) << 4;
        let w7 = (u32::from(q[7]) & 0x0F) | lo;
        q[7..11].copy_from_slice(&w7.to_le_bytes());
        q[11] = (psr_q >> 28) as u8;
        q[12..16].copy_from_slice(&((raw * 256.0).round() as i32).to_le_bytes());
        q[17] = 9;
        q[18] = 0x40; // lock time 2 s
        q[20..22].copy_from_slice(&(25u16 << 5).to_le_bytes()); // C/N0 45
        let mut body = 1u32.to_le_bytes().to_vec();
        body.extend(q);
        let r = decode(&frame(140, 2100, 5000, &body), "X");
        let g9 = &r.obs.sats[&Sat::new(GnssSystem::Gps, 9)];
        let get = |s: &str| {
            g9.iter()
                .find(|(c, _)| c.to_string() == s)
                .unwrap()
                .1
                .values[0]
        };
        assert!((get("C1C") - psr).abs() < 0.01);
        assert!(
            (get("L1C") - cycles).abs() < 0.01,
            "{} {}",
            get("L1C"),
            cycles
        );
        assert!((get("D1C") + 1000.0).abs() < 0.01);
        assert_eq!(get("S1C"), 45.0);
    }
}
