//! u-blox UBX logs: the raw measurement messages RXM-RAW (0x02 0x10, GPS L1 of older receivers) and
//! RXM-RAWX (0x02 0x15, every system and signal of M8 / F9 receivers) into the engine's observations.
//!
//! The frame (sync `B5 62`, class, id, little-endian length, payload, 8-bit Fletcher checksum over class to
//! payload end) and the message layouts follow u-blox's public protocol descriptions; the field use (which
//! tracking-status bits make a range or a phase valid, the signal-id table, the half-cycle and lock-time
//! rules for loss of lock, the half-cycle shift of BeiDou GEO phases) follows RTKLIB's `ublox.c`
//! (T. Takasu, BSD-2-Clause). No RTKLIB code is copied; the logic is restated here.
//!
//! Both messages carry the GPS week, so no outside date is needed. Navigation messages are not decoded.

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::consts::WEEK;
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};

use crate::epochs::{Epochs, stem_and_name, write_obs};
use crate::{Converted, DecodeError};

/// The UBX frames of a stream: (class, id, payload) with a good checksum, and the number of bad ones.
pub(crate) fn frames(data: &[u8]) -> (Vec<(u8, u8, &[u8])>, usize) {
    let mut out = Vec::new();
    let mut bad = 0;
    let mut i = 0;
    while i + 8 <= data.len() {
        if data[i] != 0xB5 || data[i + 1] != 0x62 {
            i += 1;
            continue;
        }
        let len = usize::from(u16::from_le_bytes([data[i + 4], data[i + 5]]));
        let end = i + 6 + len + 2;
        if end > data.len() {
            i += 1;
            continue;
        }
        let (mut a, mut b) = (0u8, 0u8);
        for &x in &data[i + 2..end - 2] {
            a = a.wrapping_add(x);
            b = b.wrapping_add(a);
        }
        if a == data[end - 2] && b == data[end - 1] {
            out.push((data[i + 2], data[i + 3], &data[i + 6..end - 2]));
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
fn r4(d: &[u8], o: usize) -> f64 {
    f64::from(f32::from_bits(u4(d, o)))
}
fn r8(d: &[u8], o: usize) -> f64 {
    let mut b = [0u8; 8];
    b.copy_from_slice(&d[o..o + 8]);
    f64::from_le_bytes(b)
}

/// The RINEX code tail of a RAWX signal (`sigId`, protocol 27.11 and later), or `None` when unknown.
fn rawx_tail(sys: GnssSystem, sig: u8, ver: u8) -> Option<&'static str> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Qzss, Sbas};
    if ver == 0 {
        // older firmware: one signal per system
        return Some(match sys {
            BeiDou => "2I",
            Galileo => "1X",
            _ => "1C",
        });
    }
    Some(match (sys, sig) {
        (Gps | Glonass | Galileo | Qzss | Sbas, 0) => "1C",
        (Gps, 3) | (Qzss, 5) => "2L",
        (Gps | Qzss, 4) => "2S",
        (Glonass, 2) => "2C",
        (Galileo, 1) => "1B",
        (Galileo, 5) => "7I",
        (Galileo, 6) => "7Q",
        (Qzss, 1) => "1Z",
        (BeiDou, 0 | 1) => "2I",
        (BeiDou, 2 | 3) => "7I",
        _ => return None,
    })
}

/// The satellite of a RAWX (`gnssId`, `svId`).
fn rawx_sat(gnss: u8, sv: u8) -> Option<Sat> {
    let sys = match gnss {
        0 => GnssSystem::Gps,
        1 => GnssSystem::Sbas,
        2 => GnssSystem::Galileo,
        3 => GnssSystem::BeiDou,
        5 => GnssSystem::Qzss,
        6 => GnssSystem::Glonass,
        _ => return None,
    };
    let prn = match sys {
        GnssSystem::Sbas if sv >= 120 => sv - 100,
        GnssSystem::Qzss if sv >= 193 => sv - 192,
        GnssSystem::Sbas => return None,
        _ => sv,
    };
    (prn > 0 && prn < 100).then(|| Sat::new(sys, prn))
}

/// What [`decode`] found.
#[derive(Debug)]
pub struct Ubx {
    /// The observations.
    pub obs: StationObs,
    /// GLONASS frequency channels from RAWX.
    pub glonass_fcn: BTreeMap<Sat, i32>,
    /// Messages seen, by (class, id).
    pub messages: BTreeMap<(u8, u8), usize>,
    /// Things the user must see.
    pub warnings: Vec<String>,
}

/// Decodes a UBX log.
#[must_use]
pub fn decode(data: &[u8], name: &str) -> Ubx {
    let (frames, bad) = frames(data);
    let mut messages: BTreeMap<(u8, u8), usize> = BTreeMap::new();
    let mut epochs = Epochs::default();
    let mut fcn: BTreeMap<Sat, i32> = BTreeMap::new();
    // lock time (ms) and half-cycle state per satellite and signal
    let mut lock: BTreeMap<(Sat, &'static str), (u16, bool)> = BTreeMap::new();
    let mut unknown = 0usize;
    for (class, id, p) in frames {
        *messages.entry((class, id)).or_insert(0) += 1;
        match (class, id) {
            (0x02, 0x10) if p.len() >= 8 => {
                let week = u2(p, 4);
                let n = usize::from(p[6]);
                if week == 0 || p.len() < 8 + 24 * n {
                    continue;
                }
                let t = f64::from(week) * WEEK + f64::from(u4(p, 0)) * 0.001;
                for k in 0..n {
                    let q = &p[8 + 24 * k..8 + 24 * (k + 1)];
                    let sv = q[20];
                    let sat = match sv {
                        1..=32 => Sat::new(GnssSystem::Gps, sv),
                        120..=158 => Sat::new(GnssSystem::Sbas, sv - 100),
                        _ => continue,
                    };
                    let c = epochs.cell(t, sat, "1C");
                    c.phase = r8(q, 0);
                    c.code = r8(q, 8);
                    c.doppler = r4(q, 16);
                    c.snr = f64::from(q[22] as i8);
                    c.lli = i8::from(q[23] & 1 != 0);
                }
            }
            (0x02, 0x15) if p.len() >= 16 => {
                let week = u2(p, 8);
                let n = usize::from(p[11]);
                let ver = p[13];
                if week == 0 || p.len() < 16 + 32 * n {
                    continue;
                }
                let t = f64::from(week) * WEEK + r8(p, 0);
                // the epoch on the receiver's millisecond, as RINEX keeps it
                let t = (t * 1000.0).round() * 0.001;
                for k in 0..n {
                    let q = &p[16 + 32 * k..16 + 32 * (k + 1)];
                    let Some(sat) = rawx_sat(q[20], q[21]) else {
                        continue;
                    };
                    let Some(tail) = rawx_tail(sat.sys, q[22], ver) else {
                        unknown += 1;
                        continue;
                    };
                    let trk = q[30];
                    let cpstd = q[28] & 15;
                    let lockt = u2(q, 24);
                    let f = if sat.sys == GnssSystem::Glonass {
                        let ch = i32::from(q[23]) - 7;
                        fcn.insert(sat, ch);
                        Some(ch)
                    } else {
                        None
                    };
                    let band = tail.as_bytes()[0] - b'0';
                    let fr = freq(sat, band, f).ok();
                    let pr = r8(q, 0);
                    let mut cp = r8(q, 8);
                    let halfv = trk & 4 != 0;
                    let halfc = trk & 8 != 0;
                    let prev = lock.insert((sat, tail), (lockt, halfc));
                    let slip = lockt == 0 || prev.is_some_and(|(l, h)| lockt < l || h != halfc);
                    let c = epochs.cell(t, sat, tail);
                    if trk & 1 != 0 {
                        c.code = pr;
                    }
                    if trk & 2 != 0 && cp != -0.5 && cpstd <= 5 {
                        if sat.sys == GnssSystem::BeiDou && (sat.prn <= 5 || sat.prn >= 59) {
                            cp += 0.5; // BeiDou GEO phases are half a cycle off
                        }
                        c.phase = cp;
                        c.lli = i8::from(slip) | (i8::from(!halfv) << 1);
                    }
                    if fr.is_some() {
                        c.doppler = r4(q, 16);
                    }
                    c.snr = f64::from(q[26]);
                }
            }
            _ => {}
        }
    }
    let mut warnings = Vec::new();
    if bad > 0 {
        warnings.push(format!(
            "{bad} UBX frames failed their checksum and were left out"
        ));
    }
    if unknown > 0 {
        warnings.push(format!(
            "{unknown} RXM-RAWX measurements of signals not known here were left out"
        ));
    }
    let obs = epochs.into_obs(name);
    fcn.retain(|s, _| obs.sats.contains_key(s));
    Ubx {
        obs,
        glonass_fcn: fcn,
        messages,
        warnings,
    }
}

/// Whether a stream reads as UBX: two good frames in its head.
#[must_use]
pub fn recognises(head: &[u8]) -> bool {
    frames(head).0.len() >= 2
}

/// Decodes a UBX log into `outdir` as RINEX 3.04 `{stem}.obs`.
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
            "{}: no u-blox raw measurement messages (RXM-RAW / RXM-RAWX) found; turn them on in the \
             receiver to log data that can be processed",
            raw.display()
        )));
    }
    let list: Vec<String> = r
        .messages
        .iter()
        .map(|((c, i), n)| format!("{c:02X}-{i:02X} x{n}"))
        .collect();
    log(&format!(
        "u-blox: {} epochs, {} satellites; messages {}",
        r.obs.t.len(),
        r.obs.sats.len(),
        list.join(", ")
    ));
    let info = serde_json::json!({ "receiver": "u-blox", "epochs": r.obs.t.len() });
    write_obs(
        raw,
        outdir,
        &r.obs,
        r.glonass_fcn,
        "u-blox UBX",
        r.warnings,
        info,
    )
}

/// A UBX frame around `payload` (for tests and for writing logs).
#[must_use]
pub fn frame(class: u8, id: u8, payload: &[u8]) -> Vec<u8> {
    let mut out = vec![0xB5, 0x62, class, id];
    out.extend((payload.len() as u16).to_le_bytes());
    out.extend_from_slice(payload);
    let (mut a, mut b) = (0u8, 0u8);
    for &x in &out[2..] {
        a = a.wrapping_add(x);
        b = b.wrapping_add(a);
    }
    out.extend([a, b]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use geomaticist_gnss_core::consts::CLIGHT as C;

    fn rawx(week: u16, tow: f64, meas: &[(u8, u8, u8, f64, f64, u16, u8)]) -> Vec<u8> {
        let mut p = vec![0u8; 16];
        p[0..8].copy_from_slice(&tow.to_le_bytes());
        p[8..10].copy_from_slice(&week.to_le_bytes());
        p[11] = meas.len() as u8;
        p[13] = 1;
        for &(g, sv, sig, pr, cp, lockt, trk) in meas {
            let mut q = vec![0u8; 32];
            q[0..8].copy_from_slice(&pr.to_le_bytes());
            q[8..16].copy_from_slice(&cp.to_le_bytes());
            q[16..20].copy_from_slice(&(-1234.5f32).to_le_bytes());
            q[20] = g;
            q[21] = sv;
            q[22] = sig;
            q[23] = 7 + 1;
            q[24..26].copy_from_slice(&lockt.to_le_bytes());
            q[26] = 45;
            q[28] = 2;
            q[30] = trk;
            p.extend(q);
        }
        frame(0x02, 0x15, &p)
    }

    #[test]
    fn rawx_signals_and_slips() {
        let mut d = b"junk".to_vec();
        d.extend(rawx(
            2300,
            345_600.000_4,
            &[
                (0, 5, 0, 2.1e7, 1.1e8, 5000, 0x0F),
                (0, 5, 3, 2.1e7 + 3.0, 8.6e7, 5000, 0x0F),
                (6, 3, 0, 1.9e7, 1.0e8, 500, 0x07),
                (3, 2, 0, 3.6e7, 1.9e8, 500, 0x0F),
                (2, 11, 9, 2.4e7, 1.2e8, 500, 0x0F),
            ],
        ));
        d.extend(rawx(
            2300,
            345_601.0,
            &[(0, 5, 0, 2.1e7 + C * 1e-9, 1.1e8 + 5.0, 100, 0x0F)],
        ));
        let r = decode(&d, "TEST");
        assert_eq!(
            r.obs.t,
            vec![2300.0 * WEEK + 345_600.0, 2300.0 * WEEK + 345_601.0]
        );
        let g5 = Sat::new(GnssSystem::Gps, 5);
        let codes: Vec<String> = r.obs.sats[&g5].keys().map(ToString::to_string).collect();
        assert!(codes.contains(&"L2L".to_owned()), "{codes:?}");
        let l1 = r.obs.sats[&g5]
            .iter()
            .find(|(c, _)| c.to_string() == "L1C")
            .unwrap()
            .1;
        assert_eq!(l1.values[0], 1.1e8);
        // the lock time went down: a slip
        assert_eq!(l1.lli.as_ref().unwrap()[1] & 1, 1);
        assert_eq!(r.glonass_fcn[&Sat::new(GnssSystem::Glonass, 3)], 1);
        // BeiDou GEO (C02) phase moved half a cycle
        let c2 = &r.obs.sats[&Sat::new(GnssSystem::BeiDou, 2)];
        let l = c2.iter().find(|(c, _)| c.to_string() == "L2I").unwrap().1;
        assert_eq!(l.values[0], 1.9e8 + 0.5);
        assert!(r.warnings.iter().any(|w| w.contains("1 RXM-RAWX")));
    }

    #[test]
    fn raw_gps_l1() {
        let mut p = vec![0u8; 8];
        p[0..4].copy_from_slice(&1_000_000u32.to_le_bytes());
        p[4..6].copy_from_slice(&2000u16.to_le_bytes());
        p[6] = 1;
        let mut q = vec![0u8; 24];
        q[0..8].copy_from_slice(&1.05e8f64.to_le_bytes());
        q[8..16].copy_from_slice(&2.0e7f64.to_le_bytes());
        q[16..20].copy_from_slice(&500.0f32.to_le_bytes());
        q[20] = 12;
        q[22] = 40;
        q[23] = 1;
        p.extend(q);
        let d = [frame(0x02, 0x10, &p), frame(0x02, 0x10, &p)].concat();
        assert!(recognises(&d));
        let r = decode(&d, "X");
        assert_eq!(r.obs.t, vec![2000.0 * WEEK + 1000.0]);
        let s = &r.obs.sats[&Sat::new(GnssSystem::Gps, 12)];
        assert_eq!(s.len(), 4);
    }
}
