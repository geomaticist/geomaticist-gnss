//! The built-in reader of packed Trimble logs (.T01): the file to Trimble's unpacked record chain, with no
//! outside program.
//!
//! Three layers: a bit stream that codes bytes with a fixed prefix code ([`code`]); packets in that byte stream
//! ([`packet`]); and their contents: measurement epochs in full or as differences ([`epoch`]), navigation
//! data, receiver information and positions ([`other`]). The output is a chain of unpacked records in the
//! layout of a .TGD file, so the RINEX writers read it as they read runpkr00's output ([`super::runpkr`]).
//!
//! A clean-room decoder written from sample receiver logs (R8 Model 2 firmware 4.61, R8 GNSS firmware 4.64;
//! GPS and GLONASS). A file of a layout it doesn't know goes to runpkr00 when that is installed.

pub(crate) mod code;
pub(crate) mod epoch;
pub(crate) mod other;
pub(crate) mod packet;

use crate::DecodeError;

use epoch::{Doppler, Meas};
use other::Position;
use packet::{T_EPOCH_DELTA, T_EPOCH_FULL, T_POS_DELTA, T_POS_FULL};

/// What the built-in reader made of a packed log.
#[derive(Debug, Clone, Default)]
pub struct Unpacked {
    /// The unpacked records, in the layout of runpkr00's .TGD (measurements as type 0x1B), with a short position
    /// record (type 0x0B) per position fix.
    pub records: Vec<u8>,
    /// Measurement epochs decoded.
    pub epochs: usize,
    /// Measurement epochs in the file that could not be decoded: damage, the epochs after it up to the next
    /// full epoch, or a layout this reader doesn't know.
    pub lost_epochs: usize,
    /// Bytes of the decoded stream that belonged to no packet (damage).
    pub skipped_bytes: usize,
    /// What went wrong, for the log (a few lines at most).
    pub notes: Vec<String>,
}

impl Unpacked {
    /// Whether everything in the file was read.
    #[must_use]
    pub fn complete(&self) -> bool {
        self.lost_epochs == 0 && self.skipped_bytes == 0 && self.notes.is_empty()
    }

    fn note(&mut self, s: String) {
        if self.notes.len() < 8 {
            self.notes.push(s);
        }
    }
}

/// Whether `data` (the file, or its first kilobytes) starts like a packed log this reader knows: its first
/// packets pass their checks.
#[must_use]
pub fn recognises(data: &[u8]) -> bool {
    let Some(coded) = data.get(code::HEAD_LEN..) else {
        return false;
    };
    let stream = code::decode(&coded[..coded.len().min(8192)]);
    let (p, _) = packet::packets(&stream);
    p.iter().filter(|x| x.ty < 0x20).count() >= 2
}

/// Unpacks a packed log. An error means the file is not one this reader knows (the caller can then try
/// runpkr00); damage inside a known file is reported in the result instead.
pub fn unpack(data: &[u8]) -> Result<Unpacked, DecodeError> {
    let unknown = || DecodeError::Format("not a packed log the built-in reader knows".to_owned());
    let stream = code::decode(data.get(code::HEAD_LEN..).ok_or_else(unknown)?);
    let (packets, skipped) = packet::packets(&stream);
    if packets.is_empty() {
        return Err(unknown());
    }
    let mut u = Unpacked {
        skipped_bytes: skipped,
        ..Unpacked::default()
    };
    if skipped > 0 {
        u.note(format!("{skipped} bytes of the log are damaged"));
    }
    let mut meas: Option<Meas> = None;
    let mut dt = 0u16;
    let mut doppler = Doppler::default();
    let mut pos: Option<Position> = None;
    for p in &packets {
        if p.after_gap {
            // a position packet may be missing; the positions after it can't be checked, so they wait for the
            // next full one. (A missing epoch packet shows in the check of the next epoch.)
            pos = None;
        }
        match p.ty {
            T_EPOCH_FULL | T_EPOCH_DELTA => {
                let res = if p.ty == T_EPOCH_FULL {
                    Meas::full(p.payload, dt).map(|m| {
                        meas = Some(m);
                    })
                } else {
                    match meas.as_mut() {
                        Some(m) => m.apply_delta(p.payload, p.check),
                        None => Err("no full epoch to build on".to_owned()),
                    }
                };
                match (res, meas.as_ref()) {
                    (Ok(()), Some(m)) => {
                        dt = m.dt;
                        let d = doppler.next(m);
                        u.records.extend(epoch::record(m, &d));
                        u.epochs += 1;
                    }
                    (Err(e), _) => {
                        if meas.is_some() || p.ty == T_EPOCH_FULL {
                            u.note(format!(
                                "epoch packet {}: {e}",
                                u.epochs + u.lost_epochs + 1
                            ));
                        }
                        meas = None;
                        u.lost_epochs += 1;
                    }
                    (Ok(()), None) => {}
                }
            }
            T_POS_FULL => match Position::full(p.payload) {
                Ok(x) => {
                    u.records.extend(x.record());
                    pos = Some(x);
                }
                Err(e) => {
                    pos = None;
                    u.note(e);
                }
            },
            T_POS_DELTA => {
                if let Some(x) = pos.as_mut() {
                    match x.apply_delta(p.payload) {
                        Ok(()) => u.records.extend(x.record()),
                        Err(e) => {
                            pos = None;
                            u.note(e);
                        }
                    }
                }
            }
            0x0f | 0x15 | 0x1c | 0x1f => {
                let r = match p.ty {
                    0x0f => other::ion_utc(p.payload),
                    0x15 => other::gps_ephemeris(p.payload),
                    0x1c => other::glonass_ephemeris(p.payload),
                    _ => other::settings(p.payload),
                };
                match r {
                    Some(r) => u.records.extend(r),
                    None => u.note(format!(
                        "a packet {:02x} of {} bytes (another layout)",
                        p.ty,
                        p.payload.len()
                    )),
                }
            }
            0x10 => {
                // text and point records: sub-type, then the body as it is
                if let Some((&sub, body)) = p.payload.split_first()
                    && body.len() <= 251
                {
                    u.records.extend(other::rec(0x10, sub, body));
                }
            }
            0x1e if p.payload.first() == Some(&1) => match other::receiver(p.payload) {
                Some(rs) => rs.iter().for_each(|r| u.records.extend(r)),
                None => u.note("the receiver information has another layout".to_owned()),
            },
            _ => {}
        }
    }
    if u.epochs == 0 {
        return Err(DecodeError::Format(match u.notes.first() {
            Some(n) => format!("no measurement epoch could be read ({n})"),
            None => "the log has no measurement epochs".to_owned(),
        }));
    }
    Ok(u)
}

/// A small packed log made for tests: receiver information, text and point records, navigation data, a
/// position, and measurement epochs (the first in full, the rest as differences) of two satellites.
#[cfg(test)]
pub(crate) mod synth {
    use super::epoch::{Meas, Sat, Sig, encode_delta};
    use super::{code, other, packet};

    /// The pieces of a synthetic log.
    pub(crate) struct Log {
        /// The decoded byte stream (0x03, then packets).
        pub stream: Vec<u8>,
        /// The measurement epochs it holds.
        pub epochs: Vec<Meas>,
    }

    fn sig(j: usize, band: u8, code: u8, snr: u16, rng: i64, ph: i64) -> Sig {
        Sig {
            len: if j == 0 { 0x11 } else { 0x0f },
            band,
            code,
            snr,
            rng,
            ph,
            slip: 3,
            fl: 0x2b,
            d_rng: 0,
            d_ph: 0,
        }
    }

    /// The first epoch: GPS 6 (L1 C/A, L2 P(Y)) and GLONASS 9 on channel -2 (L1, L2).
    pub(crate) fn first_epoch() -> Meas {
        let sat = |prn, sys, ch, sigs| Sat {
            prn,
            sys,
            ch,
            elev: 36,
            az: 30,
            b7: 2,
            sigs,
        };
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

    /// The epoch 15 s after `m`: the satellites move on, the receiver clock drifts.
    pub(crate) fn next_epoch(m: &Meas, k: i64) -> Meas {
        let mut n = m.clone();
        n.tow += 15_000;
        n.dt = 15_000;
        n.clk -= 1232 + k as i32;
        for (i, s) in n.sats.iter_mut().enumerate() {
            let rate = [-8000i64, 2500][i % 2] - 100 * k; // cycles per epoch on L1
            for (j, x) in s.sigs.iter_mut().enumerate() {
                x.fl = 0x23;
                x.snr = (i64::from(x.snr) + (k % 3) - 1) as u16;
                if j == 0 {
                    x.rng -= rate * 128 * 19 / 100;
                    x.ph += rate * 32768;
                } else {
                    x.rng += (k % 5) - 2;
                    x.ph += rate * 32768 * 60 / 77;
                }
            }
        }
        n
    }

    fn text(s: &str) -> Vec<u8> {
        let mut p = vec![0xc7];
        p.extend_from_slice(s.as_bytes());
        packet::classic(0x10, &p)
    }

    /// A log of `n` epochs.
    pub(crate) fn log(n: usize) -> Log {
        let mut stream = vec![0x03];
        // receiver information, page 1
        let mut p = vec![1, 4, 6, 0x11, 0x22, 0x33, 0x44, 0x00, 0x2c, 11];
        p.extend_from_slice(b"1234567890\0\x05R8-2\0");
        p.extend_from_slice(&[6, 0x01, 0xcd, 0x01, 0xcd, 0x01, 0x9b]);
        stream.extend(packet::classic(0x1e, &p));
        stream.extend(text("SessionMeasIntervalMsecs:15000"));
        stream.extend(text("AntennaType:186"));
        // settings with the antenna code 186, a point record (name, antenna height text at 48)
        stream.extend(packet::classic(
            0x1f,
            &[0x0a, 0x4b, 0, 3, 1, 0x1b, 0xba, 0x00, 0],
        ));
        let mut point = vec![0x03];
        point.extend_from_slice(format!("{:<48}{:<16}", "TEST", "1.5000").as_bytes());
        stream.extend(packet::classic(0x10, &point));
        // ionosphere / UTC, one GPS and one GLONASS ephemeris (field values only need to survive the trip)
        let mut ion = Vec::new();
        for k in 0..11 {
            ion.extend_from_slice(&(f64::from(k) * 1e-8).to_be_bytes());
        }
        ion.extend_from_slice(&[0, 17, 0, 17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        stream.extend(packet::classic(0x0f, &ion));
        let mut eph = vec![6];
        eph.extend_from_slice(&2000u16.to_be_bytes());
        eph.extend_from_slice(&[0, 79, 0, 79]);
        for t in [99_000u32, 100_800, 100_800] {
            eph.extend_from_slice(&t.to_be_bytes());
        }
        for k in 0..19 {
            eph.extend_from_slice(&(f64::from(k) * 1e-3).to_be_bytes());
        }
        eph.extend_from_slice(&[0; 14]);
        stream.extend(packet::classic(0x15, &eph));
        let mut glo = vec![2, 9];
        glo.extend_from_slice(&2000u16.to_be_bytes());
        glo.extend_from_slice(&100_000u32.to_be_bytes());
        glo.extend_from_slice(&[
            1, 0x23, 0x45, 0x01, 0x02, 0, 0, 0x12, 0x34, 0, 0, 0xfe, 0, 1, 3,
        ]);
        for k in 0..14 {
            glo.extend_from_slice(&(f64::from(k) * 1e3).to_be_bytes());
        }
        stream.extend(packet::classic(0x1c, &glo));
        // a position
        let pos = other::encode_position(45.25, -110.5, 123.5);
        let check = packet::record_check(packet::sum16(&pos), pos.len() + 6, 0x1D);
        stream.extend(packet::long(packet::T_POS_FULL, &pos, Some(check)));
        // the epochs
        let mut epochs = vec![first_epoch()];
        stream.extend(packet::long(
            packet::T_EPOCH_FULL,
            &epochs[0].full_payload(),
            None,
        ));
        for k in 1..n {
            let prev = &epochs[k - 1];
            let mut next = next_epoch(prev, k as i64);
            let (payload, check) = encode_delta(prev, &mut next);
            stream.extend(packet::long(packet::T_EPOCH_DELTA, &payload, Some(check)));
            epochs.push(next);
        }
        Log { stream, epochs }
    }

    /// The file of a byte stream: a 13-byte head, then the coded stream with a flush in the middle.
    pub(crate) fn file(stream: &[u8]) -> Vec<u8> {
        let mut f = vec![0xa5; code::HEAD_LEN];
        f.extend(code::encode(stream, &[stream.len() / 2]));
        f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::t01::records;

    #[test]
    fn a_synthetic_log_unpacks_to_its_records() {
        let log = synth::log(6);
        let u = unpack(&synth::file(&log.stream)).unwrap();
        assert!(u.complete(), "{:?}", u.notes);
        assert_eq!((u.epochs, u.lost_epochs), (6, 0));
        // the record chain parses, in the order of the file
        let recs = records::records(&u.records).unwrap();
        let types: Vec<u8> = recs.iter().map(|r| r.1).collect();
        assert_eq!(
            types,
            [
                0x05, 0x1e, 0x10, 0x10, 0x1f, 0x10, 0x0f, 0x15, 0x1c, 0x0b, 0x1b, 0x1b, 0x1b, 0x1b,
                0x1b, 0x1b
            ]
        );
        // the epochs are the ones the log was made from
        let eps = records::epochs(&u.records).unwrap();
        assert_eq!(eps.len(), 6);
        for (e, m) in eps.iter().zip(&log.epochs) {
            assert_eq!(
                (e.week, e.tow),
                (i64::from(m.week), f64::from(m.tow) / 1000.0)
            );
            assert_eq!(e.sats.len(), m.sats.len());
            for (s, t) in e.sats.iter().zip(&m.sats) {
                assert_eq!(
                    (s.prn, s.sysid, s.signals.len()),
                    (t.prn, t.sys, t.sigs.len())
                );
                assert_eq!(s.signals[0].range, t.sigs[0].rng as f64 / 128.0);
                assert_eq!(s.signals[1].phase, -(t.sigs[1].ph as f64) / 32768.0);
            }
        }
        assert_eq!(eps[0].sats[1].channel, -2);
        // flag 0x04 (a Doppler value from the epoch before) from the second epoch on
        assert_eq!(eps[0].sats[0].signals[0].flags, 0x2b);
        assert_eq!(eps[1].sats[0].signals[0].flags, 0x27);
        // the header the writers read
        let h = records::read_header(&[&u.records]).unwrap();
        assert_eq!(
            (h.serial.as_str(), h.model.as_str()),
            ("1234567890", "R8 Model 2")
        );
        assert_eq!((h.firmware.as_str(), h.marker.as_str()), ("4.61", "TEST"));
        assert_eq!(
            (h.arp_height, h.antenna_code, h.interval_s),
            (Some(1.5), Some(186), Some(15.0))
        );
        let (lat, lon, height) = h.approx_llh.unwrap();
        assert!((lat - 45.25).abs() < 1e-12 && (lon + 110.5).abs() < 1e-12);
        assert_eq!(height, 123.5);
    }

    #[test]
    fn files_that_are_not_packed_logs_are_refused() {
        assert!(unpack(&[]).is_err());
        assert!(!recognises(&[]) && !recognises(&[0; 4000]));
        let log = synth::log(2);
        let file = synth::file(&log.stream);
        assert!(recognises(&file) && recognises(&file[..200]));
        assert!(unpack(&[0; 12]).is_err());
        assert!(unpack(&[0; 4000]).is_err());
        assert!(
            unpack(
                b"     3.05           OBSERVATION DATA    M                   RINEX VERSION / TYPE"
            )
            .is_err()
        );
        // a log without epochs
        let log = synth::log(1);
        let cut = log.stream.len() - log.epochs[0].full_payload().len() - 6;
        assert!(unpack(&synth::file(&log.stream[..cut])).is_err());
    }

    #[test]
    fn a_cut_log_gives_the_epochs_before_the_cut() {
        let log = synth::log(8);
        let file = synth::file(&log.stream);
        let whole = unpack(&file).unwrap();
        let u = unpack(&file[..file.len() - 40]).unwrap();
        assert!(u.epochs >= 5 && u.epochs < 8, "{}", u.epochs);
        assert_eq!(u.lost_epochs, 0);
        assert!(whole.records.starts_with(&u.records));
    }

    #[test]
    fn damage_loses_the_epochs_up_to_the_next_full_one() {
        // two sessions of five epochs, each starting with a full epoch
        let a = synth::log(5);
        let mut stream = a.stream.clone();
        let m = synth::next_epoch(&a.epochs[4], 9);
        stream.extend(packet::long(T_EPOCH_FULL, &m.full_payload(), None));
        let good = unpack(&synth::file(&stream)).unwrap();
        assert_eq!((good.epochs, good.lost_epochs), (6, 0));
        // damage inside the second packet of differences: it and the two after it are lost
        let (packets, _) = packet::packets(&stream);
        let deltas: Vec<usize> = packets
            .iter()
            .filter(|p| p.ty == T_EPOCH_DELTA)
            .map(|p| p.payload.as_ptr() as usize - stream.as_ptr() as usize)
            .collect();
        let mut bad = stream.clone();
        bad[deltas[1] + 5] ^= 0x10;
        let u = unpack(&synth::file(&bad)).unwrap();
        assert_eq!((u.epochs, u.lost_epochs), (3, 3));
        assert!(!u.complete() && !u.notes.is_empty());
        let eps = records::epochs(&u.records).unwrap();
        let tows: Vec<f64> = eps.iter().map(|e| e.tow).collect();
        assert_eq!(tows, [100_000.0, 100_015.0, 100_075.0]);
        // a damaged text packet between epochs: the epoch after it still builds on the one before
        let mut stream = vec![0x03];
        stream.extend(packet::long(
            T_EPOCH_FULL,
            &a.epochs[0].full_payload(),
            None,
        ));
        let mut text = packet::classic(0x10, b"\xc7Note:1");
        text[6] ^= 1;
        stream.extend(text);
        let mut next = synth::next_epoch(&a.epochs[0], 1);
        let (payload, check) = epoch::encode_delta(&a.epochs[0], &mut next);
        stream.extend(packet::long(T_EPOCH_DELTA, &payload, Some(check)));
        let u = unpack(&synth::file(&stream)).unwrap();
        assert_eq!((u.epochs, u.lost_epochs, u.skipped_bytes), (2, 0, 13));
        // the same with the first byte of the stream damaged
        stream[0] = 0x55;
        let u = unpack(&synth::file(&stream)).unwrap();
        assert_eq!((u.epochs, u.lost_epochs, u.skipped_bytes), (2, 0, 14));
    }
}
