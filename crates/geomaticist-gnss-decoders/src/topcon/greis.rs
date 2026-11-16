//! Standard GREIS measurement messages (JAVAD / Topcon `.jps` logs, and `.tps` logs without `rD`): `~~`
//! receiver time, `RD` date, `SI` / `NN` satellite list, `R?` / `r?` / `?R` / `?r` pseudoranges, `P?` / `p?` /
//! `?P` / `?p` carrier phases, `D?` / `?d` Doppler, `E?` / `?E` C/N0 and `TC` lock time.
//!
//! The message layouts and scale factors follow the GREIS reference guide; the order of decoding and the
//! loss-of-lock rule follow RTKLIB `javad.c` (BSD-2-Clause, T. Takasu), restated here, not copied.

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::time::{DateTime, gps_seconds, utc_to_gps_seconds};
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};

use super::frame::Messages;
use crate::epochs::{Epochs, near_from_mtime, stem_and_name, write_obs};
use crate::{Converted, DecodeError};

const DAY: f64 = 86_400.0;

/// The decoded log.
#[derive(Debug)]
pub struct Greis {
    /// Observations.
    pub obs: StationObs,
    /// GLONASS channel numbers seen.
    pub glonass_fcn: BTreeMap<Sat, i32>,
    /// Message counts by ID.
    pub messages: BTreeMap<String, usize>,
    /// Things worth telling the user.
    pub warnings: Vec<String>,
}

/// A satellite of the `SI` list: the satellite (None where not known yet) and its GLONASS channel.
#[derive(Debug, Clone, Copy)]
struct Slot {
    sat: Option<Sat>,
    sys: GnssSystem,
    fcn: Option<i32>,
}

fn slot_of(usi: u8) -> Option<Slot> {
    let (sys, prn, fcn) = match usi {
        1..=37 => (GnssSystem::Gps, usi, None),
        38..=70 => (GnssSystem::Glonass, 0, Some(i32::from(usi) - 45)),
        71..=119 => (GnssSystem::Galileo, usi - 70, None),
        120..=142 => (GnssSystem::Sbas, usi - 100, None),
        193..=197 => (GnssSystem::Qzss, usi - 192, None),
        211..=240 => (GnssSystem::BeiDou, usi - 210, None),
        _ => return None,
    };
    Some(Slot {
        sat: (prn > 0).then(|| Sat::new(sys, prn)),
        sys,
        fcn,
    })
}

/// The RINEX code tail of a GREIS signal character on a system.
fn tail(sys: GnssSystem, sig: u8) -> Option<&'static str> {
    use GnssSystem::{BeiDou, Galileo, Glonass, Gps, Navic, Qzss, Sbas};
    let i = match sig {
        b'c' | b'C' => 0,
        b'1' => 1,
        b'2' => 2,
        b'3' => 3,
        b'5' => 4,
        b'l' => 5,
        _ => return None,
    };
    let row: [&str; 6] = match sys {
        Gps => ["1C", "1W", "2W", "2X", "5X", "1X"],
        Qzss => ["1C", "1Z", "6X", "2X", "5X", "1X"],
        Sbas => ["1C", "", "", "", "5X", ""],
        Galileo => ["1X", "8X", "7X", "6X", "5X", ""],
        Glonass => ["1C", "1P", "2P", "2C", "3X", ""],
        BeiDou => ["2I", "", "7I", "6I", "5X", "1X"],
        Navic => return None,
    };
    let t = row[i];
    (!t.is_empty()).then_some(t)
}

fn band(tail: &str) -> u8 {
    tail.as_bytes()[0] - b'0'
}

fn rd_i(b: &[u8], i: usize, n: usize) -> Option<&[u8]> {
    b.get(i * n..i * n + n)
}

struct Decoder {
    near: f64,
    epochs: Epochs,
    glonass_fcn: BTreeMap<Sat, i32>,
    /// GPS seconds of 00:00 of the receiver-time day, from `RD` (None: settle from `near`).
    day0: Option<f64>,
    /// Current epoch (GPS s).
    t: Option<f64>,
    tod_prev: Option<f64>,
    slots: Vec<Slot>,
    pr_ca: Vec<f64>,
    dp_ca: Vec<f64>,
    lock: BTreeMap<Sat, u16>,
    messages: BTreeMap<String, usize>,
    bad_len: usize,
}

impl Decoder {
    fn on_time(&mut self, b: &[u8]) {
        let Some(ms) = b.get(..4) else { return };
        let tod = f64::from(u32::from_le_bytes([ms[0], ms[1], ms[2], ms[3]])) * 1e-3;
        let day0 = if let Some(d) = self.day0 {
            // the day rolls over before the next RD
            if self.tod_prev.is_some_and(|p| tod < p - DAY / 2.0) {
                self.day0 = Some(d + DAY);
                d + DAY
            } else {
                d
            }
        } else {
            let mut d = (self.near / DAY).floor() * DAY;
            let t = d + tod;
            if t > self.near + DAY / 2.0 {
                d -= DAY;
            } else if t < self.near - DAY / 2.0 {
                d += DAY;
            }
            d
        };
        self.tod_prev = Some(tod);
        self.t = Some(day0 + tod);
        self.pr_ca = vec![0.0; self.slots.len()];
        self.dp_ca = vec![0.0; self.slots.len()];
    }

    fn on_date(&mut self, b: &[u8]) {
        if b.len() < 5 {
            return;
        }
        let y = i64::from(u16::from_le_bytes([b[0], b[1]]));
        let Some(dt) = DateTime::new(y, i64::from(b[2]), i64::from(b[3]), 0, 0, 0) else {
            return;
        };
        // base 0 GPS, 1 UTC (USNO), 2 GLONASS (UTC + 3 h), 3 UTC (SU)
        let d = match b[4] {
            0 => gps_seconds(&dt),
            2 => utc_to_gps_seconds(&dt) - 3.0 * 3600.0,
            _ => utc_to_gps_seconds(&dt),
        };
        self.day0 = Some(d);
        self.tod_prev = None;
    }

    fn on_sats(&mut self, b: &[u8]) {
        self.slots = b
            .iter()
            .map(|&u| {
                slot_of(u).unwrap_or(Slot {
                    sat: None,
                    sys: GnssSystem::Gps,
                    fcn: None,
                })
            })
            .collect();
        self.pr_ca = vec![0.0; self.slots.len()];
        self.dp_ca = vec![0.0; self.slots.len()];
    }

    fn on_slots(&mut self, b: &[u8]) {
        let mut k = 0;
        for s in &mut self.slots {
            if s.sys != GnssSystem::Glonass {
                continue;
            }
            let Some(&n) = b.get(k) else { break };
            k += 1;
            if (1..=32).contains(&n) {
                let sat = Sat::new(GnssSystem::Glonass, n);
                s.sat = Some(sat);
                if let Some(f) = s.fcn {
                    self.glonass_fcn.insert(sat, f);
                }
            }
        }
    }

    /// Per satellite of the list: (index, satellite, fcn, tail, its bytes) when the body has `size` bytes each.
    fn rows<'b>(
        &mut self,
        b: &'b [u8],
        size: usize,
        sig: u8,
    ) -> Vec<(usize, Sat, Option<i32>, &'static str, &'b [u8])> {
        if self.t.is_none() || self.slots.is_empty() {
            return Vec::new();
        }
        if b.len() != self.slots.len() * size {
            self.bad_len += 1;
            return Vec::new();
        }
        let mut out = Vec::new();
        for (i, s) in self.slots.iter().enumerate() {
            let (Some(sat), Some(t)) = (s.sat, tail(s.sys, sig)) else {
                continue;
            };
            if let Some(v) = rd_i(b, i, size) {
                out.push((i, sat, s.fcn, t, v));
            }
        }
        out
    }

    fn put(&mut self, sat: Sat, tail: &'static str, f: impl FnOnce(&mut crate::epochs::Cell)) {
        if let Some(t) = self.t {
            f(self.epochs.cell(t, sat, tail));
        }
    }

    fn message(&mut self, id: [u8; 2], b: &[u8]) {
        let (a, c) = (id[0], id[1]);
        match &id {
            b"~~" => return self.on_time(b),
            b"RD" => return self.on_date(b),
            b"SI" => return self.on_sats(b),
            b"NN" => return self.on_slots(b),
            b"TC" => return self.on_lock(b),
            _ if !is_obs_id(id) => return,
            _ => {}
        }
        match (a, c) {
            (b'R', s) => {
                for (i, sat, _, t, v) in self.rows(b, 8, s) {
                    let x = f64::from_le_bytes(v.try_into().unwrap_or([0; 8]));
                    if x == 0.0 {
                        continue;
                    }
                    let pr = x * CLIGHT;
                    if s == b'C' {
                        self.pr_ca[i] = pr;
                    }
                    self.put(sat, t, |cell| cell.code = pr);
                }
            }
            (b'r', s) => {
                for (i, sat, _, t, v) in self.rows(b, 4, s) {
                    let x = i32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    if x == 0x7FFF_FFFF {
                        continue;
                    }
                    let (k, add) = match sat.sys {
                        GnssSystem::Sbas => (1e-11, 0.125),
                        GnssSystem::Qzss => (2e-11, 0.125),
                        GnssSystem::BeiDou => (2e-11, 0.105),
                        GnssSystem::Galileo => (2e-11, 0.085),
                        _ => (1e-11, 0.075),
                    };
                    let pr = (f64::from(x) * k + add) * CLIGHT;
                    if s == b'c' {
                        self.pr_ca[i] = pr;
                    }
                    self.put(sat, t, |cell| cell.code = pr);
                }
            }
            (s, b'R') => {
                for (i, sat, _, t, v) in self.rows(b, 4, s) {
                    let x = f32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    if x == 0.0 || self.pr_ca[i] == 0.0 {
                        continue;
                    }
                    let pr = f64::from(x) * CLIGHT + self.pr_ca[i];
                    self.put(sat, t, |cell| cell.code = pr);
                }
            }
            (s, b'r') => {
                for (i, sat, _, t, v) in self.rows(b, 2, s) {
                    let x = i16::from_le_bytes(v.try_into().unwrap_or([0; 2]));
                    if x == 0x7FFF || self.pr_ca[i] == 0.0 {
                        continue;
                    }
                    let pr = (f64::from(x) * 1e-11 + 2e-7) * CLIGHT + self.pr_ca[i];
                    self.put(sat, t, |cell| cell.code = pr);
                }
            }
            (b'P', s) => {
                for (_, sat, _, t, v) in self.rows(b, 8, s) {
                    let x = f64::from_le_bytes(v.try_into().unwrap_or([0; 8]));
                    if x != 0.0 {
                        self.put(sat, t, |cell| cell.phase = x);
                    }
                }
            }
            (b'p', s) => {
                for (_, sat, _, t, v) in self.rows(b, 4, s) {
                    let x = u32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    if x != 0xFFFF_FFFF {
                        self.put(sat, t, |cell| cell.phase = f64::from(x) / 1024.0);
                    }
                }
            }
            (s, b'P') => {
                for (i, sat, fcn, t, v) in self.rows(b, 4, s) {
                    let x = f32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    let (Ok(f), true) = (freq(sat, band(t), fcn), self.pr_ca[i] != 0.0) else {
                        continue;
                    };
                    if x != 0.0 {
                        let cp = (f64::from(x) + self.pr_ca[i] / CLIGHT) * f;
                        self.put(sat, t, |cell| cell.phase = cp);
                    }
                }
            }
            (s, b'p') => {
                for (i, sat, fcn, t, v) in self.rows(b, 4, s) {
                    let x = i32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    let (Ok(f), true) = (freq(sat, band(t), fcn), self.pr_ca[i] != 0.0) else {
                        continue;
                    };
                    if x != 0x7FFF_FFFF {
                        let cp = (f64::from(x) * 2f64.powi(-40) + self.pr_ca[i] / CLIGHT) * f;
                        self.put(sat, t, |cell| cell.phase = cp);
                    }
                }
            }
            (b'D', s) => {
                for (i, sat, _, t, v) in self.rows(b, 4, s) {
                    let x = i32::from_le_bytes(v.try_into().unwrap_or([0; 4]));
                    if x == 0x7FFF_FFFF {
                        continue;
                    }
                    let dop = -f64::from(x) * 1e-4;
                    if s == b'C' {
                        self.dp_ca[i] = dop;
                    }
                    self.put(sat, t, |cell| cell.doppler = dop);
                }
            }
            (s, b'd') => {
                for (i, sat, fcn, t, v) in self.rows(b, 2, s) {
                    let x = i16::from_le_bytes(v.try_into().unwrap_or([0; 2]));
                    if x == 0x7FFF || self.dp_ca[i] == 0.0 {
                        continue;
                    }
                    let (Ok(f1), Ok(fnn)) = (freq(sat, 1, fcn), freq(sat, band(t), fcn)) else {
                        continue;
                    };
                    let dop = (-f64::from(x) + self.dp_ca[i] * 1e4) * fnn / f1 * 1e-4;
                    self.put(sat, t, |cell| cell.doppler = dop);
                }
            }
            (b'E', s) => {
                for (_, sat, _, t, v) in self.rows(b, 1, s) {
                    if v[0] != 255 {
                        self.put(sat, t, |cell| cell.snr = f64::from(v[0]));
                    }
                }
            }
            (s, b'E') => {
                for (_, sat, _, t, v) in self.rows(b, 1, s) {
                    if v[0] != 255 {
                        self.put(sat, t, |cell| cell.snr = f64::from(v[0]) * 0.25);
                    }
                }
            }
            _ => {}
        }
    }

    /// `TC`: C/A lock time (s) per satellite; it restarting or running backwards is a loss of lock.
    fn on_lock(&mut self, b: &[u8]) {
        for (_, sat, _, t, v) in self.rows(b, 2, b'c') {
            let tt = u16::from_le_bytes([v[0], v[1]]);
            if tt == 0xFFFF {
                continue;
            }
            let prev = self.lock.insert(sat, tt);
            if tt == 0 || prev.is_some_and(|p| tt < p) {
                self.put(sat, t, |cell| cell.lli |= 1);
            }
        }
    }
}

/// Whether the ID is a measurement or satellite-list message this decoder reads.
fn is_obs_id(id: [u8; 2]) -> bool {
    let meas = |c: u8| matches!(c, b'c' | b'C' | b'1' | b'2' | b'3' | b'5' | b'l');
    match id {
        [b'R' | b'r' | b'P' | b'p' | b'D' | b'E', s] => meas(s),
        [s, b'R' | b'r' | b'P' | b'p' | b'd' | b'E'] => meas(s),
        _ => false,
    }
}

/// Decodes the standard GREIS messages of `data`; `near` is an approximate GPS time (s) for logs without `RD`.
#[must_use]
pub fn decode(data: &[u8], name: &str, near: f64) -> Greis {
    let mut d = Decoder {
        near,
        epochs: Epochs::default(),
        glonass_fcn: BTreeMap::new(),
        day0: None,
        t: None,
        tod_prev: None,
        slots: Vec::new(),
        pr_ca: Vec::new(),
        dp_ca: Vec::new(),
        lock: BTreeMap::new(),
        messages: BTreeMap::new(),
        bad_len: 0,
    };
    let mut saw_date = false;
    for m in Messages::new(data) {
        *d.messages.entry(m.id_str().to_owned()).or_default() += 1;
        saw_date |= &m.id == b"RD";
        d.message(m.id, m.body);
    }
    let mut warnings = Vec::new();
    if !saw_date {
        warnings.push(
            "the log has no receiver date (RD); the day was taken from the file's date".to_owned(),
        );
    }
    if d.bad_len > 0 {
        warnings.push(format!(
            "{} measurement messages did not match the satellite list and were left out",
            d.bad_len
        ));
    }
    Greis {
        obs: d.epochs.into_obs(name),
        glonass_fcn: d.glonass_fcn,
        messages: d.messages,
        warnings,
    }
}

/// Whether the start of a file frames as GREIS with a satellite list and standard measurement messages.
#[must_use]
pub fn recognises(head: &[u8]) -> bool {
    let mut si = false;
    let mut meas = false;
    for m in Messages::new(head).take(4000) {
        si |= &m.id == b"SI";
        meas |= is_obs_id(m.id);
        if si && meas {
            return true;
        }
    }
    false
}

/// Decodes `raw` into `outdir` as RINEX 3.04 `{stem}.obs`; `near` is an approximate GPS time (s), else the
/// file's modification time.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    near: Option<f64>,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let (_, name) = stem_and_name(raw);
    let near = near.unwrap_or_else(|| near_from_mtime(raw));
    let r = decode(&data, &name, near);
    if r.obs.t.is_empty() {
        return Err(DecodeError::Format(format!(
            "{}: no GREIS measurement messages (SI with pseudoranges or phases) found",
            raw.display()
        )));
    }
    log(&format!(
        "GREIS: {} epochs, {} satellites, {} message types",
        r.obs.t.len(),
        r.obs.sats.len(),
        r.messages.len()
    ));
    let info = serde_json::json!({ "receiver": "JAVAD / Topcon GREIS", "epochs": r.obs.t.len() });
    write_obs(
        raw,
        outdir,
        &r.obs,
        r.glonass_fcn,
        "GREIS",
        r.warnings,
        info,
    )
}

/// A binary GREIS message: ID, 3-hex-digit length, body and its cs8 check.
#[must_use]
pub fn message(id: &[u8; 2], body: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(body.len() + 6);
    out.extend_from_slice(id);
    out.extend_from_slice(format!("{:03X}", body.len() + 1).as_bytes());
    out.extend_from_slice(body);
    out.push(super::frame::cs8(&out));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_an_epoch_and_a_slip() {
        let mut log = Vec::new();
        // 2011-01-15, GPS time base
        log.extend(message(b"RD", &[0xDB, 0x07, 1, 15, 0]));
        log.extend(message(b"SI", &[5, 46])); // G05 and a GLONASS satellite on channel 1
        log.extend(message(b"NN", &[7])); // R07
        for (k, tt) in [(0u32, 10u16), (1, 0)] {
            let tod = 3_600_000 + k * 1000;
            log.extend(message(b"~~", &tod.to_le_bytes()));
            let mut rc = Vec::new();
            rc.extend(1_000_000_000i32.to_le_bytes()); // (1e9 * 1e-11 + 0.075) s
            rc.extend(1_000_000_000i32.to_le_bytes());
            log.extend(message(b"rc", &rc));
            let mut cp = Vec::new();
            cp.extend(0i32.to_le_bytes());
            cp.extend(0x7FFF_FFFFi32.to_le_bytes());
            log.extend(message(b"cp", &cp));
            log.extend(message(b"CE", &[180, 255]));
            let mut tc = Vec::new();
            tc.extend(tt.to_le_bytes());
            tc.extend(0xFFFFu16.to_le_bytes());
            log.extend(message(b"TC", &tc));
        }
        assert!(recognises(&log));
        let g = decode(&log, "TEST", 0.0);
        let dt = DateTime::new(2011, 1, 15, 1, 0, 0).unwrap();
        assert_eq!(g.obs.t, vec![gps_seconds(&dt), gps_seconds(&dt) + 1.0]);
        let g5 = Sat::new(GnssSystem::Gps, 5);
        let r7 = Sat::new(GnssSystem::Glonass, 7);
        assert_eq!(g.glonass_fcn[&r7], 1);
        let pr = (1e9 * 1e-11 + 0.075) * CLIGHT;
        let c1 = g.obs.sats[&g5]
            .iter()
            .find(|(c, _)| c.to_string() == "C1C")
            .unwrap()
            .1;
        assert!((c1.values[0] - pr).abs() < 1e-6);
        let l1 = g.obs.sats[&g5]
            .iter()
            .find(|(c, _)| c.to_string() == "L1C")
            .unwrap()
            .1;
        let f = freq(g5, 1, None).unwrap();
        assert!((l1.values[0] - pr / CLIGHT * f).abs() < 1e-3);
        assert_eq!(l1.lli.as_ref().unwrap(), &vec![0, 1]);
        let s1 = g.obs.sats[&g5]
            .iter()
            .find(|(c, _)| c.to_string() == "S1C")
            .unwrap()
            .1;
        assert_eq!(s1.values[0], 45.0);
        // R07 has a code but no phase
        assert!(
            g.obs.sats[&r7]
                .keys()
                .all(|c| !c.to_string().starts_with('L'))
        );
    }
}
