//! Trimble .DAT / .TGD records, as Trimble's runpkr00 writes them from a .T01.
//!
//! A file is a chain of records: 't' (0x74), length low byte, sub-type (or the length's high byte for the long
//! types 0x11 and 0x1B), type. Measurements are type 0x1B (the .TGD from `runpkr00 -g`), plain little-endian
//! doubles.

use std::collections::HashMap;

use crate::util::{latin1, le};

const LONG_TYPES: [u8; 2] = [0x11, 0x1B];

/// (system, band, code) -> RINEX 2.11 (code obs, phase obs).
#[must_use]
pub fn signal(sysid: u8, band: u8, code: u8) -> Option<(&'static str, &'static str)> {
    match (sysid, band, code) {
        (0 | 2, 0, 0) => Some(("C1", "L1")),
        (0, 1, 2) | (2, 1, 1) => Some(("P2", "L2")),
        (0, 1, 5) => Some(("C2", "L2")),
        (0, 2, 8) => Some(("C5", "L5")),
        _ => None,
    }
}

/// TGD system ids seen.
#[must_use]
pub fn system(sysid: u8) -> char {
    match sysid {
        0 => 'G',
        2 => 'R',
        _ => '?',
    }
}

/// `(offset, type, sub, length)` of every record.
pub fn records(data: &[u8]) -> Result<Vec<(usize, u8, u8, usize)>, String> {
    let mut out = Vec::new();
    let (mut i, n) = (0usize, data.len());
    while i + 4 <= n {
        if data[i] != 0x74 {
            return Err(format!("no record start at byte {}", crate::t04::group3(i)));
        }
        let t = data[i + 3];
        let ln = if LONG_TYPES.contains(&t) {
            usize::from(data[i + 1]) | (usize::from(data[i + 2]) << 8)
        } else {
            usize::from(data[i + 1])
        };
        if ln < 4 || i + ln > n {
            return Err(format!(
                "bad record length at byte {}",
                crate::t04::group3(i)
            ));
        }
        out.push((i, t, data[i + 2], ln));
        i += ln;
    }
    Ok(out)
}

/// One signal.
#[derive(Debug, Clone, PartialEq)]
pub struct Signal {
    /// Band.
    pub band: u8,
    /// Code.
    pub code: u8,
    /// dB-Hz.
    pub snr: f64,
    /// Range, m.
    pub range: f64,
    /// Phase, cycles, RINEX sign.
    pub phase: f64,
    /// 0x20 phase valid, 0x04 tracked since the last epoch, 0x08 new arc.
    pub flags: u8,
}

/// One satellite.
#[derive(Debug, Clone, PartialEq)]
pub struct Sat {
    /// System letter.
    pub sys: char,
    /// System id.
    pub sysid: u8,
    /// PRN / slot.
    pub prn: u8,
    /// GLONASS frequency channel (signed); 0 for GPS.
    pub channel: i64,
    /// Elevation.
    pub elev: u8,
    /// Azimuth.
    pub azim: i64,
    /// Signals.
    pub signals: Vec<Signal>,
}

impl Sat {
    /// "G05".
    #[must_use]
    pub fn id(&self) -> String {
        format!("{}{:02}", self.sys, self.prn)
    }
}

/// One epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct Epoch {
    /// GPS week.
    pub week: i64,
    /// Time of week, s.
    pub tow: f64,
    /// Satellites.
    pub sats: Vec<Sat>,
}

/// Header fields from the runpkr00 outputs.
#[derive(Debug, Clone, Default)]
pub struct Header {
    /// Serial.
    pub serial: String,
    /// Model ("R8 Model 2").
    pub model: String,
    /// Model code ("R8-2").
    pub model_code: String,
    /// Firmware.
    pub firmware: String,
    /// Marker.
    pub marker: String,
    /// Antenna height (as the point record's text).
    pub arp_height: Option<f64>,
    /// Antenna code.
    pub antenna_code: Option<i64>,
    /// The last position fix (deg, deg, m).
    pub approx_llh: Option<(f64, f64, f64)>,
    /// Interval, s.
    pub interval_s: Option<f64>,
    /// Text records (first value of each key).
    pub text: HashMap<String, String>,
}

fn byte(r: &[u8], i: usize) -> Result<u8, String> {
    r.get(i)
        .copied()
        .ok_or_else(|| "index out of range".to_owned())
}

/// One type-0x1B record (without its 4-byte record header).
pub fn parse_epoch(r: &[u8]) -> Result<Epoch, String> {
    let hl = usize::from(byte(r, 0)?);
    if r.len() < 11 {
        return Err("unpack requires a buffer of 11 bytes".into());
    }
    let week = le::u16(r, 1);
    let tow = le::f64(r, 3) / 1000.0; // ms of the GPS week
    let mut ep = Epoch {
        week,
        tow,
        sats: Vec::new(),
    };
    let mut p = hl;
    while p < r.len() {
        let sl = usize::from(byte(r, p)?);
        let (prn, sysid, fl, nsig, el) = (
            byte(r, p + 1)?,
            byte(r, p + 2)?,
            byte(r, p + 3)?,
            byte(r, p + 4)?,
            byte(r, p + 5)?,
        );
        byte(r, p + 8)?;
        let az = le::i16(r, p + 7);
        if sl == 0 {
            return Err("satellite block of length 0".into());
        }
        p += sl;
        let mut sat = Sat {
            sys: system(sysid),
            sysid,
            prn,
            channel: if sysid == 2 { i64::from(fl as i8) } else { 0 },
            elev: el,
            azim: az,
            signals: Vec::new(),
        };
        for _ in 0..nsig {
            let gl = usize::from(byte(r, p)?);
            byte(r, p + 20)?;
            sat.signals.push(Signal {
                band: byte(r, p + 1)?,
                code: byte(r, p + 2)?,
                snr: le::u16(r, p + 3) as f64 / 10.0,
                range: le::f64(r, p + 5),
                phase: -le::f64(r, p + 13),
                flags: if gl > 26 { byte(r, p + 26)? } else { 0 },
            });
            p += gl;
        }
        ep.sats.push(sat);
    }
    Ok(ep)
}

/// Text with trailing NUL and space removed (`decode("latin-1").rstrip("\x00 ")`).
fn text(r: &[u8]) -> String {
    latin1(r).trim_end_matches(['\0', ' ']).to_owned()
}

fn find_digits10(r: &[u8]) -> Option<String> {
    r.windows(10)
        .find(|w| w.iter().all(u8::is_ascii_digit))
        .map(latin1)
}

/// `re.search(rb"\x00\x0b(\d{10})\x00\x05([\x20-\x7e]{1,4})\x00", r)`.
fn find_model(r: &[u8]) -> Option<(String, String)> {
    for i in 0..r.len() {
        let rest = &r[i..];
        if rest.len() < 2 + 10 + 2 + 1 + 1 || rest[0] != 0 || rest[1] != 0x0b {
            continue;
        }
        let serial = &rest[2..12];
        if !serial.iter().all(u8::is_ascii_digit) || rest[12] != 0 || rest[13] != 5 {
            continue;
        }
        // greedy [\x20-\x7e]{1,4} followed by \x00, with backtracking
        let printable = rest[14..]
            .iter()
            .take(4)
            .take_while(|c| (0x20..=0x7e).contains(*c))
            .count();
        for k in (1..=printable).rev() {
            if rest.get(14 + k) == Some(&0) {
                return Some((latin1(serial), latin1(&rest[14..14 + k])));
            }
        }
    }
    None
}

/// `re.match(r"\s*([-+]?\d+(\.\d+)?)", s)` -> the number.
fn leading_number(s: &str) -> Option<f64> {
    let t = s.trim_start_matches(crate::util::py_space);
    let b = t.as_bytes();
    let mut i = 0;
    if matches!(b.first(), Some(b'-' | b'+')) {
        i += 1;
    }
    let d0 = i;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == d0 {
        return None;
    }
    if i + 1 < b.len() && b[i] == b'.' && b[i + 1].is_ascii_digit() {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
    }
    t[..i].parse().ok()
}

/// Header fields from the runpkr00 outputs (.DAT and / or .TGD).
#[allow(clippy::cast_precision_loss)]
pub fn read_header(blobs: &[&[u8]]) -> Result<Header, String> {
    let mut h = Header::default();
    for data in blobs {
        let mut settings: Vec<u8> = Vec::new();
        for (i, t, s, ln) in records(data)? {
            let r = &data[i..i + ln];
            if t == 0x05 && s == 0x02 {
                if let Some(m) = r.get(16..).and_then(find_digits10)
                    && h.serial.is_empty()
                {
                    h.serial = m;
                }
                if h.firmware.is_empty() {
                    h.firmware = format!("{:.2}", le::u16(r, 6) as f64 / 100.0);
                }
            } else if t == 0x1E && s == 0x01 {
                if let Some((serial, code)) = find_model(r) {
                    if h.serial.is_empty() {
                        h.serial = serial;
                    }
                    h.model_code = code;
                }
            } else if t == 0x10 && s == 0xC7 {
                let txt = text(&r[4..]);
                let (k, v) = txt.split_once(':').unwrap_or((&txt, ""));
                h.text.entry(k.to_owned()).or_insert_with(|| v.to_owned());
            } else if t == 0x10 && s == 0x03 {
                // point record: name [0:16], height text [48:64]
                let b = &r[4..];
                h.marker = text(&b[..b.len().min(16)]);
                if let Some(v) = leading_number(&latin1(&b[48.min(b.len())..64.min(b.len())])) {
                    h.arp_height = Some(v);
                }
            } else if t == 0x0B {
                // position fix: lat, lon (semicircles), height
                if r.len() < 28 {
                    return Err("unpack_from requires a buffer of at least 28 bytes".into());
                }
                h.approx_llh = Some((
                    le::f64(r, 4) * 180.0,
                    le::f64(r, 12) * 180.0,
                    le::f64(r, 20),
                ));
            } else if t == 0x1F {
                settings.extend_from_slice(&r[4..]);
            }
        }
        if !settings.is_empty() && h.antenna_code.is_none() {
            // re.search(rb"\x03[\x00\x01]\x1b(..)", settings, re.S)
            h.antenna_code = settings
                .windows(5)
                .find(|w| w[0] == 3 && w[1] <= 1 && w[2] == 0x1b)
                .map(|w| i64::from(u16::from_le_bytes([w[3], w[4]])));
        }
    }
    h.model = match h.model_code.as_str() {
        "R8-2" => "R8 Model 2".to_owned(), // the receiver's model code (record 0x1E/1) -> REC TYPE
        other => other.to_owned(),
    };
    if h.antenna_code.is_none() && h.model_code == "R8-2" {
        h.antenna_code = Some(186); // integrated antenna by model, when the settings don't name one (assumed)
    }
    let iv = h
        .text
        .get("SessionMeasIntervalMsecs")
        .map_or("", String::as_str);
    let iv = iv.trim_start_matches(crate::util::py_space);
    let digits: String = iv.chars().take_while(char::is_ascii_digit).collect();
    if !digits.is_empty() {
        h.interval_s = digits.parse::<i64>().ok().map(|v| v as f64 / 1000.0);
    }
    Ok(h)
}

/// The measurement epochs of a .TGD.
pub fn epochs(tgd: &[u8]) -> Result<Vec<Epoch>, String> {
    let mut out = Vec::new();
    for (i, t, _s, ln) in records(tgd)? {
        if t == 0x1B {
            out.push(parse_epoch(&tgd[i + 4..i + ln])?);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regexes() {
        let r = b"xx\x00\x0b1234567890\x00\x05R8-2\x00yy";
        assert_eq!(find_model(r), Some(("1234567890".into(), "R8-2".into())));
        assert_eq!(leading_number("  1.5000xyz"), Some(1.5));
        assert_eq!(leading_number("-2."), Some(-2.0));
        assert_eq!(leading_number("abc"), None);
    }

    #[test]
    fn record_chain() {
        let d = [0x74, 4, 0, 0x05, 0x74, 5, 0, 0x1B, 9];
        let r = records(&d).unwrap();
        assert_eq!(r, vec![(0, 0x05, 0, 4), (4, 0x1B, 0, 5)]);
        assert!(records(&[0x75, 4, 0, 0]).is_err());
    }
}
