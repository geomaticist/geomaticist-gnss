//! Header facts from a T04 inner stream: everything the RINEX header needs
//!.
//!
//! - type 0x10 sub 0xC7: ASCII "Key:Value" records (ReceiverId, RefStationName, AntennaArpHeight, ...).
//! - type 0x10 sub 0x03: fixed-width point record: name [0:16], number [67:83] after the 3-byte header
//!   -> RINEX MARKER NAME / MARKER NUMBER. The last one in the file has the measured antenna height as text at
//!   [48:64].
//! - type 0x1e sub 0x01: serial, model, firmware as big-endian u16 x100 (0x026d = 621 -> "6.21").
//! - type 0x16 sub 0x07: antenna setup (215 bytes): BE doubles lat, lon (semicircles), height (m) of the
//!   reference position, then the measured antenna height (m) at +27; measure-up method index at +53; start /
//!   end minute of day at +69 / +71. One is written when logging starts and one when it ends.
//! - rec35 sub 2: first position fix (lat raw/2^40 deg, lon raw/2^39 deg, hgt raw/2^12 m) -> APPROX POSITION.

use std::collections::{BTreeMap, HashMap};

use super::T04Error;
use super::avro;
use super::container::{Inner, Records};
use super::rec35::{self, KIND_FULL, Messages};
use crate::antenna_lib;
use crate::util::{be, latin1, py_float, py_int, py_strip};

/// One antenna setup record (type 0x16 sub 7).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Setup {
    /// Reference latitude, degrees.
    pub lat: f64,
    /// Reference longitude, degrees.
    pub lon: f64,
    /// Reference height, m.
    pub h: f64,
    /// Measured antenna height, m.
    pub height: f64,
    /// Measure-up method index.
    pub method: u8,
    /// Start minute of day.
    pub start_min: i64,
    /// End minute of day.
    pub end_min: i64,
}

/// Header facts.
#[derive(Debug, Clone, Default)]
pub struct HeaderInfo {
    /// Text records, first value of each key.
    pub text: HashMap<String, String>,
    /// Receiver serial.
    pub serial: String,
    /// Receiver model.
    pub model: String,
    /// Firmware ("6.21").
    pub firmware: String,
    /// Marker name.
    pub marker_name: String,
    /// Marker number.
    pub marker_number: String,
    /// Antenna code.
    pub antenna_code: Option<i64>,
    /// IGS antenna name (or `UNKNOWN<code>`).
    pub antenna_type: String,
    /// Radome.
    pub antenna_dome: String,
    /// Antenna height to the ARP, m.
    pub antenna_height_arp: f64,
    /// Logging interval, s.
    pub interval_s: Option<f64>,
    /// First position fix (deg, deg, m).
    pub approx_llh: Option<(f64, f64, f64)>,
    /// Antenna setup records, file order.
    pub setups: Vec<Setup>,
    /// AntennaArpHeight as the text record has it.
    pub header_arp_height: Option<f64>,
    /// Measured height of the setup in force (m).
    pub setup_height: Option<f64>,
    /// Measured height in the last point record's text.
    pub point_height: Option<f64>,
    /// Its measure-up method index.
    pub setup_method: Option<u8>,
    /// Why the height isn't the header's, if it isn't.
    pub height_note: String,
    /// GLONASS slot -> frequency number k, from the ephemerides (first of each).
    pub glo_nav_k: BTreeMap<i64, i64>,
}

impl HeaderInfo {
    /// APPROX POSITION XYZ (zeros without a fix).
    #[must_use]
    pub fn approx_xyz(&self) -> (f64, f64, f64) {
        self.approx_llh.map_or((0.0, 0.0, 0.0), |(la, lo, h)| {
            crate::util::llh_deg_to_xyz(la, lo, h)
        })
    }
}

const SETUP_SIG: &[u8] = b"\xd7\x16\x07"; // length 215, type 0x16, sub 7
const POINT_SIG: &[u8] = b"\x56\x10\x03"; // length 86, type 0x10, sub 3: the point record

/// The measured antenna height a point record (0x10 sub 3, 86 bytes) holds as text at [48:64] after the 3-byte
/// header: a plain number like "1.5000000". Another form ("1.8000V123") or none -> `None`.
#[must_use]
pub fn point_record_height(b: &[u8]) -> Option<f64> {
    if b.len() != 86 || !b.starts_with(POINT_SIG) {
        return None;
    }
    let txt = latin1(&b[3 + 48..3 + 64]);
    let t = py_strip(&txt);
    let (a, c) = t.split_once('.')?;
    if a.is_empty() || c.is_empty() || !a.bytes().chain(c.bytes()).all(|x| x.is_ascii_digit()) {
        return None;
    }
    t.parse().ok()
}

/// Measured heights of the point records found by their signature in a byte range.
#[must_use]
pub fn scan_point_heights(buf: &[u8]) -> Vec<f64> {
    scan(buf, POINT_SIG, 86, point_record_height)
}

fn scan<T>(buf: &[u8], sig: &[u8], len: usize, f: impl Fn(&[u8]) -> Option<T>) -> Vec<T> {
    let mut out = Vec::new();
    let mut i = 0;
    while i + sig.len() <= buf.len() {
        if &buf[i..i + sig.len()] == sig
            && let Some(v) = f(&buf[i..(i + len).min(buf.len())])
        {
            out.push(v);
        }
        i += 1;
    }
    out
}

/// One antenna setup record (215 bytes) -> [`Setup`], or `None` if it doesn't look like one.
#[must_use]
pub fn parse_setup(b: &[u8]) -> Option<Setup> {
    if b.len() < 215 || !b.starts_with(SETUP_SIG) {
        return None;
    }
    let (lat, lon, hgt, ant) = (
        be::f64(b, 3),
        be::f64(b, 11),
        be::f64(b, 19),
        be::f64(b, 27),
    );
    if !(ant.is_finite()
        && (0.0..10.0).contains(&ant)
        && lat.abs() <= 0.5
        && lon.abs() <= 1.0
        && b[53] < 16)
    {
        return None;
    }
    Some(Setup {
        lat: lat * 180.0,
        lon: lon * 180.0,
        h: hgt,
        height: ant,
        method: b[53],
        start_min: be::u16(b, 69),
        end_min: be::u16(b, 71),
    })
}

/// Antenna setup records found by their signature in a byte range that needn't start at a record.
#[must_use]
pub fn scan_setups(buf: &[u8]) -> Vec<Setup> {
    scan(buf, SETUP_SIG, 215, parse_setup)
}

/// Which antenna height this log really carries.
///
/// The AntennaArpHeight text record is what the receiver held when logging started; a setup entered afterwards
/// is in the setup record written when logging ends. The last point record's text has the measured
/// height too, and it wins when the setup record agrees with it.
pub fn reconcile_height(h: &mut HeaderInfo) {
    h.header_arp_height = Some(h.antenna_height_arp);
    let mut point_note = String::new();
    if let (Some(ph), Some(last), Some(code)) =
        (h.point_height, h.setups.last().copied(), h.antenna_code)
    {
        let mm = antenna_lib::measure_methods(code);
        let known = mm.is_some_and(|(_, ms)| usize::from(last.method) < ms.len());
        let red = if known {
            antenna_lib::reduce_height(code, usize::from(last.method), ph)
        } else {
            None
        };
        let same = (last.height - ph).abs() < 1e-5;
        let at_apc = red.is_some_and(|(_, apc)| (last.height - apc).abs() < 1e-4);
        if known && !(same || at_apc) {
            point_note = format!(
                "; the point record's {ph:.4} m doesn't match the setup record's {:.4} m — not used, check the field book",
                last.height
            );
        }
        if let (true, Some((arp, _))) = (known && (same || at_apc), red) {
            h.setup_height = Some(ph);
            h.setup_method = Some(last.method);
            if (arp - h.antenna_height_arp).abs() < 0.0005 {
                return; // the header already has this setup's height
            }
            let name = mm.map_or("", |(_, ms)| ms[usize::from(last.method)].0);
            h.height_note = format!(
                "the header's antenna height {:.4} m (ARP) is from an earlier setup; the log's point record has \
                 {ph:.4} m to {name} = {arp:.4} m to the ARP — used{}",
                h.antenna_height_arp,
                if same {
                    String::new()
                } else {
                    format!(
                        " (its setup record had {:.4} m, the height at the phase center)",
                        last.height
                    )
                }
            );
            h.antenna_height_arp = geomaticist_gnss_core::text::round_to(arp, 6);
            return;
        }
    }
    setup_record_rule(h);
    if !point_note.is_empty() {
        if h.height_note.is_empty() {
            h.height_note = format!(
                "the header's antenna height {:.4} m (ARP) is kept{point_note}",
                h.antenna_height_arp
            );
        } else {
            h.height_note.push_str(&point_note);
        }
    }
}

/// The setup record written at the end of logging, its height read in its measure-up method's units.
fn setup_record_rule(h: &mut HeaderInfo) {
    if h.setups.len() < 2 {
        return;
    }
    let Some(code) = h.antenna_code else { return };
    let (first, last) = (h.setups[0], h.setups[h.setups.len() - 1]);
    if (last.height - first.height).abs() < 1e-6 && (last.h - first.h).abs() < 1e-3 {
        return; // no setup arrived during logging: the header stands
    }
    let Some((_, ms)) = antenna_lib::measure_methods(code) else {
        return;
    };
    if usize::from(last.method) >= ms.len() {
        return;
    }
    let Some((arp, _)) = antenna_lib::reduce_height(code, usize::from(last.method), last.height)
    else {
        return;
    };
    h.setup_height = Some(last.height);
    h.setup_method = Some(last.method);
    if (arp - h.antenna_height_arp).abs() < 0.0005 {
        return;
    }
    let name = ms[usize::from(last.method)].0;
    h.height_note = format!(
        "the header's antenna height {:.4} m (ARP) is from an earlier setup; the setup record written at the end \
         of logging has {:.4} m to {name} = {arp:.4} m to the ARP — used",
        h.antenna_height_arp, last.height
    );
    h.antenna_height_arp = geomaticist_gnss_core::text::round_to(arp, 6);
}

fn firmware_of(b: &[u8]) -> Option<String> {
    // ... serial \0 [len] model \0 [len] fw u16-BE ...   (len bytes observed as 0x06)
    let nul: Vec<usize> = b
        .iter()
        .enumerate()
        .filter(|(_, c)| **c == 0)
        .map(|(i, _)| i)
        .collect();
    for w in nul.windows(2) {
        let (i, j) = (w[0], w[1]);
        if j - i > 2 {
            let s = &b[i + 2..j];
            if s.is_ascii()
                && [&b"R"[..], b"SPS", b"NetR", b"Alloy"]
                    .iter()
                    .any(|p| s.starts_with(p))
            {
                let fw = b.get(j + 2..(j + 4).min(b.len())).unwrap_or(&[]);
                let v = fw.iter().fold(0i64, |a, &x| (a << 8) | i64::from(x));
                #[allow(clippy::cast_precision_loss)]
                return Some(format!("{:.2}", v as f64 / 100.0));
            }
        }
    }
    None
}

/// Reads the header facts.
pub fn read_header(inner: &Inner) -> Result<HeaderInfo, T04Error> {
    let mut h = HeaderInfo {
        antenna_dome: "NONE".into(),
        ..HeaderInfo::default()
    };
    for r in Records::new(inner) {
        let r = r?;
        let b = r.raw;
        if r.rtype == 0x16 && b.len() == 215 && b[2] == 0x07 {
            if let Some(s) = parse_setup(b) {
                h.setups.push(s);
            }
            continue;
        }
        if r.rtype == 0x10 && b.len() > 3 && b[2] == 0xC7 {
            let txt = latin1(&b[3..]);
            let txt = txt.trim_end_matches('\0');
            let (key, val) = txt.split_once(':').unwrap_or((txt, ""));
            h.text
                .entry(key.to_owned())
                .or_insert_with(|| val.to_owned());
        } else if r.rtype == 0x10 && b.len() >= 83 && b[2] == 0x03 {
            if let Some(v) = point_record_height(b) {
                h.point_height = Some(v); // the last point record's counts
            }
            if h.marker_name.is_empty() {
                let txt: Vec<char> = latin1(&b[3..]).chars().collect();
                let sl = |a: usize, z: usize| -> String {
                    txt.get(a.min(txt.len())..z.min(txt.len()))
                        .unwrap_or(&[])
                        .iter()
                        .collect()
                };
                h.marker_name = py_strip(&sl(0, 16)).to_owned();
                h.marker_number = py_strip(&sl(67, 83)).to_owned(); // 16 chars at +67
            }
        } else if r.rtype == 0x1C && b.len() == 137 && b[2] == 0x02 {
            // GLONASS ephemeris
            let k = i64::from(b[21] as i8);
            h.glo_nav_k.entry(i64::from(b[3])).or_insert(k);
        } else if r.rtype == 0x1E
            && b.len() > 3
            && b[2] == 0x01
            && h.firmware.is_empty()
            && let Some(fw) = firmware_of(b)
        {
            h.firmware = fw;
        }
    }
    let rid: Vec<&str> = h
        .text
        .get("ReceiverId")
        .map_or("", String::as_str)
        .split(',')
        .collect();
    if rid.len() >= 3 {
        h.model = py_strip(rid[1]).to_owned();
        h.serial = py_strip(rid[2]).to_owned();
    }
    h.antenna_code = py_int(h.text.get("AntennaType").map_or("", String::as_str));
    if let Some((igs, dome, _, _)) = h.antenna_code.and_then(antenna_lib::antenna) {
        h.antenna_type = igs.to_owned();
        h.antenna_dome = dome.to_owned();
    } else {
        h.antenna_type = format!(
            "UNKNOWN{}",
            h.antenna_code
                .map_or_else(|| "None".to_owned(), |c| c.to_string())
        );
        h.antenna_dome = "NONE".into();
    }
    if let Some(v) = py_float(h.text.get("AntennaArpHeight").map_or("0", String::as_str)) {
        h.antenna_height_arp = v;
    }
    if let Some(ms) = h
        .text
        .get("SessionMeasIntervalMsecs")
        .and_then(|s| py_int(s))
    {
        #[allow(clippy::cast_precision_loss)]
        {
            h.interval_s = Some(ms as f64 / 1000.0);
        }
    }
    h.approx_llh = first_fix(inner)?;
    reconcile_height(&mut h);
    Ok(h)
}

#[allow(clippy::cast_precision_loss)]
fn first_fix(inner: &Inner) -> Result<Option<(f64, f64, f64)>, T04Error> {
    let schemas = rec35::schemas(inner)?;
    let Some(sch) = schemas.get(&2) else {
        return Ok(None);
    };
    for m in Messages::new(inner, |s| s == 2) {
        let m = m?;
        if m.subtype == 2 && m.kind == KIND_FULL {
            let r = avro::decode_full(sch, &m.body).map_err(rec35::Rec35Error::from)?;
            if let Some(p) = r
                .get("pos")
                .and_then(avro::Value::as_array)
                .and_then(<[_]>::first)
            {
                let g = |n: &str| p.get(n).and_then(avro::Value::as_long).unwrap_or(0) as f64;
                return Ok(Some((
                    g("latPos") / 2f64.powi(40),
                    g("lonPos") / 2f64.powi(39),
                    g("hgtPos") / 2f64.powi(12),
                )));
            }
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_height_forms() {
        let mut b = vec![0u8; 86];
        b[..3].copy_from_slice(POINT_SIG);
        b[51..67].copy_from_slice(b"1.5000000       ");
        assert_eq!(point_record_height(&b), Some(1.5));
        b[51..67].copy_from_slice(b"1.8000V123      ");
        assert_eq!(point_record_height(&b), None);
    }

    #[test]
    fn firmware_field() {
        let b = b"\x40\x1e\x01\x00\x061234567890\x00\x06R10-2\x00\x06\x02\x6d";
        assert_eq!(firmware_of(b).as_deref(), Some("6.21"));
    }
}
