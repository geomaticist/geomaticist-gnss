//! Decoder for Leica MDB raw logs (.m00, GS14 / GS15).
//!
//! A clean-room decoder written from sample receiver logs. The framing: `[9C AE][u16 n] ... [u16 check]`.

pub mod epochs;
pub mod motion;
pub mod nav;
pub mod obs_writer;

use std::path::Path;

use serde_json::{Map, Value, json};

use crate::util::{gps_week_fsec, le};
use crate::{Converted, DecodeError};

/// PGM of the files this decoder writes.
pub const PROGRAM: &str = "geomaticist-gnss MDB";
/// RUN BY.
pub const RUN_BY: &str = "geomaticist-gnss";

const SYNC: &[u8] = b"\x9c\xae";
/// Header record type.
pub const T_HEADER: u8 = 0x65;
/// Receiver record type.
pub const T_RECEIVER: u8 = 0x68;
/// Antenna record type.
pub const T_ANTENNA: u8 = 0x6c;
/// Position record type.
pub const T_POSITION: u8 = 0x6e;
/// Point record type.
pub const T_POINT: u8 = 0x71;

/// True for an extension the Leica maker owns (`.m00`).
#[must_use]
pub fn is_mdb_ext(ext_lower: &str) -> bool {
    ext_lower == "m00"
}

/// `(offset, type, body)` of each record, stopping at the first break in the framing.
/// The body runs from the type byte to before the check.
#[must_use]
pub fn records(data: &[u8]) -> Vec<(usize, u8, &[u8])> {
    let mut out = Vec::new();
    let mut pos = 0;
    while pos + 4 <= data.len() && &data[pos..pos + 2] == SYNC {
        let n = usize::from(u16::from_le_bytes([data[pos + 2], data[pos + 3]]));
        if n < 4 || pos + 2 + n > data.len() {
            break;
        }
        out.push((pos, data[pos + 4], &data[pos + 4..pos + n]));
        pos += 2 + n;
    }
    out
}

/// Every observation record of the log, decoded; a record that doesn't decode is skipped and its offset and
/// reason go to `errors`.
pub fn iter_epochs(data: &[u8], errors: &mut Vec<(usize, String)>) -> Vec<epochs::Epoch> {
    let mut out = Vec::new();
    for (pos, typ, body) in records(data) {
        if typ != epochs::T_EPOCH {
            continue;
        }
        match epochs::decode(body) {
            Ok(e) => out.push(e),
            Err(e) => errors.push((pos, e)),
        }
    }
    out
}

/// The log holds observation records in its first MiB (a Captivate job's own .m00 holds only points).
#[must_use]
pub fn has_epochs(data: &[u8]) -> bool {
    let head = &data[..data.len().min(1 << 20)];
    records(head).iter().any(|r| r.1 == epochs::T_EPOCH)
}

/// `b.decode("utf-16le", "ignore").split("\0")[0]`.
#[must_use]
pub fn utf16(b: &[u8]) -> String {
    let units = b
        .as_chunks::<2>()
        .0
        .iter()
        .map(|c| u16::from_le_bytes([c[0], c[1]]));
    let s: String = char::decode_utf16(units).filter_map(Result::ok).collect();
    s.split('\0').next().unwrap_or("").to_owned()
}

/// What the RINEX header needs, from the log's first records:
/// serial, model, firmware (receiver 0x68), antenna and its height-hook offset and phase-center offsets
/// (0x6c), the first point's name and antenna height (0x71), the first position (0x6e).
/// `arp_height` = the point's height + the hook offset (mark -> ARP, as written in ANTENNA: DELTA H).
#[must_use]
pub fn header_info(data: &[u8]) -> Map<String, Value> {
    let mut out = Map::new();
    let head = &data[..data.len().min(1 << 20)];
    for (_pos, typ, body) in records(head) {
        if typ == T_RECEIVER && !out.contains_key("serial") && body.len() >= 58 {
            out.insert("serial".into(), json!(le::u32(body, 54).to_string()));
            out.insert(
                "firmware".into(),
                json!(format!("{}.{:02}", body[11], body[13])),
            );
            if let Some(i) = find(body, b"\x02\x00\x00", 7)
                && i > 0
            {
                let end = (i + 23).min(body.len());
                out.insert(
                    "short_serial".into(),
                    json!(utf16(&body[(i + 3).min(end)..end])),
                );
            }
        } else if typ == T_ANTENNA && !out.contains_key("setup") && body.len() >= 148 {
            let setup = utf16(&body[2..66]);
            let model = if setup.is_empty() {
                String::new()
            } else {
                format!("LEICA {}", setup.split_whitespace().next().unwrap_or(""))
            };
            let ant_raw = &body[0x61..0x61 + 20];
            let ant_raw = ant_raw.split(|&c| c == 0).next().unwrap_or(&[]);
            let antenna = crate::util::py_strip(&crate::util::latin1(ant_raw)).to_owned();
            out.insert("setup".into(), json!(setup));
            out.insert("model".into(), json!(model));
            out.insert("antenna".into(), json!(antenna));
            out.insert("hook_offset".into(), json!(le::f32(body, 131)));
            out.insert(
                "pco".into(),
                json!([le::f32(body, 140), le::f32(body, 144)]),
            );
        } else if typ == T_POSITION && !out.contains_key("llh") && body.len() >= 32 {
            let (lat, lon, h) = (le::f64(body, 8), le::f64(body, 16), le::f64(body, 24));
            if lat.abs() <= std::f64::consts::PI / 2.0 && lon.abs() <= 2.0 * std::f64::consts::PI {
                out.insert("llh".into(), json!([lat, lon, h]));
            }
        } else if typ == T_POINT && !out.contains_key("point") && body.len() >= 67 {
            out.insert("point".into(), json!(utf16(&body[25..57])));
            out.insert("height".into(), json!(le::f32(body, 63)));
        }
    }
    if let Some(h) = out.get("height").and_then(Value::as_f64) {
        let hook = out
            .get("hook_offset")
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
        out.insert(
            "arp_height".into(),
            json!(geomaticist_gnss_core::text::round_to(h + hook, 4)),
        );
    }
    out
}

fn find(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    hay.get(from..)?
        .windows(needle.len())
        .position(|w| w == needle)
        .map(|i| i + from)
}

/// `GPS0 + timedelta(seconds=t)`.
#[must_use]
pub fn t_of(sec: f64) -> geomaticist_gnss_core::time::DateTime {
    gps_week_fsec(0, sec)
}

/// Decodes one .m00 to `outdir/<stem>.<yy>o`, `.<yy>n` (GPS), `.<yy>g` (GLONASS).
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    convert_with(raw, outdir, log, &crate::now_utc_text())
}

/// [`convert`] with the `PGM / RUN BY / DATE` date text given.
#[allow(clippy::cast_possible_truncation)]
pub fn convert_with(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
) -> Result<Converted, DecodeError> {
    let name = raw
        .file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let stem = raw
        .file_stem()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    std::fs::create_dir_all(outdir).map_err(|e| DecodeError::io(outdir, &e))?;
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    if !has_epochs(&data) {
        return Err(DecodeError::Format(format!(
            "{name} holds no observation records (a Captivate job's own .m00 holds only its points)"
        )));
    }
    let mut text = Vec::new();
    let res = obs_writer::write_obs(
        &data, &name, PROGRAM, RUN_BY, date, None, true, "\n", &mut text,
    )
    .map_err(|e| DecodeError::Format(format!("the built-in decoder couldn't read {name} ({e})")))?;
    let obs = outdir.join(format!("{stem}.{:02}o", res.first.year.rem_euclid(100)));
    crate::write_atomic(&obs, &text)?;
    let week = crate::util::days_since_gps0(&res.first).div_euclid(7);
    let navs = nav::write_nav(
        &data,
        outdir,
        &stem,
        res.first.year,
        PROGRAM,
        RUN_BY,
        date,
        Some(week),
    )?;
    let mut info = res.info;
    let mut warnings = Vec::new();
    let mot = motion::motion(&data);
    if mot.is_kinematic() {
        let n = mot.stops().count();
        info.insert("kinematic".into(), json!(true));
        info.insert("stops".into(), json!(n));
        info.insert("moving_s".into(), json!(mot.moving_seconds().round()));
        warnings.push(format!(
            "{name} moves: {n} static occupation(s) and {:.0} min moving (the receiver's own event records). The RINEX marks them with epoch flags 3 (occupation, with its point name and antenna height) and 2 (antenna starts moving). Only the occupations can be processed as static.",
            mot.moving_seconds() / 60.0
        ));
    }
    if navs.is_empty() {
        warnings.push(format!(
            "{name} holds no broadcast ephemerides: put a navigation file (brdc) beside it."
        ));
    }
    if !res.skipped.is_empty() {
        let n = res.skipped.len();
        warnings.push(format!(
            "{name}: {n} observation record{} couldn't be decoded and {} left out (first at byte {}).",
            if n > 1 { "s" } else { "" },
            if n > 1 { "were" } else { "was" },
            crate::t04::group3(res.skipped[0].0)
        ));
    }
    if !info.contains_key("arp_height") {
        warnings.push(format!(
            "{name} has no point record, so the RINEX antenna height is 0: set the antenna height before processing."
        ));
    }
    let g = |k: &str| {
        info.get(k)
            .and_then(Value::as_str)
            .unwrap_or("?")
            .to_owned()
    };
    log(&format!(
        "Built-in Leica decoder: {} #{} fw {}, antenna {}, ARP height {:.4} m; wrote {} ({} epochs){}",
        g("model"),
        g("serial"),
        g("firmware"),
        g("antenna"),
        info.get("arp_height")
            .and_then(Value::as_f64)
            .unwrap_or(0.0),
        obs.file_name()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned()),
        res.epochs,
        if navs.is_empty() {
            String::new()
        } else {
            format!(" and {} navigation file(s)", navs.len())
        }
    ));
    Ok(Converted {
        obs,
        nav: navs,
        warnings,
        info: Value::Object(info),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn framing_stops_at_break() {
        let mut d = vec![0x9c, 0xae, 6, 0, 0x65, 1, 2, 3];
        d.extend([0x9c, 0xae, 2, 0]);
        let r = records(&d);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].1, 0x65);
        assert_eq!(r[0].2, &[0x65, 1]);
    }

    #[test]
    fn utf16_text() {
        assert_eq!(utf16(b"G\0S\x001\x004\0\0\0x\0"), "GS14");
    }
}
