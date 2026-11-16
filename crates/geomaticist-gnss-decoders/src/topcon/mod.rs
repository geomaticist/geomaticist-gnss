//! Topcon / Javad GREIS receiver logs (`.tps`) to RINEX 3.05.
//!
//! A `.tps` file is the receiver's GREIS message stream written to a file ([`frame`]). This decoder reads:
//!
//! - the compact measurement messages `rD` (all observations of an epoch) and `rP` (their signal names), into
//!   observations ([`obs`]); firmware that logs no `rP` gets its signals inferred from the Doppler ratios;
//! - the broadcast ephemerides `GE` / `QE` / `NE` / `LE` / `DE`, GPS UTC `UO` and ionosphere `IO` ([`nav`]);
//! - the text messages for the header: `JP` (file identification with the receiver model), `PM` (receiver
//!   ID, firmware, serial, the reference antenna and its offset), `RD` (date), `PV` (position for the
//!   approximate position).
//!
//! The observation file is written in the layout office software expects ([`rinex`]). A clean-room decoder
//! written from sample receiver logs, with the public GREIS framing.

pub mod frame;
pub mod greis;
pub mod nav;
pub mod obs;
pub mod rinex;

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::time::{DateTime, TICKS};
use serde_json::json;

use crate::util::le;
use crate::{Converted, DecodeError};
use frame::Messages;
use nav::{Eph, GpsUtc, NavContext};
use obs::{Dropped, Epoch, Inference};

/// PGM of the files this decoder writes.
pub const PROGRAM: &str = "geomaticist-gnss TPS";
/// RUN BY.
pub const RUN_BY: &str = "geomaticist-gnss";

/// File extensions of GREIS logs.
pub const EXTS: [&str; 1] = ["tps"];

/// Whether `data` (the start of a file is enough, 64 KiB or more) is a GREIS log this decoder reads: its first
/// message frames and passes its check, and an `rD` measurement message with a good CRC follows.
#[must_use]
pub fn recognises(data: &[u8]) -> bool {
    let start = data
        .iter()
        .position(|&c| c != b'\r' && c != b'\n')
        .unwrap_or(0);
    if frame::message_at(data, start).is_none() {
        return false;
    }
    let head = &data[..data.len().min(256 * 1024)];
    Messages::new(head).any(|m| &m.id == b"rD" && obs::split_rd(m.body).is_ok())
}

/// Whether `raw` is a `.tps` log this decoder reads (reads the first 256 KiB).
#[must_use]
pub fn can_convert(raw: &Path) -> bool {
    use std::io::Read as _;
    let mut head = Vec::new();
    std::fs::File::open(raw)
        .and_then(|f| f.take(256 * 1024).read_to_end(&mut head))
        .is_ok()
        && recognises(&head)
}

/// Header facts from the log's text messages.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LogInfo {
    /// Receiver model from the `JP` identification ("HiPer VR").
    pub model: String,
    /// Receiver ID (`rcv/id`).
    pub receiver_id: String,
    /// Firmware (`rcv/ver/main`).
    pub firmware: String,
    /// Serial number (`rcv/sn`).
    pub serial: String,
    /// Antenna name with radome, as the receiver's reference settings name it (20 characters).
    pub antenna: String,
    /// Antenna offset: up, east, north (m), from the same settings (assumed: the second of their three triples).
    pub antenna_offset: Option<[f64; 3]>,
    /// Date of the `RD` message (year, month, day).
    pub date: Option<(i64, i64, i64)>,
    /// GPS - UTC from `UO`.
    pub leap: Option<i64>,
}

fn quoted(text: &str, key: &str) -> Option<String> {
    let at = text.find(key)? + key.len();
    let rest = text[at..].strip_prefix("=\"")?;
    let end = rest.find('"')?;
    Some(rest[..end].to_owned())
}

/// The model in a `JP` text: "RLOGF JPS HiPer VR Receiver Log File" -> "HiPer VR".
fn jp_model(body: &[u8]) -> String {
    let t = String::from_utf8_lossy(body);
    let t = t.trim();
    let t = t.strip_prefix("RLOGF").unwrap_or(t).trim_start();
    let t = t
        .strip_prefix("JPS")
        .or_else(|| t.strip_prefix("TPS"))
        .unwrap_or(t)
        .trim_start();
    let end = t.find(" Receiver").unwrap_or(t.len());
    t[..end].trim().to_owned()
}

fn parse_triple(s: &str) -> Option<([f64; 3], usize)> {
    let s2 = s.strip_prefix('{')?;
    let end = s2.find('}')?;
    let v: Vec<f64> = s2[..end]
        .split(',')
        .map(|x| x.trim().parse::<f64>().ok())
        .collect::<Option<_>>()?;
    (v.len() == 3).then(|| ([v[0], v[1], v[2]], end + 2))
}

/// The antenna name and offset in the `/par={ref=...}` settings text: the first `{{a,b,c},{d,e,f},{g,h,i},NAME,`
/// (assumed: the second triple is the reference antenna's offset from the mark, east, north, up).
fn ref_antenna(text: &str) -> Option<(String, [f64; 3])> {
    let mut from = 0;
    while let Some(k) = text[from..].find("{{") {
        let at = from + k + 1;
        from = at;
        let mut p = at;
        let mut triples = Vec::new();
        while triples.len() < 3 {
            let Some((t, n)) = parse_triple(&text[p..]) else {
                break;
            };
            triples.push(t);
            p += n;
            if text[p..].starts_with(',') {
                p += 1;
            } else {
                break;
            }
        }
        if triples.len() < 3 {
            continue;
        }
        let rest = &text[p..];
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        let name = &rest[..end];
        if name.is_empty() || name.starts_with('{') {
            continue;
        }
        let t = triples[1];
        return Some((name.to_owned(), [t[2], t[0], t[1]]));
    }
    None
}

/// Joins the `PM` text messages into whole settings texts (a message starting with spaces continues the last).
fn pm_texts(pms: &[&[u8]]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for b in pms {
        let t = String::from_utf8_lossy(b);
        match out.last_mut() {
            Some(last) if t.starts_with(' ') => last.push_str(t.trim()),
            _ => out.push(t.into_owned()),
        }
    }
    out
}

impl LogInfo {
    fn read_pm(&mut self, texts: &[String]) {
        for t in texts {
            if let Some(v) = quoted(t, "rcv/id") {
                self.receiver_id = v;
            }
            if let Some(v) = quoted(t, "rcv/ver/main") {
                self.firmware = v;
            }
            if let Some(v) = quoted(t, "rcv/sn") {
                self.serial = v;
            }
            if self.antenna.is_empty()
                && t.starts_with("/par={ref=")
                && let Some((name, off)) = ref_antenna(t)
            {
                self.antenna = name;
                self.antenna_offset = Some(off);
            }
        }
    }
}

/// A navigation message kept for decoding after the log is read (the leap seconds may come later).
#[derive(Debug, Clone)]
struct NavMsg {
    id: [u8; 2],
    body: Vec<u8>,
    t_gps: Option<f64>,
}

/// What a pass over a log found besides the epochs.
#[derive(Debug, Clone, Default)]
pub struct Scan {
    /// Header facts.
    pub info: LogInfo,
    /// Navigation records, logged repeats removed.
    pub nav: Vec<Eph>,
    /// GPS ionosphere.
    pub gps_iono: Option<([f64; 4], [f64; 4])>,
    /// BeiDou ionosphere.
    pub bds_iono: Option<([f64; 4], [f64; 4])>,
    /// GPS UTC parameters.
    pub gps_utc: Option<GpsUtc>,
    /// GLONASS - UTC (tau_c, s) and the GLONASS reception time (GPS week, s of week) it was logged at.
    pub glut: Option<(f64, i64, i64)>,
    /// Mean of the receiver's own positions (`PV` with a solution).
    pub approx: Option<[f64; 3]>,
    /// What decoding the epochs left out.
    pub dropped: Dropped,
    /// Epochs that could not be decoded.
    pub bad_epochs: usize,
    /// Epochs.
    pub epochs: usize,
    /// Bytes skipped (damage or unknown framing).
    pub skipped: usize,
    /// Places bytes were skipped.
    pub gaps: usize,
    /// Messages per ID.
    pub ids: BTreeMap<String, usize>,
    /// Whether the log has `rP` signal lists.
    pub has_rp: bool,
}

/// GPS week of a civil date.
fn week_of_date(y: i64, m: i64, d: i64) -> Option<i64> {
    Some(crate::util::date_us(y, m, d)?.div_euclid(crate::util::WEEK_US))
}

/// Reads a log: calls `on_epoch` with every epoch (time in ticks since 1980-01-06, GPS time) in file order, and
/// returns the rest. `ref_week` is the GPS week to resolve the receiver's 10-bit week with when the log has
/// no `RD` date (the current week).
pub fn scan(data: &[u8], ref_week: i64, on_epoch: &mut dyn FnMut(i64, &Epoch)) -> Scan {
    scan_with(data, ref_week, false, on_epoch)
}

/// [`scan`], optionally ignoring the `rP` signal lists (`ignore_rp`: the signals are inferred as for a log
/// without them; for checking the inference against logs that have them).
pub fn scan_with(
    data: &[u8],
    ref_week: i64,
    ignore_rp: bool,
    on_epoch: &mut dyn FnMut(i64, &Epoch),
) -> Scan {
    let mut sc = Scan::default();
    let mut pms: Vec<&[u8]> = Vec::new();
    let mut navs: Vec<NavMsg> = Vec::new();
    let mut inf = Inference::default();
    let mut pending: Option<&[u8]> = None;
    let mut t_last: Option<f64> = None;
    let mut week_ref: Option<i64> = None;
    let mut pos_sum = [0.0; 3];
    let mut pos_n = 0usize;
    let mut it = Messages::new(data);
    let mut emit = |rd: &[u8],
                    rp: Option<&[u8]>,
                    sc: &mut Scan,
                    inf: &mut Inference,
                    week_ref: Option<i64>,
                    t_last: &mut Option<f64>| {
        match obs::decode_epoch(rd, rp, inf, &mut sc.dropped) {
            Ok(e) => {
                let week =
                    nav::expand_week(i64::from(e.week10), 1024, week_ref.unwrap_or(ref_week));
                let ms = week * 604_800_000 + i64::from(e.tow_ms);
                let ticks = ms * (TICKS / 1000);
                *t_last = Some(ms as f64 / 1000.0);
                sc.epochs += 1;
                on_epoch(ticks, &e);
            }
            Err(_) => sc.bad_epochs += 1,
        }
    };
    for m in it.by_ref() {
        *sc.ids.entry(m.id_str().to_owned()).or_default() += 1;
        if &m.id == b"rP" {
            sc.has_rp = true;
            if let Some(rd) = pending.take() {
                let rp = (!ignore_rp).then_some(m.body);
                emit(rd, rp, &mut sc, &mut inf, week_ref, &mut t_last);
            }
            continue;
        }
        if let Some(rd) = pending.take() {
            emit(rd, None, &mut sc, &mut inf, week_ref, &mut t_last);
        }
        match &m.id {
            b"rD" => pending = Some(m.body),
            b"JP" => sc.info.model = jp_model(m.body),
            b"PM" => pms.push(m.body),
            b"RD" if m.body.len() >= 5 && sc.info.date.is_none() => {
                let (y, mo, d) = (
                    le::u16(m.body, 0),
                    i64::from(m.body[2]),
                    i64::from(m.body[3]),
                );
                sc.info.date = Some((y, mo, d));
                week_ref = week_of_date(y, mo, d);
            }
            b"PV" if m.body.len() >= 45 => {
                if m.body[44] > 0 {
                    for (k, s) in pos_sum.iter_mut().enumerate() {
                        *s += le::f64(m.body, 8 * k);
                    }
                    pos_n += 1;
                }
            }
            b"GE" | b"QE" | b"NE" | b"LE" | b"DE" | b"UO" | b"IO" => navs.push(NavMsg {
                id: m.id,
                body: m.body.to_vec(),
                t_gps: t_last,
            }),
            _ => {}
        }
    }
    if let Some(rd) = pending.take() {
        emit(rd, None, &mut sc, &mut inf, week_ref, &mut t_last);
    }
    sc.skipped = it.skipped;
    sc.gaps = it.gaps;
    let texts = pm_texts(&pms);
    sc.info.read_pm(&texts);
    if pos_n > 0 {
        sc.approx = Some(pos_sum.map(|s| s / pos_n as f64));
    }
    let week = week_ref.unwrap_or(ref_week);
    for n in &navs {
        if &n.id == b"UO"
            && let Some(u) = nav::gps_utc(&n.body, week)
        {
            sc.info.leap.get_or_insert(u.dtls);
            sc.gps_utc = Some(u);
        }
    }
    let leap = sc.info.leap.unwrap_or(18);
    let t0 = week as f64 * 604_800.0;
    let mut seen = std::collections::HashSet::new();
    for n in &navs {
        let t_gps = n.t_gps.unwrap_or(t0);
        let ctx = NavContext {
            week: (t_gps / 604_800.0).floor() as i64,
            t_gps,
            leap,
        };
        let eph = match &n.id {
            b"GE" => nav::gps_like(&n.body, false, &ctx),
            b"QE" => nav::gps_like(&n.body, true, &ctx),
            b"NE" => {
                if n.body.len() >= 94 {
                    let tau_c = le::f64(&n.body, 63);
                    let tow = (t_gps - ctx.week as f64 * 604_800.0).floor() as i64;
                    sc.glut = Some((tau_c, tow, ctx.week));
                }
                nav::glonass(&n.body, &ctx)
            }
            b"LE" => nav::galileo(&n.body, &ctx),
            b"DE" => {
                if let Some(io) = nav::beidou_iono(&n.body) {
                    sc.bds_iono = Some(io);
                }
                nav::beidou(&n.body)
            }
            b"IO" => {
                if let Some(io) = nav::gps_iono(&n.body) {
                    sc.gps_iono = Some(io);
                }
                None
            }
            _ => None,
        };
        if let Some(e) = eph
            && seen.insert(e.key())
        {
            sc.nav.push(e);
        }
    }
    sc
}

/// Decodes `raw` into `outdir`: `<stem>.<yy>o` (RINEX 3.05 observations) and `<stem>.<yy>p` (mixed
/// navigation), with today's date in the headers.
pub fn convert(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
) -> Result<Converted, DecodeError> {
    convert_with(raw, outdir, log, &crate::now_utc_text(), None)
}

/// [`convert`] with the `PGM / RUN BY / DATE` date text given, and the GPS week to resolve the receiver's
/// 10-bit week with when the log has no date (`None`: the current week).
pub fn convert_with(
    raw: &Path,
    outdir: &Path,
    log: &mut dyn FnMut(&str),
    date: &str,
    ref_week: Option<i64>,
) -> Result<Converted, DecodeError> {
    let data = std::fs::read(raw).map_err(|e| DecodeError::io(raw, &e))?;
    let name = raw
        .file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let stem = raw
        .file_stem()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    if !recognises(&data) {
        return Err(DecodeError::Format(format!(
            "{name}: not a GREIS (.tps) log with rD measurements"
        )));
    }
    let ref_week = ref_week.unwrap_or_else(|| {
        crate::now_utc()
            .micros_since_gps0()
            .div_euclid(crate::util::WEEK_US)
    });
    log(&format!("{name}: reading the log"));
    let mut st = rinex::Stats::default();
    let sc = scan(&data, ref_week, &mut |t, e| st.add(e, t));
    let Some(first) = st.first else {
        return Err(DecodeError::Format(format!(
            "{name}: no epoch could be decoded"
        )));
    };
    let t0 = DateTime::from_micros_since_gps0(first / (TICKS / 1_000_000));
    let yy = t0.year.rem_euclid(100);
    let info = &sc.info;
    let mut comments = vec![format!("Decoded from {name} by {RUN_BY}")];
    if !sc.has_rp {
        comments.push("Signals inferred from Doppler ratios (the log has no rP)".to_owned());
    }
    let hdr = rinex::ObsHeader {
        program: PROGRAM.into(),
        run_by: RUN_BY.into(),
        date: date.into(),
        comments,
        marker: stem.clone(),
        receiver_number: info.receiver_id.clone(),
        receiver_type: if info.model.is_empty() {
            "TPS".into()
        } else {
            format!("TPS {}", info.model)
        },
        receiver_version: info.firmware.clone(),
        antenna: info.antenna.clone(),
        approx: sc.approx.unwrap_or([0.0; 3]),
        delta_hen: info.antenna_offset.unwrap_or([0.0; 3]),
        leap: info.leap,
    };
    let mut text = rinex::header_text(&st, &hdr);
    let mut lli = rinex::Lli::default();
    let _ = scan(&data, ref_week, &mut |t, e| {
        rinex::epoch_text(&mut text, &st, e, t, true, &mut lli);
    });
    let obs = outdir.join(format!("{stem}.{yy:02}o"));
    crate::write_atomic(&obs, text.as_bytes())?;
    let mut nav_paths = Vec::new();
    if !sc.nav.is_empty() {
        let nh = rinex::NavHeader {
            program: PROGRAM.into(),
            run_by: RUN_BY.into(),
            date: date.into(),
            gps_iono: sc.gps_iono,
            bds_iono: sc.bds_iono,
            gps_utc: sc.gps_utc,
            glut: sc.glut,
            leap: info.leap,
        };
        let p = outdir.join(format!("{stem}.{yy:02}p"));
        crate::write_atomic(&p, rinex::nav_text(&sc.nav, &nh).as_bytes())?;
        nav_paths.push(p);
    }
    let mut warnings = Vec::new();
    if sc.skipped > 0 {
        warnings.push(format!(
            "{name}: {} bytes in {} places could not be read and were skipped",
            sc.skipped, sc.gaps
        ));
    }
    if sc.bad_epochs > 0 {
        warnings.push(format!(
            "{name}: {} measurement epochs could not be decoded",
            sc.bad_epochs
        ));
    }
    let d = &sc.dropped;
    if d.unknown_sats + d.unknown_signals + d.glonass_no_slot + d.rp_mismatch > 0 {
        warnings.push(format!(
            "{name}: left out {} satellites of unknown index, {} GLONASS without slot, {} unknown signals; {} epochs with an unusable rP",
            d.unknown_sats, d.glonass_no_slot, d.unknown_signals, d.rp_mismatch
        ));
    }
    log(&format!(
        "{name}: {} {} fw {}, antenna {}, {} epochs, {} navigation records",
        hdr.receiver_type,
        info.receiver_id,
        info.firmware,
        if info.antenna.is_empty() {
            "?"
        } else {
            info.antenna.trim()
        },
        sc.epochs,
        sc.nav.len()
    ));
    let info_json = json!({
        "receiver": hdr.receiver_type,
        "receiver_id": info.receiver_id,
        "serial": info.serial,
        "firmware": info.firmware,
        "antenna": info.antenna.trim(),
        "antenna_offset": info.antenna_offset,
        "epochs": sc.epochs,
        "interval": st.interval(),
        "signals_from": if sc.has_rp { "rP" } else { "Doppler ratios" },
        "skipped_bytes": sc.skipped,
        "bad_epochs": sc.bad_epochs,
        "dropped": {
            "unknown_sats": d.unknown_sats,
            "glonass_no_slot": d.glonass_no_slot,
            "unknown_signals": d.unknown_signals,
            "rp_mismatch": d.rp_mismatch,
        },
        "messages": sc.ids,
    });
    Ok(Converted {
        obs,
        nav: nav_paths,
        warnings,
        info: info_json,
    })
}

#[cfg(test)]
mod tests;
