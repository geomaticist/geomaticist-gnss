//! RINEX 3.04 observations from a Leica MDB log, with the values office software expects
//!.
//!
//! The rules:
//! - phase only when the signal's phase-valid flag is set; LLI 1 when its loss-of-lock flag is;
//! - the signal-strength digit on the phase only: int((dB-Hz + 0.5) / 6), 1 ... 9;
//! - 3 decimals, half away from zero, except a Doppler that falls exactly on a half: the direction then
//!   depends on the value's power of two (`DOPPLER_TIE_TOWARD_ZERO`).
//!
//! Header: receiver, antenna, ARP height and first position from the log's own records.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use geomaticist_gnss_core::time::DateTime;
use serde_json::{Map, Value};

use super::epochs::{self, Epoch, Sat, Signal};
use crate::util::{hline, latin1_bytes, sec_f};

const TYPE_ORDER: &str = "CLDS";
const SYSTEM_ORDER: &str = "GRECJ";
/// floor(log2 |Doppler Hz|) for which an exact half is rounded toward zero.
const DOPPLER_TIE_TOWARD_ZERO: [u64; 5] = [2, 4, 7, 8, 9];

/// GPS - UTC leap seconds (the writer's table).
#[must_use]
pub fn leap_seconds(t: &DateTime) -> i64 {
    let k = (t.year, t.month, t.day);
    if k >= (2017, 1, 1) {
        18
    } else if k >= (2015, 7, 1) {
        17
    } else if k >= (2012, 7, 1) {
        16
    } else {
        15
    }
}

/// Thousandths as text, with floor division (`k // 1000`, then `k mod 1000` in 3 digits).
fn milli(k: i64, neg: bool) -> String {
    format!(
        "{}{}.{:03}",
        if neg && k != 0 { "-" } else { "" },
        k.div_euclid(1000),
        k.rem_euclid(1000)
    )
}

/// raw * 2^-10 Hz to 3 decimals: thousandths = raw * 125 / 128, half away from zero; with `lgo3`, an exact half
/// goes by `DOPPLER_TIE_TOWARD_ZERO` (the RINEX 3 rule).
#[must_use]
pub fn doppler_text(raw: i64, lgo3: bool) -> String {
    let a = raw.unsigned_abs();
    let q = i128::from(a) * 125;
    let mut k = (q + 64).div_euclid(128);
    let r = (q + 64).rem_euclid(128);
    let bits = u64::from(64 - a.leading_zeros());
    if lgo3 && r == 0 && bits >= 11 && DOPPLER_TIE_TOWARD_ZERO.contains(&(bits - 11)) {
        k -= 1;
    }
    milli(i64::try_from(k).unwrap_or(0), raw < 0)
}

/// Phase, cycles: (range + (phase - range)) / wavelength in doubles.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn phase_text(range_raw: i64, phase_raw: i64, freq: f64) -> String {
    format!(
        "{:.3}",
        (range_raw as f64 * 0.02 + phase_raw as f64 * 0.0005) / (epochs::C_LIGHT as f64 / freq)
    )
}

/// `{type3: (text, lli, ssi)}` of one satellite, in insertion order with dict overwrite.
type SatObs = Vec<(String, (String, i64, i64))>;

fn put(o: &mut SatObs, k: String, v: (String, i64, i64)) {
    match o.iter_mut().find(|e| e.0 == k) {
        Some(e) => e.1 = v,
        None => o.push((k, v)),
    }
}

/// `{sat name: {type3: (text, lli, ssi)}}` of one epoch, as written.
#[allow(clippy::cast_precision_loss)]
fn observations(ep: &Epoch, lgo3: bool) -> Result<Vec<(String, SatObs)>, String> {
    let mut out: Vec<(String, SatObs)> = Vec::new();
    for sat in &ep.sats {
        let mut obs = SatObs::new();
        for sig in &sat.signals {
            let Some(rinex) = sig.rinex else { continue };
            put(
                &mut obs,
                format!("C{rinex}"),
                (milli(sig.range_raw * 20, false), 0, 0),
            );
            if sig.has_phase() {
                let f = freq(sat, sig, rinex)?;
                put(
                    &mut obs,
                    format!("L{rinex}"),
                    (
                        phase_text(sig.range_raw, sig.phase_raw, f),
                        i64::from(sig.lli()),
                        ((sig.snr_raw + 10).div_euclid(120)).clamp(1, 9),
                    ),
                );
            }
            put(
                &mut obs,
                format!("D{rinex}"),
                (doppler_text(sig.doppler_raw, lgo3), 0, 0),
            );
            put(
                &mut obs,
                format!("S{rinex}"),
                (milli(sig.snr_raw * 50, false), 0, 0),
            );
        }
        if !obs.is_empty() {
            let name = sat.name();
            match out.iter_mut().find(|e| e.0 == name) {
                Some(e) => e.1 = obs,
                None => out.push((name, obs)),
            }
        }
    }
    Ok(out)
}

#[allow(clippy::cast_precision_loss)]
fn freq(sat: &Sat, _sig: &Signal, rinex: &str) -> Result<f64, String> {
    epochs::frequency(sat.system, rinex, sat.channel)
        .map(|f| f as f64)
        .ok_or_else(|| format!("no frequency for {} {rinex}", sat.system))
}

fn sort_types(types: &BTreeSet<String>) -> Vec<String> {
    let mut v: Vec<String> = types.iter().cloned().collect();
    v.sort_by_key(|t| {
        let b = t.as_bytes();
        (
            b.get(1).copied(),
            TYPE_ORDER.find(char::from(b[0])).unwrap_or(9),
            b.get(2).copied(),
        )
    });
    v
}

fn sat_key(s: &str) -> (usize, String) {
    (
        SYSTEM_ORDER
            .find(s.chars().next().unwrap_or('?'))
            .unwrap_or(99),
        s.to_owned(),
    )
}

/// What the writer reports.
#[derive(Debug, Clone)]
pub struct ObsResult {
    /// Epochs decoded.
    pub epochs: usize,
    /// Records that didn't decode: (offset, reason).
    pub skipped: Vec<(usize, String)>,
    /// Header facts.
    pub info: Map<String, Value>,
    /// First epoch.
    pub first: DateTime,
    /// Last epoch.
    pub last: DateTime,
}

/// Decodes the log (header first: types, counts, first / last epoch; then the body) and writes RINEX 3.04
/// observations into `out`.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn write_obs(
    data: &[u8],
    raw_name: &str,
    program: &str,
    run_by: &str,
    date: &str,
    marker: Option<&str>,
    lgo3: bool,
    line_end: &str,
    out: &mut Vec<u8>,
) -> Result<ObsResult, String> {
    let info = super::header_info(data);
    let mut errors = Vec::new();
    let mut n_epochs = 0usize;
    let mut present: BTreeMap<char, BTreeSet<String>> = BTreeMap::new();
    let mut counts: BTreeMap<String, BTreeMap<String, i64>> = BTreeMap::new();
    let mut glo: BTreeMap<String, i64> = BTreeMap::new();
    let mut first: Option<f64> = None;
    let mut last: Option<f64> = None;
    let mut steps: Vec<f64> = Vec::new();
    for (pos, typ, body) in super::records(data) {
        if typ != epochs::T_EPOCH {
            continue;
        }
        let ep = match epochs::decode(body) {
            Ok(e) => e,
            Err(e) => {
                errors.push((pos, e));
                continue;
            }
        };
        n_epochs += 1;
        let ep = &ep;
        let t = ep.gps_seconds();
        if let Some(l) = last
            && t > l
        {
            steps.push(t - l);
        }
        if first.is_none() {
            first = Some(t);
        }
        last = Some(t);
        for sat in &ep.sats {
            if sat.system == 'R' {
                glo.entry(sat.name()).or_insert(sat.channel.unwrap_or(0));
            }
        }
        let obs = observations(ep, lgo3)?;
        for (s, o) in &obs {
            for (ty, _) in o {
                present
                    .entry(s.chars().next().unwrap_or('?'))
                    .or_default()
                    .insert(ty.clone());
                *counts
                    .entry(s.clone())
                    .or_default()
                    .entry(ty.clone())
                    .or_insert(0) += 1;
            }
        }
    }
    let (Some(first), Some(last)) = (first, last) else {
        return Err(format!("{raw_name} holds no observation records"));
    };
    let types: Vec<(char, Vec<String>)> = SYSTEM_ORDER
        .chars()
        .filter_map(|s| {
            present
                .get(&s)
                .filter(|p| !p.is_empty())
                .map(|p| (s, sort_types(p)))
        })
        .collect();
    let interval = if steps.is_empty() {
        None
    } else {
        steps.sort_by(f64::total_cmp);
        Some(steps[steps.len() / 2])
    };
    let (t0, t1) = (super::t_of(first), super::t_of(last));
    let gs = |k: &str| info.get(k).and_then(Value::as_str).unwrap_or("").to_owned();
    let stem = raw_name.rsplit_once('.').map_or(raw_name, |(a, _)| a);
    let marker_name = marker
        .map(str::to_owned)
        .filter(|m| !m.is_empty())
        .or_else(|| Some(gs("point")).filter(|p| !p.is_empty()))
        .unwrap_or_else(|| stem.to_owned());
    // a log that moves (a rover's job log, or a base that logged a while before it was set): epoch flag 2
    // where the antenna starts moving, flag 3 + MARKER NAME + antenna height where an occupation starts
    // (RINEX 3.04 5.x); a static log is written as before
    let mot = super::motion::motion(data);
    let kin = mot.is_kinematic();
    let rover = kin && (mot.stops().count() > 1 || mot.moving_seconds() > 0.5 * span_s(&mot));
    let hook = info
        .get("hook_offset")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let events = if kin {
        motion_events(&mot, hook, line_end)
    } else {
        Vec::new()
    };
    let mut h = vec![
        hline(
            &format!(
                "{:>9}{:11}{:<20}{:<20}",
                "3.04", "", "OBSERVATION DATA", "M: MIXED"
            ),
            "RINEX VERSION / TYPE",
        ),
        hline(
            &format!("{program:<20.20}{run_by:<20.20}{date}"),
            "PGM / RUN BY / DATE",
        ),
        hline(
            &format!("Decoded from {raw_name} by an open clean-room decoder"),
            "COMMENT",
        ),
    ];
    if kin {
        let n = mot.stops().count();
        h.push(hline(
            &format!("Kinematic log: {n} static occupation(s), moving between"),
            "COMMENT",
        ));
        h.push(hline(
            "Epoch flag 2 = antenna starts moving, 3 = occupation",
            "COMMENT",
        ));
    }
    h.extend([
        hline(&marker_name, "MARKER NAME"),
        // a pole carried between points (RINEX 3.04 MARKER TYPE HUMAN); a base that only moved before it was set
        // stays GEODETIC
        hline(if rover { "HUMAN" } else { "GEODETIC" }, "MARKER TYPE"),
        hline("", "OBSERVER / AGENCY"),
        hline(
            &format!(
                "{:<20.20}{:<20.20}{:<20.20}",
                gs("serial"),
                gs("model"),
                gs("firmware")
            ),
            "REC # / TYPE / VERS",
        ),
        hline(
            &format!("{:<20}{:<16.16}{:<4}", "", gs("antenna"), "NONE"),
            "ANT # / TYPE",
        ),
    ]);
    if let Some(llh) = info.get("llh").and_then(Value::as_array) {
        let g = |i: usize| llh.get(i).and_then(Value::as_f64).unwrap_or(0.0);
        let (x, y, z) = crate::util::llh_rad_to_xyz(g(0), g(1), g(2));
        h.push(hline(
            &format!("{x:14.4}{y:14.4}{z:14.4}"),
            "APPROX POSITION XYZ",
        ));
    }
    let arp = info
        .get("arp_height")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    h.push(hline(
        &format!("{arp:14.4}{:14.4}{:14.4}", 0.0, 0.0),
        "ANTENNA: DELTA H/E/N",
    ));
    for (s, ty) in &types {
        for (k, chunk) in ty.chunks(13).enumerate() {
            let lead = if k == 0 {
                format!("{s}  {:3}", ty.len())
            } else {
                " ".repeat(6)
            };
            let body: String = chunk.iter().map(|t| format!(" {t}")).collect();
            h.push(hline(&(lead + &body), "SYS / # / OBS TYPES"));
        }
    }
    h.push(hline("DBHZ", "SIGNAL STRENGTH UNIT"));
    if let Some(iv) = interval.filter(|v| *v != 0.0) {
        h.push(hline(&format!("{iv:10.3}"), "INTERVAL"));
    }
    for (t, lab) in [(&t0, "TIME OF FIRST OBS"), (&t1, "TIME OF LAST OBS")] {
        h.push(hline(
            &format!(
                "{:6}{:6}{:6}{:6}{:6}{:13.7}     GPS",
                t.year,
                t.month,
                t.day,
                t.hour,
                t.minute,
                sec_f(t)
            ),
            lab,
        ));
    }
    h.push(hline(&format!("{:6}", 0), "RCV CLOCK OFFS APPL"));
    for (s, _) in &types {
        h.push(hline(&s.to_string(), "SYS / PHASE SHIFT"));
    }
    if !glo.is_empty() {
        let items: Vec<(&String, &i64)> = glo.iter().collect();
        for (k, chunk) in items.chunks(8).enumerate() {
            let body: String = chunk.iter().map(|(s, kk)| format!(" {s} {kk:2}")).collect();
            let lead = if k == 0 {
                format!("{:3}", items.len())
            } else {
                "   ".into()
            };
            h.push(hline(&(lead + &body), "GLONASS SLOT / FRQ #"));
        }
        h.push(hline("", "GLONASS COD/PHS/BIS"));
    }
    h.push(hline(&format!("{:6}", leap_seconds(&t0)), "LEAP SECONDS"));
    h.push(hline(&format!("{:6}", counts.len()), "# OF SATELLITES"));
    let mut sats: Vec<&String> = counts.keys().collect();
    sats.sort_by_key(|s| sat_key(s));
    for sat in sats {
        let ty = types
            .iter()
            .find(|t| sat.starts_with(t.0))
            .map_or(&[][..], |t| t.1.as_slice());
        let vals: Vec<i64> = ty
            .iter()
            .map(|t| counts[sat].get(t).copied().unwrap_or(0))
            .collect();
        for (k, chunk) in vals.chunks(9).enumerate() {
            let lead = if k == 0 {
                format!("   {sat}")
            } else {
                " ".repeat(6)
            };
            let body: String = chunk.iter().map(|v| format!("{v:6}")).collect();
            h.push(hline(&(lead + &body), "PRN / # OF OBS"));
        }
    }
    h.push(hline("", "END OF HEADER"));
    let mut head = h.join(line_end);
    head.push_str(line_end);
    out.extend_from_slice(&latin1_bytes(&head));
    let mut line = String::with_capacity(256);
    let mut next_event = 0usize;
    // pass 2: the epochs decoded again, written as they come
    for (_pos, typ, body) in super::records(data) {
        if typ != epochs::T_EPOCH {
            continue;
        }
        let Ok(ep) = epochs::decode(body) else {
            continue;
        };
        let mut obs = observations(&ep, lgo3)?;
        while let Some((te, text)) = events.get(next_event)
            && *te <= ep.gps_seconds() + 1e-6
        {
            out.extend_from_slice(&latin1_bytes(text));
            next_event += 1;
        }
        let t = super::t_of(ep.gps_seconds());
        line.clear();
        let _ = write!(
            line,
            "> {:4} {:02} {:02} {:02} {:02}{:11.7}  0{:3}{line_end}",
            t.year,
            t.month,
            t.day,
            t.hour,
            t.minute,
            sec_f(&t),
            obs.len()
        );
        out.extend_from_slice(line.as_bytes());
        obs.sort_by_key(|a| sat_key(&a.0));
        for (sat, o) in &obs {
            line.clear();
            line.push_str(sat);
            let ty = types
                .iter()
                .find(|t| sat.starts_with(t.0))
                .map_or(&[][..], |t| t.1.as_slice());
            for t in ty {
                match o.iter().find(|e| &e.0 == t) {
                    Some((_, (txt, lli, ssi))) if txt.len() <= 14 => {
                        let _ = write!(line, "{txt:>14}");
                        line.push(if *lli != 0 {
                            char::from(b'0' + u8::try_from(*lli).unwrap_or(0))
                        } else {
                            ' '
                        });
                        line.push(if *ssi != 0 {
                            char::from(b'0' + u8::try_from(*ssi).unwrap_or(0))
                        } else {
                            ' '
                        });
                    }
                    _ => line.push_str("                "),
                }
            }
            let trimmed = line.trim_end_matches(crate::util::py_space);
            out.extend_from_slice(&latin1_bytes(trimmed));
            out.extend_from_slice(line_end.as_bytes());
        }
    }
    Ok(ObsResult {
        epochs: n_epochs,
        skipped: errors,
        info,
        first: t0,
        last: t1,
    })
}

fn span_s(m: &super::motion::LogMotion) -> f64 {
    m.span.map_or(0.0, |(a, b)| b - a)
}

/// The event records of a log that moves, in time order: (GPS seconds, the event's lines).
fn motion_events(m: &super::motion::LogMotion, hook: f64, line_end: &str) -> Vec<(f64, String)> {
    use super::motion::SegmentKind;
    let mut out = Vec::new();
    let mut moved = false;
    for seg in &m.segments {
        let t = super::t_of(seg.start);
        let stamp = format!(
            "> {:4} {:02} {:02} {:02} {:02}{:11.7}",
            t.year,
            t.month,
            t.day,
            t.hour,
            t.minute,
            sec_f(&t)
        );
        match seg.kind {
            SegmentKind::Kinematic => {
                moved = true;
                out.push((seg.start, format!("{stamp}  2  0{line_end}")));
            }
            SegmentKind::Static if moved => {
                let dh = seg.antenna_height().unwrap_or(0.0) + hook;
                out.push((
                    seg.start,
                    format!(
                        "{stamp}  3  2{line_end}{}{line_end}{}{line_end}",
                        hline(&seg.name(), "MARKER NAME"),
                        hline(
                            &format!("{dh:14.4}{:14.4}{:14.4}", 0.0, 0.0),
                            "ANTENNA: DELTA H/E/N"
                        )
                    ),
                ));
            }
            SegmentKind::Static => {}
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doppler_ties() {
        assert_eq!(doppler_text(-283_620, true), "-276.973");
        assert_eq!(doppler_text(1, true), "0.001");
        assert_eq!(doppler_text(0, true), "0.000");
        // 2^11 * 2.5 / ... exact half at bit length 13 (11 + 2): toward zero
        let raw = 4100; // 4100 * 125 = 512500; +64 = 512564 = 4004 * 128 + 52 -> not a tie
        assert_eq!(doppler_text(raw, true), "4.004");
    }

    #[test]
    fn kinematic_log_events() {
        use crate::mdb::motion::synth::job_log;
        let t0 = 1_411_084_800.0;
        let d = job_log(t0, &[(t0 + 10.0, 4.0, "2001"), (t0 + 30.0, 5.0, "")], 5.0);
        let mut out = Vec::new();
        let r = write_obs(
            &d, "job.m00", "P", "R", "D", None, true, "
", &mut out,
        );
        assert!(r.is_ok());
        let text = String::from_utf8_lossy(&out);
        assert!(text.contains("HUMAN               "));
        assert!(text.contains("Kinematic log: 2 static occupation(s)"));
        let lines: Vec<&str> = text.lines().collect();
        let occ: Vec<usize> = (0..lines.len())
            .filter(|&i| lines[i].ends_with("  3  2"))
            .collect();
        assert_eq!(occ.len(), 2);
        assert_eq!(lines.iter().filter(|l| l.ends_with("  2  0")).count(), 3);
        assert!(lines[occ[0]].starts_with("> 2024 09 23 00 00 10.0000000"));
        assert!(
            lines[occ[0] + 1].starts_with("2001 ")
                && lines[occ[0] + 1].trim_end().ends_with("MARKER NAME")
        );
        assert!(lines[occ[0] + 2].starts_with("        1.8000"));
        assert!(lines[occ[1] + 1].starts_with("STOP-267-"));
        // the event comes before the epoch of its time
        assert!(lines[occ[0] + 3].starts_with("> 2024 09 23 00 00 10.0000000  0"));
    }

    #[test]
    fn static_log_has_no_events() {
        use crate::mdb::motion::synth::{epoch, event, frame};
        let t0 = 1_411_084_800.0;
        let mut d = frame(&[crate::mdb::T_HEADER, 0, 0, 0, 0, 0, 0]);
        let at = u32::try_from(d.len()).unwrap_or(0);
        d.extend(frame(&event(
            t0,
            at,
            crate::mdb::motion::EV_STATIC_START,
            t0,
        )));
        for k in 0..20 {
            d.extend(frame(&epoch(t0 + f64::from(k))));
        }
        let mut out = Vec::new();
        assert!(
            write_obs(
                &d, "base.m00", "P", "R", "D", None, true, "
", &mut out
            )
            .is_ok()
        );
        let text = String::from_utf8_lossy(&out);
        assert!(text.contains("GEODETIC            "));
        assert!(!text.contains("Kinematic"));
        assert!(
            !text
                .lines()
                .any(|l| l.ends_with("  3  2") || l.ends_with("  2  0"))
        );
    }

    #[test]
    fn milli_text() {
        assert_eq!(milli(20_416_285_680, false), "20416285.680");
        assert_eq!(milli(-5, true), "--1.995"); // sign prefix on top of floor division
        assert_eq!(milli(-5, false), "-1.995");
    }
}
