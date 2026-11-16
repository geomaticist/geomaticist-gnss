//! RINEX 3.05 observation writer for decoded T04 data.
//!
//! Writes the RINEX office software expects from a T04. The options: Doppler, SNR, clock offsets (epoch-line receiver clock),
//! continuous clock ("fix msec steps").
//!
//! Trimble-compatible filtering (`trimble_filter`, default):
//! - a signal is written only if measFlags has bits 0x001 and 0x400 set;
//! - its phase is written only if measFlags bit 0x008 is set (else C/D/S with a blank L);
//! - a satellite is skipped for the epoch if its FIRST signal block is present but invalid.
//!
//! LLI: bit 0 = 1 if measFlags bit 0x100 is set, or the signal's own phase was not valid at the
//! previous epoch. SSI = min(9, max(1, floor(dBHz / 6))) on C, L and D. Pseudoranges are printed from 15
//! significant digits, then 3 decimals with ties away from zero; phase, Doppler and the clock are exact
//! fractions rounded with ties away from zero.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;

use geomaticist_gnss_core::time::DateTime;

use super::container::Inner;
use super::header::HeaderInfo;
use super::measurements::{self, Obs};
use super::{SatId, T04Error};
use crate::decimal::{round_float_3, round_ratio};
use crate::util::{gps_week_ms, hline, sec_f};

/// Header order of SYS / # / OBS TYPES.
const SYSTEM_ORDER: &str = "CRGJES";
/// Systems with line-padding memory (G and C; not R).
const PAD_SYSTEMS: &str = "GC";

/// Writer options.
#[derive(Debug, Clone)]
pub struct ObsOptions {
    /// "Log Doppler observations".
    pub doppler: bool,
    /// "Log SNR observations".
    pub snr: bool,
    /// "Log clock offsets" (epoch-line receiver clock).
    pub clock_offsets: bool,
    /// "Log continuous clock (fix msec steps)".
    pub continuous_clock: bool,
    /// Trimble-compatible write rules.
    pub trimble_filter: bool,
    /// PGM.
    pub program: String,
    /// RUN BY.
    pub run_by: String,
    /// The `PGM / RUN BY / DATE` date text (`"20261104 093000 UTC"`).
    pub date: String,
    /// OBSERVER.
    pub observer: String,
    /// AGENCY.
    pub agency: String,
    /// MARKER NAME override.
    pub marker_name: Option<String>,
    /// MARKER NUMBER override.
    pub marker_number: Option<String>,
    /// Antenna height override, m.
    pub antenna_height: Option<f64>,
    /// Line end.
    pub line_end: String,
}

impl Default for ObsOptions {
    fn default() -> Self {
        Self {
            doppler: true,
            snr: true,
            clock_offsets: false,
            continuous_clock: false,
            trimble_filter: true,
            program: "t04dec 0.3".into(),
            run_by: "t04dec".into(),
            date: crate::now_utc_text(),
            observer: String::new(),
            agency: String::new(),
            marker_name: None,
            marker_number: None,
            antenna_height: None,
            line_end: "\r\n".into(),
        }
    }
}

/// GPS - UTC leap seconds at a GPS date (the writer's own table).
#[must_use]
pub fn leap_seconds(t: &DateTime) -> i64 {
    let k = (t.year, t.month, t.day);
    if k >= (2017, 1, 1) {
        18
    } else if k >= (2015, 7, 1) {
        17
    } else if k >= (2012, 7, 1) {
        16
    } else if k >= (2009, 1, 1) {
        15
    } else if k >= (2006, 1, 1) {
        14
    } else {
        13
    }
}

/// SSI from dB-Hz.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn ssi(dbhz: f64) -> i64 {
    ((dbhz / 6.0) as i64).clamp(1, 9)
}

/// The kept signals of one satellite: `(code, payload, phase written)`, last one of a code wins on lookup.
type Kept<T> = Vec<(&'static str, T, bool)>;

/// Trimble-compatible write rules over the known-signal blocks `(sat, code, measFlags, payload)` in receiver order.
fn select<T: Copy>(
    blocks: impl Iterator<Item = (SatId, &'static str, i64, T)>,
    trimble_filter: bool,
    out: &mut Vec<(SatId, Kept<T>)>,
) {
    out.clear();
    let mut by_sat: Vec<(SatId, Vec<(&'static str, i64, T)>)> = Vec::new();
    for (sat, code, f, o) in blocks {
        match by_sat.iter_mut().find(|e| e.0 == sat) {
            Some(e) => e.1.push((code, f, o)),
            None => by_sat.push((sat, vec![(code, f, o)])),
        }
    }
    for (sat, obs) in by_sat {
        if trimble_filter && (obs[0].1 & 0x401) != 0x401 {
            continue;
        }
        let keep: Kept<T> = obs
            .iter()
            .filter(|(_, f, _)| !trimble_filter || (f & 0x401) == 0x401)
            .map(|&(code, f, o)| (code, o, !trimble_filter || f & 0x8 != 0))
            .collect();
        if !keep.is_empty() {
            out.push((sat, keep));
        }
    }
}

/// The codes of a satellite: first-occurrence order, the last occurrence's value.
fn unique<T>(kept: &Kept<T>) -> Vec<(&'static str, bool)> {
    let mut out: Vec<(&'static str, bool)> = Vec::new();
    for k in kept {
        match out.iter_mut().find(|e| e.0 == k.0) {
            Some(e) => e.1 = k.2,
            None => out.push((k.0, k.2)),
        }
    }
    out
}

fn lookup<'a, T>(kept: &'a Kept<T>, code: &str) -> Option<&'a (&'static str, T, bool)> {
    kept.iter().rev().find(|k| k.0 == code)
}

/// `fmt_val` for an already-rounded value text.
fn push_val(line: &mut String, txt: &str, lli: i64, s: Option<i64>) {
    if txt.len() > 14 {
        // RINEX field overflow: leave blank rather than corrupt the line
        line.push_str("                ");
        return;
    }
    let _ = write!(line, "{txt:>14}");
    line.push(if lli != 0 {
        char::from(b'0' + u8::try_from(lli % 10).unwrap_or(0))
    } else {
        ' '
    });
    match s {
        Some(s) => line.push(char::from(
            b'0' + u8::try_from(s.rem_euclid(10)).unwrap_or(0),
        )),
        None => line.push(' '),
    }
}

struct Pass1 {
    obs_types: Vec<(char, Vec<String>)>,
    counts: BTreeMap<SatId, HashMap<String, i64>>,
    glo_k: BTreeMap<String, i64>,
    first: (i64, i64),
    last: (i64, i64),
    n_epochs: usize,
}

fn pass1(inner: &Inner, hi: &HeaderInfo, opt: &ObsOptions) -> Result<Pass1, T04Error> {
    let types_wanted: Vec<char> = [
        Some('C'),
        opt.doppler.then_some('D'),
        Some('L'),
        opt.snr.then_some('S'),
    ]
    .into_iter()
    .flatten()
    .collect();
    let mut present: Vec<(char, BTreeSet<String>)> = Vec::new();
    let mut counts: BTreeMap<SatId, HashMap<String, i64>> = BTreeMap::new();
    let mut glo_k: BTreeMap<String, i64> = BTreeMap::new();
    let mut first = None;
    let mut last = (0, 0);
    let mut n_epochs = 0usize;
    let mut sigs = Vec::new();
    let mut sel = Vec::new();
    measurements::for_each_raw_epoch(inner, |_kind, r| {
        n_epochs += 1;
        let wt = (r.week()?, r.receive_time()?);
        if first.is_none() {
            first = Some(wt);
        }
        last = wt;
        let sattype = r.sattype()?;
        let fdma = r.fdma()?;
        for (i, &sv) in r.sv_id()?.iter().enumerate() {
            if sattype.get(i).copied() == Some(2) {
                let k = fdma.get(i).copied().ok_or_else(|| {
                    T04Error::Decode("IndexError: list index out of range".into())
                })?;
                glo_k.entry(format!("R{sv:02}")).or_insert(k);
            }
        }
        measurements::epoch_signals(r, &mut sigs)?;
        select(
            sigs.iter()
                .filter_map(|&(sat, code, f)| code.map(|c| (sat, c, f, ()))),
            opt.trimble_filter,
            &mut sel,
        );
        for (sat, kept) in &sel {
            for (code, ph) in unique(kept) {
                for &t in &types_wanted {
                    if t == 'L' && !ph {
                        continue;
                    }
                    let name = format!("{t}{code}");
                    let sys = sat.sys();
                    match present.iter_mut().find(|p| p.0 == sys) {
                        Some(p) => {
                            p.1.insert(name.clone());
                        }
                        None => present.push((sys, BTreeSet::from([name.clone()]))),
                    }
                    *counts.entry(*sat).or_default().entry(name).or_insert(0) += 1;
                }
            }
        }
        Ok(())
    })?;
    let Some(first) = first else {
        return Err(T04Error::Decode("no measurement epochs in the file".into()));
    };
    let mut obs_types: Vec<(char, Vec<String>)> = Vec::new();
    for s in SYSTEM_ORDER.chars() {
        if let Some(p) = present.iter().find(|p| p.0 == s) {
            obs_types.push((s, p.1.iter().cloned().collect()));
        }
    }
    for p in &present {
        if !obs_types.iter().any(|o| o.0 == p.0) {
            obs_types.push((p.0, p.1.iter().cloned().collect()));
        }
    }
    let count_names: HashSet<String> = counts.keys().map(ToString::to_string).collect();
    glo_k.retain(|k, _| count_names.contains(k));
    for (&slot, &k) in &hi.glo_nav_k {
        // also list a satellite that has an ephemeris but no observations written
        glo_k.entry(format!("R{slot:02}")).or_insert(k);
    }
    Ok(Pass1 {
        obs_types,
        counts,
        glo_k,
        first,
        last,
        n_epochs,
    })
}

fn header_text(raw_name: &str, hi: &HeaderInfo, opt: &ObsOptions, p: &Pass1) -> (String, DateTime) {
    let first = gps_week_ms(p.first.0, p.first.1); // uncorrected times
    let last = gps_week_ms(p.last.0, p.last.1);
    let mut l = vec![
        hline(
            &format!(
                "{:>9}{:11}{:<20}{:<20}",
                "3.05", "", "OBSERVATION DATA", "Mixed(MIXED)"
            ),
            "RINEX VERSION / TYPE",
        ),
        hline(
            &format!("{:<20.20}{:<20.20}{}", opt.program, opt.run_by, opt.date),
            "PGM / RUN BY / DATE",
        ),
        hline(
            &format!("Decoded from {raw_name} by an open clean-room decoder"),
            "COMMENT",
        ),
        hline(
            opt.marker_name.as_deref().unwrap_or(&hi.marker_name),
            "MARKER NAME",
        ),
        hline(
            opt.marker_number.as_deref().unwrap_or(&hi.marker_number),
            "MARKER NUMBER",
        ),
        hline("GEODETIC", "MARKER TYPE"),
        hline(
            &format!("{:<20.20}{:<40.40}", opt.observer, opt.agency),
            "OBSERVER / AGENCY",
        ),
        hline(
            &format!(
                "{:<20.20}{:<20.20}{:<20.20}",
                hi.serial, hi.model, hi.firmware
            ),
            "REC # / TYPE / VERS",
        ),
        hline(
            &format!(
                "{:<20}{:<16.16}{:<4.4}",
                "", hi.antenna_type, hi.antenna_dome
            ),
            "ANT # / TYPE",
        ),
    ];
    let (x, y, z) = hi.approx_xyz();
    l.push(hline(
        &format!("{x:14.4}{y:14.4}{z:14.4}"),
        "APPROX POSITION XYZ",
    ));
    let dh = opt.antenna_height.unwrap_or(hi.antenna_height_arp);
    l.push(hline(
        &format!("{dh:14.4}{:14.4}{:14.4}", 0.0, 0.0),
        "ANTENNA: DELTA H/E/N",
    ));
    for (s, ty) in &p.obs_types {
        for (k, chunk) in ty.chunks(13).enumerate() {
            let body: String = chunk.iter().map(|t| format!(" {t}")).collect();
            let lead = if k == 0 {
                format!("{s}  {:3}", ty.len())
            } else {
                " ".repeat(6)
            };
            l.push(hline(&(lead + &body), "SYS / # / OBS TYPES"));
        }
    }
    if let Some(iv) = hi.interval_s.filter(|v| *v != 0.0) {
        l.push(hline(&format!("{iv:10.3}"), "INTERVAL"));
    }
    for (t, lab) in [(&first, "TIME OF FIRST OBS"), (&last, "TIME OF LAST OBS")] {
        l.push(hline(
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
    l.push(hline(&format!("{:6}", 0), "RCV CLOCK OFFS APPL"));
    if !p.glo_k.is_empty() {
        let items: Vec<(&String, &i64)> = p.glo_k.iter().collect();
        for (k, chunk) in items.chunks(8).enumerate() {
            let body: String = chunk.iter().map(|(s, kk)| format!(" {s} {kk:2}")).collect();
            let lead = if k == 0 {
                format!("{:3}", items.len())
            } else {
                "   ".into()
            };
            l.push(hline(&(lead + &body), "GLONASS SLOT / FRQ #"));
        }
        l.push(hline("", "GLONASS COD/PHS/BIS"));
    }
    for (s, code, v) in [('R', "L2C", -0.25), ('G', "L2X", -0.25), ('J', "L2X", 0.25)] {
        if p.obs_types
            .iter()
            .any(|(sys, ty)| *sys == s && ty.iter().any(|t| t == code))
        {
            l.push(hline(&format!("{s} {code} {v:+8.5}"), "SYS / PHASE SHIFT"));
        }
    }
    l.push(hline(
        &format!("{:6}", leap_seconds(&first)),
        "LEAP SECONDS",
    ));
    l.push(hline(&format!("{:6}", p.counts.len()), "# OF SATELLITES"));
    for (sat, c) in &p.counts {
        let ty = p
            .obs_types
            .iter()
            .find(|o| o.0 == sat.sys())
            .map_or(&[][..], |o| o.1.as_slice());
        let vals: Vec<i64> = ty.iter().map(|t| c.get(t).copied().unwrap_or(0)).collect();
        for (k, chunk) in vals.chunks(9).enumerate() {
            let lead = if k == 0 {
                format!("   {sat}")
            } else {
                " ".repeat(6)
            };
            let body: String = chunk.iter().map(|v| format!("{v:6}")).collect();
            l.push(hline(&(lead + &body), "PRN / # OF OBS"));
        }
    }
    l.push(hline("DBHZ", "SIGNAL STRENGTH UNIT")); // written with or without SNR
    l.push(hline("", "END OF HEADER"));
    let mut head = l.join(&opt.line_end);
    head.push_str(&opt.line_end);
    (head, first)
}

/// One system's column layout for the body.
struct SysCols {
    sys: char,
    /// (type letter, code) per column.
    cols: Vec<(char, String)>,
    /// Last column index of a GPS ...2X type, if any.
    two_x: Option<usize>,
}

/// Writes the observation file text into `out`; returns the header text and the first
/// epoch's time.
pub fn write_obs(
    inner: &Inner,
    hi: &HeaderInfo,
    raw_name: &str,
    opt: &ObsOptions,
    out: &mut Vec<u8>,
) -> Result<(String, DateTime), T04Error> {
    // pass 1, over the raw epoch records: what will be written, so header types/counts are exact; the GLONASS
    // frequency numbers as the epochs list them (each satellite's first); the first and last epoch
    let p = pass1(inner, hi, opt)?;
    let (head, first) = header_text(raw_name, hi, opt, &p);
    out.reserve(p.n_epochs * 2500);
    out.extend_from_slice(&crate::util::latin1_bytes(&head));
    // pass 2: the epochs, decoded again and written as they come
    let sys_cols: Vec<SysCols> = p
        .obs_types
        .iter()
        .map(|(s, ty)| SysCols {
            sys: *s,
            cols: ty
                .iter()
                .map(|t| (t.chars().next().unwrap_or(' '), t[1..].to_owned()))
                .collect(),
            two_x: ty.iter().rposition(|t| t.ends_with("2X")),
        })
        .collect();
    let cols_of = |s: char| sys_cols.iter().find(|c| c.sys == s);
    let le = opt.line_end.as_str();
    let mut prev_phase: HashSet<(SatId, &'static str)> = HashSet::new();
    let mut cur_phase: HashSet<(SatId, &'static str)> = HashSet::new();
    let mut pass_cols: HashMap<SatId, usize> = HashMap::new();
    let mut seen: HashSet<SatId> = HashSet::new();
    let mut sel: Vec<(SatId, Kept<usize>)> = Vec::new();
    let mut line = String::with_capacity(512);
    let mut fields: Vec<(usize, usize)> = Vec::new(); // (start, end) of each field in `line`
    measurements::for_each_epoch(inner, opt.continuous_clock, |ep| {
        let obs: &[Obs] = &ep.obs;
        select(
            obs.iter()
                .enumerate()
                .filter(|(_, o)| o.c.is_some())
                .filter_map(|(i, o)| o.code.map(|c| (o.sat, c, o.flags(), i))),
            opt.trimble_filter,
            &mut sel,
        );
        sel.sort_by_key(|a| a.0);
        let t = gps_week_ms(ep.week, ep.tow_ms);
        line.clear();
        let _ = write!(
            line,
            "> {:4} {:02} {:02} {:02} {:02}{:11.7}  0{:3}",
            t.year,
            t.month,
            t.day,
            t.hour,
            t.minute,
            sec_f(&t),
            sel.len()
        );
        if opt.clock_offsets || opt.continuous_clock {
            let ns = round_ratio(i128::from(ep.clk_raw), 524_288 * 1000, 9);
            let _ = write!(line, "{:6}{:>15}", "", format!("{ns}000"));
        }
        line.push_str(le);
        out.extend_from_slice(line.as_bytes());
        // LLI "not continuous": the signal's own phase was not valid at the previous epoch. Not "not written": a
        // satellite whose first signal is invalid is left out whole, and none of its others is then flagged
        cur_phase.clear();
        seen.clear();
        for o in obs {
            seen.insert(o.sat);
            let fl = o.flags();
            if o.c.is_some()
                && (!opt.trimble_filter || ((fl & 0x401) == 0x401 && fl & 0x8 != 0))
                && let Some(c) = o.code
            {
                cur_phase.insert((o.sat, c));
            }
            if (fl & 0x401) != 0x401 {
                continue;
            }
            let Some(sc) = cols_of(o.sat.sys()) else {
                continue;
            };
            if !PAD_SYSTEMS.contains(o.sat.sys()) {
                continue;
            }
            let Some(code) = o.code else { continue };
            let cc = if o.sat.sys() == 'G' && code == "2W" && sc.two_x.is_some() {
                "2X"
            } else {
                code
            };
            if let Some(last) = sc.cols.iter().rposition(|(_, c)| c == cc) {
                let e = pass_cols.entry(o.sat).or_insert(0);
                *e = (*e).max(last + 1);
            }
        }
        pass_cols.retain(|k, _| seen.contains(k));
        for (sat, codes) in &sel {
            let Some(sc) = cols_of(sat.sys()) else {
                continue;
            };
            line.clear();
            let _ = write!(line, "{sat}");
            fields.clear();
            for (t_, c) in &sc.cols {
                let start = line.len();
                match lookup(codes, c) {
                    None => line.push_str("                "),
                    Some(&(code, idx, ph_ok)) => {
                        let o = &obs[idx];
                        let q = Some(ssi(o.s()));
                        match t_ {
                            'C' => push_val(&mut line, &round_float_3(o.c.unwrap_or(0.0)), 0, q),
                            'L' => {
                                if ph_ok {
                                    let lli = i64::from(
                                        o.flags() & 0x100 != 0
                                            || !prev_phase.contains(&(*sat, code)),
                                    );
                                    // L = -phase / 32768 (+ 1/4 for QZSS L2X)
                                    let mut num = -i128::from(o.phase);
                                    if sat.sys() == 'J' && c == "2X" {
                                        num += 8192;
                                    }
                                    push_val(&mut line, &round_ratio(num, 32768, 3), lli, q);
                                } else {
                                    line.push_str("                ");
                                }
                            }
                            'D' => {
                                let d = o.doppler.ok_or_else(|| {
                                    T04Error::Decode("KeyError: 'sv__meas__doppler'".into())
                                })?;
                                push_val(&mut line, &round_ratio(i128::from(d), 256, 3), 0, q);
                            }
                            'S' => {
                                let _ = write!(line, "{:14.3}  ", o.s());
                            }
                            _ => line.push_str("                "),
                        }
                    }
                }
                fields.push((start, line.len()));
            }
            // each line is padded out to the furthest column any code-valid signal of this satellite reached
            // during its current pass (reset when the satellite vanishes from the receiver data); GPS 2W occupies
            // the 2X slot. GPS lines always reach the 2X column (R12i G13).
            let mut keep_to = pass_cols.get(sat).copied().unwrap_or(0);
            if sat.sys() == 'G'
                && let Some(i) = sc.two_x
            {
                keep_to = keep_to.max(i + 1);
            }
            let mut n = fields.len();
            while n > keep_to && line[fields[n - 1].0..fields[n - 1].1].trim().is_empty() {
                n -= 1;
            }
            let end = if n == 0 {
                fields.first().map_or(line.len(), |f| f.0)
            } else {
                fields[n - 1].1
            };
            out.extend_from_slice(&crate::util::latin1_bytes(&line[..end]));
            out.extend_from_slice(le.as_bytes());
        }
        std::mem::swap(&mut prev_phase, &mut cur_phase);
        Ok(())
    })?;
    Ok((head, first))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ssi_steps() {
        assert_eq!(ssi(5.9), 1);
        assert_eq!(ssi(-3.0), 1);
        assert_eq!(ssi(42.0), 7);
        assert_eq!(ssi(70.0), 9);
    }

    #[test]
    fn select_rules() {
        let g1 = SatId::new('G', 1);
        let g2 = SatId::new('G', 2);
        let blocks = vec![
            (g1, "1C", 0x409, 0usize),
            (g1, "2W", 0x401, 1),
            (g2, "1C", 0x001, 2), // first block invalid: whole satellite skipped
            (g2, "2W", 0x409, 3),
        ];
        let mut out = Vec::new();
        select(blocks.into_iter(), true, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1, vec![("1C", 0, true), ("2W", 1, false)]);
    }

    #[test]
    fn value_field() {
        let mut l = String::new();
        push_val(&mut l, "123.456", 1, Some(7));
        assert_eq!(l, "       123.45617");
        let mut l = String::new();
        push_val(&mut l, "123456789012.456", 0, Some(7));
        assert_eq!(l, " ".repeat(16));
    }
}
