//! Writing a RINEX observation file: any supported version, a time window, a new interval, an
//! edited header.
//!
//! The body is built first (in memory) while counting observations, so the header carries exact TIME OF FIRST/LAST OBS, # OF SATELLITES and
//! PRN / # OF OBS. When epochs are dropped (trimming, a longer interval), a loss-of-lock flag
//! (LLI bit 0) seen on a dropped epoch is carried to the satellite's next written epoch, and a
//! power-failure flag to the next written epoch.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use geomaticist_gnss_core::system::SYSTEM_ORDER;
use geomaticist_gnss_core::text::{
    format_g, py_float, py_isspace, py_rstrip, py_strip, sl, to_latin1,
};
use geomaticist_gnss_core::time::{
    DateTime, TICKS, fmt, header_time, v2_epoch_time, v3_epoch_time,
};
use geomaticist_gnss_core::{ObsCode, Sat};

use crate::RinexError;
use crate::header::{
    HeaderEdits, HeaderLine, RinexHeader, f14, obs_types_lines_v2, obs_types_lines_v3,
    prn_obs_lines,
};
use crate::obsmap::{Plan, make_plan, system_name};
use crate::reader::{BLANK, Field, RinexObs};
use crate::stats::{check_interval, effective_interval, select_epochs};

/// Program name written into PGM / RUN BY / DATE.
pub const PROGRAM: &str = "geomaticist-gnss";

/// Header records that exist only in RINEX 3 and later.
pub const V3_ONLY: [&str; 15] = [
    "MARKER TYPE",
    "SYS / # / OBS TYPES",
    "SYS / PHASE SHIFT",
    "GLONASS SLOT / FRQ #",
    "GLONASS COD/PHS/BIS",
    "SYS / DCBS APPLIED",
    "SYS / PCVS APPLIED",
    "SYS / SCALE FACTOR",
    "SIGNAL STRENGTH UNIT",
    "ANTENNA: DELTA X/Y/Z",
    "ANTENNA: PHASECENTER",
    "ANTENNA: B.SIGHT XYZ",
    "ANTENNA: ZERODIR AZI",
    "ANTENNA: ZERODIR XYZ",
    "CENTER OF MASS: XYZ",
];
/// RINEX 4 only.
pub const V4_ONLY: [&str; 3] = ["DOI", "LICENSE OF USE", "STATION INFORMATION"];
/// RINEX 2 only.
pub const V2_ONLY: [&str; 2] = ["# / TYPES OF OBSERV", "WAVELENGTH FACT L1/2"];
/// Records the writer builds itself.
pub const REBUILT: [&str; 25] = [
    "RINEX VERSION / TYPE",
    "PGM / RUN BY / DATE",
    "COMMENT",
    "MARKER NAME",
    "MARKER NUMBER",
    "MARKER TYPE",
    "OBSERVER / AGENCY",
    "REC # / TYPE / VERS",
    "ANT # / TYPE",
    "APPROX POSITION XYZ",
    "ANTENNA: DELTA H/E/N",
    "SYS / # / OBS TYPES",
    "# / TYPES OF OBSERV",
    "INTERVAL",
    "TIME OF FIRST OBS",
    "TIME OF LAST OBS",
    "# OF SATELLITES",
    "PRN / # OF OBS",
    "END OF HEADER",
    "WAVELENGTH FACT L1/2",
    "SYS / PHASE SHIFT",
    "GLONASS SLOT / FRQ #",
    "GLONASS COD/PHS/BIS",
    "LEAP SECONDS",
    "RCV CLOCK OFFS APPL",
];
/// RINEX 3/4 antenna records, copied from the source just before the observation types.
pub const ANTENNA_GROUP: [&str; 6] = [
    "ANTENNA: DELTA X/Y/Z",
    "ANTENNA: PHASECENTER",
    "ANTENNA: B.SIGHT XYZ",
    "ANTENNA: ZERODIR AZI",
    "ANTENNA: ZERODIR XYZ",
    "CENTER OF MASS: XYZ",
];
/// RINEX 3/4 records copied from the source just after the observation types.
pub const AFTER_TYPES: [&str; 4] = [
    "SIGNAL STRENGTH UNIT",
    "SYS / DCBS APPLIED",
    "SYS / PCVS APPLIED",
    "SYS / SCALE FACTOR",
];

/// What to write.
#[derive(Debug, Clone, PartialEq)]
pub struct WriteOptions {
    /// Target version.
    pub version: String,
    /// First epoch (ticks, inclusive).
    pub start: Option<i64>,
    /// Last epoch (ticks, inclusive).
    pub end: Option<i64>,
    /// New interval (ticks); `None` keeps every epoch.
    pub interval: Option<i64>,
    /// Observation kinds to keep (C code, L phase, D Doppler, S SNR).
    pub kinds: String,
    /// Systems to keep.
    pub systems: String,
    /// Header corrections.
    pub edits: HeaderEdits,
    /// Line ending; `None` = the source's.
    pub eol: Option<String>,
    /// Note what was done in COMMENT lines.
    pub add_comments: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            version: "3.04".into(),
            start: None,
            end: None,
            interval: None,
            kinds: "CL".into(),
            systems: SYSTEM_ORDER.into(),
            edits: HeaderEdits::default(),
            eol: None,
            add_comments: true,
        }
    }
}

/// What was written.
#[derive(Debug, Clone, PartialEq)]
pub struct WriteResult {
    /// Output file (empty for in-memory writes).
    pub path: PathBuf,
    /// Epochs written.
    pub epochs: usize,
    /// First epoch written (ticks).
    pub first: Option<i64>,
    /// Last epoch written (ticks).
    pub last: Option<i64>,
    /// Satellites written.
    pub satellites: usize,
    /// Per satellite: values written per output code.
    pub counts: BTreeMap<Sat, Vec<usize>>,
    /// Remarks for the person.
    pub notes: Vec<String>,
    /// The header written, lines right-stripped.
    pub header_lines: Vec<Vec<u8>>,
    /// The plan used.
    pub plan: Plan,
}

fn fmt3(v: [f64; 3]) -> String {
    v.iter().map(|&x| f14(x)).collect()
}

/// `s[:n]` by characters, then latin-1 bytes left-justified to `n`.
fn cut(s: &str, n: usize) -> Vec<u8> {
    let mut b = to_latin1(&s.chars().take(n).collect::<String>());
    b.resize(b.len().max(n), b' ');
    b
}

fn edited(e: Option<&String>, from_header: String) -> String {
    e.cloned().unwrap_or(from_header)
}

/// The source header with the edited fields replaced (same version).
pub fn header_with_edits(hdr: &RinexHeader, e: &HeaderEdits) -> Result<RinexHeader, RinexError> {
    let mut lines = hdr.lines.clone();
    let mut put = |label: &str, content: Vec<u8>| {
        if let Some(l) = lines.iter_mut().find(|l| l.is(label)) {
            l.content = content;
        } else {
            let at = lines.len().saturating_sub(1);
            lines.insert(at, HeaderLine::new(content, label));
        }
    };
    let trim = |s: &str, n: usize| to_latin1(&s.chars().take(n).collect::<String>());
    if let Some(v) = &e.marker_name {
        put("MARKER NAME", trim(v, 60));
    }
    if let Some(v) = &e.marker_number {
        put("MARKER NUMBER", trim(v, 20));
    }
    if let Some(v) = &e.marker_type
        && hdr.major() >= 3
    {
        put("MARKER TYPE", trim(v, 20));
    }
    if e.observer.is_some() || e.agency.is_some() {
        let mut c = cut(&edited(e.observer.as_ref(), hdr.observer()), 20);
        c.extend(cut(&edited(e.agency.as_ref(), hdr.agency()), 40));
        put("OBSERVER / AGENCY", c);
    }
    if e.rec_number.is_some() || e.rec_type.is_some() || e.rec_version.is_some() {
        let mut c = cut(&edited(e.rec_number.as_ref(), hdr.rec_number()), 20);
        c.extend(cut(&edited(e.rec_type.as_ref(), hdr.rec_type()), 20));
        c.extend(cut(&edited(e.rec_version.as_ref(), hdr.rec_version()), 20));
        put("REC # / TYPE / VERS", c);
    }
    if e.ant_number.is_some() || e.ant_type.is_some() || e.ant_radome.is_some() {
        let mut c = cut(&edited(e.ant_number.as_ref(), hdr.ant_number()), 20);
        c.extend(cut(&edited(e.ant_type.as_ref(), hdr.ant_type()), 16));
        c.extend(cut(&edited(e.ant_radome.as_ref(), hdr.ant_radome()), 4));
        put("ANT # / TYPE", c);
    }
    if let Some(v) = e.approx_xyz {
        put("APPROX POSITION XYZ", fmt3(v).into_bytes());
    }
    if let Some(v) = e.delta_hen {
        put("ANTENNA: DELTA H/E/N", fmt3(v).into_bytes());
    }
    RinexHeader::new(lines)
}

fn hl(content: impl AsRef<[u8]>, label: &str) -> HeaderLine {
    HeaderLine::new(content, label)
}

fn label_in(l: &HeaderLine, set: &[&str]) -> bool {
    let lab = py_strip(&l.label);
    set.iter().any(|s| s.as_bytes() == lab)
}

/// The output header.
#[allow(clippy::too_many_arguments)]
pub fn build_header(
    src: &RinexHeader,
    plan: &Plan,
    opts: &WriteOptions,
    first: Option<i64>,
    last: Option<i64>,
    counts: &BTreeMap<Sat, Vec<usize>>,
    interval_s: Option<f64>,
    source_name: &str,
    now: &DateTime,
    decimated: bool,
) -> Result<Vec<HeaderLine>, RinexError> {
    let tgt = &plan.target;
    let major = plan.major();
    let hdr = header_with_edits(src, &opts.edits)?;
    let systems: Vec<u8> = SYSTEM_ORDER
        .bytes()
        .filter(|s| plan.out_types.contains_key(s))
        .collect();
    let mut out: Vec<HeaderLine> = Vec::new();
    let sysc = if systems.len() == 1 { systems[0] } else { b'M' };
    let desc = if major >= 3 {
        if sysc == b'M' {
            "M: Mixed".to_owned()
        } else {
            format!("{}: {}", char::from(sysc), system_name(sysc))
        }
    } else {
        match sysc {
            b'G' => "G (GPS)",
            b'R' => "R (GLONASS)",
            b'E' => "E (GALILEO)",
            b'S' => "S (GEO)",
            _ => "M (MIXED)",
        }
        .to_owned()
    };
    out.push(hl(
        format!("{tgt:>9}{:11}{:<20}{desc:<20}", "", "OBSERVATION DATA"),
        "RINEX VERSION / TYPE",
    ));
    let run_by = opts.edits.run_by.clone().unwrap_or_else(|| src.run_by());
    let mut pgm = cut(PROGRAM, 20);
    pgm.extend(cut(&run_by, 20));
    pgm.extend(
        format!(
            "{:04}{:02}{:02} {:02}{:02}{:02} UTC",
            now.year, now.month, now.day, now.hour, now.minute, now.second
        )
        .bytes(),
    );
    out.push(hl(pgm, "PGM / RUN BY / DATE"));
    if major >= 4 {
        out.extend(src.lines.iter().filter(|l| label_in(l, &V4_ONLY)).cloned());
    }
    for c in src.comments() {
        out.push(hl(sl(&c, 0, 60), "COMMENT"));
    }
    if opts.add_comments {
        let c60 = |s: String| s.chars().take(60).collect::<String>();
        let mut done = vec![c60(format!("Converted by {PROGRAM} from {source_name}"))];
        let program = src.program();
        if !program.is_empty() {
            done.push(c60(format!(
                "Source created by {program} (RINEX {})",
                src.version_str
            )));
        }
        if opts.start.is_some() || opts.end.is_some() {
            // option, else first epoch, else 0: a zero tick counts as unset
            let pick = |a: Option<i64>, b: Option<i64>| {
                a.filter(|&x| x != 0).or(b.filter(|&x| x != 0)).unwrap_or(0)
            };
            done.push(c60(format!(
                "Time window {} - {}",
                fmt(pick(opts.start, first), 1),
                fmt(pick(opts.end, last), 1)
            )));
        }
        if decimated && let Some(iv) = interval_s.filter(|&x| x != 0.0) {
            done.push(c60(format!(
                "Decimated to {} s (nearest epoch, no averaging)",
                format_g(iv)
            )));
        }
        let changed: Vec<&str> = opts
            .edits
            .active()
            .into_iter()
            .filter(|k| *k != "comments" && *k != "run_by")
            .collect();
        if !changed.is_empty() {
            done.push(c60(format!("Header edited: {}", changed.join(", "))));
        }
        out.extend(done.iter().map(|c| hl(to_latin1(c), "COMMENT")));
    }
    for c in &opts.edits.comments {
        out.push(hl(
            to_latin1(&c.chars().take(60).collect::<String>()),
            "COMMENT",
        ));
    }

    let copy = |label: &str| -> Vec<HeaderLine> {
        hdr.lines.iter().filter(|l| l.is(label)).cloned().collect()
    };
    let or = |v: Vec<HeaderLine>, d: HeaderLine| if v.is_empty() { vec![d] } else { v };
    out.extend(or(copy("MARKER NAME"), hl("", "MARKER NAME")));
    out.extend(copy("MARKER NUMBER"));
    if major >= 3 {
        out.extend(copy("MARKER TYPE"));
    }
    out.extend(or(copy("OBSERVER / AGENCY"), hl("", "OBSERVER / AGENCY")));
    out.extend(or(
        copy("REC # / TYPE / VERS"),
        hl("", "REC # / TYPE / VERS"),
    ));
    out.extend(or(copy("ANT # / TYPE"), hl("", "ANT # / TYPE")));
    out.extend(or(
        copy("APPROX POSITION XYZ"),
        hl(fmt3([0.0; 3]), "APPROX POSITION XYZ"),
    ));
    out.extend(or(
        copy("ANTENNA: DELTA H/E/N"),
        hl(fmt3([0.0; 3]), "ANTENNA: DELTA H/E/N"),
    ));
    if major >= 3 {
        for lab in ANTENNA_GROUP {
            out.extend(copy(lab));
        }
        out.extend(obs_types_lines_v3(&plan.out_types));
        for lab in AFTER_TYPES {
            out.extend(copy(lab));
        }
    } else {
        let wl = if src.major() < 3 {
            copy("WAVELENGTH FACT L1/2")
        } else {
            Vec::new()
        };
        // RINEX 3 phase is always full-cycle
        out.extend(or(
            wl,
            hl(format!("{:6}{:6}", 1, 1), "WAVELENGTH FACT L1/2"),
        ));
        out.extend(obs_types_lines_v2(&plan.v2_types));
    }
    if let Some(iv) = interval_s.filter(|&x| x != 0.0) {
        out.push(hl(format!("{iv:10.3}"), "INTERVAL"));
    }
    let ts = src.time_system();
    if let Some(f) = first {
        out.push(hl(header_time(f, &ts), "TIME OF FIRST OBS"));
    }
    if let Some(l) = last {
        out.push(hl(header_time(l, &ts), "TIME OF LAST OBS"));
    }
    out.extend(copy("RCV CLOCK OFFS APPL"));
    if major == 3 {
        // unknown phase shifts (e.g. from RINEX 2) are left out, not guessed
        let mut keep_line = false;
        for l in copy("SYS / PHASE SHIFT") {
            let c = &l.content;
            if !py_strip(sl(c, 0, 1)).is_empty() {
                // a new system / code line; a blank column 1 continues the previous one
                let code = py_strip(sl(c, 2, 5));
                keep_line = plan.out_types.get(&c[0]).is_some_and(|kept| {
                    code.is_empty() || kept.iter().any(|k| k.as_bytes() == code)
                });
            }
            if keep_line {
                out.push(l);
            }
        }
    }
    if major >= 3 && plan.out_types.contains_key(&b'R') {
        out.extend(copy("GLONASS SLOT / FRQ #"));
        if major == 3 {
            out.extend(or(
                copy("GLONASS COD/PHS/BIS"),
                hl("", "GLONASS COD/PHS/BIS"),
            ));
        }
    }
    if let Some(leap) = src.leap_seconds() {
        out.push(hl(
            if major >= 3 { leap } else { sl(leap, 0, 6) },
            "LEAP SECONDS",
        ));
    }
    for l in &hdr.lines {
        // anything else the source had, if this version knows it
        if label_in(l, &REBUILT)
            || label_in(l, &V4_ONLY)
            || (label_in(l, &V3_ONLY) && major < 3)
            || (label_in(l, &V2_ONLY) && major >= 3)
            || label_in(l, &ANTENNA_GROUP)
            || label_in(l, &AFTER_TYPES)
        {
            continue;
        }
        out.push(l.clone());
    }
    out.push(hl(format!("{:6}", counts.len()), "# OF SATELLITES"));
    out.extend(prn_obs_lines(counts, major));
    out.push(hl("", "END OF HEADER"));
    Ok(out)
}

/// A RINEX 3/4 epoch line.
#[must_use]
pub fn epoch_line_v3(tick: Option<i64>, flag: u8, n: usize, clock: &[u8]) -> Vec<u8> {
    let t = tick.map_or_else(|| " ".repeat(27), v3_epoch_time);
    let mut line = format!("> {t}  {flag}{n:3}").into_bytes();
    let c = py_strip(clock);
    if !c.is_empty() {
        line.extend_from_slice(&[b' '; 6]);
        line.extend(std::iter::repeat_n(b' ', 15usize.saturating_sub(c.len())));
        line.extend_from_slice(c);
    }
    line
}

/// The receiver clock offset reformatted for the other family (F15.12 / F12.9); empty when
/// blank or unreadable.
#[must_use]
pub fn clock_to(clock: &[u8], major: i64) -> Vec<u8> {
    if py_strip(clock).is_empty() {
        return Vec::new();
    }
    let Some(v) = py_float(clock) else {
        return Vec::new();
    };
    if major >= 3 {
        format!("{v:15.12}").into_bytes()
    } else {
        format!("{v:12.9}").into_bytes()
    }
}

/// RINEX 2 epoch line(s): 12 satellites per line, the clock after the first twelve.
#[must_use]
pub fn epoch_lines_v2(
    tick: Option<i64>,
    flag: u8,
    sats: &[Sat],
    n: usize,
    clock: &[u8],
) -> Vec<Vec<u8>> {
    let t = tick.map_or_else(|| " ".repeat(26), v2_epoch_time);
    let head = format!("{t}  {flag}{n:3}").into_bytes();
    if sats.is_empty() {
        return vec![head];
    }
    let ids: Vec<String> = sats
        .iter()
        .map(|s| format!("{}{:2}", s.sys.letter(), s.prn))
        .collect();
    let mut lines = Vec::new();
    for (k, chunk) in ids.chunks(12).enumerate() {
        let chunk = chunk.concat();
        if k == 0 {
            let mut line = head.clone();
            line.extend(format!("{chunk:<36}").bytes());
            if !clock.is_empty() {
                line.extend(std::iter::repeat_n(
                    b' ',
                    12usize.saturating_sub(clock.len()),
                ));
                line.extend_from_slice(clock);
            }
            lines.push(py_rstrip(&line).to_vec());
        } else {
            let line = format!("{}{chunk}", " ".repeat(32));
            lines.push(py_rstrip(line.as_bytes()).to_vec());
        }
    }
    lines
}

fn obs_types_of(hdr: &RinexHeader) -> BTreeMap<u8, Vec<ObsCode>> {
    if hdr.major() >= 3 {
        hdr.obs_types.clone()
    } else {
        hdr.obs_types
            .keys()
            .map(|s| (*s, hdr.v2_types.clone()))
            .collect()
    }
}

fn sat_codes(obs: &RinexObs) -> BTreeMap<Sat, BTreeSet<ObsCode>> {
    obs.sats
        .iter()
        .map(|(sat, st)| {
            (
                *sat,
                st.counts
                    .iter()
                    .filter(|(_, n)| **n > 0)
                    .map(|(c, _)| *c)
                    .collect(),
            )
        })
        .collect()
}

/// The plan for writing `obs` with `opts`.
#[must_use]
pub fn plan_for(obs: &RinexObs, opts: &WriteOptions) -> Plan {
    let src = &obs.header;
    make_plan(
        &obs_types_of(src),
        src.version,
        &opts.version,
        &opts.kinds,
        Some(&opts.systems),
        Some(&sat_codes(obs)),
    )
}

/// After a mid-file change of observation types: points every output code at its field in the new
/// layout (by code name). Returns the new codes with no place in the output.
fn remap(
    plan: &mut Plan,
    obs: &RinexObs,
    opts: &WriteOptions,
    hdr: &RinexHeader,
) -> BTreeSet<String> {
    let new = make_plan(
        &obs_types_of(hdr),
        hdr.version,
        &opts.version,
        &opts.kinds,
        Some(&opts.systems),
        Some(&sat_codes(obs)),
    );
    let mut lost = BTreeSet::new();
    let systems: Vec<u8> = plan.field_map.keys().copied().collect();
    for s in systems {
        let old_codes = plan.out_types.get(&s).cloned().unwrap_or_default();
        let empty = Vec::new();
        let new_codes = new.out_types.get(&s).unwrap_or(&empty);
        let nmap = new.field_map.get(&s).cloned().unwrap_or_default();
        let fm = old_codes
            .iter()
            .map(|c| {
                new_codes
                    .iter()
                    .position(|x| x == c)
                    .and_then(|i| nmap.get(i).copied().flatten())
            })
            .collect();
        plan.field_map.insert(s, fm);
        for c in new_codes {
            if !old_codes.contains(c) {
                lost.insert(format!("{} {c}", char::from(s)));
            }
        }
    }
    let sats: Vec<Sat> = plan.sat_map.keys().copied().collect();
    for sat in sats {
        if let Some(nm) = new.sat_map.get(&sat) {
            plan.sat_map.insert(sat, nm.clone()); // RINEX 2 output: same codes, new source fields
        } else {
            plan.sat_map.remove(&sat);
        }
    }
    lost
}

/// The value of an LLI character: a digit, else an "invalid literal" error.
fn lli_digit(c: u8) -> Result<u8, RinexError> {
    if c.is_ascii_digit() {
        Ok(c - b'0')
    } else {
        Err(RinexError::new(format!(
            "invalid literal for int() with base 10: '{}'",
            char::from(c)
        )))
    }
}

fn field_map(plan: &Plan, sat: Sat) -> Option<&[Option<usize>]> {
    // `plan.sat_map.get(sat) or plan.field_map.get(sat[0])`
    plan.sat_map
        .get(&sat)
        .filter(|v| !v.is_empty())
        .or_else(|| plan.field_map.get(&(sat.sys.letter() as u8)))
        .map(Vec::as_slice)
        .filter(|v| !v.is_empty())
}

/// The current UTC time (for PGM / RUN BY / DATE).
#[must_use]
pub fn utc_now() -> DateTime {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let base = geomaticist_gnss_core::time::ordinal(1970, 1, 1).unwrap_or(719_163);
    let (y, m, d) = geomaticist_gnss_core::time::from_ordinal(base + days);
    DateTime {
        year: y,
        month: m,
        day: d,
        hour: rem / 3600,
        minute: rem % 3600 / 60,
        second: rem % 60,
        micro: 0,
    }
}

/// Writes `obs` with `opts` and returns the file's bytes.
#[allow(clippy::too_many_lines, clippy::cast_precision_loss)]
pub fn write_obs_bytes(
    obs: &RinexObs,
    opts: &WriteOptions,
    now: Option<&DateTime>,
) -> Result<(Vec<u8>, WriteResult), RinexError> {
    let src = &obs.header;
    let mut plan = plan_for(obs, opts);
    let mut notes = plan.notes.clone();
    if plan.out_types.is_empty() {
        return Err(RinexError::new(
            "nothing left to write - every system or observation type was filtered out",
        ));
    }
    if obs.out_of_order > 0 {
        notes.push(format!(
            "{} epoch(s) in the source are out of time order - written as found.",
            obs.out_of_order
        ));
    }
    let sample = obs.sample_ticks();
    check_interval(opts.interval, sample)?;
    let eff_iv = effective_interval(opts.interval, sample);
    let keep_idx = select_epochs(&obs.epoch_ticks, opts.start, opts.end, eff_iv, sample);
    if keep_idx.is_empty() {
        return Err(RinexError::new("no epochs in the chosen time window"));
    }
    let mut keep = vec![false; obs.epoch_ticks.len()];
    for &i in &keep_idx {
        keep[i] = true;
    }
    let win_lo = opts.start.unwrap_or(-(1 << 62));
    let win_hi = opts.end.unwrap_or(1 << 62);
    let major = plan.major();
    let same_family = (src.major() >= 3) == (major >= 3);
    let eol: Vec<u8> = opts.eol.as_deref().unwrap_or(&obs.eol).as_bytes().to_vec();
    let mut body: Vec<u8> = Vec::with_capacity(obs.data.len());
    let mut counts: BTreeMap<Sat, Vec<usize>> = BTreeMap::new();
    let mut pending: BTreeMap<Sat, BTreeSet<usize>> = BTreeMap::new();
    let (mut first, mut last): (Option<i64>, Option<i64>) = (None, None);
    let mut written = 0usize;
    let mut dropped_events = 0usize;
    let mut i_obs: i64 = -1;
    let mut lost_codes: BTreeSet<String> = BTreeSet::new();
    let mut flag1_pending = false;
    let mut flag6 = 0usize;
    let mut rd = obs.epochs();
    let mut types_now = rd.header().types_generation;
    let push = |body: &mut Vec<u8>, line: &[u8]| {
        body.extend_from_slice(line);
        body.extend_from_slice(&eol);
    };
    while let Some(ep) = rd.next() {
        let ep = ep?;
        if rd.header().types_generation != types_now {
            // a flag-4 event changed the observation types: map the new layout onto the output
            types_now = rd.header().types_generation;
            let h = rd.header().clone();
            lost_codes.extend(remap(&mut plan, obs, opts, &h));
        }
        if ep.flag == 6 {
            flag6 += 1; // cycle-slip records repeat data; not carried over
            continue;
        }
        if ep.flag > 1 {
            if ep.tick.is_some_and(|t| !(win_lo <= t && t <= win_hi)) {
                continue;
            }
            if ep.flag == 4 {
                dropped_events += 1; // header records: the output header stands for them
                continue;
            }
            let n = ep.event_lines.len();
            let line = if major >= 3 {
                epoch_line_v3(ep.tick, ep.flag, n, b"")
            } else {
                epoch_lines_v2(ep.tick, ep.flag, &[], n, b"").swap_remove(0)
            };
            push(&mut body, &line);
            for l in &ep.event_lines {
                push(&mut body, l);
            }
            continue;
        }
        let Some(tick) = ep.tick else { continue };
        i_obs += 1;
        let in_window = win_lo <= tick && tick <= win_hi;
        if !keep[i_obs as usize] {
            if in_window && ep.flag == 1 {
                flag1_pending = true; // power failure: tell the next written epoch
            }
            if in_window {
                // remember slips on dropped epochs
                for (sat, flds) in ep.records() {
                    let Some(fm) = field_map(&plan, sat) else {
                        continue;
                    };
                    let Some(po) = plan.phase_out.get(&(sat.sys.letter() as u8)) else {
                        continue;
                    };
                    for &j in po {
                        if let Some(Some(si)) = fm.get(j).copied()
                            && let Some(f) = flds.get(si)
                            && !py_isspace(f[14])
                            && lli_digit(f[14])? & 1 == 1
                        {
                            pending.entry(sat).or_default().insert(j);
                        }
                    }
                }
            }
            continue;
        }
        let mut recs: Vec<(Sat, Vec<Field>)> = Vec::new();
        for (sat, flds) in ep.records() {
            let Some(fm) = field_map(&plan, sat) else {
                continue;
            };
            let mut outf: Vec<Field> = fm
                .iter()
                .map(|si| si.and_then(|i| flds.get(i).copied()).unwrap_or(BLANK))
                .collect();
            if outf.iter().all(|f| py_strip(&f[..14]).is_empty()) {
                continue;
            }
            if let Some(p) = pending.get_mut(&sat)
                && !p.is_empty()
            {
                let js: Vec<usize> = p.iter().copied().collect();
                for j in js {
                    let f = &mut outf[j];
                    if !py_strip(&f[..14]).is_empty() {
                        let old = if py_isspace(f[14]) {
                            0
                        } else {
                            lli_digit(f[14])?
                        };
                        f[14] = b'0' + (old | 1);
                        p.remove(&j);
                    }
                }
            }
            let c = counts.entry(sat).or_insert_with(|| vec![0; fm.len()]);
            for (j, f) in outf.iter().enumerate() {
                if !py_strip(&f[..14]).is_empty() {
                    c[j] += 1;
                }
            }
            recs.push((sat, outf));
        }
        if recs.is_empty() {
            continue;
        }
        let mut flag = ep.flag;
        if flag1_pending {
            flag = 1;
            flag1_pending = false;
        }
        first = Some(first.map_or(tick, |f| f.min(tick)));
        last = Some(last.map_or(tick, |l| l.max(tick)));
        written += 1;
        if major >= 3 {
            let clock = if same_family {
                ep.clock.clone()
            } else {
                clock_to(&ep.clock, major)
            };
            push(
                &mut body,
                &epoch_line_v3(Some(tick), flag, recs.len(), &clock),
            );
            for (sat, outf) in &recs {
                let mut line = sat.id_bytes().to_vec();
                for f in outf {
                    line.extend_from_slice(f);
                }
                push(&mut body, py_rstrip(&line));
            }
        } else {
            let clock = if same_family {
                ep.clock.clone()
            } else {
                clock_to(&ep.clock, 2)
            };
            let sats: Vec<Sat> = recs.iter().map(|r| r.0).collect();
            for l in epoch_lines_v2(Some(tick), flag, &sats, recs.len(), py_strip(&clock)) {
                push(&mut body, &l);
            }
            for (_, outf) in &recs {
                for chunk in outf.chunks(5) {
                    let line = chunk.concat();
                    push(&mut body, py_rstrip(&line));
                }
            }
        }
    }
    if dropped_events > 0 {
        notes.push(format!(
            "{dropped_events} header-change record(s) (event flag 4) folded into the output header."
        ));
    }
    if flag6 > 0 {
        notes.push(format!(
            "{flag6} cycle-slip record(s) (event flag 6) left out - they repeat earlier data."
        ));
    }
    if !lost_codes.is_empty() {
        notes.push(format!(
            "Observation types added mid-file ({}) aren't in the output header - left out.",
            lost_codes.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    let interval_s = eff_iv.unwrap_or(sample) as f64 / TICKS as f64;
    let source_name = obs
        .source_path
        .file_name()
        .map_or_default(|n| n.to_string_lossy().into_owned());
    let now = now.copied().unwrap_or_else(utc_now);
    let header = build_header(
        src,
        &plan,
        opts,
        first,
        last,
        &counts,
        Some(interval_s),
        &source_name,
        &now,
        eff_iv.is_some(),
    )?;
    let mut out = Vec::with_capacity(body.len() + 80 * header.len());
    let mut header_lines = Vec::with_capacity(header.len());
    for l in &header {
        let t = l.text();
        let t = py_rstrip(&t);
        out.extend_from_slice(t);
        out.extend_from_slice(&eol);
        header_lines.push(t.to_vec());
    }
    out.extend_from_slice(&body);
    let res = WriteResult {
        path: PathBuf::new(),
        epochs: written,
        first,
        last,
        satellites: counts.len(),
        counts,
        notes,
        header_lines,
        plan,
    };
    Ok((out, res))
}

/// Writes `obs` to `out_path` (through `<name>.part`, then renamed).
pub fn write_obs(
    obs: &RinexObs,
    out_path: &Path,
    opts: &WriteOptions,
    now: Option<&DateTime>,
) -> Result<WriteResult, RinexError> {
    let (bytes, mut res) = write_obs_bytes(obs, opts, now)?;
    if let Some(parent) = out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    let mut part = out_path.as_os_str().to_owned();
    part.push(".part");
    let part = PathBuf::from(part);
    std::fs::write(&part, &bytes)?;
    std::fs::rename(&part, out_path)?;
    res.path = out_path.to_path_buf();
    Ok(res)
}

/// The RINEX normalisation used to compare files: columns 41-60 of every PGM / RUN BY / DATE line become
/// `DATE-NORMALIZED     `; nothing else changes.
#[must_use]
pub fn normalize_rinex(data: &[u8]) -> Vec<u8> {
    const LABEL: &[u8] = b"PGM / RUN BY / DATE";
    const MASK: &[u8] = b"DATE-NORMALIZED     ";
    let mut out = Vec::with_capacity(data.len());
    let mut i = 0;
    while i < data.len() {
        // bytes.splitlines(keepends=True): \n, \r, \r\n
        let mut j = i;
        while j < data.len() && data[j] != b'\n' && data[j] != b'\r' {
            j += 1;
        }
        let mut k = j;
        if k < data.len() {
            k += if data[k] == b'\r' && data.get(k + 1) == Some(&b'\n') {
                2
            } else {
                1
            };
        }
        let line = &data[i..j];
        // body = line.rstrip(b"\r\n"): the terminator only, as the line holds no other \r / \n
        if sl(line, 60, line.len()).starts_with(LABEL) {
            let mut b = sl(line, 0, 40).to_vec();
            b.resize(40, b' ');
            b.extend_from_slice(MASK);
            b.extend_from_slice(sl(line, 60, line.len()));
            out.extend_from_slice(&b);
        } else {
            out.extend_from_slice(line);
        }
        out.extend_from_slice(&data[j..k]);
        i = k;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn epoch_lines() {
        let g = |s: &str| s.parse::<Sat>().unwrap();
        let l = epoch_line_v3(Some(0), 0, 2, b"  -0.000123456789");
        assert_eq!(
            String::from_utf8(l).unwrap(),
            "> 1980 01 06 00 00  0.0000000  0  2      -0.000123456789"
        );
        let sats: Vec<Sat> = (1..=13).map(|i| g(&format!("G{i:02}"))).collect();
        let v2 = epoch_lines_v2(Some(0), 0, &sats, 13, b"-0.000123456");
        assert_eq!(v2.len(), 2);
        assert_eq!(v2[0].len(), 80);
        assert_eq!(
            String::from_utf8_lossy(&v2[1]),
            format!("{}G13", " ".repeat(32))
        );
        assert_eq!(clock_to(b"  -.5", 3), b"-0.500000000000".to_vec());
        assert_eq!(clock_to(b"x", 3), Vec::<u8>::new());
    }

    #[test]
    fn normalises_pgm_line() {
        let l = format!(
            "{:<20}{:<20}{:<20}PGM / RUN BY / DATE\r\nX\n",
            "geomaticist-gnss", "me", "20260101 000000 UTC"
        );
        let n = normalize_rinex(l.as_bytes());
        assert_eq!(
            String::from_utf8(n).unwrap(),
            format!(
                "{:<20}{:<20}DATE-NORMALIZED     PGM / RUN BY / DATE\r\nX\n",
                "geomaticist-gnss", "me"
            )
        );
    }
}
