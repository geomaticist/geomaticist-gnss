//! Reading RINEX observation files (2.10/2.11, 3.0x, 4.0x) as a stream of epochs, and the indexed
//! file the summary, trim and convert tools work on.
//!
//! Observation fields are kept as their original 16-character text (value F14.3 + LLI + signal
//! strength), so a value is never re-rounded on its way through the writer.
//!
//! A satellite id that is not a satellite (`"G-1"`, `"Gab"`, an unknown system letter) is
//! dropped, its observation records still consumed; [`EpochReader::dropped_sats`] counts them.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use geomaticist_gnss_core::system::SYSTEM_ORDER;
use geomaticist_gnss_core::text::{is_blank, lines_lf, py_strip, sl};
use geomaticist_gnss_core::time::{TICKS, parse_seconds, py_int, to_ticks};
use geomaticist_gnss_core::{GnssSystem, ObsCode, Sat};
use serde::Serialize;

use crate::RinexError;
use crate::header::{HeaderLine, RinexHeader, read_header_lines};

/// A 16-character observation field (value F14.3, LLI, signal strength), blank-padded.
pub type Field = [u8; 16];

/// Blank field.
pub const BLANK: Field = [b' '; 16];

/// One epoch record.
#[derive(Debug, Clone, PartialEq)]
pub struct Epoch {
    /// Ticks (file time system), `None` when the epoch time is unreadable.
    pub tick: Option<i64>,
    /// Event flag.
    pub flag: u8,
    /// Satellites, in record order.
    pub sats: Vec<Sat>,
    fields: Vec<Field>,
    offsets: Vec<usize>,
    /// Receiver clock offset text, right-stripped, as written.
    pub clock: Vec<u8>,
    /// Flags 2-6 (RINEX 2: 2-5): the records that follow.
    pub event_lines: Vec<Vec<u8>>,
}

impl Epoch {
    fn new(tick: Option<i64>, flag: u8, clock: &[u8]) -> Self {
        let mut end = clock.len();
        while end > 0 && geomaticist_gnss_core::text::py_isspace(clock[end - 1]) {
            end -= 1;
        }
        Self {
            tick,
            flag,
            sats: Vec::new(),
            fields: Vec::new(),
            offsets: vec![0],
            clock: clock[..end].to_vec(),
            event_lines: Vec::new(),
        }
    }

    fn push_sat(&mut self, sat: Sat) {
        self.sats.push(sat);
        self.offsets.push(self.fields.len());
    }

    /// The fields of satellite `i` (as many as the header had types, or the record held).
    #[must_use]
    pub fn fields(&self, i: usize) -> &[Field] {
        &self.fields[self.offsets[i]..self.offsets[i + 1]]
    }

    /// (satellite, fields) pairs.
    pub fn records(&self) -> impl Iterator<Item = (Sat, &[Field])> {
        self.sats
            .iter()
            .enumerate()
            .map(|(i, s)| (*s, self.fields(i)))
    }
}

/// A three-character satellite id with a default system. `None` when it is not a satellite
/// (`"G-1"`, `"Gab"`, an unknown system letter).
#[must_use]
pub fn parse_sat(s: &[u8], default_sys: u8) -> Option<Sat> {
    let c = py_strip(sl(s, 0, 1));
    let sys = c.first().map_or(default_sys, u8::to_ascii_uppercase);
    let sys = GnssSystem::from_letter(char::from(sys))?;
    let num = sl(s, 1, 3);
    let prn = match py_int(num) {
        Some(n) => n,
        None if py_strip(num).is_empty() => 0, // blank number: PRN 0 ("G00")
        None => return None,
    };
    let prn = u8::try_from(prn).ok().filter(|p| *p <= 99)?;
    Some(Sat::new(sys, prn))
}

fn year2(yy: i64) -> i64 {
    if yy < 80 { 2000 + yy } else { 1900 + yy }
}

fn int_or0(s: &[u8]) -> Result<i64, RinexError> {
    if s.is_empty() {
        return Ok(0);
    }
    py_int(s).ok_or_else(|| {
        RinexError::new(format!(
            "invalid literal for int(): {:?}",
            String::from_utf8_lossy(s)
        ))
    })
}

fn field16(src: &[u8]) -> Field {
    let mut f = BLANK;
    let n = src.len().min(16);
    f[..n].copy_from_slice(&src[..n]);
    f
}

/// Streams the epochs of a file's bytes; the working header it owns
/// is updated by flag-4 header events.
pub struct EpochReader<'a> {
    lines: Box<dyn Iterator<Item = &'a [u8]> + 'a>,
    header: RinexHeader,
    v3: bool,
    default_sys: u8,
    /// Satellite ids dropped (non-satellite strings).
    pub dropped_sats: usize,
    failed: bool,
}

impl<'a> EpochReader<'a> {
    /// A reader over the whole file (header included; it is skipped).
    #[must_use]
    pub fn new(data: &'a [u8], header: RinexHeader) -> Self {
        let mut lines = lines_lf(data);
        for line in lines.by_ref() {
            if py_strip(sl(line, 60, 80)) == b"END OF HEADER" {
                break;
            }
        }
        let v3 = header.major() >= 3;
        let default_sys = if header.system == b'M' {
            b'G'
        } else {
            header.system
        };
        Self {
            lines: Box::new(lines),
            header,
            v3,
            default_sys,
            dropped_sats: 0,
            failed: false,
        }
    }

    /// The working header (observation types as changed by flag-4 events so far).
    #[must_use]
    pub fn header(&self) -> &RinexHeader {
        &self.header
    }

    fn nxt(&mut self) -> Result<&'a [u8], RinexError> {
        self.lines
            .next()
            .ok_or_else(|| RinexError::new("the file ends in the middle of an epoch"))
    }

    fn event_header(&mut self, ep: &Epoch) -> Result<(), RinexError> {
        if ep.flag == 4 {
            let lines: Vec<HeaderLine> = ep
                .event_lines
                .iter()
                .map(|r| HeaderLine::from_raw(r))
                .collect();
            self.header.update_from_event(&lines)?;
        }
        Ok(())
    }

    fn next_v3(&mut self, line: &'a [u8]) -> Result<Option<Epoch>, RinexError> {
        if !line.starts_with(b">") {
            return Ok(None); // stray line; RINEX 3 epochs always start with '>'
        }
        let flag = int_or0(sl(line, 31, 32))?;
        let n = int_or0(sl(line, 32, 35))?;
        let tick = (|| {
            to_ticks(
                py_int(sl(line, 2, 6))?,
                py_int(sl(line, 7, 9))?,
                py_int(sl(line, 10, 12))?,
                py_int(sl(line, 13, 15))?,
                py_int(sl(line, 16, 18))?,
                parse_seconds(sl(line, 18, 29))?,
            )
        })();
        let flag = u8::try_from(flag).unwrap_or(u8::MAX);
        let mut ep = Epoch::new(tick, flag, sl(line, 41, 56));
        if flag <= 1 {
            for _ in 0..n {
                let rec = self.nxt()?;
                let sat = parse_sat(sl(rec, 0, 3), b'G');
                let sys_letter = sat.map_or_else(
                    || {
                        py_strip(sl(rec, 0, 1))
                            .first()
                            .map_or(b'G', u8::to_ascii_uppercase)
                    },
                    |s| s.sys.letter() as u8,
                );
                let mut nt = self.header.types_for(sys_letter).len();
                if nt == 0 {
                    nt = rec.len().saturating_sub(3).div_ceil(16).max(1);
                }
                let Some(sat) = sat else {
                    self.dropped_sats += 1;
                    continue;
                };
                let body = sl(rec, 3, rec.len());
                ep.push_sat(sat);
                for i in 0..nt {
                    ep.fields.push(field16(sl(body, 16 * i, 16 * i + 16)));
                }
                *ep.offsets.last_mut().unwrap_or(&mut 0) = ep.fields.len();
            }
        } else {
            for _ in 0..n {
                let rec = self.nxt()?;
                ep.event_lines.push(rec.to_vec());
            }
            self.event_header(&ep)?;
        }
        Ok(Some(ep))
    }

    fn next_v2(&mut self, line: &'a [u8]) -> Result<Option<Epoch>, RinexError> {
        if line.len() < 29 {
            return Ok(None);
        }
        let flag = int_or0(sl(line, 28, 29))?;
        let n = int_or0(sl(line, 29, 32))?;
        let tick = (|| {
            to_ticks(
                year2(py_int(sl(line, 1, 3))?),
                py_int(sl(line, 4, 6))?,
                py_int(sl(line, 7, 9))?,
                py_int(sl(line, 10, 12))?,
                py_int(sl(line, 13, 15))?,
                parse_seconds(sl(line, 15, 26))?,
            )
        })();
        let flag = u8::try_from(flag).unwrap_or(u8::MAX);
        let mut ep = Epoch::new(tick, flag, sl(line, 68, 80));
        if matches!(flag, 0 | 1 | 6) {
            // flag 6 (cycle-slip records) has the observation layout
            let n = usize::try_from(n).unwrap_or(0);
            let mut ids: Vec<Option<Sat>> = Vec::with_capacity(n);
            let mut cur = line;
            loop {
                for k in 0..12 {
                    if ids.len() >= n {
                        break;
                    }
                    ids.push(parse_sat(sl(cur, 32 + 3 * k, 35 + 3 * k), self.default_sys));
                }
                if ids.len() >= n {
                    break;
                }
                cur = self.nxt()?;
            }
            let nt = self.header.v2_types.len();
            let nlines = nt.div_ceil(5).max(1);
            for sat in ids {
                let start = ep.fields.len();
                let mut got = 0;
                for _ in 0..nlines {
                    let rec = self.nxt()?;
                    for i in 0..5 {
                        if got < nt {
                            ep.fields.push(field16(sl(rec, 16 * i, 16 * i + 16)));
                            got += 1;
                        }
                    }
                }
                if let Some(s) = sat {
                    ep.push_sat(s);
                    *ep.offsets.last_mut().unwrap_or(&mut 0) = ep.fields.len();
                } else {
                    self.dropped_sats += 1;
                    ep.fields.truncate(start);
                }
            }
        } else {
            for _ in 0..n {
                let rec = self.nxt()?;
                ep.event_lines.push(rec.to_vec());
            }
            self.event_header(&ep)?;
        }
        Ok(Some(ep))
    }
}

impl Iterator for EpochReader<'_> {
    type Item = Result<Epoch, RinexError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        loop {
            let line = self.lines.next()?;
            if is_blank(line) {
                continue;
            }
            let r = if self.v3 {
                self.next_v3(line)
            } else {
                self.next_v2(line)
            };
            match r {
                Ok(Some(ep)) => return Some(Ok(ep)),
                Ok(None) => {}
                Err(e) => {
                    self.failed = true;
                    return Some(Err(e));
                }
            }
        }
    }
}

/// Reads the header of a file's bytes and the line ending of its first line.
pub fn read_header(data: &[u8]) -> Result<(RinexHeader, String), RinexError> {
    let first = data.split(|&b| b == b'\n').next().unwrap_or(&[]);
    let eol = if first.ends_with(b"\r") && data.len() > first.len() {
        "\r\n"
    } else {
        "\n"
    };
    Ok((RinexHeader::new(read_header_lines(data)?)?, eol.into()))
}

/// Per-satellite statistics of the indexed file.
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct SatStats {
    /// Epochs with at least one observation.
    pub ticks: Vec<i64>,
    /// Epochs with a loss-of-lock flag on any phase.
    pub slips: Vec<i64>,
    /// Observation code -> number of values.
    pub counts: BTreeMap<ObsCode, usize>,
    /// Phase code -> number of loss-of-lock flags.
    pub lli: BTreeMap<ObsCode, usize>,
    /// Signal strength indicator sums.
    pub ssi_sum: BTreeMap<ObsCode, i64>,
    /// Signal strength indicator counts.
    pub ssi_n: BTreeMap<ObsCode, usize>,
    /// S observation sums (dB-Hz).
    pub snr_sum: BTreeMap<ObsCode, f64>,
    /// S observation counts.
    pub snr_n: BTreeMap<ObsCode, usize>,
}

/// An observation file: its header plus an index of epoch times and per-satellite statistics.
#[derive(Debug, Clone)]
pub struct RinexObs {
    /// The path the user opened.
    pub source_path: PathBuf,
    /// Plain RINEX text (decompressed).
    pub data: Vec<u8>,
    /// The file's header.
    pub header: RinexHeader,
    /// Line ending of the first line.
    pub eol: String,
    /// Observation epochs (flags 0/1), in file order.
    pub epoch_ticks: Vec<i64>,
    /// Satellites per observation epoch.
    pub epoch_nsat: Vec<i64>,
    /// Per-satellite statistics.
    pub sats: BTreeMap<Sat, SatStats>,
    /// (tick, flag, number of records) of event epochs.
    pub events: Vec<(Option<i64>, u8, usize)>,
    /// Per system, the types in force when its first satellite was seen.
    pub obs_types_seen: BTreeMap<u8, Vec<ObsCode>>,
    /// Epochs earlier than the one before.
    pub out_of_order: usize,
    /// Epochs equal to the one before.
    pub duplicates: usize,
    /// Satellite ids dropped (non-satellite strings).
    pub dropped_sats: usize,
}

impl RinexObs {
    /// Opens (gzip / Hatanaka as needed) and indexes a file.
    pub fn open(path: &Path) -> Result<Self, RinexError> {
        let data = crate::open::read_plain(path)?;
        Self::from_bytes(path, data)
    }

    /// Indexes plain RINEX bytes.
    pub fn from_bytes(source: &Path, data: Vec<u8>) -> Result<Self, RinexError> {
        let (header, eol) = read_header(&data)?;
        let mut obs = Self {
            source_path: source.to_path_buf(),
            data,
            header,
            eol,
            epoch_ticks: Vec::new(),
            epoch_nsat: Vec::new(),
            sats: BTreeMap::new(),
            events: Vec::new(),
            obs_types_seen: BTreeMap::new(),
            out_of_order: 0,
            duplicates: 0,
            dropped_sats: 0,
        };
        obs.index()?;
        Ok(obs)
    }

    fn index(&mut self) -> Result<(), RinexError> {
        let data = std::mem::take(&mut self.data);
        let mut rd = EpochReader::new(&data, self.header.clone());
        let mut prev: Option<i64> = None;
        while let Some(ep) = rd.next() {
            let ep = ep?;
            if ep.flag > 1 {
                self.events.push((ep.tick, ep.flag, ep.event_lines.len()));
                continue;
            }
            let Some(tick) = ep.tick else { continue };
            if let Some(p) = prev {
                if tick < p {
                    self.out_of_order += 1;
                } else if tick == p {
                    self.duplicates += 1;
                }
            }
            prev = Some(tick);
            self.epoch_ticks.push(tick);
            self.epoch_nsat.push(ep.sats.len() as i64);
            for (sat, flds) in ep.records() {
                let sys = sat.sys.letter() as u8;
                let types = rd.header().types_for(sys);
                let st = self.sats.entry(sat).or_insert_with(|| {
                    self.obs_types_seen
                        .entry(sys)
                        .or_insert_with(|| types.to_vec());
                    SatStats::default()
                });
                let mut any_obs = false;
                let mut slip = false;
                for (code, f) in types.iter().zip(flds) {
                    let v = &f[..14];
                    if is_blank(v) {
                        continue;
                    }
                    any_obs = true;
                    *st.counts.entry(*code).or_default() += 1;
                    let kind = code.kind_char();
                    let lli = f[14];
                    if kind == 'L' && !is_blank(&[lli]) {
                        let l = py_int(&[lli]).ok_or_else(|| RinexError::new("invalid LLI"))?;
                        if l & 1 != 0 {
                            *st.lli.entry(*code).or_default() += 1;
                            slip = true;
                        }
                    }
                    let ssi = f[15];
                    if !is_blank(&[ssi]) {
                        let s = py_int(&[ssi]).ok_or_else(|| RinexError::new("invalid SSI"))?;
                        *st.ssi_sum.entry(*code).or_default() += s;
                        *st.ssi_n.entry(*code).or_default() += 1;
                    }
                    if kind == 'S'
                        && let Some(x) = geomaticist_gnss_core::text::py_float(v)
                    {
                        *st.snr_sum.entry(*code).or_default() += x;
                        *st.snr_n.entry(*code).or_default() += 1;
                    }
                }
                if any_obs {
                    st.ticks.push(tick);
                }
                if slip {
                    st.slips.push(tick);
                }
            }
        }
        self.dropped_sats = rd.dropped_sats;
        drop(rd);
        self.data = data;
        if self.epoch_ticks.is_empty() {
            return Err(RinexError::new("no observation epochs found in this file"));
        }
        Ok(())
    }

    /// True for RINEX 2.
    #[must_use]
    pub fn is_v2(&self) -> bool {
        self.header.major() < 3
    }

    /// First epoch (the smallest when out of order).
    #[must_use]
    pub fn first_tick(&self) -> i64 {
        if self.out_of_order > 0 {
            self.epoch_ticks.iter().copied().min().unwrap_or(0)
        } else {
            self.epoch_ticks[0]
        }
    }

    /// Last epoch.
    #[must_use]
    pub fn last_tick(&self) -> i64 {
        if self.out_of_order > 0 {
            self.epoch_ticks.iter().copied().max().unwrap_or(0)
        } else {
            self.epoch_ticks[self.epoch_ticks.len() - 1]
        }
    }

    /// The usual spacing between epochs (the most common positive difference, sampled every
    /// `(n-1)//20000` for huge files; ties: the first seen).
    #[must_use]
    pub fn sample_ticks(&self) -> i64 {
        let t = &self.epoch_ticks;
        let fallback = || (self.header.interval().unwrap_or(1.0) * TICKS as f64) as i64;
        if t.len() < 2 {
            return fallback();
        }
        let mut diffs: Vec<(i64, usize)> = Vec::new(); // in first-seen order
        let step = ((t.len() - 1) / 20_000).max(1);
        let mut i = 0;
        while i < t.len() - 1 {
            let d = t[i + 1] - t[i];
            if d > 0 {
                match diffs.iter_mut().find(|(k, _)| *k == d) {
                    Some(e) => e.1 += 1,
                    None => diffs.push((d, 1)),
                }
            }
            i += step;
        }
        let mut best: Option<(i64, usize)> = None;
        for (d, n) in diffs {
            if best.is_none_or(|(_, bn)| n > bn) {
                best = Some((d, n));
            }
        }
        best.map_or_else(fallback, |(d, _)| d)
    }

    /// Systems present, in RINEX order.
    #[must_use]
    pub fn systems(&self) -> Vec<GnssSystem> {
        let mut v: Vec<GnssSystem> = self.sats.keys().map(|s| s.sys).collect();
        v.sort_by_key(|s| s.order());
        v.dedup();
        debug_assert!(v.iter().all(|s| SYSTEM_ORDER.contains(s.letter())));
        v
    }

    /// The epochs of the file (a fresh reader over the plain text).
    #[must_use]
    pub fn epochs(&self) -> EpochReader<'_> {
        EpochReader::new(&self.data, self.header.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const V3: &str =
        "     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    4 C1C L1C S1C C2W                                      SYS / # / OBS TYPES
E    2 C1C L1C                                              SYS / # / OBS TYPES
    30.000                                                  INTERVAL
                                                            END OF HEADER
> 2024 09 23 00 00  0.0000000  0  2
G01  20000000.123   100000000.12345        45.000
E05  21000000.000 6
> 2024 09 23 00 00 30.0000000  0  1
G01  20000000.223   100000005.1231
> 2024 09 23 00 01  0.0000000  4  1
G    3 C1C L1C S1C                                          SYS / # / OBS TYPES
> 2024 09 23 00 01  0.0000000  0  1
G01  20000000.323   100000010.123          46.000
";

    #[test]
    fn index_v3_with_event() {
        let obs = RinexObs::from_bytes(Path::new("t.rnx"), V3.as_bytes().to_vec()).unwrap();
        assert_eq!(obs.epoch_ticks.len(), 3);
        assert_eq!(obs.events, vec![(Some(obs.epoch_ticks[2]), 4, 1)]);
        let g: Sat = "G01".parse().unwrap();
        let st = &obs.sats[&g];
        assert_eq!(st.counts[&"C1C".parse::<ObsCode>().unwrap()], 3);
        assert_eq!(st.counts[&"S1C".parse::<ObsCode>().unwrap()], 2); // types changed by the event
        assert_eq!(st.slips.len(), 1);
        assert_eq!(st.ssi_sum[&"L1C".parse::<ObsCode>().unwrap()], 5);
        assert_eq!(st.snr_sum[&"S1C".parse::<ObsCode>().unwrap()], 91.0);
        assert_eq!(obs.sample_ticks(), 300_000_000);
        assert_eq!(obs.systems(), vec![GnssSystem::Gps, GnssSystem::Galileo]);
        assert_eq!(obs.obs_types_seen[&b'E'].len(), 2);
    }

    #[test]
    fn sat_ids_like_python() {
        assert_eq!(parse_sat(b"G 5", b'G').unwrap().to_string(), "G05");
        assert_eq!(parse_sat(b" 12", b'R').unwrap().to_string(), "R12");
        assert_eq!(parse_sat(b"   ", b'G').unwrap().to_string(), "G00");
        assert_eq!(parse_sat(b"g07", b'G').unwrap().to_string(), "G07");
        assert!(parse_sat(b"Gab", b'G').is_none());
        assert!(parse_sat(b"X01", b'G').is_none());
    }

    #[test]
    fn truncated_epoch_is_an_error() {
        let text = V3.replace("E05  21000000.000 6\n", "");
        let text = &text[..text.find("> 2024 09 23 00 00 30").unwrap()];
        let r = RinexObs::from_bytes(Path::new("t.rnx"), text.as_bytes().to_vec());
        assert!(r.is_err());
    }
}
