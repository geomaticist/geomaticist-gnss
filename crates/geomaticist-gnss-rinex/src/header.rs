//! RINEX observation header: parsing, the fields a surveyor edits, and the formatting helpers.
//!
//! A header is kept as its original lines (content = columns 1-60, label = columns 61-80, both as
//! latin-1 bytes) so anything this program doesn't know is carried through untouched.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use geomaticist_gnss_core::system::SYSTEM_ORDER;
use geomaticist_gnss_core::text::{latin1, py_float, py_rstrip, py_strip, sl};
use geomaticist_gnss_core::time::py_int;
use geomaticist_gnss_core::{ObsCode, Sat};

use crate::RinexError;

/// One header record: content (columns 1-60, padded to 60) and label (columns 61-80, stripped).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeaderLine {
    /// Columns 1-60 (latin-1 bytes).
    pub content: Vec<u8>,
    /// Columns 61-80, stripped.
    pub label: Vec<u8>,
}

impl HeaderLine {
    /// A record from content and label text.
    #[must_use]
    pub fn new(content: impl AsRef<[u8]>, label: impl AsRef<[u8]>) -> Self {
        Self {
            content: content.as_ref().to_vec(),
            label: label.as_ref().to_vec(),
        }
    }

    /// From a raw line: `line[:60].ljust(60)`, `line[60:80].strip()`.
    #[must_use]
    pub fn from_raw(line: &[u8]) -> Self {
        let mut content = sl(line, 0, 60).to_vec();
        content.resize(content.len().max(60), b' ');
        Self {
            content,
            label: py_strip(sl(line, 60, 80)).to_vec(),
        }
    }

    /// True when the stripped label equals `label`.
    #[must_use]
    pub fn is(&self, label: &str) -> bool {
        py_strip(&self.label) == label.as_bytes()
    }

    /// `f"{content[:60]:<60}{label[:20]:<20}"`.
    #[must_use]
    pub fn text(&self) -> Vec<u8> {
        let mut out = sl(&self.content, 0, 60).to_vec();
        out.resize(60, b' ');
        out.extend_from_slice(sl(&self.label, 0, 20));
        out.resize(80, b' ');
        out
    }
}

/// Header fields a person may correct; `None` keeps what the file says.
#[derive(Debug, Clone, Default, PartialEq)]
#[allow(missing_docs)] // field names are the RINEX header fields
pub struct HeaderEdits {
    pub marker_name: Option<String>,
    pub marker_number: Option<String>,
    pub marker_type: Option<String>,
    pub observer: Option<String>,
    pub agency: Option<String>,
    pub rec_number: Option<String>,
    pub rec_type: Option<String>,
    pub rec_version: Option<String>,
    pub ant_number: Option<String>,
    /// Antenna type, 16 characters max (NGS / IGS name).
    pub ant_type: Option<String>,
    /// Radome, 4 characters (e.g. NONE).
    pub ant_radome: Option<String>,
    pub approx_xyz: Option<[f64; 3]>,
    pub delta_hen: Option<[f64; 3]>,
    /// Extra COMMENT lines to add.
    pub comments: Vec<String>,
    pub run_by: Option<String>,
}

impl HeaderEdits {
    /// Names of the fields that are set (not `None`, empty or `[]`), in field order.
    #[must_use]
    pub fn active(&self) -> Vec<&'static str> {
        let s = |v: &Option<String>| v.as_ref().is_some_and(|x| !x.is_empty());
        let mut out = Vec::new();
        for (name, on) in [
            ("marker_name", s(&self.marker_name)),
            ("marker_number", s(&self.marker_number)),
            ("marker_type", s(&self.marker_type)),
            ("observer", s(&self.observer)),
            ("agency", s(&self.agency)),
            ("rec_number", s(&self.rec_number)),
            ("rec_type", s(&self.rec_type)),
            ("rec_version", s(&self.rec_version)),
            ("ant_number", s(&self.ant_number)),
            ("ant_type", s(&self.ant_type)),
            ("ant_radome", s(&self.ant_radome)),
            ("approx_xyz", self.approx_xyz.is_some()),
            ("delta_hen", self.delta_hen.is_some()),
            ("comments", !self.comments.is_empty()),
            ("run_by", s(&self.run_by)),
        ] {
            if on {
                out.push(name);
            }
        }
        out
    }
}

/// A RINEX observation header.
#[derive(Debug, Clone, PartialEq)]
pub struct RinexHeader {
    /// The records, in file order, END OF HEADER last.
    pub lines: Vec<HeaderLine>,
    /// Version as written (`"3.04"`).
    pub version_str: String,
    /// Version as a number.
    pub version: f64,
    /// File type letter (`O`).
    pub file_type: u8,
    /// Satellite system letter (`G`, `M` ...).
    pub system: u8,
    /// Observation codes per system letter (RINEX 2: the one list copied to each system).
    pub obs_types: BTreeMap<u8, Vec<ObsCode>>,
    /// RINEX 2 observation codes.
    pub v2_types: Vec<ObsCode>,
    /// Bumped whenever a flag-4 event replaces the observation types, so a
    /// reader can tell the types changed.
    pub types_generation: u32,
}

impl RinexHeader {
    /// Parses the records; fails on a header that cannot be read.
    pub fn new(lines: Vec<HeaderLine>) -> Result<Self, RinexError> {
        let first = lines
            .first()
            .cloned()
            .unwrap_or_else(|| HeaderLine::new("", ""));
        if !first.is("RINEX VERSION / TYPE") {
            return Err(RinexError::new(
                "not a RINEX file: the first line isn't 'RINEX VERSION / TYPE'",
            ));
        }
        let version_str = latin1(py_strip(sl(&first.content, 0, 9)));
        let version = py_float(version_str.as_bytes())
            .ok_or_else(|| RinexError::new(format!("unreadable RINEX version {version_str:?}")))?;
        let file_type = sl(&first.content, 20, 21)
            .first()
            .map_or(0, u8::to_ascii_uppercase);
        let system = py_strip(sl(&first.content, 40, 41))
            .first()
            .map_or(b'G', u8::to_ascii_uppercase);
        if file_type != b'O' {
            return Err(RinexError::new(format!(
                "this is a RINEX '{}' file, not observation data",
                latin1(py_strip(sl(&first.content, 20, 40)))
            )));
        }
        let mut h = Self {
            lines,
            version_str,
            version,
            file_type,
            system,
            obs_types: BTreeMap::new(),
            v2_types: Vec::new(),
            types_generation: 0,
        };
        h.parse_obs_types()?;
        Ok(h)
    }

    /// `int(self.version)`.
    #[must_use]
    pub fn major(&self) -> i64 {
        self.version.trunc() as i64
    }

    /// Content of the first record with this label.
    #[must_use]
    pub fn get(&self, label: &str) -> Option<&[u8]> {
        self.lines
            .iter()
            .find(|l| l.is(label))
            .map(|l| l.content.as_slice())
    }

    /// Contents of every record with this label.
    #[must_use]
    pub fn get_all(&self, label: &str) -> Vec<&[u8]> {
        self.lines
            .iter()
            .filter(|l| l.is(label))
            .map(|l| l.content.as_slice())
            .collect()
    }

    fn parse_obs_types(&mut self) -> Result<(), RinexError> {
        let bad = || RinexError::new("unreadable observation types record");
        if self.major() >= 3 {
            let mut cur: Option<u8> = None;
            let mut want = 0usize;
            let recs: Vec<Vec<u8>> = self
                .get_all("SYS / # / OBS TYPES")
                .into_iter()
                .map(<[u8]>::to_vec)
                .collect();
            for c in recs {
                if !py_strip(sl(&c, 0, 1)).is_empty() {
                    let s = c[0].to_ascii_uppercase();
                    cur = Some(s);
                    let n = sl(&c, 3, 6);
                    want = if n.is_empty() {
                        0
                    } else {
                        usize::try_from(py_int(n).ok_or_else(bad)?).unwrap_or(0)
                    };
                    self.obs_types.insert(s, Vec::new());
                }
                let Some(s) = cur else { continue };
                let list = self.obs_types.entry(s).or_default();
                for i in 0..13 {
                    let code = py_strip(sl(&c, 7 + 4 * i, 10 + 4 * i));
                    if !code.is_empty()
                        && list.len() < want
                        && let Some(oc) = ObsCode::new(code)
                    {
                        list.push(oc);
                    }
                }
            }
        } else {
            let mut want = 0usize;
            let recs: Vec<Vec<u8>> = self
                .get_all("# / TYPES OF OBSERV")
                .into_iter()
                .map(<[u8]>::to_vec)
                .collect();
            for c in recs {
                if !py_strip(sl(&c, 0, 6)).is_empty() {
                    want = usize::try_from(py_int(sl(&c, 0, 6)).ok_or_else(bad)?).unwrap_or(0);
                }
                for i in 0..9 {
                    let code = py_strip(sl(&c, 10 + 6 * i, 12 + 6 * i));
                    if !code.is_empty()
                        && self.v2_types.len() < want
                        && let Some(oc) = ObsCode::new(code)
                    {
                        self.v2_types.push(oc);
                    }
                }
            }
            let systems: Vec<u8> = if self.system == b'M' {
                SYSTEM_ORDER.bytes().collect()
            } else {
                vec![self.system]
            };
            for s in systems {
                self.obs_types.insert(s, self.v2_types.clone());
            }
        }
        Ok(())
    }

    /// An epoch with event flag 4 can carry new header records (new observation types).
    pub fn update_from_event(&mut self, lines: &[HeaderLine]) -> Result<(), RinexError> {
        let has = |l: &str| lines.iter().any(|x| x.is(l));
        if has("SYS / # / OBS TYPES") || has("# / TYPES OF OBSERV") {
            let is_types =
                |x: &HeaderLine| x.is("SYS / # / OBS TYPES") || x.is("# / TYPES OF OBSERV");
            let keep: Vec<HeaderLine> = self
                .lines
                .iter()
                .filter(|l| !is_types(l))
                .cloned()
                .collect();
            let mut new = keep[..keep.len().saturating_sub(1)].to_vec();
            // Insert the event lines whose label is in the event's own label set: all of
            // them (comments and markers included), before END OF HEADER.
            new.extend(lines.iter().cloned());
            new.extend(keep.last().cloned());
            self.lines = new;
            self.obs_types.clear();
            self.v2_types.clear();
            self.parse_obs_types()?;
            self.types_generation += 1;
        }
        Ok(())
    }

    fn field(&self, label: &str, a: usize, b: usize) -> String {
        self.get(label)
            .map_or(String::new(), |c| latin1(py_strip(sl(c, a, b))))
    }

    /// MARKER NAME.
    #[must_use]
    pub fn marker_name(&self) -> String {
        self.field("MARKER NAME", 0, 60)
    }
    /// MARKER NUMBER.
    #[must_use]
    pub fn marker_number(&self) -> String {
        self.field("MARKER NUMBER", 0, 20)
    }
    /// MARKER TYPE.
    #[must_use]
    pub fn marker_type(&self) -> String {
        self.field("MARKER TYPE", 0, 20)
    }
    /// Observer.
    #[must_use]
    pub fn observer(&self) -> String {
        self.field("OBSERVER / AGENCY", 0, 20)
    }
    /// Agency.
    #[must_use]
    pub fn agency(&self) -> String {
        self.field("OBSERVER / AGENCY", 20, 60)
    }
    /// Receiver number.
    #[must_use]
    pub fn rec_number(&self) -> String {
        self.field("REC # / TYPE / VERS", 0, 20)
    }
    /// Receiver type.
    #[must_use]
    pub fn rec_type(&self) -> String {
        self.field("REC # / TYPE / VERS", 20, 40)
    }
    /// Receiver firmware version.
    #[must_use]
    pub fn rec_version(&self) -> String {
        self.field("REC # / TYPE / VERS", 40, 60)
    }
    /// Antenna number.
    #[must_use]
    pub fn ant_number(&self) -> String {
        self.field("ANT # / TYPE", 0, 20)
    }
    /// Antenna type.
    #[must_use]
    pub fn ant_type(&self) -> String {
        self.field("ANT # / TYPE", 20, 36)
    }
    /// Radome.
    #[must_use]
    pub fn ant_radome(&self) -> String {
        self.field("ANT # / TYPE", 36, 40)
    }
    /// Run by.
    #[must_use]
    pub fn run_by(&self) -> String {
        self.field("PGM / RUN BY / DATE", 20, 40)
    }
    /// Program.
    #[must_use]
    pub fn program(&self) -> String {
        self.field("PGM / RUN BY / DATE", 0, 20)
    }

    fn floats(&self, label: &str) -> Option<[f64; 3]> {
        let c = self.get(label)?;
        Some([
            py_float(sl(c, 0, 14))?,
            py_float(sl(c, 14, 28))?,
            py_float(sl(c, 28, 42))?,
        ])
    }

    /// APPROX POSITION XYZ (`None` when unreadable).
    #[must_use]
    pub fn approx_xyz(&self) -> Option<[f64; 3]> {
        self.floats("APPROX POSITION XYZ")
    }

    /// ANTENNA: DELTA H/E/N.
    #[must_use]
    pub fn delta_hen(&self) -> Option<[f64; 3]> {
        self.floats("ANTENNA: DELTA H/E/N")
    }

    /// INTERVAL (s).
    #[must_use]
    pub fn interval(&self) -> Option<f64> {
        self.get("INTERVAL").and_then(|c| py_float(sl(c, 0, 10)))
    }

    /// Time system of TIME OF FIRST OBS, else by the file's system.
    #[must_use]
    pub fn time_system(&self) -> String {
        let ts = self.field("TIME OF FIRST OBS", 48, 51);
        if !ts.is_empty() {
            return ts;
        }
        match self.system {
            b'R' => "GLO",
            b'E' => "GAL",
            b'C' => "BDT",
            b'J' => "QZS",
            b'I' => "IRN",
            _ => "GPS",
        }
        .into()
    }

    /// COMMENT contents, right-stripped (bytes).
    #[must_use]
    pub fn comments(&self) -> Vec<Vec<u8>> {
        self.get_all("COMMENT")
            .into_iter()
            .map(|c| py_rstrip(c).to_vec())
            .collect()
    }

    /// LEAP SECONDS content (60 characters).
    #[must_use]
    pub fn leap_seconds(&self) -> Option<&[u8]> {
        self.get("LEAP SECONDS")
    }

    /// Observation codes of a system (`obs_types.get(sys, [])`).
    #[must_use]
    pub fn types_for(&self, sys: u8) -> &[ObsCode] {
        self.obs_types.get(&sys).map_or(&[], Vec::as_slice)
    }
}

/// Header lines up to and including END OF HEADER from a file's bytes, split at `\n`, `\r` or
/// `\r\n`. Returns the lines and the byte offset
/// just after END OF HEADER's line.
pub fn read_header_lines(data: &[u8]) -> Result<Vec<HeaderLine>, RinexError> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let start = i;
        while i < data.len() && data[i] != b'\n' && data[i] != b'\r' {
            i += 1;
        }
        let line = &data[start..i];
        if i < data.len() {
            if data[i] == b'\r' && data.get(i + 1) == Some(&b'\n') {
                i += 1;
            }
            i += 1;
        }
        let hl = HeaderLine::from_raw(line);
        let end = hl.is("END OF HEADER");
        out.push(hl);
        if end {
            return Ok(out);
        }
        if out.len() > 5000 {
            break;
        }
    }
    Err(RinexError::new(
        "no END OF HEADER found - this isn't a complete RINEX header",
    ))
}

/// `f"{x:14.4f}"`.
#[must_use]
pub fn f14(x: f64) -> String {
    format!("{x:14.4}")
}

/// SYS / # / OBS TYPES records.
#[must_use]
pub fn obs_types_lines_v3(types: &BTreeMap<u8, Vec<ObsCode>>) -> Vec<HeaderLine> {
    let mut order: Vec<u8> = SYSTEM_ORDER.bytes().collect();
    order.extend(
        types
            .keys()
            .filter(|k| !SYSTEM_ORDER.as_bytes().contains(k)),
    );
    let mut out = Vec::new();
    for s in order {
        let Some(codes) = types.get(&s) else { continue };
        if codes.is_empty() {
            continue;
        }
        for (k, chunk) in codes.chunks(13).enumerate() {
            let mut line = if k == 0 {
                let mut h = vec![s];
                h.extend_from_slice(format!("  {:3}", codes.len()).as_bytes());
                h
            } else {
                vec![b' '; 6]
            };
            for c in chunk {
                line.push(b' ');
                line.extend_from_slice(format!("{:<3}", c.as_str()).as_bytes());
            }
            out.push(HeaderLine::new(line, "SYS / # / OBS TYPES"));
        }
    }
    out
}

/// # / TYPES OF OBSERV records.
#[must_use]
pub fn obs_types_lines_v2(codes: &[ObsCode]) -> Vec<HeaderLine> {
    let mut out = Vec::new();
    let n = codes.len().max(1);
    let mut i = 0;
    while i < n {
        let chunk = &codes[i.min(codes.len())..(i + 9).min(codes.len())];
        let mut line = if i == 0 {
            format!("{:6}", codes.len())
        } else {
            " ".repeat(6)
        };
        for c in chunk {
            let _ = write!(line, "    {:<2}", c.as_str());
        }
        out.push(HeaderLine::new(line, "# / TYPES OF OBSERV"));
        i += 9;
    }
    out
}

/// PRN / # OF OBS records, satellites in RINEX order.
#[must_use]
pub fn prn_obs_lines(counts: &BTreeMap<Sat, Vec<usize>>, major: i64) -> Vec<HeaderLine> {
    let mut sats: Vec<&Sat> = counts.keys().collect();
    sats.sort_by_key(|s| s.rinex_order());
    let mut out = Vec::new();
    for sat in sats {
        let n = &counts[sat];
        let sid = if major >= 3 {
            sat.to_string()
        } else {
            format!("{}{:2}", sat.sys.letter(), sat.prn)
        };
        let total = n.len().max(1);
        let mut i = 0;
        while i < total {
            let chunk = &n[i.min(n.len())..(i + 9).min(n.len())];
            let mut line = if i == 0 {
                format!("   {sid}")
            } else {
                " ".repeat(6)
            };
            for v in chunk {
                let _ = write!(line, "{v:6}");
            }
            out.push(HeaderLine::new(line, "PRN / # OF OBS"));
            i += 9;
        }
    }
    out
}

const A_GRS80: f64 = 6_378_137.0;
const F_GRS80: f64 = 1.0 / 298.257_222_101;

fn e2() -> f64 {
    F_GRS80 * (2.0 - F_GRS80)
}

/// Geodetic (degrees, m) to ECEF on GRS80.
#[must_use]
pub fn llh_to_xyz(lat_deg: f64, lon_deg: f64, h: f64) -> [f64; 3] {
    let (lat, lon) = (lat_deg.to_radians(), lon_deg.to_radians());
    let n = A_GRS80 / (1.0 - e2() * lat.sin().powi(2)).sqrt();
    [
        (n + h) * lat.cos() * lon.cos(),
        (n + h) * lat.cos() * lon.sin(),
        (n * (1.0 - e2()) + h) * lat.sin(),
    ]
}

/// ECEF to geodetic (degrees, m) on GRS80.
#[must_use]
pub fn xyz_to_llh(x: f64, y: f64, z: f64) -> [f64; 3] {
    let e2 = e2();
    let lon = y.atan2(x);
    let p = x.hypot(y);
    let mut lat = z.atan2(p * (1.0 - e2));
    for _ in 0..10 {
        let n = A_GRS80 / (1.0 - e2 * lat.sin().powi(2)).sqrt();
        let h = p / lat.cos() - n;
        let lat_new = z.atan2(p * (1.0 - e2 * n / (n + h)));
        if (lat_new - lat).abs() < 1e-14 {
            lat = lat_new;
            break;
        }
        lat = lat_new;
    }
    let n = A_GRS80 / (1.0 - e2 * lat.sin().powi(2)).sqrt();
    let h = p / lat.cos() - n;
    [lat.to_degrees(), lon.to_degrees(), h]
}

/// `64° 58' 40.80000" N`.
#[must_use]
pub fn dms(deg: f64, positive: &str, negative: &str) -> String {
    let hemi = if deg >= 0.0 { positive } else { negative };
    let a = deg.abs();
    let mut d = a.trunc() as i64;
    let mut m = ((a - a.trunc()) * 60.0).trunc() as i64;
    let mut s = ((a - a.trunc()) * 60.0 - m as f64) * 60.0;
    if s >= 59.999_995 {
        s = 0.0;
        m += 1;
    }
    if m == 60 {
        m = 0;
        d += 1;
    }
    format!("{d}\u{b0} {m:02}' {s:08.5}\" {hemi}")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hdr(lines: &[&str]) -> RinexHeader {
        let text = lines.join("\n");
        RinexHeader::new(read_header_lines(text.as_bytes()).unwrap()).unwrap()
    }

    fn line(content: &str, label: &str) -> String {
        format!("{content:<60}{label}")
    }

    #[test]
    fn parses_v3_types_and_fields() {
        let h = hdr(&[
            &line(
                "     3.04           OBSERVATION DATA    M (MIXED)",
                "RINEX VERSION / TYPE",
            ),
            &line("FAIR", "MARKER NAME"),
            &line(
                " -2281621.3290 -1453595.7660  5756961.9680",
                "APPROX POSITION XYZ",
            ),
            &line(
                "G   16 C1C L1C D1C S1C C1W L1W D1W S1W C2W L2W D2W S2W C2L",
                "SYS / # / OBS TYPES",
            ),
            &line("       L2L D2L S2L", "SYS / # / OBS TYPES"),
            &line("E    4 C1C L1C C5Q L5Q", "SYS / # / OBS TYPES"),
            &line("", "END OF HEADER"),
        ]);
        assert_eq!(h.major(), 3);
        assert_eq!(h.version_str, "3.04");
        assert_eq!(h.system, b'M');
        assert_eq!(h.types_for(b'G').len(), 16);
        assert_eq!(h.types_for(b'G')[15].as_str(), "S2L");
        assert_eq!(h.marker_name(), "FAIR");
        assert_eq!(h.approx_xyz().unwrap()[2], 5_756_961.968);
        assert_eq!(h.time_system(), "GPS");
        let back = obs_types_lines_v3(&h.obs_types);
        assert_eq!(back.len(), 3);
        assert_eq!(
            latin1(&back[0].text()).trim_end(),
            line(
                "G   16 C1C L1C D1C S1C C1W L1W D1W S1W C2W L2W D2W S2W C2L",
                "SYS / # / OBS TYPES"
            )
        );
    }

    #[test]
    fn parses_v2_types() {
        let h = hdr(&[
            &line(
                "     2.11           OBSERVATION DATA    M (MIXED)",
                "RINEX VERSION / TYPE",
            ),
            &line(
                "    10    L1    L2    C1    P1    P2    S1    S2    C5    L5",
                "# / TYPES OF OBSERV",
            ),
            &line("          S5", "# / TYPES OF OBSERV"),
            &line("", "END OF HEADER"),
        ]);
        assert_eq!(h.v2_types.len(), 10);
        assert_eq!(h.types_for(b'E').len(), 10);
        assert_eq!(h.types_for(b'S').len(), 10);
        let back = obs_types_lines_v2(&h.v2_types);
        assert_eq!(back.len(), 2);
        assert!(RinexHeader::new(read_header_lines(b"x\n").unwrap_or_default()).is_err());
    }

    #[test]
    fn geodesy_and_dms() {
        let xyz = llh_to_xyz(64.978, -147.499, 319.0);
        let llh = xyz_to_llh(xyz[0], xyz[1], xyz[2]);
        assert!((llh[0] - 64.978).abs() < 1e-11 && (llh[2] - 319.0).abs() < 1e-6);
        assert_eq!(dms(64.5, "N", "S"), "64\u{b0} 30' 00.00000\" N");
        assert_eq!(f14(1.5), "        1.5000");
    }

    #[test]
    fn edits_active() {
        let e = HeaderEdits {
            marker_name: Some("X".into()),
            approx_xyz: Some([0.0; 3]),
            run_by: Some(String::new()),
            ..HeaderEdits::default()
        };
        assert_eq!(e.active(), vec!["marker_name", "approx_xyz"]);
    }
}
