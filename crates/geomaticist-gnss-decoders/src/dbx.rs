//! Leica Captivate job database (DBX): the tables `<job>.X01` … `<job>.X23` and `<job>.XCF` in the job's
//! folder. A clean-room decoder written from our own Captivate jobs.
//!
//! A table is 1024-byte pages. Page 0 is the table header. Every other page: 12 bytes of page links, a
//! flag byte, a byte, u16 entry count, u16 free bytes, then the entries back to back. An entry is
//! `[u8 length − 1][u8 kind] …`; kind 19 holds a data record, kinds 1–6 are index keys (not read here). A data
//! record is
//!
//! ```text
//! [len − 1][19][u8 flag][u32 record number][5 bytes][u8 record type] tagged fields … [ff ff]
//! ```
//!
//! A record longer than one entry continues in the next entries with the same number, flag 0x80 on the
//! last. Each field is a u16 tag and its value; tags rise within a record and an absent field is simply
//! left out. The value's type is fixed per table, record type and tag ([`schema`]); a tag outside the
//! schema ends the read of that record (what was read is kept), so an unknown newer field never breaks an
//! import.
//!
//! X01 kind 30 ("averaged": a point's mean of its measured shots: it follows two measured records of the
//! same name, has no solution type, and its height is their mean) is named, and imported as a position only.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::util::le;

/// Page size.
pub const PAGE: usize = 1024;
/// Entry kind of a data record.
pub const KIND_DATA: u8 = 19;
/// Bytes before a data record's first field.
pub const HEADER: usize = 13;

/// A field's value type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ty {
    /// NUL-terminated Latin-1 text.
    Text,
    /// `[u16 n][n bytes]`.
    Bytes,
    /// No value: present / absent.
    Flag,
    /// u8.
    U8,
    /// u16.
    U16,
    /// i32.
    I32,
    /// f32.
    F32,
    /// f64.
    F64,
}

/// A field's value.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldValue {
    /// Text.
    Text(String),
    /// Bytes.
    Bytes(Vec<u8>),
    /// A present flag.
    Flag,
    /// An integer (u8, u16, i32).
    Int(i64),
    /// A float (f32 exactly, f64).
    Float(f64),
}

impl FieldValue {
    /// The integer, if it is one.
    #[must_use]
    pub fn int(&self) -> Option<i64> {
        match self {
            FieldValue::Int(v) => Some(*v),
            _ => None,
        }
    }
    /// The number (integer or float) as f64.
    #[must_use]
    pub fn num(&self) -> Option<f64> {
        match self {
            FieldValue::Int(v) => Some(*v as f64),
            FieldValue::Float(v) => Some(*v),
            _ => None,
        }
    }
    /// The text, if it is text.
    #[must_use]
    pub fn text(&self) -> Option<&str> {
        match self {
            FieldValue::Text(s) => Some(s),
            _ => None,
        }
    }
    /// The bytes, if they are bytes.
    #[must_use]
    pub fn bytes(&self) -> Option<&[u8]> {
        match self {
            FieldValue::Bytes(b) => Some(b),
            _ => None,
        }
    }
}

/// `(table, record type) -> [(tag, type)]`.
#[must_use]
pub fn schema(table: &str, rtype: u8) -> Option<Vec<(u16, Ty)>> {
    use Ty::{Bytes, F32, F64, Flag, I32, Text, U8, U16};
    Some(match (table, rtype) {
        // X01 points
        ("X01", 0x00) => {
            let mut v = vec![
                (0x00, Text),
                (0x01, U16),
                (0x02, I32),
                (0x04, I32),
                (0x07, F64),
                (0x08, I32),
                (0x0a, U8),
                (0x0c, U8),
                (0x0d, U8),
                (0x0f, U8),
                (0x10, U8),
                (0x11, F64),
                (0x12, F64),
                (0x13, F64),
                (0x16, U8),
            ];
            v.extend((0x17..0x1e).map(|t| (t, F32)));
            v.extend([
                (0x1e, U8),
                (0x1f, U8),
                (0x20, I32),
                (0x21, Text),
                (0x27, Flag),
                (0x28, Bytes),
            ]);
            v
        }
        // X06 codelist: groups (3) and codes (4)
        ("X06", 0x03) => vec![(0x00, Text)],
        ("X06", 0x04) => vec![
            (0x01, I32),
            (0x02, I32),
            (0x03, F64),
            (0x04, F64),
            (0x05, Text),
            (0x07, Text),
        ],
        // X14 point -> code links
        ("X14", 0x2e) => vec![
            (0x00, U8),
            (0x02, I32),
            (0x03, I32),
            (0x04, F64),
            (0x05, Flag),
        ],
        // X12 instruments: controller (0x2a), reference station as received over RTCM (0x35)
        ("X12", 0x2a) => vec![
            (0x01, I32),
            (0x02, Text),
            (0x03, Text),
            (0x04, Text),
            (0x05, Text),
            (0x06, I32),
        ],
        ("X12", 0x35) => vec![
            (0x00, I32),
            (0x01, U8),
            (0x02, U16),
            (0x04, U8),
            (0x05, U8),
            (0x06, F32),
            (0x07, F64),
            (0x08, F64),
            (0x09, F64),
            (0x0b, Text),
            (0x0d, U8),
            (0x0e, U16),
        ],
        // X18 coordinate system: system (0x10), ellipsoid (0x13), projection (0x12), geoid (0x19)
        ("X18", 0x10) => vec![
            (0x00, I32),
            (0x01, I32),
            (0x03, I32),
            (0x04, F64),
            (0x05, U8),
            (0x06, U8),
            (0x08, Text),
        ],
        ("X18", 0x13) => vec![
            (0x00, F64),
            (0x01, Text),
            (0x02, F64),
            (0x03, F64),
            (0x04, U8),
        ],
        ("X18", 0x12) => vec![(0x00, F64), (0x01, U8), (0x02, Text), (0x03, Bytes)],
        ("X18", 0x19) => vec![
            (0x00, U8),
            (0x01, U8),
            (0x02, U8),
            (0x03, Text),
            (0x04, Text),
            (0x05, F64),
        ],
        // X04 photos (0x3c), X05 point -> photo links
        ("X04", 0x3c) => vec![
            (0x00, I32),
            (0x02, I32),
            (0x05, U16),
            (0x06, I32),
            (0x0d, U16),
            (0x0e, U16),
            (0x0f, Text),
            (0x10, Text),
            (0x12, F64),
            (0x14, I32),
        ],
        ("X05", 0x3d) => vec![(0x00, U8), (0x02, I32), (0x03, Flag), (0x04, I32)],
        _ => return None,
    })
}

/// Why a Captivate job couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbxError(pub String);

impl std::fmt::Display for DbxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DbxError {}

/// One data record of a table.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// The table (`X01` …).
    pub table: String,
    /// Record number.
    pub number: u32,
    /// Record type.
    pub rtype: u8,
    /// Fields by tag.
    pub fields: BTreeMap<u16, FieldValue>,
    /// Every field read (false: stopped at a tag outside the schema).
    pub complete: bool,
}

impl Record {
    /// The field with this tag.
    #[must_use]
    pub fn get(&self, tag: u16) -> Option<&FieldValue> {
        self.fields.get(&tag)
    }
    /// The field as a number.
    #[must_use]
    pub fn num(&self, tag: u16) -> Option<f64> {
        self.get(tag).and_then(FieldValue::num)
    }
    /// The field as an integer.
    #[must_use]
    pub fn int(&self, tag: u16) -> Option<i64> {
        self.get(tag).and_then(FieldValue::int)
    }
    /// The field as text ("" when absent).
    #[must_use]
    pub fn text(&self, tag: u16) -> &str {
        self.get(tag).and_then(FieldValue::text).unwrap_or("")
    }
}

/// `(page, kind, entry bytes)` for every entry of every page after the header.
#[must_use]
pub fn entries(data: &[u8]) -> Vec<(usize, u8, &[u8])> {
    let mut out = Vec::new();
    for pg in 1..data.len() / PAGE {
        let p = &data[pg * PAGE..(pg + 1) * PAGE];
        let count = le::u16(p, 14);
        let mut pos = 18usize;
        for _ in 0..count {
            if pos >= PAGE {
                break;
            }
            let n = usize::from(p[pos]) + 1;
            let end = (pos + n).min(PAGE);
            out.push((pg, p.get(pos + 1).copied().unwrap_or(0), &p[pos..end]));
            pos += n;
        }
    }
    out
}

/// `(record number, record type, bytes from the type byte on)` for every data record, continuations joined.
/// Index leaves of kind 19 are 12-byte pointers, never data, and are skipped by their size.
#[must_use]
pub fn raw_records(data: &[u8]) -> Vec<(u32, u8, Vec<u8>)> {
    let mut out = Vec::new();
    let mut part: Vec<(u32, Vec<u8>)> = Vec::new();
    for (_pg, kind, e) in entries(data) {
        if kind != KIND_DATA || e.len() < HEADER {
            continue;
        }
        let num = u32::try_from(le::u32(e, 3)).unwrap_or(0);
        let last = e[2] & 0x80 != 0;
        let k = part.iter().position(|p| p.0 == num);
        if last && k.is_none() {
            if e.ends_with(b"\xff\xff") {
                out.push((num, e[12], e[12..].to_vec()));
            }
            continue;
        }
        // a continued record: the first piece holds the header, later pieces only the next bytes
        match k {
            None => part.push((num, e[12..].to_vec())),
            Some(k) => part[k].1.extend_from_slice(&e[7..]),
        }
        if last && let Some(k) = part.iter().position(|p| p.0 == num) {
            let whole = part.remove(k).1;
            if whole.ends_with(b"\xff\xff") && !whole.is_empty() {
                out.push((num, whole[0], whole));
            }
        }
    }
    out
}

/// Tagged fields of a record body (its first byte is the record type) → (fields, read to the end).
#[must_use]
pub fn parse_fields(body: &[u8], schema: &[(u16, Ty)]) -> (BTreeMap<u16, FieldValue>, bool) {
    let mut out = BTreeMap::new();
    let end = body.len().saturating_sub(2);
    let mut pos = 1usize;
    let mut last: i64 = -1;
    while pos < end {
        if pos + 2 > end {
            return (out, false);
        }
        let tag = u16::try_from(le::u16(body, pos)).unwrap_or(0);
        let Some(&(_, ty)) = schema.iter().find(|s| s.0 == tag) else {
            return (out, false);
        };
        if i64::from(tag) <= last {
            return (out, false);
        }
        pos += 2;
        let need = |n: usize| pos + n <= end;
        let v = match ty {
            Ty::Text => {
                let Some(z) = body[pos..end].iter().position(|&c| c == 0) else {
                    return (out, false);
                };
                let s = crate::util::latin1(&body[pos..pos + z]);
                pos += z + 1;
                FieldValue::Text(s)
            }
            Ty::Bytes => {
                if !need(2) {
                    return (out, false);
                }
                let n = usize::try_from(le::u16(body, pos)).unwrap_or(0);
                let b = body
                    .get(pos + 2..(pos + 2 + n).min(body.len()))
                    .unwrap_or(&[])
                    .to_vec();
                pos += 2 + n;
                FieldValue::Bytes(b)
            }
            Ty::Flag => FieldValue::Flag,
            Ty::U8 | Ty::U16 | Ty::I32 | Ty::F32 | Ty::F64 => {
                let n = match ty {
                    Ty::U8 => 1,
                    Ty::U16 => 2,
                    Ty::I32 | Ty::F32 => 4,
                    _ => 8,
                };
                if !need(n) {
                    return (out, false);
                }
                let v = match ty {
                    Ty::U8 => FieldValue::Int(i64::from(body[pos])),
                    Ty::U16 => FieldValue::Int(le::u16(body, pos)),
                    Ty::I32 => FieldValue::Int(i64::from(le::u32(body, pos) as u32 as i32)),
                    Ty::F32 => FieldValue::Float(le::f32(body, pos)),
                    _ => FieldValue::Float(le::f64(body, pos)),
                };
                pos += n;
                v
            }
        };
        out.insert(tag, v);
        last = i64::from(tag);
    }
    (out, pos == end)
}

/// Every data record of one table's bytes, its fields read where the table / record type has a schema.
pub fn read_table_bytes(data: &[u8], table: &str, name: &str) -> Result<Vec<Record>, DbxError> {
    if data.len() < PAGE || !data.len().is_multiple_of(PAGE) {
        return Err(DbxError(format!(
            "{name} isn't a Captivate DBX table (size {} isn't whole 1024-byte pages).",
            data.len()
        )));
    }
    Ok(raw_records(data)
        .into_iter()
        .map(|(number, rtype, body)| match schema(table, rtype) {
            None => Record {
                table: table.to_owned(),
                number,
                rtype,
                fields: BTreeMap::new(),
                complete: false,
            },
            Some(s) => {
                let (fields, complete) = parse_fields(&body, &s);
                Record {
                    table: table.to_owned(),
                    number,
                    rtype,
                    fields,
                    complete,
                }
            }
        })
        .collect())
}

/// Every data record of one table file (`table`: `X01` …; `None`: from the file's suffix).
pub fn read_table(path: &Path, table: Option<&str>) -> Result<Vec<Record>, DbxError> {
    let name = path
        .file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let t = table.map_or_else(
        || {
            path.extension()
                .map_or_else(String::new, |e| e.to_string_lossy().to_uppercase())
        },
        str::to_owned,
    );
    let data = std::fs::read(path).map_err(|e| DbxError(format!("{name}: {e}")))?;
    read_table_bytes(&data, &t, &name)
}

// ------------------------------------------------------------------ the job

/// What an X01 point record is (field 0x0a).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointKind {
    /// 15: keyed in or kept from an earlier job.
    Stored,
    /// 25: the RTCM reference (the base mark).
    Reference,
    /// 30: a mean of a point's measured shots.
    Averaged,
    /// 35: an RTK measurement.
    Measured,
    /// 40: a navigated (autonomous) fix.
    Navigated,
    /// Anything else.
    Other,
}

impl PointKind {
    /// From field 0x0a.
    #[must_use]
    pub fn of(code: Option<i64>) -> Self {
        match code {
            Some(15) => PointKind::Stored,
            Some(25) => PointKind::Reference,
            Some(30) => PointKind::Averaged,
            Some(35) => PointKind::Measured,
            Some(40) => PointKind::Navigated,
            _ => PointKind::Other,
        }
    }
    /// Its name ("stored", "reference" …).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            PointKind::Stored => "stored",
            PointKind::Reference => "reference",
            PointKind::Averaged => "averaged",
            PointKind::Measured => "measured",
            PointKind::Navigated => "navigated",
            PointKind::Other => "other",
        }
    }
}

/// An X01 point record.
#[derive(Debug, Clone, PartialEq)]
pub struct DbxPoint {
    /// Record number.
    pub number: u32,
    /// Field 0x02: the point's id, shared by its re-measurements.
    pub pid: i64,
    /// The point's name ("" when the record has none).
    pub name: String,
    /// What it is.
    pub kind: PointKind,
    /// Field 0x0a.
    pub kind_code: Option<i64>,
    /// Latitude, longitude (degrees), ellipsoid height (m).
    pub lat: f64,
    /// Longitude, degrees.
    pub lon: f64,
    /// Ellipsoid height, m.
    pub h: f64,
    /// Time stored, MJD in the controller's local time.
    pub mjd: Option<f64>,
    /// Covariance, m², (north, east, up): nn ne nu ee eu uu.
    pub cov: Option<[f64; 6]>,
    /// Field 0x08 present: replaced by a later shot or deleted in the field.
    pub superseded: bool,
    /// Field 0x1e (4 RTK fixed, 10 reference, 2 navigated).
    pub solution: Option<i64>,
    /// Field 0x20: its point record in the job's .m00.
    pub m00_offset: Option<usize>,
    /// Code (X14 → X06).
    pub code: String,
    /// The code's description.
    pub code_desc: String,
    /// Photo file names in the job's Images/ folder.
    pub photos: Vec<String>,
}

/// A reference station as received over RTCM (X12 type 0x35).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DbxReference {
    /// Station id.
    pub station: Option<i64>,
    /// Antenna text (`Tripod`).
    pub antenna: String,
    /// Antenna height as entered at the base, m.
    pub height: Option<f64>,
    /// The broadcast antenna position, ECEF m.
    pub arp_xyz: Option<[f64; 3]>,
}

/// The job's projection (X18 type 0x12): a Transverse Mercator.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DbxProjection {
    /// Its name.
    pub name: String,
    /// False easting, m.
    pub x_0: f64,
    /// False northing, m.
    pub y_0: f64,
    /// Latitude of origin, degrees.
    pub lat_0: f64,
    /// Central meridian, degrees.
    pub lon_0: f64,
    /// Scale.
    pub k_0: f64,
}

/// The job's coordinate system (X18).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DbxCs {
    /// The system's name (`UTM 32N`).
    pub name: String,
    /// Ellipsoid (name, a, f).
    pub ellipsoid: Option<(String, f64, f64)>,
    /// Projection.
    pub projection: Option<DbxProjection>,
    /// Geoid (name, file).
    pub geoid: Option<(String, String)>,
}

/// A Captivate job.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DbxJob {
    /// The job's name (its tables' stem).
    pub name: String,
    /// Its folder.
    pub folder: PathBuf,
    /// X01 points, in record order.
    pub points: Vec<DbxPoint>,
    /// X06 record number → (code, description).
    pub codes: BTreeMap<u32, (String, String)>,
    /// Controller type and serial.
    pub controller: (String, String),
    /// Reference stations received.
    pub references: Vec<DbxReference>,
    /// Coordinate system.
    pub cs: DbxCs,
    /// File names in Images/.
    pub images: Vec<String>,
    /// The job's own raw log.
    pub m00: Option<PathBuf>,
    /// Table → records read.
    pub tables: BTreeMap<String, usize>,
    /// Table → records whose fields weren't all read.
    pub unread: BTreeMap<String, usize>,
}

/// `(folder, stem)` of a Captivate job from its folder, its .XCF or any of its tables.
pub fn job_paths(path: &Path) -> Result<(PathBuf, String), DbxError> {
    if path.is_dir() {
        let mut xcf: Vec<PathBuf> = std::fs::read_dir(path)
            .map_err(|e| DbxError(format!("{}: {e}", path.display())))?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("xcf")))
            .collect();
        xcf.sort();
        let Some(x) = xcf.first() else {
            return Err(DbxError(format!(
                "{} holds no Captivate job (no .XCF file).",
                path.file_name()
                    .map_or_else(String::new, |s| s.to_string_lossy().into_owned())
            )));
        };
        return Ok((
            path.to_path_buf(),
            x.file_stem()
                .map_or_else(String::new, |s| s.to_string_lossy().into_owned()),
        ));
    }
    Ok((
        path.parent().map_or_default(Path::to_path_buf),
        path.file_stem()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned()),
    ))
}

/// A Captivate job: its folder (holding an .XCF), its .XCF or one of its .X01–.X99 tables.
#[must_use]
pub fn is_dbx(path: &Path) -> bool {
    if path.is_dir() {
        return std::fs::read_dir(path).is_ok_and(|rd| {
            rd.filter_map(Result::ok).any(|e| {
                e.path()
                    .extension()
                    .is_some_and(|x| x.eq_ignore_ascii_case("xcf"))
            })
        });
    }
    let s = path
        .extension()
        .map_or_else(String::new, |e| e.to_string_lossy().to_uppercase());
    s == "XCF" || (s.len() == 3 && s.starts_with('X') && s[1..].bytes().all(|b| b.is_ascii_digit()))
}

fn table_path(folder: &Path, stem: &str, t: &str) -> Option<PathBuf> {
    [
        folder.join(format!("{stem}.{t}")),
        folder.join(format!("{stem}.{}", t.to_lowercase())),
    ]
    .into_iter()
    .find(|c| c.is_file())
}

/// Reads a Captivate job: points (X01) with their codes (X14 → X06), instruments and the RTCM reference
/// (X12), coordinate system (X18), and the photos in Images/.
pub fn read_job(path: &Path) -> Result<DbxJob, DbxError> {
    let (folder, stem) = job_paths(path)?;
    let Some(x01) = table_path(&folder, &stem, "X01") else {
        return Err(DbxError(format!(
            "Job {stem}: no point table ({stem}.X01) beside it."
        )));
    };
    let mut job = DbxJob {
        name: stem.clone(),
        folder: folder.clone(),
        ..DbxJob::default()
    };
    let recs = read_table(&x01, Some("X01"))?;
    job.tables.insert("X01".into(), recs.len());
    for r in &recs {
        if r.rtype != 0 || r.get(0x11).is_none() {
            continue;
        }
        if !r.complete {
            *job.unread.entry("X01".into()).or_insert(0) += 1;
        }
        let name = r.text(0x00).split('\n').next().unwrap_or("").to_owned();
        let k = r.int(0x0a);
        let cov = (0x18..0x1e)
            .map(|t| r.num(t))
            .collect::<Option<Vec<f64>>>()
            .map(|v| [v[0], v[1], v[2], v[3], v[4], v[5]]);
        job.points.push(DbxPoint {
            number: r.number,
            pid: r.int(0x02).unwrap_or(i64::from(r.number)),
            name,
            kind: PointKind::of(k),
            kind_code: k,
            lat: r.num(0x11).unwrap_or(0.0).to_degrees(),
            lon: r.num(0x12).unwrap_or(0.0).to_degrees(),
            h: r.num(0x13).unwrap_or(0.0),
            mjd: r.num(0x07),
            cov,
            superseded: r.get(0x08).is_some(),
            solution: r.int(0x1e),
            m00_offset: r.int(0x20).and_then(|v| usize::try_from(v).ok()),
            code: String::new(),
            code_desc: String::new(),
            photos: Vec::new(),
        });
    }

    // codes: X06 codes, X14 links point id -> code record
    if let Some(x06) = table_path(&folder, &stem, "X06") {
        for r in read_table(&x06, Some("X06"))? {
            let key = r.text(0x07);
            if r.rtype == 0x04 && !key.is_empty() {
                let code = key
                    .chars()
                    .skip(1)
                    .collect::<String>()
                    .split('\n')
                    .next()
                    .unwrap_or("")
                    .to_owned();
                job.codes.insert(r.number, (code, r.text(0x05).to_owned()));
            }
        }
    }
    if let Some(x14) = table_path(&folder, &stem, "X14") {
        let mut links: BTreeMap<i64, Vec<i64>> = BTreeMap::new();
        for r in read_table(&x14, Some("X14"))? {
            if r.rtype == 0x2e
                && let Some(pid) = r.int(0x02)
            {
                links
                    .entry(pid)
                    .or_default()
                    .push(r.int(0x03).unwrap_or(-1));
            }
        }
        for p in &mut job.points {
            for c in links.get(&p.pid).map_or(&[][..], Vec::as_slice) {
                if let Some((code, desc)) = u32::try_from(*c).ok().and_then(|c| job.codes.get(&c)) {
                    p.code.clone_from(code);
                    p.code_desc.clone_from(desc);
                    break;
                }
            }
        }
    }

    if let Some(x12) = table_path(&folder, &stem, "X12") {
        for r in read_table(&x12, Some("X12"))? {
            if r.rtype == 0x2a {
                let serial = if r.text(0x05).is_empty() {
                    r.text(0x02)
                } else {
                    r.text(0x05)
                };
                job.controller = (r.text(0x03).to_owned(), serial.to_owned());
            } else if r.rtype == 0x35 && r.get(0x07).is_some() {
                let xyz = match (r.num(0x07), r.num(0x08), r.num(0x09)) {
                    (Some(x), Some(y), Some(z)) => Some([x, y, z]),
                    _ => None,
                };
                job.references.push(DbxReference {
                    station: r.int(0x00),
                    antenna: r.text(0x0b).to_owned(),
                    height: r
                        .num(0x06)
                        .map(|h| geomaticist_gnss_core::text::round_to(h, 4)),
                    arp_xyz: xyz,
                });
            }
        }
    }

    if let Some(x18) = table_path(&folder, &stem, "X18") {
        for r in read_table(&x18, Some("X18"))? {
            match r.rtype {
                0x10 => job.cs.name = r.text(0x08).to_owned(),
                0x13 if r.num(0x00).is_some_and(|a| a != 0.0) => {
                    job.cs.ellipsoid = Some((
                        r.text(0x01).to_owned(),
                        r.num(0x00).unwrap_or(0.0),
                        r.num(0x02).unwrap_or(0.0),
                    ));
                }
                0x12 => {
                    if let Some(b) = r.get(0x03).and_then(FieldValue::bytes)
                        && b.len() >= 48
                    {
                        let f = |k: usize| le::f64(b, 8 * k);
                        job.cs.projection = Some(DbxProjection {
                            name: r.text(0x02).to_owned(),
                            x_0: f(0),
                            y_0: f(1),
                            lat_0: f(2).to_degrees(),
                            lon_0: f(3).to_degrees(),
                            k_0: f(4),
                        });
                    }
                }
                0x19 => {
                    job.cs.geoid = Some((r.text(0x03).to_owned(), r.text(0x04).to_owned()));
                }
                _ => {}
            }
        }
    }

    let imgs = folder.join("Images");
    if imgs.is_dir()
        && let Ok(rd) = std::fs::read_dir(&imgs)
    {
        job.images = rd
            .filter_map(Result::ok)
            .filter(|e| e.path().is_file())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        job.images.sort();
    }
    if let (Some(x04), Some(x05)) = (
        table_path(&folder, &stem, "X04"),
        table_path(&folder, &stem, "X05"),
    ) {
        let files: BTreeMap<u32, String> = read_table(&x04, Some("X04"))?
            .into_iter()
            .filter(|r| r.rtype == 0x3c && !r.text(0x0f).is_empty())
            .map(|r| {
                let ext = if r.get(0x10).is_some() {
                    r.text(0x10).to_owned()
                } else {
                    "jpg".to_owned()
                };
                (r.number, format!("{}.{ext}", r.text(0x0f)))
            })
            .collect();
        let mut by_pid: BTreeMap<i64, Vec<String>> = BTreeMap::new();
        for r in read_table(&x05, Some("X05"))? {
            if r.rtype == 0x3d
                && let Some(f) = r
                    .int(0x04)
                    .and_then(|n| u32::try_from(n).ok())
                    .and_then(|n| files.get(&n))
            {
                by_pid
                    .entry(r.int(0x02).unwrap_or(-1))
                    .or_default()
                    .push(f.clone());
            }
        }
        for p in &mut job.points {
            if !p.superseded {
                p.photos = by_pid.get(&p.pid).cloned().unwrap_or_default();
            }
        }
    }
    job.m00 = table_path(&folder, &stem, "m00");
    Ok(job)
}

#[doc(hidden)]
#[allow(clippy::must_use_candidate)]
pub mod synth {
    //! Synthetic DBX tables for tests (no real data); public so the app's tests can build them too.

    /// A field: tag and its value bytes.
    pub fn field(tag: u16, val: &[u8]) -> Vec<u8> {
        let mut v = tag.to_le_bytes().to_vec();
        v.extend_from_slice(val);
        v
    }

    /// Text value bytes.
    pub fn text(s: &str) -> Vec<u8> {
        let mut v: Vec<u8> = s.bytes().collect();
        v.push(0);
        v
    }

    /// A data entry: `[len − 1][19][flag][u32 number][5 bytes][type] fields [ff ff]`.
    pub fn entry(flag: u8, number: u32, rtype: u8, fields: &[u8]) -> Vec<u8> {
        let mut e = vec![0, super::KIND_DATA, flag];
        e.extend_from_slice(&number.to_le_bytes());
        e.extend_from_slice(&[0, 1, 2, 3, 4, rtype]);
        e.extend_from_slice(fields);
        e.extend_from_slice(&[0xff, 0xff]);
        e[0] = u8::try_from(e.len() - 1).unwrap_or(0);
        e
    }

    /// A table: a header page, then the entries packed into pages.
    pub fn table(entries: &[Vec<u8>]) -> Vec<u8> {
        let mut d = vec![0u8; super::PAGE];
        d[4] = 0x56;
        let mut i = 0;
        while i < entries.len() {
            let mut p = vec![0u8; super::PAGE];
            p[4..12].copy_from_slice(&[0xff; 8]);
            let mut pos = 18;
            let mut n = 0u16;
            while i < entries.len() && pos + entries[i].len() <= super::PAGE {
                p[pos..pos + entries[i].len()].copy_from_slice(&entries[i]);
                pos += entries[i].len();
                n += 1;
                i += 1;
            }
            p[14..16].copy_from_slice(&n.to_le_bytes());
            p[16..18].copy_from_slice(&u16::try_from(super::PAGE - pos).unwrap_or(0).to_le_bytes());
            d.extend(p);
        }
        d
    }

    /// An X01 point record's fields.
    pub fn x01_point(
        name: &str,
        pid: i32,
        kind: u8,
        llh_deg: (f64, f64, f64),
        mjd: f64,
        superseded: bool,
        m00_offset: Option<i32>,
    ) -> Vec<u8> {
        let mut f = field(0x00, &text(&format!("{name}\n{kind}\n0")));
        f.extend(field(0x01, &0xc000u16.to_le_bytes()));
        f.extend(field(0x02, &pid.to_le_bytes()));
        f.extend(field(0x07, &mjd.to_le_bytes()));
        if superseded {
            f.extend(field(0x08, &1i32.to_le_bytes()));
        }
        f.extend(field(0x0a, &[kind]));
        f.extend(field(0x11, &llh_deg.0.to_radians().to_le_bytes()));
        f.extend(field(0x12, &llh_deg.1.to_radians().to_le_bytes()));
        f.extend(field(0x13, &llh_deg.2.to_le_bytes()));
        if kind == 35 {
            for (k, v) in [1e-4f32, 0.0, 0.0, 1e-4, 0.0, 4e-4].iter().enumerate() {
                f.extend(field(
                    0x18 + u16::try_from(k).unwrap_or(0),
                    &v.to_le_bytes(),
                ));
            }
            f.extend(field(0x1e, &[4]));
        }
        if let Some(o) = m00_offset {
            f.extend(field(0x20, &o.to_le_bytes()));
            f.extend(field(0x21, &text("m00")));
        }
        f
    }
}

#[cfg(test)]
mod tests {
    use super::synth::*;
    use super::*;

    #[test]
    fn pages_and_records() {
        let p1 = x01_point("A1", 7, 35, (46.0, 8.0, 10.0), 58700.5, false, Some(1234));
        // a long record split over two entries (flag 0 then 0x80), the second piece after its 7-byte head
        let p2 = x01_point("B2", 8, 15, (46.1, 8.1, 20.0), 58700.4, false, None);
        let (a, b) = p2.split_at(20);
        let mut first = entry(0x00, 2, 0, a);
        // drop the ff ff of the first piece: only the last piece ends the record
        first.truncate(first.len() - 2);
        first[0] = u8::try_from(first.len() - 1).unwrap_or(0);
        let mut second = vec![0, KIND_DATA, 0x80, 2, 0, 0, 0];
        second.extend_from_slice(b);
        second.extend_from_slice(&[0xff, 0xff]);
        second[0] = u8::try_from(second.len() - 1).unwrap_or(0);
        let t = table(&[entry(0x80, 1, 0, &p1), first, second]);
        let recs = read_table_bytes(&t, "X01", "j.X01").unwrap_or_default();
        assert_eq!(recs.len(), 2);
        assert!(recs.iter().all(|r| r.complete));
        assert_eq!(recs[0].text(0x00), "A1\n35\n0");
        assert_eq!(recs[0].int(0x20), Some(1234));
        assert!((recs[0].num(0x11).unwrap_or(0.0).to_degrees() - 46.0).abs() < 1e-12);
        assert_eq!(recs[1].number, 2);
        assert_eq!(recs[1].int(0x0a), Some(15));
        assert!(read_table_bytes(&t[..1500], "X01", "j.X01").is_err());
    }

    #[test]
    fn unknown_tag_stops_the_record() {
        let mut f = field(0x00, &text("C\n35\n0"));
        f.extend(field(0x05, &[1, 2, 3])); // 0x05 isn't in the X01 schema
        let recs =
            read_table_bytes(&table(&[entry(0x80, 3, 0, &f)]), "X01", "j").unwrap_or_default();
        assert_eq!(recs.len(), 1);
        assert!(!recs[0].complete);
        assert_eq!(recs[0].text(0x00), "C\n35\n0");
    }

    #[test]
    fn job_folder() {
        let dir = std::env::temp_dir().join(format!("geom_dbx_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap_or(());
        let w = |n: &str, d: &[u8]| std::fs::write(dir.join(n), d).unwrap_or(());
        w("J.XCF", &[0u8; 393]);
        w(
            "J.X01",
            &table(&[
                entry(
                    0x80,
                    1,
                    0,
                    &x01_point(
                        "RTCM-Ref 0000",
                        1,
                        25,
                        (46.0, 8.0, 10.0),
                        58700.4,
                        false,
                        None,
                    ),
                ),
                entry(
                    0x80,
                    2,
                    0,
                    &x01_point("2001", 2, 35, (46.001, 8.0, 11.0), 58700.5, false, Some(99)),
                ),
            ]),
        );
        let mut proj = 500_000.0f64.to_le_bytes().to_vec();
        for v in [
            0.0f64,
            0f64.to_radians(),
            9f64.to_radians(),
            0.9996,
            4f64.to_radians(),
        ] {
            proj.extend_from_slice(&v.to_le_bytes());
        }
        let mut pf = field(0x02, &text("Zone X"));
        pf.extend(field(0x03, &{
            let mut b = u16::try_from(proj.len())
                .unwrap_or(0)
                .to_le_bytes()
                .to_vec();
            b.extend_from_slice(&proj);
            b
        }));
        w("J.X18", &table(&[entry(0x80, 1, 0x12, &pf)]));
        assert!(is_dbx(&dir));
        assert!(is_dbx(Path::new("x/J.X01")));
        assert!(!is_dbx(Path::new("x/J.m00")));
        let j = read_job(&dir.join("J.XCF")).unwrap_or_default();
        assert_eq!(j.name, "J");
        assert_eq!(j.points.len(), 2);
        assert_eq!(j.points[0].kind, PointKind::Reference);
        assert_eq!(j.points[1].kind, PointKind::Measured);
        assert_eq!(j.points[1].m00_offset, Some(99));
        assert!(j.points[1].cov.is_some());
        let pr = j.cs.projection.clone().unwrap_or_default();
        assert_eq!(pr.x_0, 500_000.0);
        assert!(pr.lat_0.abs() < 1e-9 && (pr.lon_0 - 9.0).abs() < 1e-9);
        assert_eq!(read_job(&dir).map(|j| j.name), Ok("J".to_owned()));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
