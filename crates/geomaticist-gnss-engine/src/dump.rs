//! Per-stage debug dumps: `NN_stage.npz` (a STORED zip of `.npy` arrays, keys sorted, fixed entry
//! dates) and `NN_stage.json` (sorted keys, indent 1, NaN as null, ±inf as the strings "inf" / "-inf"),
//! so two runs can be compared file by file.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use serde_json::{Map, Value};

use crate::arr::{B2, F2, I2, Vec3Arr};
use crate::linalg::DMat;
use crate::stages::{
    ArcsStage, FloatStage, NarrowLaneStage, ObservationsStage, SigmaStage, Stages, ValidationStage,
    WideLaneStage,
};

/// One array to write: shape and typed data.
#[derive(Debug, Clone, PartialEq)]
pub enum NpyOut {
    /// `<f8`.
    F64(Vec<usize>, Vec<f64>),
    /// `<i8`.
    I64(Vec<usize>, Vec<i64>),
    /// `|b1`.
    Bool(Vec<usize>, Vec<bool>),
    /// `<U{n}` (fixed-width UCS-4).
    Str(Vec<usize>, Vec<String>),
}

fn npy_bytes(a: &NpyOut) -> Vec<u8> {
    let (descr, shape, body): (String, &Vec<usize>, Vec<u8>) = match a {
        NpyOut::F64(s, v) => (
            "<f8".into(),
            s,
            v.iter().flat_map(|x| x.to_le_bytes()).collect(),
        ),
        NpyOut::I64(s, v) => (
            "<i8".into(),
            s,
            v.iter().flat_map(|x| x.to_le_bytes()).collect(),
        ),
        NpyOut::Bool(s, v) => ("|b1".into(), s, v.iter().map(|&b| u8::from(b)).collect()),
        NpyOut::Str(s, v) => {
            let w = v
                .iter()
                .map(|x| x.chars().count())
                .max()
                .unwrap_or(1)
                .max(1);
            let mut b = Vec::with_capacity(v.len() * w * 4);
            for x in v {
                let mut n = 0;
                for c in x.chars() {
                    b.extend_from_slice(&(c as u32).to_le_bytes());
                    n += 1;
                }
                for _ in n..w {
                    b.extend_from_slice(&0u32.to_le_bytes());
                }
            }
            (format!("<U{w}"), s, b)
        }
    };
    let shape_s = if shape.len() == 1 {
        format!("({},)", shape[0])
    } else {
        format!(
            "({})",
            shape
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let mut header =
        format!("{{'descr': '{descr}', 'fortran_order': False, 'shape': {shape_s}, }}");
    while (10 + header.len() + 1) % 64 != 0 {
        header.push(' ');
    }
    header.push('\n');
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend_from_slice(&(header.len() as u16).to_le_bytes());
    out.extend_from_slice(header.as_bytes());
    out.extend_from_slice(&body);
    out
}

fn crc32(data: &[u8]) -> u32 {
    let mut table = [0u32; 256];
    for i in 0..256u32 {
        let mut c = i;
        for _ in 0..8 {
            c = if c & 1 != 0 {
                0xEDB8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
        }
        table[i as usize] = c;
    }
    let mut crc = 0xFFFF_FFFFu32;
    for &b in data {
        crc = table[((crc ^ u32::from(b)) & 0xFF) as usize] ^ (crc >> 8);
    }
    crc ^ 0xFFFF_FFFF
}

/// A STORED zip of `.npy` members (keys sorted), entry dates 1980-01-01.
#[must_use]
pub fn npz_bytes(arrays: &BTreeMap<String, NpyOut>) -> Vec<u8> {
    let mut out = Vec::new();
    let mut central = Vec::new();
    let mut n = 0u16;
    for (key, a) in arrays {
        let name = format!("{key}.npy");
        let data = npy_bytes(a);
        let crc = crc32(&data);
        let off = out.len() as u32;
        let size = data.len() as u32;
        // local file header
        out.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        out.extend_from_slice(&20u16.to_le_bytes()); // version needed
        out.extend_from_slice(&0u16.to_le_bytes()); // flags
        out.extend_from_slice(&0u16.to_le_bytes()); // method: stored
        out.extend_from_slice(&0u16.to_le_bytes()); // time 00:00:00
        out.extend_from_slice(&0x0021u16.to_le_bytes()); // date 1980-01-01
        out.extend_from_slice(&crc.to_le_bytes());
        out.extend_from_slice(&size.to_le_bytes());
        out.extend_from_slice(&size.to_le_bytes());
        out.extend_from_slice(&(name.len() as u16).to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(&data);
        // central directory entry
        central.extend_from_slice(&0x0201_4b50u32.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes()); // made by
        central.extend_from_slice(&20u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0x0021u16.to_le_bytes());
        central.extend_from_slice(&crc.to_le_bytes());
        central.extend_from_slice(&size.to_le_bytes());
        central.extend_from_slice(&size.to_le_bytes());
        central.extend_from_slice(&(name.len() as u16).to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes()); // extra
        central.extend_from_slice(&0u16.to_le_bytes()); // comment
        central.extend_from_slice(&0u16.to_le_bytes()); // disk
        central.extend_from_slice(&0u16.to_le_bytes()); // internal attrs
        central.extend_from_slice(&0u32.to_le_bytes()); // external attrs
        central.extend_from_slice(&off.to_le_bytes());
        central.extend_from_slice(name.as_bytes());
        n += 1;
    }
    let cd_off = out.len() as u32;
    let cd_size = central.len() as u32;
    out.extend_from_slice(&central);
    out.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&n.to_le_bytes());
    out.extend_from_slice(&cd_size.to_le_bytes());
    out.extend_from_slice(&cd_off.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out
}

/// JSON text as `json.dumps(sort_keys=True, indent=1, ensure_ascii=False)` writes it (plus the trailing
/// newline). Keys are sorted because `serde_json::Map` is a `BTreeMap` here.
#[must_use]
pub fn json_text(v: &Value) -> String {
    let mut buf = Vec::new();
    let fmt = serde_json::ser::PrettyFormatter::with_indent(b" ");
    let mut ser = serde_json::Serializer::with_formatter(&mut buf, fmt);
    serde::Serialize::serialize(v, &mut ser).unwrap_or(());
    let mut s = String::from_utf8(buf).unwrap_or_default();
    s.push('\n');
    s
}

// ---- JSON helpers ----------------------------------------------------------------

/// A float as the stage JSON writes it: NaN -> null, ±inf -> "inf" / "-inf".
#[must_use]
pub fn jf(x: f64) -> Value {
    if x.is_nan() {
        Value::Null
    } else if x.is_infinite() {
        Value::String(if x > 0.0 { "inf" } else { "-inf" }.into())
    } else {
        serde_json::Number::from_f64(x).map_or(Value::Null, Value::Number)
    }
}

/// A float as the result JSON writes it: NaN -> "NaN", ±inf -> "Infinity" /
/// "-Infinity". The fix dict is built with this convention; [`stage_nan`] turns it into the stage files' one.
#[must_use]
pub fn jfx(x: f64) -> Value {
    if x.is_nan() {
        Value::String("NaN".into())
    } else if x.is_infinite() {
        Value::String(if x > 0.0 { "Infinity" } else { "-Infinity" }.into())
    } else {
        serde_json::Number::from_f64(x).map_or(Value::Null, Value::Number)
    }
}

/// A float rounded to `nd` decimals, half to even (result convention: NaN -> "NaN").
#[must_use]
pub fn jround(x: f64, nd: usize) -> Value {
    jfx(geomaticist_gnss_core::text::round_to(x, nd))
}

/// The fix dict in the stage files' convention (`debugdump._jsonable`): the result convention's "NaN" becomes
/// null, "Infinity" / "-Infinity" become "inf" / "-inf".
#[must_use]
pub fn stage_nan(v: &Value) -> Value {
    match v {
        Value::String(s) if s == "NaN" => Value::Null,
        Value::String(s) if s == "Infinity" => Value::String("inf".into()),
        Value::String(s) if s == "-Infinity" => Value::String("-inf".into()),
        Value::Array(a) => Value::Array(a.iter().map(stage_nan).collect()),
        Value::Object(m) => {
            Value::Object(m.iter().map(|(k, x)| (k.clone(), stage_nan(x))).collect())
        }
        other => other.clone(),
    }
}

/// A float rounded half to even -> int.
#[must_use]
pub fn jround0(x: f64) -> Value {
    Value::from(crate::nputil::py_round_i(x))
}

/// A float list.
#[must_use]
pub fn jfl(v: &[f64]) -> Value {
    Value::Array(v.iter().map(|&x| jf(x)).collect())
}

/// A 3 × 3 as nested lists.
#[must_use]
pub fn jm3(m: &[[f64; 3]; 3]) -> Value {
    Value::Array(m.iter().map(|r| jfl(r)).collect())
}

/// Strings.
#[must_use]
pub fn jsl<S: AsRef<str>>(v: &[S]) -> Value {
    Value::Array(
        v.iter()
            .map(|s| Value::String(s.as_ref().to_owned()))
            .collect(),
    )
}

/// Integers.
#[must_use]
pub fn jil(v: &[usize]) -> Value {
    Value::Array(v.iter().map(|&x| Value::from(x)).collect())
}

/// Bools.
#[must_use]
pub fn jbl(v: &[bool]) -> Value {
    Value::Array(v.iter().map(|&x| Value::Bool(x)).collect())
}

/// A `DMat` as nested lists.
#[must_use]
pub fn jmat(m: &DMat) -> Value {
    Value::Array((0..m.n).map(|i| jfl(m.row(i))).collect())
}

fn f2(a: &F2) -> NpyOut {
    NpyOut::F64(vec![a.ne, a.ns], a.data.clone())
}
fn b2(a: &B2) -> NpyOut {
    NpyOut::Bool(vec![a.ne, a.ns], a.data.clone())
}
fn i2(a: &I2) -> NpyOut {
    NpyOut::I64(vec![a.ne, a.ns], a.data.clone())
}
fn v3(a: &Vec3Arr) -> NpyOut {
    NpyOut::F64(vec![a.ne, a.ns, 3], a.flat())
}
fn f1(a: &[f64]) -> NpyOut {
    NpyOut::F64(vec![a.len()], a.to_vec())
}
fn rows3(a: &[[f64; 3]]) -> NpyOut {
    NpyOut::F64(vec![a.len(), 3], a.iter().flatten().copied().collect())
}
fn mat(m: &DMat) -> NpyOut {
    NpyOut::F64(vec![m.n, m.n], m.data.clone())
}

/// The arrays and the JSON of a stage, in the dump layout.
pub type StageFiles = (BTreeMap<String, NpyOut>, Value);

fn notes(v: &[String]) -> Value {
    jsl(v)
}

/// Stage 1 files.
#[must_use]
pub fn observations_files(s: &ObservationsStage) -> StageFiles {
    let mut a = BTreeMap::new();
    let mut put = |k: &str, v: NpyOut| {
        a.insert(k.to_owned(), v);
    };
    put("t", f1(&s.t));
    put(
        "cols",
        NpyOut::Str(
            vec![s.cols.len()],
            s.cols.iter().map(ToString::to_string).collect(),
        ),
    );
    for (k, v) in [
        ("phA", &s.ph_a),
        ("phB", &s.ph_b),
        ("cdA", &s.cd_a),
        ("cdB", &s.cd_b),
        ("gfA", &s.gf_a),
        ("gfB", &s.gf_b),
        ("mwA", &s.mw_a),
        ("mwB", &s.mw_b),
        ("elA", &s.el_a),
        ("elB", &s.el_b),
        ("azA", &s.az_a),
        ("azB", &s.az_b),
        ("rhoA", &s.rho_a),
        ("rhoB", &s.rho_b),
        ("clkA", &s.clk_a),
        ("clkB", &s.clk_b),
        ("dryA", &s.dry_a),
        ("dryB", &s.dry_b),
        ("mapA", &s.map_a),
        ("mapB", &s.map_b),
        ("antA", &s.ant_a),
        ("antB", &s.ant_b),
    ] {
        put(k, f2(v));
    }
    for (k, v) in [
        ("okA", &s.ok_a),
        ("okB", &s.ok_b),
        ("valid", &s.valid),
        ("lli", &s.lli),
        ("hard", &s.hard),
        ("halfA", &s.half_a),
        ("halfB", &s.half_b),
        ("msA", &s.ms_a),
        ("msB", &s.ms_b),
    ] {
        put(k, b2(v));
    }
    put("satA", v3(&s.sat_a));
    put("satB", v3(&s.sat_b));
    put("tideA", rows3(&s.tide_a));
    put("tideB", rows3(&s.tide_b));
    put("lam1", f1(&s.lam1));
    put("lam2", f1(&s.lam2));
    put("aif", f1(&s.aif));
    put("bif", f1(&s.bif));
    let mut m = Map::new();
    m.insert("base".into(), Value::String(s.base.clone()));
    m.insert("rover".into(), Value::String(s.rover.clone()));
    m.insert("base_xyz".into(), jfl(&s.base_xyz));
    m.insert("rover_xyz0".into(), jfl(&s.rover_xyz0));
    let pairs: Map<String, Value> = s
        .pairs
        .iter()
        .map(|(k, (b1, b2))| {
            (
                k.to_string(),
                Value::Array(vec![Value::from(*b1), Value::from(*b2)]),
            )
        })
        .collect();
    m.insert("pairs".into(), Value::Object(pairs));
    m.insert("signals_base".into(), s.signals_base.clone());
    m.insert("signals_rover".into(), s.signals_rover.clone());
    let age: Map<String, Value> = s
        .eph_age_min
        .iter()
        .map(|(k, v)| (k.to_string(), jf(*v)))
        .collect();
    m.insert("eph_age_min".into(), Value::Object(age));
    m.insert(
        "events_base".into(),
        serde_json::to_value(s.events_base).unwrap_or(Value::Null),
    );
    m.insert(
        "events_rover".into(),
        serde_json::to_value(s.events_rover).unwrap_or(Value::Null),
    );
    m.insert("notes".into(), notes(&s.notes));
    (a, Value::Object(m))
}

/// Stage 2 files.
#[must_use]
pub fn arcs_files(s: &ArcsStage) -> StageFiles {
    let mut a = BTreeMap::new();
    a.insert("arc".to_owned(), i2(&s.arc));
    a.insert("starts".to_owned(), b2(&s.starts));
    a.insert("valid".to_owned(), b2(&s.valid));
    a.insert("corr".to_owned(), f2(&s.corr));
    a.insert("corr_mw".to_owned(), f2(&s.corr_mw));
    a.insert("corr_gf".to_owned(), f2(&s.corr_gf));
    a.insert("half".to_owned(), b2(&s.half));
    a.insert("ms".to_owned(), b2(&s.ms));
    let mut m = Map::new();
    m.insert("n_amb".into(), Value::from(s.n_amb));
    m.insert("n_arc_raw".into(), Value::from(s.n_arc_raw));
    m.insert(
        "wl_links".into(),
        Value::Array(
            s.wl_links
                .iter()
                .map(|(x, y)| Value::Array(vec![Value::from(*x), Value::from(*y)]))
                .collect(),
        ),
    );
    m.insert(
        "slip_events".into(),
        serde_json::to_value(s.slip_events).unwrap_or(Value::Null),
    );
    let ss: Map<String, Value> = s
        .sat_slips
        .iter()
        .map(|(k, (r, p))| {
            (
                k.to_string(),
                Value::Array(vec![Value::from(*r), Value::from(*p)]),
            )
        })
        .collect();
    m.insert("sat_slips".into(), Value::Object(ss));
    m.insert("notes".into(), notes(&s.notes));
    (a, Value::Object(m))
}

/// Stage 3 files.
#[must_use]
pub fn float_files(s: &FloatStage) -> StageFiles {
    let mut a = BTreeMap::new();
    a.insert("x".to_owned(), f1(&s.x));
    a.insert("N".to_owned(), mat(&s.n));
    a.insert("Q".to_owned(), mat(&s.q));
    a.insert("arc".to_owned(), i2(&s.arc));
    a.insert("use".to_owned(), b2(&s.use_));
    a.insert("wf_ph".to_owned(), f2(&s.wf_ph));
    a.insert("wf_cd".to_owned(), f2(&s.wf_cd));
    a.insert("res_ph".to_owned(), f2(&s.res_ph));
    a.insert("res_cd".to_owned(), f2(&s.res_cd));
    a.insert("cf_ph".to_owned(), f2(&s.cf_ph));
    a.insert("var_ph".to_owned(), f2(&s.var_ph));
    if let Some(v) = &s.var_cd {
        a.insert("var_cd".to_owned(), f2(v));
    }
    let mut m = Map::new();
    m.insert("s0".into(), jf(s.s0));
    m.insert("vector".into(), jfl(&s.vector));
    m.insert("cov".into(), jm3(&s.cov));
    m.insert("n_amb".into(), Value::from(s.n_amb));
    m.insert("iAmb".into(), Value::from(s.i_amb));
    m.insert("iTB".into(), Value::from(s.i_tb));
    m.insert("iTA".into(), Value::from(s.i_ta));
    m.insert("nn".into(), Value::from(s.nn));
    m.insert("nTA".into(), Value::from(s.n_ta));
    m.insert("P".into(), Value::from(s.p));
    m.insert("acc".into(), acc_json(&s.acc));
    let b: Map<String, Value> = s
        .biased
        .iter()
        .map(|(k, v)| (k.to_string(), jf(*v)))
        .collect();
    m.insert("biased".into(), Value::Object(b));
    m.insert("notes".into(), notes(&s.notes));
    (a, Value::Object(m))
}

/// The accumulators as a JSON object (ints as ints).
#[must_use]
pub fn acc_json(acc: &crate::stages::Acc) -> Value {
    let mut m = Map::new();
    m.insert("vtpv".into(), jf(acc.vtpv));
    m.insert("nrow".into(), Value::from(acc.nrow));
    m.insert("ss_ph".into(), jf(acc.ss_ph));
    m.insert("n_ph".into(), Value::from(acc.n_ph));
    m.insert("ss_cd".into(), jf(acc.ss_cd));
    m.insert("n_cd".into(), Value::from(acc.n_cd));
    Value::Object(m)
}

/// Stage 4 files.
#[must_use]
pub fn wide_lane_files(s: &WideLaneStage) -> StageFiles {
    let mut a = BTreeMap::new();
    let mut m = Map::new();
    if let Some(w) = &s.wl {
        a.insert("nw".to_owned(), NpyOut::I64(vec![w.nw.len()], w.nw.clone()));
        a.insert(
            "ok".to_owned(),
            NpyOut::Bool(vec![w.ok.len()], w.ok.clone()),
        );
        a.insert("span".to_owned(), f1(&w.span));
        a.insert("cnt".to_owned(), f1(&w.cnt));
        a.insert("jj".to_owned(), NpyOut::I64(vec![w.jj.len()], w.jj.clone()));
        a.insert(
            "sys".to_owned(),
            NpyOut::Str(vec![w.sys.len()], w.sys.clone()),
        );
        a.insert("mw".to_owned(), f1(&w.mw));
        a.insert("sig".to_owned(), f1(&w.sig));
        a.insert(
            "cand".to_owned(),
            NpyOut::Bool(vec![w.cand.len()], w.cand.clone()),
        );
        let bias: Map<String, Value> = w.bias.iter().map(|(k, v)| (k.clone(), jf(*v))).collect();
        m.insert("bias".into(), Value::Object(bias));
        m.insert("n".into(), jil(&w.n));
        m.insert("groups".into(), Value::from(w.groups));
    }
    m.insert("keep_float".into(), jbl(&s.keep_float));
    m.insert(
        "half_arcs".into(),
        s.half_arcs.as_ref().map_or(Value::Null, |h| jbl(h)),
    );
    m.insert(
        "w_arc".into(),
        s.w_arc.as_ref().map_or(Value::Null, |w| jfl(w)),
    );
    m.insert("pc_arc".into(), jfl(&s.pc_arc));
    (a, Value::Object(m))
}

/// Stage 5 files.
#[must_use]
pub fn narrow_lane_files(s: &NarrowLaneStage) -> StageFiles {
    let mut a = BTreeMap::new();
    if let (Some(x), Some(c)) = (&s.x, &s.cov) {
        a.insert("x".to_owned(), f1(x));
        a.insert("cov".to_owned(), mat(c));
    }
    (a, stage_nan(&Value::Object(s.info.clone())))
}

/// Stage 6 files.
#[must_use]
pub fn validation_files(s: &ValidationStage) -> StageFiles {
    let mut a = BTreeMap::new();
    a.insert("x".to_owned(), f1(&s.x));
    a.insert("cov".to_owned(), mat(&s.cov));
    let mut m = Map::new();
    m.insert("fix".into(), stage_nan(&Value::Object(s.fix.clone())));
    m.insert("failed".into(), Value::String(s.failed.clone()));
    m.insert("s0".into(), jf(s.s0));
    m.insert("vector".into(), jfl(&s.vector));
    m.insert("notes".into(), notes(&s.notes));
    (a, Value::Object(m))
}

/// Stage 7 files (no arrays).
#[must_use]
pub fn sigma_files(s: &SigmaStage) -> StageFiles {
    let mut m = Map::new();
    m.insert("cov_formal".into(), jm3(&s.cov_formal));
    m.insert("cov".into(), jm3(&s.cov));
    m.insert("sigma_q".into(), s.sigma_q.map_or(Value::Null, jf));
    m.insert("frag".into(), s.frag.map_or(Value::Null, jf));
    m.insert("rej".into(), s.rej.map_or(Value::Null, jf));
    m.insert("partial".into(), s.partial.map_or(Value::Null, jf));
    m.insert("fdist".into(), jf(s.fdist));
    m.insert("s0".into(), jf(s.s0));
    m.insert("hours".into(), jf(s.hours));
    m.insert("length_m".into(), jf(s.length_m));
    m.insert("n_rejected".into(), Value::from(s.n_rejected));
    m.insert("n_dd".into(), Value::from(s.n_dd));
    m.insert("notes".into(), notes(&s.notes));
    (BTreeMap::new(), Value::Object(m))
}

/// The stage files of every stage present, in order, with their stems (`01_observations` ...).
#[must_use]
pub fn stage_files(st: &Stages) -> Vec<(String, StageFiles)> {
    let mut out = Vec::new();
    let mut n = 0;
    let mut push = |name: &str, f: StageFiles| {
        n += 1;
        out.push((format!("{n:02}_{name}"), f));
    };
    if let Some(s) = &st.observations {
        push("observations", observations_files(s));
    }
    if let Some(s) = &st.arcs {
        push("arcs", arcs_files(s));
    }
    if let Some(s) = &st.float {
        push("float", float_files(s));
    }
    if let Some(s) = &st.wide_lane {
        push("wide_lane", wide_lane_files(s));
    }
    if let Some(s) = &st.narrow_lane {
        push("narrow_lane", narrow_lane_files(s));
    }
    if let Some(s) = &st.validation {
        push("validation", validation_files(s));
    }
    if let Some(s) = &st.sigma {
        push("sigma", sigma_files(s));
    }
    out
}

/// Writes every stage into `folder` (made if needed).
pub fn write_stages(folder: &Path, st: &Stages) -> std::io::Result<Vec<PathBuf>> {
    std::fs::create_dir_all(folder)?;
    let mut written = Vec::new();
    for (stem, (arrays, meta)) in stage_files(st) {
        if !arrays.is_empty() {
            let p = folder.join(format!("{stem}.npz"));
            std::fs::File::create(&p)?.write_all(&npz_bytes(&arrays))?;
            written.push(p);
        }
        let p = folder.join(format!("{stem}.json"));
        std::fs::File::create(&p)?.write_all(json_text(&meta).as_bytes())?;
        written.push(p);
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn npy_header_and_zip_shape() {
        let mut a = BTreeMap::new();
        a.insert(
            "x".to_owned(),
            NpyOut::F64(vec![2, 2], vec![1.0, 2.0, 3.0, 4.0]),
        );
        a.insert(
            "s".to_owned(),
            NpyOut::Str(vec![2], vec!["G05".into(), "R1".into()]),
        );
        let z = npz_bytes(&a);
        assert_eq!(&z[..4], &0x0403_4b50u32.to_le_bytes());
        let npy = npy_bytes(&NpyOut::Bool(vec![3], vec![true, false, true]));
        assert_eq!(&npy[..6], b"\x93NUMPY");
        let hlen = u16::from_le_bytes([npy[8], npy[9]]) as usize;
        assert_eq!((10 + hlen) % 64, 0);
        assert_eq!(&npy[10 + hlen..], &[1, 0, 1]);
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        let v = serde_json::json!({"b": [1.0, 2], "a": {"z": null}});
        assert_eq!(
            json_text(&v),
            "{\n \"a\": {\n  \"z\": null\n },\n \"b\": [\n  1.0,\n  2\n ]\n}\n"
        );
        assert_eq!(jround(26.25, 1), serde_json::json!(26.2));
        assert_eq!(jf(f64::NAN), Value::Null);
    }
}
