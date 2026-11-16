//! The PPP stage files (`NN_stage.npz` of sorted `.npy` members, `NN_stage.json` with sorted keys,
//! indent 1, NaN as null), built with the engine's writers.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use geomaticist_gnss_engine::arr::{B2, F2, I2, Vec3Arr};
use geomaticist_gnss_engine::dump::{
    NpyOut, StageFiles, jbl, jf, jfl, jil, jm3, jsl, json_text, npz_bytes,
};
use geomaticist_gnss_engine::linalg::DMat;
use serde_json::{Map, Value};

use crate::options::antenna_json;
use crate::stages::{
    ArcsStage, FloatStage, NarrowLaneStage, ObservationsStage, SigmaStage, Stages, ValidationStage,
    WideLaneStage,
};

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
fn chars(v: &[char]) -> Value {
    Value::Array(v.iter().map(|c| Value::String(c.to_string())).collect())
}

/// Stage 1 files.
#[must_use]
pub fn observations_files(s: &ObservationsStage) -> StageFiles {
    let mut a = BTreeMap::new();
    let mut put = |k: &str, v: NpyOut| {
        a.insert(k.to_owned(), v);
    };
    put("t", f1(&s.t));
    put("t_obs", f1(&s.t_obs));
    put(
        "cols",
        NpyOut::Str(
            vec![s.cols.len()],
            s.cols.iter().map(ToString::to_string).collect(),
        ),
    );
    put("phase", f2(&s.phase));
    put("code", f2(&s.code));
    put("gf", f2(&s.gf));
    put("mw", f2(&s.mw));
    put("lli", b2(&s.lli));
    put("sxyz", v3(&s.sxyz));
    put("sclk", f2(&s.sclk));
    put("lamif", f2(&s.lamif));
    put("ok", b2(&s.ok));
    put("ok_raw", b2(&s.ok_raw));
    put("shadow", b2(&s.shadow));
    put("pcvs", f2(&s.pcvs));
    put("az", f2(&s.az));
    put("el", f2(&s.el));
    put("windup", f2(&s.windup));
    put("disp", rows3(&s.disp));
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
    let mut m = Map::new();
    m.insert("station".into(), Value::String(s.station.clone()));
    m.insert("x0".into(), jfl(&s.x0));
    m.insert("pairs".into(), Value::Object(pairs));
    m.insert("signals".into(), jsl(&s.signals));
    m.insert("a_if".into(), jfl(&s.a_if));
    m.insert("b_if".into(), jfl(&s.b_if));
    m.insert("lam1".into(), jfl(&s.lam1));
    m.insert("lam2".into(), jfl(&s.lam2));
    m.insert("osb".into(), jbl(&s.osb));
    m.insert("osb_hits".into(), Value::from(s.osb_hits));
    m.insert("osb_miss".into(), Value::from(s.osb_miss));
    m.insert("notes".into(), jsl(&s.notes));
    (a, Value::Object(m))
}

/// Stage 2 files.
#[must_use]
pub fn arcs_files(s: &ArcsStage) -> StageFiles {
    let mut a = BTreeMap::new();
    a.insert("arc".to_owned(), i2(&s.arc));
    a.insert("amb0".to_owned(), f2(&s.amb0));
    a.insert("ok".to_owned(), b2(&s.ok));
    a.insert("phase_ok".to_owned(), b2(&s.phase_ok));
    a.insert("code_ok".to_owned(), b2(&s.code_ok));
    let mut m = Map::new();
    m.insert("n_arcs".into(), Value::from(s.n_arcs));
    m.insert("antenna".into(), antenna_json(&s.antenna));
    m.insert("off_sys".into(), chars(&s.off_sys));
    m.insert("notes".into(), jsl(&s.notes));
    (a, Value::Object(m))
}

/// Stage 3 files.
#[must_use]
pub fn float_files(s: &FloatStage) -> StageFiles {
    let mut a = BTreeMap::new();
    let mut put = |k: &str, v: NpyOut| {
        a.insert(k.to_owned(), v);
    };
    put("x", f1(&s.x));
    put("N", mat(&s.n));
    put("Ninv", mat(&s.ninv));
    put("sol", f1(&s.sol));
    put("resp", f2(&s.resp));
    put("resc", f2(&s.resc));
    put("clk", f1(&s.clk));
    put("arc", i2(&s.arc));
    put("use_p", b2(&s.use_p));
    put("use_c", b2(&s.use_c));
    put("rej_phase", b2(&s.rej_phase));
    put("rej_code", b2(&s.rej_code));
    put("amb0", f2(&s.amb0));
    let mut m = Map::new();
    m.insert("sigma0".into(), jf(s.sigma0));
    m.insert("cov".into(), jm3(&s.cov));
    m.insert("npar".into(), Value::from(s.npar));
    m.insert("iI".into(), Value::from(s.i_i));
    m.insert("iG".into(), Value::from(s.i_g));
    m.insert("iO".into(), Value::from(s.i_o));
    m.insert("iZ".into(), Value::from(s.i_z));
    m.insert("iT".into(), Value::from(s.i_t));
    m.insert("iA".into(), Value::from(s.i_a));
    m.insert("nz".into(), Value::from(s.nz));
    m.insert("ng".into(), Value::from(s.ng));
    m.insert("namb".into(), Value::from(s.namb));
    m.insert("isb_sys".into(), chars(&s.isb_sys));
    m.insert("ref_sys".into(), Value::String(s.ref_sys.to_string()));
    m.insert("glo_sats".into(), jil(&s.glo_sats));
    m.insert("off_sys".into(), chars(&s.off_sys));
    m.insert("n_rej".into(), Value::from(s.n_rej));
    m.insert("n_split".into(), Value::from(s.n_split));
    m.insert("nobs".into(), Value::from(s.nobs));
    m.insert("notes".into(), jsl(&s.notes));
    (a, Value::Object(m))
}

/// Stage 4 files (JSON only).
#[must_use]
pub fn wide_lane_files(s: &WideLaneStage) -> StageFiles {
    let mut m = Map::new();
    m.insert("ar".into(), Value::Bool(s.ar));
    m.insert(
        "wl".into(),
        s.wl.map_or(Value::Null, |w| jil(&[w[0], w[1]])),
    );
    m.insert(
        "why".into(),
        s.why.clone().map_or(Value::Null, Value::String),
    );
    m.insert("ok".into(), s.ok.map_or(Value::Null, Value::Bool));
    (BTreeMap::new(), Value::Object(m))
}

/// Stage 5 files (JSON only).
#[must_use]
pub fn narrow_lane_files(s: &NarrowLaneStage) -> StageFiles {
    (
        BTreeMap::new(),
        geomaticist_gnss_engine::dump::stage_nan(&Value::Object(s.ar.clone())),
    )
}

/// Stage 6 files.
#[must_use]
pub fn validation_files(s: &ValidationStage) -> StageFiles {
    let mut a = BTreeMap::new();
    a.insert("x".to_owned(), f1(&s.x));
    a.insert(
        "cov".to_owned(),
        NpyOut::F64(vec![3, 3], s.cov.iter().flatten().copied().collect()),
    );
    let mut m = Map::new();
    m.insert("fixed".into(), Value::Bool(s.fixed));
    m.insert("sigma0".into(), jf(s.sigma0));
    m.insert("n_rejected".into(), Value::from(s.n_rejected));
    m.insert("n_split".into(), Value::from(s.n_split));
    m.insert("notes".into(), jsl(&s.notes));
    (a, Value::Object(m))
}

/// Stage 7 files (JSON only).
#[must_use]
pub fn sigma_files(s: &SigmaStage) -> StageFiles {
    let mut m = Map::new();
    m.insert("cov_formal".into(), jm3(&s.cov_formal));
    m.insert("cov_enu".into(), jm3(&s.cov_enu));
    m.insert("cov_realistic_enu".into(), jm3(&s.cov_realistic_enu));
    m.insert("sigma0".into(), jf(s.sigma0));
    m.insert("hours".into(), jf(s.hours));
    m.insert("rms_phase".into(), jf(s.rms_phase));
    m.insert("rms_code".into(), jf(s.rms_code));
    m.insert("n_obs_phase".into(), Value::from(s.n_obs_phase));
    m.insert("n_obs_code".into(), Value::from(s.n_obs_code));
    m.insert("n_amb".into(), Value::from(s.n_amb));
    m.insert("n_rejected".into(), Value::from(s.n_rejected));
    m.insert("notes".into(), jsl(&s.notes));
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
    fn wide_lane_json_has_the_expected_keys() {
        let s = WideLaneStage {
            ar: false,
            wl: None,
            why: None,
            ok: None,
        };
        let (a, m) = wide_lane_files(&s);
        assert!(a.is_empty());
        assert_eq!(
            json_text(&m),
            "{\n \"ar\": false,\n \"ok\": null,\n \"why\": null,\n \"wl\": null\n}\n"
        );
        let st = Stages {
            wide_lane: Some(s),
            ..Stages::default()
        };
        assert_eq!(stage_files(&st)[0].0, "01_wide_lane");
    }
}
