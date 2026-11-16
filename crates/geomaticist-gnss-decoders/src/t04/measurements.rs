//! rec35 sub19 "Measurement Record" -> per-epoch observations.
//!
//! Keyframes (kind 1) carry absolute values. Delta records (kind 2) carry differences, applied element by
//! element by array index against the previous epoch's arrays (missing old elements count as 0; an empty delta
//! array means all-zero deltas).
//!
//! Each field has a delta ORDER: order 1: `x[n] = x[n-1] + d[n]`; order 2: `v[n] = v[n-1] + d[n]`,
//! `x[n] = x[n-1] + v[n]` (v reset to 0 at keyframes). Orders are listed in [`delta_order`]; the others are 1.

use std::rc::Rc;

use super::avro::{self, Flat};
use super::container::{GAP, Inner};
use super::rec35::{KIND_DELTA, KIND_FULL, Messages, Rec35Error};
use super::signals;
use super::{SatId, T04Error};

/// The Measurement Record's subtype.
pub const SUBTYPE: i64 = 19;

/// clkOffset units per millisecond (2^19).
pub const MS: i64 = 524_288;

/// Delta order of a field. Confirmed against RINEX: phase (2), range (1); receiveTime (2) by epoch times.
#[must_use]
pub fn delta_order(name: &str) -> u8 {
    match name {
        "epoch__receiveTime"
        | "sv__meas__phase"
        | "clock__clkOffset"
        | "clock__steerClkOffsetM"
        | "clock__system__systemClockOffset" => 2,
        _ => 1,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Group {
    Meas,
    ClockSystem,
    Sv,
    None,
}

fn group_of(name: &str) -> Group {
    if name.starts_with("sv__meas__") {
        Group::Meas
    } else if name.starts_with("clock__system__") {
        Group::ClockSystem
    } else if name.starts_with("sv__") || name == "_Nb_Na" {
        Group::Sv
    } else {
        Group::None
    }
}

/// Indexes of the fields the decoder reads (by name in the file's schema).
#[derive(Debug, Clone)]
pub struct Fields {
    names: Vec<Rc<str>>,
    receive_time: Option<usize>,
    week: Option<usize>,
    clk_offset: Option<usize>,
    na: Option<usize>,
    nb: Option<usize>,
    nb_na: Option<usize>,
    sv_id: Option<usize>,
    sattype: Option<usize>,
    fdma: Option<usize>,
    block: Option<usize>,
    track: Option<usize>,
    flags: Option<usize>,
    snr: Option<usize>,
    doppler: Option<usize>,
    phase: Option<usize>,
    range: Option<usize>,
}

impl Fields {
    fn new(layout: &[(Rc<str>, bool)]) -> Self {
        let names: Vec<Rc<str>> = layout.iter().map(|l| l.0.clone()).collect();
        let f = |n: &str| names.iter().position(|x| &**x == n);
        Self {
            receive_time: f("epoch__receiveTime"),
            week: f("epoch__weekNum"),
            clk_offset: f("clock__clkOffset"),
            na: f("_Na"),
            nb: f("_Nb"),
            nb_na: f("_Nb_Na"),
            sv_id: f("sv__svId"),
            sattype: f("sv__rt_sattype"),
            fdma: f("sv__fdmaChannel"),
            block: f("sv__meas__blockType"),
            track: f("sv__meas__trackType"),
            flags: f("sv__meas__measFlags"),
            snr: f("sv__meas__snr"),
            doppler: f("sv__meas__doppler"),
            phase: f("sv__meas__phase"),
            range: f("sv__meas__range"),
            names,
        }
    }

    fn key_err(&self, i: Option<usize>, what: &str) -> T04Error {
        let _ = (self, i);
        T04Error::Decode(format!("KeyError: '{what}'"))
    }
}

fn scalar(r: &[Flat], i: Option<usize>, name: &str) -> Result<i64, T04Error> {
    match i.and_then(|i| r.get(i)) {
        Some(Flat::Scalar(v)) => Ok(*v),
        Some(Flat::Array(_)) => Err(T04Error::Decode(format!("TypeError: '{name}' is a list"))),
        None => Err(T04Error::Decode(format!("KeyError: '{name}'"))),
    }
}

fn array<'a>(r: &'a [Flat], i: Option<usize>, name: &str) -> Result<&'a [i64], T04Error> {
    match i.and_then(|i| r.get(i)) {
        Some(Flat::Array(v)) => Ok(v),
        Some(Flat::Scalar(_)) => Err(T04Error::Decode(format!("TypeError: '{name}' is an int"))),
        None => Err(T04Error::Decode(format!("KeyError: '{name}'"))),
    }
}

fn at(a: &[i64], j: usize) -> Result<i64, T04Error> {
    a.get(j)
        .copied()
        .ok_or_else(|| T04Error::Decode("IndexError: list index out of range".into()))
}

/// One raw epoch record (absolute values), with the schema's field indexes.
#[derive(Debug)]
pub struct RawEpoch<'a> {
    /// Field values in schema order.
    pub vals: &'a [Flat],
    /// The field indexes.
    pub fields: &'a Fields,
}

impl RawEpoch<'_> {
    /// `epoch__weekNum`.
    pub fn week(&self) -> Result<i64, T04Error> {
        scalar(self.vals, self.fields.week, "epoch__weekNum")
    }
    /// `epoch__receiveTime` (ms of week).
    pub fn receive_time(&self) -> Result<i64, T04Error> {
        scalar(self.vals, self.fields.receive_time, "epoch__receiveTime")
    }
    /// `clock__clkOffset` (2^-19 ms).
    pub fn clk_offset(&self) -> Result<i64, T04Error> {
        scalar(self.vals, self.fields.clk_offset, "clock__clkOffset")
    }
    /// `sv__svId`.
    pub fn sv_id(&self) -> Result<&[i64], T04Error> {
        array(self.vals, self.fields.sv_id, "sv__svId")
    }
    /// `sv__rt_sattype`.
    pub fn sattype(&self) -> Result<&[i64], T04Error> {
        array(self.vals, self.fields.sattype, "sv__rt_sattype")
    }
    /// `sv__fdmaChannel`.
    pub fn fdma(&self) -> Result<&[i64], T04Error> {
        array(self.vals, self.fields.fdma, "sv__fdmaChannel")
    }
    fn arr(&self, i: Option<usize>, name: &str) -> Result<&[i64], T04Error> {
        array(self.vals, i, name)
    }
}

struct State {
    x: Vec<Flat>,
    v: Vec<Flat>,
    absent: Vec<bool>,
    have: bool,
}

fn pad_to(v: &mut Vec<i64>, n: usize) {
    v.truncate(n);
    v.resize(n, 0);
}

/// Calls `f(kind, epoch)` for every measurement epoch with absolute values.
/// A damaged spot drops the state: the next epoch used is the next keyframe.
pub fn for_each_raw_epoch<F>(inner: &Inner, mut f: F) -> Result<(), T04Error>
where
    F: FnMut(i32, &RawEpoch<'_>) -> Result<(), T04Error>,
{
    let schemas = super::rec35::schemas(inner)?;
    let sch = schemas.get(&SUBTYPE).ok_or_else(|| {
        T04Error::Decode("no Measurement Record schema (rec35 sub19) in file".into())
    })?;
    let layout = avro::flat_layout(sch).ok_or_else(|| {
        T04Error::Decode(
            "unsupported Measurement Record schema (fields not long / array of long)".into(),
        )
    })?;
    let fields = Fields::new(&layout);
    let order: Vec<u8> = layout.iter().map(|l| delta_order(&l.0)).collect();
    let groups: Vec<Group> = layout.iter().map(|l| group_of(&l.0)).collect();
    let mut st = State {
        x: Vec::new(),
        v: Vec::new(),
        absent: vec![false; layout.len()],
        have: false,
    };
    let mut d: Vec<Flat> = Vec::new();
    for m in Messages::new(inner, |s| s == SUBTYPE) {
        let m = m?;
        if m.subtype == i64::from(GAP) {
            st.have = false; // damaged data: wait for the next keyframe
            continue;
        }
        if m.subtype != SUBTYPE || (m.kind != KIND_FULL && m.kind != KIND_DELTA) {
            continue;
        }
        avro::decode_flat(&layout, &m.body, &mut d).map_err(Rec35Error::from)?;
        if m.kind == KIND_FULL {
            st.x.clone_from(&d);
            st.v = d
                .iter()
                .map(|v| match v {
                    Flat::Scalar(_) => Flat::Scalar(0),
                    Flat::Array(a) => Flat::Array(vec![0; a.len()]),
                })
                .collect();
            // a field this keyframe leaves empty while its group has members isn't logged and stays empty; one
            // empty only because its group is (fw 5.41 just after power-on, tracking nothing) fills in
            let nb0 = scalar(&d, fields.nb, "_Nb")?;
            let na0 = scalar(&d, fields.na, "_Na")?;
            let nmeas0: i64 = array(&d, fields.nb_na, "_Nb_Na")?.iter().sum();
            for (k, val) in d.iter().enumerate() {
                let g = match groups[k] {
                    Group::Meas => nmeas0,
                    Group::ClockSystem => na0,
                    Group::Sv => nb0,
                    Group::None => 0,
                };
                st.absent[k] = matches!(val, Flat::Array(a) if a.is_empty()) && g > 0;
            }
            st.have = true;
            f(
                m.kind,
                &RawEpoch {
                    vals: &st.x,
                    fields: &fields,
                },
            )?;
            continue;
        }
        if !st.have {
            continue; // delta before the first keyframe: cannot resolve
        }
        // array group sizes for this epoch (schema field order puts the counts first)
        let nb = scalar(&st.x, fields.nb, "_Nb")? + scalar(&d, fields.nb, "_Nb")?;
        let na = scalar(&st.x, fields.na, "_Na")? + scalar(&d, fields.na, "_Na")?;
        let nbu = usize::try_from(nb).unwrap_or(0);
        let xa = array(&st.x, fields.nb_na, "_Nb_Na")?;
        let da = array(&d, fields.nb_na, "_Nb_Na")?;
        let nmeas: i64 = (0..nbu)
            .map(|i| xa.get(i).copied().unwrap_or(0) + da.get(i).copied().unwrap_or(0))
            .sum();
        for k in 0..d.len() {
            let o = order[k];
            match (&d[k], &mut st.x[k], &mut st.v[k]) {
                (Flat::Array(dd), Flat::Array(xv), Flat::Array(vv)) => {
                    let g = match groups[k] {
                        Group::Meas => Some(nmeas),
                        Group::ClockSystem => Some(na),
                        Group::Sv => Some(nb),
                        Group::None => None,
                    };
                    let n = if !dd.is_empty() {
                        st.absent[k] = false;
                        dd.len()
                    } else if let Some(g) = g {
                        if st.absent[k] {
                            0
                        } else {
                            usize::try_from(g).unwrap_or(0)
                        }
                    } else {
                        xv.len()
                    };
                    pad_to(xv, n);
                    if o == 2 {
                        pad_to(vv, n);
                        for i in 0..n {
                            vv[i] += dd.get(i).copied().unwrap_or(0);
                            xv[i] += vv[i];
                        }
                    } else {
                        for (x, dv) in xv.iter_mut().zip(dd.iter()) {
                            *x += dv;
                        }
                    }
                }
                (Flat::Scalar(dd), Flat::Scalar(xv), Flat::Scalar(vv)) => {
                    if o == 2 {
                        *vv += dd;
                        *xv += *vv;
                    } else {
                        *xv += dd;
                    }
                }
                _ => {
                    return Err(T04Error::Decode(format!(
                        "TypeError: field {} changed type",
                        fields.names[k]
                    )));
                }
            }
        }
        f(
            m.kind,
            &RawEpoch {
                vals: &st.x,
                fields: &fields,
            },
        )?;
    }
    Ok(())
}

/// One decoded measurement block.
#[derive(Debug, Clone)]
pub struct Obs {
    /// Satellite.
    pub sat: SatId,
    /// RINEX code suffix ("1C"), `None` for a signal not in the table.
    pub code: Option<&'static str>,
    /// Pseudorange, m (`None` for an unknown signal).
    pub c: Option<f64>,
    /// Raw phase (cycles = -phase / 32768).
    pub phase: i64,
    /// Raw Doppler (Hz = doppler / 256), if logged.
    pub doppler: Option<i64>,
    /// Raw SNR (dB-Hz = snr / 10), if logged.
    pub snr: Option<i64>,
    /// measFlags, if logged.
    pub flags: Option<i64>,
}

impl Obs {
    /// measFlags or 0 (`raw.get("sv__meas__measFlags", 0)`).
    #[must_use]
    pub fn flags(&self) -> i64 {
        self.flags.unwrap_or(0)
    }

    /// SNR in dB-Hz.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn s(&self) -> f64 {
        self.snr.unwrap_or(0) as f64 / 10.0
    }
}

/// The observations of a raw epoch, into `out`.
#[allow(clippy::cast_precision_loss)]
pub fn epoch_obs(r: &RawEpoch<'_>, out: &mut Vec<Obs>) -> Result<(), T04Error> {
    out.clear();
    let f = r.fields;
    let sv_id = r.sv_id()?;
    let sattype = r.sattype()?;
    let fdma = r.fdma()?;
    let nb_na = r.arr(f.nb_na, "_Nb_Na")?;
    let block = r.arr(f.block, "sv__meas__blockType")?;
    let track = r.arr(f.track, "sv__meas__trackType")?;
    let phase = r.arr(f.phase, "sv__meas__phase")?;
    let range = r.arr(f.range, "sv__meas__range")?;
    let snr = r.arr(f.snr, "sv__meas__snr").ok();
    let doppler = r.arr(f.doppler, "sv__meas__doppler").ok();
    let flags = r.arr(f.flags, "sv__meas__measFlags").ok();
    if snr.is_none() || doppler.is_none() || flags.is_none() {
        // the measurements are built from these names: a schema without them cannot be decoded
        return Err(f.key_err(f.snr, "sv__meas__snr"));
    }
    let mut j = 0usize;
    for (i, &sv) in sv_id.iter().enumerate() {
        let st = at(sattype, i)?;
        let k = fdma.get(i).copied().unwrap_or(0);
        let (sysl, off) = signals::sattype(st);
        let sat = SatId::new(sysl, sv + off);
        for _ in 0..at(nb_na, i)?.max(0) {
            let b = at(block, j)?;
            let t = at(track, j)?;
            let ph = at(phase, j)?;
            let rng = at(range, j)?;
            let get = |a: Option<&[i64]>| a.and_then(|a| a.get(j).copied());
            let mut o = Obs {
                sat,
                code: None,
                c: None,
                phase: ph,
                doppler: get(doppler),
                snr: get(snr),
                flags: get(flags),
            };
            j += 1;
            if let Some((code, freq)) = signals::signal(st, b, t) {
                let lam = signals::wavelength(freq, k);
                o.code = Some(code);
                o.c = Some((-ph as f64 / 32768.0) * lam + rng as f64 / 128.0);
            }
            out.push(o);
        }
    }
    Ok(())
}

/// `(sat, code, measFlags)` of every measurement block of a raw epoch, in [`epoch_obs`]'s order, code `None`
/// for a signal not in the table.
pub fn epoch_signals(
    r: &RawEpoch<'_>,
    out: &mut Vec<(SatId, Option<&'static str>, i64)>,
) -> Result<(), T04Error> {
    out.clear();
    let f = r.fields;
    let sv_id = r.sv_id()?;
    let sattype = r.sattype()?;
    let nb_na = r.arr(f.nb_na, "_Nb_Na")?;
    let block = r.arr(f.block, "sv__meas__blockType")?;
    let track = r.arr(f.track, "sv__meas__trackType")?;
    let flags = r.arr(f.flags, "sv__meas__measFlags")?;
    let mut j = 0usize;
    for (i, &sv) in sv_id.iter().enumerate() {
        let st = at(sattype, i)?;
        let (sysl, off) = signals::sattype(st);
        let sat = SatId::new(sysl, sv + off);
        for _ in 0..at(nb_na, i)?.max(0) {
            let sig = signals::signal(st, at(block, j)?, at(track, j)?);
            out.push((sat, sig.map(|s| s.0), flags.get(j).copied().unwrap_or(0)));
            j += 1;
        }
    }
    Ok(())
}

/// A decoded epoch.
#[derive(Debug, Clone, Default)]
pub struct Epoch {
    /// GPS week.
    pub week: i64,
    /// Time of week, ms (continuous-clock corrected when asked).
    pub tow_ms: i64,
    /// 1 keyframe, 2 delta.
    pub kind: i32,
    /// Observations.
    pub obs: Vec<Obs>,
    /// clkOffset, units of 2^-19 ms.
    pub clk_raw: i64,
    /// Receiver time before any continuous-clock correction.
    pub tow_raw_ms: i64,
}

/// Calls `f(epoch)` for every decoded epoch.
///
/// `continuous_clock` is the "continuous clock (fix msec steps)" option: each 1 ms
/// receiver clock step seen in `clock__clkOffset` is accumulated (n ms), and the epoch time and every
/// pseudorange are shifted by -n ms (-n * c * 1e-3 m). Phase, Doppler, SNR and the clock offset are unchanged.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
pub fn for_each_epoch<F>(inner: &Inner, continuous_clock: bool, mut f: F) -> Result<(), T04Error>
where
    F: FnMut(&Epoch) -> Result<(), T04Error>,
{
    let mut n_ms: i64 = 0;
    let mut prev: Option<i64> = None;
    let mut ep = Epoch::default();
    for_each_raw_epoch(inner, |kind, r| {
        let clk = r.clk_offset()?;
        if let Some(p) = prev {
            n_ms += ((clk - p) as f64 / MS as f64).round_ties_even() as i64;
        }
        prev = Some(clk);
        epoch_obs(r, &mut ep.obs)?;
        let mut tow = r.receive_time()?;
        if continuous_clock && n_ms != 0 {
            tow -= n_ms;
            let shift = n_ms as f64 * signals::C_LIGHT * 1e-3;
            for o in &mut ep.obs {
                if let Some(c) = o.c.as_mut() {
                    *c -= shift;
                }
            }
        }
        ep.week = r.week()?;
        ep.tow_ms = tow;
        ep.kind = kind;
        ep.clk_raw = clk;
        ep.tow_raw_ms = r.receive_time()?;
        f(&ep)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orders_and_groups() {
        assert_eq!(delta_order("sv__meas__phase"), 2);
        assert_eq!(delta_order("sv__meas__range"), 1);
        assert_eq!(group_of("sv__meas__snr"), Group::Meas);
        assert_eq!(group_of("_Nb_Na"), Group::Sv);
        assert_eq!(group_of("clock__system__rt_clktype"), Group::ClockSystem);
        assert_eq!(group_of("_Nb"), Group::None);
    }

    #[test]
    fn padding() {
        let mut v = vec![1, 2, 3];
        pad_to(&mut v, 5);
        assert_eq!(v, vec![1, 2, 3, 0, 0]);
        pad_to(&mut v, 2);
        assert_eq!(v, vec![1, 2]);
    }
}
