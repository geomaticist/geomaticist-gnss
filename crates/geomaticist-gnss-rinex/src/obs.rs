//! The engine's observation reader: a RINEX 2/3/4 file (plain, gzip or Hatanaka) into a
//! [`StationObs`].
//!
//! Rules: epochs with flags 0/1 and a readable time, inside `[t0 - 1 ms, t1 + 1 ms]` (reading stops at the
//! first epoch past `t1`), systems in `systems`, kinds C/L/S with a band; a value is the field's
//! first 14 characters as a float, blank, unreadable or exactly 0.0 values are missing; the LLI is
//! the 15th character when it is a digit. Epochs are sorted (stable) and duplicates keep the first.
//! The observation-type plan comes from the file's header, not from flag-4 changes (the epoch
//! reader itself follows the event).

use std::collections::HashMap;
use std::path::Path;

use geomaticist_gnss_core::obs::{Signal, norm_antenna};
use geomaticist_gnss_core::signal::{KINDS, SignalKey, code_key};
use geomaticist_gnss_core::text::{is_blank, latin1, py_float, sl};
use geomaticist_gnss_core::time::{TICKS, ticks_to_seconds};
use geomaticist_gnss_core::{ObsCode, Sat, StationObs};

use crate::RinexError;
use crate::header::RinexHeader;
use crate::reader::{EpochReader, read_header};

/// Reads an observation file between GPS seconds `t0` and `t1` (inclusive) for the systems in
/// `systems` (usually `"GRECJ"`).
pub fn read_obs(
    path: &Path,
    name: &str,
    t0: Option<f64>,
    t1: Option<f64>,
    systems: &str,
) -> Result<StationObs, RinexError> {
    let data = crate::open::read_plain(path)?;
    read_obs_bytes(path, &data, name, t0, t1, systems)
}

struct Row {
    idx: Vec<u32>,
    vals: Vec<f64>,
    llis: Vec<i8>,
    key: SignalKey,
    order: u32,
}

/// [`read_obs`] on plain RINEX bytes; `path` names the station when the header has no marker.
pub fn read_obs_bytes(
    path: &Path,
    data: &[u8],
    name: &str,
    t0: Option<f64>,
    t1: Option<f64>,
    systems: &str,
) -> Result<StationObs, RinexError> {
    read_obs_bytes_kinds(path, data, name, t0, t1, systems, KINDS)
}

/// [`read_obs`] keeping the observation kinds in `kinds` (`C` code, `L` phase, `S` signal strength,
/// `D` Doppler). [`read_obs`] keeps `"CLS"`; kinematic processing reads `"CLSD"` so that Doppler is
/// available for velocities and cycle-slip checks. A Doppler value is in Hz, positive when the range
/// decreases (the RINEX sign).
pub fn read_obs_kinds(
    path: &Path,
    name: &str,
    t0: Option<f64>,
    t1: Option<f64>,
    systems: &str,
    kinds: &str,
) -> Result<StationObs, RinexError> {
    let data = crate::open::read_plain(path)?;
    read_obs_bytes_kinds(path, &data, name, t0, t1, systems, kinds)
}

/// [`read_obs_kinds`] on plain RINEX bytes.
pub fn read_obs_bytes_kinds(
    path: &Path,
    data: &[u8],
    name: &str,
    t0: Option<f64>,
    t1: Option<f64>,
    systems: &str,
    kinds: &str,
) -> Result<StationObs, RinexError> {
    let (hdr, _eol) = read_header(data)?;
    let stem = path
        .file_stem()
        .map_or_default(|s| s.to_string_lossy().into_owned());
    let first4 = |s: &str| s.chars().take(4).collect::<String>().to_uppercase();
    let st_name = if !name.is_empty() {
        name.to_owned()
    } else if !hdr.marker_name().is_empty() {
        first4(&hdr.marker_name())
    } else {
        first4(&stem)
    };
    let mut st = StationObs::new(st_name, path);
    st.antenna = hdr
        .get("ANT # / TYPE")
        .map_or(String::new(), |c| norm_antenna(&latin1(sl(c, 20, 40))));
    st.receiver = hdr
        .get("REC # / TYPE / VERS")
        .map_or(String::new(), |c| latin1(sl(c, 20, 40)).trim().to_owned());
    st.approx = hdr.approx_xyz().unwrap_or([0.0; 3]);
    st.delta_hen = hdr.delta_hen().unwrap_or([0.0; 3]);
    let major = hdr.major();
    let ver = hdr.version;
    // (t0 - 1 ms) in ticks, rounded half to even
    let lo = t0.map(|t| ((t - 1e-3) * TICKS as f64).round_ties_even() as i64);
    let hi = t1.map(|t| ((t + 1e-3) * TICKS as f64).round_ties_even() as i64);

    // plan per system letter: (field index, code, key), from the file's header (`hdr`)
    let mut plans: HashMap<u8, Vec<(usize, ObsCode, SignalKey)>> = HashMap::new();
    let plan_for = |sysc: u8, hdr: &RinexHeader| -> Vec<(usize, ObsCode, SignalKey)> {
        let types = match hdr.obs_types.get(&sysc) {
            Some(t) if !t.is_empty() => t.as_slice(),
            _ => hdr.v2_types.as_slice(),
        };
        types
            .iter()
            .enumerate()
            .filter_map(|(i, c)| {
                let k = code_key(c.as_bytes(), char::from(sysc), major, ver);
                (kinds.contains(k.kind) && k.band != 0).then_some((i, *c, k))
            })
            .collect()
    };

    let mut ticks: Vec<i64> = Vec::new();
    let mut pf: Vec<bool> = Vec::new();
    let mut rows: HashMap<(Sat, ObsCode), Row> = HashMap::new();
    let mut next_order: u32 = 0;
    let rd = EpochReader::new(data, hdr.clone());
    for ep in rd {
        let ep = ep?;
        let Some(tick) = ep.tick else { continue };
        if ep.flag > 1 {
            continue;
        }
        if lo.is_some_and(|lo| tick < lo) {
            continue;
        }
        if hi.is_some_and(|hi| tick > hi) {
            break;
        }
        let ei = u32::try_from(ticks.len()).map_err(|_| RinexError::new("too many epochs"))?;
        ticks.push(tick);
        pf.push(ep.flag == 1);
        for (sat, flds) in ep.records() {
            let sysc = sat.sys.letter();
            if !systems.contains(sysc) {
                continue;
            }
            let plan = plans
                .entry(sysc as u8)
                .or_insert_with(|| plan_for(sysc as u8, &hdr));
            for (i, c, k) in plan.iter() {
                let Some(f) = flds.get(*i) else { continue };
                let v = &f[..14];
                if is_blank(v) {
                    continue;
                }
                let Some(x) = py_float(v) else { continue };
                if x == 0.0 {
                    continue;
                }
                let r = rows.entry((sat, *c)).or_insert_with(|| {
                    next_order += 1;
                    Row {
                        idx: Vec::new(),
                        vals: Vec::new(),
                        llis: Vec::new(),
                        key: *k,
                        order: next_order - 1,
                    }
                });
                r.idx.push(ei);
                r.vals.push(x);
                if k.kind == 'L' {
                    let ll = f[14];
                    r.llis.push(if ll.is_ascii_digit() {
                        (ll - b'0') as i8
                    } else {
                        0
                    });
                }
            }
        }
    }

    // sort (stable), duplicates keep the first
    let n = ticks.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| ticks[i]);
    let mut remap = vec![-1i64; n];
    let mut t_out: Vec<i64> = Vec::with_capacity(n);
    let mut pf_out: Vec<bool> = Vec::with_capacity(n);
    for (k, &i) in order.iter().enumerate() {
        if k > 0 && ticks[i] == ticks[order[k - 1]] {
            continue;
        }
        remap[i] = t_out.len() as i64;
        t_out.push(ticks[i]);
        pf_out.push(pf[i]);
    }
    let m = t_out.len();
    st.t = t_out.iter().map(|&t| ticks_to_seconds(t)).collect();
    st.power_fail = pf_out;
    for ((sat, c), r) in rows {
        let mut vals = vec![f64::NAN; m];
        let mut lli = (r.key.kind == 'L').then(|| vec![0i8; m]);
        for (j, &e) in r.idx.iter().enumerate() {
            let ii = remap[e as usize];
            if ii < 0 {
                continue;
            }
            let ii = ii as usize;
            // numpy fancy assignment: the last write to an index wins
            vals[ii] = r.vals[j];
            if let Some(l) = lli.as_mut() {
                l[ii] = r.llis[j];
            }
        }
        st.sats.entry(sat).or_default().insert(
            c,
            Signal {
                key: r.key,
                values: vals,
                lli,
                first_seen: r.order,
            },
        );
    }
    Ok(st)
}

#[cfg(test)]
mod tests {
    use super::*;

    const V2: &str =
        "     2.11           OBSERVATION DATA    G (GPS)             RINEX VERSION / TYPE
TEST                                                        MARKER NAME
                    TRM59800.99     SCIT                    ANT # / TYPE
     4    L1    L2    C1    P2                              # / TYPES OF OBSERV
                                                            END OF HEADER
 24  9 23  0  0 30.0000000  0  2G01G02
 100000000.1236   78000000.4565   20000000.000    20000004.000
 110000000.123    88000000.456    21000000.000           0.000
 24  9 23  0  0  0.0000000  0  1G01
  99999999.1231   77999999.456    19999999.000    20000003.000
 24  9 23  0  0 30.0000000  0  1G01
         1.000            1.000            1.000           1.000
";

    const V3D: &str =
        "     3.04           OBSERVATION DATA    G                   RINEX VERSION / TYPE
G    3 C1C L1C D1C                                          SYS / # / OBS TYPES
                                                            END OF HEADER
> 2024 09 23 00 00 30.0000000  0  1
G01  20000000.000   100000000.123       -1234.567
";

    #[test]
    fn reads_sorts_and_dedups() {
        let st = read_obs_bytes(
            Path::new("x/abcd2670.24o"),
            V2.as_bytes(),
            "",
            None,
            None,
            "GRECJ",
        )
        .unwrap();
        assert_eq!(st.name, "TEST");
        assert_eq!(st.antenna, "TRM59800.99     SCIT");
        assert_eq!(st.t.len(), 2);
        let g1: Sat = "G01".parse().unwrap();
        let g2: Sat = "G02".parse().unwrap();
        let l1: ObsCode = "L1".parse().unwrap();
        let p2: ObsCode = "P2".parse().unwrap();
        let s = &st.sats[&g1][&l1];
        assert_eq!(s.values, vec![99_999_999.123, 100_000_000.123]); // duplicate: first kept
        assert_eq!(s.lli.as_ref().unwrap(), &vec![1, 6]);
        assert!(!st.sats[&g2].contains_key(&p2)); // 0.0 is missing
        assert_eq!(st.sats[&g1][&p2].key.attr, Some('P'));
        let d = read_obs_bytes_kinds(
            Path::new("x.o"),
            V3D.as_bytes(),
            "N",
            None,
            None,
            "G",
            "CLSD",
        )
        .unwrap();
        let d1: ObsCode = "D1C".parse().unwrap();
        assert_eq!(d.sats[&g1][&d1].values, vec![-1234.567]);
        assert_eq!(d.sats[&g1][&d1].key.kind, 'D');
        assert!(d.sats[&g1][&d1].lli.is_none());
        let nod = read_obs_bytes(Path::new("x.o"), V3D.as_bytes(), "N", None, None, "G").unwrap();
        assert!(!nod.sats[&g1].contains_key(&d1));
        let w = read_obs_bytes(
            Path::new("x.o"),
            V2.as_bytes(),
            "N",
            Some(1_411_084_830.0),
            None,
            "G",
        )
        .unwrap();
        assert_eq!(w.t, vec![1_411_084_830.0]);
        let none = read_obs_bytes(Path::new("x.o"), V2.as_bytes(), "N", None, None, "E").unwrap();
        assert!(none.sats.is_empty());
    }
}
