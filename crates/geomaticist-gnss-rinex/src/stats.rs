//! Choosing epochs (time window and new interval) and describing a file or a selection.
//!
//! Decimation keeps, for every multiple of
//! the new interval counted from GPS time zero, the one epoch nearest to it within half the source
//! interval. Nothing is averaged or interpolated.

use std::collections::BTreeMap;

use geomaticist_gnss_core::system::SYSTEM_ORDER;
use geomaticist_gnss_core::text::{format_g, round_to};
use geomaticist_gnss_core::time::{TICKS, fmt, fmt_duration, from_ticks, ordinal};
use serde_json::{Map, Value, json};

use crate::RinexError;
use crate::header::{RinexHeader, xyz_to_llh};
use crate::obsmap::system_name;
use crate::reader::RinexObs;

/// Sampling rates OPUS accepts (s).
pub const OPUS_RATES: [f64; 7] = [1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 30.0];

const FAR: i64 = 1 << 62;

/// The decimation interval actually applied: `None` when it would drop nothing.
#[must_use]
pub fn effective_interval(interval: Option<i64>, sample: i64) -> Option<i64> {
    interval.filter(|&i| i != 0 && i > sample)
}

/// Fails when a longer interval is not a multiple of the data's.
pub fn check_interval(interval: Option<i64>, sample: i64) -> Result<(), RinexError> {
    if let Some(i) = interval
        && i != 0
        && i > sample
        && i % sample != 0
    {
        return Err(RinexError::new(format!(
            "{} s isn't a multiple of the data's {} s interval - the epochs would be unevenly \
             spaced. Choose a multiple.",
            format_g(i as f64 / 1e7),
            format_g(sample as f64 / 1e7)
        )));
    }
    Ok(())
}

fn is_sorted(t: &[i64]) -> bool {
    t.windows(2).all(|w| w[0] <= w[1])
}

/// `(t + interval // 2) // interval` and the distance to that multiple.
fn nearest(t: i64, interval: i64) -> (i64, i64) {
    let k = (t + interval.div_euclid(2)).div_euclid(interval);
    (k, (t - k * interval).abs())
}

/// Indexes (into `ticks`) of the epochs to keep, in file order; `start` / `end` inclusive.
/// Files with out-of-order epochs take the slow path.
#[must_use]
pub fn select_epochs(
    ticks: &[i64],
    start: Option<i64>,
    end: Option<i64>,
    interval: Option<i64>,
    sample: i64,
) -> Vec<usize> {
    let decimate = interval.filter(|&i| i != 0 && i > sample);
    if !is_sorted(ticks) {
        let lo_t = start.unwrap_or(-FAR);
        let hi_t = end.unwrap_or(FAR);
        let inside: Vec<usize> = (0..ticks.len())
            .filter(|&i| lo_t <= ticks[i] && ticks[i] <= hi_t)
            .collect();
        let Some(iv) = decimate else { return inside };
        let tol = (sample.min(iv) / 2).max(1);
        let mut best: BTreeMap<i64, (usize, i64)> = BTreeMap::new();
        for i in inside {
            let (k, d) = nearest(ticks[i], iv);
            if d <= tol && best.get(&k).is_none_or(|b| d < b.1) {
                best.insert(k, (i, d));
            }
        }
        let mut out: Vec<usize> = best.values().map(|v| v.0).collect();
        out.sort_unstable();
        return out;
    }
    let lo = start.map_or(0, |s| ticks.partition_point(|&t| t < s));
    let hi = end.map_or(ticks.len(), |e| ticks.partition_point(|&t| t <= e));
    let Some(iv) = decimate else {
        return (lo..hi.max(lo)).collect();
    };
    let tol = (sample.min(iv) / 2).max(1);
    let mut keep = Vec::new();
    let mut cur: Option<(i64, usize, i64)> = None; // (k, best index, best distance)
    for (i, &t) in ticks.iter().enumerate().take(hi).skip(lo) {
        let (k, d) = nearest(t, iv);
        if d > tol {
            continue;
        }
        match cur {
            Some((ck, _, bd)) if ck == k => {
                if d < bd {
                    cur = Some((k, i, d));
                }
            }
            _ => {
                if let Some((_, bi, _)) = cur {
                    keep.push(bi);
                }
                cur = Some((k, i, d));
            }
        }
    }
    if let Some((_, bi, _)) = cur {
        keep.push(bi);
    }
    keep
}

/// Gaps longer than 1.5 sampling intervals (at most `limit`).
#[must_use]
pub fn gaps(ticks: &[i64], sample: i64, limit: usize) -> Vec<(i64, i64)> {
    let lim = sample * 3 / 2;
    let mut out = Vec::new();
    for w in ticks.windows(2) {
        if w[1] - w[0] > lim {
            out.push((w[0], w[1]));
            if out.len() >= limit {
                break;
            }
        }
    }
    out
}

fn opt_xyz(v: Option<[f64; 3]>) -> Value {
    v.map_or(Value::Null, |x| json!(x))
}

/// `RinexHeader.summary()`.
#[must_use]
pub fn header_summary(h: &RinexHeader) -> Value {
    let rec: Vec<String> = [h.rec_type(), h.rec_number()]
        .into_iter()
        .filter(|x| !x.is_empty())
        .collect();
    let types: Map<String, Value> = h
        .obs_types
        .iter()
        .map(|(s, c)| {
            (
                char::from(*s).to_string(),
                json!(
                    c.iter()
                        .map(geomaticist_gnss_core::ObsCode::as_str)
                        .collect::<Vec<_>>()
                ),
            )
        })
        .collect();
    json!({
        "version": h.version_str, "system": char::from(h.system).to_string(),
        "marker_name": h.marker_name(), "marker_number": h.marker_number(),
        "marker_type": h.marker_type(), "observer": h.observer(), "agency": h.agency(),
        "receiver": rec.join(" "), "receiver_version": h.rec_version(),
        "antenna_type": h.ant_type(), "antenna_radome": h.ant_radome(),
        "antenna_number": h.ant_number(), "approx_xyz": opt_xyz(h.approx_xyz()),
        "delta_hen": opt_xyz(h.delta_hen()), "interval": h.interval(),
        "program": h.program(), "obs_types": types,
    })
}

#[derive(Default)]
struct SysSum {
    sats: Vec<String>,
    epochs: usize,
    slips: usize,
    obs: BTreeMap<String, usize>,
    ssi: BTreeMap<String, (i64, usize)>,
    snr: BTreeMap<String, (f64, usize)>,
}

/// Everything `rinex info` and `.stats.json` show, for the whole file or a selection;
/// `header` is the header as edited (default the file's own).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn summarize(
    obs: &RinexObs,
    start: Option<i64>,
    end: Option<i64>,
    interval: Option<i64>,
    header: Option<&RinexHeader>,
) -> Value {
    let hdr = header.unwrap_or(&obs.header);
    let sample = obs.sample_ticks();
    let ticks = &obs.epoch_ticks;
    let idx = select_epochs(ticks, start, end, interval, sample);
    let mut sel: Vec<i64> = idx.iter().map(|&i| ticks[i]).collect();
    sel.sort_unstable();
    let decimating = interval.is_some_and(|i| i != 0 && i > sample);
    let eff = if decimating {
        interval.unwrap_or(sample)
    } else {
        sample
    };
    let (first, last) = (sel.first().copied(), sel.last().copied());
    let duration = match (first, last) {
        (Some(a), Some(b)) => b - a,
        _ => 0,
    };
    let expected = if sel.is_empty() {
        0
    } else {
        duration.div_euclid(eff) + 1
    };
    let whole = start.is_none() && end.is_none() && !decimating;

    let mut systems: BTreeMap<u8, SysSum> = BTreeMap::new();
    for (sat, st) in &obs.sats {
        let (mut n_ep, n_slip) = if whole {
            (st.ticks.len(), st.slips.len())
        } else {
            let (lo, hi) = (first.unwrap_or(0), last.unwrap_or(-1));
            let cnt = |v: &[i64]| {
                v.partition_point(|&t| t <= hi)
                    .saturating_sub(v.partition_point(|&t| t < lo))
            };
            (cnt(&st.ticks), cnt(&st.slips))
        };
        if !whole && decimating && n_ep > 0 {
            let iv = interval.unwrap_or(sample) as f64;
            n_ep = (geomaticist_gnss_core::text::round_half_even(n_ep as f64 * sample as f64 / iv)
                as usize)
                .max(1);
        }
        if n_ep == 0 {
            continue;
        }
        let d = systems.entry(sat.sys.letter() as u8).or_default();
        d.sats.push(sat.to_string());
        d.epochs += n_ep;
        d.slips += n_slip;
        for (c, n) in &st.counts {
            *d.obs.entry(c.to_string()).or_default() += n;
        }
        for (c, n) in &st.ssi_n {
            let a = d.ssi.entry(c.to_string()).or_default();
            a.0 += st.ssi_sum.get(c).copied().unwrap_or(0);
            a.1 += n;
        }
        for (c, n) in &st.snr_n {
            let a = d.snr.entry(c.to_string()).or_default();
            a.0 += st.snr_sum.get(c).copied().unwrap_or(0.0);
            a.1 += n;
        }
    }
    let order = SYSTEM_ORDER.as_bytes();
    let mut keys: Vec<u8> = order
        .iter()
        .copied()
        .filter(|s| systems.contains_key(s))
        .collect();
    keys.extend(systems.keys().filter(|s| !order.contains(s)));
    let mut sys_json = Map::new();
    let (mut n_sats, mut n_slips) = (0, 0);
    for s in keys {
        let d = &systems[&s];
        let avg1 = |m: Vec<(String, f64, usize)>| -> Map<String, Value> {
            m.into_iter()
                .filter(|x| x.2 > 0)
                .map(|(c, a, n)| (c, json!(round_to(a / n as f64, 1))))
                .collect()
        };
        n_sats += d.sats.len();
        n_slips += d.slips;
        sys_json.insert(
            char::from(s).to_string(),
            json!({
                "name": system_name(s), "satellites": d.sats, "epochs": d.epochs,
                "slips": d.slips, "obs_counts": d.obs,
                "ssi": avg1(d.ssi.iter().map(|(c, v)| (c.clone(), v.0 as f64, v.1)).collect()),
                "snr": avg1(d.snr.iter().map(|(c, v)| (c.clone(), v.0, v.1)).collect()),
                "count": d.sats.len(),
            }),
        );
    }
    let nsat: Vec<i64> = idx.iter().map(|&i| obs.epoch_nsat[i]).collect();
    let llh = hdr
        .approx_xyz()
        .filter(|x| x.iter().any(|v| *v != 0.0))
        .map(|x| xyz_to_llh(x[0], x[1], x[2]));
    let tf = |t: Option<i64>| t.map_or(Value::Null, |t| json!(fmt(t, 7)));
    let mut out = json!({
        "file": obs.source_path.display().to_string(),
        "rinex_version": hdr.version_str,
        "selection": {"start": tf(first), "end": tf(last), "whole_file": whole},
        "first_epoch": tf(first),
        "last_epoch": tf(last),
        "duration_s": duration as f64 / TICKS as f64,
        "duration": fmt_duration(duration),
        "interval_header_s": hdr.interval(),
        "interval_actual_s": sample as f64 / TICKS as f64,
        "interval_output_s": eff as f64 / TICKS as f64,
        "epochs_in_file": ticks.len(),
        "epochs": sel.len(),
        "epochs_expected": expected,
        "completeness_pct": if expected > 0 {
            round_to(100.0 * sel.len() as f64 / expected as f64, 1)
        } else {
            0.0
        },
        "gaps": gaps(&sel, eff, 100).iter().map(|(a, b)| json!({
            "from": fmt(*a, 1), "to": fmt(*b, 1), "seconds": (b - a) as f64 / TICKS as f64,
        })).collect::<Vec<_>>(),
        "satellites": n_sats,
        "mean_sats_per_epoch": if nsat.is_empty() {
            json!(0)
        } else {
            json!(round_to(nsat.iter().sum::<i64>() as f64 / nsat.len() as f64, 1))
        },
        "min_sats_per_epoch": nsat.iter().copied().min().unwrap_or(0),
        "systems": sys_json,
        "cycle_slips": n_slips,
        "events": obs.events.len(),
        "out_of_order_epochs": obs.out_of_order,
        "duplicate_epochs": obs.duplicates,
        "header": header_summary(hdr),
        "approx_llh": llh.map_or(Value::Null, |v| json!(v)),
    });
    let checks = opus_checks(obs, &out, first, last, eff, hdr);
    out["opus"] = Value::Array(checks);
    if !whole {
        out["included_epochs"] = if sel.len() <= 5000 {
            json!(sel.iter().map(|&t| fmt(t, 7)).collect::<Vec<_>>())
        } else {
            json!(format!("{} epochs (list omitted above 5000)", sel.len()))
        };
    }
    out
}

/// NGS OPUS's published requirements checked against a selection.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn opus_checks(
    obs: &RinexObs,
    s: &Value,
    first: Option<i64>,
    last: Option<i64>,
    interval: i64,
    hdr: &RinexHeader,
) -> Vec<Value> {
    let mut out = Vec::new();
    let mut add = |level: &str, text: String| out.push(json!({"level": level, "text": text}));
    let dur_h = s["duration_s"].as_f64().unwrap_or(0.0) / 3600.0;
    let dur = s["duration"].as_str().unwrap_or("");
    if dur_h < 0.25 {
        add(
            "fail",
            format!("Session is {dur} - OPUS needs at least 15 minutes."),
        );
    } else if dur_h < 2.0 {
        add(
            "ok",
            format!("{dur}: OPUS-RS (rapid static, 15 min to 2 h)."),
        );
    } else if dur_h <= 48.0 {
        add("ok", format!("{dur}: OPUS-S (static, 2 to 48 h)."));
    } else {
        add(
            "fail",
            format!("Session is {dur} - OPUS accepts at most 48 hours."),
        );
    }
    let secs = interval as f64 / TICKS as f64;
    if OPUS_RATES.iter().any(|r| (secs - r).abs() < 1e-6) {
        add(
            "ok",
            format!(
                "Interval {} s is one OPUS accepts (it decimates everything to 30 s).",
                format_g(secs)
            ),
        );
    } else {
        add(
            "warn",
            format!(
                "Interval {} s - OPUS accepts 1, 2, 3, 5, 10, 15 or 30 s.",
                format_g(secs)
            ),
        );
    }
    if let (Some(f), Some(l)) = (first, last) {
        let (a, b) = (from_ticks(f), from_ticks(l));
        let da = ordinal(a.year, a.month, a.day).unwrap_or(0);
        let db = ordinal(b.year, b.month, b.day).unwrap_or(0);
        let crossings = db - da;
        if crossings > 1 {
            add(
                "fail",
                format!("The session crosses UTC midnight {crossings} times - OPUS allows once."),
            );
        }
    }
    let g: Vec<String> = obs
        .obs_types_seen
        .get(&b'G')
        .filter(|v| !v.is_empty())
        .or_else(|| obs.header.obs_types.get(&b'G'))
        .map_or_default(|v| v.iter().map(ToString::to_string).collect());
    let has = |c: &str| g.iter().any(|x| x == c);
    let v2 = obs.is_v2();
    let has_l1 = g.iter().any(|c| c.starts_with("L1"));
    let has_l2 = g.iter().any(|c| c.starts_with("L2"));
    if g.is_empty() {
        add("fail", "No GPS observations - OPUS uses GPS only.".into());
    } else if !(has_l1 && has_l2) {
        add(
            "fail",
            "GPS needs both L1 and L2 carrier phase for OPUS.".into(),
        );
    } else {
        add("ok", "Dual-frequency GPS carrier phase (L1 + L2).".into());
        if dur_h < 2.0 {
            let p2 = if v2 {
                has("P2")
            } else {
                ["C2W", "C2P", "C2Y"].iter().any(|c| has(c))
            };
            let p1 = if v2 {
                ["P1", "C1"].iter().any(|c| has(c))
            } else {
                ["C1C", "C1W", "C1P"].iter().any(|c| has(c))
            };
            if !(p1 && p2) {
                add(
                    "fail",
                    "Under 2 hours OPUS needs P2 and P1 or C1 (RINEX 3: C2W and C1C/C1W).".into(),
                );
            }
        }
    }
    let ant = hdr.ant_type();
    if ant.is_empty() || ["UNKNOWN", "NONE"].contains(&ant.to_uppercase().as_str()) {
        add(
            "warn",
            "No antenna type in the header - set the NGS antenna name (Header page).".into(),
        );
    } else {
        let t = format!("Antenna {ant} {}", hdr.ant_radome());
        add(
            "ok",
            format!("{} - check it matches the NGS name.", t.trim_end()),
        );
    }
    match hdr.delta_hen() {
        Some(h) if h[0] != 0.0 => add("ok", format!("Antenna height (to ARP) {:.4} m.", h[0])),
        _ => add(
            "warn",
            "ANTENNA: DELTA H is 0 - you'll need to enter the ARP height in OPUS.".into(),
        ),
    }
    if hdr.approx_xyz().is_none_or(|x| x.iter().all(|v| *v == 0.0)) {
        add("warn", "No approximate position in the header.".into());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn select_and_decimate() {
        let t: Vec<i64> = (0..10).map(|i| i * 300_000_000 + 10_000).collect(); // 30 s, 1 ms late
        assert_eq!(
            select_epochs(&t, None, None, None, 300_000_000),
            (0..10).collect::<Vec<_>>()
        );
        assert_eq!(
            select_epochs(&t, None, None, Some(600_000_000), 300_000_000),
            vec![0, 2, 4, 6, 8]
        );
        assert_eq!(
            select_epochs(&t, Some(t[3]), Some(t[5]), None, 300_000_000),
            vec![3, 4, 5]
        );
        let mut u = t.clone();
        u.swap(1, 2);
        assert_eq!(
            select_epochs(&u, None, None, Some(600_000_000), 300_000_000),
            vec![0, 1, 4, 6, 8]
        );
        assert!(check_interval(Some(450_000_000), 300_000_000).is_err());
        assert!(check_interval(Some(600_000_000), 300_000_000).is_ok());
        assert_eq!(effective_interval(Some(300_000_000), 300_000_000), None);
        assert_eq!(gaps(&[0, 30, 60, 150], 30, 100), vec![(60, 150)]);
    }
}
