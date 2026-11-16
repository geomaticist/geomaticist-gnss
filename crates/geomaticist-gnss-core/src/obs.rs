//! Observations as the engine consumes them: one station's epochs, per satellite and code a value
//! series (NaN where missing) and the LLI of phases; signal choice per satellite and the
//! ionosphere-free / geometry-free / Melbourne-Wübbena combinations.
//!
//! The RINEX reading (`read_obs`) is in `geomaticist-gnss-rinex`.
//!
//! Layout: epoch times are one shared `Vec<f64>` of GPS seconds;
//! each satellite holds its signals in a `BTreeMap` by code (lookups by best code never depend on
//! map order except at an exact tie, which is broken by the order the code was first seen in the
//! file).

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::consts::{CLIGHT, pairs};
use crate::signal::{
    FreqError, ObsCode, SignalKey, attr_priority, attr_score, class_priority, freq, signal_class,
};
use crate::system::Sat;

/// Per band, the (phase code, code code) to use; either may be `None` (then the best at the
/// station). The shape of [`common_signals`]'s result per station.
pub type BandSignals = BTreeMap<u8, (Option<ObsCode>, Option<ObsCode>)>;

/// One code of one satellite at one station on the station's epoch grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Signal {
    /// What the code measures ([`crate::signal::code_key`]).
    pub key: SignalKey,
    /// Values per epoch: phase in cycles, code in metres, SNR in dB-Hz; NaN where missing.
    pub values: Vec<f64>,
    /// Loss-of-lock indicator per epoch (phases only; 0 where missing or blank).
    pub lli: Option<Vec<i8>>,
    /// Order in which the reader first saw this (satellite, code),
    /// used only to break exact ties in [`StationObs::best_code`].
    pub first_seen: u32,
}

impl Signal {
    /// Number of finite values.
    #[must_use]
    pub fn count(&self) -> usize {
        self.values.iter().filter(|v| v.is_finite()).count()
    }
}

/// The signals of one satellite, by code.
pub type SatSignals = BTreeMap<ObsCode, Signal>;

/// One station's observations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StationObs {
    /// Station name (four characters by default).
    pub name: String,
    /// The file it was read from.
    pub path: PathBuf,
    /// IGS antenna name + radome, 16 + 4 characters (`norm_antenna`).
    pub antenna: String,
    /// Receiver type.
    pub receiver: String,
    /// APPROX POSITION XYZ (m).
    pub approx: [f64; 3],
    /// ANTENNA: DELTA H/E/N (m).
    pub delta_hen: [f64; 3],
    /// Epochs, GPS seconds, sorted and unique.
    pub t: Vec<f64>,
    /// Epoch flag 1 (power failure) per epoch.
    pub power_fail: Vec<bool>,
    /// Satellites with at least one value, each with its signals.
    pub sats: BTreeMap<Sat, SatSignals>,
}

impl StationObs {
    /// An empty station.
    #[must_use]
    pub fn new(name: impl Into<String>, path: impl Into<PathBuf>) -> Self {
        Self {
            name: name.into(),
            path: path.into(),
            antenna: String::new(),
            receiver: String::new(),
            approx: [0.0; 3],
            delta_hen: [0.0; 3],
            t: Vec::new(),
            power_fail: Vec::new(),
            sats: BTreeMap::new(),
        }
    }

    /// Satellites, sorted.
    #[must_use]
    pub fn sat_list(&self) -> Vec<Sat> {
        self.sats.keys().copied().collect()
    }

    /// A copy holding epochs `t0 <= t <= t1` only (1e-6 s slack); satellites with no finite value
    /// left are dropped.
    #[must_use]
    pub fn window(&self, t0: f64, t1: f64) -> StationObs {
        let m: Vec<bool> = self
            .t
            .iter()
            .map(|&t| t >= t0 - 1e-6 && t <= t1 + 1e-6)
            .collect();
        let pick_f = |v: &[f64]| -> Vec<f64> {
            v.iter()
                .zip(&m)
                .filter(|(_, k)| **k)
                .map(|(x, _)| *x)
                .collect()
        };
        let pick_i = |v: &[i8]| -> Vec<i8> {
            v.iter()
                .zip(&m)
                .filter(|(_, k)| **k)
                .map(|(x, _)| *x)
                .collect()
        };
        let mut out = StationObs {
            name: self.name.clone(),
            path: self.path.clone(),
            antenna: self.antenna.clone(),
            receiver: self.receiver.clone(),
            approx: self.approx,
            delta_hen: self.delta_hen,
            t: pick_f(&self.t),
            power_fail: self
                .power_fail
                .iter()
                .zip(&m)
                .filter(|(_, k)| **k)
                .map(|(x, _)| *x)
                .collect(),
            sats: BTreeMap::new(),
        };
        for (sat, sigs) in &self.sats {
            let dd: SatSignals = sigs
                .iter()
                .map(|(c, s)| {
                    (
                        *c,
                        Signal {
                            key: s.key,
                            values: pick_f(&s.values),
                            lli: s.lli.as_deref().map(pick_i),
                            first_seen: s.first_seen,
                        },
                    )
                })
                .collect();
            if dd.values().any(|s| s.values.iter().any(|v| v.is_finite())) {
                out.sats.insert(*sat, dd);
            }
        }
        out
    }

    /// Joins the observations of one receiver session that was cut into several files (a drone that
    /// starts a new file per flight while its receiver keeps tracking): epochs are merged in time order,
    /// and an epoch present in more than one part keeps the first part's values. The header facts (name,
    /// antenna, approximate position, antenna height) are those of the first part.
    ///
    /// Returns `None` for an empty list.
    #[must_use]
    pub fn concat(parts: &[StationObs]) -> Option<StationObs> {
        let first = parts.first()?;
        if parts.len() == 1 {
            return Some(first.clone());
        }
        // (time, part, index) in time order; ties keep the earlier part
        let mut order: Vec<(f64, usize, usize)> = Vec::new();
        for (p, st) in parts.iter().enumerate() {
            order.extend(st.t.iter().enumerate().map(|(i, &t)| (t, p, i)));
        }
        order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        order.dedup_by(|b, a| (a.0 - b.0).abs() < 1e-7);
        let n = order.len();
        let mut out = StationObs {
            name: first.name.clone(),
            path: first.path.clone(),
            antenna: first.antenna.clone(),
            receiver: first.receiver.clone(),
            approx: first.approx,
            delta_hen: first.delta_hen,
            t: order.iter().map(|o| o.0).collect(),
            power_fail: order
                .iter()
                .map(|&(_, p, i)| parts[p].power_fail.get(i).copied().unwrap_or(false))
                .collect(),
            sats: BTreeMap::new(),
        };
        let mut seen: u32 = 0;
        for st in parts {
            for (sat, sigs) in &st.sats {
                for (code, sig) in sigs {
                    let slot = out.sats.entry(*sat).or_default();
                    slot.entry(*code).or_insert_with(|| {
                        seen += 1;
                        Signal {
                            key: sig.key,
                            values: vec![f64::NAN; n],
                            lli: sig.lli.as_ref().map(|_| vec![0; n]),
                            first_seen: seen - 1,
                        }
                    });
                }
            }
        }
        for (k, &(_, p, i)) in order.iter().enumerate() {
            for (sat, sigs) in &parts[p].sats {
                let Some(slot) = out.sats.get_mut(sat) else {
                    continue;
                };
                for (code, sig) in sigs {
                    let Some(dst) = slot.get_mut(code) else {
                        continue;
                    };
                    dst.values[k] = sig.values.get(i).copied().unwrap_or(f64::NAN);
                    if let (Some(d), Some(s)) = (dst.lli.as_mut(), sig.lli.as_ref()) {
                        d[k] = s.get(i).copied().unwrap_or(0);
                    }
                }
            }
        }
        Some(out)
    }

    /// The code of `kind` on `band` this satellite carries most, ties by attribute priority (then by
    /// first appearance).
    #[must_use]
    pub fn best_code(&self, sat: Sat, kind: char, band: u8) -> Option<ObsCode> {
        let sigs = self.sats.get(&sat)?;
        let pri = attr_priority(sat.sys.letter(), band);
        sigs.iter()
            .filter(|(_, s)| s.key.kind == kind && s.key.band == band)
            .min_by_key(|(_, s)| {
                (
                    std::cmp::Reverse(s.count()),
                    attr_score(pri, s.key.attr),
                    s.first_seen,
                )
            })
            .map(|(c, _)| *c)
    }

    /// Bands with both phase and code on this satellite.
    #[must_use]
    pub fn bands(&self, sat: Sat) -> BTreeSet<u8> {
        let Some(sigs) = self.sats.get(&sat) else {
            return BTreeSet::new();
        };
        sigs.values()
            .map(|s| s.key.band)
            .collect::<BTreeSet<u8>>()
            .into_iter()
            .filter(|&b| {
                self.best_code(sat, 'L', b).is_some() && self.best_code(sat, 'C', b).is_some()
            })
            .collect()
    }

    fn signal(&self, sat: Sat, code: ObsCode) -> Option<&Signal> {
        self.sats.get(&sat)?.get(&code)
    }

    /// Phase (cycles) and LLI of a band (best code, or `code`). NaN / 0 when absent.
    #[must_use]
    pub fn phase(&self, sat: Sat, band: u8, code: Option<ObsCode>) -> (Vec<f64>, Vec<i8>) {
        let n = self.t.len();
        let c = code.or_else(|| self.best_code(sat, 'L', band));
        match c.and_then(|c| self.signal(sat, c)) {
            Some(s) => (
                s.values.clone(),
                s.lli.clone().unwrap_or_else(|| vec![0; n]),
            ),
            None => (vec![f64::NAN; n], vec![0; n]),
        }
    }

    /// Code (m) of a band (best code, or `code`); NaN when absent.
    #[must_use]
    pub fn code(&self, sat: Sat, band: u8, code: Option<ObsCode>) -> Vec<f64> {
        let c = code.or_else(|| self.best_code(sat, 'C', band));
        c.and_then(|c| self.signal(sat, c))
            .map_or_else(|| vec![f64::NAN; self.t.len()], |s| s.values.clone())
    }

    /// Signal strength (dB-Hz) of a band; NaN when absent.
    #[must_use]
    pub fn snr(&self, sat: Sat, band: u8) -> Vec<f64> {
        self.best_code(sat, 'S', band)
            .and_then(|c| self.signal(sat, c))
            .map_or_else(|| vec![f64::NAN; self.t.len()], |s| s.values.clone())
    }
}

/// `'TRM59800.99     SCIT'` (as RINEX writes it, 20 characters) -> IGS 16 + 4 (`NONE` when the
/// radome is blank).
#[must_use]
pub fn norm_antenna(ant: &str) -> String {
    let chars: Vec<char> = ant.chars().collect();
    let words: Vec<&str> = ant.split_whitespace().collect();
    let a: String = if chars.len() >= 16 {
        chars[..16].iter().collect::<String>().trim().to_owned()
    } else {
        words.first().map_or(String::new(), |w| (*w).to_owned())
    };
    let rad: String = if chars.len() > 16 {
        chars[16..chars.len().min(20)]
            .iter()
            .collect::<String>()
            .trim()
            .to_owned()
    } else {
        words.get(1).map_or(String::new(), |w| (*w).to_owned())
    };
    if a.is_empty() {
        return String::new();
    }
    let rad = if rad.is_empty() {
        "NONE".to_owned()
    } else {
        rad
    };
    format!("{a:<16}{rad:>4}")
}

/// The dual-frequency pair for `sat`, common to all `stations`: the first in `PAIRS` order that has
/// at least 80 % of the best pair's epochs with both phases and the first code at every station.
#[must_use]
pub fn pick_pair(sat: Sat, stations: &[&StationObs]) -> Option<(u8, u8)> {
    // Share of the best pair's epochs a pair needs.
    const PAIR_SHARE_MIN: f64 = 0.8;
    let cands = pairs(sat.sys);
    let mut best = 0usize;
    let mut counts = Vec::with_capacity(cands.len());
    for &(b1, b2) in cands {
        let mut n: Option<usize> = None;
        for st in stations {
            if !st.sats.contains_key(&sat) {
                n = Some(0);
                break;
            }
            let (l1, _) = st.phase(sat, b1, None);
            let (l2, _) = st.phase(sat, b2, None);
            let c1 = st.code(sat, b1, None);
            let k = l1
                .iter()
                .zip(&l2)
                .zip(&c1)
                .filter(|((a, b), c)| a.is_finite() && b.is_finite() && c.is_finite())
                .count();
            n = Some(n.map_or(k, |n| n.min(k)));
        }
        let n = n.unwrap_or(0);
        counts.push(n);
        best = best.max(n);
    }
    if best == 0 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    cands
        .iter()
        .zip(counts)
        .find(|(_, n)| *n as f64 >= PAIR_SHARE_MIN * best as f64)
        .map(|(p, _)| *p)
}

/// Per station, `{band: (phase code, code code)}` of the same signal class at every station (either
/// may be `None`: the caller then takes that station's best), with
/// `frac` (default 0.5).
#[must_use]
pub fn common_signals(
    sat: Sat,
    pair: (u8, u8),
    stations: &[&StationObs],
    frac: f64,
) -> Vec<BandSignals> {
    let mut out: Vec<BandSignals> = vec![BTreeMap::new(); stations.len()];
    let sysc = sat.sys.letter();
    for band in [pair.0, pair.1] {
        for kind in ['L', 'C'] {
            let Some(pri) = class_priority(sysc, kind, band) else {
                continue;
            };
            // per station: code -> (class, count)
            let per: Vec<Vec<(ObsCode, String, usize)>> = stations
                .iter()
                .map(|st| {
                    st.sats.get(&sat).map_or_else(Vec::new, |sigs| {
                        sigs.iter()
                            .filter(|(_, s)| s.key.kind == kind && s.key.band == band)
                            .map(|(c, s)| (*c, signal_class(sysc, band, s.key.attr), s.count()))
                            .collect()
                    })
                })
                .collect();
            let mut pick: Option<Vec<ObsCode>> = None;
            for cl in pri.chars() {
                let mut names = Vec::new();
                let mut ok = true;
                for cnt in &per {
                    let best = cnt.iter().map(|x| x.2).max().unwrap_or(0);
                    #[allow(clippy::cast_precision_loss)]
                    let have = cnt
                        .iter()
                        .filter(|(_, k, n)| {
                            k.chars().eq(std::iter::once(cl))
                                && *n as f64 >= frac * best as f64
                                && *n > 0
                        })
                        .max_by(|a, b| (a.2, a.0).cmp(&(b.2, b.0)));
                    if let Some(h) = have {
                        names.push(h.0);
                    } else {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    pick = Some(names);
                    break;
                }
            }
            let Some(pick) = pick else { continue };
            for (i, c) in pick.into_iter().enumerate() {
                let e = out[i].entry(band).or_insert((None, None));
                if kind == 'L' {
                    e.0 = Some(c);
                } else {
                    e.1 = Some(c);
                }
            }
        }
    }
    out
}

/// Ionosphere-free combination `X_IF = a X1 - b X2`: (a, b) = (f1², f2²) / (f1² - f2²).
#[must_use]
pub fn if_coeffs(f1: f64, f2: f64) -> (f64, f64) {
    let d = f1 * f1 - f2 * f2;
    (f1 * f1 / d, f2 * f2 / d)
}

/// One satellite at one station on the station's epochs, for one signal pair.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Combo {
    /// Ionosphere-free phase (m).
    pub phase_if: Vec<f64>,
    /// Ionosphere-free code (m; NaN where a code is missing).
    pub code_if: Vec<f64>,
    /// L1 - L2 in metres (geometry-free).
    pub gf: Vec<f64>,
    /// Melbourne-Wübbena wide lane (cycles).
    pub mw: Vec<f64>,
    /// Loss of lock on either phase (LLI bit 0).
    pub lli: Vec<bool>,
    /// First-band code (m): transmit time.
    pub code1: Vec<f64>,
    /// First frequency (Hz).
    pub f1: f64,
    /// Second frequency (Hz).
    pub f2: f64,
    /// Half-cycle ambiguity on either phase (LLI bit 1).
    pub half: Vec<bool>,
}

/// The combinations of `pair` for one satellite. `signals` (`{band: (phase code, code code)}`, from
/// [`common_signals`]) picks the observations; otherwise each band's best at this station.
pub fn combine(
    st: &StationObs,
    sat: Sat,
    pair: (u8, u8),
    fcn: Option<i32>,
    signals: Option<&BandSignals>,
) -> Result<Combo, FreqError> {
    let (b1, b2) = pair;
    let sg = |b: u8| {
        signals
            .and_then(|s| s.get(&b))
            .copied()
            .unwrap_or((None, None))
    };
    let f1 = freq(sat, b1, fcn)?;
    let f2 = freq(sat, b2, fcn)?;
    let (l1, ll1) = st.phase(sat, b1, sg(b1).0);
    let (l2, ll2) = st.phase(sat, b2, sg(b2).0);
    let p1 = st.code(sat, b1, sg(b1).1);
    let p2 = st.code(sat, b2, sg(b2).1);
    let (lam1, lam2) = (CLIGHT / f1, CLIGHT / f2);
    let (a, b) = if_coeffs(f1, f2);
    let lw = CLIGHT / (f1 - f2);
    let n = st.t.len();
    let mut c = Combo {
        phase_if: Vec::with_capacity(n),
        code_if: Vec::with_capacity(n),
        gf: Vec::with_capacity(n),
        mw: Vec::with_capacity(n),
        lli: Vec::with_capacity(n),
        code1: p1.clone(),
        f1,
        f2,
        half: Vec::with_capacity(n),
    };
    for i in 0..n {
        let big_l1 = l1[i] * lam1;
        let big_l2 = l2[i] * lam2;
        c.phase_if.push(a * big_l1 - b * big_l2);
        c.code_if.push(a * p1[i] - b * p2[i]);
        c.gf.push(big_l1 - big_l2);
        c.mw.push(
            ((f1 * big_l1 - f2 * big_l2) / (f1 - f2) - (f1 * p1[i] + f2 * p2[i]) / (f1 + f2)) / lw,
        );
        c.lli.push((ll1[i] & 1) > 0 || (ll2[i] & 1) > 0);
        c.half.push((ll1[i] & 2) > 0 || (ll2[i] & 2) > 0);
    }
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signal::code_key;

    fn station() -> StationObs {
        let mut st = StationObs::new("TEST", "x.obs");
        st.t = vec![0.0, 30.0, 60.0, 90.0];
        st.power_fail = vec![false; 4];
        let g: Sat = "G05".parse().unwrap();
        let mut sigs = SatSignals::new();
        let nan = f64::NAN;
        let mut add = |code: &str, vals: Vec<f64>, order: u32| {
            let c: ObsCode = code.parse().unwrap();
            let key = code_key(code.as_bytes(), 'G', 3, 3.04);
            let lli = (key.kind == 'L').then(|| vec![0, 1, 0, 2]);
            sigs.insert(
                c,
                Signal {
                    key,
                    values: vals,
                    lli,
                    first_seen: order,
                },
            );
        };
        add("C1C", vec![2.0e7, 2.0e7 + 1.0, nan, 2.0e7 + 3.0], 0);
        add("C1W", vec![2.0e7, 2.0e7 + 1.0, 2.0e7 + 2.0, nan], 1);
        add(
            "L1C",
            vec![1.0e8, 1.0e8 + 5.0, 1.0e8 + 10.0, 1.0e8 + 15.0],
            2,
        );
        add(
            "C2W",
            vec![2.0e7 + 4.0, 2.0e7 + 5.0, 2.0e7 + 6.0, 2.0e7 + 7.0],
            3,
        );
        add(
            "L2W",
            vec![7.8e7, 7.8e7 + 4.0, 7.8e7 + 8.0, 7.8e7 + 12.0],
            4,
        );
        st.sats.insert(g, sigs);
        st
    }

    #[test]
    fn best_code_and_bands() {
        let st = station();
        let g: Sat = "G05".parse().unwrap();
        // C1C and C1W both have 3 values; W scores 1 in "PWYCMLSXN", C scores 3.
        assert_eq!(st.best_code(g, 'C', 1).unwrap().as_str(), "C1W");
        assert_eq!(st.bands(g).into_iter().collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(pick_pair(g, &[&st]), Some((1, 2)));
        let w = st.window(30.0, 60.0);
        assert_eq!(w.t, vec![30.0, 60.0]);
        assert_eq!(
            w.sats[&g][&"L1C".parse::<ObsCode>().unwrap()].lli,
            Some(vec![1, 0])
        );
    }

    #[test]
    fn combinations() {
        let st = station();
        let g: Sat = "G05".parse().unwrap();
        let c = combine(&st, g, (1, 2), None, None).unwrap();
        let (a, b) = if_coeffs(1575.42e6, 1227.60e6);
        let lam1 = CLIGHT / 1575.42e6;
        let lam2 = CLIGHT / 1227.60e6;
        assert_eq!(c.phase_if[0], a * (1.0e8 * lam1) - b * (7.8e7 * lam2));
        assert_eq!(c.lli, vec![false, true, false, false]);
        assert_eq!(c.half, vec![false, false, false, true]);
        assert!(c.code_if[3].is_nan()); // C1W missing there
        let cs = common_signals(g, (1, 2), &[&st, &st], 0.5);
        assert_eq!(cs[0][&1].0.unwrap().as_str(), "L1C");
        assert_eq!(cs[0][&1].1.unwrap().as_str(), "C1W");
        assert_eq!(cs[1][&2].0.unwrap().as_str(), "L2W");
    }

    #[test]
    fn antenna_names() {
        assert_eq!(norm_antenna("TRM59800.99     SCIT"), "TRM59800.99     SCIT");
        assert_eq!(norm_antenna("TRM59800.99         "), "TRM59800.99     NONE");
        assert_eq!(norm_antenna("ASH701945G_M"), "ASH701945G_M    NONE");
        assert_eq!(norm_antenna(""), "");
    }
}
