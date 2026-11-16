//! Which signals of each satellite the two receivers share, and how they group for differencing.
//!
//! A double difference only cancels the receivers' hardware delays when both satellites are tracked the
//! same way at each receiver. Signals are therefore grouped by system, frequency band and the pair of
//! tracking attributes (rover / base), and differences are formed inside a group only. Two receivers that
//! track a band differently (`L1C` at one, `L1X` at the other) still form a consistent group, because the
//! offset between the two tracking modes is the same for every satellite in it.

use std::collections::BTreeMap;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::signal::freq;
use geomaticist_gnss_core::{GnssSystem, ObsCode, Sat, StationObs};

use super::options::PpkOptions;
use crate::orbits::Orbits;

/// One frequency of one satellite, as both receivers observe it.
#[derive(Debug, Clone, PartialEq)]
pub struct BandPlan {
    /// RINEX band number.
    pub band: u8,
    /// Rover phase code.
    pub rover_l: ObsCode,
    /// Rover pseudorange code.
    pub rover_c: ObsCode,
    /// Rover Doppler code, when the file has one.
    pub rover_d: Option<ObsCode>,
    /// Rover signal-strength code, when the file has one.
    pub rover_s: Option<ObsCode>,
    /// Base phase code.
    pub base_l: ObsCode,
    /// Base pseudorange code.
    pub base_c: ObsCode,
    /// Carrier frequency, Hz.
    pub freq: f64,
    /// Wavelength, m.
    pub lam: f64,
    /// Index into [`Plan::groups`].
    pub group: usize,
    /// Index of this ambiguity among all of the session.
    pub key: usize,
}

/// The frequencies of one satellite.
#[derive(Debug, Clone, PartialEq)]
pub struct SatPlan {
    /// The satellite.
    pub sat: Sat,
    /// Its shared frequencies, lowest band number first.
    pub bands: Vec<BandPlan>,
}

/// Signals that may be differenced with each other.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    /// System.
    pub sys: GnssSystem,
    /// RINEX band number.
    pub band: u8,
    /// Rover phase code.
    pub rover_l: ObsCode,
    /// Base phase code.
    pub base_l: ObsCode,
    /// Whether the ambiguities of the group are integers that may be fixed.
    pub fixable: bool,
    /// Satellites in the group.
    pub n_sats: usize,
}

/// The signal plan of a session.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Plan {
    /// The satellites both receivers observed, with their shared frequencies.
    pub sats: Vec<SatPlan>,
    /// The groups.
    pub groups: Vec<Group>,
    /// Number of (satellite, frequency) ambiguities.
    pub n_keys: usize,
    /// Satellites both receivers observed but that were left out, with the reason.
    pub left_out: Vec<(Sat, String)>,
}

fn count(st: &StationObs, sat: Sat, code: ObsCode) -> usize {
    st.sats
        .get(&sat)
        .and_then(|s| s.get(&code))
        .map_or(0, geomaticist_gnss_core::Signal::count)
}

/// The codes of `kind` on `band` with data, most data first.
fn codes(st: &StationObs, sat: Sat, kind: char, band: u8) -> Vec<(ObsCode, Option<char>, usize)> {
    let Some(sigs) = st.sats.get(&sat) else {
        return Vec::new();
    };
    let mut v: Vec<(ObsCode, Option<char>, usize, u32)> = sigs
        .iter()
        .filter(|(_, s)| s.key.kind == kind && s.key.band == band)
        .map(|(c, s)| (*c, s.key.attr, s.count(), s.first_seen))
        .filter(|x| x.2 > 0)
        .collect();
    v.sort_by(|a, b| b.2.cmp(&a.2).then(a.3.cmp(&b.3)));
    v.into_iter().map(|x| (x.0, x.1, x.2)).collect()
}

/// The code of `kind` on `band` with attribute `attr`, else the receiver's best.
fn code_like(
    st: &StationObs,
    sat: Sat,
    kind: char,
    band: u8,
    attr: Option<char>,
) -> Option<ObsCode> {
    let all = codes(st, sat, kind, band);
    all.iter()
        .find(|c| c.1 == attr)
        .or_else(|| all.first())
        .map(|c| c.0)
}

impl Plan {
    /// Builds the plan of a rover / base pair.
    #[must_use]
    pub fn build(
        rover: &StationObs,
        base: &StationObs,
        orb: &Orbits,
        fcn: &BTreeMap<Sat, i32>,
        opts: &PpkOptions,
    ) -> Plan {
        let mut plan = Plan::default();
        let mut group_of: BTreeMap<(char, u8, ObsCode, ObsCode), usize> = BTreeMap::new();
        for sat in rover.sats.keys().filter(|s| base.sats.contains_key(*s)) {
            let sat = *sat;
            if !opts.systems.contains(sat.sys.letter()) || opts.exclude.contains(&sat) {
                continue;
            }
            if !orb.has(sat) {
                plan.left_out.push((sat, "no orbit".into()));
                continue;
            }
            let chan = fcn.get(&sat).copied();
            let mut bands: Vec<u8> = rover
                .bands(sat)
                .intersection(&base.bands(sat))
                .copied()
                .collect();
            bands.sort_unstable();
            let mut sp = SatPlan {
                sat,
                bands: Vec::new(),
            };
            for band in bands {
                if sp.bands.len() >= opts.max_bands {
                    break;
                }
                let Ok(f) = freq(sat, band, chan) else {
                    continue;
                };
                // the pair of phase codes: the same attribute at both ends when there is one, then most data
                let rl = codes(rover, sat, 'L', band);
                let bl = codes(base, sat, 'L', band);
                let mut best: Option<(bool, usize, ObsCode, ObsCode)> = None;
                for r in &rl {
                    for b in &bl {
                        let cand = (r.1 == b.1, r.2.min(b.2), r.0, b.0);
                        let better = best.is_none_or(|x| (cand.0, cand.1) > (x.0, x.1));
                        if better {
                            best = Some(cand);
                        }
                    }
                }
                let Some((_, _, rover_l, base_l)) = best else {
                    continue;
                };
                let attr_of = |st: &StationObs, c: ObsCode| {
                    st.sats
                        .get(&sat)
                        .and_then(|s| s.get(&c))
                        .and_then(|s| s.key.attr)
                };
                let (ra, ba) = (attr_of(rover, rover_l), attr_of(base, base_l));
                let (Some(rover_c), Some(base_c)) = (
                    code_like(rover, sat, 'C', band, ra),
                    code_like(base, sat, 'C', band, ba),
                ) else {
                    continue;
                };
                let rover_d = code_like(rover, sat, 'D', band, ra);
                let rover_s = code_like(rover, sat, 'S', band, ra);
                let gkey = (sat.sys.letter(), band, rover_l, base_l);
                let ng = group_of.len();
                let group = *group_of.entry(gkey).or_insert(ng);
                if group == plan.groups.len() {
                    plan.groups.push(Group {
                        sys: sat.sys,
                        band,
                        rover_l,
                        base_l,
                        fixable: sat.sys != GnssSystem::Glonass || opts.ar_glonass,
                        n_sats: 0,
                    });
                }
                plan.groups[group].n_sats += 1;
                sp.bands.push(BandPlan {
                    band,
                    rover_l,
                    rover_c,
                    rover_d,
                    rover_s,
                    base_l,
                    base_c,
                    freq: f,
                    lam: CLIGHT / f,
                    group,
                    key: plan.n_keys,
                });
                plan.n_keys += 1;
            }
            if sp.bands.is_empty() {
                let why = if sat.sys == GnssSystem::Glonass && chan.is_none() {
                    "GLONASS frequency channel unknown"
                } else {
                    "no frequency with code and phase at both receivers"
                };
                plan.left_out.push((sat, why.into()));
            } else {
                plan.sats.push(sp);
            }
        }
        plan
    }

    /// One line per system: its frequencies as rover / base phase codes and the satellite count.
    #[must_use]
    pub fn describe(&self) -> Vec<String> {
        let mut by: BTreeMap<char, Vec<String>> = BTreeMap::new();
        for g in &self.groups {
            by.entry(g.sys.letter()).or_default().push(format!(
                "{}/{} ({} sats{})",
                g.rover_l.as_str(),
                g.base_l.as_str(),
                g.n_sats,
                if g.fixable { "" } else { ", float" }
            ));
        }
        by.into_iter()
            .map(|(s, v)| format!("{s} {}", v.join(" ")))
            .collect()
    }

    /// Number of satellites.
    #[must_use]
    pub fn n_sats(&self) -> usize {
        self.sats.len()
    }

    /// How many (satellite, frequency) observations `rover` and `base` could share at best.
    #[must_use]
    pub fn data_count(&self, rover: &StationObs) -> usize {
        self.sats
            .iter()
            .map(|s| {
                s.bands
                    .iter()
                    .map(|b| count(rover, s.sat, b.rover_l))
                    .sum::<usize>()
            })
            .sum()
    }
}
