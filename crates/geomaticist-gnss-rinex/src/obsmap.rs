//! Which observation types go into a written file, and where each comes from in the source.
//!
//! Within RINEX 3 and 4 the 3-character codes are the same, except that RINEX 3.02 wrote BeiDou B1I as band 1 (3.03 and later: band 2). RINEX 2
//! has 2-character codes in one list for all systems; the tables follow the RINEX 3 specification's
//! own 2 <-> 3 correspondence.

use std::collections::{BTreeMap, BTreeSet};

use geomaticist_gnss_core::system::SYSTEM_ORDER;
use geomaticist_gnss_core::{GnssSystem, ObsCode, Sat};

/// Versions the writer writes.
pub const VERSIONS: [&str; 8] = [
    "2.11", "3.02", "3.03", "3.04", "3.05", "4.00", "4.01", "4.02",
];

/// Observation kinds and their names.
pub const KINDS: [(char, &str); 4] = [
    ('C', "Code (pseudorange)"),
    ('L', "Carrier phase"),
    ('D', "Doppler"),
    ('S', "Signal strength (SNR)"),
];

/// `float(version)` of a version string, 0 when unreadable.
#[must_use]
pub fn version_number(version: &str) -> f64 {
    geomaticist_gnss_core::text::py_float(version.as_bytes()).unwrap_or(0.0)
}

/// Systems a version can hold.
#[must_use]
pub fn systems_allowed(version: &str) -> &'static str {
    let v = version_number(version);
    if v < 3.0 {
        "GRES"
    } else if v < 3.03 {
        "GRECJS"
    } else {
        SYSTEM_ORDER
    }
}

/// System name for notes (`SYSTEM_NAMES.get(s, s)`).
#[must_use]
pub fn system_name(s: u8) -> String {
    GnssSystem::from_letter(char::from(s))
        .map_or_else(|| char::from(s).to_string(), |g| g.name().to_owned())
}

type V2Table = &'static [(&'static str, &'static [&'static str])];

/// RINEX 2 code -> RINEX 3 codes to take it from, best first, per system
/// (the first matching candidate wins).
#[must_use]
pub fn v2_from_v3(sys: u8) -> V2Table {
    match sys {
        b'G' => &[
            ("C1", &["C1C"]),
            ("P1", &["C1W", "C1P", "C1Y"]),
            ("C2", &["C2L", "C2X", "C2S", "C2C"]),
            ("P2", &["C2W", "C2P", "C2Y", "C2D"]),
            ("C5", &["C5X", "C5Q", "C5I"]),
            ("L1", &["L1C", "L1W", "L1P", "L1Y"]),
            (
                "L2",
                &["L2W", "L2P", "L2Y", "L2D", "L2X", "L2L", "L2S", "L2C"],
            ),
            ("L5", &["L5X", "L5Q", "L5I"]),
            ("D1", &["D1C", "D1W", "D1P"]),
            ("D2", &["D2W", "D2P", "D2X", "D2L", "D2S", "D2C"]),
            ("D5", &["D5X", "D5Q", "D5I"]),
            ("S1", &["S1C", "S1W", "S1P"]),
            ("S2", &["S2W", "S2P", "S2X", "S2L", "S2S", "S2C"]),
            ("S5", &["S5X", "S5Q", "S5I"]),
        ],
        b'R' => &[
            ("C1", &["C1C"]),
            ("P1", &["C1P"]),
            ("C2", &["C2C"]),
            ("P2", &["C2P"]),
            ("L1", &["L1C", "L1P"]),
            ("L2", &["L2P", "L2C"]),
            ("D1", &["D1C", "D1P"]),
            ("D2", &["D2P", "D2C"]),
            ("S1", &["S1C", "S1P"]),
            ("S2", &["S2P", "S2C"]),
        ],
        b'E' => &[
            ("C1", &["C1X", "C1C", "C1B"]),
            ("C5", &["C5X", "C5Q", "C5I"]),
            ("C7", &["C7X", "C7Q", "C7I"]),
            ("C8", &["C8X", "C8Q", "C8I"]),
            ("C6", &["C6X", "C6C", "C6B"]),
            ("L1", &["L1X", "L1C", "L1B"]),
            ("L5", &["L5X", "L5Q", "L5I"]),
            ("L7", &["L7X", "L7Q", "L7I"]),
            ("L8", &["L8X", "L8Q", "L8I"]),
            ("L6", &["L6X", "L6C", "L6B"]),
            ("D1", &["D1X", "D1C", "D1B"]),
            ("D5", &["D5X", "D5Q", "D5I"]),
            ("D7", &["D7X", "D7Q", "D7I"]),
            ("D8", &["D8X", "D8Q", "D8I"]),
            ("D6", &["D6X", "D6C", "D6B"]),
            ("S1", &["S1X", "S1C", "S1B"]),
            ("S5", &["S5X", "S5Q", "S5I"]),
            ("S7", &["S7X", "S7Q", "S7I"]),
            ("S8", &["S8X", "S8Q", "S8I"]),
            ("S6", &["S6X", "S6C", "S6B"]),
        ],
        b'S' => &[
            ("C1", &["C1C"]),
            ("C5", &["C5I", "C5Q", "C5X"]),
            ("L1", &["L1C"]),
            ("L5", &["L5I", "L5Q", "L5X"]),
            ("D1", &["D1C"]),
            ("D5", &["D5I", "D5Q", "D5X"]),
            ("S1", &["S1C"]),
            ("S5", &["S5I", "S5Q", "S5X"]),
        ],
        _ => &[],
    }
}

/// Order of RINEX 2 codes in a written `# / TYPES OF OBSERV`.
pub const V2_ORDER: [&str; 26] = [
    "L1", "L2", "L5", "L7", "L8", "L6", "C1", "P1", "C2", "P2", "C5", "C7", "C8", "C6", "D1", "D2",
    "D5", "D7", "D8", "D6", "S1", "S2", "S5", "S7", "S8", "S6",
];

fn oc(s: &str) -> ObsCode {
    ObsCode::new(s.as_bytes()).unwrap_or_default()
}

fn code3(kind: u8, band: u8, attr: u8) -> Option<ObsCode> {
    ObsCode::new(&[kind, band, attr])
}

/// The RINEX 3 code for a RINEX 2 code.
#[must_use]
pub fn v3_from_v2(sys: u8, code: ObsCode, v2_types: &[ObsCode]) -> Option<ObsCode> {
    let b = code.as_bytes();
    let (kind, band) = (*b.first()?, *b.get(1)?);
    let has_p2 = v2_types.iter().any(|c| c.as_bytes() == b"P2");
    if kind == b'P' {
        return match (sys, band) {
            (b'G', b'1') => Some(oc("C1W")),
            (b'G', b'2') => Some(oc("C2W")),
            (b'R', b'1') => Some(oc("C1P")),
            (b'R', b'2') => Some(oc("C2P")),
            _ => None,
        };
    }
    match sys {
        b'G' => match band {
            b'1' => code3(kind, b'1', b'C'),
            b'2' if kind == b'C' => Some(oc("C2X")),
            b'2' => code3(kind, b'2', if has_p2 { b'W' } else { b'X' }),
            b'5' => code3(kind, b'5', b'X'),
            _ => None,
        },
        b'R' => match band {
            b'1' => code3(kind, b'1', b'C'),
            b'2' if kind == b'C' => Some(oc("C2C")),
            b'2' => code3(kind, b'2', if has_p2 { b'P' } else { b'C' }),
            _ => None,
        },
        b'E' if b"15786".contains(&band) => code3(kind, band, b'X'),
        b'S' => match band {
            b'1' => code3(kind, b'1', b'C'),
            b'5' => code3(kind, b'5', b'I'),
            _ => None,
        },
        _ => None,
    }
}

/// BeiDou B1I: band 1 in RINEX 3.02, band 2 from 3.03.
#[must_use]
pub fn bds_band(code: ObsCode, to_302: bool) -> Option<ObsCode> {
    let b = code.as_bytes();
    let band = b.get(1).copied();
    let attr = b.get(2).copied();
    if to_302 {
        return match band {
            Some(b'2') => code3(b[0], b'1', attr.unwrap_or(b' ')),
            Some(b'1') => None, // B1C (BDS-3) doesn't exist in 3.02
            _ => Some(code),
        };
    }
    if band == Some(b'1') && attr.is_some_and(|a| b"IQX".contains(&a)) {
        return code3(b[0], b'2', attr.unwrap_or(b' '));
    }
    Some(code)
}

/// What to write and from where.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Plan {
    /// Target version (`"3.04"`).
    pub target: String,
    /// Output codes per system (RINEX 2: the shared list, per system).
    pub out_types: BTreeMap<u8, Vec<ObsCode>>,
    /// RINEX 2: the one list.
    pub v2_types: Vec<ObsCode>,
    /// Per system: source field index for each output code.
    pub field_map: BTreeMap<u8, Vec<Option<usize>>>,
    /// RINEX 2 output from RINEX 3: per satellite, the source field for each output code.
    pub sat_map: BTreeMap<Sat, Vec<Option<usize>>>,
    /// Per system: output indexes of the phase codes.
    pub phase_out: BTreeMap<u8, Vec<usize>>,
    /// Systems written, in order.
    pub systems: String,
    /// Remarks for the person (systems or codes left out).
    pub notes: Vec<String>,
}

impl Plan {
    /// `int(float(target))`.
    #[must_use]
    pub fn major(&self) -> i64 {
        version_number(&self.target) as i64
    }
}

fn kind_of(c: ObsCode) -> u8 {
    c.as_bytes().first().copied().unwrap_or(b' ')
}

fn join_codes(cs: &[ObsCode]) -> String {
    cs.iter()
        .map(ObsCode::as_str)
        .collect::<Vec<_>>()
        .join(", ")
}

/// Builds the plan. `src_types`: per system, the source's own codes (2 characters for RINEX 2,
/// 3 for 3/4). `sat_codes`: per satellite, the codes it has data for; going to RINEX 2 the source
/// code is then chosen per satellite.
#[must_use]
pub fn make_plan(
    src_types: &BTreeMap<u8, Vec<ObsCode>>,
    src_version: f64,
    target: &str,
    kinds: &str,
    systems: Option<&str>,
    sat_codes: Option<&BTreeMap<Sat, BTreeSet<ObsCode>>>,
) -> Plan {
    let mut plan = Plan {
        target: target.to_owned(),
        ..Plan::default()
    };
    let tgt = version_number(target);
    let allowed = systems_allowed(target).as_bytes();
    let want: &[u8] = match systems {
        Some(s) if !s.is_empty() => s.as_bytes(),
        _ => SYSTEM_ORDER.as_bytes(),
    };
    let order = SYSTEM_ORDER.as_bytes();
    let mut present: Vec<u8> = order
        .iter()
        .copied()
        .filter(|s| src_types.contains_key(s))
        .collect();
    present.extend(src_types.keys().filter(|s| !order.contains(s))); // sorted (BTreeMap)
    let dropped: Vec<u8> = present
        .iter()
        .copied()
        .filter(|s| want.contains(s) && !allowed.contains(s))
        .collect();
    if !dropped.is_empty() {
        let names: Vec<String> = dropped.iter().map(|&s| system_name(s)).collect();
        plan.notes.push(format!(
            "RINEX {target} can't hold {} - left out.",
            names.join(", ")
        ));
    }
    let use_: Vec<u8> = present
        .iter()
        .copied()
        .filter(|s| want.contains(s) && allowed.contains(s))
        .collect();
    plan.systems = use_.iter().map(|&s| char::from(s)).collect();
    let src_v2 = src_version < 3.0;
    let kinds = kinds.as_bytes();

    if tgt >= 3.0 {
        for &s in &use_ {
            let codes = &src_types[&s];
            let mut out: Vec<ObsCode> = Vec::new();
            let mut fmap: Vec<Option<usize>> = Vec::new();
            for (i, &c) in codes.iter().enumerate() {
                let mut c3 = if src_v2 {
                    if let Some(x) = v3_from_v2(s, c, codes) {
                        x
                    } else {
                        plan.notes.push(format!(
                            "{} {c}: no RINEX 3 equivalent - left out.",
                            char::from(s)
                        ));
                        continue;
                    }
                } else {
                    c
                };
                if s == b'C' {
                    let conv = if tgt < 3.03 && (src_v2 || src_version >= 3.03) {
                        bds_band(c3, true)
                    } else if tgt >= 3.03 && !src_v2 && src_version < 3.03 {
                        bds_band(c3, false)
                    } else {
                        Some(c3)
                    };
                    let Some(x) = conv else {
                        plan.notes
                            .push(format!("BeiDou {c}: not in RINEX {target} - left out."));
                        continue;
                    };
                    c3 = x;
                }
                if !kinds.contains(&kind_of(c3)) || out.contains(&c3) {
                    continue;
                }
                out.push(c3);
                fmap.push(Some(i));
            }
            if !out.is_empty() {
                plan.phase_out.insert(
                    s,
                    out.iter()
                        .enumerate()
                        .filter(|(_, c)| kind_of(**c) == b'L')
                        .map(|(j, _)| j)
                        .collect(),
                );
                plan.out_types.insert(s, out);
                plan.field_map.insert(s, fmap);
            }
        }
        return plan;
    }

    // ---- target RINEX 2
    let mut chosen: Vec<(u8, Vec<(ObsCode, usize)>)> = Vec::new(); // insertion-ordered dicts
    for &s in &use_ {
        let codes = &src_types[&s];
        let mut pick: Vec<(ObsCode, usize)> = Vec::new();
        if src_v2 {
            for (i, &c) in codes.iter().enumerate() {
                if !pick.iter().any(|(k, _)| *k == c) {
                    pick.push((c, i));
                }
            }
        } else {
            let table = v2_from_v3(s);
            for (v2c, cands) in table {
                for cand in *cands {
                    if let Some(i) = codes.iter().position(|c| c.as_str() == *cand) {
                        pick.push((oc(v2c), i));
                        break;
                    }
                }
            }
            let unmapped: Vec<ObsCode> = codes
                .iter()
                .copied()
                .filter(|c| !table.iter().any(|(_, cands)| cands.contains(&c.as_str())))
                .collect();
            if !unmapped.is_empty() {
                plan.notes.push(format!(
                    "{}: {} have no RINEX 2 code - left out.",
                    char::from(s),
                    join_codes(&unmapped)
                ));
            }
        }
        pick.retain(|(c, _)| {
            let k = kind_of(*c);
            kinds.contains(&if k == b'P' { b'C' } else { k })
        });
        if !pick.is_empty() {
            chosen.push((s, pick));
        }
    }
    let in_chosen = |c: &str| {
        chosen
            .iter()
            .any(|(_, p)| p.iter().any(|(k, _)| k.as_str() == c))
    };
    let mut v2: Vec<ObsCode> = V2_ORDER
        .iter()
        .filter(|c| in_chosen(c))
        .map(|c| oc(c))
        .collect();
    let extra: BTreeSet<ObsCode> = chosen
        .iter()
        .flat_map(|(_, p)| p.iter().map(|(k, _)| *k))
        .filter(|k| !v2.contains(k))
        .collect();
    v2.extend(extra);
    plan.v2_types.clone_from(&v2);
    for (s, pick) in &chosen {
        let get = |c: &ObsCode| pick.iter().find(|(k, _)| k == c).map(|(_, i)| *i);
        plan.field_map.insert(*s, v2.iter().map(get).collect());
        plan.out_types.insert(*s, v2.clone());
        plan.phase_out.insert(
            *s,
            v2.iter()
                .enumerate()
                .filter(|(_, c)| kind_of(**c) == b'L')
                .map(|(j, _)| j)
                .collect(),
        );
    }
    if let Some(sat_codes) = sat_codes.filter(|_| !src_v2) {
        for (sat, have) in sat_codes {
            let s = sat.sys.letter() as u8;
            let Some((_, pick)) = chosen.iter().find(|(x, _)| *x == s) else {
                continue;
            };
            let codes = &src_types[&s];
            let table = v2_from_v3(s);
            let m = v2
                .iter()
                .map(|c| {
                    let cands = table
                        .iter()
                        .find(|(k, _)| *k == c.as_str())
                        .map_or(&[][..], |(_, v)| *v);
                    let idx = cands.iter().find_map(|cand| {
                        let cc = oc(cand);
                        if have.contains(&cc) {
                            codes.iter().position(|x| *x == cc)
                        } else {
                            None
                        }
                    });
                    idx.filter(|_| pick.iter().any(|(k, _)| k == c))
                })
                .collect();
            plan.sat_map.insert(*sat, m);
        }
    }
    plan.systems = use_
        .iter()
        .filter(|s| chosen.iter().any(|(x, _)| x == *s))
        .map(|&s| char::from(s))
        .collect();
    plan
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codes(s: &str) -> Vec<ObsCode> {
        s.split_whitespace().map(oc).collect()
    }

    #[test]
    fn v2_to_v3_and_back() {
        let mut src = BTreeMap::new();
        src.insert(b'G', codes("L1 C1 L2 P2"));
        let p = make_plan(&src, 2.10, "3.04", "CL", None, None);
        assert_eq!(p.out_types[&b'G'], codes("L1C C1C L2W C2W"));
        assert_eq!(p.phase_out[&b'G'], vec![0, 2]);
        let mut src3 = BTreeMap::new();
        src3.insert(b'G', codes("C1C L1C C2W L2W C2L L2L S1C"));
        src3.insert(b'C', codes("C2I L2I"));
        let p2 = make_plan(&src3, 3.04, "2.11", "CL", None, None);
        assert_eq!(p2.v2_types, codes("L1 L2 C1 C2 P2"));
        assert_eq!(
            p2.field_map[&b'G'],
            vec![Some(1), Some(3), Some(0), Some(4), Some(2)]
        );
        assert_eq!(p2.systems, "G");
        assert!(p2.notes[0].starts_with("RINEX 2.11 can't hold BeiDou"));
        let p302 = make_plan(&src3, 3.04, "3.02", "CL", Some("C"), None);
        assert_eq!(p302.out_types[&b'C'], codes("C1I L1I"));
    }
}
