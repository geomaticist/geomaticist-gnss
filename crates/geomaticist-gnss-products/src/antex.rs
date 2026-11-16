//! ANTEX 1.4: receiver and satellite antenna calibrations (PCO, PCV) and the receiver range
//! correction.
//!
//! Only the entries asked for are parsed from the (large) file.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock};

use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::text::{latin1, py_rstrip, py_splitlines, py_strip, sl};
use geomaticist_gnss_core::time::{DateTime, gps_seconds};

use crate::util::{
    ProductError, is_alpha_latin1, np_arange, np_divmod, np_interp, pf, pi, split_ws,
};

/// Carrier frequency of `const.FREQ[sys][band]` (`None` where the table has no entry, which
/// includes every NavIC and SBAS band).
fn freq_table(sys: u8, band: u8) -> Option<f64> {
    let s = GnssSystem::from_letter(char::from(sys))?;
    if matches!(s, GnssSystem::Navic | GnssSystem::Sbas) {
        return None;
    }
    geomaticist_gnss_core::consts::base_freq(s, band)
}

/// One ANTEX entry.
#[derive(Debug, Clone, PartialEq)]
pub struct AntennaCal {
    /// Type + radome (receivers) or antenna type (satellites).
    pub name: String,
    /// Serial number field (satellites: the PRN / SVN columns).
    pub serial: String,
    /// Zenith (receivers) or nadir (satellites) grid, degrees.
    pub zen: Vec<f64>,
    /// Azimuth step (degrees), 0 = no azimuth rows.
    pub dazi: f64,
    /// Frequency key (`"G01"`) -> PCO (m): receivers N, E, U; satellites x, y, z.
    pub pco: BTreeMap<String, [f64; 3]>,
    /// Frequency key -> NOAZI PCV row (m).
    pub noazi: BTreeMap<String, Vec<f64>>,
    /// Frequency key -> PCV rows per azimuth (m).
    pub grid: BTreeMap<String, Vec<Vec<f64>>>,
    /// Valid from (GPS seconds).
    pub valid_from: f64,
    /// Valid until (GPS seconds).
    pub valid_until: f64,
}

impl AntennaCal {
    fn empty() -> Self {
        Self {
            name: String::new(),
            serial: String::new(),
            zen: Vec::new(),
            dazi: 0.0,
            pco: BTreeMap::new(),
            noazi: BTreeMap::new(),
            grid: BTreeMap::new(),
            valid_from: -1e30,
            valid_until: 1e30,
        }
    }

    /// The calibration used for `band` of `sys`: its own, else GPS's of the same band, else the
    /// nearest GPS one.
    #[must_use]
    pub fn freq_key(&self, sys: u8, band: u8) -> Option<String> {
        let k = format!("{}{band:02}", char::from(sys));
        if self.pco.contains_key(&k) {
            return Some(k);
        }
        let g = format!("G{band:02}");
        if self.pco.contains_key(&g) {
            return Some(g);
        }
        let mut alt = match band {
            5 | 7 | 8 | 6 | 3 => "G02",
            _ => "G01",
        };
        if sys == b'R' && band == 3 {
            alt = if self.pco.contains_key("R02") {
                "R02"
            } else {
                "G02"
            };
        }
        self.pco.contains_key(alt).then(|| alt.to_owned())
    }

    /// How [`Self::freq_key`] models a band: (`calibrated` | `same frequency` | `substituted` |
    /// `none`, key).
    #[must_use]
    pub fn freq_status(&self, sys: u8, band: u8) -> (&'static str, Option<String>) {
        let Some(key) = self.freq_key(sys, band) else {
            return ("none", None);
        };
        if key == format!("{}{band:02}", char::from(sys)) {
            return ("calibrated", Some(key));
        }
        let f = freq_table(sys, band);
        let kb = key.as_bytes();
        let fk = std::str::from_utf8(&kb[1..])
            .ok()
            .and_then(|s| s.parse::<u8>().ok())
            .and_then(|b| freq_table(kb[0], b));
        if let (Some(f), Some(fk)) = (f, fk)
            && f != 0.0
            && fk != 0.0
            && (f - fk).abs() < 1e6
        {
            return ("same frequency", Some(key));
        }
        ("substituted", Some(key))
    }

    /// PCV (m) at zenith / nadir angles (degrees), from the NOAZI row unless azimuths (degrees) are
    /// given and the entry has azimuth rows.
    #[must_use]
    pub fn pcv(&self, key: &str, zen_deg: &[f64], az_deg: Option<&[f64]>) -> Vec<f64> {
        let Some(noazi) = self.noazi.get(key) else {
            return vec![0.0; zen_deg.len()];
        };
        let (Some(&z0), Some(&z1)) = (self.zen.first(), self.zen.last()) else {
            return vec![0.0; zen_deg.len()]; // no ZEN1 / ZEN2 / DZEN (no correction)
        };
        let clip = |z: f64| z.max(z0).min(z1);
        let grid = self.grid.get(key);
        match (az_deg, grid) {
            (Some(az), Some(g)) if self.dazi > 0.0 => zen_deg
                .iter()
                .zip(az)
                .map(|(&z, &a)| {
                    let z = clip(z);
                    let (_, a) = np_divmod(a, 360.0);
                    let (q, _) = np_divmod(a, self.dazi);
                    let i = (q as i64).clamp(0, g.len() as i64 - 2) as usize;
                    let w = (a - i as f64 * self.dazi) / self.dazi;
                    let lo = np_interp(z, &self.zen, &g[i]);
                    let hi = np_interp(z, &self.zen, &g[i + 1]);
                    (1.0 - w) * lo + w * hi
                })
                .collect(),
            _ => zen_deg
                .iter()
                .map(|&z| np_interp(clip(z), &self.zen, noazi))
                .collect(),
        }
    }
}

/// `"YYYY MM DD HH MM SS.S"` -> GPS seconds.
fn antex_date(s: &[u8]) -> Option<f64> {
    let p = split_ws(s);
    if p.len() < 6 {
        return None;
    }
    let n: Vec<i64> = p[..5].iter().map(|v| pi(v)).collect::<Option<_>>()?;
    let dt = DateTime::new(n[0], n[1], n[2], n[3], n[4], 0)?;
    Some(gps_seconds(&dt) + pf(p[5])?)
}

fn find(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if from > hay.len() {
        return None;
    }
    hay[from..]
        .windows(needle.len())
        .position(|w| w == needle)
        .map(|i| i + from)
}

fn rfind(hay: &[u8], needle: &[u8], end: usize) -> Option<usize> {
    let end = end.min(hay.len());
    hay[..end].windows(needle.len()).rposition(|w| w == needle)
}

/// One entry's text -> calibration.
fn parse_block(txt: &[u8]) -> AntennaCal {
    let mut cal = AntennaCal::empty();
    let mut cur: Option<String> = None;
    let mut rows: Vec<Vec<f64>> = Vec::new();
    let f = |ln: &[u8], a: usize, b: usize| pf(sl(ln, a, b));
    for ln in py_splitlines(txt) {
        let lab_raw = sl(ln, 60, 80);
        let lab: &[u8] = if lab_raw.iter().any(|&c| is_alpha_latin1(c)) {
            py_strip(lab_raw)
        } else {
            b""
        };
        match lab {
            b"TYPE / SERIAL NO" => {
                cal.name = latin1(py_rstrip(sl(ln, 0, 20)));
                cal.serial = latin1(py_strip(sl(ln, 20, 40)));
            }
            b"DAZI" => {
                if let Some(v) = f(ln, 2, 8) {
                    cal.dazi = v;
                }
            }
            _ if lab.starts_with(b"ZEN1 / ZEN2 / DZEN") => {
                if let (Some(z1), Some(z2), Some(dz)) = (f(ln, 2, 8), f(ln, 8, 14), f(ln, 14, 20)) {
                    cal.zen = np_arange(z1, z2 + dz / 2.0, dz);
                }
            }
            b"VALID FROM" => {
                if let Some(v) = antex_date(sl(ln, 0, 43)) {
                    cal.valid_from = v;
                }
            }
            b"VALID UNTIL" => {
                if let Some(v) = antex_date(sl(ln, 0, 43)) {
                    cal.valid_until = v;
                }
            }
            b"START OF FREQUENCY" => {
                cur = Some(latin1(py_strip(sl(ln, 3, 6))));
                rows.clear();
            }
            b"NORTH / EAST / UP" if cur.as_ref().is_some_and(|c| !c.is_empty()) => {
                if let (Some(a), Some(b), Some(c)) = (f(ln, 0, 10), f(ln, 10, 20), f(ln, 20, 30)) {
                    let k = cur.clone().unwrap_or_default();
                    cal.pco.insert(k, [a * 1e-3, b * 1e-3, c * 1e-3]);
                }
            }
            b"END OF FREQUENCY" if cur.as_ref().is_some_and(|c| !c.is_empty()) => {
                if !rows.is_empty() {
                    let k = cur.clone().unwrap_or_default();
                    cal.grid.insert(k, std::mem::take(&mut rows));
                }
                cur = None;
            }
            b"" if cur.as_ref().is_some_and(|c| !c.is_empty()) && !py_strip(ln).is_empty() => {
                let nz = cal.zen.len();
                let k = cur.clone().unwrap_or_default();
                if sl(ln, 3, 8) == b"NOAZI" {
                    let vals: Option<Vec<f64>> = split_ws(sl(ln, 8, ln.len()))
                        .into_iter()
                        .take(nz)
                        .map(|v| pf(v).map(|x| x * 1e-3))
                        .collect();
                    if let Some(v) = vals {
                        cal.noazi.insert(k, v);
                    }
                } else {
                    let vals = split_ws(ln);
                    let row: Option<Vec<f64>> = vals
                        .iter()
                        .skip(1)
                        .take(nz)
                        .map(|v| pf(v).map(|x| x * 1e-3))
                        .collect();
                    if let Some(r) = row {
                        rows.push(r);
                    }
                }
            }
            _ => {}
        }
    }
    cal
}

/// An ANTEX file held in memory; entries are parsed on request.
///
/// Clones share the bytes and the parsed entries (the type list, the receiver entries asked for, the
/// satellite table): an application that keeps one `Antex` per file and hands each solve its own clone
/// (several at once) parses the 25 MB of `igs20.atx` once, not once per baseline.
#[derive(Debug, Clone)]
pub struct Antex {
    raw: Arc<Vec<u8>>,
    rec: Arc<Mutex<BTreeMap<String, Option<AntennaCal>>>>,
    sat: Arc<OnceLock<BTreeMap<String, Vec<AntennaCal>>>>,
    types: Arc<OnceLock<BTreeSet<String>>>,
}

impl Antex {
    /// Reads the file.
    pub fn open(path: &Path) -> Result<Self, ProductError> {
        let raw =
            std::fs::read(path).map_err(|e| ProductError(format!("{}: {e}", path.display())))?;
        Ok(Self::from_bytes(raw))
    }

    /// From the file's bytes.
    #[must_use]
    pub fn from_bytes(raw: Vec<u8>) -> Self {
        Self {
            raw: Arc::new(raw),
            rec: Arc::new(Mutex::new(BTreeMap::new())),
            sat: Arc::new(OnceLock::new()),
            types: Arc::new(OnceLock::new()),
        }
    }

    /// Every receiver antenna type + radome, right-stripped (computed once per file).
    #[must_use]
    pub fn types(&self) -> BTreeSet<String> {
        self.types.get_or_init(|| self.scan_types()).clone()
    }

    fn scan_types(&self) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        // re.findall(rb"^(.{60})TYPE / SERIAL NO", raw, re.M): 60 bytes without \n at a line start
        let raw: &[u8] = &self.raw;
        let mut start = 0;
        while start <= raw.len() {
            let end = find(raw, b"\n", start).unwrap_or(raw.len());
            let line = &raw[start..end];
            if line.len() >= 60 + 16
                && &line[60..76] == b"TYPE / SERIAL NO"
                && py_strip_bytes(&line[20..40]).is_empty()
            {
                out.insert(latin1(py_rstrip(&line[..20])));
            }
            start = end + 1;
        }
        out
    }

    /// The calibration of an IGS antenna name + radome; a missing radome falls back to NONE.
    pub fn receiver(&mut self, name: &str) -> Option<AntennaCal> {
        let chars: Vec<char> = name.chars().collect();
        let part = |a: usize, b: usize| -> String {
            chars[a.min(chars.len())..b.min(chars.len())]
                .iter()
                .collect::<String>()
                .trim()
                .to_owned()
        };
        let radome = part(16, 20);
        let name = format!(
            "{:<16}{:>4}",
            part(0, 16),
            if radome.is_empty() { "NONE" } else { &radome }
        );
        if let Some(c) = self.rec.lock().ok().and_then(|m| m.get(&name).cloned()) {
            return c;
        }
        let mut cal = None;
        let head: String = name.chars().take(16).collect();
        for cand in [name.clone(), format!("{head}NONE")] {
            let cb = cand.as_bytes();
            let raw: &[u8] = &self.raw;
            let mut i = 0;
            while let Some(j) = find(raw, cb, i) {
                let ls = rfind(raw, b"\n", j).map_or(0, |p| p + 1);
                let le = find(raw, b"\n", j);
                let line = &raw[ls..le.unwrap_or(raw.len().saturating_sub(1)).max(ls)];
                if sl(line, 60, line.len()).starts_with(b"TYPE / SERIAL NO")
                    && py_rstrip(sl(line, 0, 20)) == py_rstrip(cb)
                    && py_strip(sl(line, 20, 40)).is_empty()
                {
                    let s = rfind(raw, b"START OF ANTENNA", ls).unwrap_or(0);
                    let e = le
                        .and_then(|le| find(raw, b"END OF ANTENNA", le))
                        .unwrap_or(raw.len());
                    cal = Some(parse_block(&raw[s..e]));
                    break;
                }
                i = j + 1;
            }
            if let Some(c) = cal.as_mut() {
                c.name = cand;
                break;
            }
        }
        if let Ok(mut m) = self.rec.lock() {
            m.insert(name, cal.clone());
        }
        cal
    }

    /// Satellite (`"G05"`) -> calibrations by validity period.
    pub fn satellites(&mut self) -> &BTreeMap<String, Vec<AntennaCal>> {
        self.sat.get_or_init(|| {
            let mut out: BTreeMap<String, Vec<AntennaCal>> = BTreeMap::new();
            let raw: &[u8] = &self.raw;
            let mut from = 0;
            while let Some(m) = find(raw, b"START OF ANTENNA", from) {
                let body = m + b"START OF ANTENNA".len();
                let Some(end) = find(raw, b"END OF ANTENNA", body) else {
                    break;
                };
                let txt = &raw[body..end];
                from = body;
                let mut parts = txt.splitn(3, |&b| b == b'\n');
                let _ = parts.next();
                let first = if txt.contains(&b'\n') {
                    parts.next().unwrap_or(b"")
                } else {
                    b""
                };
                if !sl(first, 60, first.len()).starts_with(b"TYPE / SERIAL NO") {
                    continue;
                }
                let prn = py_strip(sl(first, 20, 23));
                if prn.len() == 3
                    && b"GRECJ".contains(&prn[0])
                    && prn[1..].iter().all(u8::is_ascii_digit)
                {
                    out.entry(latin1(prn)).or_default().push(parse_block(txt));
                }
            }
            out
        })
    }

    /// The calibration valid at `t_gps`.
    pub fn satellite(&mut self, sat: &str, t_gps: f64) -> Option<AntennaCal> {
        self.satellites()
            .get(sat)?
            .iter()
            .find(|c| c.valid_from <= t_gps && t_gps <= c.valid_until)
            .cloned()
    }
}

/// `bytes.strip()` (ASCII whitespace only).
fn py_strip_bytes(s: &[u8]) -> &[u8] {
    let ws = |c: &u8| matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c);
    let a = s.iter().position(|c| !ws(c)).unwrap_or(s.len());
    let b = s.iter().rposition(|c| !ws(c)).map_or(a, |i| i + 1);
    &s[a..b]
}

/// Range correction (m) of the antenna phase centre against its ARP: `-PCO·e + PCV`, at
/// azimuths / elevations (radians).
#[must_use]
pub fn receiver_correction(
    cal: Option<&AntennaCal>,
    sys: u8,
    band: u8,
    az: &[f64],
    el: &[f64],
    use_azimuth: bool,
) -> Vec<f64> {
    let Some(cal) = cal else {
        return vec![0.0; el.len()];
    };
    let Some(key) = cal.freq_key(sys, band) else {
        return vec![0.0; el.len()];
    };
    let [n, e, u] = cal.pco[&key];
    let zen: Vec<f64> = el.iter().map(|x| 90.0 - x.to_degrees()).collect();
    let azd: Vec<f64> = az.iter().map(|x| x.to_degrees()).collect();
    let pcv = cal.pcv(&key, &zen, use_azimuth.then_some(azd.as_slice()));
    az.iter()
        .zip(el)
        .zip(pcv)
        .map(|((&a, &l), p)| {
            let ce = l.cos();
            let los = e * a.sin() * ce + n * a.cos() * ce + u * l.sin();
            -los + p
        })
        .collect()
}

/// Satellite PCV of the ionosphere-free combination (m), at nadir angles (degrees).
#[must_use]
pub fn sat_pcv_if(
    cal: Option<&AntennaCal>,
    sys: u8,
    b1: u8,
    b2: u8,
    a: f64,
    b: f64,
    nadir_deg: &[f64],
) -> Vec<f64> {
    let keys = cal.and_then(|c| Some((c, c.freq_key(sys, b1)?, c.freq_key(sys, b2)?)));
    let Some((cal, k1, k2)) = keys else {
        return vec![0.0; nadir_deg.len()];
    };
    let p1 = cal.pcv(&k1, nadir_deg, None);
    let p2 = cal.pcv(&k2, nadir_deg, None);
    p1.iter().zip(&p2).map(|(x, y)| a * x - b * y).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const ATX: &str =
        "     1.4            M                                       ANTEX VERSION / SYST
                                                            END OF HEADER
                                                            START OF ANTENNA
TEST1           NONE                                        TYPE / SERIAL NO
   180.0                                                    DAZI
     0.0  10.0   5.0                                        ZEN1 / ZEN2 / DZEN
   G01                                                      START OF FREQUENCY
      1.00      2.00     60.00                              NORTH / EAST / UP
   NOAZI    0.00   -1.00   -2.00
     0.0    0.00   -1.00   -2.00
   180.0    0.00   -3.00   -4.00
   360.0    0.00   -1.00   -2.00
   G01                                                      END OF FREQUENCY
                                                            END OF ANTENNA
";

    #[test]
    fn parses_and_corrects() {
        let mut a = Antex::from_bytes(ATX.as_bytes().to_vec());
        assert_eq!(
            a.types().into_iter().collect::<Vec<_>>(),
            vec!["TEST1           NONE"]
        );
        let c = a.receiver("TEST1").unwrap();
        assert_eq!(c.name, "TEST1           NONE");
        assert_eq!(c.zen, vec![0.0, 5.0, 10.0]);
        assert_eq!(c.pco["G01"], [1e-3, 2e-3, 60e-3]);
        assert_eq!(c.grid["G01"].len(), 3);
        assert_eq!(
            c.freq_status(b'E', 1),
            ("same frequency", Some("G01".into()))
        );
        assert_eq!(c.freq_status(b'E', 5), ("none", None));
        assert_eq!(c.pcv("G01", &[7.5], None), vec![-1.5e-3]);
        assert_eq!(c.pcv("G01", &[5.0], Some(&[90.0])), vec![-2e-3]);
        let r = receiver_correction(
            Some(&c),
            b'G',
            1,
            &[0.0],
            &[std::f64::consts::FRAC_PI_2],
            false,
        );
        assert!((r[0] + 60e-3).abs() < 1e-15);
        assert!(a.receiver("NOPE").is_none());
    }
}
