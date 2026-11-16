//! The epoch table the streaming decoders (RTCM 3, u-blox, NovAtel, GREIS) fill: per epoch, per satellite
//! and RINEX code tail (`1C`, `2W` ...), one value set; then the observations as series and the RINEX file.

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::obs::Signal;
use geomaticist_gnss_core::signal::code_key;
use geomaticist_gnss_core::{ObsCode, Sat, StationObs};

use crate::{Converted, DecodeError};

/// One value set of a (satellite, signal) at an epoch; NaN where the message has no value.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Cell {
    pub code: f64,
    pub phase: f64,
    pub doppler: f64,
    pub snr: f64,
    pub lli: i8,
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            code: f64::NAN,
            phase: f64::NAN,
            doppler: f64::NAN,
            snr: f64::NAN,
            lli: 0,
        }
    }
}

/// Epochs keyed by GPS milliseconds.
#[derive(Debug, Default)]
pub(crate) struct Epochs {
    pub map: BTreeMap<i64, BTreeMap<(Sat, &'static str), Cell>>,
}

impl Epochs {
    /// The cell of a satellite's signal at GPS time `t` (s), made empty when new.
    pub fn cell(&mut self, t: f64, sat: Sat, tail: &'static str) -> &mut Cell {
        let key = (t * 1000.0).round() as i64;
        self.map
            .entry(key)
            .or_default()
            .entry((sat, tail))
            .or_default()
    }

    /// The observations as series.
    pub fn into_obs(self, name: &str) -> StationObs {
        let n = self.map.len();
        let mut st = StationObs::new(name, "");
        st.t = self.map.keys().map(|k| *k as f64 * 0.001).collect();
        st.power_fail = vec![false; n];
        let mut order = 0u32;
        for (i, cells) in self.map.values().enumerate() {
            for ((sat, tail), cell) in cells {
                let sysc = sat.sys.letter();
                for (kind, v) in [
                    ('C', cell.code),
                    ('L', cell.phase),
                    ('D', cell.doppler),
                    ('S', cell.snr),
                ] {
                    if !v.is_finite() {
                        continue;
                    }
                    let text = format!("{kind}{tail}");
                    let Ok(code) = text.parse::<ObsCode>() else {
                        continue;
                    };
                    let sig = st
                        .sats
                        .entry(*sat)
                        .or_default()
                        .entry(code)
                        .or_insert_with(|| {
                            order += 1;
                            Signal {
                                key: code_key(text.as_bytes(), sysc, 3, 3.04),
                                values: vec![f64::NAN; n],
                                lli: (kind == 'L').then(|| vec![0; n]),
                                first_seen: order - 1,
                            }
                        });
                    sig.values[i] = v;
                    if let Some(l) = sig.lli.as_mut() {
                        l[i] = cell.lli;
                    }
                }
            }
        }
        st
    }
}

/// An approximate GPS time of a log from its modification time (settles the week of time-of-week streams).
pub(crate) fn near_from_mtime(path: &Path) -> f64 {
    let unix = std::fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0.0, |d| d.as_secs_f64());
    // 1970-01-01 is 315 964 800 s before 1980-01-06; GPS runs ahead of UTC by the leap seconds
    unix - 315_964_800.0 + 18.0
}

/// The file stem and the 4-character marker name made from it.
pub(crate) fn stem_and_name(raw: &Path) -> (String, String) {
    let stem = raw
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let name: String = stem.chars().take(4).collect::<String>().to_uppercase();
    (stem, name)
}

/// Writes decoded observations as RINEX 3.04 `{stem}.obs` in `outdir`.
pub(crate) fn write_obs(
    raw: &Path,
    outdir: &Path,
    obs: &StationObs,
    glonass_fcn: BTreeMap<Sat, i32>,
    what: &str,
    warnings: Vec<String>,
    info: serde_json::Value,
) -> Result<Converted, DecodeError> {
    let (stem, _) = stem_and_name(raw);
    let header = geomaticist_gnss_rinex::station::StationHeader {
        program: "geomaticist-gnss".to_owned(),
        run_by: String::new(),
        date: crate::now_utc_text(),
        comments: vec![format!(
            "Decoded from the {what} log {}",
            raw.file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_default()
        )],
        glonass_fcn,
    };
    let text = geomaticist_gnss_rinex::station::station_to_rinex3(obs, &header);
    let path = outdir.join(format!("{stem}.obs"));
    crate::write_atomic(&path, text.as_bytes())?;
    Ok(Converted {
        obs: path,
        nav: Vec::new(),
        warnings,
        info,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use geomaticist_gnss_core::GnssSystem;

    #[test]
    fn cells_into_series() {
        let mut e = Epochs::default();
        let g1 = Sat::new(GnssSystem::Gps, 1);
        e.cell(100.0, g1, "1C").code = 2.0e7;
        e.cell(101.0, g1, "1C").phase = 1.0e8;
        e.cell(101.0, g1, "1C").lli = 1;
        let st = e.into_obs("TEST");
        assert_eq!(st.t, vec![100.0, 101.0]);
        let sigs = &st.sats[&g1];
        assert_eq!(sigs.len(), 2);
        let l = sigs.values().find(|s| s.lli.is_some()).unwrap();
        assert!(l.values[0].is_nan());
        assert_eq!(l.values[1], 1.0e8);
        assert_eq!(l.lli.as_ref().unwrap()[1], 1);
    }
}
