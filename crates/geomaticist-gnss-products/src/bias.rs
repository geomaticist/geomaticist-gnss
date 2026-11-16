//! Bias-SINEX observable-specific biases (OSB) and IGS Earth rotation parameters (ERP), and
//! Onsala ocean-loading coefficients (BLQ).

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::text::{latin1, py_splitlines, py_strip, sl};

use crate::util::{ProductError, np_interp, np_mean, pf, read_text, split_ws};

/// Satellite OSB from Bias-SINEX files: (satellite, observable) -> metres, the mean over the files'
/// records; `ns` × c, `m` as is, other units (`cyc`) left out.
pub fn read_osb<P: AsRef<Path>>(
    paths: &[P],
) -> Result<BTreeMap<(String, String), f64>, ProductError> {
    let mut acc: Vec<((String, String), Vec<f64>)> = Vec::new();
    for p in paths {
        let txt = read_text(p.as_ref())?;
        let mut on = false;
        for ln in py_splitlines(&txt) {
            if ln.starts_with(b"+BIAS/SOLUTION") {
                on = true;
                continue;
            }
            if ln.starts_with(b"-BIAS/SOLUTION") {
                on = false;
                continue;
            }
            if !on || !ln.starts_with(b" OSB ") {
                continue;
            }
            let prn = py_strip(sl(ln, 11, 14));
            let sta = py_strip(sl(ln, 15, 24));
            if !sta.is_empty() || prn.len() != 3 {
                continue;
            }
            let obs = latin1(py_strip(sl(ln, 25, 29)));
            let unit = py_strip(sl(ln, 65, 69));
            let Some(mut v) = pf(sl(ln, 70, 91)) else {
                continue;
            };
            if unit == b"ns" {
                v *= 1e-9 * CLIGHT;
            } else if unit != b"m" {
                continue;
            }
            let key = (latin1(prn), obs);
            match acc.iter_mut().find(|(k, _)| *k == key) {
                Some((_, vs)) => vs.push(v),
                None => acc.push((key, vec![v])),
            }
        }
    }
    Ok(acc.into_iter().map(|(k, v)| (k, np_mean(&v))).collect())
}

/// IGS ERP (version 2): MJD, x pole and y pole (arcsec), sorted and de-duplicated.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Erp {
    /// Modified Julian dates.
    pub mjd: Vec<f64>,
    /// x pole (arcsec).
    pub xp: Vec<f64>,
    /// y pole (arcsec).
    pub yp: Vec<f64>,
}

/// `\d{5}\.\d+` (ASCII digits).
fn is_mjd(s: &[u8]) -> bool {
    s.len() > 6
        && s[..5].iter().all(u8::is_ascii_digit)
        && s[5] == b'.'
        && s[6..].iter().all(u8::is_ascii_digit)
}

/// Reads ERP files.
pub fn read_erp<P: AsRef<Path>>(paths: &[P]) -> Result<Erp, ProductError> {
    let mut rows: Vec<(f64, f64, f64)> = Vec::new();
    for p in paths {
        let txt = read_text(p.as_ref())?;
        for ln in py_splitlines(&txt) {
            let f = split_ws(ln);
            if f.len() >= 3
                && is_mjd(f[0])
                && let (Some(a), Some(b), Some(c)) = (pf(f[0]), pf(f[1]), pf(f[2]))
            {
                rows.push((a, b * 1e-6, c * 1e-6));
            }
        }
    }
    // sorted(set(rows))
    rows.sort_by(|a, b| {
        a.0.total_cmp(&b.0)
            .then(a.1.total_cmp(&b.1))
            .then(a.2.total_cmp(&b.2))
    });
    rows.dedup();
    Ok(Erp {
        mjd: rows.iter().map(|r| r.0).collect(),
        xp: rows.iter().map(|r| r.1).collect(),
        yp: rows.iter().map(|r| r.2).collect(),
    })
}

/// Pole coordinates (arcsec) at a GPS time, linear in MJD; `None` without data.
#[must_use]
pub fn pole_at(erp: &Erp, t_gps: f64) -> Option<(f64, f64)> {
    if erp.mjd.is_empty() {
        return None;
    }
    let m = 44244.0 + t_gps / 86400.0;
    Some((
        np_interp(m, &erp.mjd, &erp.xp),
        np_interp(m, &erp.mjd, &erp.yp),
    ))
}

/// BLQ text -> station name (upper case) -> 6 x 11 rows: amplitudes (m) radial, west, south, then
/// phases (degrees, Greenwich lag).
#[must_use]
pub fn read_blq(text: &[u8]) -> BTreeMap<String, [[f64; 11]; 6]> {
    let lines: Vec<&[u8]> = py_splitlines(text)
        .into_iter()
        .filter(|l| !l.starts_with(b"$$") && !py_strip(l).is_empty())
        .collect();
    let mut out = BTreeMap::new();
    let mut i = 0;
    while i < lines.len() {
        let name = split_ws(lines[i])
            .first()
            .map(|n| latin1(n).to_uppercase())
            .unwrap_or_default();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        let mut j = i + 1;
        while j < lines.len() && rows.len() < 6 {
            let vals: Option<Vec<f64>> = split_ws(lines[j]).iter().map(|v| pf(v)).collect();
            let Some(mut vals) = vals else { break };
            vals.truncate(11);
            rows.push(vals);
            j += 1;
        }
        if rows.len() == 6 && rows.iter().all(|r| r.len() == 11) {
            let mut a = [[0.0; 11]; 6];
            for (k, r) in rows.iter().enumerate() {
                a[k].copy_from_slice(r);
            }
            out.insert(name, a);
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blq_and_erp_rules() {
        let mut t = String::from("$$ header\n  fair\n");
        for k in 0..6 {
            t.push_str(&vec![format!("{:.5}", 0.001 * f64::from(k)); 11].join(" "));
            t.push('\n');
        }
        t.push_str("  BAD\n 1 2 3\n");
        let b = read_blq(t.as_bytes());
        assert_eq!(b.keys().collect::<Vec<_>>(), vec!["FAIR"]);
        assert_eq!(b["FAIR"][5][10], 0.005);
        assert!(is_mjd(b"60576.50") && !is_mjd(b"6057.50") && !is_mjd(b"60576."));
        let e = Erp {
            mjd: vec![60576.5, 60577.5],
            xp: vec![0.1, 0.3],
            yp: vec![0.2, 0.2],
        };
        let t_mid = (60577.0 - 44244.0) * 86400.0;
        let (x, y) = pole_at(&e, t_mid).unwrap();
        assert!((x - 0.2).abs() < 1e-12 && (y - 0.2).abs() < 1e-15);
    }
}
