//! Precise orbits and clocks: SP3 (a, c, d) and RINEX clock files, and positions by 10th-order
//! Lagrange interpolation in an Earth-fixed frame rotated to the epoch.
//!
//! Reading rules: an epoch within 1 ms of the previous one is not added (later
//! records go to the slot NumPy's `searchsorted` gives), all-zero positions are left out, SP3 clocks
//! of 999999 or more are missing; positions km -> m, clocks µs -> s. Satellite ids with an unknown
//! system letter are left out.

use std::collections::BTreeMap;
use std::path::Path;

use geomaticist_gnss_core::Sat;
use geomaticist_gnss_core::consts::OMGE;
use geomaticist_gnss_core::time::{DateTime, gps_seconds};

use crate::util::{
    ProductError, lines, np_sum, pf, pi, read_text, searchsorted, sort_paths, split_ws,
};

/// Lagrange order (11 points; RTKLIB uses 10th order too).
pub const ORDER: usize = 10;

/// `ln[a:b].replace(" ", "0")` parsed as a satellite.
fn sat_id(s: &[u8]) -> Option<Sat> {
    let t: String = s
        .iter()
        .map(|&b| if b == b' ' { '0' } else { char::from(b) })
        .collect();
    t.parse().ok()
}

fn epoch_seconds(y: i64, mo: i64, d: i64, h: i64, mi: i64, sec: f64) -> Option<f64> {
    Some(gps_seconds(&DateTime::new(y, mo, d, h, mi, 0)?) + sec)
}

/// Insertion-ordered `dict[float, V]` (a repeated key keeps its place, the value is replaced).
#[derive(Debug, Default)]
struct KeyedRows<V> {
    rows: Vec<(f64, V)>,
}

impl<V> KeyedRows<V> {
    fn set(&mut self, k: f64, v: V) {
        if let Some(r) = self.rows.iter_mut().rev().find(|r| r.0 == k) {
            r.1 = v;
        } else {
            self.rows.push((k, v));
        }
    }
}

/// SP3 orbits (centre of mass), optionally with RINEX clock files.
#[derive(Debug, Clone, Default)]
pub struct Precise {
    /// SP3 epochs (GPS seconds), increasing.
    pub t: Vec<f64>,
    /// Positions (m) per satellite at the epochs, NaN where missing.
    pub pos: BTreeMap<Sat, Vec<[f64; 3]>>,
    /// SP3 clocks (s), NaN where missing.
    pub clk_sp3: BTreeMap<Sat, Vec<f64>>,
    /// RINEX clock records per satellite: (times, clocks), sorted.
    pub clk: ClockSeries,
}

impl Precise {
    /// Reads SP3 files (sorted by path) and clock files.
    pub fn open<P: AsRef<Path>, Q: AsRef<Path>>(
        sp3_paths: &[P],
        clk_paths: &[Q],
    ) -> Result<Self, ProductError> {
        let mut epochs: Vec<f64> = Vec::new();
        let mut pos: BTreeMap<Sat, KeyedRows<[f64; 3]>> = BTreeMap::new();
        let mut clk: BTreeMap<Sat, KeyedRows<f64>> = BTreeMap::new();
        for p in sort_paths(sp3_paths) {
            let text = read_text(&p)?;
            let bad = |what: &str| ProductError(format!("{}: {what}", p.display()));
            let mut cur: Option<f64> = None;
            for ln in lines(&text) {
                if ln.first() == Some(&b'*') {
                    let parts = split_ws(&ln[1..]);
                    if parts.len() < 6 {
                        return Err(bad("short epoch line"));
                    }
                    let n: Option<Vec<i64>> = parts[..5].iter().map(|v| pi(v)).collect();
                    let (Some(n), Some(sec)) = (n, pf(parts[5])) else {
                        return Err(bad("unreadable epoch line"));
                    };
                    let c = epoch_seconds(n[0], n[1], n[2], n[3], n[4], sec)
                        .ok_or_else(|| bad("invalid epoch date"))?;
                    cur = Some(c);
                    let last = epochs.last().copied();
                    if last.is_some_and(|l| (c - l).abs() < 1e-3) {
                        // the same epoch again
                    } else if last.is_none_or(|l| c > l) {
                        epochs.push(c);
                    }
                } else if ln.first() == Some(&b'P')
                    && let Some(c) = cur
                {
                    let f = |a: usize, b: usize| pf(geomaticist_gnss_core::text::sl(ln, a, b));
                    let (Some(x), Some(y), Some(z)) = (f(4, 18), f(18, 32), f(32, 46)) else {
                        continue;
                    };
                    let ck = geomaticist_gnss_core::text::sl(ln, 46, 60);
                    let cv = if ln.len() > 46 && !geomaticist_gnss_core::text::is_blank(ck) {
                        match pf(ck) {
                            Some(v) => v,
                            None => continue,
                        }
                    } else {
                        999_999.999_999
                    };
                    if x == 0.0 && y == 0.0 && z == 0.0 {
                        continue;
                    }
                    let Some(sat) = sat_id(geomaticist_gnss_core::text::sl(ln, 1, 4)) else {
                        continue;
                    };
                    pos.entry(sat)
                        .or_default()
                        .set(c, [x * 1e3, y * 1e3, z * 1e3]);
                    if cv.abs() < 999_999.0 {
                        clk.entry(sat).or_default().set(c, cv * 1e-6);
                    }
                }
            }
        }
        let n = epochs.len();
        let mut out = Self {
            t: epochs,
            ..Self::default()
        };
        for (sat, rows) in pos {
            let mut arr = vec![[f64::NAN; 3]; n];
            for (k, v) in rows.rows {
                let i = searchsorted(&out.t, k);
                *arr.get_mut(i).ok_or_else(|| {
                    ProductError(format!("SP3 record of {sat} past the last epoch"))
                })? = v;
            }
            out.pos.insert(sat, arr);
            let mut c = vec![f64::NAN; n];
            if let Some(rows) = clk.get(&sat) {
                for (k, v) in &rows.rows {
                    let i = searchsorted(&out.t, *k);
                    *c.get_mut(i).ok_or_else(|| {
                        ProductError(format!("SP3 clock of {sat} past the last epoch"))
                    })? = *v;
                }
            }
            out.clk_sp3.insert(sat, c);
        }
        if !clk_paths.is_empty() {
            out.clk = read_clk(clk_paths)?;
        }
        Ok(out)
    }

    /// True when the SP3 has the satellite.
    #[must_use]
    pub fn has(&self, sat: Sat) -> bool {
        self.pos.contains_key(&sat)
    }

    /// Position (m, ECEF, centre of mass), clock (s) and validity at GPS times `t`.
    #[must_use]
    pub fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        let n = t.len();
        let nt = self.t.len();
        let Some(p) = self.pos.get(&sat).filter(|_| nt > ORDER) else {
            return (vec![[0.0; 3]; n], vec![0.0; n], vec![false; n]);
        };
        let clk = self.clock(sat, t);
        let mut xyz_out = Vec::with_capacity(n);
        let mut clk_out = Vec::with_capacity(n);
        let mut ok_out = Vec::with_capacity(n);
        let np1 = ORDER + 1;
        for (i, &ti) in t.iter().enumerate() {
            let k = searchsorted(&self.t, ti) as i64;
            let i0 = (k - (np1 / 2) as i64).clamp(0, (nt - np1) as i64) as usize;
            let tn = &self.t[i0..i0 + np1];
            let pn = &p[i0..i0 + np1];
            let w = lagrange_weights(tn, ti);
            let mut wx = [0.0; ORDER + 1];
            let mut wy = [0.0; ORDER + 1];
            let mut wz = [0.0; ORDER + 1];
            for j in 0..np1 {
                // Earth rotation between the node and t, so the interpolation is in one frame
                let ang = OMGE * (tn[j] - ti);
                let (sa, ca) = (ang.sin(), ang.cos());
                let px = pn[j][0] * ca - pn[j][1] * sa;
                let py = pn[j][0] * sa + pn[j][1] * ca;
                wx[j] = w[j] * px;
                wy[j] = w[j] * py;
                wz[j] = w[j] * pn[j][2];
            }
            let xyz = [np_sum(&wx), np_sum(&wy), np_sum(&wz)];
            let finite = xyz.iter().all(|v| v.is_finite());
            let inside = ti >= self.t[0] - 1.0 && ti <= self.t[nt - 1] + 1.0 && finite;
            let c = clk[i];
            ok_out.push(inside && c.is_finite());
            xyz_out.push(xyz.map(|v| if v.is_finite() { v } else { 0.0 }));
            clk_out.push(if c.is_finite() { c } else { 0.0 });
        }
        (xyz_out, clk_out, ok_out)
    }

    fn clock(&self, sat: Sat, t: &[f64]) -> Vec<f64> {
        if let Some((ct, cv)) = self.clk.get(&sat) {
            return linear(ct, cv, t);
        }
        match self.clk_sp3.get(&sat) {
            Some(c) => linear(&self.t, c, t),
            None => vec![f64::NAN; t.len()],
        }
    }
}

/// Lagrange basis weights at `t` for the nodes `tn`, multiplied in a fixed order (node by node).
#[must_use]
pub fn lagrange_weights(tn: &[f64], t: f64) -> Vec<f64> {
    let n = tn.len();
    let mut w = vec![1.0; n];
    for j in 0..n {
        for m in 0..n {
            if m != j {
                w[j] *= (t - tn[m]) / (tn[j] - tn[m]);
            }
        }
    }
    w
}

/// Linear interpolation of a clock series: NaN across gaps over 900 s and more than 300 s
/// outside the data.
#[must_use]
pub fn linear(ts: &[f64], vs: &[f64], t: &[f64]) -> Vec<f64> {
    let (ts, vs): (Vec<f64>, Vec<f64>) = ts
        .iter()
        .zip(vs)
        .filter(|(_, v)| v.is_finite())
        .map(|(a, b)| (*a, *b))
        .unzip();
    if ts.len() < 2 {
        return vec![f64::NAN; t.len()];
    }
    let m = ts.len();
    t.iter()
        .map(|&x| {
            let k = (searchsorted(&ts, x) as i64 - 1).clamp(0, m as i64 - 2) as usize;
            let (t0, t1) = (ts[k], ts[k + 1]);
            let out = vs[k] + (vs[k + 1] - vs[k]) * (x - t0) / (t1 - t0);
            if (t1 - t0) > 900.0 || x < ts[0] - 300.0 || x > ts[m - 1] + 300.0 {
                f64::NAN
            } else {
                out
            }
        })
        .collect()
}

/// Per satellite: (times, clocks).
pub type ClockSeries = BTreeMap<Sat, (Vec<f64>, Vec<f64>)>;

/// RINEX clock files: satellite (`AS`) records -> (times, clocks) sorted by time then value.
pub fn read_clk<P: AsRef<Path>>(paths: &[P]) -> Result<ClockSeries, ProductError> {
    let mut d: BTreeMap<Sat, Vec<(f64, f64)>> = BTreeMap::new();
    for p in sort_paths(paths) {
        let text = read_text(&p)?;
        for ln in lines(&text) {
            if !ln.starts_with(b"AS ") {
                continue;
            }
            let sl = geomaticist_gnss_core::text::sl;
            let parts = split_ws(sl(ln, 8, 27));
            if parts.len() < 5 {
                continue;
            }
            let n: Option<Vec<i64>> = parts[..5].iter().map(|v| pi(v)).collect();
            let (Some(n), Some(sec), Some(v)) = (n, pf(sl(ln, 27, 34)), pf(sl(ln, 40, 59))) else {
                continue;
            };
            let Some(t) = epoch_seconds(n[0], n[1], n[2], n[3], n[4], sec) else {
                continue; // an invalid date: the record is skipped
            };
            let Some(sat) = sat_id(sl(ln, 3, 6)) else {
                continue;
            };
            d.entry(sat).or_default().push((t, v));
        }
    }
    Ok(d.into_iter()
        .map(|(s, mut rows)| {
            rows.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
            (s, rows.into_iter().unzip())
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lagrange_reproduces_a_polynomial() {
        let tn: Vec<f64> = (0..11).map(|i| f64::from(i) * 900.0).collect();
        let w = lagrange_weights(&tn, 1234.5);
        let f = |x: f64| 3.0 + 2e-3 * x - 1e-7 * x * x;
        let v: f64 = tn.iter().zip(&w).map(|(t, w)| w * f(*t)).sum();
        assert!((v - f(1234.5)).abs() < 1e-9);
        let c = linear(
            &[0.0, 300.0, 2000.0],
            &[0.0, 3.0, 5.0],
            &[150.0, 1000.0, -400.0],
        );
        assert_eq!(c[0], 1.5);
        assert!(c[1].is_nan() && c[2].is_nan());
    }
}
