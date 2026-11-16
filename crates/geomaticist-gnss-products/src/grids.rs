//! TU Wien 1° grids: GPT3 (`gpt3_1.grd`), daily VMF3 (`VMF3_YYYYMMDD.Hhh`) and the orography
//! (`orography_ell_1x1`), with the GPT3 site evaluation, the four-cell bilinear weights and the
//! VMF3 site values.
//!
//! GPT3 is kept in single precision (f32); every value is widened back to f64 before use. Row order: latitude 89.5 down to -89.5, longitude 0.5 to 359.5.

use std::path::Path;

use crate::util::{ProductError, py_mod, read_text, split_ws};

/// Rows of every TU Wien 1° grid.
pub const ROWS: usize = 64_800;
/// GPT3 columns.
pub const GPT3_COLS: usize = 64;

/// Standard gravity (m/s²).
pub const G_M: f64 = 9.80665;
/// Molar mass of dry air (kg/mol).
pub const M_DRY: f64 = 28.965e-3;
/// Gas constant (J/K/mol).
pub const R_GAS: f64 = 8.3143;

/// A numeric text grid: whitespace-separated numbers, `c` to the end of a line ignored,
/// blank lines skipped; every row must have `cols` values (0 = any). Returns the values row-major.
fn loadtxt(
    path: &Path,
    comment: Option<u8>,
    cols: usize,
) -> Result<(Vec<f64>, usize), ProductError> {
    let text = read_text(path)?;
    let mut out = Vec::new();
    let mut nrows = 0;
    let mut width = cols;
    for (k, line) in text.split(|&b| b == b'\n').enumerate() {
        let line = match comment.and_then(|c| line.iter().position(|&b| b == c)) {
            Some(p) => &line[..p],
            None => line,
        };
        let parts = split_ws(line);
        if parts.is_empty() {
            continue;
        }
        if width == 0 {
            width = parts.len();
        }
        if parts.len() != width {
            return Err(ProductError(format!(
                "{}: line {} has {} values, expected {width}",
                path.display(),
                k + 1,
                parts.len()
            )));
        }
        for p in parts {
            let s = std::str::from_utf8(p).map_err(|_| bad_value(path, k))?;
            out.push(s.parse::<f64>().map_err(|_| bad_value(path, k))?);
        }
        nrows += 1;
    }
    Ok((out, nrows))
}

fn bad_value(path: &Path, k: usize) -> ProductError {
    ProductError(format!(
        "{}: unreadable number on line {}",
        path.display(),
        k + 1
    ))
}

/// The GPT3 1° grid, 64800 x 64, in single precision.
#[derive(Debug, Clone)]
pub struct Gpt3Grid {
    data: Vec<f32>,
}

impl Gpt3Grid {
    /// Reads `gpt3_1.grd`.
    pub fn open(path: &Path) -> Result<Self, ProductError> {
        let (v, rows) = loadtxt(path, Some(b'%'), 0)?;
        if rows != ROWS || v.len() != ROWS * GPT3_COLS {
            return Err(ProductError(format!(
                "{}: ({rows}, {}), expected 64800 x 64",
                path.display(),
                v.len() / rows.max(1)
            )));
        }
        Ok(Self {
            data: v.into_iter().map(|x| x as f32).collect(),
        })
    }

    /// Row `r`, widened to f64.
    #[must_use]
    pub fn row(&self, r: usize) -> [f64; GPT3_COLS] {
        let mut out = [0.0; GPT3_COLS];
        for (o, v) in out
            .iter_mut()
            .zip(&self.data[r * GPT3_COLS..(r + 1) * GPT3_COLS])
        {
            *o = f64::from(*v);
        }
        out
    }
}

/// A daily VMF3 grid (lat lon ah aw zhd zwd), 64800 x 6.
pub fn load_vmf3(path: &Path) -> Result<Vec<[f64; 6]>, ProductError> {
    let (v, rows) = loadtxt(path, Some(b'!'), 0)?;
    if rows != ROWS || v.len() != ROWS * 6 {
        return Err(ProductError(format!(
            "{}: ({rows}, {}), expected 64800 x 6",
            path.display(),
            v.len() / rows.max(1)
        )));
    }
    Ok(v.as_chunks::<6>()
        .0
        .iter()
        .map(|c| [c[0], c[1], c[2], c[3], c[4], c[5]])
        .collect())
}

/// The VMF3 orography (64800 heights, m).
pub fn load_orography(path: &Path) -> Result<Vec<f64>, ProductError> {
    let text = read_text(path)?;
    let mut out = Vec::with_capacity(ROWS);
    for p in split_ws(&text) {
        let s = std::str::from_utf8(p).map_err(|_| bad_value(path, 0))?;
        out.push(s.parse::<f64>().map_err(|_| bad_value(path, 0))?);
    }
    if out.len() != ROWS {
        return Err(ProductError(format!(
            "{}: {} values, expected 64800",
            path.display(),
            out.len()
        )));
    }
    Ok(out)
}

/// `VMF3_YYYYMMDD.Hhh`.
#[must_use]
pub fn vmf3_name(year: i64, month: i64, day: i64, hour: i64) -> String {
    format!("VMF3_{year:04}{month:02}{day:02}.H{hour:02}")
}

/// Rows of the four cells around (lat, lon) (radians) and their bilinear weights; one row near the
/// poles.
#[must_use]
pub fn cells(lat: f64, lon: f64) -> Vec<(usize, f64)> {
    let pod = 90.0 - lat.to_degrees();
    let plon = py_mod(lon.to_degrees(), 360.0);
    let ip = (pod.floor() as i64).min(179);
    let il = (plon.floor() as i64).rem_euclid(360);
    if pod <= 0.5 || pod >= 179.5 {
        return vec![((ip * 360 + il) as usize, 1.0)];
    }
    let dp = pod - (ip as f64 + 0.5);
    let dl = plon - (il as f64 + 0.5);
    let sign = |x: f64| {
        if x > 0.0 {
            1
        } else if x < 0.0 {
            -1
        } else {
            0
        }
    };
    let ip1 = (ip + sign(dp)).clamp(0, 179);
    let il1 = (il + sign(dl)).rem_euclid(360);
    let (a, b) = (dp.abs(), dl.abs());
    vec![
        ((ip * 360 + il) as usize, (1.0 - a) * (1.0 - b)),
        ((ip1 * 360 + il) as usize, a * (1.0 - b)),
        ((ip * 360 + il1) as usize, (1.0 - a) * b),
        ((ip1 * 360 + il1) as usize, a * b),
    ]
}

/// The annual and semi-annual terms `[1, cos w, sin w, cos 2w, sin 2w]`.
#[must_use]
pub fn season(doy: f64) -> [f64; 5] {
    let w = 2.0 * std::f64::consts::PI * doy / 365.25;
    [1.0, w.cos(), w.sin(), (2.0 * w).cos(), (2.0 * w).sin()]
}

/// GPT3 at one site.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(non_snake_case)] // the usual symbol names (T, dT, Tm)
pub struct Gpt3 {
    /// Pressure (hPa).
    pub p: f64,
    /// Temperature (°C).
    pub T: f64,
    /// Lapse rate (K/km).
    pub dT: f64,
    /// Water-vapour weighted mean temperature (K).
    pub Tm: f64,
    /// Water vapour pressure (hPa).
    pub e: f64,
    /// VMF3 hydrostatic a at zero height.
    pub ah: f64,
    /// VMF3 wet a.
    pub aw: f64,
    /// Water vapour decrease factor.
    pub la: f64,
    /// Geoid undulation of the grid (m).
    pub undu: f64,
    /// Hydrostatic north gradient (m).
    pub gn_h: f64,
    /// Hydrostatic east gradient (m).
    pub ge_h: f64,
    /// Wet north gradient (m).
    pub gn_w: f64,
    /// Wet east gradient (m).
    pub ge_w: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// GPT3 at lat / lon (radians), ellipsoidal height `h` (m) and day of year (with fraction).
#[must_use]
pub fn gpt3(grid: &Gpt3Grid, lat: f64, lon: f64, h: f64, doy: f64) -> Gpt3 {
    let cs = cells(lat, lon);
    let w: Vec<f64> = cs.iter().map(|c| c.1).collect();
    let g: Vec<[f64; GPT3_COLS]> = cs.iter().map(|c| grid.row(c.0)).collect();
    let s = season(doy);
    let sea = |c0: usize, scale: f64| -> Vec<f64> {
        g.iter().map(|r| dot(&r[c0..c0 + 5], &s) / scale).collect()
    };
    let (p0, t0, q, dt) = (sea(2, 1.0), sea(7, 1.0), sea(12, 1e3), sea(17, 1e3));
    let n = g.len();
    let mut tl = vec![0.0; n];
    let mut pl = vec![0.0; n];
    let mut el = vec![0.0; n];
    let la = sea(34, 1.0);
    let undu: Vec<f64> = g.iter().map(|r| r[22]).collect();
    for k in 0..n {
        let redh = (h - undu[k]) - g[k][23];
        tl[k] = t0[k] + dt[k] * redh - 273.15;
        let tv = t0[k] * (1.0 + 0.6077 * q[k]);
        pl[k] = p0[k] * (-G_M * M_DRY / (R_GAS * tv) * redh).exp() / 100.0;
        let e0 = q[k] * p0[k] / (0.622 + 0.378 * q[k]) / 100.0;
        el[k] = e0 * (100.0 * pl[k] / p0[k]).powf(la[k] + 1.0);
    }
    Gpt3 {
        p: dot(&w, &pl),
        T: dot(&w, &tl),
        dT: dot(&w, &dt) * 1e3,
        Tm: dot(&w, &sea(39, 1.0)),
        e: dot(&w, &el),
        ah: dot(&w, &sea(24, 1e3)),
        aw: dot(&w, &sea(29, 1e3)),
        la: dot(&w, &la),
        undu: dot(&w, &undu),
        gn_h: dot(&w, &sea(44, 1e5)),
        ge_h: dot(&w, &sea(49, 1e5)),
        gn_w: dot(&w, &sea(54, 1e5)),
        ge_w: dot(&w, &sea(59, 1e5)),
    }
}

/// Zenith hydrostatic delay (m) from the site pressure (hPa): Saastamoinen (Davis et al. 1985).
#[must_use]
pub fn zhd_from_pressure(p: f64, lat: f64, h: f64) -> f64 {
    0.002_276_8 * p / (1.0 - 0.002_66 * (2.0 * lat).cos() - 0.28e-6 * h)
}

/// A-priori zenith wet delay (m), Askne & Nordius (1987).
#[must_use]
pub fn zwd_askne(e: f64, tm: f64, la: f64) -> f64 {
    let (k1, k2, k3) = (77.604, 64.79, 377_600.0);
    let k2p = k2 - k1 * 18.0152 / 28.9644;
    1e-6 * (k2p + k3 / tm) * (R_GAS / M_DRY) / ((la + 1.0) * G_M) * e
}

/// ah, aw, zhd, zwd of one VMF3 grid at a site (delays lifted to `h`, then bilinear).
#[must_use]
pub fn vmf3_site(grid: &[[f64; 6]], orog: &[f64], lat: f64, lon: f64, h: f64) -> [f64; 4] {
    let cs = cells(lat, lon);
    let c = 1.0 - 0.002_66 * (2.0 * lat).cos();
    let mut out = [0.0; 4];
    for &(r, w) in &cs {
        let g = grid[r];
        let ho = orog[r];
        let p = g[4] / 0.002_276_8 * (c - 0.28e-6 * ho);
        let p = p * (1.0 - 0.000_022_6 * (h - ho)).powf(5.225);
        let zhd = 0.002_276_8 * p / (c - 0.28e-6 * h);
        let zwd = g[5] * (-(h - ho) / 2000.0).exp();
        out[0] += w * g[2];
        out[1] += w * g[3];
        out[2] += w * zhd;
        out[3] += w * zwd;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cells_and_weights() {
        let c = cells(64.978_f64.to_radians(), (-147.499_f64).to_radians());
        assert_eq!(c.len(), 4);
        let s: f64 = c.iter().map(|x| x.1).sum();
        assert!((s - 1.0).abs() < 1e-12);
        assert_eq!(c[0].0, 25 * 360 + 212);
        assert_eq!(cells(89.9_f64.to_radians(), 0.0).len(), 1);
        assert_eq!(vmf3_name(2024, 9, 23, 6), "VMF3_20240923.H06");
        let z = zhd_from_pressure(1013.25, 0.0, 0.0);
        assert!((z - 2.3125).abs() < 1e-3);
    }
}
