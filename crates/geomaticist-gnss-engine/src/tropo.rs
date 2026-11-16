//! Troposphere: the standard atmosphere with Niell mapping, GPT3 and VMF3 (TU Wien grids) with VMF3's
//! empirical b and c, and the choice between them with its fall-backs.
//!
//! The grid readers, `gpt3`, `zhd_from_pressure`, `zwd_askne`
//! and `vmf3_site` are in `geomaticist_gnss_products::grids`. No downloads here (rule E1): the grids the
//! caller has are handed in as [`TropoInputs`]; a missing grid falls back vmf3 -> gpt3 -> niell with a
//! note.

use geomaticist_gnss_core::time::{DateTime, gps_datetime, gps_seconds};
use geomaticist_gnss_products::grids::{
    Gpt3Grid, gpt3, season, vmf3_site, zhd_from_pressure, zwd_askne,
};

use crate::arr::F2;
use crate::models::{NMF_HT, mapf, niell_at, niell_coeffs, saastamoinen_zenith};
use crate::nputil::np_interp;
use crate::vmf3_coeffs as vc;

/// The troposphere models.
pub const MODELS: [&str; 3] = ["niell", "gpt3", "vmf3"];

/// One VMF3 6-hourly grid with its epoch (GPS seconds of the 0/6/12/18 h UTC boundary as
/// `(ep - gps_datetime(0)).total_seconds()` gives it).
#[derive(Debug, Clone)]
pub struct Vmf3Grid {
    /// Epoch (GPS seconds of the grid's datetime).
    pub t: f64,
    /// The 64800 rows (lat lon ah aw zhd zwd).
    pub rows: Vec<[f64; 6]>,
}

/// The grids a caller has (loaded by the app or the CLI). `None` means "not available".
#[derive(Debug, Clone, Default)]
pub struct TropoInputs<'a> {
    /// GPT3 1° grid.
    pub gpt3: Option<&'a Gpt3Grid>,
    /// VMF3 orography (64800 heights).
    pub orography: Option<&'a [f64]>,
    /// VMF3 grids by epoch.
    pub vmf3: &'a [Vmf3Grid],
}

/// The 6-hourly grid epochs (GPS seconds of the UTC-less datetimes, computed on GPS
/// datetimes) that bracket `[start, end]`.
#[must_use]
pub fn vmf3_epochs(start_gps: f64, end_gps: f64) -> Vec<(DateTime, f64)> {
    let t0 = gps_datetime(start_gps);
    let t1 = gps_datetime(end_gps);
    let mut e = DateTime {
        year: t0.year,
        month: t0.month,
        day: t0.day,
        hour: t0.hour / 6 * 6,
        minute: 0,
        second: 0,
        micro: 0,
    };
    let mut out = vec![(e, gps_seconds(&e))];
    while e < t1 {
        e = DateTime::from_micros_since_gps0(e.micros_since_gps0() + 6 * 3600 * 1_000_000);
        out.push((e, gps_seconds(&e)));
    }
    out
}

/// The grid file name of an epoch.
#[must_use]
pub fn vmf3_name(ep: &DateTime) -> String {
    format!(
        "VMF3_{:04}{:02}{:02}.H{:02}",
        ep.year, ep.month, ep.day, ep.hour
    )
}

/// Day of year with its fraction (microseconds ignored).
#[must_use]
pub fn doy_of(t_gps: f64) -> f64 {
    let d = gps_datetime(t_gps);
    d.yday() as f64 + (d.hour as f64 + d.minute as f64 / 60.0 + d.second as f64 / 3600.0) / 24.0
}

const NMAX: usize = 12;

/// Unnormalised associated Legendre terms V (cos part) and W (sin part), n, m ≤ 12, in the tables' order
/// (`(n, m)` for n in 0..=12, m in 0..=n).
#[must_use]
pub fn legendre(lat: f64, lon: f64) -> (Vec<f64>, Vec<f64>) {
    let pd = std::f64::consts::FRAC_PI_2 - lat;
    let (x, y, z) = (pd.sin() * lon.cos(), pd.sin() * lon.sin(), pd.cos());
    let mut v = [[0.0f64; NMAX + 1]; NMAX + 1];
    let mut w = [[0.0f64; NMAX + 1]; NMAX + 1];
    v[0][0] = 1.0;
    v[1][0] = z;
    for n in 2..=NMAX {
        v[n][0] = ((2 * n - 1) as f64 * z * v[n - 1][0] - (n - 1) as f64 * v[n - 2][0]) / n as f64;
    }
    for m in 1..=NMAX {
        v[m][m] = (2 * m - 1) as f64 * (x * v[m - 1][m - 1] - y * w[m - 1][m - 1]);
        w[m][m] = (2 * m - 1) as f64 * (x * w[m - 1][m - 1] + y * v[m - 1][m - 1]);
        if m < NMAX {
            v[m + 1][m] = (2 * m + 1) as f64 * z * v[m][m];
            w[m + 1][m] = (2 * m + 1) as f64 * z * w[m][m];
        }
        for n in m + 2..=NMAX {
            v[n][m] = ((2 * n - 1) as f64 * z * v[n - 1][m] - (n + m - 1) as f64 * v[n - 2][m])
                / (n - m) as f64;
            w[n][m] = ((2 * n - 1) as f64 * z * w[n - 1][m] - (n + m - 1) as f64 * w[n - 2][m])
                / (n - m) as f64;
        }
    }
    let mut vv = Vec::with_capacity(91);
    let mut ww = Vec::with_capacity(91);
    for n in 0..=NMAX {
        for m in 0..=n {
            vv.push(v[n][m]);
            ww.push(w[n][m]);
        }
    }
    (vv, ww)
}

fn sh_eval(v: &[f64], w: &[f64], a: &[[f64; 5]; 91], b: &[[f64; 5]; 91], s: &[f64; 5]) -> f64 {
    // (v @ a + w @ b) @ s
    let mut col = [0.0; 5];
    for c in 0..5 {
        let mut va = 0.0;
        let mut wb = 0.0;
        for k in 0..91 {
            va += v[k] * a[k][c];
            wb += w[k] * b[k][c];
        }
        col[c] = va + wb;
    }
    col.iter().zip(s).map(|(x, y)| x * y).sum()
}

/// VMF3's empirical `(bh, bw, ch, cw)` at a site and day.
#[must_use]
pub fn vmf3_bc(lat: f64, lon: f64, doy: f64) -> (f64, f64, f64, f64) {
    let (v, w) = legendre(lat, lon);
    let s = season(doy);
    (
        sh_eval(&v, &w, &vc::ANM_BH, &vc::BNM_BH, &s),
        sh_eval(&v, &w, &vc::ANM_BW, &vc::BNM_BW, &s),
        sh_eval(&v, &w, &vc::ANM_CH, &vc::BNM_CH, &s),
        sh_eval(&v, &w, &vc::ANM_CW, &vc::BNM_CW, &s),
    )
}

/// VMF3 hydrostatic (lifted to height `h` by Niell's correction) and wet mapping at one elevation (rad).
#[inline]
#[must_use]
pub fn vmf3_map_at(el: f64, ah: f64, aw: f64, bc: (f64, f64, f64, f64), h: f64) -> (f64, f64) {
    let (bh, bw, ch, cw) = bc;
    let el = el.max(1e-3);
    let mh =
        mapf(el, ah, bh, ch) + (1.0 / el.sin() - mapf(el, NMF_HT.0, NMF_HT.1, NMF_HT.2)) * h / 1e3;
    (mh, mapf(el, aw, bw, cw))
}

/// [`vmf3_map_at`] over a slice of elevations.
#[must_use]
pub fn vmf3_map(
    el: &[f64],
    ah: f64,
    aw: f64,
    bc: (f64, f64, f64, f64),
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut mh = Vec::with_capacity(el.len());
    let mut mw = Vec::with_capacity(el.len());
    for &e in el {
        let (a, b) = vmf3_map_at(e, ah, aw, bc, h);
        mh.push(a);
        mw.push(b);
    }
    (mh, mw)
}

/// Hydrostatic delay and mapping functions for one site.
#[derive(Debug, Clone, PartialEq)]
pub struct Tropo {
    /// Zenith hydrostatic delay (m) per epoch (one value repeated for niell / gpt3; `zhd_per_epoch`
    /// says whether it varies).
    pub zhd: Vec<f64>,
    /// Whether `zhd` is a per-epoch series (vmf3) or a scalar (niell, gpt3).
    pub zhd_per_epoch: bool,
    /// Hydrostatic mapping, the shape of `el`.
    pub mh: F2,
    /// Wet mapping.
    pub mw: F2,
    /// A-priori zenith wet delay (m), 0 for niell.
    pub zwd0: f64,
    /// A-priori total gradients (north, east), m.
    pub grad: (f64, f64),
    /// The model actually used.
    pub used: &'static str,
    /// Why another model than asked.
    pub note: String,
}

/// The VMF3 series `(ah, aw, zhd, zwd)` per epoch, linear in time between the grids; `None` when any grid
/// of the span is missing.
fn vmf3_series(t: &[f64], lat: f64, lon: f64, h: f64, inp: &TropoInputs) -> Option<[Vec<f64>; 4]> {
    let orog = inp.orography?;
    let eps = vmf3_epochs(t[0], t[t.len() - 1]);
    let mut vals: Vec<[f64; 4]> = Vec::with_capacity(eps.len());
    let mut tt: Vec<f64> = Vec::with_capacity(eps.len());
    for (_ep, te) in &eps {
        let g = inp.vmf3.iter().find(|g| (g.t - te).abs() < 1e-6)?;
        vals.push(vmf3_site(&g.rows, orog, lat, lon, h));
        tt.push(*te);
    }
    let n = t.len();
    let mut out: [Vec<f64>; 4] = [vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    for k in 0..4 {
        if tt.len() == 1 {
            out[k] = vec![vals[0][k]; n];
        } else {
            let fp: Vec<f64> = vals.iter().map(|v| v[k]).collect();
            out[k] = t.iter().map(|&x| np_interp(x, &tt, &fp)).collect();
        }
    }
    Some(out)
}

/// Hydrostatic delay and mapping functions for one site: lat / lon (rad), h (m), epochs `t`, elevations
/// `el` `(ne, ns)`. `model`: `niell | gpt3 | vmf3` (falls back vmf3 -> gpt3 -> niell).
#[must_use]
pub fn troposphere(
    model: &str,
    lat: f64,
    lon: f64,
    h: f64,
    t: &[f64],
    el: &F2,
    inp: &TropoInputs,
) -> Tropo {
    let (ne, ns) = (el.ne, el.ns);
    let mut note = String::new();
    let doy = doy_of(0.5 * (t[0] + t[t.len() - 1]));
    let mut model = model;
    if model == "vmf3" {
        if let Some([ah, aw, zhd, zwd]) = vmf3_series(t, lat, lon, h, inp) {
            let bc = vmf3_bc(lat, lon, doy);
            let mut mh = F2::zeros(ne, ns);
            let mut mw = F2::zeros(ne, ns);
            for e in 0..ne {
                for j in 0..ns {
                    let (a, b) = vmf3_map_at(el.at(e, j), ah[e], aw[e], bc, h);
                    mh.set(e, j, a);
                    mw.set(e, j, b);
                }
            }
            let grad = inp.gpt3.map_or((0.0, 0.0), |g| {
                let g3 = gpt3(g, lat, lon, h, doy);
                (g3.gn_h + g3.gn_w, g3.ge_h + g3.ge_w)
            });
            return Tropo {
                zhd,
                zhd_per_epoch: true,
                mh,
                mw,
                zwd0: crate::nputil::mean(&zwd),
                grad,
                used: "vmf3",
                note,
            };
        }
        note = "VMF3 grids missing for these days: GPT3 used".to_owned();
        model = "gpt3";
    }
    if model == "gpt3" {
        if let Some(grid) = inp.gpt3 {
            let g3 = gpt3(grid, lat, lon, h, doy);
            let bc = vmf3_bc(lat, lon, doy);
            let mut mh = F2::zeros(ne, ns);
            let mut mw = F2::zeros(ne, ns);
            for e in 0..ne {
                for j in 0..ns {
                    let (a, b) = vmf3_map_at(el.at(e, j), g3.ah, g3.aw, bc, h);
                    mh.set(e, j, a);
                    mw.set(e, j, b);
                }
            }
            return Tropo {
                zhd: vec![zhd_from_pressure(g3.p, lat, h); ne],
                zhd_per_epoch: false,
                mh,
                mw,
                zwd0: zwd_askne(g3.e, g3.Tm, g3.la),
                grad: (g3.gn_h + g3.gn_w, g3.ge_h + g3.ge_w),
                used: "gpt3",
                note,
            };
        }
        if !note.is_empty() {
            note.push_str("; ");
        }
        note.push_str("GPT3 grid missing: standard atmosphere and Niell used");
    }
    let (zhd, _) = saastamoinen_zenith(lat, h, 0.0);
    let (ah, aw) = niell_coeffs(lat, doy.trunc());
    let mut mh = F2::zeros(ne, ns);
    let mut mw = F2::zeros(ne, ns);
    for e in 0..ne {
        for j in 0..ns {
            let (a, b) = niell_at(el.at(e, j), &ah, &aw, h);
            mh.set(e, j, a);
            mw.set(e, j, b);
        }
    }
    Tropo {
        zhd: vec![zhd; ne],
        zhd_per_epoch: false,
        mh,
        mw,
        zwd0: 0.0,
        grad: (0.0, 0.0),
        used: "niell",
        note,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vmf3_coefficients_are_sane() {
        let (bh, bw, ch, cw) = vmf3_bc(1.1, -2.5, 267.0);
        // reference values at the same site and day: 0.002732529980686709, 0.0013134027544609063,
        // 0.05857145391694364, 0.0384059208694693
        assert!((bh - 0.002_732_529_980_686_709).abs() < 1e-15, "{bh}");
        assert!((bw - 0.001_313_402_754_460_906_3).abs() < 1e-15, "{bw}");
        assert!((ch - 0.058_571_453_916_943_64).abs() < 1e-14, "{ch}");
        assert!((cw - 0.038_405_920_869_469_3).abs() < 1e-14, "{cw}");
        let (mh, mw) = vmf3_map(&[0.2, 1.0], 0.00124, 0.00058, (bh, bw, ch, cw), 300.0);
        assert!(mh[0] > mh[1] && mw[0] > mw[1]);
        let eps = vmf3_epochs(1_411_149_600.0, 1_411_153_199.0);
        assert_eq!(eps.len(), 2);
        assert_eq!(vmf3_name(&eps[0].0), "VMF3_20240923.H18");
        assert_eq!(vmf3_name(&eps[1].0), "VMF3_20240924.H00");
        let t: Vec<f64> = (0..12)
            .map(|k| 1_411_149_600.0 + 300.0 * k as f64)
            .collect();
        let el = F2::filled(12, 2, 0.5);
        let tr = troposphere(
            "vmf3",
            1.134,
            -2.574,
            319.0,
            &t,
            &el,
            &TropoInputs::default(),
        );
        assert_eq!(tr.used, "niell");
        assert!(tr.note.contains("GPT3 grid missing"));
        assert!((doy_of(1_411_149_600.0) - (267.0 + 0.75)).abs() < 1e-9);
    }
}
