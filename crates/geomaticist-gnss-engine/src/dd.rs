//! Double differences: reference satellite per system and epoch, the DD operator and its covariance, and
//! the clock-free residual the screening uses.
//!
//! The engine itself eliminates the clock per (epoch, system) instead of forming DDs; these are kept for the
//! PPP crate and the tests.

/// Reference satellite among `cols` (one system, one epoch): the highest; with `prev` and `keep_deg > 0`,
/// `prev` is kept while within `keep_deg` of the highest.
#[must_use]
pub fn reference(cols: &[usize], el: &[f64], prev: Option<usize>, keep_deg: f64) -> usize {
    let mut best = cols[0];
    let mut be = f64::NEG_INFINITY;
    for &c in cols {
        if el[c] > be {
            be = el[c];
            best = c;
        }
    }
    if let Some(p) = prev
        && keep_deg > 0.0
        && cols.contains(&p)
        && el[p] >= el[best] - keep_deg.to_radians()
    {
        return p;
    }
    best
}

/// `D (n-1, n)`: `DD = D · SD`, rows for every satellite but the reference.
#[must_use]
pub fn operator(n: usize, ref_pos: usize) -> Vec<Vec<f64>> {
    (0..n)
        .filter(|&i| i != ref_pos)
        .map(|i| {
            let mut row = vec![0.0; n];
            row[i] = 1.0;
            row[ref_pos] = -1.0;
            row
        })
        .collect()
}

/// DD covariance `D diag(var) Dᵀ`.
#[must_use]
pub fn cov(var_sd: &[f64], ref_pos: usize) -> Vec<Vec<f64>> {
    let v: Vec<f64> = var_sd
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != ref_pos)
        .map(|(_, x)| *x)
        .collect();
    let vr = var_sd[ref_pos];
    (0..v.len())
        .map(|i| {
            (0..v.len())
                .map(|j| if i == j { v[i] + vr } else { vr })
                .collect()
        })
        .collect()
}

/// Inverse of the DD covariance (Sherman-Morrison).
#[must_use]
pub fn weight(var_sd: &[f64], ref_pos: usize) -> Vec<Vec<f64>> {
    let v: Vec<f64> = var_sd
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != ref_pos)
        .map(|(_, x)| *x)
        .collect();
    let vr = var_sd[ref_pos];
    let di: Vec<f64> = v.iter().map(|x| 1.0 / x).collect();
    let s: f64 = di.iter().sum();
    let f = vr / (1.0 + vr * s);
    (0..v.len())
        .map(|i| {
            (0..v.len())
                .map(|j| if i == j { di[i] } else { 0.0 } - di[i] * di[j] * f)
                .collect()
        })
        .collect()
}

/// SD residuals with the epoch's weighted mean removed.
#[must_use]
pub fn clock_free(res_sd: &[f64], var_sd: &[f64]) -> Vec<f64> {
    let w: Vec<f64> = var_sd.iter().map(|v| 1.0 / v).collect();
    let sw: f64 = w.iter().sum();
    let m: f64 = w.iter().zip(res_sd).map(|(a, b)| a * b).sum::<f64>() / sw;
    res_sd.iter().map(|r| r - m).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn operator_cov_weight() {
        let var = [1.0, 2.0, 4.0];
        let d = operator(3, 1);
        assert_eq!(d, vec![vec![1.0, -1.0, 0.0], vec![0.0, -1.0, 1.0]]);
        let c = cov(&var, 1);
        let w = weight(&var, 1);
        for i in 0..2 {
            for j in 0..2 {
                let v: f64 = (0..2).map(|k| c[i][k] * w[k][j]).sum();
                assert!((v - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        assert_eq!(reference(&[0, 1, 2], &[0.3, 0.9, 0.5], Some(2), 30.0), 2);
        assert_eq!(reference(&[0, 1, 2], &[0.3, 0.9, 0.5], Some(2), 1.0), 1);
        let cf = clock_free(&[1.0, 2.0, 3.0], &[1.0, 1.0, 1.0]);
        assert!((cf[0] + 1.0).abs() < 1e-15);
    }
}
