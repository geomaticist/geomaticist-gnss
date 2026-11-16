//! NumPy's reductions with NumPy's operation order (pairwise sums), so the engine reproduces NumPy's
//! means, medians and standard deviations to the last bit where that is possible.
//!
//! `np_sum` is the products crate's pairwise sum; `median` is `np.median` (the mean of the two middle
//! values for an even count), `nanmedian` ignores NaN, `std_ddof` is `np.std(ddof=…)` as `_var` computes
//! it (mean, squared deviations, pairwise sum, divide), `var` likewise with ddof 0.

pub use geomaticist_gnss_products::util::np_sum;

/// `np.mean` of a 1-D array (NaN for an empty one, as NumPy warns and returns NaN).
#[must_use]
pub fn mean(a: &[f64]) -> f64 {
    if a.is_empty() {
        return f64::NAN;
    }
    np_sum(a) / a.len() as f64
}

/// `np.var(a, ddof)`: squared deviations from the mean, pairwise-summed, over `n - ddof`.
#[must_use]
pub fn var_ddof(a: &[f64], ddof: usize) -> f64 {
    let n = a.len();
    if n <= ddof {
        return f64::NAN;
    }
    let m = mean(a);
    let sq: Vec<f64> = a
        .iter()
        .map(|&x| {
            let d = x - m;
            d * d
        })
        .collect();
    np_sum(&sq) / (n - ddof) as f64
}

/// `np.std(a, ddof)`.
#[must_use]
pub fn std_ddof(a: &[f64], ddof: usize) -> f64 {
    var_ddof(a, ddof).sqrt()
}

/// `np.median` of a 1-D array (NaN when empty; NaN values propagate as NumPy's sort puts them last and the
/// median then reads NaN — callers filter first).
#[must_use]
pub fn median(a: &[f64]) -> f64 {
    let n = a.len();
    if n == 0 {
        return f64::NAN;
    }
    if a.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }
    let mut v = a.to_vec();
    v.sort_by(f64::total_cmp);
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// `np.nanmedian`: the median of the finite values (NaN when there are none).
#[must_use]
pub fn nanmedian(a: &[f64]) -> f64 {
    let v: Vec<f64> = a.iter().copied().filter(|x| !x.is_nan()).collect();
    median(&v)
}

/// `1.4826 * median(|a - median(a)|)` (the robust scatter the slip detector and the screens use).
#[must_use]
pub fn mad_scatter(a: &[f64]) -> f64 {
    let m = median(a);
    let d: Vec<f64> = a.iter().map(|x| (x - m).abs()).collect();
    1.4826 * median(&d)
}

/// `np.round(x, decimals)`: `rint(x * 10^d) / 10^d`, half to even.
#[must_use]
pub fn np_round(x: f64, decimals: i32) -> f64 {
    if decimals == 0 {
        return x.round_ties_even();
    }
    let s = 10f64.powi(decimals);
    (x * s).round_ties_even() / s
}

/// `x` rounded half to even, as an `i64`; 0 for non-finite input.
#[must_use]
pub fn py_round_i(x: f64) -> i64 {
    if x.is_finite() {
        x.round_ties_even() as i64
    } else {
        0
    }
}

/// `np.searchsorted(a, v)` (left): first index with `a[i] >= v`.
#[must_use]
pub fn searchsorted(a: &[f64], v: f64) -> usize {
    a.partition_point(|&x| x < v)
}

/// `np.interp` on one value (ends clamped).
pub use geomaticist_gnss_products::util::np_interp;

/// `np.hypot`.
#[must_use]
pub fn hypot(a: f64, b: f64) -> f64 {
    a.hypot(b)
}

/// The `np.argmax` of a slice (first maximum; NaN never wins unless all are NaN, where NumPy returns the
/// first NaN: here index 0).
#[must_use]
pub fn argmax(a: &[f64]) -> usize {
    let mut k = 0;
    let mut best = f64::NEG_INFINITY;
    let mut found = false;
    for (i, &x) in a.iter().enumerate() {
        if x.is_nan() {
            // NumPy: the first NaN wins argmax
            return i;
        }
        if !found || x > best {
            best = x;
            k = i;
            found = true;
        }
    }
    k
}

/// `np.argmin` (first minimum).
#[must_use]
pub fn argmin(a: &[f64]) -> usize {
    let mut k = 0;
    let mut best = f64::INFINITY;
    let mut found = false;
    for (i, &x) in a.iter().enumerate() {
        if x.is_nan() {
            return i;
        }
        if !found || x < best {
            best = x;
            k = i;
            found = true;
        }
    }
    k
}

/// `np.linspace(0, n - 1, k).astype(int)` made unique: the sample rows the orbit checks use.
#[must_use]
pub fn sample_rows(n: usize, k: usize) -> Vec<usize> {
    let k = k.min(n);
    if k == 0 {
        return Vec::new();
    }
    let mut out: Vec<usize> = Vec::with_capacity(k);
    // np.linspace: arange(num) * step + start with step = (stop - start) / (num - 1); the last value is set
    // to stop exactly; astype(int) truncates; np.unique drops repeats (the values are increasing).
    let step = if k == 1 {
        0.0
    } else {
        (n - 1) as f64 / (k - 1) as f64
    };
    for i in 0..k {
        let idx = if i + 1 == k {
            n - 1
        } else {
            (step * i as f64).trunc() as usize
        };
        if out.last() != Some(&idx) {
            out.push(idx);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reductions() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[3.0, 1.0, 2.0, 4.0]), 2.5);
        assert_eq!(nanmedian(&[f64::NAN, 1.0, 3.0]), 2.0);
        assert!((std_ddof(&[1.0, 2.0, 3.0, 4.0], 1) - 1.290_994_448_735_805_6).abs() < 1e-15);
        assert_eq!(np_round(2.5, 0), 2.0);
        assert_eq!(np_round(1.0005, 3), 1.0);
        assert_eq!(py_round_i(-2.5), -2);
        assert_eq!(argmax(&[1.0, 3.0, 3.0]), 1);
        assert_eq!(sample_rows(10, 240).len(), 10);
        assert_eq!(sample_rows(1000, 240).len(), 240);
        assert_eq!(sample_rows(5, 3), vec![0, 2, 4]);
    }
}
