//! Shared helpers: reading text files, and NumPy-style numerics where the result
//! depends on the exact operation order (`np.interp`, `np.sum`, `np.arange`, `np.mod`, `//`).

use std::io::Read;
use std::path::Path;

use geomaticist_gnss_core::text::{py_float, py_isspace};

/// A product reading error (file and reason).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductError(pub String);

impl std::fmt::Display for ProductError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ProductError {}

impl ProductError {
    /// An error with this message.
    pub fn new(msg: impl Into<String>) -> Self {
        Self(msg.into())
    }
}

/// A file's bytes, gunzipped when the name ends in `.gz`, with
/// universal newlines (`\r\n` and `\r` become `\n`).
pub fn read_text(path: &Path) -> Result<Vec<u8>, ProductError> {
    let raw = std::fs::read(path).map_err(|e| ProductError(format!("{}: {e}", path.display())))?;
    let gz = path
        .file_name()
        .is_some_and(|n| n.to_string_lossy().to_ascii_lowercase().ends_with(".gz"));
    let data = if gz {
        let mut out = Vec::with_capacity(raw.len() * 4);
        flate2::read::MultiGzDecoder::new(raw.as_slice())
            .read_to_end(&mut out)
            .map_err(|e| ProductError(format!("{}: {e}", path.display())))?;
        out
    } else {
        raw
    };
    Ok(universal_newlines(&data))
}

/// `\r\n` and lone `\r` to `\n`.
#[must_use]
pub fn universal_newlines(data: &[u8]) -> Vec<u8> {
    if !data.contains(&b'\r') {
        return data.to_vec();
    }
    let mut out = Vec::with_capacity(data.len());
    let mut i = 0;
    while i < data.len() {
        if data[i] == b'\r' {
            out.push(b'\n');
            if data.get(i + 1) == Some(&b'\n') {
                i += 1;
            }
        } else {
            out.push(data[i]);
        }
        i += 1;
    }
    out
}

/// Lines split at `\n` (text already in universal-newline form), without the `\n`.
pub fn lines(data: &[u8]) -> impl Iterator<Item = &[u8]> {
    let end = data.strip_suffix(b"\n").unwrap_or(data);
    end.split(|&b| b == b'\n')
        .take(if data.is_empty() { 0 } else { usize::MAX })
}

/// Splits latin-1 bytes on runs of whitespace (empty parts dropped).
#[must_use]
pub fn split_ws(s: &[u8]) -> Vec<&[u8]> {
    s.split(|&b| py_isspace(b))
        .filter(|p| !p.is_empty())
        .collect()
}

/// A float from text, `None` when it is not a valid number.
#[must_use]
pub fn pf(s: &[u8]) -> Option<f64> {
    py_float(s)
}

/// An integer from text (whitespace around, optional sign, digits).
#[must_use]
pub fn pi(s: &[u8]) -> Option<i64> {
    geomaticist_gnss_core::time::py_int(s)
}

/// `np.interp(x, xp, fp)` for one `x` (`xp` increasing), as NumPy computes it:
/// `slope * (x - xp[j]) + fp[j]`, ends clamped.
#[must_use]
pub fn np_interp(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return fp[0];
    }
    if x.is_nan() {
        return x;
    }
    if x < xp[0] {
        return fp[0];
    }
    if x > xp[n - 1] {
        return fp[n - 1];
    }
    if x == xp[n - 1] {
        return fp[n - 1];
    }
    // j with xp[j] <= x < xp[j + 1]
    let j = xp.partition_point(|&v| v <= x).saturating_sub(1).min(n - 2);
    if xp[j] == x {
        return fp[j]; // NumPy: "avoid potential non-finite interpolation"
    }
    let slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    let r = slope * (x - xp[j]) + fp[j];
    if r.is_nan() {
        let r2 = slope * (x - xp[j + 1]) + fp[j + 1];
        if r2.is_nan() && fp[j] == fp[j + 1] {
            return fp[j];
        }
        return r2;
    }
    r
}

/// NumPy's pairwise summation (`np.add.reduce` over a contiguous axis).
#[must_use]
pub fn np_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - n % 8 {
            for k in 0..8 {
                r[k] += a[i + k];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_sum(&a[..n2]) + np_sum(&a[n2..])
    }
}

/// `np.mean` of a 1-D array.
#[must_use]
pub fn np_mean(a: &[f64]) -> f64 {
    np_sum(a) / a.len() as f64
}

/// `np.arange(start, stop, step)` for floats: `ceil((stop - start) / step)` values
/// `start + i * delta` with `delta = (start + step) - start` (NumPy's `DOUBLE_fill`).
#[must_use]
pub fn np_arange(start: f64, stop: f64, step: f64) -> Vec<f64> {
    let len = ((stop - start) / step).ceil();
    if !(len.is_finite() && len > 0.0) {
        return Vec::new();
    }
    let len = len as usize;
    let mut out = Vec::with_capacity(len);
    out.push(start);
    if len > 1 {
        let delta = (start + step) - start;
        out.push(start + delta);
        for i in 2..len {
            out.push(start + i as f64 * delta);
        }
    }
    out
}

/// NumPy's `divmod` for floats: (floor quotient, modulo with the divisor's sign).
#[must_use]
pub fn np_divmod(a: f64, b: f64) -> (f64, f64) {
    let mut m = a % b;
    if b == 0.0 {
        return (a / b, m);
    }
    let mut div = (a - m) / b;
    if m == 0.0 {
        m = 0.0_f64.copysign(b);
    } else if (b < 0.0) != (m < 0.0) {
        m += b;
        div -= 1.0;
    }
    let floordiv = if div == 0.0 {
        0.0_f64.copysign(a / b)
    } else {
        let mut f = div.floor();
        if div - f > 0.5 {
            f += 1.0;
        }
        f
    };
    (floordiv, m)
}

/// Float modulo, floored (result has the divisor's sign).
#[must_use]
pub fn py_mod(a: f64, b: f64) -> f64 {
    let mut m = a % b;
    if m == 0.0 {
        m = 0.0_f64.copysign(b);
    } else if (b < 0.0) != (m < 0.0) {
        m += b;
    }
    m
}

/// `np.searchsorted(a, v)` (side `left`).
#[must_use]
pub fn searchsorted(a: &[f64], v: f64) -> usize {
    a.partition_point(|&x| x < v)
}

/// Whether one latin-1 character is a letter.
#[must_use]
pub fn is_alpha_latin1(c: u8) -> bool {
    c.is_ascii_alphabetic()
        || matches!(c, 0xAA | 0xB5 | 0xBA | 0xC0..=0xD6 | 0xD8..=0xF6 | 0xF8..=0xFF)
}

/// Sorts paths in the platform's path order (Windows: case-insensitive,
/// `/` as `\`).
pub fn sort_paths<P: AsRef<Path>>(paths: &[P]) -> Vec<std::path::PathBuf> {
    let mut v: Vec<std::path::PathBuf> = paths.iter().map(|p| p.as_ref().to_path_buf()).collect();
    if cfg!(windows) {
        v.sort_by_key(|p| p.to_string_lossy().replace('/', "\\").to_lowercase());
    } else {
        v.sort();
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numpy_numerics() {
        let xp = [0.0, 5.0, 10.0];
        let fp = [1.0, 2.0, 4.0];
        assert_eq!(np_interp(2.5, &xp, &fp), 1.5);
        assert_eq!(np_interp(-1.0, &xp, &fp), 1.0);
        assert_eq!(np_interp(10.0, &xp, &fp), 4.0);
        assert_eq!(np_interp(7.5, &xp, &fp), 3.0);
        assert_eq!(np_sum(&[1.0; 11]), 11.0);
        assert_eq!(
            np_sum(&(0..300).map(f64::from).collect::<Vec<_>>()),
            44_850.0
        );
        assert_eq!(np_arange(0.0, 92.5, 5.0).len(), 19);
        assert_eq!(np_divmod(-30.0, 360.0), (-1.0, 330.0));
        assert_eq!(np_divmod(725.0, 360.0), (2.0, 5.0));
        assert_eq!(py_mod(-147.72, 360.0), 360.0 - 147.72);
        let t = universal_newlines(b"a\r\nb\rc\n");
        assert_eq!(lines(&t).collect::<Vec<_>>(), vec![&b"a"[..], b"b", b"c"]);
        assert_eq!(lines(b"").count(), 0);
    }
}
