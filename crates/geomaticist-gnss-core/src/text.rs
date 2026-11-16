//! Text rules the readers and writers depend on, over latin-1 bytes.
//!
//! Files are read as latin-1 (one byte = one character), so byte offsets are character offsets
//! and a column range is a clipped byte slice ([`sl`]). Stripping removes Unicode whitespace,
//! which in latin-1 includes `\x1c`..`\x1f`, `\x85` and `\xa0` ([`py_isspace`]); number fields may
//! have surrounding whitespace ([`py_float`]); rounding is half to even on the exact binary value
//! (Rust's `{:.N}` does the same); `%g`-style formatting has no Rust equivalent ([`format_g`]).

/// Unicode whitespace test for a latin-1 character.
#[must_use]
pub const fn py_isspace(b: u8) -> bool {
    matches!(b, 9..=13 | 28..=32 | 0x85 | 0xa0)
}

/// Strips whitespace ([`py_isspace`]) from both ends of latin-1 bytes.
#[must_use]
pub fn py_strip(s: &[u8]) -> &[u8] {
    let start = s.iter().position(|&b| !py_isspace(b)).unwrap_or(s.len());
    let end = s
        .iter()
        .rposition(|&b| !py_isspace(b))
        .map_or(start, |p| p + 1);
    &s[start..end.max(start)]
}

/// Strips whitespace ([`py_isspace`]) from the end.
#[must_use]
pub fn py_rstrip(s: &[u8]) -> &[u8] {
    let end = s.iter().rposition(|&b| !py_isspace(b)).map_or(0, |p| p + 1);
    &s[..end]
}

/// True when `s` is empty after stripping.
#[must_use]
pub fn is_blank(s: &[u8]) -> bool {
    s.iter().all(|&b| py_isspace(b))
}

/// Bytes `a..b`, clipped to the slice, never panics.
#[must_use]
pub fn sl(s: &[u8], a: usize, b: usize) -> &[u8] {
    let n = s.len();
    let a = a.min(n);
    let b = b.min(n).max(a);
    &s[a..b]
}

/// Latin-1 bytes to a `String` (each byte one char).
#[must_use]
pub fn latin1(s: &[u8]) -> String {
    s.iter().map(|&b| char::from(b)).collect()
}

/// A `String` back to latin-1 bytes (chars above U+00FF become `?`).
#[must_use]
pub fn to_latin1(s: &str) -> Vec<u8> {
    s.chars()
        .map(|c| u8::try_from(u32::from(c)).unwrap_or(b'?'))
        .collect()
}

/// A float from a text field: whitespace around, then a decimal or `inf` / `nan`.
/// `None` when the field is not a number. Underscores between digits are not accepted; they
/// never occur in RINEX.
#[must_use]
pub fn py_float(s: &[u8]) -> Option<f64> {
    let t = py_strip(s);
    if t.is_empty() {
        return None;
    }
    let st = std::str::from_utf8(t).ok()?;
    st.parse::<f64>().ok()
}

/// `x` rounded to an integer, half to even, as an `f64`.
#[must_use]
pub fn round_half_even(x: f64) -> f64 {
    x.round_ties_even()
}

/// `x` rounded to `ndigits` decimals for display values (half to even on the exact binary value,
/// via correctly rounded decimal formatting).
#[must_use]
pub fn round_to(x: f64, ndigits: usize) -> f64 {
    if !x.is_finite() {
        return x;
    }
    format!("{x:.ndigits$}").parse().unwrap_or(x)
}

/// `%g` formatting (6 significant digits, trailing zeros removed, exponent form below
/// 1e-4 and from 1e6).
#[must_use]
pub fn format_g(x: f64) -> String {
    format_g_prec(x, 6)
}

/// `%g` formatting with `p` significant digits.
#[must_use]
pub fn format_g_prec(x: f64, p: usize) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.into();
    }
    let p = p.max(1);
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let e = format!("{:.*e}", p - 1, x);
    let (mant, exp) = e.split_once('e').unwrap_or((&e, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let p_i = i32::try_from(p).unwrap_or(6);
    if exp < -4 || exp >= p_i {
        let m = strip_zeros(mant);
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    } else {
        let decimals = usize::try_from(p_i - 1 - exp).unwrap_or(0);
        strip_zeros(&format!("{x:.decimals$}")).to_owned()
    }
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// Lines of latin-1 bytes, with Unicode line boundaries: breaks at `\n`, `\r`, `\r\n`, `\x0b`, `\x0c`,
/// `\x1c`, `\x1d`, `\x1e`, `\x85`.
#[must_use]
pub fn py_splitlines(s: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i < s.len() {
        let b = s[i];
        if matches!(b, b'\n' | b'\r' | 0x0b | 0x0c | 0x1c | 0x1d | 0x1e | 0x85) {
            out.push(&s[start..i]);
            if b == b'\r' && s.get(i + 1) == Some(&b'\n') {
                i += 1;
            }
            start = i + 1;
        }
        i += 1;
    }
    if start < s.len() {
        out.push(&s[start..]);
    }
    out
}

/// Lines split after `\n`, each with trailing `\r`/`\n` characters removed.
pub fn lines_lf(s: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut rest = s;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let (line, next) = match rest.iter().position(|&b| b == b'\n') {
            Some(p) => (&rest[..p], &rest[p + 1..]),
            None => (rest, &b""[..]),
        };
        rest = next;
        let mut end = line.len();
        while end > 0 && (line[end - 1] == b'\r' || line[end - 1] == b'\n') {
            end -= 1;
        }
        Some(&line[..end])
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_and_slices() {
        assert_eq!(py_strip(b" \xa0ab \x85"), b"ab");
        assert_eq!(py_rstrip(b"ab  "), b"ab");
        assert_eq!(sl(b"abc", 2, 10), b"c");
        assert_eq!(sl(b"abc", 5, 10), b"");
        assert!(is_blank(b"  \t"));
        assert_eq!(latin1(b"\xe9"), "\u{e9}");
        assert_eq!(to_latin1("\u{e9}x"), b"\xe9x");
    }

    #[test]
    fn floats_like_python() {
        assert_eq!(py_float(b"  12345678.123"), Some(12_345_678.123));
        assert_eq!(py_float(b"-.5"), Some(-0.5));
        assert_eq!(py_float(b"  "), None);
        assert_eq!(py_float(b"1.0D+01"), None);
        assert_eq!(round_half_even(2.5), 2.0);
        assert_eq!(round_half_even(3.5), 4.0);
        assert_eq!(round_half_even(-2.5), -2.0);
        assert_eq!(round_half_even(2.4), 2.0);
        assert_eq!(round_to(0.25, 1), 0.2);
        assert_eq!(round_to(33.35, 1), 33.4); // 33.35 is above the tie in binary
    }

    #[test]
    fn g_format() {
        for (x, s) in [
            (30.0, "30"),
            (0.2, "0.2"),
            (1e-5, "1e-05"),
            (1_234_567.0, "1.23457e+06"),
            (123_456.0, "123456"),
            (0.0001, "0.0001"),
            (-2.5, "-2.5"),
            (1.0 / 3.0, "0.333333"),
        ] {
            assert_eq!(format_g(x), s, "{x}");
        }
    }

    #[test]
    fn line_splitting() {
        assert_eq!(py_splitlines(b"a\r\nb\rc\n"), vec![&b"a"[..], b"b", b"c"]);
        assert_eq!(py_splitlines(b"a\n\nb"), vec![&b"a"[..], b"", b"b"]);
        let v: Vec<&[u8]> = lines_lf(b"a\r\nb\rc\n\nd").collect();
        assert_eq!(v, vec![&b"a"[..], b"b\rc", b"", b"d"]);
    }
}
