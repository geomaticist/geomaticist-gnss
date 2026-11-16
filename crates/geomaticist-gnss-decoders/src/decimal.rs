//! Exact decimal formatting of doubles and rationals.
//!
//! Every RINEX value the decoders print goes through one of these rules:
//! - the exact binary value of the double, rounded with ties away from zero. [`Dec::exact`] gives the exact
//!   value, [`Dec::round_at`] the rounding.
//! - 15 significant digits (`%.15g`): Rust's correctly rounded `{:e}` formatting (ties to even on the exact
//!   value); [`Dec::parse`] reads it back.
//! - rationals with a power-of-two denominator are rounded in integers ([`round_ratio`]).

use std::fmt::Write as _;

/// A decimal number `±digits × 10^exp` (digits without leading zeros; `"0"` for zero).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Dec {
    /// Sign (kept for zero, as `Decimal` keeps `-0`).
    pub neg: bool,
    /// ASCII digits, most significant first.
    pub digits: Vec<u8>,
    /// Power of ten of the last digit.
    pub exp: i32,
}

impl Dec {
    /// The exact value of a finite double.
    #[must_use]
    pub fn exact(v: f64) -> Self {
        let neg = v.is_sign_negative();
        let bits = v.abs().to_bits();
        let e_bits = i32::try_from((bits >> 52) & 0x7ff).unwrap_or(0);
        let frac = bits & ((1u64 << 52) - 1);
        let (mut m, mut e) = if e_bits == 0 {
            (frac, -1074)
        } else {
            (frac | (1u64 << 52), e_bits - 1075)
        };
        if m == 0 {
            return Self {
                neg,
                digits: b"0".to_vec(),
                exp: 0,
            };
        }
        while m & 1 == 0 {
            m >>= 1;
            e += 1;
        }
        // value = m * 2^e
        let mut big = Big::from_u64(m);

        let exp10 = if e >= 0 {
            for _ in 0..e {
                big.mul_small(2);
            }
            0
        } else {
            let k = -e;
            let mut left = k;
            while left >= 13 {
                big.mul_small(1_220_703_125); // 5^13
                left -= 13;
            }
            if left > 0 {
                big.mul_small(5u64.pow(u32::try_from(left).unwrap_or(0)));
            }
            -k
        };
        let mut digits = big.to_digits();
        let mut exp = exp10;
        while digits.len() > 1 && digits.last() == Some(&b'0') {
            digits.pop();
            exp += 1;
        }
        Self { neg, digits, exp }
    }

    /// Parses a plain or exponent decimal literal (`"-12.5"`, `"1.5e+07"`).
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        let (neg, body) = match s.as_bytes().first()? {
            b'-' => (true, &s[1..]),
            b'+' => (false, &s[1..]),
            _ => (false, s),
        };
        let (mant, e) = match body.find(['e', 'E']) {
            Some(i) => (&body[..i], body[i + 1..].parse::<i32>().ok()?),
            None => (body, 0),
        };
        let (int, frac) = mant.split_once('.').unwrap_or((mant, ""));
        let mut digits: Vec<u8> = int.bytes().chain(frac.bytes()).collect();
        if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
            return None;
        }
        let exp = e - i32::try_from(frac.len()).ok()?;
        let lead = digits.iter().take_while(|&&d| d == b'0').count();
        digits.drain(..lead.min(digits.len() - 1));
        Some(Self { neg, digits, exp })
    }

    /// True for zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.digits == b"0"
    }

    /// `Decimal.adjusted()`: the power of ten of the most significant digit.
    #[must_use]
    pub fn adjusted(&self) -> i32 {
        self.exp + i32::try_from(self.digits.len()).unwrap_or(0) - 1
    }

    /// Rounded (ties away from zero, `ROUND_HALF_UP`) so the last digit is at 10^`q`
    /// (`quantize(Decimal(10) ** q, ROUND_HALF_UP)`).
    #[must_use]
    pub fn round_at(&self, q: i32) -> Self {
        let mut digits = self.digits.clone();
        if self.exp >= q {
            if !self.is_zero() {
                digits.extend(std::iter::repeat_n(b'0', (self.exp - q) as usize));
            }
            return Self {
                neg: self.neg,
                digits,
                exp: q,
            };
        }
        let drop = (q - self.exp) as usize;

        let up = if drop > digits.len() {
            digits.clear();
            false
        } else {
            let first = digits[digits.len() - drop];
            digits.truncate(digits.len() - drop);
            first >= b'5'
        };
        if up {
            let mut i = digits.len();
            loop {
                if i == 0 {
                    digits.insert(0, b'1');
                    break;
                }
                i -= 1;
                if digits[i] == b'9' {
                    digits[i] = b'0';
                } else {
                    digits[i] += 1;
                    break;
                }
            }
        }
        let lead = digits.iter().take_while(|&&d| d == b'0').count();
        digits.drain(..lead);
        if digits.is_empty() {
            digits.push(b'0');
        }
        Self {
            neg: self.neg,
            digits,
            exp: q,
        }
    }

    /// Rounded to `n` significant digits, ties away from zero.
    #[must_use]
    pub fn round_sig(&self, n: i32) -> Self {
        self.round_at(self.adjusted() - (n - 1))
    }

    /// The absolute value printed with exactly `-exp` decimals (`exp <= 0`), as `f"{abs(d):.Nf}"` prints a
    /// quantized `Decimal`; for `exp > 0` the zeros are written out.
    #[must_use]
    pub fn abs_fixed(&self) -> String {
        let s: String = self.digits.iter().map(|&b| b as char).collect();
        if self.exp >= 0 {
            let mut out = s;
            if !self.is_zero() {
                out.extend(std::iter::repeat_n('0', self.exp as usize));
            }
            return out;
        }
        let places = (-self.exp) as usize;
        if s.len() > places {
            let (a, b) = s.split_at(s.len() - places);
            format!("{a}.{b}")
        } else {
            format!("0.{}{s}", "0".repeat(places - s.len()))
        }
    }

    /// Shifted by 10^`k` (`scaleb`).
    #[must_use]
    pub fn scaleb(&self, k: i32) -> Self {
        let mut d = self.clone();
        if !d.is_zero() {
            d.exp += k;
        }
        d
    }
}

/// Little-endian base-10^9 big integer, just enough for [`Dec::exact`].
struct Big(Vec<u64>);

impl Big {
    const BASE: u64 = 1_000_000_000;

    fn from_u64(mut v: u64) -> Self {
        let mut limbs = Vec::new();
        while v > 0 {
            limbs.push(v % Self::BASE);
            v /= Self::BASE;
        }
        Self(limbs)
    }

    fn mul_small(&mut self, k: u64) {
        let mut carry = 0u64;
        for l in &mut self.0 {
            let x = *l * k + carry;
            *l = x % Self::BASE;
            carry = x / Self::BASE;
        }
        while carry > 0 {
            self.0.push(carry % Self::BASE);
            carry /= Self::BASE;
        }
    }

    fn to_digits(&self) -> Vec<u8> {
        let mut s = String::new();
        for (i, l) in self.0.iter().rev().enumerate() {
            if i == 0 {
                let _ = write!(s, "{l}");
            } else {
                let _ = write!(s, "{l:09}");
            }
        }
        if s.is_empty() {
            s.push('0');
        }
        s.into_bytes()
    }
}

/// `|num / den|` rounded to `places` decimals with ties away from zero, as text with the sign
/// (`-` only when the rounded value is not zero).
#[must_use]
pub fn round_ratio(num: i128, den: i128, places: u32) -> String {
    let scale = 10i128.pow(places);
    let n = num.abs() * scale;
    let d = den.abs();
    // k = int(n/d + 1/2) = floor((2n + d) / 2d)
    let k = (2 * n + d) / (2 * d);
    let neg = ((num < 0) != (den < 0)) && k != 0;
    let whole = k / scale;
    let frac = k % scale;
    let w = places as usize;
    if places == 0 {
        format!("{}{whole}", if neg { "-" } else { "" })
    } else {
        format!("{}{whole}.{frac:0w$}", if neg { "-" } else { "" })
    }
}

/// `v` at 15 significant digits (`%.15g`) read back as a decimal (the value the writers start from).
#[must_use]
pub fn g15(v: f64) -> Dec {
    // `{:.14e}` is correctly rounded with ties to even, as '%.15g' is; the plain/exponent layout of
    // '%g' does not change the value.
    let s = format!("{v:.14e}");
    Dec::parse(&s).unwrap_or_else(|| Dec::exact(v))
}

/// `str(Decimal(f"{v:.15g}").quantize(Decimal("0.001"), ROUND_HALF_UP))`: the T04 pseudorange text.
#[must_use]
pub fn round_float_3(v: f64) -> String {
    let d = g15(v).round_at(-3);
    let body = d.abs_fixed();
    if d.neg { format!("-{body}") } else { body }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_values() {
        assert_eq!(Dec::exact(0.5).abs_fixed(), "0.5");
        assert_eq!(Dec::exact(0.1).digits.len(), 55);
        assert_eq!(Dec::exact(1e22).abs_fixed(), "10000000000000000000000");
        assert_eq!(Dec::exact(3.0).abs_fixed(), "3");
        assert!(Dec::exact(-0.0).neg);
    }

    #[test]
    fn rounding_half_up() {
        let d = Dec::parse("2.5").unwrap();
        assert_eq!(d.round_at(0).abs_fixed(), "3");
        let d = Dec::parse("9.9996").unwrap();
        assert_eq!(d.round_at(-3).abs_fixed(), "10.000");
        let d = Dec::parse("0.0004").unwrap();
        assert_eq!(d.round_at(-3).abs_fixed(), "0.000");
        assert_eq!(Dec::exact(1234.5678).round_sig(3).abs_fixed(), "1230");
    }

    #[test]
    fn ratios() {
        assert_eq!(round_ratio(-1, 2000, 3), "-0.001"); // -0.0005 -> away from zero
        assert_eq!(round_ratio(-1, 4000, 3), "0.000");
        assert_eq!(round_ratio(12345, 256, 3), "48.223");
    }

    #[test]
    fn pseudorange_text() {
        assert_eq!(round_float_3(21_234_567.123_449_9), "21234567.123");
        assert_eq!(round_float_3(20_000_000.000_5), "20000000.001");
        // ties to even at 15 digits first ('%.15g'), then half up
        assert_eq!(format!("{:.0}", 2.5_f64), "2");
    }
}
