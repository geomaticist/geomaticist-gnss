//! Hatanaka compact RINEX (CRINEX 1.0 for RINEX 2, CRINEX 3.0 for RINEX 3/4) decoder.
//!
//! Written from the published format description: Y. Hatanaka (2008), "A Compression Format and
//! Tools for GNSS Observation Data", Bulletin of the Geographical Survey Institute 55, and the
//! format notes distributed with RNXCMP. This decoder gives the same observation text as the C
//! `crx2rnx` (checked byte for byte on CORS files).
//!
//! Format, as decoded here:
//! - two CRINEX header lines (`CRINEX VERS   / TYPE`, `CRINEX PROG / DATE`), then the RINEX header
//!   unchanged up to END OF HEADER;
//! - per epoch: the epoch line as a text difference against the previous one (a space keeps the old
//!   character, `&` writes a space, anything else replaces it; a line starting with `&` in 1.0 or
//!   `>` in 3.0 is a full line), with every satellite id appended (1.0: from column 33 without
//!   continuation lines; 3.0: from column 42);
//! - event epochs (flag 2-5): the epoch line and its records as written;
//! - the receiver clock line: empty, `n&value` (a new arc of difference order n) or a difference;
//! - per satellite one line: the observation fields separated by one space, each empty (missing),
//!   `n&value` (new arc) or an n-th order difference of the value × 1000; after the last field the
//!   LLI / signal strength characters as a text difference against the satellite's previous ones.

use std::collections::HashMap;

use geomaticist_gnss_core::text::{lines_lf, py_strip, sl};

use crate::RinexError;
use crate::header::{HeaderLine, RinexHeader};

/// A differencing arc: `u[0]` is the latest value, `u[k]` its k-th difference.
#[derive(Debug, Clone)]
struct Arc {
    order: usize,
    max_order: usize,
    u: [i64; 10],
}

impl Arc {
    fn init(max_order: usize, value: i64) -> Self {
        let mut u = [0; 10];
        u[0] = value;
        Self {
            order: 0,
            max_order: max_order.min(9),
            u,
        }
    }

    fn step(&mut self, d: i64) -> i64 {
        if self.order < self.max_order {
            self.order += 1;
        }
        let k = self.order;
        self.u[k] = d;
        for j in (1..=k).rev() {
            self.u[j - 1] = self.u[j - 1].wrapping_add(self.u[j]);
        }
        self.u[0]
    }
}

fn err(msg: impl Into<String>) -> RinexError {
    RinexError::new(format!("CRINEX: {}", msg.into()))
}

/// Parses one compressed field: `None` = blank, `Some(Ok((Some(order), v)))` = new arc,
/// `Some(Ok((None, d)))` = difference.
fn parse_field(f: &[u8]) -> Result<Option<(Option<usize>, i64)>, RinexError> {
    if f.is_empty() {
        return Ok(None);
    }
    if let Some(p) = f.iter().position(|&b| b == b'&') {
        let order = parse_int(&f[..p])
            .and_then(|o| usize::try_from(o).ok())
            .ok_or_else(|| err("bad arc order"))?;
        let v = parse_int(&f[p + 1..]).ok_or_else(|| bad_number(&f[p + 1..]))?;
        Ok(Some((Some(order), v)))
    } else {
        Ok(Some((None, parse_int(f).ok_or_else(|| bad_number(f))?)))
    }
}

fn bad_number(s: &[u8]) -> RinexError {
    err(format!("bad number {:?}", String::from_utf8_lossy(s)))
}

/// A signed decimal integer with optional surrounding blanks (no other characters).
fn parse_int(s: &[u8]) -> Option<i64> {
    let s = py_strip(s);
    let (neg, digits) = match s.first()? {
        b'-' => (true, &s[1..]),
        b'+' => (false, &s[1..]),
        _ => (false, s),
    };
    if digits.is_empty() || digits.len() > 18 {
        return None;
    }
    let mut v: i64 = 0;
    for &c in digits {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v * 10 + i64::from(c - b'0');
    }
    Some(if neg { -v } else { v })
}

/// Text difference: a space keeps the old character, `&` writes a space, others replace; a longer
/// new text extends the old.
fn repair(old: &mut Vec<u8>, diff: &[u8]) {
    for (i, &d) in diff.iter().enumerate() {
        if i < old.len() {
            match d {
                b' ' => {}
                b'&' => old[i] = b' ',
                c => old[i] = c,
            }
        } else {
            old.push(if d == b'&' { b' ' } else { d });
        }
    }
}

/// `value / 10^dec` as fixed-point text with `dec` decimals, right-justified to `width`, as
/// crx2rnx prints it: no "0" before the point (`.621`, `-.000123456789`).
fn fixed(v: i64, dec: u32, width: usize) -> String {
    let mut out = Vec::with_capacity(width);
    push_fixed(&mut out, v, dec, width);
    String::from_utf8_lossy(&out).into_owned()
}

/// [`fixed`] appended to a buffer, without allocating.
fn push_fixed(out: &mut Vec<u8>, v: i64, dec: u32, width: usize) {
    let mut tmp = [0u8; 24];
    let mut i = tmp.len();
    let mut a = v.unsigned_abs();
    for _ in 0..dec {
        i -= 1;
        tmp[i] = b'0' + (a % 10) as u8;
        a /= 10;
    }
    i -= 1;
    tmp[i] = b'.';
    while a > 0 {
        i -= 1;
        tmp[i] = b'0' + (a % 10) as u8;
        a /= 10;
    }
    if v < 0 {
        i -= 1;
        tmp[i] = b'-';
    }
    let n = tmp.len() - i;
    out.extend(std::iter::repeat_n(b' ', width.saturating_sub(n)));
    out.extend_from_slice(&tmp[i..]);
}

/// One satellite of an epoch: id, values (x 1000) and LLI / SSI characters.
type SatRecord = ([u8; 3], Vec<Option<i64>>, Vec<u8>);

#[derive(Debug, Clone, Default)]
struct SatState {
    fields: Vec<Option<Arc>>,
    flags: Vec<u8>,
}

fn push_line(out: &mut Vec<u8>, line: &[u8]) {
    let mut end = line.len();
    while end > 0 && line[end - 1] == b' ' {
        end -= 1;
    }
    out.extend_from_slice(&line[..end]);
    out.push(b'\n');
}

/// Decodes CRINEX text to plain RINEX text (LF line ends, trailing blanks trimmed).
pub fn decompress(input: &[u8]) -> Result<Vec<u8>, RinexError> {
    let mut lines = lines_lf(input);
    let l1 = lines.next().ok_or_else(|| err("empty file"))?;
    if !sl(l1, 60, 80).starts_with(b"CRINEX VERS") {
        return Err(err("not a CRINEX file (no 'CRINEX VERS / TYPE' line)"));
    }
    let crx_ver: f64 = std::str::from_utf8(py_strip(sl(l1, 0, 9)))
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| err("unreadable CRINEX version"))?;
    let v3 = crx_ver >= 3.0;
    lines
        .next()
        .ok_or_else(|| err("no CRINEX PROG / DATE line"))?;

    let mut out: Vec<u8> = Vec::with_capacity(input.len() * 4);
    let mut hlines = Vec::new();
    for line in lines.by_ref() {
        out.extend_from_slice(line);
        out.push(b'\n');
        let hl = HeaderLine::from_raw(line);
        let end = hl.is("END OF HEADER");
        hlines.push(hl);
        if end {
            break;
        }
    }
    let hdr = RinexHeader::new(hlines)?;
    let ntypes = |sys: u8| -> usize {
        if v3 {
            hdr.types_for(sys).len()
        } else {
            hdr.v2_types.len()
        }
    };

    let mut epoch: Vec<u8> = Vec::new();
    let mut clock: Option<Arc> = None;
    let mut sats: HashMap<[u8; 3], SatState> = HashMap::new();
    let (flag_col, n_col, sat_col) = if v3 {
        (31, 32..35, 41)
    } else {
        (28, 29..32, 32)
    };

    while let Some(line) = lines.next() {
        if line.is_empty() && epoch.is_empty() {
            continue;
        }
        let init = if v3 {
            line.first() == Some(&b'>')
        } else {
            line.first() == Some(&b'&')
        };
        if init {
            epoch = line.to_vec();
            if !v3 {
                epoch[0] = b' ';
            }
        } else {
            repair(&mut epoch, line);
        }
        let flag = sl(&epoch, flag_col, flag_col + 1)
            .first()
            .copied()
            .unwrap_or(b'0');
        let n: usize = std::str::from_utf8(sl(&epoch, n_col.start, n_col.end))
            .ok()
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0);
        if (b'2'..=b'5').contains(&flag) {
            // event: the epoch line (without satellites) and its records as written
            push_line(&mut out, sl(&epoch, 0, if v3 { 35 } else { 32 }));
            for _ in 0..n {
                let rec = lines
                    .next()
                    .ok_or_else(|| err("file ends inside an event"))?;
                push_line(&mut out, rec);
            }
            continue;
        }
        let ids: Vec<[u8; 3]> = (0..n)
            .map(|k| {
                let s = sl(&epoch, sat_col + 3 * k, sat_col + 3 * k + 3);
                let mut id = [b' '; 3];
                id[..s.len()].copy_from_slice(s);
                id
            })
            .collect();
        // clock line
        let cl = lines
            .next()
            .ok_or_else(|| err("file ends before the clock line"))?;
        let clock_val = match parse_field(py_strip(cl))? {
            None => {
                clock = None;
                None
            }
            Some((Some(order), v)) => {
                clock = Some(Arc::init(order, v));
                Some(v)
            }
            Some((None, d)) => match clock.as_mut() {
                Some(a) => Some(a.step(d)),
                None => return Err(err("clock difference without an arc")),
            },
        };
        let mut new_sats: HashMap<[u8; 3], SatState> = HashMap::with_capacity(n);
        let mut recs: Vec<SatRecord> = Vec::with_capacity(n);
        for id in &ids {
            let dl = lines
                .next()
                .ok_or_else(|| err("file ends inside an epoch"))?;
            let nt = ntypes(id[0]).max(usize::from(!v3));
            let mut st = sats.remove(id).unwrap_or_default();
            st.fields.resize(nt, None);
            let mut vals = Vec::with_capacity(nt);
            let mut p = 0usize;
            for j in 0..nt {
                let start = p.min(dl.len());
                let end = dl[start..]
                    .iter()
                    .position(|&b| b == b' ')
                    .map_or(dl.len(), |q| start + q);
                let f = &dl[start..end];
                p = if end < dl.len() { end + 1 } else { dl.len() };
                let v = match parse_field(f)? {
                    None => {
                        st.fields[j] = None;
                        None
                    }
                    Some((Some(order), v)) => {
                        st.fields[j] = Some(Arc::init(order, v));
                        Some(v)
                    }
                    Some((None, d)) => match st.fields[j].as_mut() {
                        Some(a) => Some(a.step(d)),
                        None => {
                            return Err(err(format!(
                                "difference without an arc ({})",
                                String::from_utf8_lossy(id)
                            )));
                        }
                    },
                };
                vals.push(v);
            }
            let flags_diff = &dl[p.min(dl.len())..];
            repair(&mut st.flags, flags_diff);
            let flags = st.flags.clone();
            // crx2rnx blanks the stored flags of a missing field, so the next diff is taken
            // against blanks (a stale LLI never comes back)
            for (j, v) in vals.iter().enumerate() {
                if v.is_none() {
                    for k in [2 * j, 2 * j + 1] {
                        if let Some(c) = st.flags.get_mut(k) {
                            *c = b' ';
                        }
                    }
                }
            }
            new_sats.insert(*id, st);
            recs.push((*id, vals, flags));
        }
        sats = new_sats;
        // RINEX text of the epoch
        if v3 {
            let mut el = sl(&epoch, 0, 35).to_vec();
            el.resize(35, b' ');
            if let Some(c) = clock_val {
                el.extend_from_slice(b"      ");
                el.extend_from_slice(fixed(c, 12, 15).as_bytes());
            }
            push_line(&mut out, &el);
            for (id, vals, flags) in &recs {
                let mut l = id.to_vec();
                for (j, v) in vals.iter().enumerate() {
                    push_field(&mut l, *v, flags, j);
                }
                push_line(&mut out, &l);
            }
        } else {
            let mut first = sl(&epoch, 0, 32).to_vec();
            first.resize(32, b' ');
            for (k, chunk) in ids.chunks(12).enumerate() {
                let mut l = if k == 0 {
                    first.clone()
                } else {
                    vec![b' '; 32]
                };
                for id in chunk {
                    l.extend_from_slice(id);
                }
                if k == 0
                    && let Some(c) = clock_val
                {
                    l.resize(68, b' ');
                    l.extend_from_slice(fixed(c, 9, 12).as_bytes());
                }
                push_line(&mut out, &l);
            }
            if ids.is_empty() {
                push_line(&mut out, &first);
            }
            for (_id, vals, flags) in &recs {
                let mut l = Vec::with_capacity(80);
                for (j, v) in vals.iter().enumerate() {
                    push_field(&mut l, *v, flags, j);
                    if j % 5 == 4 {
                        push_line(&mut out, &l);
                        l.clear();
                    }
                }
                if vals.len() % 5 != 0 || vals.is_empty() {
                    push_line(&mut out, &l);
                }
            }
        }
    }
    Ok(out)
}

fn push_field(l: &mut Vec<u8>, v: Option<i64>, flags: &[u8], j: usize) {
    // crx2rnx: a value prints without a leading "0" (`.621`, `-.985`), and a missing value
    // prints as 16 blanks whatever its (stale) flags hold
    match v {
        Some(x) => {
            push_fixed(l, x, 3, 14);
            l.push(*flags.get(2 * j).unwrap_or(&b' '));
            l.push(*flags.get(2 * j + 1).unwrap_or(&b' '));
        }
        None => l.extend_from_slice(&[b' '; 16]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arcs_reconstruct() {
        let xs = [100i64, 107, 115, 126, 140, 150];
        let mut a = Arc::init(3, xs[0]);
        // differences as the compressor writes them: order ramps up to 3
        let mut hist = vec![xs[0]];
        for &x in &xs[1..] {
            hist.push(x);
            let k = (hist.len() - 1).min(3);
            let mut d: Vec<i64> = hist[hist.len() - 1 - k..].to_vec();
            for _ in 0..k {
                d = d.windows(2).map(|w| w[1] - w[0]).collect();
            }
            assert_eq!(a.step(d[0]), x);
        }
    }

    #[test]
    fn text_repair_and_fixed() {
        let mut s = b"> 2024 09 23 00 00  0.0000000  0  3".to_vec();
        let diff = format!("{}30{}2", " ".repeat(19), " ".repeat(13));
        repair(&mut s, diff.as_bytes());
        assert_eq!(&s, b"> 2024 09 23 00 00 30.0000000  0  2");
        assert_eq!(fixed(-500, 3, 14), "         -.500");
        assert_eq!(fixed(123_456_789, 3, 14), "    123456.789");
        assert_eq!(fixed(-1_234, 12, 15), " -.000000001234");
        assert_eq!(fixed(621, 3, 14), "          .621");
        assert_eq!(fixed(-985, 3, 14), "         -.985");
        assert_eq!(fixed(-1_327, 3, 14), "        -1.327");
    }

    #[test]
    fn decodes_small_crx3() {
        let crx = "\
3.0                 COMPACT RINEX FORMAT                    CRINEX VERS   / TYPE
RNX2CRX ver.4.1.0                       04-Nov-26 00:00     CRINEX PROG / DATE
     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    2 C1C L1C                                              SYS / # / OBS TYPES
                                                            END OF HEADER
> 2024 09 23 00 00  0.0000000  0  1      G01

3&20000000123 3&100000000123 1 6
                   30
3&5
100 5000 &
";
        let txt = decompress(crx.as_bytes()).unwrap();
        let s = String::from_utf8(txt).unwrap();
        let body: Vec<&str> = s.lines().skip(3).collect();
        assert_eq!(body[0], "> 2024 09 23 00 00  0.0000000  0  1");
        assert_eq!(body[1], "G01  20000000.1231  100000000.1236");
        assert_eq!(
            body[2],
            "> 2024 09 23 00 00 30.0000000  0  1        .000000000005"
        );
        assert_eq!(body[3], "G01  20000000.223   100000005.1236");
    }
}
