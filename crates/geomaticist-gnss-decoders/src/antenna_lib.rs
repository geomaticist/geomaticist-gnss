//! Antennas: Trimble antenna code -> IGS name and radome, and each antenna's measure-up methods.
//!
//! A method's offsets: horizontal (for lever / slant marks) and vertical from the measured point to the phase
//! center (APC). The "ARP method" (usually *Bottom of antenna mount*) is the antenna reference point, so
//! height to ARP = height to APC - that method's vertical offset.

use crate::antennas::{ANTENNAS, MEASURE_METHODS, Method};

/// `(IGS name, dome, manufacturer, model name)` of a code.
#[must_use]
pub fn antenna(code: i64) -> Option<(&'static str, &'static str, &'static str, &'static str)> {
    ANTENNAS
        .binary_search_by_key(&code, |e| e.0)
        .ok()
        .map(|i| ANTENNAS[i].1)
}

/// `(index of the ARP method, methods)` of a code.
#[must_use]
pub fn measure_methods(code: i64) -> Option<(usize, &'static [Method])> {
    MEASURE_METHODS
        .binary_search_by_key(&code, |e| e.0)
        .ok()
        .map(|i| MEASURE_METHODS[i].1)
}

/// The measure-up methods of a code, empty if unknown.
#[must_use]
pub fn methods(code: i64) -> &'static [Method] {
    measure_methods(code).map_or(&[], |m| m.1)
}

/// Their names.
#[must_use]
pub fn method_names(code: i64) -> Vec<&'static str> {
    methods(code).iter().map(|m| m.0).collect()
}

/// Model name (`"R10-2 Internal"`), empty if unknown.
#[must_use]
pub fn model_name(code: i64) -> &'static str {
    antenna(code).map_or("", |a| a.3)
}

/// 20-character RINEX antenna name: IGS name (16) + radome (4), e.g. `"TRMR10-2        NONE"`; empty if unknown.
#[must_use]
pub fn igs_name(code: i64) -> String {
    antenna(code).map_or_else(String::new, |(igs, dome, _, _)| {
        format!(
            "{igs:<16}{:>4}",
            if dome.is_empty() { "NONE" } else { dome }
        )
    })
}

/// Antenna code from a code (`"549"`, `"code 549"`), a Trimble model name (`"R10-2 Internal"`), an IGS name, or
/// `"IGS DOME"`.
#[must_use]
pub fn code_for(antenna_text: &str) -> Option<i64> {
    let s = antenna_text.trim();
    if s.is_empty() {
        return None;
    }
    let lower = s.to_lowercase();
    let num = lower
        .strip_prefix("code")
        .map_or(lower.as_str(), str::trim_start);
    if !num.is_empty() && num.bytes().all(|b| b.is_ascii_digit()) {
        let c: i64 = num.parse().ok()?;
        return antenna(c).map(|_| c);
    }
    let collapsed = lower.split_whitespace().collect::<Vec<_>>().join(" ");
    for key in [lower.as_str(), collapsed.as_str()] {
        for &(code, (igs, dome, _mf, name)) in ANTENNAS {
            let cands = [name.to_owned(), igs.to_owned(), format!("{igs} {dome}")];
            if cands
                .iter()
                .any(|k| !k.is_empty() && k.trim().to_lowercase() == key)
            {
                return Some(code);
            }
        }
    }
    None
}

/// Measured height with a measure-up method -> (height to ARP, height to APC), metres.
///
/// Vertical-to-mark = sqrt((measured + tape)^2 - horizontal^2) for slant/lever methods, else measured;
/// APC = that + the method's vertical offset; ARP = APC - the ARP method's vertical offset.
/// `None` when the code or method is unknown.
#[must_use]
pub fn reduce_height(code: i64, method: usize, measured: f64) -> Option<(f64, f64)> {
    let (arp_i, ms) = measure_methods(code)?;
    let &(_name, h, v, tape) = ms.get(method)?;
    let m = measured + tape;
    let vert = if h == 0.0 {
        m
    } else {
        (m * m - h * h).max(0.0).sqrt()
    };
    let apc = vert + v;
    Some((apc - ms.get(arp_i)?.2, apc))
}

/// The same with the method given by name (case-insensitive).
#[must_use]
pub fn reduce_height_named(code: i64, method: &str, measured: f64) -> Option<(f64, f64)> {
    let want = method.trim().to_lowercase();
    let i = methods(code)
        .iter()
        .position(|m| m.0.to_lowercase() == want)?;
    reduce_height(code, i, measured)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_lookups() {
        assert_eq!(antenna(549).map(|a| a.0), Some("TRMR10-2"));
        assert_eq!(igs_name(549), "TRMR10-2        NONE");
        assert_eq!(code_for("R10-2 Internal"), Some(549));
        assert_eq!(code_for("code 549"), Some(549));
        assert_eq!(code_for("trmr10-2"), Some(549));
        assert_eq!(code_for("nothing"), None);
        assert_eq!(method_names(549)[0], "Bottom of antenna mount");
    }

    #[test]
    fn lever_height() {
        // R10-2, lever of R10 extension: sqrt(h^2 - 0.12^2) + 0.3491 - 0.1491
        let (arp, apc) = reduce_height(549, 2, 1.5).unwrap();
        let vert = (1.5f64 * 1.5 - 0.12 * 0.12).sqrt();
        assert_eq!(apc, vert + 0.3491);
        assert_eq!(arp, vert + 0.3491 - 0.1491);
        assert_eq!(
            reduce_height_named(549, "lever of r10 extension", 1.5),
            Some((arp, apc))
        );
    }
}
