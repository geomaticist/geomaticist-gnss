//! Opening RINEX files: plain text, gzip, Hatanaka (CRINEX). Everything is read into memory.
//!
//! Gzip is unpacked in memory and CRINEX (`.crx` / `.??d`) goes through our own decoder
//! ([`crate::crinex`]).
//!
//! Besides the name rule, a file whose first line says `CRINEX VERS` is decoded too.

use std::io::Read;
use std::path::Path;

use crate::RinexError;

/// Gunzips when the bytes are a gzip stream (multi-member streams included).
pub fn gunzip_if(data: Vec<u8>) -> Result<Vec<u8>, RinexError> {
    if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        let mut out = Vec::with_capacity(data.len() * 4);
        flate2::read::MultiGzDecoder::new(data.as_slice())
            .read_to_end(&mut out)
            .map_err(|e| RinexError::new(format!("gzip: {e}")))?;
        Ok(out)
    } else {
        Ok(data)
    }
}

/// The Hatanaka name rule on the (gzip-stripped) lower-case name: `.crx`, or `.??d` with two
/// digits (`.24d`).
#[must_use]
#[allow(clippy::case_sensitive_file_extension_comparisons)] // the name is lower-cased first
pub fn is_hatanaka_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    let b = n.as_bytes();
    n.ends_with(".crx")
        || (b.len() > 4
            && b[b.len() - 4] == b'.'
            && b[b.len() - 1] == b'd'
            && b[b.len() - 3].is_ascii_digit()
            && b[b.len() - 2].is_ascii_digit())
}

/// True when the text starts with a CRINEX header line.
#[must_use]
pub fn looks_like_crinex(data: &[u8]) -> bool {
    let first = data.split(|&b| b == b'\n').next().unwrap_or(&[]);
    first.len() >= 70 && first[60..].starts_with(b"CRINEX VERS")
}

/// The plain RINEX text of a file: gunzipped when the name ends in `.gz` (or the bytes are gzip),
/// Hatanaka-decoded when the name says so or the text is CRINEX.
#[allow(clippy::case_sensitive_file_extension_comparisons)] // the name is lower-cased first
pub fn read_plain(path: &Path) -> Result<Vec<u8>, RinexError> {
    let raw =
        std::fs::read(path).map_err(|e| RinexError::new(format!("{}: {e}", path.display())))?;
    let name = path
        .file_name()
        .map_or_default(|n| n.to_string_lossy().into_owned());
    let lower = name.to_ascii_lowercase();
    let (data, inner) = if lower.ends_with(".gz") {
        (gunzip_if(raw)?, lower[..lower.len() - 3].to_owned())
    } else {
        (gunzip_if(raw)?, lower)
    };
    if is_hatanaka_name(&inner) || looks_like_crinex(&data) {
        return crate::crinex::decompress(&data);
    }
    Ok(data)
}

/// The text of a navigation / product file, gunzipped when it is a gzip stream
/// (decided by content, not by name).
pub fn read_text(path: &Path) -> Result<Vec<u8>, RinexError> {
    let raw =
        std::fs::read(path).map_err(|e| RinexError::new(format!("{}: {e}", path.display())))?;
    gunzip_if(raw)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names() {
        assert!(is_hatanaka_name("FAIR00USA_R_20242670000_01D_30S_MO.crx"));
        assert!(is_hatanaka_name("fair2670.24d"));
        assert!(!is_hatanaka_name("fair2670.24o"));
        assert!(!is_hatanaka_name("x.d"));
    }

    #[test]
    fn gzip_round_trip() {
        use std::io::Write;
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(b"hello").unwrap();
        let z = e.finish().unwrap();
        assert_eq!(gunzip_if(z).unwrap(), b"hello");
        assert_eq!(gunzip_if(b"plain".to_vec()).unwrap(), b"plain");
    }
}
