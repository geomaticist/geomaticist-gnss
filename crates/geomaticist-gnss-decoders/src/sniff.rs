//! Which receiver log a file is, by its content: the formats decoded here, and the ones only recognised so
//! the user is told plainly what the file is (`Unsupported`), not that it is unreadable.
//!
//! The framing signatures (sync bytes, checksums) are the published ones of each maker; the list of formats
//! follows the receivers RTKLIB reads (BSD-2-Clause), with no code taken from it.

use crate::RawKind;
use crate::topcon;

/// What the start of a file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sniffed {
    /// A log decoded here.
    Supported(RawKind),
    /// A receiver log that is recognised but not decoded (its name, for the message).
    Unsupported(&'static str),
    /// Not a receiver log known here.
    Unknown,
}

fn count(head: &[u8], pat: &[u8]) -> usize {
    head.windows(pat.len()).filter(|w| *w == pat).count()
}

/// SkyTraq binary: `A0 A1`, a big-endian payload length, the payload, an xor check and `0D 0A`.
fn skytraq_frames(head: &[u8]) -> usize {
    let mut n = 0;
    let mut i = 0;
    while i + 7 <= head.len() {
        if head[i] == 0xA0 && head[i + 1] == 0xA1 {
            let len = usize::from(u16::from_be_bytes([head[i + 2], head[i + 3]]));
            let end = i + 4 + len;
            if len > 0 && end + 3 <= head.len() {
                let cs = head[i + 4..end].iter().fold(0u8, |a, b| a ^ b);
                if cs == head[end] && head[end + 1] == 0x0D && head[end + 2] == 0x0A {
                    n += 1;
                    i = end + 3;
                    continue;
                }
            }
        }
        i += 1;
    }
    n
}

/// NovAtel Superstar II: `01`, the message ID, its complement, a length, the payload and a 16-bit sum.
fn superstar_frames(head: &[u8]) -> usize {
    let mut n = 0;
    let mut i = 0;
    while i + 6 <= head.len() {
        if head[i] == 0x01 && head[i + 1] ^ head[i + 2] == 0xFF {
            let len = usize::from(head[i + 3]);
            let end = i + 4 + len;
            if end + 2 <= head.len() {
                let sum = head[i..end]
                    .iter()
                    .fold(0u16, |a, &b| a.wrapping_add(u16::from(b)));
                if sum == u16::from_le_bytes([head[end], head[end + 1]]) {
                    n += 1;
                    i = end + 2;
                    continue;
                }
            }
        }
        i += 1;
    }
    n
}

/// RTCM 2: 30-bit words carried 6 bits per byte with `01` in the top two bits, so a stretch of the stream is
/// almost all bytes 0x40..=0x7F (which plain text, with its spaces, digits and line ends, never is).
fn looks_rtcm2(head: &[u8]) -> bool {
    head.chunks(1024).filter(|c| c.len() == 1024).any(|c| {
        let six = c.iter().filter(|&&b| (0x40..=0x7F).contains(&b)).count();
        six * 100 >= c.len() * 97
    })
}

/// What `head` (the first 256 KiB of a file, or all of a smaller one) is.
#[must_use]
pub fn sniff(head: &[u8]) -> Sniffed {
    if crate::ubx::recognises(head) {
        return Sniffed::Supported(RawKind::Ubx);
    }
    if crate::novatel::recognises(head) {
        return Sniffed::Supported(RawKind::NovAtel);
    }
    if crate::sbf::recognises(head) {
        return Sniffed::Supported(RawKind::Sbf);
    }
    if crate::rtcm3::recognises(head) {
        return Sniffed::Supported(RawKind::Rtcm3);
    }
    if topcon::recognises(head) {
        return Sniffed::Supported(RawKind::Tps);
    }
    if topcon::greis::recognises(head) {
        return Sniffed::Supported(RawKind::Greis);
    }
    if crate::novatel::is_oem3(head) {
        return Sniffed::Unsupported("NovAtel OEM3");
    }
    if count(head, b"$BIN") >= 2 {
        return Sniffed::Unsupported("Hemisphere (Crescent) binary");
    }
    if skytraq_frames(head) >= 2 {
        return Sniffed::Unsupported("SkyTraq binary");
    }
    if superstar_frames(head) >= 3 {
        return Sniffed::Unsupported("NovAtel Superstar II");
    }
    if head.first() == Some(&0x8B) && count(&head[..head.len().min(4096)], &[0x8B]) >= 16 {
        return Sniffed::Unsupported("Furuno GW-10");
    }
    if looks_rtcm2(head) {
        return Sniffed::Unsupported("RTCM 2");
    }
    Sniffed::Unknown
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_the_formats() {
        let mut ubx = crate::ubx::frame(0x01, 0x07, &[0; 8]);
        ubx.extend(crate::ubx::frame(0x01, 0x07, &[0; 8]));
        assert_eq!(sniff(&ubx), Sniffed::Supported(RawKind::Ubx));

        let mut sbf = crate::sbf::block(5914, 0, 1000, 2300, &[0; 10]);
        sbf.extend(crate::sbf::block(5914, 0, 2000, 2300, &[0; 10]));
        assert_eq!(sniff(&sbf), Sniffed::Supported(RawKind::Sbf));

        let mut hemi = b"$BIN\x01\x00\x04\x00abcd".to_vec();
        hemi.extend_from_slice(b"$BIN\x01\x00\x04\x00abcd");
        assert_eq!(
            sniff(&hemi),
            Sniffed::Unsupported("Hemisphere (Crescent) binary")
        );

        let mut sky = Vec::new();
        for _ in 0..2 {
            sky.extend([0xA0, 0xA1, 0x00, 0x02, 0x10, 0x22, 0x32, 0x0D, 0x0A]);
        }
        assert_eq!(sniff(&sky), Sniffed::Unsupported("SkyTraq binary"));

        let mut ss2 = Vec::new();
        for _ in 0..3 {
            let mut m = vec![0x01, 0x17, 0xE8, 0x02, 0x05, 0x06];
            let s = m.iter().fold(0u16, |a, &b| a + u16::from(b));
            m.extend(s.to_le_bytes());
            ss2.extend(m);
        }
        assert_eq!(sniff(&ss2), Sniffed::Unsupported("NovAtel Superstar II"));

        let rtcm2: Vec<u8> = (0..2048u32).map(|i| 0x40 | (i * 7 % 64) as u8).collect();
        assert_eq!(sniff(&rtcm2), Sniffed::Unsupported("RTCM 2"));

        let text = b"PNT,1,1000.000,2000.000,50.000,CP\r\n".repeat(100);
        assert_eq!(sniff(&text), Sniffed::Unknown);
    }
}
