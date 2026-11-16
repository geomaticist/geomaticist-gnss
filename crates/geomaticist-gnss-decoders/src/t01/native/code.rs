//! Layer 1 of a packed log: the bit stream.
//!
//! After a 13-byte head the file is one stream of bits, most significant bit of each byte first, that codes
//! bytes with a fixed prefix code: 256 codewords of 3 to 11 bits, one per byte value, the same in every file.
//! The one 11-bit pattern no byte uses is a flush mark: the writer pads with zero bits to the next byte
//! boundary after it and goes on from there.

use std::sync::OnceLock;

/// Bytes at the start of the file that are not part of the coded stream.
pub(crate) const HEAD_LEN: usize = 13;

/// The longest codeword, in bits.
const MAX_BITS: usize = 11;

/// The flush mark: `10111001111`.
const FLUSH: (u16, u8) = (0x5cf, 11);

/// The symbol the lookup table gives for the flush mark.
const FLUSH_SYM: u16 = 256;

/// `(bits, length)` of the codeword of each byte value, the bits right-aligned.
#[rustfmt::skip]
const CODES: [(u16, u8); 256] = [
    (0x004, 3), (0x006, 6), (0x00b, 5), (0x00d, 5), (0x01d, 6), (0x068, 7), (0x05a, 7), (0x004, 7), // 0x00
    (0x03f, 6), (0x001, 8), (0x028, 8), (0x03a, 6), (0x020, 8), (0x024, 8), (0x07c, 7), (0x030, 8), // 0x08
    (0x003, 8), (0x02e, 8), (0x072, 7), (0x027, 8), (0x036, 8), (0x050, 8), (0x0a4, 8), (0x0b2, 8), // 0x10
    (0x0e2, 8), (0x029, 7), (0x056, 9), (0x019, 9), (0x017, 9), (0x0fb, 8), (0x0df, 8), (0x005, 8), // 0x18
    (0x002, 5), (0x0cf, 8), (0x004, 4), (0x030, 6), (0x03d, 7), (0x023, 8), (0x02d, 8), (0x0cb, 8), // 0x20
    (0x0cd, 8), (0x0c3, 9), (0x03f, 7), (0x142, 9), (0x175, 9), (0x1a4, 9), (0x1a5, 9), (0x17e, 9), // 0x28
    (0x036, 6), (0x1e2, 9), (0x1c6, 9), (0x199, 9), (0x172, 9), (0x1dc, 9), (0x193, 9), (0x1ae, 9), // 0x30
    (0x14f, 9), (0x177, 9), (0x198, 9), (0x01a, 9), (0x1bc, 9), (0x1ba, 9), (0x1d8, 9), (0x0d5, 8), // 0x38
    (0x031, 6), (0x1df, 9), (0x02b, 7), (0x1ef, 9), (0x039, 9), (0x1e3, 9), (0x1d9, 9), (0x166, 9), // 0x40
    (0x079, 7), (0x1ed, 9), (0x070, 10), (0x02a, 10), (0x010, 10), (0x1db, 9), (0x1f5, 9), (0x1da, 9), // 0x48
    (0x05b, 7), (0x1dd, 9), (0x018, 10), (0x011, 10), (0x019, 10), (0x1ee, 9), (0x1ec, 9), (0x1bb, 9), // 0x50
    (0x17f, 9), (0x076, 10), (0x088, 10), (0x02b, 10), (0x075, 10), (0x037, 10), (0x077, 10), (0x00d, 9), // 0x58
    (0x031, 8), (0x0e5, 9), (0x005, 9), (0x15d, 9), (0x0a4, 10), (0x0bf, 10), (0x02a, 8), (0x15a, 9), // 0x60
    (0x04a, 9), (0x0d0, 10), (0x190, 10), (0x144, 10), (0x180, 10), (0x0a7, 10), (0x18b, 10), (0x147, 10), // 0x68
    (0x0a6, 8), (0x0ae, 10), (0x0c8, 10), (0x189, 10), (0x198, 10), (0x1e0, 10), (0x19d, 10), (0x181, 10), // 0x70
    (0x191, 10), (0x074, 10), (0x097, 10), (0x15c, 9), (0x182, 10), (0x1ce, 10), (0x1cf, 10), (0x04c, 9), // 0x78
    (0x038, 7), (0x05e, 7), (0x02c, 8), (0x007, 8), (0x037, 8), (0x055, 8), (0x007, 7), (0x021, 8), // 0x80
    (0x070, 7), (0x07d, 8), (0x0a8, 8), (0x0ac, 8), (0x0b8, 8), (0x0be, 8), (0x0c8, 8), (0x0ca, 8), // 0x88
    (0x033, 8), (0x0e7, 8), (0x000, 9), (0x0d3, 8), (0x045, 9), (0x014, 9), (0x04d, 9), (0x016, 9), // 0x90
    (0x0c2, 9), (0x065, 9), (0x0cb, 9), (0x0a8, 9), (0x15e, 9), (0x14b, 9), (0x15f, 9), (0x0d6, 8), // 0x98
    (0x03c, 9), (0x192, 9), (0x1a8, 9), (0x176, 9), (0x013, 10), (0x1f4, 9), (0x002, 10), (0x1c7, 9), // 0xa0
    (0x071, 10), (0x036, 10), (0x07a, 7), (0x1de, 9), (0x089, 10), (0x096, 10), (0x0a6, 10), (0x1a9, 9), // 0xa8
    (0x0f1, 9), (0x003, 10), (0x0a5, 10), (0x030, 10), (0x0bd, 10), (0x145, 10), (0x0d1, 10), (0x058, 7), // 0xb0
    (0x0be, 10), (0x0ce, 8), (0x1cd, 9), (0x035, 8), (0x153, 10), (0x146, 10), (0x152, 10), (0x153, 9), // 0xb8
    (0x01f, 8), (0x0d2, 10), (0x0c9, 10), (0x199, 10), (0x19e, 10), (0x183, 10), (0x19c, 10), (0x19f, 10), // 0xc0
    (0x0d3, 10), (0x194, 10), (0x188, 10), (0x283, 10), (0x1f3, 10), (0x1c8, 10), (0x29d, 10), (0x1e4, 10), // 0xc8
    (0x1cc, 9), (0x19b, 10), (0x2a4, 10), (0x1cc, 10), (0x1f1, 10), (0x282, 10), (0x29c, 10), (0x1e7, 10), // 0xd0
    (0x07a, 10), (0x287, 10), (0x1f0, 10), (0x1e5, 10), (0x19a, 10), (0x1e1, 10), (0x2b6, 10), (0x055, 7), // 0xd8
    (0x1af, 9), (0x280, 10), (0x2e6, 10), (0x2a5, 10), (0x295, 10), (0x286, 10), (0x2cf, 10), (0x281, 10), // 0xe0
    (0x18a, 10), (0x1f2, 10), (0x5ce, 11), (0x294, 10), (0x1e6, 10), (0x1cd, 10), (0x2b7, 10), (0x195, 10), // 0xe8
    (0x012, 10), (0x0f0, 8), (0x2ce, 10), (0x0af, 10), (0x174, 9), (0x07b, 10), (0x1c9, 10), (0x0bc, 10), // 0xf0
    (0x1bd, 9), (0x0dc, 8), (0x031, 10), (0x0c9, 9), (0x063, 8), (0x051, 7), (0x004, 9), (0x007, 5), // 0xf8
];

/// For every 11-bit window: `length << 9 | symbol` of the codeword at its start.
fn lut() -> &'static [u16; 1 << MAX_BITS] {
    static LUT: OnceLock<[u16; 1 << MAX_BITS]> = OnceLock::new();
    LUT.get_or_init(|| {
        let mut t = [0u16; 1 << MAX_BITS];
        let mut put = |bits: u16, len: u8, sym: u16| {
            let free = MAX_BITS - usize::from(len);
            let first = usize::from(bits) << free;
            for e in &mut t[first..first + (1 << free)] {
                *e = u16::from(len) << 9 | sym;
            }
        };
        for (sym, &(bits, len)) in CODES.iter().enumerate() {
            put(bits, len, sym as u16);
        }
        put(FLUSH.0, FLUSH.1, FLUSH_SYM);
        t
    })
}

/// The bytes a coded stream holds (`coded` is the file without its head). A codeword cut by the end of the
/// data is dropped.
#[must_use]
pub(crate) fn decode(coded: &[u8]) -> Vec<u8> {
    let lut = lut();
    let total = coded.len() * 8;
    let mut out = Vec::with_capacity(coded.len());
    let mut pos = 0usize;
    while pos < total {
        // the next 11 bits, zero-filled past the end
        let (byte, shift) = (pos / 8, pos % 8);
        let mut w = 0u32;
        for k in 0..3 {
            w = w << 8 | u32::from(coded.get(byte + k).copied().unwrap_or(0));
        }
        let window = (w >> (24 - MAX_BITS - shift)) as usize & ((1 << MAX_BITS) - 1);
        let e = lut[window];
        let len = usize::from(e >> 9);
        if len == 0 || pos + len > total {
            break;
        }
        pos += len;
        let sym = e & 0x1ff;
        if sym == FLUSH_SYM {
            pos = pos.div_ceil(8) * 8;
        } else {
            out.push(sym as u8);
        }
    }
    out
}

/// Bytes to a coded stream (the inverse of [`decode`]), with a flush mark after the byte counts in `flush_after`
/// and at the end.
#[cfg(test)]
#[must_use]
pub(crate) fn encode(bytes: &[u8], flush_after: &[usize]) -> Vec<u8> {
    let mut out = Vec::new();
    let (mut acc, mut n) = (0u32, 0u32);
    let put = |out: &mut Vec<u8>, acc: &mut u32, n: &mut u32, bits: u16, len: u8| {
        *acc = *acc << len | u32::from(bits);
        *n += u32::from(len);
        while *n >= 8 {
            out.push((*acc >> (*n - 8)) as u8);
            *n -= 8;
        }
        *acc &= (1 << *n) - 1;
    };
    let flush = |out: &mut Vec<u8>, acc: &mut u32, n: &mut u32| {
        if *n > 0 {
            out.push((*acc << (8 - *n)) as u8);
            *acc = 0;
            *n = 0;
        }
    };
    for (i, &b) in bytes.iter().enumerate() {
        let (bits, len) = CODES[usize::from(b)];
        put(&mut out, &mut acc, &mut n, bits, len);
        if flush_after.contains(&(i + 1)) {
            put(&mut out, &mut acc, &mut n, FLUSH.0, FLUSH.1);
            flush(&mut out, &mut acc, &mut n);
        }
    }
    put(&mut out, &mut acc, &mut n, FLUSH.0, FLUSH.1);
    flush(&mut out, &mut acc, &mut n);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_code_is_complete_and_prefix_free() {
        // Kraft sum of the 256 codewords and the flush mark is exactly 1 ...
        let sum: u32 = CODES
            .iter()
            .chain(std::iter::once(&FLUSH))
            .map(|&(_, len)| 1u32 << (MAX_BITS - usize::from(len)))
            .sum();
        assert_eq!(sum, 1 << MAX_BITS);
        // ... and no window is claimed twice (so no codeword is a prefix of another)
        let mut seen = [false; 1 << MAX_BITS];
        for &(bits, len) in CODES.iter().chain(std::iter::once(&FLUSH)) {
            assert!(u32::from(bits) < 1 << len);
            let free = MAX_BITS - usize::from(len);
            for w in 0..1usize << free {
                let i = (usize::from(bits) << free) + w;
                assert!(!seen[i]);
                seen[i] = true;
            }
        }
        assert!(lut().iter().all(|&e| e >> 9 != 0));
    }

    #[test]
    fn known_codewords() {
        // 0x00 is `100`, a space `00010`, the record start 03 02 the ten bits `0110101011`
        assert_eq!(CODES[0], (0b100, 3));
        assert_eq!(CODES[0x20], (0b00010, 5));
        let (etx, stx) = (CODES[3], CODES[2]);
        assert_eq!(
            (u32::from(etx.0) << stx.1 | u32::from(stx.0), etx.1 + stx.1),
            (0b01_1010_1011, 10)
        );
    }

    #[test]
    fn round_trip_with_flushes() {
        let bytes: Vec<u8> = (0..=255u8)
            .chain((0..=255u8).rev())
            .chain([0; 40])
            .collect();
        for flush in [vec![], vec![1], vec![7, 300, 511], vec![552]] {
            let coded = encode(&bytes, &flush);
            assert_eq!(decode(&coded), bytes, "flush after {flush:?}");
        }
        // a stream cut in the middle of a codeword gives the whole bytes before the cut
        let coded = encode(&bytes, &[]);
        let cut = decode(&coded[..coded.len() / 2]);
        assert!(!cut.is_empty() && bytes.starts_with(&cut));
        assert_eq!(decode(&[]).len(), 0);
    }
}
