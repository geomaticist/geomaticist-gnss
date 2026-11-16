//! Layer 2 of a packed log: the packets of the decoded byte stream.
//!
//! The stream is a 0x03 followed by packets back to back, each starting with 0x02:
//!
//! - types below 0x20: `02 LEN TYPE payload SUM16 03`, `LEN` the whole packet, `SUM16` the big-endian sum of
//!   every byte before it;
//! - 0x26 (full epoch): `02 LENlo 26 LENhi payload SUM16`, the same sum, no closing 0x03;
//! - 0x25 (epoch as differences): the same frame, but its two last bytes check the epoch after it is rebuilt
//!   (see [`super::epoch`]);
//! - 0x39 (full position) and 0x38 (position as differences): the frame of 0x26; 0x39's check is the sum of its
//!   payload, of the two bytes of its length + 3, and of 0x1D;
//! - 0x29 and any other type from 0x20 up: `02 LEN TYPE ...`, `LEN` the whole packet (contents not read here).

/// Epoch as differences from the epoch before.
pub(crate) const T_EPOCH_DELTA: u8 = 0x25;
/// Full epoch.
pub(crate) const T_EPOCH_FULL: u8 = 0x26;
/// Position as differences from the position before.
pub(crate) const T_POS_DELTA: u8 = 0x38;
/// Full position.
pub(crate) const T_POS_FULL: u8 = 0x39;

/// One packet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Packet<'a> {
    /// Type.
    pub ty: u8,
    /// The bytes between the header and the two check bytes.
    pub payload: &'a [u8],
    /// The two check bytes, big-endian.
    pub check: u16,
    /// Whether the bytes before this packet had to be skipped to find it (damage in the stream).
    pub after_gap: bool,
}

/// The 16-bit sum of bytes.
#[must_use]
pub(crate) fn sum16(b: &[u8]) -> u16 {
    b.iter().fold(0u16, |s, &x| s.wrapping_add(u16::from(x)))
}

/// The check of a record whose own type is `kind` and whose packet is `len` bytes long, over `body`: the sum
/// of the body, of the two bytes of `len + 3`, and of `kind` (0x1B for an epoch, 0x1D for a position).
#[must_use]
pub(crate) fn record_check(body_sum: u16, len: usize, kind: u8) -> u16 {
    let l = len + 3;
    body_sum
        .wrapping_add((l & 0xff) as u16)
        .wrapping_add(((l >> 8) & 0xff) as u16)
        .wrapping_add(u16::from(kind))
}

fn long_frame(ty: u8) -> bool {
    matches!(ty, T_EPOCH_DELTA | T_EPOCH_FULL | T_POS_DELTA | T_POS_FULL)
}

/// The packet at `b[i..]`, when one is there whole and passes the checks its type allows. `in_damage`: the bytes
/// before it belonged to no packet, so one that has no check of its own must be followed by another packet (not
/// by the end of the data) to be believed; an epoch of differences is taken either way, since its content is
/// checked when it is rebuilt.
fn packet_at(b: &[u8], i: usize, in_damage: bool) -> Option<(Packet<'_>, usize)> {
    let rest = &b[i..];
    if rest.len() < 6 || rest[0] != 0x02 {
        return None;
    }
    let ty = rest[2];
    let (len, head, tail) = if long_frame(ty) {
        (usize::from(rest[1]) | usize::from(rest[3]) << 8, 4, 2)
    } else if ty >= 0x20 {
        (usize::from(rest[1]), 3, 2)
    } else {
        (usize::from(rest[1]), 3, 3)
    };
    if len < head + tail || len > rest.len() {
        return None;
    }
    let raw = &rest[..len];
    let check = u16::from(raw[len - tail]) << 8 | u16::from(raw[len - tail + 1]);
    let ok = if ty < 0x20 {
        raw[len - 1] == 0x03 && sum16(&raw[..len - 3]) == check
    } else if ty == T_EPOCH_FULL {
        sum16(&raw[..len - 2]) == check
    } else if ty == T_POS_FULL {
        record_check(sum16(&raw[4..len - 2]), len, 0x1D) == check
    } else {
        // no check of its own bytes: the next packet has to start right after it
        match rest.get(len) {
            Some(&x) => x == 0x02,
            None => !in_damage || ty == T_EPOCH_DELTA,
        }
    };
    ok.then_some((
        Packet {
            ty,
            payload: &raw[head..len - tail],
            check,
            after_gap: false,
        },
        len,
    ))
}

/// The packets of a decoded stream, and how many bytes belonged to none (damage; a packet cut by the end of
/// the file is not counted).
#[must_use]
pub(crate) fn packets(b: &[u8]) -> (Vec<Packet<'_>>, usize) {
    let mut out = Vec::new();
    let mut i = usize::from(b.first() == Some(&0x03));
    let (mut skipped, mut pending) = (0usize, 0usize);
    while i < b.len() {
        if let Some((mut p, len)) = packet_at(b, i, pending > 0) {
            p.after_gap = pending > 0;
            skipped += pending;
            pending = 0;
            out.push(p);
            i += len;
        } else {
            pending += 1;
            i += 1;
        }
    }
    (out, skipped)
}

/// A packet of a type below 0x20.
#[cfg(test)]
#[must_use]
pub(crate) fn classic(ty: u8, payload: &[u8]) -> Vec<u8> {
    assert!(ty < 0x20 && payload.len() + 6 <= 255);
    let mut p = vec![0x02, (payload.len() + 6) as u8, ty];
    p.extend_from_slice(payload);
    let s = sum16(&p);
    p.extend_from_slice(&[(s >> 8) as u8, s as u8, 0x03]);
    p
}

/// A packet with the two-byte length (0x25, 0x26, 0x38, 0x39), its check given or, for `None`, the plain sum.
#[cfg(test)]
#[must_use]
pub(crate) fn long(ty: u8, payload: &[u8], check: Option<u16>) -> Vec<u8> {
    assert!(long_frame(ty));
    let len = payload.len() + 6;
    let mut p = vec![0x02, (len & 0xff) as u8, ty, (len >> 8) as u8];
    p.extend_from_slice(payload);
    let s = check.unwrap_or_else(|| sum16(&p));
    p.extend_from_slice(&[(s >> 8) as u8, s as u8]);
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames() {
        let mut b = vec![0x03];
        b.extend(classic(0x10, b"\xc7Key:Value"));
        b.extend(long(T_EPOCH_FULL, &[1, 2, 3], None));
        b.extend(long(T_EPOCH_DELTA, &[9; 300], Some(0xbeef)));
        let pos = [0x0f, 7, 8];
        let check = record_check(sum16(&pos), pos.len() + 6, 0x1D);
        b.extend(long(T_POS_FULL, &pos, Some(check)));
        b.extend([0x02, 0x07, 0x29, 0x01, 0x20, 0xaa, 0xbb]); // a type that is not read
        b.extend(classic(0x0f, &[0; 104]));
        let (p, skipped) = packets(&b);
        assert_eq!(skipped, 0);
        let types: Vec<u8> = p.iter().map(|x| x.ty).collect();
        assert_eq!(types, [0x10, 0x26, 0x25, 0x39, 0x29, 0x0f]);
        assert_eq!(p[0].payload, b"\xc7Key:Value");
        assert_eq!(p[1].payload, [1, 2, 3]);
        assert_eq!((p[2].payload.len(), p[2].check), (300, 0xbeef));
        assert_eq!(p[3].payload, pos);
        assert_eq!(p[4].payload, [0x01, 0x20]);
        assert!(p.iter().all(|x| !x.after_gap));
    }

    #[test]
    fn damage_and_a_cut_end() {
        let good = classic(0x10, b"\xc7A:1");
        let mut b = vec![0x03];
        b.extend(&good);
        let mut bad = classic(0x10, b"\xc7B:2");
        bad[5] ^= 0x40; // the sum no longer fits
        b.extend(&bad);
        b.extend(&good);
        let cut = classic(0x15, &[5; 185]);
        b.extend(&cut[..100]);
        let (p, skipped) = packets(&b);
        assert_eq!(p.len(), 2);
        assert!(!p[0].after_gap && p[1].after_gap);
        // the damaged packet is skipped; so is the cut one, which nothing follows
        assert_eq!(skipped, bad.len());
        // bytes in a cut tail that look like a packet without a check of its own are not taken for one
        let mut b = vec![0x03];
        b.extend(&good);
        b.extend([0x77, 0x02, 0x06, 0x29, 1, 2, 3]);
        let (p, skipped) = packets(&b);
        assert_eq!((p.len(), skipped), (1, 0));
        // ... but such a packet is the last of a whole file
        b.remove(good.len() + 1);
        assert_eq!(packets(&b).0.len(), 2);
        assert_eq!(packets(&[]).0.len(), 0);
        assert_eq!(packets(&[0x03]).0.len(), 0);
    }
}
