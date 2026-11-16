//! The packets other than epochs, as the unpacked records the writers read (layer 3):
//! navigation data and receiver information are the file's bytes with their multi-byte fields turned from
//! big- to little-endian; positions are scaled integers, in full (0x39) or as differences (0x38).

/// An unpacked record: `'t'`, its length, sub-type, type, body.
#[must_use]
pub(crate) fn rec(ty: u8, sub: u8, body: &[u8]) -> Vec<u8> {
    let mut r = Vec::with_capacity(body.len() + 4);
    r.extend_from_slice(&[0x74, ((body.len() + 4) & 0xff) as u8, sub, ty]);
    r.extend_from_slice(body);
    r
}

/// `b` with each field of the given widths byte-reversed.
fn swapped(b: &[u8], widths: &[usize]) -> Vec<u8> {
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    for &w in widths {
        out.extend(b[i..i + w].iter().rev());
        i += w;
    }
    debug_assert_eq!(i, b.len());
    out
}

fn be_i16(b: &[u8]) -> i16 {
    i16::from_be_bytes([b[0], b[1]])
}

/// Ionosphere / UTC (packet 0x0F, 104 bytes): eleven doubles (alpha 0-3, beta 0-3, A1, A0, tot), two 16-bit
/// leap-second counts and a 32-bit count that the unpacked record holds as doubles, eight bytes as they are.
pub(crate) fn ion_utc(p: &[u8]) -> Option<Vec<u8>> {
    if p.len() != 104 {
        return None;
    }
    let mut o = swapped(&p[..88], &[8; 11]);
    o.extend_from_slice(&f64::from(be_i16(&p[88..90])).to_le_bytes());
    o.extend_from_slice(&f64::from(be_i16(&p[90..92])).to_le_bytes());
    let n = u32::from_be_bytes([p[92], p[93], p[94], p[95]]);
    o.extend_from_slice(&f64::from(n).to_le_bytes());
    o.extend_from_slice(&p[96..104]);
    Some(rec(0x0f, 0, &o))
}

/// GPS ephemeris (packet 0x15, 185 bytes): PRN, week, IODC, two bytes, three 32-bit times, nineteen doubles, a
/// 32-bit field, ten bytes (zero in every sample; copied as they are).
pub(crate) fn gps_ephemeris(p: &[u8]) -> Option<Vec<u8>> {
    if p.len() != 185 {
        return None;
    }
    let mut o = vec![p[0], 0];
    o.extend(swapped(&p[1..5], &[2, 2]));
    o.extend_from_slice(&p[5..7]);
    o.extend(swapped(&p[7..19], &[4, 4, 4]));
    o.extend(swapped(&p[19..171], &[8; 19]));
    o.extend(swapped(&p[171..175], &[4]));
    o.extend_from_slice(&p[175..185]);
    Some(rec(0x15, 0, &o))
}

/// GLONASS ephemeris (packet 0x1C, 135 bytes): system, slot, week, time, day, three bytes, a 32-bit field, six
/// bytes, fourteen doubles.
pub(crate) fn glonass_ephemeris(p: &[u8]) -> Option<Vec<u8>> {
    if p.len() != 135 {
        return None;
    }
    let mut o = p[..2].to_vec();
    o.extend(swapped(&p[2..10], &[2, 4, 2]));
    o.extend_from_slice(&p[10..13]);
    o.extend(swapped(&p[13..17], &[4]));
    o.extend_from_slice(&p[17..23]);
    o.extend(swapped(&p[23..135], &[8; 14]));
    Some(rec(0x1c, 0, &o))
}

/// Receiver settings (packet 0x1F): the first two bytes change places.
pub(crate) fn settings(p: &[u8]) -> Option<Vec<u8>> {
    if p.len() < 2 || p.len() > 251 {
        return None;
    }
    let mut o = vec![p[1], p[0]];
    o.extend_from_slice(&p[2..]);
    Some(rec(0x1f, 0, &o))
}

/// Receiver information, page 1 of packet 0x1E: a count, then fields each led by its length. The first holds a
/// 32-bit and a 16-bit number, the second the serial number, the third the model code, the fourth three 16-bit
/// versions (the first is the firmware, times 100). Gives the unpacked 0x1E/1 record and a 0x05/2 record with
/// what the writers read from it (firmware and serial; the rest of Trimble's 0x05/2 is not rebuilt).
pub(crate) fn receiver(p: &[u8]) -> Option<Vec<Vec<u8>>> {
    if p.first() != Some(&1) || p.len() < 2 {
        return None;
    }
    let mut body = vec![p[1]];
    let (mut i, mut k) = (2, 0);
    let (mut serial, mut versions) = (None, None);
    while i < p.len() {
        let len = usize::from(p[i]);
        let f = p.get(i + 1..i + 1 + len)?;
        i += 1 + len;
        body.push(len as u8);
        match (k, len) {
            (0, 6) => body.extend(swapped(f, &[4, 2])),
            (3, 6) => {
                versions = Some([f[1], f[0], f[3], f[2]]);
                body.extend(swapped(f, &[2, 2, 2]));
            }
            _ => {
                if k == 1 {
                    serial = Some(f);
                }
                body.extend_from_slice(f);
            }
        }
        k += 1;
    }
    if body.len() > 251 {
        return None;
    }
    let mut out = Vec::new();
    if let (Some(serial), Some(v)) = (serial, versions) {
        let mut b = [0u8; 60];
        b[2..6].copy_from_slice(&v);
        let n = serial.len().min(30);
        b[26..26 + n].copy_from_slice(&serial[..n]);
        out.push(rec(0x05, 2, &b));
    }
    out.push(rec(0x1e, 1, &body));
    Some(out)
}

fn int(b: &[u8]) -> i64 {
    let v = b.iter().fold(0i64, |v, &x| v << 8 | i64::from(x));
    let s = 64 - 8 * b.len() as u32;
    (v << s) >> s
}

/// The receiver's position: latitude in 2^-40 degree, longitude in 2^-39 degree, height in 1/4096 m, and the
/// last difference of each.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct Position {
    /// Latitude, longitude, height.
    pub v: [i64; 3],
    /// Their last differences.
    pub d: [i64; 3],
}

impl Position {
    /// A full position (the payload of a packet 0x39): a 15-byte header (week, time, counts), the length of
    /// the block that follows, then latitude and longitude in six bytes each and the height in four.
    pub(crate) fn full(p: &[u8]) -> Result<Self, String> {
        if p.len() < 32 || p[0] != 0x0f || p[15] != 0x35 {
            return Err("a position packet of another layout".to_owned());
        }
        Ok(Position {
            v: [int(&p[16..22]), int(&p[22..28]), int(&p[28..32])],
            d: [0; 3],
        })
    }

    /// The position after this one, from the payload of a packet 0x38: two flag bytes, week and time, up to
    /// two counts and six header bytes named by the second flag byte, a block length, then a nibble per
    /// coordinate (8: the value itself in six bytes; else 4 for a negative second difference and its size in
    /// bytes) and the coordinates. Only the three coordinates are read.
    pub(crate) fn apply_delta(&mut self, p: &[u8]) -> Result<(), String> {
        let short = || "the position packet is too short".to_owned();
        let (f, g) = (*p.first().ok_or_else(short)?, *p.get(1).ok_or_else(short)?);
        if !matches!(f, 0x08 | 0x09 | 0x0a | 0x10) || g & 0x80 == 0 {
            return Err(format!("position flags {f:02x} {g:02x}"));
        }
        let mut i = 8 + usize::from(g >> 4 & 7 == 4) + usize::from(g >> 1 & 7 == 4);
        if g & 1 != 0 {
            i += 6;
        }
        i += 1;
        let sizes = *p.get(i).ok_or_else(short)?;
        i += 1;
        let mut nibble = [sizes >> 4, sizes & 15, 0];
        for k in 0..3 {
            if k == 2 {
                nibble[2] = *p.get(i).ok_or_else(short)? >> 4;
                i += 1;
            }
            let nib = nibble[k];
            if nib & 8 != 0 {
                if k == 2 || nib & 3 != 0 {
                    return Err(format!("a position coordinate described by {nib:x}"));
                }
                let v = int(p.get(i..i + 6).ok_or_else(short)?);
                i += 6;
                self.d[k] = v - self.v[k];
                self.v[k] = v;
            } else {
                let n = usize::from(nib & 3);
                let m = p
                    .get(i..i + n)
                    .ok_or_else(short)?
                    .iter()
                    .fold(0i64, |v, &x| v << 8 | i64::from(x));
                i += n;
                self.d[k] += if nib & 4 != 0 { -m } else { m };
                self.v[k] += self.d[k];
            }
        }
        Ok(())
    }

    /// Latitude and longitude in degrees, height in metres.
    #[must_use]
    pub(crate) fn llh(&self) -> (f64, f64, f64) {
        (
            self.v[0] as f64 / 2f64.powi(40),
            self.v[1] as f64 / 2f64.powi(39),
            self.v[2] as f64 / 4096.0,
        )
    }

    /// The start of an unpacked position record (type 0x0B): latitude and longitude in semicircles and the
    /// height, as doubles (the rest of Trimble's record is not rebuilt).
    #[must_use]
    pub(crate) fn record(&self) -> Vec<u8> {
        let (lat, lon, h) = self.llh();
        let mut b = Vec::with_capacity(24);
        b.extend_from_slice(&(lat / 180.0).to_le_bytes());
        b.extend_from_slice(&(lon / 180.0).to_le_bytes());
        b.extend_from_slice(&h.to_le_bytes());
        rec(0x0b, 0, &b)
    }
}

/// The payload of a full position packet for `lat`, `lon` (degrees) and `h` (m).
#[cfg(test)]
pub(crate) fn encode_position(lat: f64, lon: f64, h: f64) -> Vec<u8> {
    let mut p = vec![
        0x0f, 0x07, 0xd0, 0, 0, 0, 0, 1, 9, 9, 0x60, 6, 5, 0, 0, 0x35,
    ];
    p.extend_from_slice(&((lat * 2f64.powi(40)).round() as i64).to_be_bytes()[2..]);
    p.extend_from_slice(&((lon * 2f64.powi(39)).round() as i64).to_be_bytes()[2..]);
    p.extend_from_slice(&((h * 4096.0).round() as i32).to_be_bytes());
    p.extend_from_slice(&[0; 43]);
    p
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(n: usize) -> Vec<u8> {
        (0..n).map(|i| (i * 7 + 3) as u8).collect()
    }

    #[test]
    fn navigation_records_swap_their_fields() {
        let mut p = seq(104);
        p[..8].copy_from_slice(&1.5e-8f64.to_be_bytes());
        p[88..92].copy_from_slice(&[0, 17, 0xff, 0xfe]);
        p[92..96].copy_from_slice(&1_127_669_296u32.to_be_bytes());
        let r = ion_utc(&p).unwrap();
        assert_eq!((r.len(), &r[..4]), (124, &[0x74, 124, 0, 0x0f][..]));
        let f = |i: usize| f64::from_le_bytes(r[i..i + 8].try_into().unwrap());
        assert_eq!(
            (f(4), f(92), f(100), f(108)),
            (1.5e-8, 17.0, -2.0, 1_127_669_296.0)
        );
        assert_eq!(r[116..], p[96..]);
        assert!(ion_utc(&p[..103]).is_none());

        let mut p = seq(185);
        p[0] = 6;
        p[1..3].copy_from_slice(&2000u16.to_be_bytes());
        p[19..27].copy_from_slice(&(-4.2e-9f64).to_be_bytes());
        let r = gps_ephemeris(&p).unwrap();
        assert_eq!(r.len(), 190);
        assert_eq!((r[4], r[5], u16::from_le_bytes([r[6], r[7]])), (6, 0, 2000));
        assert_eq!(r[10..12], p[5..7]);
        assert_eq!(
            u32::from_le_bytes(r[12..16].try_into().unwrap()),
            u32::from_be_bytes(p[7..11].try_into().unwrap())
        );
        assert_eq!(f64::from_le_bytes(r[24..32].try_into().unwrap()), -4.2e-9);
        assert_eq!(r[180..], p[175..]);
        assert!(gps_ephemeris(&p[1..]).is_none());

        let mut p = seq(135);
        p[23..31].copy_from_slice(&(-1.234e7f64).to_be_bytes());
        let r = glonass_ephemeris(&p).unwrap();
        assert_eq!(r.len(), 139);
        assert_eq!(r[4..6], p[..2]);
        assert_eq!(
            u16::from_le_bytes([r[6], r[7]]),
            u16::from_be_bytes([p[2], p[3]])
        );
        assert_eq!(r[21..27], p[17..23]);
        assert_eq!(f64::from_le_bytes(r[27..35].try_into().unwrap()), -1.234e7);
        assert!(glonass_ephemeris(&p[..10]).is_none());

        assert_eq!(
            settings(&[0x0a, 0x4b, 1, 2]).unwrap(),
            [0x74, 8, 0, 0x1f, 0x4b, 0x0a, 1, 2]
        );
        assert!(settings(&[1]).is_none());
    }

    #[test]
    fn receiver_information() {
        let mut p = vec![1, 4, 6, 0x11, 0x22, 0x33, 0x44, 0x00, 0x2c, 11];
        p.extend_from_slice(b"1234567890\0");
        p.push(5);
        p.extend_from_slice(b"R8-2\0");
        p.extend_from_slice(&[6, 0x01, 0xcd, 0x01, 0xcc, 0x01, 0x9b]);
        let r = receiver(&p).unwrap();
        assert_eq!(r.len(), 2);
        // 0x05/2: firmware at +6, the serial further on
        assert_eq!((r[0][2], r[0][3]), (2, 0x05));
        assert_eq!(u16::from_le_bytes([r[0][6], r[0][7]]), 461);
        assert_eq!(&r[0][30..40], b"1234567890");
        // 0x1E/1: the numbers turned, the texts as they are
        let mut want = vec![
            0x74, 0, 1, 0x1e, 4, 6, 0x44, 0x33, 0x22, 0x11, 0x2c, 0x00, 11,
        ];
        want.extend_from_slice(b"1234567890\0\x05R8-2\0");
        want.extend_from_slice(&[6, 0xcd, 0x01, 0xcc, 0x01, 0x9b, 0x01]);
        want[1] = want.len() as u8;
        assert_eq!(r[1], want);
        // other pages, and a field that runs past the end
        assert!(receiver(&[2, 5, 3, 0, 1]).is_none());
        assert!(receiver(&[1, 1, 9, 0]).is_none());
    }

    #[test]
    fn positions() {
        let mut pos = Position::full(&encode_position(45.5, -110.25, 8.5)).unwrap();
        assert_eq!(pos.llh(), (45.5, -110.25, 8.5));
        let r = pos.record();
        assert_eq!((r.len(), r[3]), (28, 0x0b));
        assert_eq!(
            f64::from_le_bytes(r[4..12].try_into().unwrap()) * 180.0,
            45.5
        );
        // differences: latitude -3 (one byte), longitude +0x0102 (two bytes), height +5 (one byte)
        let head = [0x08, 0x80, 0x07, 0xd0, 1, 2, 3, 4, 0x20];
        let mut p = head.to_vec();
        p.extend_from_slice(&[0x52, 3, 1, 2, 0x10, 5]);
        pos.apply_delta(&p).unwrap();
        assert_eq!(pos.d, [-3, 0x102, 5]);
        // again: second differences add to the first; here none, so the position keeps moving
        let mut p = head.to_vec();
        p.extend_from_slice(&[0x00, 0x00]);
        pos.apply_delta(&p).unwrap();
        let v0 = Position::full(&encode_position(45.5, -110.25, 8.5))
            .unwrap()
            .v;
        assert_eq!(pos.v, [v0[0] - 6, v0[1] + 0x204, v0[2] + 10]);
        // a coordinate given in full, both counts given in full, six more header bytes
        let mut p = vec![
            0x0a, 0xc9, 0x07, 0xd0, 1, 2, 3, 4, 9, 9, 1, 0x60, 6, 5, 0, 0, 0x20, 0x80,
        ];
        p.extend_from_slice(&v0[0].to_be_bytes()[2..]);
        p.extend_from_slice(&[0x40]);
        pos.apply_delta(&p).unwrap();
        assert_eq!((pos.v[0], pos.d[0]), (v0[0], 6));
        assert_eq!(pos.v[1], v0[1] + 0x306);
        // refused: other flags, a height given in full, a cut packet
        assert!(
            pos.clone()
                .apply_delta(&[0x40, 0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                .is_err()
        );
        let mut p = head.to_vec();
        p.extend_from_slice(&[0x00, 0x80, 0, 0, 0, 0]);
        assert!(pos.clone().apply_delta(&p).is_err());
        assert!(pos.clone().apply_delta(&head).is_err());
        assert!(Position::full(&[0; 40]).is_err());
    }
}
