//! Layer 1 of a Topcon `.tps` log: the message stream.
//!
//! A `.tps` file is the receiver's output stream written to a file: one message after the other, each
//! `ID ID L L L body`, where the two-character ID is printable ASCII (`0` to `~`) and `LLL` is the body length
//! in upper-case hexadecimal. Carriage returns and line feeds between messages are separators. This is the
//! public GREIS framing; the checks are:
//!
//! - binary messages: the last body byte is the 8-bit GREIS checksum (rotate left by 2, xor) of everything
//!   before it, from the first ID character on;
//! - text messages (`PM`, `MF`, `CL`, ...): the body ends in `@` and two hexadecimal digits, the same checksum of
//!   everything before the two digits;
//! - the compact measurement messages `rD` and `rP`: the last two body bytes are a CRC-16 (polynomial 0x1021,
//!   initial value 0, no reflection; "XMODEM"), little-endian, of everything before them;
//! - the message format `MF` ("JP010008F": ID, versions, byte order) ends in the checksum as two hexadecimal
//!   digits, without `@`;
//! - the file identification `JP` at the start of the file has no check (elsewhere it is not accepted).
//!
//! A binary body may end in `@` and two hexadecimal digits by chance; its own checksum is tried first.
//!
//! A message whose header can't be read or whose check fails is skipped one byte at a time until a message
//! passes again (resynchronisation); the skipped bytes are counted.

/// The 8-bit GREIS checksum of `data`: for each byte, rotate left by two bits and xor the byte; rotate once
/// more at the end.
#[must_use]
pub fn cs8(data: &[u8]) -> u8 {
    let mut r = 0u8;
    for &b in data {
        r = r.rotate_left(2) ^ b;
    }
    r.rotate_left(2)
}

/// CRC-16/XMODEM of `data` (polynomial 0x1021, initial value 0, no reflection, no final xor).
#[must_use]
pub fn crc16(data: &[u8]) -> u16 {
    let mut crc: u16 = 0;
    for &b in data {
        crc ^= u16::from(b) << 8;
        for _ in 0..8 {
            crc = if crc & 0x8000 != 0 {
                (crc << 1) ^ 0x1021
            } else {
                crc << 1
            };
        }
    }
    crc
}

/// How a message is checked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckKind {
    /// Last byte: [`cs8`].
    Binary,
    /// Ends in `@HH`: [`cs8`] in hexadecimal.
    Text,
    /// Last two bytes: [`crc16`], little-endian.
    Crc16,
    /// Ends in two hexadecimal digits, [`cs8`] of everything before them, without `@` (`MF`).
    HexTail,
    /// No check (`JP`).
    None,
}

/// The check a message with this ID carries.
#[must_use]
pub fn check_kind(id: [u8; 2], body: &[u8]) -> CheckKind {
    match &id {
        b"rD" | b"rP" => CheckKind::Crc16,
        b"JP" => CheckKind::None,
        b"MF" => CheckKind::HexTail,
        _ if is_text_tail(body) => CheckKind::Text,
        _ => CheckKind::Binary,
    }
}

fn is_hex(c: u8) -> bool {
    c.is_ascii_digit() || (b'A'..=b'F').contains(&c)
}

fn is_text_tail(body: &[u8]) -> bool {
    let n = body.len();
    n >= 3 && body[n - 3] == b'@' && is_hex(body[n - 2]) && is_hex(body[n - 1])
}

/// One message.
#[derive(Debug, Clone, Copy)]
pub struct Message<'a> {
    /// Byte offset of its first ID character in the file.
    pub offset: usize,
    /// The two ID characters.
    pub id: [u8; 2],
    /// The body without its check bytes (text messages: without `@HH`).
    pub body: &'a [u8],
}

impl Message<'_> {
    /// The ID as text.
    #[must_use]
    pub fn id_str(&self) -> &str {
        std::str::from_utf8(&self.id).unwrap_or("??")
    }
}

fn hex_len(h: &[u8]) -> Option<usize> {
    let mut n = 0usize;
    for &c in h {
        let d = match c {
            b'0'..=b'9' => c - b'0',
            b'A'..=b'F' => c - b'A' + 10,
            _ => return None,
        };
        n = n * 16 + usize::from(d);
    }
    Some(n)
}

/// Reads the message at `pos`: `Some((message, next position))` if its header is well formed and it passes its
/// check.
#[must_use]
pub fn message_at(data: &[u8], pos: usize) -> Option<(Message<'_>, usize)> {
    let h = data.get(pos..pos + 5)?;
    if !(0x30..=0x7e).contains(&h[0]) || !(0x30..=0x7e).contains(&h[1]) {
        return None;
    }
    let len = hex_len(&h[2..5])?;
    let end = pos + 5 + len;
    let raw = data.get(pos + 5..end)?;
    let id = [h[0], h[1]];
    if &id == b"JP" && pos != 0 {
        // the unchecked file identification only counts at the start of the file
        return None;
    }
    let all = &data[pos..end];
    let binary = || -> Option<&[u8]> {
        let (&last, rest) = raw.split_last()?;
        (cs8(&all[..all.len() - 1]) == last).then_some(rest)
    };
    let hex_check = || {
        let n = raw.len();
        let want = hex_len(raw.get(n.checked_sub(2)?..)?)?;
        (usize::from(cs8(&all[..all.len() - 2])) == want).then_some(())
    };
    let body = match check_kind(id, raw) {
        CheckKind::Binary => binary()?,
        // a binary body can end in '@' and two hex digits by chance: its own checksum decides first
        CheckKind::Text => {
            if let Some(b) = binary() {
                b
            } else {
                hex_check()?;
                &raw[..raw.len() - 3]
            }
        }
        CheckKind::HexTail => {
            hex_check()?;
            &raw[..raw.len() - 2]
        }
        CheckKind::Crc16 => {
            let n = raw.len();
            if n < 2 {
                return None;
            }
            let want = u16::from_le_bytes([raw[n - 2], raw[n - 1]]);
            (crc16(&all[..all.len() - 2]) == want).then_some(&raw[..n - 2])?
        }
        CheckKind::None => raw,
    };
    Some((
        Message {
            offset: pos,
            id,
            body,
        },
        end,
    ))
}

/// The messages of a file, in order. Bytes that belong to no message that passes its check are skipped and
/// counted in [`Messages::skipped`].
#[derive(Debug, Clone)]
pub struct Messages<'a> {
    data: &'a [u8],
    pos: usize,
    /// Bytes skipped so far (damage, or what the reader doesn't know), separators not included.
    pub skipped: usize,
    /// Number of places where bytes were skipped.
    pub gaps: usize,
    /// Offset of the first skipped byte, if any.
    pub first_gap: Option<usize>,
}

impl<'a> Messages<'a> {
    /// Messages of `data` from its start.
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            skipped: 0,
            gaps: 0,
            first_gap: None,
        }
    }
}

impl<'a> Iterator for Messages<'a> {
    type Item = Message<'a>;

    fn next(&mut self) -> Option<Message<'a>> {
        let mut in_gap = false;
        while self.pos < self.data.len() {
            let c = self.data[self.pos];
            if c == b'\r' || c == b'\n' {
                self.pos += 1;
                continue;
            }
            if let Some((m, next)) = message_at(self.data, self.pos) {
                self.pos = next;
                return Some(m);
            }
            if !in_gap {
                in_gap = true;
                self.gaps += 1;
                self.first_gap.get_or_insert(self.pos);
            }
            self.skipped += 1;
            self.pos += 1;
        }
        None
    }
}

/// Appends one message with its check (for tests and for writing logs): binary messages get [`cs8`], `rD` /
/// `rP` the CRC-16, text messages `@HH`.
pub fn push_message(out: &mut Vec<u8>, id: [u8; 2], body: &[u8], text: bool) {
    let start = out.len();
    let extra = match &id {
        b"rD" | b"rP" => 2,
        b"JP" => 0,
        b"MF" => 4,
        _ if text => 3,
        _ => 1,
    };
    out.extend_from_slice(&id);
    let tail = if extra == 4 { 2 } else { extra };
    out.extend_from_slice(format!("{:03X}", body.len() + tail).as_bytes());
    out.extend_from_slice(body);
    match extra {
        2 => {
            let c = crc16(&out[start..]);
            out.extend_from_slice(&c.to_le_bytes());
        }
        3 => {
            out.push(b'@');
            let c = cs8(&out[start..]);
            out.extend_from_slice(format!("{c:02X}").as_bytes());
        }
        1 => {
            let c = cs8(&out[start..]);
            out.push(c);
        }
        4 => {
            let c = cs8(&out[start..]);
            out.extend_from_slice(format!("{c:02X}").as_bytes());
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_values() {
        // the standard check value of CRC-16/XMODEM
        assert_eq!(crc16(b"123456789"), 0x31C3);
        assert_eq!(cs8(b""), 0);
        assert_eq!(cs8(&[1]), 4);
        assert_eq!(cs8(&[0x80]), 0x02);
    }

    #[test]
    fn frames_and_resyncs() {
        let mut f = Vec::new();
        push_message(&mut f, *b"JP", b"RLOGF test", false);
        f.extend_from_slice(b"\r\n");
        push_message(&mut f, *b"~~", &[1, 2, 3, 4], false);
        f.push(b'\n');
        push_message(&mut f, *b"PM", b"rcv/id=\"X\",", true);
        push_message(&mut f, *b"rD", &[9; 30], false);
        let ids: Vec<String> = Messages::new(&f).map(|m| m.id_str().to_owned()).collect();
        assert_eq!(ids, ["JP", "~~", "PM", "rD"]);
        // the message format, and a binary body that ends like a text one
        let mut g = Vec::new();
        push_message(&mut g, *b"MF", b"JP01000", false);
        assert_eq!(&g[..], b"MF009JP010008F");
        let mut body = vec![7u8; 10];
        body.extend_from_slice(b"@2D");
        push_message(&mut g, *b"GE", &body, false);
        let m: Vec<Message> = Messages::new(&g).collect();
        assert_eq!(m.len(), 2);
        assert_eq!(m[0].body, b"JP01000");
        assert_eq!(m[1].body, &body[..]);
        let m: Vec<Message> = Messages::new(&f).collect();
        assert_eq!(m[1].body, &[1, 2, 3, 4]);
        assert_eq!(m[2].body, b"rcv/id=\"X\",");
        assert_eq!(m[3].body, &[9; 30]);
        // damage one byte of the ~~ body: that message is lost, the reader finds the next one
        let mut g = f.clone();
        let at = m[1].offset + 6;
        g[at] ^= 0x40;
        let mut it = Messages::new(&g);
        let ids: Vec<String> = it.by_ref().map(|m| m.id_str().to_owned()).collect();
        assert_eq!(ids, ["JP", "PM", "rD"]);
        assert_eq!(it.gaps, 1);
        assert!(it.skipped >= 10);
        // a cut file: the last message is incomplete and dropped
        let cut = &f[..f.len() - 3];
        assert_eq!(Messages::new(cut).count(), 3);
    }
}
