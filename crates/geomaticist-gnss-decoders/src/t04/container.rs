//! T04 container layers.
//!
//! Layer 1 (outer): the .T04 file is a sequence of chunks
//! `[type:u8][length:u24 big-endian][payload: length bytes]`:
//! type 0x00 file header (13 bytes, identical in every file so far), 0x05 bzip2 stream ("BZh3..."),
//! 0x08 Zstandard frame (magic 28 B5 2F FD). Decompressed payloads are concatenated into one continuous
//! "inner stream" (records cross chunk boundaries).
//!
//! Layer 2 (inner): the inner stream is a sequence of records `[len_lo:u8][type:u8][...]`. For most types the
//! length is `len_lo` (it counts the whole record, header included). For type 0x23 ("rec35") byte 2 is
//! `(frag_flags << 4) | len_hi` and the length is `len_hi << 8 | len_lo` (12 bits). The chain ends exactly at
//! the end of the stream in every file tested.
//!
//! A damaged chunk (bad block, broken length) is skipped: decoding resumes at the
//! next good chunk and the records cut in half are dropped (the stream keeps `gaps` where data was cut out).

use std::io::Read;

/// Chunk type: file header.
pub const CHUNK_HEADER: u8 = 0x00;
/// Chunk type: bzip2 stream.
pub const CHUNK_BZIP2: u8 = 0x05;
/// Chunk type: Zstandard frame.
pub const CHUNK_ZSTD: u8 = 0x08;
/// Inner record type of the Avro messages.
pub const REC35: u8 = 0x23;

const ZSTD_MAGIC: &[u8] = b"\x28\xb5\x2f\xfd";

/// A damaged spot in the file, for reporting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Damage {
    /// File offset of the bad chunk header.
    pub file_offset: usize,
    /// Bytes skipped up to the next good chunk (or the end of the file).
    pub bytes_skipped: usize,
    /// Why (`"damaged bzip2 block (OSError)"` ...).
    pub reason: String,
}

/// The joined inner stream. `gaps` are stream offsets where damaged data was cut out (the record chain stays
/// valid across them); `damage` describes each damaged spot in the file.
#[derive(Debug, Clone, Default)]
pub struct Inner {
    /// The decompressed records.
    pub data: Vec<u8>,
    /// Offsets in `data` where damaged data was cut out.
    pub gaps: Vec<usize>,
    /// The damaged spots.
    pub damage: Vec<Damage>,
}

/// Errors of the container layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerError(pub String);

impl std::fmt::Display for ContainerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ContainerError {}

fn length24(data: &[u8], p: usize) -> usize {
    (usize::from(data[p + 1]) << 16) | (usize::from(data[p + 2]) << 8) | usize::from(data[p + 3])
}

/// Decompresses one chunk payload; `Err(reason)` for a corrupt block.
fn decompress(ctype: u8, payload: &[u8]) -> Result<Vec<u8>, String> {
    match ctype {
        CHUNK_HEADER => Ok(Vec::new()),
        CHUNK_ZSTD => zstd(payload).map_err(|e| format!("damaged zstd block ({e})")),
        CHUNK_BZIP2 => bzip2(payload).map_err(|e| format!("damaged bzip2 block ({e})")),
        t => Err(format!("unknown chunk type {t:#04x}")),
    }
}

/// A Zstandard frame (a corrupt frame is an error).
fn zstd(payload: &[u8]) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(payload.len() * 4);
    let mut src = payload;
    let mut dec = ruzstd::decoding::StreamingDecoder::new(&mut src).map_err(|_| "ZstdError")?;
    dec.read_to_end(&mut out).map_err(|_| "ZstdError")?;
    Ok(out)
}

/// A bzip2 stream (one or more concatenated streams, as `bz2.decompress`).
fn bzip2(payload: &[u8]) -> Result<Vec<u8>, String> {
    if !payload.starts_with(b"BZh") {
        return Err("OSError".into());
    }
    let mut out = Vec::with_capacity(payload.len() * 4);
    let mut r = bzip2_rs::DecoderReader::new(payload);
    r.read_to_end(&mut out).map_err(|_| "OSError")?;
    Ok(out)
}

fn magic_ok(data: &[u8], q: usize) -> bool {
    let t = data[q];
    (t == CHUNK_ZSTD && data.get(q + 4..q + 8) == Some(ZSTD_MAGIC))
        || (t == CHUNK_BZIP2 && data.get(q + 4..q + 7) == Some(b"BZh".as_slice()))
}

fn find(hay: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if from >= hay.len() {
        return None;
    }
    hay[from..]
        .windows(needle.len())
        .position(|w| w == needle)
        .map(|i| i + from)
}

/// Next offset > p where a chunk header starts and the chain continues from it (next header or EOF).
fn resync(data: &[u8], p: usize) -> Option<usize> {
    let n = data.len();
    let mut q = p + 1;
    while q + 8 <= n {
        let i = [find(data, b"BZh", q + 4), find(data, ZSTD_MAGIC, q + 4)]
            .into_iter()
            .flatten()
            .min()?;
        if i >= 4 {
            let q2 = i - 4;
            if q2 > p && magic_ok(data, q2) {
                let nxt = q2 + 4 + length24(data, q2);
                if nxt == n || (nxt + 8 <= n && magic_ok(data, nxt)) {
                    return Some(q2);
                }
            }
        }
        q = i.saturating_sub(3).max(q + 1);
    }
    None
}

type Segment = (usize, Vec<u8>, Option<Damage>);

/// Decompressed pieces of the file, split where damage was skipped: `(file_offset, bytes, damage_before)`.
fn segments(data: &[u8]) -> Vec<Segment> {
    let mut segs: Vec<Segment> = Vec::new();
    let mut cur: Vec<Vec<u8>> = Vec::new();
    let mut start = 0;
    let mut dmg: Option<Damage> = None;
    let (mut p, n) = (0usize, data.len());
    while p < n {
        let bad: Option<String>;
        let mut out = Vec::new();
        if p + 4 > n {
            bad = Some("truncated chunk header".into());
        } else {
            let ctype = data[p];
            let length = length24(data, p);
            if p + 4 + length > n {
                bad = Some("chunk runs past the end of the file".into());
            } else if ![CHUNK_HEADER, CHUNK_BZIP2, CHUNK_ZSTD].contains(&ctype) {
                bad = Some(format!("unknown chunk type {ctype:#04x}"));
            } else {
                match decompress(ctype, &data[p + 4..p + 4 + length]) {
                    Ok(o) => {
                        out = o;
                        bad = None;
                    }
                    Err(e) => bad = Some(e),
                }
            }
        }
        let Some(reason) = bad else {
            if cur.is_empty() {
                start = p;
            }
            cur.push(out);
            p += 4 + length24(data, p);
            continue;
        };
        let q = resync(data, p);
        if !cur.is_empty() {
            segs.push((start, cur.concat(), dmg.take()));
            cur.clear();
        }
        dmg = Some(Damage {
            file_offset: p,
            bytes_skipped: q.unwrap_or(n) - p,
            reason,
        });
        match q {
            None => {
                segs.push((n, Vec::new(), dmg.take()));
                break;
            }
            Some(q) => p = q,
        }
    }
    if !cur.is_empty() {
        segs.push((start, cur.concat(), dmg.take()));
    }
    segs
}

fn rec_len(seg: &[u8], p: usize, rtype: u8) -> usize {
    let lo = usize::from(seg[p]);
    if rtype == REC35 && p + 2 < seg.len() {
        lo | (usize::from(seg[p + 2] & 0x0F) << 8)
    } else {
        lo
    }
}

/// Offset where the record chain starting at p stops fitting inside seg.
fn chain_end(seg: &[u8], mut p: usize) -> usize {
    let n = seg.len();
    while p + 2 <= n {
        let length = rec_len(seg, p, seg[p + 1]);
        if length < 2 || p + length > n {
            break;
        }
        p += length;
    }
    p
}

/// First offset in a segment that begins after damage from which `need` records parse with known types.
fn record_start(seg: &[u8], types: &[bool; 256], need: usize) -> Option<usize> {
    let n = seg.len();
    for o in 0..n.min(1 << 16) {
        let (mut p, mut k) = (o, 0);
        while k < need && p + 3 <= n {
            let t = seg[p + 1];
            if !types[usize::from(t)] {
                break;
            }
            let length = rec_len(seg, p, t);
            if length < 2 || p + length > n {
                break;
            }
            p += length;
            k += 1;
        }
        if k >= need || (k > 0 && p == n) {
            return Some(o);
        }
    }
    None
}

/// All chunk payloads joined, damage skipped.
#[must_use]
pub fn inner_stream(data: &[u8]) -> Inner {
    let segs = segments(data);
    if segs.iter().all(|s| s.2.is_none()) {
        let mut out = Vec::with_capacity(segs.iter().map(|s| s.1.len()).sum());
        for s in &segs {
            out.extend_from_slice(&s.1);
        }
        return Inner {
            data: out,
            ..Inner::default()
        };
    }
    let mut types = [false; 256];
    let mut ntypes = 0;
    let first: &[u8] = match segs.first() {
        Some((_, s, None)) => s,
        _ => &[],
    };
    let mut p = 0;
    while p + 2 <= first.len() && ntypes < 64 {
        let t = first[p + 1];
        if !types[usize::from(t)] {
            types[usize::from(t)] = true;
            ntypes += 1;
        }
        let length = rec_len(first, p, t);
        if length < 2 {
            break;
        }
        p += length;
    }
    let mut out = Inner::default();
    let mut pos = 0;
    let mut have_parts = false;
    for (_off, seg, dmg) in segs {
        let mut o = 0;
        if let Some(d) = dmg {
            out.damage.push(d);
            let Some(start) = (if seg.is_empty() {
                None
            } else {
                record_start(&seg, &types, 24)
            }) else {
                continue;
            };
            o = start;
            if have_parts {
                out.gaps.push(pos);
            }
        }
        let end = chain_end(&seg, o);
        out.data.extend_from_slice(&seg[o..end]);
        have_parts = true;
        pos += end - o;
    }
    out
}

/// The strict form: every chunk must decompress.
pub fn inner_stream_strict(data: &[u8]) -> Result<Vec<u8>, ContainerError> {
    let mut out = Vec::new();
    let mut p = 0;
    while p < data.len() {
        if p + 4 > data.len() {
            return Err(ContainerError(format!("truncated chunk header at {p:#x}")));
        }
        let length = length24(data, p);
        let Some(payload) = data.get(p + 4..p + 4 + length) else {
            return Err(ContainerError(format!(
                "truncated chunk at {p:#x}: want {length}, have {}",
                data.len() - p - 4
            )));
        };
        out.extend(decompress(data[p], payload).map_err(ContainerError)?);
        p += 4 + length;
    }
    Ok(out)
}

/// The first chunk's payload (the 13-byte file header).
pub fn file_header(data: &[u8]) -> Result<&[u8], ContainerError> {
    if data.len() < 4 || data[0] != CHUNK_HEADER {
        return Err(ContainerError("first chunk is not a header chunk".into()));
    }
    data.get(4..4 + length24(data, 0))
        .ok_or_else(|| ContainerError("truncated chunk at 0x0".into()))
}

/// One inner record.
#[derive(Debug, Clone, Copy)]
pub struct Record<'a> {
    /// Offset in the inner stream.
    pub offset: usize,
    /// Record type ([`GAP`] for a cut-out spot).
    pub rtype: i32,
    /// The whole record, header included (empty for a gap).
    pub raw: &'a [u8],
}

/// Pseudo record type: damaged data was cut out here.
pub const GAP: i32 = -1;

/// The records of an inner stream, gaps included.
pub struct Records<'a> {
    inner: &'a Inner,
    p: usize,
    gap: usize,
    failed: bool,
}

impl<'a> Records<'a> {
    /// Iterator over `inner`.
    #[must_use]
    pub fn new(inner: &'a Inner) -> Self {
        Self {
            inner,
            p: 0,
            gap: 0,
            failed: false,
        }
    }
}

impl<'a> Iterator for Records<'a> {
    type Item = Result<Record<'a>, ContainerError>;

    fn next(&mut self) -> Option<Self::Item> {
        let data = &self.inner.data;
        let n = data.len();
        if self.failed || self.p >= n {
            return None;
        }
        if let Some(&g) = self.inner.gaps.get(self.gap)
            && g <= self.p
        {
            self.gap += 1;
            return Some(Ok(Record {
                offset: g,
                rtype: GAP,
                raw: &[],
            }));
        }
        let p = self.p;
        if p + 2 > n {
            self.failed = true;
            return Some(Err(ContainerError(format!(
                "truncated record header at {p:#x}"
            ))));
        }
        let rtype = data[p + 1];
        let mut length = usize::from(data[p]);
        if rtype == REC35 {
            length |= usize::from(data.get(p + 2).copied().unwrap_or(0) & 0x0F) << 8;
        }
        if length < 2 || p + length > n {
            self.failed = true;
            return Some(Err(ContainerError(format!(
                "bad record length {length} type {rtype:#04x} at {p:#x}"
            ))));
        }
        self.p += length;
        Some(Ok(Record {
            offset: p,
            rtype: i32::from(rtype),
            raw: &data[p..p + length],
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_walk_and_gaps() {
        let inner = Inner {
            data: vec![3, 0x10, 0xC7, 2, 0x15, 4, 0x23, 0x00, 9],
            gaps: vec![3],
            damage: vec![],
        };
        let recs: Vec<_> = Records::new(&inner).map(Result::unwrap).collect();
        let types: Vec<i32> = recs.iter().map(|r| r.rtype).collect();
        assert_eq!(types, vec![0x10, GAP, 0x15, 0x23]);
        assert_eq!(recs[3].raw.len(), 4);
    }

    #[test]
    fn bad_length_is_an_error() {
        let inner = Inner {
            data: vec![9, 0x10, 0],
            ..Inner::default()
        };
        assert!(Records::new(&inner).next().unwrap().is_err());
    }

    #[test]
    fn chain_end_stops_at_partial_record() {
        assert_eq!(chain_end(&[2, 1, 4, 1, 0], 0), 2);
        assert_eq!(chain_end(&[2, 1, 3, 1, 0], 0), 5);
    }
}
