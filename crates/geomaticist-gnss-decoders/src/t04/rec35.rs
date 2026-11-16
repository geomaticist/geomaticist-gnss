//! rec35 (inner record type 0x23): Avro-encoded messages with schemas embedded in the file.
//!
//! Record layout:
//! `[len_lo][0x23][flags<<4 | len_hi][subtype: LEB128 varint][kind: u8][body ...]`;
//! flags: 0 = complete, 1 = first fragment, 2 = middle, 3 = last;
//! kind: 0 = Avro schema (JSON text), 1 = Avro datum (full), 2 = datum of a second, more frequent form
//! (the delta form). Fragments of one message share the subtype; bodies are concatenated in order.

use std::collections::HashMap;

use super::avro::{self, Schema};
use super::container::{ContainerError, GAP, Inner, REC35, Records};

/// Kind: Avro schema.
pub const KIND_SCHEMA: i32 = 0;
/// Kind: full datum (keyframe).
pub const KIND_FULL: i32 = 1;
/// Kind: delta datum.
pub const KIND_DELTA: i32 = 2;

/// One message.
#[derive(Debug, Clone)]
pub struct Message {
    /// Inner-stream offset of the first fragment.
    pub offset: usize,
    /// Subtype ([`GAP`] for a cut-out spot).
    pub subtype: i64,
    /// Kind ([`GAP`] for a cut-out spot).
    pub kind: i32,
    /// Joined body.
    pub body: Vec<u8>,
    /// Number of fragments.
    pub nfrag: usize,
}

/// Errors while reading messages.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rec35Error(pub String);

impl std::fmt::Display for Rec35Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Rec35Error {}

impl From<ContainerError> for Rec35Error {
    fn from(e: ContainerError) -> Self {
        Self(e.0)
    }
}

impl From<avro::AvroError> for Rec35Error {
    fn from(e: avro::AvroError) -> Self {
        Self(e.0)
    }
}

struct Pending {
    offset: usize,
    kind: i32,
    body: Vec<u8>,
    nfrag: usize,
}

/// Messages in order. Where damaged data was cut out a message with subtype [`GAP`]
/// is yielded and unfinished fragments are dropped. Only the bodies of subtypes `keep` accepts are collected;
/// the others are tracked (so a continuation without a start is still an error) but yielded with an empty body.
pub struct Messages<'a, F: Fn(i64) -> bool> {
    recs: Records<'a>,
    pending: HashMap<i64, Pending>,
    after_gap: bool,
    keep: F,
    done: bool,
}

impl<'a, F: Fn(i64) -> bool> Messages<'a, F> {
    /// Iterator over the messages of `inner`, keeping the bodies of the subtypes `keep` accepts.
    pub fn new(inner: &'a Inner, keep: F) -> Self {
        Self {
            recs: Records::new(inner),
            pending: HashMap::new(),
            after_gap: false,
            keep,
            done: false,
        }
    }
}

impl<F: Fn(i64) -> bool> Iterator for Messages<'_, F> {
    type Item = Result<Message, Rec35Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            let r = match self.recs.next()? {
                Ok(r) => r,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e.into()));
                }
            };
            if r.rtype == GAP {
                self.pending.clear();
                self.after_gap = true;
                return Some(Ok(Message {
                    offset: r.offset,
                    subtype: i64::from(GAP),
                    kind: GAP,
                    body: Vec::new(),
                    nfrag: 0,
                }));
            }
            if r.rtype != i32::from(REC35) {
                continue;
            }
            let raw = r.raw;
            let flags = raw.get(2).copied().unwrap_or(0) >> 4;
            let parsed = avro::read_uvarint(raw, 3).and_then(|(sub, q)| {
                raw.get(q)
                    .map(|&_k| (sub, q))
                    .ok_or_else(|| avro::AvroError("index out of range".into()))
            });
            let (sub, q) = match parsed {
                Ok(v) => v,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e.into()));
                }
            };
            #[allow(clippy::cast_possible_wrap)]
            let sub = sub as i64;
            let kind = i32::from(raw[q]);
            let body = &raw[q + 1..];
            let keep = (self.keep)(sub);
            if flags == 0 || flags == 1 {
                self.pending.insert(
                    sub,
                    Pending {
                        offset: r.offset,
                        kind,
                        body: if keep { body.to_vec() } else { Vec::new() },
                        nfrag: 1,
                    },
                );
            } else if let Some(p) = self.pending.get_mut(&sub) {
                if keep {
                    p.body.extend_from_slice(body);
                }
                p.nfrag += 1;
            } else {
                if self.after_gap {
                    continue; // tail of a message whose start was in the damaged part
                }
                self.done = true;
                return Some(Err(Rec35Error(format!(
                    "continuation without start, sub {sub} at {:#x}",
                    r.offset
                ))));
            }
            if (flags == 0 || flags == 3)
                && let Some(p) = self.pending.remove(&sub)
            {
                return Some(Ok(Message {
                    offset: p.offset,
                    subtype: sub,
                    kind: p.kind,
                    body: p.body,
                    nfrag: p.nfrag,
                }));
            }
        }
    }
}

/// The schema: the JSON text with trailing `\0`, space and newline padding removed.
pub fn parse_schema(body: &[u8]) -> Result<serde_json::Value, Rec35Error> {
    let mut end = body.len();
    while end > 0 && matches!(body[end - 1], 0 | b' ' | b'\n') {
        end -= 1;
    }
    let text = std::str::from_utf8(&body[..end]).map_err(|e| Rec35Error(e.to_string()))?;
    serde_json::from_str(text).map_err(|e| Rec35Error(e.to_string()))
}

/// The last schema of each subtype in the file, compiled.
pub fn schemas(inner: &Inner) -> Result<HashMap<i64, Schema>, Rec35Error> {
    let mut json = HashMap::new();
    for m in Messages::new(inner, |_| true) {
        let m = m?;
        if m.kind == KIND_SCHEMA && m.body.len() > 2 {
            json.insert(m.subtype, parse_schema(&m.body)?);
        }
    }
    let mut out = HashMap::new();
    for (k, j) in json {
        // a schema this decoder can't compile is only an error when its messages are decoded
        if let Ok(s) = avro::compile(&j) {
            out.insert(k, s);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(flags: u8, sub: u8, kind: u8, body: &[u8]) -> Vec<u8> {
        let mut r = vec![0, REC35, flags << 4, sub, kind];
        r.extend_from_slice(body);
        r[0] = u8::try_from(r.len()).unwrap();
        r
    }

    #[test]
    fn fragments_join() {
        let mut data = rec(1, 19, 1, b"ab");
        data.extend(rec(2, 19, 1, b"cd"));
        data.extend(rec(0, 5, 2, b"x"));
        data.extend(rec(3, 19, 1, b"e"));
        let inner = Inner {
            data,
            ..Inner::default()
        };
        let msgs: Vec<Message> = Messages::new(&inner, |_| true)
            .map(Result::unwrap)
            .collect();
        assert_eq!(msgs.len(), 2);
        assert_eq!(msgs[0].subtype, 5);
        assert_eq!(msgs[1].body, b"abcde");
        assert_eq!(msgs[1].nfrag, 3);
    }

    #[test]
    fn continuation_without_start() {
        let inner = Inner {
            data: rec(2, 19, 1, b"cd"),
            ..Inner::default()
        };
        assert!(Messages::new(&inner, |_| true).next().unwrap().is_err());
    }

    #[test]
    fn schema_padding() {
        let v = parse_schema(b"{\"type\": \"long\"}\n  \x00").unwrap();
        assert_eq!(v["type"], "long");
    }
}
