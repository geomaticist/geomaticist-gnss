//! Minimal Apache Avro binary decoder (public spec: avro.apache.org/docs/current/specification).
//!
//! Covers what T04 schemas use: long/int (zigzag varint), arrays, records, plus the common primitives.
//! A schema (JSON, embedded in the file) is compiled once into [`Schema`]; named types are resolved when the
//! schema is compiled.

use std::collections::HashMap;
use std::rc::Rc;

use serde_json::Value as Json;

/// Decoding error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AvroError(pub String);

impl std::fmt::Display for AvroError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for AvroError {}

type R<T> = Result<T, AvroError>;

fn eof() -> AvroError {
    AvroError("index out of range".into())
}

/// Unsigned LEB128 varint at `p` -> (value, next offset).
pub fn read_uvarint(b: &[u8], mut p: usize) -> R<(u64, usize)> {
    let mut v: u64 = 0;
    let mut shift = 0u32;
    loop {
        let byte = *b.get(p).ok_or_else(eof)?;
        p += 1;
        if shift < 64 {
            v |= u64::from(byte & 0x7F) << shift;
        }
        shift += 7;
        if byte < 0x80 {
            return Ok((v, p));
        }
    }
}

/// Avro long/int: zigzag-encoded varint.
#[inline]
pub fn read_long(b: &[u8], p: usize) -> R<(i64, usize)> {
    // fast path: one byte
    if let Some(&byte) = b.get(p)
        && byte < 0x80
    {
        let v = i64::from(byte);
        return Ok(((v >> 1) ^ -(v & 1), p + 1));
    }
    let (v, p) = read_uvarint(b, p)?;
    #[allow(clippy::cast_possible_wrap)]
    let s = ((v >> 1) as i64) ^ -((v & 1) as i64);
    Ok((s, p))
}

/// A compiled schema.
#[derive(Debug, Clone)]
pub enum Schema {
    /// long or int.
    Long,
    /// null.
    Null,
    /// boolean.
    Boolean,
    /// float.
    Float,
    /// double.
    Double,
    /// bytes.
    Bytes,
    /// string.
    Str,
    /// array of items.
    Array(Rc<Schema>),
    /// map of string to values.
    Map(Rc<Schema>),
    /// record: field names and types, in order.
    Record(Rc<str>, Rc<[(Rc<str>, Schema)]>),
    /// enum symbols.
    Enum(Rc<[String]>),
    /// fixed size.
    Fixed(usize),
    /// union of branches.
    Union(Rc<[Schema]>),
}

/// A decoded value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// long / int.
    Long(i64),
    /// null.
    Null,
    /// boolean.
    Bool(bool),
    /// float or double.
    Double(f64),
    /// bytes, string (raw), fixed.
    Bytes(Vec<u8>),
    /// enum symbol.
    Symbol(String),
    /// array.
    Array(Vec<Value>),
    /// map.
    Map(Vec<(String, Value)>),
    /// record fields in schema order.
    Record(Vec<(Rc<str>, Value)>),
}

impl Value {
    /// Field of a record.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Value> {
        match self {
            Value::Record(f) => f.iter().find(|(k, _)| &**k == name).map(|(_, v)| v),
            _ => None,
        }
    }

    /// The long value.
    #[must_use]
    pub fn as_long(&self) -> Option<i64> {
        match self {
            Value::Long(v) => Some(*v),
            _ => None,
        }
    }

    /// The array items.
    #[must_use]
    pub fn as_array(&self) -> Option<&[Value]> {
        match self {
            Value::Array(v) => Some(v),
            _ => None,
        }
    }
}

/// Compiles a schema from its JSON.
pub fn compile(json: &Json) -> R<Schema> {
    let mut named = HashMap::new();
    compile_inner(json, &mut named)
}

fn compile_inner(json: &Json, named: &mut HashMap<String, Schema>) -> R<Schema> {
    match json {
        Json::String(s) => Ok(match s.as_str() {
            "long" | "int" => Schema::Long,
            "null" => Schema::Null,
            "boolean" => Schema::Boolean,
            "float" => Schema::Float,
            "double" => Schema::Double,
            "bytes" => Schema::Bytes,
            "string" => Schema::Str,
            other => named
                .get(other)
                .cloned()
                .ok_or_else(|| AvroError(format!("unsupported type '{other}'")))?,
        }),
        Json::Array(branches) => Ok(Schema::Union(
            branches
                .iter()
                .map(|b| compile_inner(b, named))
                .collect::<R<Vec<_>>>()?
                .into(),
        )),
        Json::Object(o) => {
            let t = o
                .get("type")
                .ok_or_else(|| AvroError("schema without type".into()))?;
            match t.as_str() {
                Some("array") => Ok(Schema::Array(Rc::new(compile_inner(
                    o.get("items").unwrap_or(&Json::Null),
                    named,
                )?))),
                Some("map") => Ok(Schema::Map(Rc::new(compile_inner(
                    o.get("values").unwrap_or(&Json::Null),
                    named,
                )?))),
                Some("record") => {
                    let name = o.get("name").and_then(Json::as_str).unwrap_or("");
                    let mut fields = Vec::new();
                    for f in o
                        .get("fields")
                        .and_then(Json::as_array)
                        .map_or(&[][..], Vec::as_slice)
                    {
                        let fname = f.get("name").and_then(Json::as_str).unwrap_or("");
                        let ft = compile_inner(f.get("type").unwrap_or(&Json::Null), named)?;
                        fields.push((Rc::from(fname), ft));
                    }
                    let s = Schema::Record(Rc::from(name), fields.into());
                    if !name.is_empty() {
                        named.insert(name.to_owned(), s.clone());
                    }
                    Ok(s)
                }
                Some("enum") => Ok(Schema::Enum(
                    o.get("symbols")
                        .and_then(Json::as_array)
                        .map_or(&[][..], Vec::as_slice)
                        .iter()
                        .map(|s| s.as_str().unwrap_or("").to_owned())
                        .collect::<Vec<_>>()
                        .into(),
                )),
                Some("fixed") => Ok(Schema::Fixed(
                    o.get("size")
                        .and_then(Json::as_u64)
                        .and_then(|v| usize::try_from(v).ok())
                        .unwrap_or(0),
                )),
                _ => compile_inner(t, named),
            }
        }
        _ => Err(AvroError(format!("unsupported type {json}"))),
    }
}

fn read_bytes(b: &[u8], p: usize) -> R<(Vec<u8>, usize)> {
    let (n, p) = read_long(b, p)?;
    let n = usize::try_from(n).map_err(|_| eof())?;
    let raw = b.get(p..p + n).ok_or_else(eof)?;
    Ok((raw.to_vec(), p + n))
}

/// Decodes one value at `p` -> (value, next offset).
pub fn read(schema: &Schema, b: &[u8], p: usize) -> R<(Value, usize)> {
    Ok(match schema {
        Schema::Long => {
            let (v, p) = read_long(b, p)?;
            (Value::Long(v), p)
        }
        Schema::Null => (Value::Null, p),
        Schema::Boolean => (Value::Bool(*b.get(p).ok_or_else(eof)? != 0), p + 1),
        Schema::Float => {
            let x = b.get(p..p + 4).ok_or_else(eof)?;
            (
                Value::Double(f64::from(f32::from_le_bytes([x[0], x[1], x[2], x[3]]))),
                p + 4,
            )
        }
        Schema::Double => {
            let x = b.get(p..p + 8).ok_or_else(eof)?;
            let mut a = [0u8; 8];
            a.copy_from_slice(x);
            (Value::Double(f64::from_le_bytes(a)), p + 8)
        }
        Schema::Bytes | Schema::Str => {
            let (v, p) = read_bytes(b, p)?;
            (Value::Bytes(v), p)
        }
        Schema::Array(items) => {
            let mut out = Vec::new();
            let mut p = p;
            loop {
                let (n, q) = read_long(b, p)?;
                p = q;
                if n == 0 {
                    break;
                }
                let n = if n < 0 {
                    let (_, q) = read_long(b, p)?;
                    p = q;
                    -n
                } else {
                    n
                };
                for _ in 0..n {
                    let (v, q) = read(items, b, p)?;
                    p = q;
                    out.push(v);
                }
            }
            (Value::Array(out), p)
        }
        Schema::Map(values) => {
            let mut out = Vec::new();
            let mut p = p;
            loop {
                let (n, q) = read_long(b, p)?;
                p = q;
                if n == 0 {
                    break;
                }
                let n = if n < 0 {
                    let (_, q) = read_long(b, p)?;
                    p = q;
                    -n
                } else {
                    n
                };
                for _ in 0..n {
                    let (k, q) = read_bytes(b, p)?;
                    let (v, q) = read(values, b, q)?;
                    p = q;
                    out.push((String::from_utf8_lossy(&k).into_owned(), v));
                }
            }
            (Value::Map(out), p)
        }
        Schema::Record(_, fields) => {
            let mut out = Vec::with_capacity(fields.len());
            let mut p = p;
            for (name, t) in fields.iter() {
                let (v, q) = read(t, b, p)?;
                p = q;
                out.push((name.clone(), v));
            }
            (Value::Record(out), p)
        }
        Schema::Enum(symbols) => {
            let (i, p) = read_long(b, p)?;
            let s = usize::try_from(i)
                .ok()
                .and_then(|i| symbols.get(i))
                .ok_or_else(eof)?;
            (Value::Symbol(s.clone()), p)
        }
        Schema::Fixed(n) => (
            Value::Bytes(b.get(p..p + n).ok_or_else(eof)?.to_vec()),
            p + n,
        ),
        Schema::Union(branches) => {
            let (i, p) = read_long(b, p)?;
            let t = usize::try_from(i)
                .ok()
                .and_then(|i| branches.get(i))
                .ok_or_else(eof)?;
            read(t, b, p)?
        }
    })
}

/// Decodes a whole message body, which must be consumed exactly.
pub fn decode_full(schema: &Schema, body: &[u8]) -> R<Value> {
    let (v, p) = read(schema, body, 0)?;
    if p != body.len() {
        let name = match schema {
            Schema::Record(n, _) => n.to_string(),
            _ => "None".into(),
        };
        return Err(AvroError(format!(
            "{name}: consumed {p} of {} bytes",
            body.len()
        )));
    }
    Ok(v)
}

/// A field of a flat record: a long or an array of longs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Flat {
    /// long.
    Scalar(i64),
    /// array of long.
    Array(Vec<i64>),
}

/// For a record whose fields are all `long` or arrays of `long` (the Measurement Record): the field names,
/// and whether each is an array. `None` for any other schema.
#[must_use]
pub fn flat_layout(schema: &Schema) -> Option<Vec<(Rc<str>, bool)>> {
    let Schema::Record(_, fields) = schema else {
        return None;
    };
    fields
        .iter()
        .map(|(n, t)| match t {
            Schema::Long => Some((n.clone(), false)),
            Schema::Array(items) if matches!(**items, Schema::Long) => Some((n.clone(), true)),
            _ => None,
        })
        .collect()
}

/// Decodes a flat record body (layout from [`flat_layout`]) into `out`, reusing its buffers; the body must
/// be consumed exactly.
pub fn decode_flat(layout: &[(Rc<str>, bool)], body: &[u8], out: &mut Vec<Flat>) -> R<()> {
    out.resize_with(layout.len(), || Flat::Scalar(0));
    let mut p = 0;
    for (slot, (_, is_array)) in out.iter_mut().zip(layout) {
        if *is_array {
            let mut v = match std::mem::replace(slot, Flat::Scalar(0)) {
                Flat::Array(mut v) => {
                    v.clear();
                    v
                }
                Flat::Scalar(_) => Vec::new(),
            };
            loop {
                let (n, q) = read_long(body, p)?;
                p = q;
                if n == 0 {
                    break;
                }
                let n = if n < 0 {
                    let (_, q) = read_long(body, p)?;
                    p = q;
                    -n
                } else {
                    n
                };
                for _ in 0..n {
                    let (x, q) = read_long(body, p)?;
                    p = q;
                    v.push(x);
                }
            }
            *slot = Flat::Array(v);
        } else {
            let (x, q) = read_long(body, p)?;
            p = q;
            *slot = Flat::Scalar(x);
        }
    }
    if p != body.len() {
        return Err(AvroError(format!(
            "rec35_sub19: consumed {p} of {} bytes",
            body.len()
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zigzag() {
        assert_eq!(read_long(&[0x00], 0).unwrap().0, 0);
        assert_eq!(read_long(&[0x01], 0).unwrap().0, -1);
        assert_eq!(read_long(&[0x02], 0).unwrap().0, 1);
        assert_eq!(read_long(&[0xAC, 0x02], 0).unwrap(), (150, 2));
    }

    #[test]
    fn record_with_array() {
        let schema = compile(
            &serde_json::from_str(
                r#"{"type":"record","name":"r","fields":[{"name":"a","type":"long"},
                {"name":"b","type":{"type":"array","items":"long"}}]}"#,
            )
            .unwrap(),
        )
        .unwrap();
        let body = [0x04, 0x04, 0x02, 0x03, 0x00];
        let v = decode_full(&schema, &body).unwrap();
        assert_eq!(v.get("a"), Some(&Value::Long(2)));
        assert_eq!(
            v.get("b"),
            Some(&Value::Array(vec![Value::Long(1), Value::Long(-2)]))
        );
        let layout = flat_layout(&schema).unwrap();
        let mut out = Vec::new();
        decode_flat(&layout, &body, &mut out).unwrap();
        assert_eq!(out, vec![Flat::Scalar(2), Flat::Array(vec![1, -2])]);
        assert!(decode_full(&schema, &body[..4]).is_err());
    }
}
