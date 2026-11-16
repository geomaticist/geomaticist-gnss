//! `geomaticist-gnss-rinex`: RINEX 2/3/4 observation and navigation reading, a
//! byte-exact writer, observation maps and statistics, gzip and Hatanaka (CRINEX) decompression.
//!
//! Modules:
//! - [`header`]: observation header records and edits
//! - [`open`] (plain text, gzip, Hatanaka) and [`reader`]: indexing an observation file
//! - [`obsmap`]: which codes a version can hold
//! - [`stats`]: summaries and checks
//! - [`writer`]: writing an observation file
//! - [`obs`]: reading observations for the engine
//! - [`nav`]: reading navigation files (with the B1 record-type dispatch)
//! - [`crinex`]: the Hatanaka compact RINEX format (Hatanaka 2008, CRINEX 1.0 / 3.0), written
//!   from the published format description.

// Numeric code: casts between integer ticks, counts and f64 are deliberate and bounded; exact
// float comparisons are intended.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]

pub mod crinex;
pub mod header;
pub mod nav;
pub mod obs;
pub mod obsmap;
pub mod open;
pub mod reader;
pub mod station;
pub mod stats;
pub mod writer;

pub use header::{HeaderEdits, HeaderLine, RinexHeader};
pub use reader::{Epoch, RinexObs};

/// An error reading or writing RINEX.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RinexError(pub String);

impl RinexError {
    /// An error with a message.
    pub fn new(msg: impl Into<String>) -> Self {
        Self(msg.into())
    }
}

impl std::fmt::Display for RinexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RinexError {}

impl From<std::io::Error> for RinexError {
    fn from(e: std::io::Error) -> Self {
        Self(e.to_string())
    }
}
