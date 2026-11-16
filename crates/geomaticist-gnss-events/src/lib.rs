//! `geomaticist-gnss-events`: where event times come from, and going back the other way.
//!
//! The kinematic solver (`geomaticist-gnss-engine::ppk`) positions a plain list of
//! [`Event`](geomaticist_gnss_core::event::Event)s and never creates one. This crate holds the optional
//! readers that produce such a list, so that the engine stays a GNSS engine with no knowledge of any
//! camera or drone maker:
//!
//! - [`mrk`]: DJI `Timestamp.MRK` files (one line per photo: time, antenna-to-camera offset, the
//!   aircraft's own position);
//! - [`rinex`]: event records in a RINEX observation file (epoch flag 5), stop-and-go occupations (epoch
//!   flags 2 and 3), a check of the file's epoch times for known recorder glitches, and a writer that
//!   inserts event records into an observation file so that other software can read the same events.
//!
//! Every reader returns its warnings beside the data; nothing is printed and nothing is fixed silently.

// times are f64 GPS seconds and tick counts; the casts between them are bounded by the calendar
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

pub mod mrk;
pub mod rinex;

/// An error reading or writing an event file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EventError(pub String);

impl std::fmt::Display for EventError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for EventError {}
