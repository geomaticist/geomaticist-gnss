//! `geomaticist-gnss-core`: GNSS constants, satellite systems and signals, GNSS time and the
//! leap-second table, and the shared observation and navigation data types the engine consumes.
//!
//! Modules:
//! - [`consts`]: physical constants and carrier frequencies ([`signal::freq`] adds GLONASS channels)
//! - [`signal`], [`obs`]: observation codes, signal choice, combinations
//! - [`nav`]: navigation data and record layouts
//! - [`time`]: epoch ticks, GPS seconds, leap seconds
//! - [`system`]: satellite systems, their names and RINEX order
//!
//! Bottom of the engine: depends on `std` and `serde` only.

// Numeric code: casts between integer ticks, counts and f64 are deliberate and bounded
// (ticks < 2^63, counts < 2^53); exact float comparisons are intended.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names
)]

pub mod consts;
pub mod event;
pub mod nav;
pub mod obs;
pub mod signal;
pub mod system;
pub mod text;
pub mod time;

pub use nav::{GloEph, KeplerEph, NavData};
pub use obs::{Combo, Signal, StationObs};
pub use signal::{ObsCode, SignalKey};
pub use system::{GnssSystem, Sat};
