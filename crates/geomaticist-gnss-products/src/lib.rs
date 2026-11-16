//! `geomaticist-gnss-products`: parsers for SP3, CLK, ANTEX, Bias-SINEX, ERP, BLQ and the
//! GPT3 / VMF3 grids.
//!
//! Parsers only, with no downloads (rule E1). Modules:
//! - precise orbits and clocks (SP3, CLK, Lagrange): [`precise`];
//! - antenna calibrations (ANTEX) and the satellite PCV of the ionosphere-free combination: [`antex`];
//! - Bias-SINEX OSB, ERP and BLQ readers: [`bias`];
//! - the GPT3 / VMF3 grids and GPT3: [`grids`].

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names
)]

pub mod antex;
pub mod bias;
pub mod grids;
pub mod precise;
pub mod util;

pub use antex::{AntennaCal, Antex, receiver_correction, sat_pcv_if};
pub use bias::{Erp, pole_at, read_blq, read_erp, read_osb};
pub use grids::{Gpt3, Gpt3Grid, gpt3};
pub use precise::Precise;
pub use util::ProductError;
