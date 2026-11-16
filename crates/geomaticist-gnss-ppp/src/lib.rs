//! `geomaticist-gnss-ppp`: static precise point positioning (float) and PPP-AR on top of the engine.
//!
//! Modules:
//! - [`options`] (`PppOptions`, `PppResult`, `realistic_cov`, `obs_name`), [`solve`] (the batch
//!   ionosphere-free solve: epochs on the clock grid, satellite states and corrections, arcs, the float
//!   estimate with per-epoch clock elimination, outlier rounds and arc splits), [`resolve`] (PPP-AR by
//!   wide lane / narrow lane with the analysis centre's OSB)
//! - the models (full IERS tide, pole tide, ocean loading, wind-up, Shapiro, relativity) are the engine's
//!   `tides` and `attitude` modules; the OSB / ERP / BLQ readers are in `geomaticist-gnss-products`
//! - [`stages`], [`dump`] (the seven `NN_stage.npz / .json`)
//!
//! Options that change results:
//! - the normal equations are solved by Cholesky each iteration and inverted once at the end (same math
//!   as inverting every iteration, last-bit differences);
//! - ocean tide loading: `PppOptions::otl` applies BLQ coefficients handed to [`solve::solve`]; off by
//!   default.
//! - the wide-lane σ rule of PPP-AR: `PppOptions::wl_estimator` `Ar1` (one independent Melbourne-Wübbena
//!   value per `ar_wl_corr_s`) or `Legacy` (`std / sqrt(n / 10)`). PPP-AR itself stays off by default
//!   (`ar: false`).
//!
//! The crate takes parsed products handed to it; it never downloads (rule E1).

// Numeric code: casts between indices, counts and f64 are deliberate and bounded; exact float
// comparisons are intended; the short names follow the formulas (single letters and all).
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::too_many_arguments,
    clippy::struct_excessive_bools,
    clippy::needless_range_loop,
    clippy::doc_markdown,
    clippy::unreadable_literal,
    clippy::inconsistent_digit_grouping,
    clippy::must_use_candidate,
    clippy::needless_pass_by_value,
    clippy::type_complexity,
    clippy::items_after_statements,
    clippy::single_match_else,
    clippy::map_unwrap_or,
    clippy::option_if_let_else,
    clippy::redundant_closure_for_method_calls,
    clippy::if_not_else,
    clippy::manual_let_else,
    clippy::cognitive_complexity,
    clippy::cast_lossless,
    clippy::range_plus_one,
    clippy::explicit_iter_loop,
    clippy::implicit_clone,
    clippy::return_self_not_must_use,
    clippy::uninlined_format_args,
    clippy::manual_midpoint,
    clippy::match_same_arms,
    clippy::needless_late_init,
    // `!(el >= mask)` on purpose: a NaN elevation fails the mask
    clippy::neg_cmp_op_on_partial_ord
)]

pub mod dump;
pub mod kinematic;
pub mod options;
pub mod resolve;
pub mod solve;
pub mod stages;

pub use kinematic::{KinematicInput, KinematicOptions};
pub use options::{PppOptions, PppResult, WlEstimator, realistic_cov};
pub use solve::{PppError, solve};
pub use stages::Stages;

/// Crate name, as Cargo knows it.
pub const CRATE: &str = env!("CARGO_PKG_NAME");
