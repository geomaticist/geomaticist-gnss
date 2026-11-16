//! `geomaticist-gnss-engine`: Static double-difference baseline engine (float, wide lane, narrow lane, realistic sigma).
//!
//! Modules:
//! - [`models`] (geometry, Sun / Moon, degree-2 tide, antennas)
//! - [`tides`], [`attitude`]
//! - [`orbits`] (broadcast, and `Orbits` over broadcast and precise sources)
//! - [`tropo`]
//! - [`slips`]
//! - [`dd`]
//! - [`ambiguity`]
//! - [`options`], [`solve`] (the staged solver), [`stages`]
//! - [`dump`]
//!
//! [`ppk`]: kinematic post-processing (a moving receiver against a static base, epoch
//! by epoch, with positions at events). It reuses the models above and has its own documentation
//! (`docs/PPK.md`, `docs/PPK_API.md`).
//!
//! The static solve is split into typed stages (observations, arcs, float, wide lane, narrow-lane search,
//! validation, sigma) so each stage can be tested on its own. The crate takes parsed
//! products and paths handed to it; it never downloads (rule E1).
//!
//! Compatibility flags (`options::Compat`) keep the original behaviour where the engine deliberately
//! improves on it: the solid tide model and the vᵀPv form. They default to the original behaviour.

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
    // style lints that the formula-shaped closures trip over
    clippy::must_use_candidate,
    clippy::needless_pass_by_value,
    clippy::assigning_clones,
    clippy::type_complexity,
    clippy::items_after_statements,
    clippy::single_match_else,
    clippy::needless_late_init,
    clippy::map_unwrap_or,
    clippy::option_if_let_else,
    clippy::redundant_closure_for_method_calls,
    clippy::if_not_else,
    clippy::manual_let_else,
    clippy::cognitive_complexity,
    clippy::manual_midpoint,
    clippy::cast_lossless,
    clippy::range_plus_one,
    clippy::explicit_iter_loop,
    clippy::implicit_clone,
    clippy::return_self_not_must_use,
    clippy::uninlined_format_args
)]

pub mod ambiguity;
pub mod arr;
pub mod attitude;
pub mod dd;
pub mod dump;
pub mod fixing;
pub mod linalg;
pub mod models;
pub mod nputil;
pub mod options;
pub mod orbits;
pub mod ppk;
pub mod slips;
pub mod solve;
pub mod stages;
pub mod tides;
pub mod tropo;
pub mod vmf3_coeffs;

pub use options::{Compat, Options, SlipOptions, StaticResult, TideModel};
pub use orbits::{Broadcast, ExternalOrbit, OrbitSource, Orbits};
pub use solve::{EngineError, solve};
pub use stages::Stages;
