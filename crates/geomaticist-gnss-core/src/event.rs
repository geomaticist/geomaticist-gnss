//! Events: times at which a position is wanted (a camera exposure, a marker button, an external pulse).
//!
//! The processing engine never creates events. It takes a plain list of [`Event`]s and returns a position
//! for each one. Readers for the files that carry such times (RINEX event records, camera timestamp files)
//! live outside the solver, in `geomaticist-gnss-events`.

use serde::{Deserialize, Serialize};

/// A vector from the antenna phase centre to the point that is wanted at an event (a camera's
/// projection centre, the tip of a pole), in the local level frame at the event.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EventOffset {
    /// Metres towards north.
    pub north: f64,
    /// Metres towards east.
    pub east: f64,
    /// Metres **down** (positive when the point is below the antenna).
    pub down: f64,
}

impl EventOffset {
    /// The offset as east, north, up metres.
    #[must_use]
    pub fn enu(&self) -> [f64; 3] {
        [self.east, self.north, -self.down]
    }
}

/// One time to position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Event {
    /// GPS seconds since 1980-01-06 00:00:00 (GPS time has no leap seconds).
    pub t: f64,
    /// A name for the event (a sequence number, a photo name); may be empty.
    pub label: String,
    /// Antenna-to-point offset to apply to the antenna position, when known.
    pub offset: Option<EventOffset>,
}

impl Event {
    /// An event at GPS seconds `t` with no label and no offset.
    #[must_use]
    pub fn at(t: f64) -> Self {
        Self {
            t,
            label: String::new(),
            offset: None,
        }
    }

    /// An event from a GPS week and seconds of week.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn from_week_tow(week: i64, tow: f64) -> Self {
        Self::at(week as f64 * crate::consts::WEEK + tow)
    }
}

/// A stop-and-go occupation: the receiver stood still on a point from `t0` to `t1`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Occupation {
    /// The point's name.
    pub name: String,
    /// Start, GPS seconds.
    pub t0: f64,
    /// End, GPS seconds.
    pub t1: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn week_and_offset() {
        let e = Event::from_week_tow(2300, 345_600.250_000);
        assert!((e.t - (2300.0 * 604_800.0 + 345_600.25)).abs() < 1e-6);
        let o = EventOffset {
            north: -0.441,
            east: 0.348,
            down: 0.303,
        };
        assert_eq!(o.enu(), [0.348, -0.441, -0.303]);
    }
}
