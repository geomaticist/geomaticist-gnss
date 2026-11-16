# Kinematic post-processing (PPK): library entry points

This note is the contract between the engine and a program that calls it (a desktop application, a script,
another library). It lists what goes in, what comes out, and which items are stable. User-facing documentation
(how to run the command line, what the options and flags mean, how the method works) is in [PPK.md](PPK.md).

Everything here is plain Rust data: no file is read or written by the solver itself, nothing is downloaded,
and no global state is kept. The same entry points are behind the `gnss ppk` command and the Python bindings.

## Crates involved

| crate | what a caller uses from it |
|---|---|
| `geomaticist-gnss-core` | `StationObs` (one receiver's observations; `StationObs::concat` joins files), `Sat`, `event::{Event, EventOffset, Occupation}` (the plain event list) |
| `geomaticist-gnss-rinex` | `obs::read_obs_kinds` (RINEX to `StationObs`, with Doppler), `nav::read_nav` |
| `geomaticist-gnss-products` | `Antex`, `Precise` (optional precise orbits) |
| `geomaticist-gnss-engine` | `ppk::*`: the solver, options, results, event interpolation, `.pos` / CSV text |
| `geomaticist-gnss-events` | optional input readers and writers around events: DJI `Timestamp.MRK`, RINEX event records (flag 5), stop-and-go occupations (flags 2 / 3), a check of a file's epoch times, RINEX with events inserted |

The solver core (`geomaticist-gnss-engine::ppk`) never creates events and knows nothing about any camera or
drone maker: it takes a list of times. The readers that produce such a list live in `geomaticist-gnss-events`.
Raw receiver logs are decoded to RINEX by `geomaticist-gnss-decoders` first (as for static baselines).

## Inputs

```rust
use geomaticist_gnss_engine::ppk::{self, BaseStation, PpkInput, PpkOptions};

let mut input = PpkInput {
    base: &base_obs,            // &StationObs  — the static base
    rover: &rover_obs,          // &StationObs  — the moving receiver (several files: StationObs::concat)
    orbits: &orbits,            // &Orbits      — broadcast and / or precise (first source that has a satellite wins)
    fcn: &fcn,                  // &BTreeMap<Sat, i32> — GLONASS frequency channels (NavData::glo_fcn)
    base_station: BaseStation { // where the base is
        marker_xyz,             // [f64; 3] ECEF metres of the mark (the result is in the frame of this coordinate)
        delta_hen: None,        // Option<[f64; 3]>: antenna height / east / north above the mark; None = the RINEX header's
    },
    antex: Some(&mut antex),    // Option<&mut Antex>: receiver antenna calibrations (base, and rover when its type is known)
    tropo: &tropo_inputs,       // &TropoInputs: GPT3 / VMF3 grids; TropoInputs::default() = standard atmosphere
    events: &events,            // &[Event]: times to position (may be empty)
    occupations: &occupations,  // &[Occupation]: stop-and-go segments to average (may be empty)
};
let result = ppk::solve(&mut input, &PpkOptions::default(), &mut |p| { /* progress */ true })?;
```

- **Observations** are `StationObs` values. Read them with `geomaticist_gnss_rinex::obs::read_obs_kinds(path,
  name, t0, t1, systems, "CLSD")` (the `D` keeps Doppler, which the velocity output and the `hermite2` event
  option use; the solver works without it). A session split over several files (a file per flight while the
  receiver keeps tracking) is joined with `StationObs::concat`, so ambiguities carry across the file
  boundaries. The base only needs to cover the rover's span.
- **Base coordinate**: ECEF of the mark. The solver adds the antenna height (`delta_hen`, default: the base
  file's `ANTENNA: DELTA H/E/N`) and the antenna calibration. Rover positions come out in the frame and at
  the epoch of this coordinate.
- **Events**: `geomaticist_gnss_core::event::Event { t, label, offset }`. `t` is GPS seconds since
  1980-01-06 (the engine's time scale, no leap seconds). `offset` is optional: a vector from the antenna
  phase centre to the point wanted (a camera's sensor centre, say) as `EventOffset { north, east, down }` in
  metres, in the local level frame at the event.
- **Options**: `PpkOptions` (serde; every field has a default; `PpkOptions::default().with_json(text)`
  overrides some fields and rejects unknown ones). The ones a user interface will show: `systems`, `elmask`,
  `direction` (`forward`, `backward`, `combined`), `ar` (`off`, `continuous`), `ar_ratio`, `event_interp`
  (`linear2`, `cubic4`, `hermite2`), `event_max_missing`, `base_max_gap`. The full list with defaults is in
  [PPK.md](PPK.md#options).

## Outputs

`ppk::solve` returns `Result<PpkResult, PpkError>` (`PpkError::Cancelled`, `PpkError::Input(String)`).

```rust
pub struct PpkResult {
    pub epochs: Vec<EpochSolution>,        // the track asked for (options.direction), one entry per rover epoch
    pub forward: Vec<EpochSolution>,       // the two passes, kept when options.keep_passes
    pub backward: Vec<EpochSolution>,
    pub events: Vec<EventResult>,          // one per input event, same order
    pub occupations: Vec<OccupationResult>,// one per input occupation
    pub summary: PpkSummary,               // counts per status, fix share, forward / backward agreement, baseline range, signals, slips
    pub warnings: Vec<String>,             // things the user must see
    pub notes: Vec<String>,                // what the solver did (signals left out, antenna handling)
    pub base_arp: [f64; 3],                // the base antenna reference point used
    pub base_marker: [f64; 3],             // the base mark as given
}

pub struct EpochSolution {
    pub t: f64,                 // GPS seconds: the receive time (the epoch tag less the receiver clock offset)
    pub t_tag: f64,             // the epoch tag as in the file
    pub status: FixStatus,      // Fixed, Float, Single, None
    pub xyz: [f64; 3],          // ECEF metres (NaN when status is None)
    pub llh: [f64; 3],          // latitude, longitude (degrees), ellipsoidal height (m)
    pub cov_enu: [[f64; 3]; 3], // realistic covariance, east / north / up, m²
    pub n_sats: u16,            // satellites in the solution
    pub ratio: f32,             // ratio of this epoch's independent ambiguity search (0 when none ran)
    pub age: f32,               // seconds between the rover epoch and the nearest base epoch used
    pub vel_enu: [f32; 3],      // Doppler velocity, m/s (NaN without Doppler)
}

pub struct EventResult {
    pub index: usize,           // position in the input list
    pub label: String,
    pub t: f64,
    pub status: EventStatus,    // Ok, Gap (interpolated across missing epochs), NoSolution
    pub why: String,            // empty, or why there is no position
    pub fix: FixStatus,         // the worst status of the epochs used
    pub xyz: [f64; 3],          // the point wanted: antenna phase centre plus the event's offset
    pub llh: [f64; 3],
    pub antenna_xyz: [f64; 3],  // the antenna at the event
    pub sigma_enu: [f64; 3],
    pub epochs_used: Vec<usize>,// indices into `epochs`
    pub epochs_fix: Vec<FixStatus>,
    pub method: EventInterp,    // what was actually used (falls back to Linear2 where the choice cannot be formed)
    pub gap_s: f64,             // spacing of the two epochs around the event
    pub missing: usize,         // epochs missing between them
}
```

`FixStatus::q()` gives the RTKLIB quality flag: `Fixed` 1, `Float` 2, `Single` 5, `None` 0 (not written).
All result types are serde (`Serialize`, `Deserialize`); NaN serialises as `null` in JSON.

Text outputs are built from the result by pure functions (the caller writes the file):

- `ppk::pos::epochs_pos(&result, &PosHeader) -> String` — RTKLIB-compatible `.pos`
- `ppk::pos::track_pos(&result, &result.forward, &PosHeader)` — the same for any track
- `ppk::pos::events_pos(&result, &PosHeader) -> String` — the same layout for the events
- `ppk::pos::events_csv(&result) -> String` — one line per event with everything in `EventResult`

`PosHeader::new(&options)` fills the header from the options; set `inputs` (file names for the `% inp file`
lines), `time_decimals` and `velocity` as wanted.

## Progress and cancel

The third argument of `solve` is `&mut dyn FnMut(Progress) -> bool`. It is called a few times per second with
`Progress { stage, done, total }` (`stage`: `Prepare`, `Forward`, `Backward`, `Combine`, `Events`). Return
`false` to stop: `solve` then returns `Err(PpkError::Cancelled)`. The callback runs on the calling thread.
With `direction: combined` the backward pass runs on a second thread at the same time as the forward one,
and progress is reported for the forward one.

`solve` is a pure function of its inputs: calling it from several threads at once, on different data, is
fine. It uses two threads at most and no shared pool.

## Epoch by epoch

`solve` is a loop over three pieces:

- `ppk::prepare::prepare_all` turns every rover epoch and the base epochs around it into a `PreparedEpoch`
  (satellite positions at the transmit time, base corrections brought to the rover's epoch, cycle-slip
  flags, a code position, a Doppler velocity). Internally it steps through the epochs in time order and
  keeps only a memory of the previous one; the satellite positions are computed for all epochs at once
  because that is much faster on recorded data.
- `ppk::filter::PpkFilter::step(&PreparedEpoch) -> EpochSolution` advances the estimator by one epoch. The
  filter holds its own state and no reference to the data set; the backward pass is the same filter fed the
  prepared epochs in reverse.
- `ppk::events::interpolate(&track, &events, method, max_missing)` works on any slice of `EpochSolution`.

A real-time mode is not part of this version. What it needs is a preparation step fed one epoch at a time
from a stream; the filter and the event interpolation are already in that shape. `PreparedEpoch`, `Setup`
and `PpkFilter` are public so that such a front end can be added without touching them, but they are not yet
covered by the stability promise below.

## Event sources (`geomaticist-gnss-events`)

```rust
use geomaticist_gnss_events::{mrk, rinex};

let mut m = mrk::read(path)?;                       // DJI Timestamp.MRK
m.resolve_weeks(rover_obs.t[0]);                    // only matters for files without a GPS week column
let events: Vec<Event> = m.events();                // labels = photo numbers, offsets = antenna-to-camera vectors
let warnings = &m.warnings;                         // week changes, unreadable lines, missing offsets

let bytes = geomaticist_gnss_rinex::open::read_plain(obs_path)?;
let ev = rinex::events_from_bytes(&bytes)?;         // RINEX flag-5 records -> ev.events, ev.warnings
let occ = rinex::occupations_from_bytes(&bytes)?;   // flags 2 / 3 (stop and go) -> Vec<Occupation>
let chk = rinex::check_times(&bytes)?;              // epoch-time glitches -> chk.warnings (nothing is changed)
let ins = rinex::insert_events(&bytes, &events)?;   // ins.data: the file with flag-5 records added
```

Each reader returns its warnings beside the data; nothing is printed and nothing is fixed silently.

## Stability

The names and fields in "Inputs", "Outputs", "Progress and cancel" and "Event sources" are what an
application should build against. Additions (new option fields with defaults, new result fields, new enum
variants such as another `EventInterp`) are expected; removals or changes of meaning are not made without a
note here. Match `EventInterp`, `EventStatus`, `FixStatus` and `Stage` with a wildcard arm.

The static-baseline and PPP entry points (`geomaticist_gnss_engine::solve`, `geomaticist_gnss_ppp::solve`) and
the existing readers (`read_obs`, `read_nav`) are unchanged by the kinematic work. Additions outside
`engine::ppk` made for it: `core::event`, `StationObs::concat`, `rinex::obs::read_obs_kinds` /
`read_obs_bytes_kinds`, the crate `geomaticist-gnss-events`.
