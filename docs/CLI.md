# The `gnss` command line

```
gnss rinex info | trim | convert     RINEX observation files
gnss convert                         a receiver log or an RTCM 3 log to RINEX
gnss baseline                        one static baseline
gnss ppp [--kinematic]               precise point positioning, static or one position per epoch
gnss ppk                             a moving receiver against a static base, with events
gnss events list | check | insert    events and epoch-time checks
```

Common to all commands:

- Every command prints text; `--json` prints the same facts as JSON for programs. A failure prints
  `error: ...` on standard error and exits with a non-zero status.
- Observation files may be plain, gzip-compressed or Hatanaka-compressed (`.crx`, `.??d`); the kind is
  recognised from the content.
- Times on the command line are GPS time: `2024-09-23 18:00:00`.
- Positions are ECEF metres `X,Y,Z` unless an option says otherwise.
- Nothing is downloaded. Orbits, clocks, antenna files and grids are files you give.

`gnss COMMAND --help` lists every option; this page explains what the commands do.

## `gnss rinex info FILE`

Describes an observation file: version, first and last epoch, interval, completeness, satellites per system,
cycle slips flagged by the receiver, event records, and a set of checks on what an online processing service
would accept (session length, interval, dual-frequency data, antenna type and height in the header).

`--start`, `--end`, `--interval` describe a selection instead of the whole file.

## `gnss rinex trim FILE -o OUT`

Cuts a time window (`--start`, `--end`) and optionally thins the epochs (`--interval SECONDS`, a multiple of
the file's own), keeping the RINEX version and every observation. The values are copied as text, never
re-rounded.

## `gnss rinex convert FILE -o OUT`

Writes another RINEX version: `--rinex-version` 2.11, 3.02 to 3.05, 4.00 to 4.02 (default 3.04). Options:
`--kinds` (C code, L phase, D Doppler, S signal strength; default `CL`), `--systems` (default all), the same
window and interval options as `trim`, `--eol lf|crlf`, `--no-comments`.

Observation codes are mapped between the two-character names of RINEX 2 and the three-character names of
RINEX 3 and 4 by the tracking attribute; what cannot be expressed in the target version is left out and
listed in the output.

## `gnss convert RAW -o OUT`

A receiver's raw log to RINEX. `OUT` is the observation file to write (navigation files go beside it) or an
existing folder.

| input | output |
|---|---|
| Trimble `.T04` | RINEX 3.05 observation and navigation files |
| Trimble `.T01` `.T00` `.T02` `.R00` | RINEX 2.11. Packed `.T01` logs are read directly. A file the built-in reader doesn't know (it knows `.T01` logs of R8 receivers) goes through the maker's unpacking program when that is installed (`--runpkr00 PATH`, or on `PATH`); the program is not part of this software |
| Leica `.m00` (MDB) | RINEX 3.04; a log that moves between points gets "start moving" and "new occupation" records with the point names |
| RTCM 3 log `.rtcm3` `.rtcm` `.rtc` `.rt3` | RINEX 3.04 observation file |
| Septentrio `.sbf` | RINEX 3.04 observation file and RINEX 3.05 navigation file |

The command runs the separate program `gnss-decode`, which must be beside `gnss` or on `PATH`.

**RTCM 3 logs.** Decoded messages: MSM4, MSM5, MSM6 and MSM7 of GPS, GLONASS, Galileo, QZSS and BeiDou
(observations), 1005 / 1006 (station position and antenna height), 1007 / 1008 / 1033 (antenna and receiver
names). Two things a stream does not carry have to come from outside:

- the **week**: `--date YYYY-MM-DD`, any date within three days of the data. Without it the file's
  modification time is used, and a warning says so;
- the **GLONASS frequency channels**, when the stream uses MSM4 or MSM6: `--nav FILE` (a GLONASS navigation
  file of the day). Without them the GLONASS phases are left out, with a warning. MSM5 and MSM7 carry the
  channels.

Navigation messages in the stream are not decoded: use a broadcast navigation file for the day. A cycle
slip is flagged where a signal's lock time goes down.

**Septentrio SBF logs.** Decoded blocks: MeasEpoch (GPS, GLONASS, Galileo, BeiDou and QZSS observations of
the main antenna), GPSNav, GLONav, GALNav, BDSNav and QZSNav (navigation), GPSIon, GPSUtc and BDSIon (navigation
header), ReceiverSetup (marker, antenna, antenna deltas, receiver), PVTCartesian or PVTGeodetic (approximate
position) and ReceiverTime (leap seconds). Every block's CRC is checked; blocks that fail it are left out and
counted in a warning. A cycle slip is flagged where a signal's lock time drops or restarts, and a half-cycle
ambiguity is flagged from the measurement's own flag.

Limits: QZSS L6 and L1C/B, the L-band and the second and third antennas are left out; SBAS and NavIC
observations are kept but their navigation is not decoded; the navigation header has no GLONASS time
corrections; the receiver serial and firmware are not written to the RINEX header. The decoder has been
tested on synthetic data only, not on a real receiver's log. The GPS L1P phase and Doppler are written, and
a phase with a half-cycle flag is kept and marked in the loss-of-lock indicator.

## `gnss baseline BASE ROVER --base-xyz X,Y,Z --nav FILE...`

One static baseline: the vector from the base's antenna reference point (held at `--base-xyz`) to the
rover's, with a covariance meant to be realistic rather than formal, and the record of how the ambiguities
were resolved.

| option | meaning |
|---|---|
| `--nav` (repeatable) | broadcast navigation files, merged |
| `--sp3`, `--clk` | precise orbits and clocks, used before the broadcast ones |
| `--antex` | antenna calibrations; the antenna types come from the file headers |
| `--rover-xyz0` | start position of the rover (default: the header's approximate position, else a code solution) |
| `--start`, `--end`, `--systems` | selection (default systems `GRECJ`) |
| `--options FILE` | a JSON object of solver options; fields left out keep their defaults, unknown fields are an error |
| `--gpt3`, `--orography`, `--vmf3` | troposphere grids for the `gpt3` and `vmf3` models |
| `--dump-dir DIR` | write the intermediate stages of the solver (arrays and JSON) for inspection |

Output: status (`Fixed`, `Float`, `Failed`), the vector in XYZ and in east / north / up, standard
deviations, baseline length, phase RMS, rejected observations, the ratio of the ambiguity search, and notes
(antennas not found, systems left out, why a fix was refused). A receiver that moved during the session is
refused: use `gnss ppk`.

## `gnss ppp OBS --sp3 FILE...`

Precise point positioning of one receiver with an analysis centre's products.

| option | meaning |
|---|---|
| `--sp3`, `--clk` | precise orbits (one file per day of the session) and clocks |
| `--antex` | satellite and receiver antenna calibrations |
| `--osb` | Bias-SINEX observable-specific biases: code biases, and phase biases for ambiguity resolution |
| `--erp` | Earth rotation parameters (pole tide) |
| `--nav` | a broadcast navigation file, for the GLONASS frequency channels only |
| `--blq`, `--blq-station`, `--otl` | ocean tide loading coefficients, and the switch that applies them |
| `--ar` | resolve the ambiguities (needs phase biases in `--osb`) |
| `--xyz0`, `--name`, `--start`, `--end`, `--systems` | start position, station name, selection (default systems `GRE`) |
| `--options FILE` | a JSON object of solver options |
| `--gpt3`, `--orography`, `--vmf3`, `--dump-dir` | as for `baseline` |

Static (the default): one position for the session, at the antenna reference point, in the frame of the
orbits at the epoch of the data, with a covariance fitted to how such solutions actually scatter, the
zenith delays, the per-system offsets and notes.

`--kinematic`: one position per epoch for a moving receiver, with events and `.pos` outputs as for
`gnss ppk` (`--mrk`, `--rinex-events`, `--events`, `--interp`, `--max-missing`, `--direction`, `-o`,
`--keep-passes`, `--kinematic-options`). What to expect of it is in
[PPK.md](PPK.md#without-a-base-kinematic-ppp).

## `gnss ppk BASE ROVER... `

A moving receiver against a static base, epoch by epoch, with positions at events. Documented in full in
[PPK.md](PPK.md): inputs, the base coordinate and its reference frame, antennas, the `.pos` and CSV
outputs, the event options, every solver option, the method and its limits.

## `gnss events`

- `gnss events list FILE`: the events of a DJI `Timestamp.MRK` file (count, time span, the quality of the
  aircraft's own positions, warnings) or the event records (epoch flag 5) of a RINEX observation file.
- `gnss events check FILE`: looks at the epoch times of an observation file for what recorders are known to
  get wrong (epochs dated 1980, runs of epochs with the date some days off, repeated epochs, a header first
  epoch that is not the first epoch) and reports; nothing is changed.
- `gnss events insert OBS --mrk FILE... | --events FILE -o OUT`: a copy of the observation file with an event
  record for each event, so that other software reads the same events.

See [PPK.md](PPK.md#events).

## Python

`geomaticist_gnss_py` (built with maturin from `crates/geomaticist-gnss-py`) offers `read_obs`, `read_nav`,
`read_sp3`, `raw_kind`, `decode_raw`, `rinex_info`, `static_solve`, `ppk_solve`, `read_mrk`, `rinex_events`,
`rinex_occupations`, `check_times`, `insert_events`, `mark_kinematic`, `ecef_to_llh`, `llh_to_ecef` and
`synth_session` (simulated test data). The functions take file paths or handles to data read in Rust and
return JSON text; the solvers run without the interpreter lock. Worked examples for each are in
[examples/python](../examples/python/README.md).
