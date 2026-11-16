# Architecture

## Crates

| crate | contents | depends on |
|---|---|---|
| `geomaticist-gnss-core` | constants, GNSS time and leap seconds, systems, satellites, signals, the observation type (`StationObs`), broadcast navigation records, the event types | serde |
| `geomaticist-gnss-rinex` | RINEX 2 / 3 / 4 observation and navigation readers, a byte-exact observation writer (trim, decimate, convert), a writer for any `StationObs`, Hatanaka and gzip | core |
| `geomaticist-gnss-products` | SP3, RINEX clock, ANTEX, Bias-SINEX, ERP, BLQ, GPT3 / VMF3 grids | core |
| `geomaticist-gnss-engine` | models (geometry, tides, troposphere, antennas), orbits, cycle slips, LAMBDA, the static baseline solver, the kinematic solver (`ppk`) | core, rinex, products |
| `geomaticist-gnss-ppp` | static PPP with optional ambiguity resolution, kinematic PPP | core, rinex, products, engine |
| `geomaticist-gnss-events` | event sources: DJI `Timestamp.MRK`, RINEX event records, stop-and-go occupations, epoch-time checks, RINEX with events inserted | core, rinex |
| `geomaticist-gnss-decoders` | receiver formats to RINEX: Trimble T04, Leica MDB / DBX, Topcon / Javad GREIS, u-blox UBX, NovAtel, RTCM 3 MSM logs, Septentrio SBF, Trimble T01 (built-in reader, the maker's unpacking program as fallback); antenna name table; the `gnss-decode` program | core, rinex |
| `geomaticist-gnss-cli` | the `gnss` command line | core, rinex, products, engine, ppp, events |
| `geomaticist-gnss-py` | Python bindings (PyO3, abi3) | core, rinex, products, engine, events, decoders |
| `xtask` | `cargo xtask ci`, `cargo xtask check-deps` | cargo_metadata |

```
cli ──► ppp ──► engine ──► products ──► core
          └──────► rinex ───────────────┘
events ──► rinex, core
decoders ──► rinex, core
```

## Rules

Checked by `cargo xtask check-deps` on the resolved dependency graph:

- **No heavy or networked dependencies.** No crate pulls in a network stack, a GUI or GPU toolkit, GDAL or
  PROJ, directly or not. The engine works on files and parsed data handed to it and never downloads.
- **The processing path does not depend on the receiver decoders.** RINEX in, results out: core, rinex,
  products, engine, ppp, events and cli never depend on `geomaticist-gnss-decoders`. The `gnss convert`
  command runs the separate `gnss-decode` program. (The Python bindings link both.)
- **PyO3 only in the bindings crate**, and nothing depends on it.
- Every crate has tests and is licensed MIT OR Apache-2.0 (`cargo deny` checks the licences of all
  dependencies; no GPL code).

Lints: clippy `pedantic` as warnings, `unwrap_used` and `unsafe_code` denied, `missing_docs` warned; CI runs
clippy with warnings as errors.

## Data flow

**Static baseline** (`engine::solve`): two `StationObs` and an `Orbits` (broadcast and / or precise) go
through typed stages: observations on common epochs with all corrections, cycle-slip arcs, the float
solution, wide lanes, the narrow-lane search, validation, the realistic covariance. Each stage can be dumped.

**Static PPP** (`ppp::solve`): one `StationObs` with precise orbits and clocks through the same kind of
stages: observations with satellite antennas, biases, tides and wind-up; arcs; the float solution with the
receiver clock eliminated per epoch; optional ambiguity resolution; covariance.

**Kinematic** (`engine::ppk::solve`): a signal plan for the receiver pair, every rover epoch prepared once
(satellite positions, base corrections brought to the rover epoch, cycle-slip tests and repair), then a
filter run forward and backward over the prepared epochs, the two tracks merged, events interpolated. See
[PPK.md](PPK.md) for the method and [PPK_API.md](PPK_API.md) for the entry points.

**Kinematic PPP** (`ppp::kinematic::solve`): the first two stages of the static PPP, then a forward and a
backward filter with one position per epoch; the result type, the events and the outputs are those of the
kinematic solver.

**Events** are plain data (`core::event::Event`: a time, a label, an optional antenna-to-point offset). The
solvers never create them; readers in `geomaticist-gnss-events` do.

## Types worth knowing

- **Time**: RINEX files are handled in integer ticks of 0.1 µs since 1980-01-06 (`core::time`), so a time
  is never re-rounded; the solvers use `f64` GPS seconds.
- **`StationObs`**: one receiver's observations: epochs, and per satellite and RINEX code a series of
  values (NaN where missing) with the loss-of-lock indicator for phases.
- **`Orbits`**: a list of sources (broadcast, precise, or caller-supplied); the first that has a satellite
  wins.
- **Results** are serde types; text and JSON are produced from them by the command line.

## Numerics

Dense linear algebra is `faer` run sequentially, so a result does not depend on the number of threads.
Parallelism (rayon) is over independent pieces (satellites, epochs of a batch assembly, the two kinematic
passes) with results gathered in a fixed order.

## Tests

- Unit tests in every module.
- Kinematic, kinematic PPP, events, RTCM: end-to-end on simulated observations
  (`engine::ppk::synth`: a constellation on circular orbits and a receiver on a known path, generated with
  the engine's own models), including a command-line test that goes through files.
- Decoders: unit tests on small synthetic records, and a local test that decodes RTKLIB's public receiver
  sample logs (BSD-2-Clause, not copied here). Point `GEOM_TEST_DATA_DIR` at a folder holding them in
  `rcv_samples/`; without it the test is skipped.

## Build and CI

```
cargo build --workspace
cargo test --workspace
cargo xtask ci          # fmt --check, clippy -D warnings, tests, check-deps, cargo deny check
```

`cargo deny check advisories` needs network access. Line endings are LF.
