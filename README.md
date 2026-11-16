# geomaticist-gnss

A GNSS processing engine in Rust: a library, a command line (`gnss`) and Python bindings.

| what | how |
|---|---|
| **Static baselines** | double differences, float solution, wide-lane and narrow-lane ambiguity fixing with wrong-fix checks, a realistic covariance |
| **Static PPP** | ionosphere-free float solution with precise orbits and clocks, optional ambiguity resolution with an analysis centre's biases |
| **Kinematic (PPK)** | a moving receiver against a static base, epoch by epoch, forward / backward / merged, on-the-fly ambiguity fixing, positions at events (camera exposures) |
| **Kinematic PPP** | a moving receiver on its own, with precise orbits and clocks |
| **RINEX** | read versions 2, 3 and 4 (plain, gzip, Hatanaka); describe, trim, decimate and convert between versions |
| **Receiver formats** | Trimble T04 and T01, Leica MDB / DBX, Topcon / Javad / Sokkia GREIS (`.tps`, `.jps`), u-blox UBX, NovAtel OEM4 / OEMV / OEM6, Septentrio SBF and logged RTCM 3 MSM streams to RINEX; Trimble T01 logs the built-in reader doesn't know go through the maker's unpacking program |
| **Products** | SP3, RINEX clock, ANTEX, Bias-SINEX, ERP, BLQ, GPT3 / VMF3 grids |
| **Outputs** | text and JSON; RTKLIB-compatible `.pos` for kinematic tracks and events |

GPS, GLONASS, Galileo, BeiDou and QZSS.

## Why another GNSS engine

There are good open GNSS engines already: RTKLIB and its forks, Ginan, PRIDE PPP-AR, gLAB, GREAT, and the
Rust crates of nav-solutions. This one aims at a particular gap: the surveyor's day, from the receiver's own
log to a baseline or a camera position, in one permissively licensed library.

- **From the receiver's log, not just RINEX.** Built-in readers for Trimble, Leica, Topcon / Javad,
  Septentrio, u-blox and NovAtel logs and for RTCM 3 streams. Open readers for some of these formats are
  hard to find.
- **Survey results, not just positions.** Static baselines with ambiguity fixing, wrong-fix checks and a
  covariance meant to be realistic rather than optimistic; PPK with camera events and stop-and-go
  occupations; output that office software and RTKLIB tools already read.
- **Fast.** A one-hour baseline from the sample data below solves in a few hundredths of a second, and a
  full day of 30 s data (2788 epochs, 200 ambiguities) in about a second on a desktop PC. In our own
  side-by-side runs it finished well ahead of the commercial office software we compared it with. Those
  comparisons are not published here, so try it on your own data.
- **Parallel, and repeatable.** The heavy steps (normal-equation assembly, per-satellite preparation and
  cycle-slip detection, the forward and backward passes of a kinematic solution) run on all cores, arranged
  so that the numbers are the same bit for bit whatever the thread count. Several baselines can share one
  thread pool without slowing each other down.
- **Rust, with Python.** A memory-safe core with no C libraries to install, a command line, and Python
  bindings for scripting.
- **Permissive licence.** MIT or Apache-2.0, so it can go into closed products too; just keep the credit
  (see NOTICE).

It is young software. Check its results against what you trust before you rely on them,
and please report what you find.

The engine takes files and parsed products handed to it. It never downloads anything and keeps no state
between calls.

## Quick look

```
gnss rinex info site.24o                                   # what is in a file
gnss convert base.T04 -o work/                             # a receiver log to RINEX
gnss baseline base.24o rover.24o --base-xyz X,Y,Z --nav brdc.24n --antex igs20.atx
gnss ppp site.24o --sp3 orbit.SP3 --clk clock.CLK --antex igs20.atx
gnss ppk base.24o flight.obs --base-llh LAT,LON,H --nav brdc.24n --antex igs20.atx \
         --mrk flight_Timestamp.MRK -o flight.pos
```

Every command prints text for people and, with `--json`, JSON for programs.

[data/cors](data/cors/README.md) has one hour of real data from two NOAA CORS stations to try the commands on.

### Antenna calibrations

Precise work needs an ANTEX file with the satellite and receiver antenna calibrations (`--antex`). It is not
included here (about 60-70 MB, and updated often); download the current one:

- IGS: `igs20.atx` from <https://files.igs.org/pub/station/general/>
- NGS: `ngs20.atx` from <https://geodesy.noaa.gov/ANTCAL/> (also has receiver antennas the IGS file lacks)

Use the file that matches the frame of your orbits and clocks (IGS20 for current products). The receiver
antenna name and radome in the RINEX header must match an entry in the file, or the receiver antenna is
left uncorrected and the result says so.

## Documentation

- [docs/CLI.md](docs/CLI.md): every command of `gnss`, with its options and outputs
- [docs/PPK.md](docs/PPK.md): kinematic processing and events (inputs, outputs, options, method, limits)
- [docs/PPK_API.md](docs/PPK_API.md): the kinematic library entry points for programs that call the engine
- [docs/ARCH.md](docs/ARCH.md): the crates, how they depend on each other, how to build and test
- API documentation: `cargo doc --no-deps --open`

## Build

Rust stable (1.99 or later). No C libraries.

```
cargo build --release          # the gnss and gnss-decode programs in target/release
cargo test --workspace
cargo xtask ci                 # format check, clippy, tests, dependency rules, licence check
```

The Python bindings (`crates/geomaticist-gnss-py`) build with maturin into a wheel for CPython 3.13 or later;
[examples/python](examples/python/README.md) walks through what they can do, on simulated data.

## Status

Under development; interfaces may still change. The static baseline and PPP solvers were validated against
reference results that are not published here. The tests in this repository run on simulated data generated
by the engine itself and on small synthetic files, so they run anywhere; the decoders also have tests that
run on your own receiver logs when you point an environment variable at them (see docs/ARCH.md).

## Licence

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Acknowledgements

The `.pos` output layout and its quality flags follow RTKLIB (T. Takasu, BSD 2-Clause) so that existing
tools read the files; no RTKLIB code is included. Integer ambiguity search is the LAMBDA method
(P. Teunissen and co-workers). Models follow the IERS Conventions and the IGS standards (ANTEX, SP3,
RINEX, Bias-SINEX); troposphere mapping from the Niell, GPT3 and VMF3 models.
