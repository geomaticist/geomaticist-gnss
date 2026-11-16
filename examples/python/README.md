# Python examples

Short, commented scripts that show what the GNSS engine can do from Python. Each one runs as it is on
**simulated data** the engine writes itself (no downloads, no real receiver files). Most also take your own
files as arguments.

## Build or install the bindings

The module is `geomaticist_gnss_py`, built from `crates/geomaticist-gnss-py` with
[maturin](https://www.maturin.rs). It needs CPython 3.13 or later and a Rust toolchain.

```bash
pip install maturin
cd crates/geomaticist-gnss-py
maturin develop --release          # builds and installs into the active virtual environment
```

Or build a wheel and install it anywhere:

```bash
maturin build --release            # wheel in target/wheels/
pip install ../../target/wheels/geomaticist_gnss_py-*.whl
```

Without maturin you can build with cargo and copy the library next to the scripts (or onto `PYTHONPATH`)
with the Python extension name: `cargo build --release -p geomaticist-gnss-py`, then copy
`target/release/geomaticist_gnss_py.dll` to `geomaticist_gnss_py.pyd` (Windows), or
`libgeomaticist_gnss_py.so` to `geomaticist_gnss_py.so` (Linux; `.dylib` on macOS).

Check it: `python -c "import geomaticist_gnss_py as g; print(g.version())"`.

## The scripts

Run them from this folder. They write into `out/` (git-ignored): the simulated data under `out/data/`,
results beside it.

| Script | Shows |
|---|---|
| `rinex_inspect.py [OBS]` | Header, span, interval, completeness and gaps, satellites and SNR per system, event (flag 5) and stop-and-go (flags 2/3) records, the time check, NGS OPUS requirements |
| `precise_orbits.py [SP3 [CLK..] --at X,Y,Z]` | Reading SP3 orbits / RINEX clocks; satellite position and clock at any time; azimuth / elevation from a station; ECEF ⇄ latitude/longitude/height |
| `static_baseline.py [BASE ROVER --base-xyz ..]` | One static baseline: ECEF and local vector, realistic sigmas, fixed/float and ratio, QC (RMS, rejections, satellites) |
| `baseline_network.py` | Three baselines from one session, a loop closure (mm and ppm), and a small weighted least-squares adjustment with one station held fixed (plain Python) |
| `ppk_events_to_points.py [BASE ROVER --base-xyz .. --mrk ..]` | Kinematic PPK of a drone flight with camera events (DJI `Timestamp.MRK` or RINEX flag-5) → a point CSV (name from `--prefix`/`--start` or the label, lat/lon/ellipsoid height, sigmas, fix) and the trajectory `.pos` |
| `stop_and_go.py [BASE ROVER --base-xyz ..]` | Stop-and-go: occupations from the rover file's flag-3/flag-2 records → one averaged point per stop |
| `ppp.py [OBS --sp3 .. --clk .. [--kinematic]]` | PPP: one receiver, no base. Static: one position with a realistic sigma; `--kinematic`: a track and photo positions from camera events |
| `rinex_for_tbc.py [ROVER MRK]` | Getting a drone RINEX ready for office software (e.g. Trimble Business Center): `insert_events` writes the photo events as flag-5 records, `mark_kinematic` sets `MARKER TYPE AIRBORNE` and adds a flag-2 record |
| `decode_raw_log.py LOG [OUT_DIR]` | Raw receiver log → RINEX 3 (Trimble T04/T02/T01/T00, Leica MDB, Topcon/Javad GREIS, RTCM 3, u-blox UBX, NovAtel). No sample log ships with the repo, so without a file it prints its usage |
| `run_all.py` | Smoke test: runs every script above on the simulated data |
| `_common.py` | Small helpers the scripts share (station dict, orbit sources, time formatting) |

## Inputs

- **Observations**: RINEX 2.11, 3.0x or 4 observation files, plain, gzip or Hatanaka-compressed. Read with
  `read_obs(path, name="", t0=None, t1=None, systems="GRECJ", kinds="CLS")` (`kinds`: C code, L phase,
  S SNR, D Doppler; the kinematic filter uses Doppler, so the PPK scripts read the rover with `"CLSD"`).
- **Orbits**: broadcast navigation files (`read_nav`) and/or precise SP3 orbits with optional RINEX clocks
  (`read_sp3`). The solvers take a list of sources, e.g. `[("precise_handle", sp3), ("broadcast_handle",
  nav)]`; see `_common.orbit_sources`.
- **Base coordinates**: the base **mark** in ECEF metres. The antenna height and offsets come from the
  file's `ANTENNA: DELTA H/E/N` or the station dict's `delta_hen`.
- **Antenna models** (optional): an ANTEX file (`antex=`, e.g. the IGS `igs20.atx`). The antenna type in
  the RINEX header, or the station dict's `antenna`, selects the model. The simulated files name no
  antenna, so the engine warns that no phase-centre correction is applied; that warning is expected here.
- **Events**: a DJI `Timestamp.MRK` (`read_mrk`), the RINEX flag-5 records (`rinex_events`), or your own
  list of `(gps_seconds, label, None or (north, east, down))`, the last being the lever arm in metres from
  the antenna to the camera.
- **Occupations**: the RINEX flag-3/flag-2 records (`rinex_occupations`) or your own list of
  `(name, gps_t0, gps_t1)`.

Most functions return JSON text; `json.loads` turns it into dicts and lists. The solvers raise
`RuntimeError` for an input they cannot use and `ValueError` for a malformed argument.

## Simulated data

`synth_session(out_dir, kind)` writes a small data set and returns its paths, the base coordinates and the
true positions, which the scripts compare their results with:

- `"static"`: one hour at 15 s from three stations about 1.4 to 2.9 km apart (BASE, PT01, PT02), plus an
  SP3 file. GPS L1/L2, Galileo E1 and BeiDou B1 from simulated constellations; no real satellites.
- `"flight"`: a two-minute drone flight at 5 Hz with a base, and an MRK with 11 photos.
- `"stop_and_go"`: four minutes at 1 Hz with three one-minute stops marked by flag-3/flag-2 records.

Simulated observations are cleaner than field data (no multipath, no real atmosphere), so the errors
against the truth, a millimetre or so, are far better than real data will give.

## The `.pos` file

`ppk_solve(..., outputs="pos")` returns the trajectory as RTKLIB-style text, one epoch per line:

| Column | Meaning |
|---|---|
| `GPST` | Date and time, `YYYY/MM/DD hh:mm:ss.sss`, **GPS time** |
| `latitude(deg)`, `longitude(deg)` | WGS 84 (the frame of the orbits; ITRF with precise orbits) |
| `height(m)` | Ellipsoidal height of the rover antenna, not the ground (no geoid applied) |
| `Q` | 1 fixed, 2 float, 5 single (code only) |
| `ns` | Satellites used |
| `sdn(m)`, `sde(m)`, `sdu(m)` | Standard deviations north, east, up |
| `sdne(m)`, `sdeu(m)`, `sdun(m)` | Signed square roots of the covariances |
| `age(s)` | Age of the base data at the rover epoch |
| `ratio` | Ambiguity validation ratio (higher is more certain) |

The events CSV (`outputs="events_csv"`) and the event `.pos` (`outputs="events_pos"`) use the same
conventions; `ppk_events_to_points.py` writes its own simpler CSV from the `events` list.

## Times

- Every time the bindings take or return is **GPS seconds since 1980-01-06 00:00:00 GPST**, as a float.
  GPS week and seconds of week: `week = t // 604800`, `tow = t % 604800` (`_common.week_tow`).
- GPS time runs ahead of UTC by the leap seconds: **18 s** since 2017-01-01. RINEX observation files,
  SP3 and the `.pos` output are in GPS time; DJI MRK files give GPS week and seconds of week. Convert
  when comparing with a clock or log in UTC or local time.
- MRK lines may omit the GPS week; `read_mrk(path, near=t)` takes it from a nearby time (the rover's first
  epoch).
- `check_times` flags epochs before 1990, which usually mean a receiver hit by the GPS week rollover.

## What the bindings do not cover yet

The engine also has ANTEX lookups, antenna height reduction (slant or to-bottom-of-mount heights to ARP)
and the RINEX trim/convert writer, but these have no Python function yet. The RINEX writer is reachable
through the `gnss` command line (`gnss rinex trim`, `gnss rinex convert`). There is no network adjustment, geoid model or datum
transformation in the engine; `baseline_network.py` shows the adjustment in plain Python.
