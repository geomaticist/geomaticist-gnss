# Kinematic post-processing (PPK) and events

This is the user documentation of the kinematic part of the engine: a moving receiver (a drone, a vehicle, a
pole) positioned epoch by epoch against a static base receiver, and positions at arbitrary times ("events",
typically camera exposures) interpolated on the resulting track.

- [Quick start](#quick-start)
- [What goes in](#what-goes-in)
- [What comes out](#what-comes-out)
- [Events](#events)
- [Options](#options)
- [How it works](#how-it-works)
- [Without a base: kinematic PPP](#without-a-base-kinematic-ppp)
- [Checks on the input files](#checks-on-the-input-files)
- [Bases from a correction stream](#bases-from-a-correction-stream)
- [Python](#python)
- [Limits](#limits)
- [Credits](#credits)

Programmers who call the library directly: see [PPK_API.md](PPK_API.md).

## Quick start

```
# 1. raw receiver logs to RINEX, when needed (RINEX files are used as they are)
gnss convert base.T04 -o work/

# 2. process: base, then the rover file(s)
gnss ppk work/base.24o flight1.obs \
    --base-llh 45.123456789,7.123456789,512.345 \
    --nav work/base.24n --nav work/base.24g --nav work/base.24l --nav work/base.24c \
    --antex igs20.atx \
    --mrk flight1_Timestamp.MRK \
    -o out/flight1.pos
```

This writes

| file | content |
|---|---|
| `out/flight1.pos` | the track: one line per epoch, in RTKLIB's `.pos` layout |
| `out/flight1_events.pos` | one line per event (photo), same layout, at the event's time |
| `out/flight1_events.csv` | one line per event with everything known about it |

and prints a summary: the share of fixed epochs, the agreement of the forward and backward passes, the events
with and without a position, and every warning.

## What goes in

### Observations

RINEX 2, 3 or 4 observation files, plain, gzip or Hatanaka-compressed. Receiver raw logs are converted first
with `gnss convert`.

- **Base**: a receiver that did not move, logging over the whole rover session. Its interval may be longer
  than the rover's: its observations are brought to the rover's epochs by interpolation (see
  [How it works](#how-it-works)). Up to `base_max_gap` seconds between base epochs are bridged (default 2.5 s,
  enough for a 1 Hz base under a 5 or 10 Hz rover). A base logging every 15 or 30 s can be used by raising
  `base_max_gap`, at some cost in accuracy.
- **Rover**: one file, or several. Some recorders start a new file for each flight while the receiver keeps
  tracking; give all the files of such a session and they are joined, so that the ambiguities found in one
  flight carry into the next. A short file on its own (a flight of a minute) has little time to resolve its
  ambiguities; joined to its neighbours it has none to resolve.

The two receivers must share satellites **and signals**. What is used is listed in the summary, for example
`G L1C/L1C (11 sats) L2W/L2X (10 sats)`: rover code / base code per frequency. The engine uses every
frequency both receivers track on a satellite (up to three), frequency by frequency, without combining them.
Where the receivers track the same frequency in different ways (`L1C` at one, `L1X` at the other), the
satellites concerned are differenced among themselves only, which removes the offset between the two tracking
modes exactly.

### The base coordinate

`--base-xyz X,Y,Z` (ECEF metres) or `--base-llh LAT,LON,H` (degrees, ellipsoidal height in metres): the
coordinate of the **mark** under the base antenna.

**Every rover position comes out in the reference frame and at the epoch of this coordinate.** A base
coordinate in a national frame gives rover positions in that frame; a coordinate from a precise point
positioning service in the current global frame gives positions in the global frame of that day. The two can
differ by more than a metre. Nothing in the processing can tell which was given: it is the user's choice, and
the most common cause of a constant offset against other data.

Without a base coordinate the approximate position in the base file's header is used and a warning says so:
the track has the right shape but sits wherever that approximate position (often metres off) puts it.

The antenna height (`ANTENNA: DELTA H/E/N` of the base file, or `--base-hen H,E,N`) is added to the mark to
get the antenna reference point.

### Antennas

- **Base**: give an ANTEX file (`--antex igs20.atx`) and have the antenna type in the base file's header. The
  phase-centre offsets and variations of the base antenna are then applied. Without them the base's phase
  centre is taken to be at its reference point, and every height is off by the antenna's phase-centre offset
  (often 5 to 15 cm). The engine warns when this happens.
- **Rover**: when the rover file names an antenna type that the ANTEX file has, it is applied and the positions
  are those of the rover's antenna reference point (less its `ANTENNA: DELTA H/E/N`, taken along the
  vertical, for a pole). Drone receivers usually name no antenna. The positions are then those of the
  **antenna's phase centre**, which is what the offset to the camera in a drone's event file starts from; the
  small difference between the phase centres of different frequencies (millimetres) is not modelled. A note
  in the output says which case applies.

### Orbits

`--nav` broadcast navigation files (several are merged; GLONASS needs its navigation file for the frequency
channels), and / or `--sp3` precise orbits with `--clk` clocks. Over the short baselines of kinematic work
broadcast orbits are as good as precise ones.

## What comes out

### The summary

```
BASE -> ROVR: 6000 epochs at 0.20 s, fixed 99.6% (5976 fixed, 24 float, 0 single, 0 none)
signals (rover/base): G L1C/L1C (10 sats) L2W/L2W (10 sats); E L1C/L1X (9 sats); ...
baseline 35.2 - 610.4 m; new ambiguities 212 (first appearances, breaks, slips); slips repaired 41
forward fixed 99.2%, backward 99.3%; both fixed at 5930 epochs: difference rms E 0.8 N 1.9 U 2.6 mm, max 21.0 mm; conflicts 0
events: 480 (479 with a position, 479 of them fixed; 0 across a gap; 1 without)
```

- **fixed / float / single / none**: how each epoch was solved (see the `Q` column below).
- **forward / backward**: the data are processed twice, forward and backward in time, and merged. Where both
  passes are fixed their positions are compared: the RMS and the largest difference are a direct check of
  the result, made from the data themselves. **Conflicts** are epochs where the passes disagreed beyond their
  uncertainty; such epochs are reported as float.
- **warnings** must be read. They name what limits the result: no antenna calibration, a base coordinate
  taken from a header, events outside the data, glitches in a file's epoch times.

### The `.pos` file

The layout of RTKLIB's solution files, so that scripts and programs made for those read these unchanged:
header lines starting with `%`, then one line per epoch:

```
%  GPST                  latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio
2024/02/08 00:00:49.200   45.123456789    7.123456789   612.3456   1  31   0.0050   0.0040   0.0100  -0.0010   0.0000   0.0020   0.00   42.3
```

| column | meaning |
|---|---|
| `GPST` | GPS time of the position. It is the time the signals were received: the epoch's time tag less the receiver's clock offset. For receivers that steer their clock the two are the same. |
| latitude, longitude, height | degrees and metres on the ellipsoid, in the frame of the base coordinate |
| `Q` | **1** fixed: integer ambiguities resolved and validated, centimetre level. **2** float: carrier-phase solution without resolved integers, decimetres, better as it converges. **5** single: a code-only position of the rover alone, metres (written when there is no base data for the epoch or too few common satellites). **6** precise point positioning (`gnss ppp --kinematic`). Epochs with no position at all are not written. |
| `ns` | satellites in the solution |
| `sdn sde sdu` | standard deviations north, east, up (m) |
| `sdne sdeu sdun` | square roots of the covariances between them, with the covariance's sign (RTKLIB's convention) |
| `age` | seconds between the rover epoch and the nearest base epoch used |
| `ratio` | the ratio of this epoch's independent ambiguity search (second-best candidate over best). See below. |

Options: `--time-decimals N` (default 3), `--pos-velocity` (adds `vn ve vu`, the Doppler velocity in m/s),
`--keep-passes` (also writes `NAME_forward.pos` and `NAME_backward.pos`).

**The standard deviations are meant to be realistic, not formal.** For a fixed epoch the formal covariance
of the least-squares solution is scaled by `sigma_fixed_scale` (2) and a floor is added (5 mm horizontal,
10 mm vertical, plus 1 ppm of the baseline), because the formal value knows nothing of multipath or of an
uncalibrated rover antenna. For a float epoch the filter's covariance is scaled by `sigma_float_scale` (3).

**About `ratio` on fixed epochs.** An epoch is fixed when its *confirmed* integers fit this epoch's phases
(see [How it works](#how-it-works)); the search that produces `ratio` is run afresh at every epoch from the
float solution alone, as an independent audit. During a manoeuvre that audit can be weak (a low ratio) while
the confirmed integers still fit to millimetres: the epoch is then fixed with a low ratio. A ratio is not
needed to trust a fixed epoch; a failed audit that persists removes the fix by itself.

### Event files

`NAME_events.pos` has the `.pos` layout with one line per event that got a position: the time is the event's
(six decimals by default), the position is the point wanted (the antenna position plus the event's offset,
when it has one), `Q` is the worst `Q` of the epochs it was interpolated between.

`NAME_events.csv` has one line for **every** event, with or without a position:

| column | meaning |
|---|---|
| `index`, `label` | position in the input list (from 1), and the event's label (photo sequence number; with several MRK files, `file:number`) |
| `gps_week`, `gps_tow_s`, `gpst` | the event's time |
| `status` | `ok`; `gap`: interpolated across missing epochs (allowed by `--max-missing`); `none`: no position, `why` says why |
| `fix`, `q` | the worst status of the epochs used (`fixed` 1, `float` 2, `single` 5) |
| `method` | the interpolation actually used (`cubic4` falls back to `linear2` in the first and last interval and beside a gap) |
| `latitude_deg`, `longitude_deg`, `height_m`, `x_m`, `y_m`, `z_m` | the point wanted (antenna plus offset) |
| `antenna_latitude_deg`, `antenna_longitude_deg`, `antenna_height_m` | the antenna at the event |
| `sd_e_m`, `sd_n_m`, `sd_u_m` | standard deviations |
| `gap_s`, `missing_epochs` | spacing of the two epochs around the event, and how many epochs are missing between them |
| `epochs_used`, `epochs_q` | indices of the epochs used (from 0, into the track) and their `Q` |
| `why` | why there is no position |

## Events

The solver positions a list of times. It never invents an event; the list comes from one or more of:

- `--mrk FILE` (repeatable): DJI `Timestamp.MRK` files;
- `--rinex-events`: the event records (epoch flag 5) of the rover file(s);
- `--events FILE`: a text file, one event per line, `WEEK SECONDS [LABEL]` or `GPS_SECONDS [LABEL]`.

All times are GPS time.

### Interpolation

An event falls between two epochs. `--interp` chooses how its position is found:

| option | uses | when |
|---|---|---|
| `linear2` | the two epochs around the event, straight line | always available; exact for constant velocity |
| `cubic4` (default) | two epochs before and two after, a cubic through the four positions | follows curves and speed changes; recommended at any rate |
| `hermite2` | the two epochs around the event **and their Doppler velocities** | low-rate data (1 Hz) with good Doppler |

Interpolation error on exact tracks, with no measurement noise (millimetres, largest over many events):

| track | rate | `linear2` | `cubic4` | `hermite2` (exact velocities) |
|---|---|---|---|---|
| straight, 15 m/s | 1 Hz | 0 | 0 | 0 |
| straight, 15 m/s | 5 Hz | 0 | 0 | 0 |
| turn, 30 m radius at 12 m/s (4.8 m/s²) | 1 Hz | 598 | 18 | 2 |
| turn, 30 m radius at 12 m/s (4.8 m/s²) | 5 Hz | 24 | 0.03 | 0.00 |
| speeding up at 2 m/s² | 1 Hz | 280 | 0 | 0 |
| speeding up at 2 m/s² | 5 Hz | 11 | 0 | 0 |

The linear error is (acceleration × interval²) / 8 in the middle of an interval; the cubic's falls with the
fourth power of the interval. At 5 Hz and above the cubic's own error is far below the noise of the positions
it interpolates, on any manoeuvre an aircraft can fly.

`hermite2` replaces the two outer positions by the velocities the receiver measured. With perfect velocities
it is the most exact of the three. Real Doppler velocities carry a few centimetres per second of noise, lag
behind in accelerations, and on some receivers are stale for an epoch after a loss of lock. Their effect on
the event grows with the interval (interval × 0.15 × velocity error at mid-interval): at 1 Hz, where the
cubic through positions cannot see what happens between epochs, the velocities help clearly; at 5 Hz the
4-point cubic is at least as good and does not depend on the Doppler. Hence the recommendation above. The
option needs Doppler in the rover file; where an epoch has no velocity the linear form is used for that
event and `method` says so.

### Rules that protect the result

- **Never extrapolated.** An event before the first or after the last epoch with a position gets none.
- **Missing epochs.** When one or more epochs are missing next to an event, the event gets no position unless
  `--max-missing N` allows that many missing epochs. Such an event has status `gap`, uses the epochs on either
  side of the hole, and its standard deviation grows with the longer span.
- The cubic needs four epochs with no more than the allowed gaps between them; where it cannot be formed the
  linear form is used and `method` says so.
- The standard deviation of an event is that of the epochs it rests on (combined as fully correlated:
  neighbouring epochs share their errors, averaging them gains nothing) plus an allowance for the
  interpolation itself, estimated from the curvature of the track around the event.

### Camera offsets (DJI `Timestamp.MRK`)

A line of the file:

```
1   345600.250000   [2300]   -120,N   45,E   193,V   45.00000001,Lat   7.00000002,Lon   512.345,Ellh   0.012000, 0.011000, 0.025000   50,Q
```

| field | meaning |
|---|---|
| `1` | photo number |
| `345600.250000`, `[2300]` | exposure time: GPS seconds of the week, GPS week |
| `-120,N` `45,E` `193,V` | the vector from the antenna phase centre to the centre of the camera sensor at that moment, in millimetres: towards north, towards east, and **downwards** (the camera hangs below the antenna, so `V` is positive) |
| `Lat` `Lon` `Ellh` | the position the aircraft computed in flight **for the camera** (antenna position plus the vector above), not for the antenna. It is of RTK quality only when `Q` is 50 |
| three numbers | its standard deviations north, east, vertical (m) |
| `Q` | how it was computed: 50 fixed RTK, 34 float RTK, 16 single point (no corrections received) |

The vector changes from photo to photo because it is the fixed lever arm between antenna and camera turned
by the aircraft's attitude and the gimbal at the exposure. Applying it needs no attitude data:

```
camera north  = antenna north + N
camera east   = antenna east  + E
camera height = antenna height − V
```

The engine does exactly this (in metres, in the local level frame at the event) and writes the camera
position in the event files, with the antenna position beside it in the CSV. `--no-offsets` leaves the
positions at the antenna.

The aircraft's own positions in the file are not used. They are useful as a check: with RTK corrections in
flight (`Q` 50) they should agree with the post-processed camera positions to a few centimetres, provided
both are in the same reference frame; without corrections (`Q` 16) they are single-point positions, good to
a metre or two, and sit in the frame of the broadcast orbits rather than in that of the base coordinate.

When the GPS week changes during a flight (Saturday to Sunday, GPS time) the seconds of the week restart at
zero. The reader handles it, also when the file's week number did not follow, and says so.

### Stop and go

`--occupations` reads the occupations marked in the rover file (epoch flag 3 "new site occupation" with its
`MARKER NAME`, until flag 2 "start moving antenna") and prints, for each, the mean of its fixed epochs (of its
float ones when none is fixed), their scatter, and a standard deviation that is not divided by the square
root of the number of epochs: the epochs of one occupation share most of their errors.

### Writing events into a RINEX file

```
gnss events insert rover.obs --mrk rover_Timestamp.MRK -o rover_events.obs
```

writes a copy of the observation file with an event record (`> 2024 02 08 00 00 49.2000000  5  0`) before the
first epoch that is later than each event, so that other software reads the same events. Everything else is
copied byte for byte. Events outside the file's epochs are not written and are counted in the output.

`gnss events list FILE` lists the events of an MRK file or of a RINEX file.

## Options

`--options FILE` takes a JSON object with any of the following fields; fields left out keep their default,
unknown fields are an error. The command-line switches `--systems`, `--elmask`, `--direction`, `--interp`,
`--max-missing` and `--no-fix` set the corresponding fields.

| field | default | meaning |
|---|---|---|
| `systems` | `"GRECJ"` | systems to use: G GPS, R GLONASS, E Galileo, C BeiDou, J QZSS |
| `exclude` | `[]` | satellites to leave out, e.g. `["G05", "E11"]` |
| `elmask` | 12 | elevation mask, degrees |
| `max_bands` | 3 | most frequencies used per satellite |
| `direction` | `"combined"` | `forward`, `backward` or `combined` |
| `keep_passes` | false | also return the forward and backward tracks |
| `err_phase`, `err_phase_el` | 0.003, 0.003 | phase noise of one receiver (m): σ² = a² + (b / sin elevation)² |
| `code_ratio` | 100 | code noise as a multiple of the phase noise |
| `code_corr_s` | 5 | how long code errors last (s); at shorter intervals the code is down-weighted accordingly |
| `efact_glonass` | 1.5 | extra noise factor for GLONASS |
| `trop_model` | `"niell"` | troposphere for both ends: `niell` (standard atmosphere), `gpt3`, `vmf3` (need grids) |
| `humidity` | 0.5 | relative humidity of the standard atmosphere, for the wet delay |
| `pcv_azimuth` | true | use the azimuth-dependent part of antenna calibrations |
| `ar` | `"continuous"` | integer ambiguity resolution; `off` gives float solutions only |
| `ar_ratio` | 3.0 | least ratio of the ambiguity search |
| `ar_elmask` | 15 | elevation above which a satellite's ambiguities are searched, degrees |
| `ar_min_lock_s` | 2 | seconds without a slip before an ambiguity is searched |
| `ar_min_dd` | 5 | least number of double-difference ambiguities for a fix |
| `ar_max_res` | 0.02 | largest RMS of the phase residuals of a fixed solution (m); a single phase may misfit by three times this |
| `ar_confirm_s` | 1.0 | seconds over which an integer must keep coming out the same before it is relied on |
| `ar_max_shift_sigma` | 5 | a fixed position farther than this many (realistic) float σ plus 5 cm from the float one is not accepted |
| `ar_partial_rounds` | 3 | rounds of the partial search when the full set fails |
| `ar_glonass` | false | fix GLONASS ambiguities (see [Limits](#limits)) |
| `amb_sigma0`, `amb_rw` | 30, 1e-4 | σ of a new ambiguity (cycles) and its random walk (cycles per √s) |
| `pos_sigma0` | 30 | σ of the position before each epoch's update (m) |
| `max_outage_s` | 1 | a satellite out of the data for longer (s) starts a new ambiguity; with data at a long interval, two and a half intervals are allowed instead |
| `slip_gf`, `slip_td` | 0.05, 0.05 | cycle-slip limits (m), see [How it works](#how-it-works) |
| `outlier_k`, `code_gross` | 5, 30 | outlier limits: standard deviations, and metres for the code |
| `base_max_gap`, `base_tol` | 2.5, 0.0005 | base interpolation: largest gap bridged, and what counts as the same epoch (s) |
| `max_baseline_km` | 20 | beyond this distance a warning says the model is being stretched |
| `min_sats` | 5 | fewest satellites for a phase solution |
| `velocity` | true | compute Doppler velocities |
| `event_interp` | `"cubic4"` | `linear2`, `cubic4`, `hermite2` |
| `event_max_missing` | 0 | missing epochs next to an event that may be skipped |
| `sigma_fixed_scale`, `sigma_fixed_h`, `sigma_fixed_v`, `sigma_ppm` | 2, 0.005, 0.010, 1 | realistic σ of a fixed epoch |
| `sigma_float_scale`, `sigma_single` | 3, 3 | realistic σ of float and code-only epochs |

## How it works

**Observations.** For every satellite and frequency both receivers track, the difference between rover and
base is formed (single difference). It removes the satellite's clock and, over a short baseline, its orbit
error and the ionosphere. The troposphere is modelled at both ends, each at its own height: a drone 100 m
above its base sees centimetres less of it. The base's observations are not interpolated as such: what is
interpolated to the rover's epoch is their *correction* (observed minus modelled), which changes slowly
except for the base receiver's clock, and that is the same for every satellite and drops out of the
differences between satellites whatever it does between two base epochs.

**States.** The rover position and one ambiguity per satellite and frequency. The position carries nothing
from one epoch to the next: before each update it is set to a code position with a loose prior. No
assumption is made about how the receiver moves, so no manoeuvre can violate one. The ambiguities are
constants until a cycle slip.

**Cycle slips.** A new ambiguity is started

- when a satellite appears, or returns after more than `max_outage_s`;
- when a receiver's geometry-free phase jumps (two frequencies; limit `slip_gf`);
- when the change of a phase from one epoch to the next does not fit the motion that all other satellites
  agree on (three position components and a clock are solved from all phases; limit `slip_td`). This works
  on a single frequency and for any motion.

A loss-of-lock flag from the receiver, or a break of a few epochs in one satellite's phase, does not by
itself start a new ambiguity. The phase is tested against the motion of the others since it was last seen:
no jump means it was whole; a jump of a whole number of cycles (to within 0.15 cycle) is **taken out** and
the ambiguity carries on; anything else starts a new one. Receivers under dynamic stress (a hard brake)
lose and regain the highest satellites within half a second; without the repair those satellites would be
missing from the solution for seconds exactly when the geometry needs them most.

**Float solution.** Each epoch, double differences against the highest satellite of every signal group are
formed with their full covariance and the states are updated. A phase that still misfits afterwards slipped
unnoticed: its ambiguity is restarted and the update repeated. Code errors last several seconds (receiver
smoothing, multipath), so at high rates the code is down-weighted to keep the float solution honest about
what it knows.

**Integers.** The search is LAMBDA with the ratio test, on all candidate double differences, afresh at
every epoch from the float solution. An integer that comes out the same at least three times over
`ar_confirm_s` is *confirmed*. From then on:

1. the **fixed position** of an epoch is computed from that epoch's phases alone with the confirmed integers
   taken as known (a small least-squares problem: three coordinates and a clock per signal group). It has no
   memory: two passes that hold the same integers give the same position;
2. every phase must fit that position (each to three times `ar_max_res`, their RMS to `ar_max_res`); a
   satellite that does not loses its integer;
3. the float solution audits the confirmed integers at every epoch: if a fresh search would prefer another
   set by the ratio a fix needs, the epoch is not fixed, and when that lasts `ar_confirm_s` all integers are
   given up and resolution starts again;
4. a new or restarted satellite gets its integer directly from the fixed position (with position and clock
   known, one phase gives its integer), and has to repeat it over `ar_confirm_s` before it counts.

The float solution is never overwritten by a fix, so a wrong set of integers cannot feed back into it.

**Forward, backward, combined.** The same filter runs over the data in both directions. Where both passes
are fixed and agree, the result is their mean weighted by their covariances; where they disagree beyond
their uncertainty the epoch is reported float with a σ widened to the disagreement; where one is fixed and
the other float, the fixed one is taken if it lies inside the float one's uncertainty. Epochs before the
first fix of the forward pass are covered by the backward pass.

**GLONASS** satellites take part with float ambiguities by default: their frequencies differ from satellite
to satellite, and their integers can only be relied on between receivers of the same make.

## Without a base: kinematic PPP

`gnss ppp OBS --kinematic` positions a moving receiver on its own, with an analysis centre's precise orbits
and clocks instead of a base receiver:

```
gnss ppp flight1.obs --kinematic \
    --sp3 COD0OPSFIN_20240390000_01D_05M_ORB.SP3 --clk COD0OPSFIN_20240390000_01D_30S_CLK.CLK \
    --osb COD0OPSFIN_20240390000_01D_01D_OSB.BIA --erp COD0OPSFIN_20240390000_01D_01D_ERP.ERP \
    --antex igs20.atx --nav brdc0390.24g \
    --mrk flight1_Timestamp.MRK -o out/flight1_ppp.pos
```

The models are those of the static `gnss ppp` (ionosphere-free phase and code on two frequencies, satellite
antennas, code and phase biases, tides, wind-up, a zenith wet delay estimated as a random walk); the estimate
is one position per epoch, forward, backward and merged. Events (`--mrk`, `--rinex-events`, `--events`,
`--interp`, `--max-missing`) and the outputs (`.pos`, `_events.pos`, `_events.csv`) work as for `gnss ppk`.
Epochs have `Q` 6. The positions are in the frame of the orbits (the current IGS frame at the day of the
data), not in the frame of any base.

**What to expect.** A float precise point solution has to converge: its ambiguities start from the code,
metres off, and are separated from the position only as the satellites move across the sky. A static
station processed this way, one position per 30 s epoch with nothing assumed about its motion, gives:

| session | merged track, 95 % of epochs within |
|---|---|
| 24 h | 6 cm horizontal, 13 cm vertical (after the first 12 minutes) |
| 2 h | 20 cm horizontal, 18 cm vertical |
| 50 min | 20 cm horizontal, 30 to 50 cm vertical |

The forward pass alone needs about half an hour to come under 30 cm and one to two hours to come under
10 cm; the merged track is far better at the start because the backward pass arrives there converged.
For a flight of some tens of minutes this is a decimetre-level method. Where a base is available,
`gnss ppk` is the centimetre-level one. Needs two frequencies on each satellite used; the receiver must
log before and after the part that matters, the longer the better.

`--kinematic-options FILE` takes a JSON object with the filter's settings (`pos_sigma0`, `clk_sigma`,
`amb_sigma0`, `isb_sigma0`, `isb_rw`, `code_corr_s`, `outlier_k`, `min_arc_s`, `min_sats`, `sigma_scale`);
`--options` sets the models, as for the static solution.

## Checks on the input files

`gnss events check FILE` (and `gnss ppk`, for every rover file) looks at the epoch times for what recorders
are known to get wrong, and reports without changing anything:

- epochs dated 1980 at the start of a file, written before the receiver knew the date;
- a run of epochs whose date is some days off while the time of day continues (seen shortly before a GPS
  week change);
- repeated epochs;
- a `TIME OF FIRST OBS` in the header that is not the first epoch.

The first two make a file unusable as written: correct the dates in a copy of the file before processing.

A rover file whose Doppler has the opposite sign of the RINEX definition is recognised from the file's own
phases, used with the sign turned, and reported.

## Bases from a correction stream

A base that was logged as an RTCM 3 stream rather than as a receiver file is converted first:

```
gnss convert base.rtcm3 -o work/ --date 2024-02-08
```

and the resulting RINEX file is the base of `gnss ppk`. The station position and antenna height of the
stream's 1005 / 1006 message go into the RINEX header (`APPROX POSITION XYZ`, `ANTENNA: DELTA H/E/N`), and the
antenna name of its 1033 message into `ANT # / TYPE`; the base coordinate is still given to `gnss ppk`
explicitly, because the frame it is in is the user's to know. See [CLI.md](CLI.md#gnss-convert-raw--o-out).

## Python

The bindings (`geomaticist_gnss_py`) expose the same functions:

```python
import json, struct
import geomaticist_gnss_py as g

def station(h):
    return {"handle": h, "name": h.name, "antenna": h.antenna, "receiver": h.receiver,
            "approx": h.approx, "delta_hen": h.delta_hen}

base = g.read_obs("base.24o")
rover = g.read_obs("flight1.obs", kinds="CLSD")            # keep the Doppler
nav = g.read_nav(["base.24n", "base.24g", "base.24l", "base.24c"])
first_epoch = struct.unpack("<d", rover.t_bytes()[:8])[0]
mrk = json.loads(g.read_mrk("flight1_Timestamp.MRK", first_epoch))

out = json.loads(g.ppk_solve(
    station(base), [station(rover)], [("broadcast_handle", nav)],
    base_xyz,                                               # ECEF metres of the base mark
    options_json='{"event_interp": "cubic4"}',
    antex="igs20.atx", fcn=nav.glo_fcn(),
    events=[tuple(e) for e in mrk["events"]],               # (gps_seconds, label, (north, east, down) or None)
    progress=lambda stage, done, total: True,               # return False to cancel
    outputs="epochs,pos,events_pos,events_csv"))

print(out["summary"]["fix_percent"], out["warnings"])
open("flight1.pos", "w", newline="").write(out["pos"])
```

`ppk_solve` returns JSON: `summary`, `warnings`, `notes`, `events` (one object per event), `occupations`,
and what `outputs` asks for: `epochs` (the track as columns `t`, `lat`, `lon`, `h`, `q`, `sd_e` ...), `passes`,
and the texts `pos`, `events_pos`, `events_csv`. Also available: `read_mrk`, `rinex_events`,
`rinex_occupations`, `check_times`, `insert_events`.

## Limits

- **Short baselines.** The model takes the ionosphere to be the same at both receivers. That holds to
  centimetres over a few kilometres and degrades with distance and ionospheric activity; beyond
  `max_baseline_km` (20) a warning is given. A solution for long baselines (ionosphere estimated) is not
  part of this version.
- **Static base.** A moving base is not supported.
- **GLONASS** ambiguities are float unless `ar_glonass` is set; set it only for two receivers of the same
  type.
- **Tides** are not applied: over these distances the base and the rover move together.
- **Memory**: every epoch is prepared once and kept for the two passes, about 5 kB per epoch (a one-hour
  flight at 10 Hz: about 200 MB).

## Credits

The `.pos` layout and the quality flags are those of RTKLIB (T. Takasu; BSD 2-Clause), reproduced from its
published format so that existing tools keep working; the combination of a forward and a backward pass
follows the same practice. No RTKLIB code is included. The integer search is the LAMBDA method (P. Teunissen
and co-workers). The troposphere, orbit, antenna and geometry models are the ones the static engine of this
repository uses.
