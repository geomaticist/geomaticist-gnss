# geomaticist-gnss-py

Python bindings for the Geomaticist GNSS engine (Rust): survey-grade GNSS post-processing from Python.

- Read RINEX 2/3/4 observation and navigation files (plain, gzip, Hatanaka), SP3 orbits and RINEX clocks.
- Decode raw receiver logs to RINEX 3: Trimble T04/T02/T01/T00, Leica MDB, Topcon/Javad GREIS, RTCM 3,
  u-blox UBX, NovAtel.
- Static baselines with ambiguity fixing and quality control.
- Kinematic PPK: forward/backward filter, camera events (DJI `Timestamp.MRK` or RINEX flag 5) with lever
  arms, stop-and-go occupations, RTKLIB-style `.pos` output.
- PPP, static and kinematic, with precise products (SP3/CLK, Bias-SINEX, ERP, ANTEX, BLQ).

```python
import json, geomaticist_gnss_py as g

obs = g.read_obs("SITE00USA_R_20240010000_01H_30S_MO.crx.gz")
out = json.loads(g.ppp_solve(
    {"handle": obs, "name": obs.name, "antenna": obs.antenna, "receiver": obs.receiver,
     "approx": obs.approx, "delta_hen": obs.delta_hen},
    sp3=["orbits.sp3"], clk=["clocks.clk"], antex="igs20.atx"))
print(out["llh_deg"])
```

Install: `pip install maturin`, then `maturin develop --release` in this folder (CPython 3.13+, Rust
toolchain). Commented examples for every function are in
[`examples/python`](../../examples/python/README.md).

Licensed under either of MIT or Apache-2.0, at your option.
