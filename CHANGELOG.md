# Changelog

All notable changes to this project are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## Unreleased

First public release.

- Engine: static baselines, static and kinematic PPP (optional PPP-AR), kinematic PPK with camera events
  and stop-and-go occupations.
- RINEX 2/3/4 reading (plain, gzip, Hatanaka), trimming, decimation and conversion; event and
  kinematic-marking tools for office software.
- Raw receiver decoders: Trimble T04/T02/T01/T00, Leica MDB, Topcon/Javad GREIS, RTCM 3, u-blox UBX,
  NovAtel, Septentrio SBF.
- The `gnss` command line and the `geomaticist_gnss_py` Python bindings, with commented examples.
- Sample data: one hour from two NOAA CORS stations (GODE, GODN) with the day's broadcast orbits.
