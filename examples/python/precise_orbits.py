"""Precise orbits (SP3): satellite positions and clocks, and a sky plot table for a station.

    python precise_orbits.py                         # synthetic SP3 and station
    python precise_orbits.py ORBITS.sp3 [CLOCKS.clk ...] --at X,Y,Z [--sat G05]

``read_sp3`` returns a ``Precise`` handle: ``span`` (first, last, count of the SP3 epochs, GPS seconds),
``sats`` and ``position(sat, [t, ...])`` -> ``[(x, y, z, clock_s, ok), ...]`` (ECEF metres, interpolated;
``ok`` is false outside the file's span). The same handle is an orbit source for the solvers:
``("precise_handle", handle)``.

Azimuth and elevation are worked out here from the station's ECEF position (``ecef_to_llh`` gives its
latitude and longitude for the local east-north-up frame).
"""

import argparse
import math

import geomaticist_gnss_py as g

from _common import gpst, synthetic, xyz_arg


def az_el(station_xyz, sat_xyz):
    """Azimuth and elevation (degrees) of a satellite seen from a station, both ECEF metres."""
    lat, lon, _ = g.ecef_to_llh(station_xyz)
    la, lo = math.radians(lat), math.radians(lon)
    d = [sat_xyz[i] - station_xyz[i] for i in range(3)]
    e = -math.sin(lo) * d[0] + math.cos(lo) * d[1]
    n = -math.sin(la) * math.cos(lo) * d[0] - math.sin(la) * math.sin(lo) * d[1] + math.cos(la) * d[2]
    u = math.cos(la) * math.cos(lo) * d[0] + math.cos(la) * math.sin(lo) * d[1] + math.sin(la) * d[2]
    return math.degrees(math.atan2(e, n)) % 360.0, math.degrees(math.atan2(u, math.hypot(e, n)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sp3", nargs="?")
    ap.add_argument("clk", nargs="*", default=[])
    ap.add_argument("--at", type=xyz_arg, help="station ECEF X,Y,Z metres")
    ap.add_argument("--sat", help="satellite to tabulate, e.g. G05 (default: the first)")
    a = ap.parse_args()
    if a.sp3 is None:
        d = synthetic("static")
        a.sp3, a.at = d["sp3"], d["base_xyz"]
    elif a.at is None:
        ap.error("give --at X,Y,Z")

    p = g.read_sp3([a.sp3], a.clk)
    first, last, count = p.span
    sats = p.sats
    print(f"{a.sp3}\n  {count} epochs, {gpst(first)} .. {gpst(last)} GPS time; {len(sats)} satellites")
    print(f"  {' '.join(sats)}")

    lat, lon, h = g.ecef_to_llh(a.at)
    print(f"\nsky at lat {lat:.5f} lon {lon:.5f} h {h:.1f} m, {gpst(first)}:")
    print("  sat   azimuth  elevation")
    for sat in sats:
        x, y, z, _, ok = p.position(sat, [first])[0]
        if ok:
            az, el = az_el(a.at, (x, y, z))
            if el > 0:
                print(f"  {sat}  {az:7.1f}  {el:9.1f}")

    sat = a.sat or sats[0]
    print(f"\n{sat} every 15 minutes (ECEF km, clock microseconds, elevation degrees):")
    ts = [first + k * 900.0 for k in range(int((last - first) // 900) + 1)]
    for t, (x, y, z, clk, ok) in zip(ts, p.position(sat, ts)):
        if ok:
            print(f"  {gpst(t, 0)}  {x / 1e3:11.3f} {y / 1e3:11.3f} {z / 1e3:11.3f}  {clk * 1e6:9.3f}"
                  f"  {az_el(a.at, (x, y, z))[1]:6.1f}")  # fmt: skip

    # the round trip between ECEF and latitude / longitude / ellipsoidal height (WGS 84)
    back = g.llh_to_ecef(lat, lon, h)
    print(f"\nllh_to_ecef(ecef_to_llh(station)) differs by {math.dist(back, a.at) * 1e3:.6f} mm")


if __name__ == "__main__":
    main()
