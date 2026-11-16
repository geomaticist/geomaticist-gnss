"""PPK with camera events: base + moving rover + event marks -> a point CSV (one row per photo).

    python ppk_events_to_points.py                                      # synthetic drone flight
    python ppk_events_to_points.py BASE.obs ROVER.obs --base-xyz X,Y,Z --nav BRDC.rnx
        [--mrk Timestamp.MRK] [--prefix P --start 1001] [--out points.csv]

The events come from a DJI ``Timestamp.MRK`` (``--mrk``; it also carries each photo's lever arm from the
antenna to the camera) or, without ``--mrk``, from the rover RINEX's flag-5 event records. Each point is
named ``<prefix><number>`` counting up from ``--start``, or by the event's own label when no prefix is
given. The engine interpolates the PPK track at each event time and applies the lever arm, so the
coordinates are at the camera (``xyz``/``llh``); ``antenna_xyz`` is the antenna phase centre.

Also written: the whole trajectory as an RTKLIB-style ``.pos`` (see README for the columns).
"""

import argparse
import csv
import json
import math
import os

import geomaticist_gnss_py as g

from _common import ensure_out, epochs, gpst, orbit_sources, station, synthetic, xyz_arg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", nargs="?")
    ap.add_argument("rover", nargs="?")
    ap.add_argument("--base-xyz", type=xyz_arg)
    ap.add_argument("--nav", nargs="*", default=[])
    ap.add_argument("--sp3", nargs="*", default=[])
    ap.add_argument("--clk", nargs="*", default=[])
    ap.add_argument("--antex")
    ap.add_argument("--mrk", help="DJI Timestamp.MRK; without it the rover's flag-5 events are used")
    ap.add_argument("--prefix", default="", help="point name prefix; empty = use the event labels")
    ap.add_argument("--start", type=int, default=1, help="first point number (with --prefix)")
    ap.add_argument("--out", help="CSV to write (default out/ppk_points.csv)")
    a = ap.parse_args()

    truth = None
    if a.base is None:
        d = synthetic("flight")
        a.base, a.rover, a.base_xyz = d["base"], d["rover"], d["base_xyz"]
        a.sp3, a.mrk, truth = [d["sp3"]], d["mrk"], d["truth"]
        a.prefix = a.prefix or "IMG"
    elif a.rover is None or a.base_xyz is None:
        ap.error("give BASE, ROVER and --base-xyz")

    base = g.read_obs(a.base)
    # kinds "CLSD" also keeps Doppler, which the kinematic filter uses for velocity
    rover = g.read_obs(a.rover, kinds="CLSD")

    # 1. the events: [GPS seconds, label, [north, east, down] lever arm in metres or None]
    if a.mrk:
        # MRK lines may lack the GPS week; `near` (the rover's first epoch) supplies it
        m = json.loads(g.read_mrk(a.mrk, epochs(rover)[0]))
        events = m["events"]
        print(f"{len(events)} events from {os.path.basename(a.mrk)}; RTK quality flags {m['quality']}")
    else:
        m = json.loads(g.rinex_events(a.rover))
        events = m["events"]
        print(f"{len(events)} flag-5 events in {os.path.basename(a.rover)}")
    for w in m["warnings"]:
        print(f"  warning: {w}")

    # 2. the solve (forward + backward Kalman filter, combined); outputs asks for the .pos text too
    sources, fcn = orbit_sources(a.nav, a.sp3, a.clk)
    out = json.loads(
        g.ppk_solve(
            station(base), station(rover), sources, a.base_xyz, "{}",
            antex=a.antex, fcn=fcn, events=events, outputs="pos",
        )
    )
    s = out["summary"]
    print(
        f"PPK: {s['n_epochs']} epochs, {s['fix_percent']:.1f} % fixed, baseline {s['baseline_m'][0]:.0f}-"
        f"{s['baseline_m'][1]:.0f} m; events ok {s['events_ok']}, failed {s['events_failed']}, "
        f"in a gap {s['events_gap']}"
    )
    for w in out["warnings"]:
        print(f"  warning: {w}")

    # 3. the point CSV
    path = a.out or os.path.join(ensure_out(), "ppk_points.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "gpst", "latitude_deg", "longitude_deg", "ellipsoid_height_m",
                    "sd_e_m", "sd_n_m", "sd_u_m", "fix", "status"])  # fmt: skip
        for k, e in enumerate(out["events"]):
            name = f"{a.prefix}{a.start + k}" if a.prefix else e["label"]
            if e["llh"] is None or e["llh"][0] is None:  # no position (gap, outside the track)
                w.writerow([name, gpst(e["t"], 3), "", "", "", "", "", "", e["fix"], f"{e['status']}: {e['why']}"])
                continue
            lat, lon, h = e["llh"]
            se, sn, su = e["sigma_enu"]
            w.writerow([name, gpst(e["t"], 3), f"{lat:.9f}", f"{lon:.9f}", f"{h:.4f}",
                        f"{se:.4f}", f"{sn:.4f}", f"{su:.4f}", e["fix"], e["status"]])  # fmt: skip
    print(f"points  -> {path}")

    pos = os.path.join(os.path.dirname(path), "ppk_track.pos")
    with open(pos, "w", newline="\n") as f:
        f.write(out["pos"])
    print(f"track   -> {pos}")

    if truth:  # synthetic data: antenna positions against the simulated ones
        worst = max(
            math.dist(e["antenna_xyz"], tr["antenna_xyz"]) for e, tr in zip(out["events"], truth)
        )
        print(f"vs truth: worst antenna position error {worst * 1000:.1f} mm")


if __name__ == "__main__":
    main()
