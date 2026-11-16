"""Precise point positioning (PPP): one receiver, no base station, precise orbits and clocks.

    python ppp.py                                         # synthetic: static PPP of one station
    python ppp.py --kinematic                             # synthetic: kinematic PPP of a drone flight
    python ppp.py OBS --sp3 ORB.sp3 --clk CLK.clk [--antex igs20.atx] [--osb BIAS.bia] [--erp ERP.erp]
        [--nav BRDC.rnx] [--blq OCEAN.blq --otl] [--ar] [--systems GRE]
        [--kinematic [--mrk Timestamp.MRK | --rinex-events] [--out track.pos]]

Static PPP gives one position for the whole session with a realistic sigma (scaled for the session
length): hours of data give centimetres, a short session decimetres. PPP-AR (``--ar``) fixes the
ambiguities and needs phase biases from the analysis centre (``--osb``).

Kinematic PPP gives one position per epoch. With camera events (``--mrk`` or the rover's flag-5 records)
it interpolates the track at each photo, as ``ppk_events_to_points.py`` does with a base station. PPP needs
no base, but it converges over tens of minutes, so give it a long static start before the flight.

Products: the IGS analysis centres publish SP3 orbits, CLK clocks, Bias-SINEX and ERP files for every day
(final products about two weeks later, rapid ones the next day); igs20.atx is the IGS antenna file.
"""

import argparse
import json
import math
import os

import geomaticist_gnss_py as g

from _common import ensure_out, epochs, station, synthetic


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("obs", nargs="?")
    ap.add_argument("--sp3", nargs="*", default=[])
    ap.add_argument("--clk", nargs="*", default=[])
    ap.add_argument("--antex")
    ap.add_argument("--osb", nargs="*", default=[])
    ap.add_argument("--erp", nargs="*", default=[])
    ap.add_argument("--nav", nargs="*", default=[], help="broadcast files, for the GLONASS channels only")
    ap.add_argument("--blq")
    ap.add_argument("--otl", action="store_true", help="apply the ocean loading of --blq")
    ap.add_argument("--ar", action="store_true", help="PPP-AR (needs phase biases in --osb)")
    ap.add_argument("--systems", default="GRE")
    ap.add_argument("--kinematic", action="store_true")
    ap.add_argument("--mrk", help="kinematic: DJI Timestamp.MRK with the photo events")
    ap.add_argument("--rinex-events", action="store_true", help="kinematic: the flag-5 records of OBS")
    ap.add_argument("--out", help="kinematic: write the track .pos here (default out/ppp_track.pos)")
    a = ap.parse_args()

    truth = None
    if a.obs is None:
        # synthetic: the simulated SP3 carries the satellite clocks too, so no CLK file is needed
        d = synthetic("flight" if a.kinematic else "static")
        if a.kinematic:
            a.obs, a.mrk, truth = d["rover"], d["mrk"], d["truth"]
        else:
            a.obs, truth = d["stations"]["PT01"], d["truth"]["PT01"]
        a.sp3 = [d["sp3"]]
        a.systems = "G"  # the simulated Galileo and BeiDou are single-frequency
    elif not a.sp3:
        ap.error("give --sp3 (and usually --clk)")

    obs = g.read_obs(a.obs, kinds="CLSD" if a.kinematic else "CLS")
    opts = {"systems": a.systems, "ar": a.ar, "otl": a.otl}

    events = None
    if a.kinematic and a.mrk:
        events = json.loads(g.read_mrk(a.mrk, epochs(obs)[0]))["events"]
    elif a.kinematic and a.rinex_events:
        events = json.loads(g.rinex_events(a.obs))["events"]

    out = json.loads(
        g.ppp_solve(
            station(obs), a.sp3, clk=a.clk, antex=a.antex, osb=a.osb, erp=a.erp, nav=a.nav, blq=a.blq,
            options_json=json.dumps(opts), kinematic=a.kinematic, events=events,
            outputs="pos,events_csv" if a.kinematic else "",
        )
    )

    if not a.kinematic:
        r = out["result"]
        lat, lon, h = out["llh_deg"]
        sig = [math.sqrt(max(out["cov_realistic_enu"][i][i], 0.0)) for i in range(3)]
        print(f"{obs.name}: {'fixed' if r['fixed'] else 'float'} PPP, {r['hours']:.2f} h, {r['n_epochs']} epochs")
        print(f"  X Y Z     {r['xyz'][0]:14.4f} {r['xyz'][1]:14.4f} {r['xyz'][2]:14.4f} m (antenna reference point)")
        print(f"  lat lon h {lat:14.9f} {lon:14.9f} {h:10.4f}")
        print(f"  sigma ENU {sig[0]:.3f} {sig[1]:.3f} {sig[2]:.3f} m (realistic)")
        print(f"  rms phase {r['rms_phase'] * 1000:.1f} mm, code {r['rms_code']:.2f} m, rejected {r['n_rejected']}")
        for n in r["notes"]:
            print(f"  note: {n}")
        if truth:
            print(f"  vs truth: {math.dist(r['xyz'], truth) * 100:.1f} cm")
        return

    s = out["summary"]
    print(f"{obs.name}: kinematic PPP, {s['n_epochs']} epochs, {s['n_float']} with a position")
    for w in out["warnings"]:
        print(f"  warning: {w}")
    path = a.out or os.path.join(ensure_out(), "ppp_track.pos")
    with open(path, "w", newline="\n") as f:
        f.write(out["pos"])
    print(f"track  -> {path}")
    if out["events"]:
        ev_path = path.rsplit(".", 1)[0] + "_events.csv"
        with open(ev_path, "w", newline="\n") as f:
            f.write(out["events_csv"])
        ok = sum(1 for e in out["events"] if e["llh"] and e["llh"][0] is not None)
        print(f"events -> {ev_path} ({ok} of {len(out['events'])} with a position)")
        if truth:
            errs = [
                math.dist(e["antenna_xyz"], t["antenna_xyz"])
                for e, t in zip(out["events"], truth)
                if e["antenna_xyz"] and e["antenna_xyz"][0] is not None
            ]
            if errs:
                print(f"vs truth: worst antenna position error {max(errs):.2f} m (PPP, not yet converged)")


if __name__ == "__main__":
    main()
