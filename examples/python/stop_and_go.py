"""Stop-and-go: one kinematic rover file with several short occupations -> one point per stop.

    python stop_and_go.py                                   # synthetic: three stops PT01..PT03
    python stop_and_go.py BASE.obs ROVER.obs --base-xyz X,Y,Z --nav BRDC.rnx [--sp3 ...]

The occupations come from the rover RINEX's special records: flag 3 (new site, with a MARKER NAME line)
starts a stop, flag 2 (start moving) ends it. Survey controllers and the office software write these.
If your file has none, build the list yourself: ``[(name, gps_t0, gps_t1), ...]``.

The engine averages the fixed epochs inside each occupation; ``sigma_enu`` is the formal sigma of that
mean and ``scatter_enu`` the spread of the epochs about it.
"""

import argparse
import json
import math

import geomaticist_gnss_py as g

from _common import gpst, orbit_sources, station, synthetic, xyz_arg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", nargs="?")
    ap.add_argument("rover", nargs="?")
    ap.add_argument("--base-xyz", type=xyz_arg)
    ap.add_argument("--nav", nargs="*", default=[])
    ap.add_argument("--sp3", nargs="*", default=[])
    ap.add_argument("--clk", nargs="*", default=[])
    ap.add_argument("--antex")
    a = ap.parse_args()

    truth = None
    if a.base is None:
        d = synthetic("stop_and_go")
        a.base, a.rover, a.base_xyz = d["base"], d["rover"], d["base_xyz"]
        a.sp3, truth = [d["sp3"]], d["truth"]
    elif a.rover is None or a.base_xyz is None:
        ap.error("give BASE, ROVER and --base-xyz")

    occ = json.loads(g.rinex_occupations(a.rover))  # [[name, t0, t1], ...]
    if not occ:
        raise SystemExit("no flag-3 occupations in the rover file; pass your own list to ppk_solve")
    for name, t0, t1 in occ:
        print(f"occupation {name:<8} {gpst(t0)} .. {gpst(t1)}  ({t1 - t0:.0f} s)")

    sources, fcn = orbit_sources(a.nav, a.sp3, a.clk)
    out = json.loads(
        g.ppk_solve(
            station(g.read_obs(a.base)), station(g.read_obs(a.rover, kinds="CLSD")), sources, a.base_xyz,
            antex=a.antex, fcn=fcn, occupations=occ, outputs="",
        )
    )
    print(f"PPK: {out['summary']['fix_percent']:.1f} % of epochs fixed\n")

    print(f"{'name':<8} {'latitude':>14} {'longitude':>14} {'h (ellips.)':>11}  sigma E/N/U mm  epochs  fix")
    for o in out["occupations"]:
        lat, lon, h = o["llh"]
        se, sn, su = (x * 1000 for x in o["sigma_enu"])
        line = f"{o['name']:<8} {lat:14.9f} {lon:14.9f} {h:11.4f}  {se:4.1f} {sn:4.1f} {su:4.1f}  {o['n']:6d}  {o['fix']}"
        if truth and o["name"] in truth:
            line += f"   vs truth {math.dist(o['xyz'], truth[o['name']]) * 1000:.1f} mm"
        print(line)


if __name__ == "__main__":
    main()
