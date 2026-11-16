"""One static baseline: base + rover RINEX -> vector, sigmas, fixed/float and a short QC summary.

    python static_baseline.py                                   # synthetic data
    python static_baseline.py BASE.obs ROVER.obs --base-xyz X,Y,Z --nav BRDC.rnx [--sp3 ORB.sp3 ...]
        [--antex igs20.atx] [--options '{"elmask": 15}']

--base-xyz is the base MARK in ECEF metres (the antenna height from the base file's header, or the
station dict's "delta_hen", is applied by the engine). The rover's coordinates are those of its mark too.
"""

import argparse
import json
import math

import geomaticist_gnss_py as g

from _common import orbit_sources, station, synthetic, xyz_arg


def solve(base_path, rover_path, base_xyz, sources, fcn, antex=None, options="{}"):
    """Runs one baseline and returns the engine's result dict (see README: static result fields)."""
    base = g.read_obs(base_path)
    rover = g.read_obs(rover_path)
    out = json.loads(
        g.static_solve(station(base), station(rover), sources, base_xyz, options, antex=antex, fcn=fcn)
    )
    return out["result"]


def report(r, name="rover"):
    """Prints the vector, its sigmas, the fix and a QC summary."""
    dx, dy, dz = r["vector"]
    de, dn, du = r["enu"]
    sig = [math.sqrt(r["cov_enu"][i][i]) for i in range(3)]
    fix = r["fix"]
    lat, lon, h = g.ecef_to_llh(r["rover_xyz"])
    print(f"baseline to {name}: {math.sqrt(dx * dx + dy * dy + dz * dz):.4f} m")
    print(f"  dX dY dZ  : {dx:12.4f} {dy:12.4f} {dz:12.4f} m (ECEF)")
    print(f"  dE dN dU  : {de:12.4f} {dn:12.4f} {du:12.4f} m (local at the base)")
    print(f"  sigma ENU : {sig[0] * 1000:.1f} {sig[1] * 1000:.1f} {sig[2] * 1000:.1f} mm (realistic, 1 sigma)")
    print(f"  {name:<9} : lat {lat:.9f}  lon {lon:.9f}  h {h:.4f} m (WGS 84 ellipsoid)")
    print(f"  solution  : {'FIXED' if r['fixed'] else 'FLOAT'} ({fix.get('status')}, ratio {fix.get('ratio')})")
    if r["failed"]:
        print(f"  !! not to be trusted: {fix.get('why')}")
    # QC: how well the observations fit, how much was thrown out
    print(
        f"  QC        : {r['n_epochs']} epochs, {r['n_dd']} double differences, {r['n_amb']} ambiguities, "
        f"sigma0 {r['sigma0']:.2f}"
    )
    print(
        f"              RMS phase {r['rms_phase'] * 1000:.1f} mm, code {r['rms_code']:.2f} m; "
        f"rejected {r['n_rejected']} phase / {r['n_rejected_code']} code"
    )
    used = ", ".join(f"{s} {len(v)}" for s, v in r["sats"].items())
    print(f"              satellites used: {used}")
    worst = sorted(r["sat_report"].items(), key=lambda kv: -(kv[1].get("rms_mm") or 0))[:3]
    for sat, rep in worst:
        print(f"              {sat}: rms {rep.get('rms_mm')} mm, arcs {rep.get('arcs')}")
    for n in r["notes"]:
        print(f"  note: {n}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", nargs="?")
    ap.add_argument("rover", nargs="?")
    ap.add_argument("--base-xyz", type=xyz_arg)
    ap.add_argument("--nav", nargs="*", default=[])
    ap.add_argument("--sp3", nargs="*", default=[])
    ap.add_argument("--clk", nargs="*", default=[])
    ap.add_argument("--antex")
    ap.add_argument("--options", default="{}", help="engine options as JSON, e.g. '{\"elmask\": 15}'")
    a = ap.parse_args()

    truth = None
    if a.base is None:
        d = synthetic("static")
        a.base, a.rover = d["stations"]["BASE"], d["stations"]["PT01"]
        a.base_xyz, a.sp3, truth = d["base_xyz"], [d["sp3"]], d["truth"]["PT01"]
    elif a.rover is None or a.base_xyz is None:
        ap.error("give BASE, ROVER and --base-xyz")

    sources, fcn = orbit_sources(a.nav, a.sp3, a.clk)
    r = solve(a.base, a.rover, a.base_xyz, sources, fcn, a.antex, a.options)
    report(r)
    if truth:  # synthetic data: compare with the simulated position
        err = [(r["rover_xyz"][i] - truth[i]) * 1000 for i in range(3)]
        print(f"  vs truth  : {err[0]:+.1f} {err[1]:+.1f} {err[2]:+.1f} mm (dX dY dZ)")


if __name__ == "__main__":
    main()
