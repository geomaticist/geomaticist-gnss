"""Several baselines in one session: a loop closure and a small least-squares network adjustment.

    python baseline_network.py

The engine solves single baselines; the network part (loop check, adjustment) is plain Python here, so it
is easy to follow and to change. The synthetic "static" set has three stations observed together:
BASE (known) and the new points PT01 and PT02. The script

1. solves BASE->PT01, BASE->PT02 and PT01->PT02 (the last from PT01's coordinates of the first solve),
2. closes the triangle (vector sum, mm and ppm of the loop length),
3. adjusts PT01 and PT02 with BASE held fixed, each vector weighted by its 3x3 covariance.

Use it on your own session by editing STATIONS / KNOWN below (or turning them into arguments).
"""

import math

from _common import orbit_sources, synthetic
from static_baseline import solve

d = synthetic("static")
STATIONS = d["stations"]  # name -> RINEX path
KNOWN = {"BASE": d["base_xyz"]}  # held fixed, ECEF metres
SP3 = [d["sp3"]]
TRUTH = d["truth"]  # synthetic data only


# ---- tiny linear algebra (no numpy needed) -------------------------------------------------------------
def solve_linear(a, b):
    """Solves a x = b (Gauss-Jordan with partial pivoting); a is n x n, b has n rows of m columns."""
    n = len(a)
    m = [list(a[i]) + list(b[i]) for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        piv = m[c][c]
        m[c] = [x / piv for x in m[c]]
        for r in range(n):
            if r != c and m[r][c] != 0.0:
                f = m[r][c]
                m[r] = [x - f * y for x, y in zip(m[r], m[c])]
    return [row[n:] for row in m]


def inverse(a):
    n = len(a)
    return solve_linear(a, [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])


def main():
    sources, fcn = orbit_sources(sp3=SP3)

    # 1. the baselines: (from, to) -> result
    plan = [("BASE", "PT01"), ("BASE", "PT02"), ("PT01", "PT02")]
    coords = dict(KNOWN)
    res = {}
    for a, b in plan:
        r = solve(STATIONS[a], STATIONS[b], coords[a], sources, fcn)
        res[(a, b)] = r
        coords.setdefault(b, r["rover_xyz"])  # first estimate of a new point
        sig = [math.sqrt(r["cov_enu"][i][i]) * 1000 for i in range(3)]
        length = math.sqrt(sum(x * x for x in r["vector"]))
        print(
            f"{a}->{b}: {length:10.4f} m  {'fixed' if r['fixed'] else 'FLOAT'}  "
            f"sigma E/N/U {sig[0]:.1f}/{sig[1]:.1f}/{sig[2]:.1f} mm"
        )

    # 2. loop closure: BASE->PT01 + PT01->PT02 - BASE->PT02 should be zero
    v1, v2, v3 = (res[k]["vector"] for k in plan)
    mis = [v1[i] + v3[i] - v2[i] for i in range(3)]
    perim = sum(math.sqrt(sum(x * x for x in v)) for v in (v1, v2, v3))
    m = math.sqrt(sum(x * x for x in mis))
    print(
        f"\nloop BASE-PT01-PT02: misclosure {mis[0] * 1000:+.1f} {mis[1] * 1000:+.1f} {mis[2] * 1000:+.1f} mm "
        f"(dX dY dZ), {m * 1000:.1f} mm over {perim:.0f} m = {m / perim * 1e6:.2f} ppm"
    )

    # 3. least squares: unknowns PT01, PT02 (6 parameters); each vector is an observation of
    #    to - from with weight inverse(covariance). The normal equations are N dx = u.
    unknown = [s for s in STATIONS if s not in KNOWN]
    idx = {s: 3 * k for k, s in enumerate(unknown)}
    n = 3 * len(unknown)
    nmat = [[0.0] * n for _ in range(n)]
    u = [[0.0] for _ in range(n)]
    for (a, b), r in res.items():
        w = inverse(r["cov"])
        # misclosure of the observation at the current coordinates: observed - computed
        comp = [coords[b][i] - coords[a][i] for i in range(3)]
        lvec = [r["vector"][i] - comp[i] for i in range(3)]
        for s, sign in ((b, 1.0), (a, -1.0)):
            if s not in idx:
                continue
            for i in range(3):
                u[idx[s] + i][0] += sign * sum(w[i][j] * lvec[j] for j in range(3))
                for s2, sign2 in ((b, 1.0), (a, -1.0)):
                    if s2 in idx:
                        for j in range(3):
                            nmat[idx[s] + i][idx[s2] + j] += sign * sign2 * w[i][j]
    dx = [x[0] for x in solve_linear(nmat, u)]
    for s in unknown:
        coords[s] = [coords[s][i] + dx[idx[s] + i] for i in range(3)]
    qxx = inverse(nmat)

    # residuals and the variance factor (degrees of freedom: 3 per vector minus the unknowns)
    vtpv = 0.0
    print("\nresiduals (adjusted - observed, mm):")
    for (a, b), r in res.items():
        v = [(coords[b][i] - coords[a][i]) - r["vector"][i] for i in range(3)]
        w = inverse(r["cov"])
        vtpv += sum(v[i] * w[i][j] * v[j] for i in range(3) for j in range(3))
        print(f"  {a}->{b}: {v[0] * 1000:+.1f} {v[1] * 1000:+.1f} {v[2] * 1000:+.1f}")
    dof = 3 * len(res) - n
    s0 = math.sqrt(vtpv / dof) if dof else float("nan")
    print(f"variance factor sigma0 = {s0:.2f} ({dof} degrees of freedom)")

    print("\nadjusted coordinates (ECEF m, sigma mm):")
    for s in unknown:
        k = idx[s]
        sig = [math.sqrt(qxx[k + i][k + i]) * 1000 for i in range(3)]
        x = coords[s]
        line = f"  {s}: {x[0]:.4f} {x[1]:.4f} {x[2]:.4f}  ({sig[0]:.1f} {sig[1]:.1f} {sig[2]:.1f})"
        if s in TRUTH:
            e = [(x[i] - TRUTH[s][i]) * 1000 for i in range(3)]
            line += f"  vs truth {e[0]:+.1f} {e[1]:+.1f} {e[2]:+.1f} mm"
        print(line)


if __name__ == "__main__":
    main()
