"""Smoke test: runs every example on its synthetic data and reports which ones failed.

    python run_all.py

Each script runs in its own Python process from this folder; a non-zero exit is a failure.
"""

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    "rinex_inspect.py",
    "precise_orbits.py",
    "static_baseline.py",
    "baseline_network.py",
    "ppk_events_to_points.py",
    "stop_and_go.py",
    "ppp.py",
    "rinex_for_tbc.py",
    "decode_raw_log.py",  # no data: prints its usage
]


def main():
    failed = []
    for s in SCRIPTS:
        t = time.time()
        p = subprocess.run([sys.executable, s], cwd=HERE, capture_output=True, text=True)
        ok = p.returncode == 0
        print(f"{'ok  ' if ok else 'FAIL'} {s:<26} {time.time() - t:6.1f} s")
        if not ok:
            failed.append(s)
            print(p.stdout[-2000:], p.stderr[-2000:], sep="\n")
    print(f"\n{len(SCRIPTS) - len(failed)} of {len(SCRIPTS)} examples ran")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
