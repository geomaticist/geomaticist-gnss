"""Prepare a drone rover RINEX for office post-processing software (e.g. Trimble Business Center).

    python rinex_for_tbc.py                                   # synthetic flight + MRK
    python rinex_for_tbc.py ROVER.obs Timestamp.MRK [--out ROVER_events.obs]

Office software reads photo events only from the RINEX itself, and treats a file as kinematic only when
the header says so. Two steps:

1. ``insert_events`` writes one flag-5 event record per MRK photo, in time order between the epochs.
2. ``mark_kinematic`` sets ``MARKER TYPE`` (default AIRBORNE) and adds a flag-2 "start moving" record
   before the first epoch, so the whole file is processed as a moving receiver. It also drops blank
   header lines that some readers reject.

The script then reads the result back to check both changes.
"""

import argparse
import json
import os

import geomaticist_gnss_py as g

from _common import epochs, ensure_out, synthetic


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rover", nargs="?")
    ap.add_argument("mrk", nargs="?")
    ap.add_argument("--out")
    ap.add_argument("--marker-type", default="AIRBORNE")
    a = ap.parse_args()
    if a.rover is None:
        d = synthetic("flight")
        a.rover, a.mrk = d["rover"], d["mrk"]
    elif a.mrk is None:
        ap.error("give ROVER and MRK")

    out = a.out or os.path.join(ensure_out(), "rover_events.obs")
    tmp = out + ".step1"

    # 1. events from the MRK (week resolved from the rover's first epoch), written as flag-5 records
    first = epochs(g.read_obs(a.rover, kinds="C"))[0]
    events = json.loads(g.read_mrk(a.mrk, first))["events"]
    r1 = json.loads(g.insert_events(a.rover, events, tmp))
    print(f"insert_events : {json.dumps(r1)}")

    # 2. kinematic header + flag 2
    r2 = json.loads(g.mark_kinematic(tmp, out, a.marker_type))
    os.remove(tmp)
    print(f"mark_kinematic: flag-2 added {r2['flag_added']}, blank header lines removed {r2['blank_removed']}")

    # 3. check
    info = json.loads(g.rinex_info(out))
    back = json.loads(g.rinex_events(out))["events"]
    print(f"result        : {out}")
    print(f"  MARKER TYPE {info['header']['marker_type']!r}, special records {info['events']}")
    print(f"  {len(back)} flag-5 events read back (of {len(events)} in the MRK)")


if __name__ == "__main__":
    main()
