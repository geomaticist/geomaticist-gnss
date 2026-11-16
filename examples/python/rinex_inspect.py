"""Look inside a RINEX observation file: header, epochs and gaps, satellites, event records, time check.

    python rinex_inspect.py [OBS]

OBS is a RINEX 2 / 3 / 4 observation file (plain, .gz or Hatanaka .crx/.d). Without it the script
inspects the synthetic stop-and-go rover, which carries flag 3 / flag 2 records.
"""

import json
import sys

import geomaticist_gnss_py as g

from _common import synthetic


def main(path):
    # rinex_info: what `gnss rinex info --json` prints (the same summary the app's Overview page shows)
    info = json.loads(g.rinex_info(path))
    h = info["header"]
    print(f"{info['file']}  (RINEX {info['rinex_version']})")
    print(f"  marker   : {h['marker_name']!r}  type {h['marker_type'] or '-'}")
    print(f"  receiver : {h.get('receiver') or '-'}   antenna: {h['antenna_type'] or '-'} {h['antenna_radome']}")
    print(f"  ARP      : delta H/E/N {h['delta_hen']} m, approx XYZ {h['approx_xyz']}")
    print(f"  span     : {info['first_epoch']} -> {info['last_epoch']}  ({info['duration']}), GPS time")
    print(f"  interval : {info['interval_actual_s']} s (header says {info['interval_header_s']})")
    print(
        f"  epochs   : {info['epochs']} of {info['epochs_expected']} expected "
        f"({info['completeness_pct']} % complete), {info['duplicate_epochs']} duplicate, "
        f"{info['out_of_order_epochs']} out of order"
    )
    print(f"  gaps     : {len(info['gaps'])}")
    for gap in info["gaps"][:10]:
        print(f"             {gap}")
    print(f"  sats     : {info['satellites']} (mean {info['mean_sats_per_epoch']} per epoch)")
    for sys_letter, s in info["systems"].items():
        snr = ", ".join(f"{k} {v} dB-Hz" for k, v in s["snr"].items())
        print(f"    {s['name']:<8} {s['count']:3d} sats, {s['slips']} receiver-flagged slips; mean SNR {snr}")
        print(f"             obs types: {' '.join(h['obs_types'].get(sys_letter, []))}")
    print(f"  special records (flags 2-5): {info['events']}")

    # rinex_events: external events (flag 5, e.g. camera shutter marks) as (GPS seconds, label, offset)
    ev = json.loads(g.rinex_events(path))
    print(f"  flag-5 events: {len(ev['events'])}")
    for t, label, _ in ev["events"][:5]:
        print(f"    {t:.6f}  {label}")
    for w in ev["warnings"]:
        print(f"    warning: {w}")

    # rinex_occupations: stop-and-go occupations from flag 3 (new site) / flag 2 (start moving) records
    occ = json.loads(g.rinex_occupations(path))
    print(f"  occupations: {len(occ)}")
    for name, t0, t1 in occ:
        print(f"    {name:<8} {t1 - t0:6.0f} s from GPS second {t0:.1f}")

    # check_times: epochs before 1990 (GPS week rollover), displaced or non-increasing epochs
    tc = json.loads(g.check_times(path))
    print(
        f"  time check: {tc['epochs']} epochs, {tc['before_1990']} before 1990, "
        f"{tc['displaced']} displaced, {tc['not_later']} not later than the one before"
    )
    for w in tc["warnings"]:
        print(f"    warning: {w}")

    # NGS OPUS's published requirements, as a quick static-session QC
    print("  OPUS checks:")
    for c in info["opus"]:
        print(f"    [{c['level']:>4}] {c['text']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else synthetic("stop_and_go")["rover"])
