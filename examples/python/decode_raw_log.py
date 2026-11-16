"""Raw receiver log -> RINEX: Trimble .T04/.T02/.T01/.T00, Leica .MDB, Topcon .tps / Javad .jps (GREIS),
RTCM 3, u-blox UBX and NovAtel OEM binary.

    python decode_raw_log.py LOGFILE [OUT_DIR]

``raw_kind`` names the format from the file's content (``None`` when it is not one the decoders read);
``decode_raw`` writes the observation file (RINEX 3) and any navigation files into OUT_DIR and returns
``(obs_path, [nav_paths], [warnings], [log_lines])``.

The repository ships no raw logs, so this example has no built-in data: give it a log of your own.
"""

import os
import sys

import geomaticist_gnss_py as g

from _common import ensure_out


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else ensure_out("decoded")
    kind = g.raw_kind(path)
    print(f"{os.path.basename(path)}: {kind or 'not a raw log the decoders know'}")
    if kind is None:
        return
    obs, navs, warnings, log = g.decode_raw(path, out_dir)
    for line in log:
        print(f"  {line}")
    for w in warnings:
        print(f"  warning: {w}")
    print(f"observations -> {obs}")
    for n in navs:
        print(f"navigation   -> {n}")


if __name__ == "__main__":
    main()
