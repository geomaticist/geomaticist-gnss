"""Small helpers shared by the examples (not part of the engine).

The engine's Python module is ``geomaticist_gnss_py``. Its functions take plain values (paths, lists,
dicts) and most return JSON text, which the examples turn into Python objects with ``json.loads``.
"""

import datetime as _dt
import json
import os
import struct

import geomaticist_gnss_py as g

HERE = os.path.dirname(os.path.abspath(__file__))
# Where the examples write the synthetic data and their outputs (git-ignored).
OUT = os.path.join(HERE, "out")

GPS_EPOCH = _dt.datetime(1980, 1, 6)
WEEK = 604800.0


def station(h):
    """The station dict the solvers take, from a ``Station`` handle returned by ``read_obs``.

    The dict lets you override what the RINEX header says: the station name, the antenna type (as it is
    written in ANTEX, e.g. ``"TRM59800.00     NONE"``), the receiver, the approximate position and the
    antenna height / east / north offsets (``delta_hen``, metres, ARP above the mark).
    """
    return {
        "handle": h,
        "name": h.name,
        "antenna": h.antenna,
        "receiver": h.receiver,
        "approx": h.approx,
        "delta_hen": h.delta_hen,
    }


def epochs(h):
    """The epochs of a ``Station`` handle as a list of GPS seconds."""
    b = h.t_bytes()
    return list(struct.unpack(f"<{len(b) // 8}d", b))


def gpst(t, decimals=1):
    """GPS seconds -> 'YYYY-MM-DD hh:mm:ss.s' in GPS time (GPST runs 18 s ahead of UTC since 2017)."""
    s = (GPS_EPOCH + _dt.timedelta(seconds=t)).strftime("%Y-%m-%d %H:%M:%S.%f")
    return s[: len(s) - 6 + decimals] if decimals else s[:-7]


def week_tow(t):
    """GPS seconds -> (GPS week, seconds of week)."""
    w = int(t // WEEK)
    return w, t - w * WEEK


def synthetic(kind):
    """Writes (once) and returns the synthetic data set ``kind``: "static", "flight" or "stop_and_go".

    The data are simulated by the engine (``synth_session``); nothing comes from a real receiver. The
    returned dict holds the file paths, ``base_xyz`` (the base mark, ECEF metres) and ``truth``.
    """
    folder = os.path.join(OUT, "data", kind)
    return json.loads(g.synth_session(folder, kind))


def orbit_sources(nav=None, sp3=None, clk=None):
    """The orbit sources list for ``static_solve`` / ``ppk_solve``, and the GLONASS channels.

    Broadcast navigation files (``nav``) and precise SP3 orbits (``sp3``, with optional RINEX clocks
    ``clk``) can be given together; the solver uses them in the order given.
    """
    sources, fcn = [], {}
    if sp3:
        sources.append(("precise_handle", g.read_sp3(list(sp3), list(clk or []))))
    if nav:
        n = g.read_nav(list(nav))
        for name, err in n.bad:
            print(f"  nav file not read: {name}: {err}")
        sources.append(("broadcast_handle", n))
        fcn = n.glo_fcn()
    if not sources:
        raise SystemExit("give navigation files (--nav) or precise orbits (--sp3)")
    return sources, fcn


def xyz_arg(text):
    """'X,Y,Z' (ECEF metres) -> [x, y, z]."""
    v = [float(x) for x in text.split(",")]
    if len(v) != 3:
        raise ValueError("expected X,Y,Z")
    return v


def ensure_out(*parts):
    """Creates and returns a folder under ``out/``."""
    p = os.path.join(OUT, *parts)
    os.makedirs(p, exist_ok=True)
    return p
