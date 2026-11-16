//! RINEX 3.05 observation and mixed navigation files from a decoded Topcon `.tps` log.
//!
//! The observation file is laid out as office software expects for these logs:
//!
//! - observation types per system in a fixed order (each signal as C, L, D, S), only the signals present;
//! - satellites of an epoch in the receiver's order; the receiver clock offset on the epoch line;
//! - the phases with these phase shifts taken off (`SYS / PHASE SHIFT`: GPS L1W +1/4, L2X -1/4, GLONASS
//!   L1P and L2P +1/4, BeiDou L7I 0.16667 cycle);
//! - signal strength truncated to three decimals (it is a sixteenth of a dB-Hz).
//!
//! It also writes loss-of-lock indicators (bit 0 only, from the lock time and the block flags, rules at
//! [`Lli::of`]), the receiver firmware, `SIGNAL STRENGTH UNIT`, and BeiDou PRN 38 and up.

use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;

use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::time::{header_time, v3_epoch_time};

use super::nav::{Eph, GpsUtc};
use super::obs::{Epoch, SatObs, Signal, carrier};
use crate::t04::nav::{fmt_d19, fmt_e};
use crate::util::hline;

/// RINEX header order of systems.
const SYSTEMS: [GnssSystem; 6] = [
    GnssSystem::Gps,
    GnssSystem::Glonass,
    GnssSystem::Galileo,
    GnssSystem::BeiDou,
    GnssSystem::Qzss,
    GnssSystem::Sbas,
];

fn sys_rank(s: GnssSystem) -> usize {
    SYSTEMS.iter().position(|&x| x == s).unwrap_or(9)
}

/// The order of the signals of a system in `SYS / # / OBS TYPES`.
fn header_order(sys: GnssSystem) -> &'static [[u8; 2]] {
    const G: [[u8; 2]; 5] = [*b"1C", *b"1W", *b"2X", *b"2W", *b"5X"];
    const R: [[u8; 2]; 5] = [*b"1C", *b"1P", *b"2C", *b"2P", *b"3X"];
    const E: [[u8; 2]; 4] = [*b"1X", *b"8X", *b"5X", *b"7X"];
    const C: [[u8; 2]; 4] = [*b"2I", *b"5P", *b"7I", *b"7D"];
    match sys {
        GnssSystem::Gps => &G,
        GnssSystem::Glonass => &R,
        GnssSystem::Galileo => &E,
        GnssSystem::BeiDou => &C,
        _ => &[],
    }
}

/// The phase shift of a signal, cycles (`None`: the system's reference signal, written without a value).
#[must_use]
pub fn phase_shift(sys: GnssSystem, sig: [u8; 2]) -> Option<f64> {
    match (sys, &sig) {
        (GnssSystem::Gps, b"1C")
        | (GnssSystem::Glonass, b"1C" | b"2C")
        | (GnssSystem::BeiDou, b"2I") => None,
        (GnssSystem::Gps, b"1W") | (GnssSystem::Glonass, b"1P" | b"2P") => Some(0.25),
        (GnssSystem::Gps, b"2X") => Some(-0.25),
        (GnssSystem::BeiDou, b"7I") => Some(0.16667),
        _ => Some(0.0),
    }
}

/// A satellite as RINEX names it, ordered for the header.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct SatKey(usize, u8, char);

fn sat_key(s: &SatObs) -> SatKey {
    SatKey(sys_rank(s.sys), s.prn, s.sys.letter())
}

fn sat_id(s: &SatObs) -> String {
    format!("{}{:02}", s.sys.letter(), s.prn)
}

/// What the header needs, gathered over all epochs (pass 1).
#[derive(Debug, Clone, Default)]
pub struct Stats {
    /// Signals seen per system, in first-seen order.
    pub signals: BTreeMap<usize, Vec<[u8; 2]>>,
    /// Values per satellite and signal: (C, L, D, S) counts.
    counts: BTreeMap<SatKey, BTreeMap<[u8; 2], [usize; 4]>>,
    /// First and last epoch, ticks.
    pub first: Option<i64>,
    /// Last epoch, ticks.
    pub last: Option<i64>,
    /// Epoch spacings (ticks) and how often each was seen.
    pub spacings: BTreeMap<i64, usize>,
    /// GLONASS slot -> frequency channel.
    pub glonass: BTreeMap<u8, i8>,
    /// Epochs.
    pub epochs: usize,
    systems: BTreeMap<usize, GnssSystem>,
}

impl Stats {
    /// Adds an epoch at `ticks`.
    pub fn add(&mut self, e: &Epoch, ticks: i64) {
        if let Some(l) = self.last {
            *self.spacings.entry(ticks - l).or_default() += 1;
        }
        self.first.get_or_insert(ticks);
        self.last = Some(ticks);
        self.epochs += 1;
        for s in &e.sats {
            let r = sys_rank(s.sys);
            self.systems.insert(r, s.sys);
            let list = self.signals.entry(r).or_default();
            let c = self.counts.entry(sat_key(s)).or_default();
            for g in &s.signals {
                if !list.contains(&g.sig) {
                    list.push(g.sig);
                }
                let n = c.entry(g.sig).or_default();
                for v in n.iter_mut() {
                    *v += 1;
                }
                if carrier(s.sys, g.sig[0] - b'0', s.fcn).is_none() {
                    n[1] -= 1;
                }
            }
            if let (GnssSystem::Glonass, Some(f)) = (s.sys, s.fcn) {
                self.glonass.insert(s.prn, f);
            }
        }
    }

    /// The signals of a system in header order.
    #[must_use]
    pub fn header_signals(&self, sys: GnssSystem) -> Vec<[u8; 2]> {
        let seen = self
            .signals
            .get(&sys_rank(sys))
            .cloned()
            .unwrap_or_default();
        let mut out: Vec<[u8; 2]> = header_order(sys)
            .iter()
            .copied()
            .filter(|s| seen.contains(s))
            .collect();
        let extra: Vec<[u8; 2]> = seen.iter().copied().filter(|s| !out.contains(s)).collect();
        out.extend(extra);
        out
    }

    /// The most common epoch spacing, s.
    #[must_use]
    pub fn interval(&self) -> Option<f64> {
        self.spacings
            .iter()
            .filter(|(d, _)| **d > 0)
            .max_by_key(|(d, n)| (**n, -**d))
            .map(|(d, _)| *d as f64 / 1e7)
    }
}

/// Header facts that don't come from the epochs.
#[derive(Debug, Clone, Default)]
pub struct ObsHeader {
    /// PGM.
    pub program: String,
    /// RUN BY.
    pub run_by: String,
    /// The `PGM / RUN BY / DATE` date text.
    pub date: String,
    /// COMMENT lines.
    pub comments: Vec<String>,
    /// MARKER NAME.
    pub marker: String,
    /// Receiver number (the receiver ID).
    pub receiver_number: String,
    /// Receiver type.
    pub receiver_type: String,
    /// Receiver firmware.
    pub receiver_version: String,
    /// Antenna type with radome (20 characters).
    pub antenna: String,
    /// APPROX POSITION XYZ.
    pub approx: [f64; 3],
    /// ANTENNA: DELTA H/E/N.
    pub delta_hen: [f64; 3],
    /// LEAP SECONDS.
    pub leap: Option<i64>,
}

fn fmt_milli(m: i64, width: usize) -> String {
    let sign = if m < 0 { "-" } else { "" };
    let a = m.unsigned_abs();
    format!("{:>width$}", format!("{sign}{}.{:03}", a / 1000, a % 1000))
}

/// The header text (LF line ends).
#[must_use]
pub fn header_text(st: &Stats, h: &ObsHeader) -> String {
    let mut o = String::new();
    let mut line = |c: &str, l: &str| {
        o.push_str(hline(c, l).trim_end());
        o.push('\n');
    };
    line(
        "     3.05           OBSERVATION DATA    M",
        "RINEX VERSION / TYPE",
    );
    line(
        &format!("{:<20.20}{:<20.20}{:<20.20}", h.program, h.run_by, h.date),
        "PGM / RUN BY / DATE",
    );
    for c in &h.comments {
        line(c, "COMMENT");
    }
    line(&h.marker, "MARKER NAME");
    line("", "MARKER TYPE");
    line("", "OBSERVER / AGENCY");
    line(
        &format!(
            "{:<20.20}{:<20.20}{:<20.20}",
            h.receiver_number, h.receiver_type, h.receiver_version
        ),
        "REC # / TYPE / VERS",
    );
    line(&format!("{:<20}{:<20.20}", "", h.antenna), "ANT # / TYPE");
    line(
        &format!(
            "{:14.4}{:14.4}{:14.4}",
            h.approx[0], h.approx[1], h.approx[2]
        ),
        "APPROX POSITION XYZ",
    );
    line(
        &format!(
            "{:14.4}{:14.4}{:14.4}",
            h.delta_hen[0], h.delta_hen[1], h.delta_hen[2]
        ),
        "ANTENNA: DELTA H/E/N",
    );
    let systems: Vec<GnssSystem> = st.systems.values().copied().collect();
    for &sys in &systems {
        let sigs = st.header_signals(sys);
        let types: Vec<String> = sigs
            .iter()
            .flat_map(|s| {
                let t = String::from_utf8_lossy(s).into_owned();
                ["C", "L", "D", "S"].map(|k| format!("{k}{t}"))
            })
            .collect();
        for (k, chunk) in types.chunks(13).enumerate() {
            let mut c = if k == 0 {
                format!("{}  {:3}", sys.letter(), types.len())
            } else {
                "      ".to_owned()
            };
            for t in chunk {
                let _ = write!(c, " {t}");
            }
            line(&c, "SYS / # / OBS TYPES");
        }
    }
    line("DBHZ", "SIGNAL STRENGTH UNIT");
    if let Some(iv) = st.interval() {
        line(&format!("{iv:10.3}"), "INTERVAL");
    }
    if let (Some(a), Some(b)) = (st.first, st.last) {
        line(&header_time(a, "GPS"), "TIME OF FIRST OBS");
        line(&header_time(b, "GPS"), "TIME OF LAST OBS");
    }
    if let Some(l) = h.leap {
        line(&format!("{l:6}"), "LEAP SECONDS");
    }
    line(&format!("{:6}", 0), "RCV CLOCK OFFS APPL");
    for &sys in &systems {
        for s in st.header_signals(sys) {
            let t = String::from_utf8_lossy(&s).into_owned();
            let c = match phase_shift(sys, s) {
                Some(v) => format!("{} L{t} {v:8.5}", sys.letter()),
                None => format!("{} L{t}", sys.letter()),
            };
            line(&c, "SYS / PHASE SHIFT");
        }
    }
    if !st.glonass.is_empty() {
        let all: Vec<(&u8, &i8)> = st.glonass.iter().collect();
        for (k, chunk) in all.chunks(8).enumerate() {
            let mut c = if k == 0 {
                format!("{:3} ", all.len())
            } else {
                "    ".to_owned()
            };
            for (slot, f) in chunk {
                let _ = write!(c, "R{slot:02} {f:2} ");
            }
            line(&c, "GLONASS SLOT / FRQ #");
        }
        line(
            " C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000",
            "GLONASS COD/PHS/BIS",
        );
    }
    line(&format!("{:6}", st.counts.len()), "# OF SATELLITES");
    for (key, sigs) in &st.counts {
        let sys = st.systems.get(&key.0).copied().unwrap_or(GnssSystem::Gps);
        let mut n = Vec::new();
        for s in st.header_signals(sys) {
            let c = sigs.get(&s).copied().unwrap_or_default();
            n.extend_from_slice(&c);
        }
        for (k, chunk) in n.chunks(9).enumerate() {
            let mut c = if k == 0 {
                format!("   {}{:02}", key.2, key.1)
            } else {
                "      ".to_owned()
            };
            for v in chunk {
                let _ = write!(c, "{v:6}");
            }
            line(&c, "PRN / # OF OBS");
        }
    }
    line("", "END OF HEADER");
    o
}

/// The loss-of-lock indicator state of one signal.
#[derive(Debug, Clone, Copy)]
struct Track {
    ticks: i64,
    settled: bool,
    half: bool,
}

/// Loss-of-lock indicators over a log (pass 2).
#[derive(Debug, Clone, Default)]
pub struct Lli {
    tracks: HashMap<(GnssSystem, u8, [u8; 2]), Track>,
}

impl Lli {
    /// The LLI of a signal at `ticks`: bit 0 (a new lock, or a possible slip) when the lock time is shorter than
    /// the time since the signal was last seen, while the block's bit 72 is clear and at the epoch it sets, or
    /// when its bit 73 changes (at the start of a track that change came with a jump of half a cycle in most
    /// cases). The first epoch of a signal gets no flag. Bit 1 is never set: bit 73 stays the
    /// same over whole tracks, so it can't be the half-cycle status RINEX means.
    pub fn of(&mut self, s: &SatObs, g: &Signal, ticks: i64) -> u8 {
        let key = (s.sys, s.prn, g.sig);
        let prev = self.tracks.insert(
            key,
            Track {
                ticks,
                settled: g.settled,
                half: g.half_cycle_known,
            },
        );
        let Some(p) = prev else {
            return 0;
        };
        let since = (ticks - p.ticks) as f64 / 1e7;
        let lock = f64::from(g.lock) * 0.2;
        u8::from(lock + 1e-6 < since || !g.settled || !p.settled || p.half != g.half_cycle_known)
    }
}

/// Writes the observation records of one epoch.
pub fn epoch_text(out: &mut String, st: &Stats, e: &Epoch, ticks: i64, clock: bool, lli: &mut Lli) {
    let n = e.sats.len();
    if clock {
        let _ = writeln!(
            out,
            "> {}  0{n:3}      {:15.12}",
            v3_epoch_time(ticks),
            e.clock_s()
        );
    } else {
        let _ = writeln!(out, "> {}  0{n:3}", v3_epoch_time(ticks));
    }
    for s in &e.sats {
        let mut rec = sat_id(s);
        let mut by_sig: HashMap<[u8; 2], (&Signal, u8)> = HashMap::new();
        for g in &s.signals {
            let l = lli.of(s, g, ticks);
            by_sig.insert(g.sig, (g, l));
        }
        for sig in st.header_signals(s.sys) {
            match by_sig.get(&sig) {
                Some((g, l)) => {
                    rec.push_str(&fmt_milli(g.range * 20, 14));
                    rec.push_str("  ");
                    match carrier(s.sys, sig[0] - b'0', s.fcn) {
                        Some(f) => {
                            let v = g.phase_cycles(f) - phase_shift(s.sys, sig).unwrap_or(0.0);
                            let _ = write!(rec, "{v:14.3}");
                            if *l > 0 {
                                let _ = write!(rec, "{l} ");
                            } else {
                                rec.push_str("  ");
                            }
                        }
                        None => rec.push_str("                "),
                    }
                    rec.push_str(&fmt_milli(i64::from(g.doppler), 14));
                    rec.push_str("  ");
                    rec.push_str(&fmt_milli(i64::from(g.snr16) * 125 / 2, 14));
                    rec.push_str("  ");
                }
                None => rec.push_str(&" ".repeat(64)),
            }
        }
        out.push_str(rec.trim_end());
        out.push('\n');
    }
}

/// Navigation header facts.
#[derive(Debug, Clone, Default)]
pub struct NavHeader {
    /// PGM.
    pub program: String,
    /// RUN BY.
    pub run_by: String,
    /// Date text.
    pub date: String,
    /// GPS ionosphere alpha, beta.
    pub gps_iono: Option<([f64; 4], [f64; 4])>,
    /// BeiDou ionosphere alpha, beta.
    pub bds_iono: Option<([f64; 4], [f64; 4])>,
    /// GPS UTC parameters.
    pub gps_utc: Option<GpsUtc>,
    /// GLONASS - UTC: tau_c (s), and the GPS time of week and week it was logged at.
    pub glut: Option<(f64, i64, i64)>,
    /// GPS - UTC, s.
    pub leap: Option<i64>,
}

fn nav_line_fields(o: &mut String, vals: &[Option<f64>]) {
    let mut s = "    ".to_owned();
    for v in vals {
        s.push_str(&fmt_d19(*v));
    }
    o.push_str(s.trim_end());
    o.push('\n');
}

/// The mixed RINEX 3.05 navigation file: records grouped by system (G R E C J S), each group by epoch then
/// satellite. LF line ends.
#[must_use]
pub fn nav_text(recs: &[Eph], h: &NavHeader) -> String {
    let mut o = String::new();
    let mut line = |c: &str, l: &str| {
        o.push_str(hline(c, l).trim_end());
        o.push('\n');
    };
    line(
        "     3.05           N: GNSS NAV DATA    M: MIXED",
        "RINEX VERSION / TYPE",
    );
    line(
        &format!("{:<20.20}{:<20.20}{:<20.20}", h.program, h.run_by, h.date),
        "PGM / RUN BY / DATE",
    );
    let iono = |name: &str, v: &[f64; 4]| {
        let mut c = format!("{name} ");
        for x in v {
            c.push_str(&fmt_e(*x, 4, 12));
        }
        c
    };
    if let Some((a, b)) = &h.gps_iono {
        line(&iono("GPSA", a), "IONOSPHERIC CORR");
        line(&iono("GPSB", b), "IONOSPHERIC CORR");
    }
    if let Some((a, b)) = &h.bds_iono {
        line(&iono("BDSA", a), "IONOSPHERIC CORR");
        line(&iono("BDSB", b), "IONOSPHERIC CORR");
    }
    if let Some(u) = &h.gps_utc {
        line(
            &format!(
                "GPUT {}{} {:6} {:4}",
                fmt_e(u.a0, 10, 17),
                fmt_e(u.a1, 9, 16),
                u.tot,
                u.wnt
            ),
            "TIME SYSTEM CORR",
        );
    }
    if let Some((tau_c, tow, week)) = h.glut {
        line(
            &format!(
                "GLUT {}{} {tow:6} {week:4}",
                fmt_e(-tau_c, 10, 17),
                fmt_e(0.0, 9, 16)
            ),
            "TIME SYSTEM CORR",
        );
    }
    if let Some(u) = &h.gps_utc {
        line(
            &format!("{:6}{:6}{:6}{:6}", u.dtls, u.dtlsf, u.wnlsf, u.dn),
            "LEAP SECONDS",
        );
    } else if let Some(l) = h.leap {
        line(&format!("{l:6}"), "LEAP SECONDS");
    }
    line("", "END OF HEADER");
    let order = "GREJCS";
    let mut v: Vec<&Eph> = recs.iter().collect();
    v.sort_by_key(|e| {
        (
            order.find(e.rec.sys).unwrap_or(9),
            e.rec.epoch.micros_since_gps0(),
            e.rec.sat.clone(),
            e.iod,
        )
    });
    for e in v {
        let r = &e.rec;
        let t = &r.epoch;
        let mut first = format!(
            "{} {:04} {:02} {:02} {:02} {:02} {:02}",
            r.sat, t.year, t.month, t.day, t.hour, t.minute, t.second
        );
        for x in r.vals.iter().take(3) {
            first.push_str(&fmt_d19(*x));
        }
        o.push_str(first.trim_end());
        o.push('\n');
        for chunk in r.vals[3..].chunks(4) {
            nav_line_fields(&mut o, chunk);
        }
    }
    o
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn milli_format() {
        assert_eq!(fmt_milli(455_317_624_800, 14), " 455317624.800");
        assert_eq!(fmt_milli(-1514, 14), "        -1.514");
        assert_eq!(fmt_milli(-999, 6), "-0.999");
        assert_eq!(fmt_milli(0, 6), " 0.000");
        // 49 + 3/16 dB-Hz prints truncated
        assert_eq!(
            fmt_milli(i64::from(49u16 * 16 + 3) * 125 / 2, 8),
            "  49.187"
        );
    }

    #[test]
    fn shifts() {
        assert_eq!(phase_shift(GnssSystem::Gps, *b"1C"), None);
        assert_eq!(phase_shift(GnssSystem::Gps, *b"2X"), Some(-0.25));
        assert_eq!(phase_shift(GnssSystem::Galileo, *b"1X"), Some(0.0));
    }
}
