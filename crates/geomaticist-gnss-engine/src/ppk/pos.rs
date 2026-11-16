//! Text outputs: the `.pos` solution format of RTKLIB for epochs and for events, and a CSV for events.
//!
//! The `.pos` layout (header lines starting with `%`, then one line per position with the columns
//! `GPST latitude longitude height Q ns sdn sde sdu sdne sdeu sdun age ratio`) is the one written by RTKLIB
//! (T. Takasu, BSD-2-Clause; `solution.c`), reproduced here from its published format so that programs and
//! scripts made for RTKLIB output read these files unchanged. `Q` is 1 for a fixed, 2 for a float and 5 for
//! a code-only (single) position; epochs without a position are not written.
//!
//! These functions return text; the caller writes the file.

use std::fmt::Write as _;

use geomaticist_gnss_core::time::{TICKS, from_ticks, gps_week_tow};

use super::events::{EventResult, EventStatus};
use super::options::{ArMode, Direction, PpkOptions};
use super::{EpochSolution, PpkResult, llh_deg};

/// What the `.pos` header says about the run.
#[derive(Debug, Clone, PartialEq)]
pub struct PosHeader {
    /// Program name and version for the `% program` line.
    pub program: String,
    /// Input files for the `% inp file` lines (rover, base, navigation ...).
    pub inputs: Vec<String>,
    /// Decimals of the seconds in the time column (RTKLIB's default is 3).
    pub time_decimals: usize,
    /// Add the velocity columns `vn ve vu` (m/s) after `ratio`.
    pub velocity: bool,
    /// The options of the run.
    pub options: PpkOptions,
    /// The `% pos mode` line (`Kinematic` for a differential solution).
    pub mode: String,
}

impl PosHeader {
    /// A header for a run with `options`.
    #[must_use]
    pub fn new(options: &PpkOptions) -> Self {
        Self {
            program: concat!("geomaticist-gnss ", env!("CARGO_PKG_VERSION")).to_owned(),
            inputs: Vec::new(),
            time_decimals: 3,
            velocity: false,
            options: options.clone(),
            mode: "Kinematic".to_owned(),
        }
    }
}

/// `2024/02/08 00:00:49.000` of GPS seconds, rounded to `decimals` (at most 7).
#[must_use]
pub fn gpst_text(t: f64, decimals: usize) -> String {
    let dec = decimals.min(7);
    let scale = 10f64.powi(dec as i32);
    let n = (t * scale).round();
    let whole = (n / scale).floor();
    let frac = (n - whole * scale).max(0.0) as u64;
    let c = from_ticks(whole as i64 * TICKS);
    let mut s = format!(
        "{:04}/{:02}/{:02} {:02}:{:02}:{:02}",
        c.year,
        c.month,
        c.day,
        c.hour,
        c.minute,
        c.sec_ticks / TICKS
    );
    if dec > 0 {
        let _ = write!(s, ".{frac:0dec$}");
    }
    s
}

/// RTKLIB's `sqvar`: the signed square root of a covariance element.
fn sqvar(c: f64) -> f64 {
    if c < 0.0 { -(-c).sqrt() } else { c.sqrt() }
}

fn systems_text(s: &str) -> String {
    let mut v = Vec::new();
    for (c, name) in [
        ('G', "GPS"),
        ('R', "GLONASS"),
        ('E', "Galileo"),
        ('J', "QZSS"),
        ('C', "BeiDou"),
    ] {
        if s.contains(c) {
            v.push(name);
        }
    }
    v.join(" ")
}

fn header(res: &PpkResult, h: &PosHeader, first: Option<f64>, last: Option<f64>) -> String {
    let o = &h.options;
    let mut s = String::new();
    let _ = writeln!(s, "% program   : {}", h.program);
    for f in &h.inputs {
        let _ = writeln!(s, "% inp file  : {f}");
    }
    for (name, t) in [("obs start", first), ("obs end  ", last)] {
        if let Some(t) = t {
            let (w, tow) = gps_week_tow(t);
            let _ = writeln!(
                s,
                "% {name} : {} GPST (week{w:04} {tow:8.1}s)",
                gpst_text(t, 1)
            );
        }
    }
    let _ = writeln!(s, "% pos mode  : {}", h.mode);
    let _ = writeln!(
        s,
        "% freqs     : {}",
        if res.summary.signals.is_empty() {
            "-".to_owned()
        } else {
            res.summary.signals.join("; ")
        }
    );
    let _ = writeln!(
        s,
        "% solution  : {}",
        match o.direction {
            Direction::Forward => "Forward",
            Direction::Backward => "Backward",
            Direction::Combined => "Combined",
        }
    );
    let _ = writeln!(s, "% elev mask : {:.1} deg", o.elmask);
    let _ = writeln!(s, "% dynamics  : off");
    let _ = writeln!(s, "% tidecorr  : off");
    let _ = writeln!(s, "% ionos opt : off (short baseline)");
    let _ = writeln!(s, "% tropo opt : {}", o.trop_model);
    let _ = writeln!(s, "% ephemeris : as given");
    let _ = writeln!(s, "% navi sys  : {}", systems_text(&o.systems));
    let _ = writeln!(
        s,
        "% amb res   : {}",
        match o.ar {
            ArMode::Off => "Off",
            ArMode::Continuous => "Continuous",
        }
    );
    let _ = writeln!(
        s,
        "% amb glo   : {}",
        if o.ar_glonass { "On" } else { "Off" }
    );
    let _ = writeln!(s, "% val thres : {:.1}", o.ar_ratio);
    let m = llh_deg(res.base_marker);
    if m[0].is_finite() {
        let _ = writeln!(s, "% ref pos   :{:14.9} {:14.9} {:10.4}", m[0], m[1], m[2]);
    }
    for w in &res.warnings {
        let _ = writeln!(s, "% warning   : {w}");
    }
    let _ = writeln!(s, "%");
    let _ = writeln!(
        s,
        "% (lat/lon/height=WGS84/ellipsoidal,Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)"
    );
    let _ = write!(
        s,
        "%  GPST{}latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  \
         sdun(m) age(s)  ratio",
        " ".repeat(15 + h.time_decimals.min(7))
    );
    if h.velocity {
        s.push_str("    vn(m/s)    ve(m/s)    vu(m/s)");
    }
    s.push('\n');
    s
}

#[allow(clippy::too_many_arguments)]
fn line(
    out: &mut String,
    t: f64,
    dec: usize,
    llh: [f64; 3],
    q: u8,
    ns: u16,
    cov: &[[f64; 3]; 3],
    age: f32,
    ratio: f32,
    vel: Option<[f32; 3]>,
) {
    // cov is east / north / up; the columns are sdn sde sdu sdne sdeu sdun
    let _ = write!(
        out,
        "{} {:14.9} {:14.9} {:10.4} {:3} {:3} {:8.4} {:8.4} {:8.4} {:8.4} {:8.4} {:8.4} {:6.2} {:6.1}",
        gpst_text(t, dec),
        llh[0],
        llh[1],
        llh[2],
        q,
        ns,
        sqvar(cov[1][1]),
        sqvar(cov[0][0]),
        sqvar(cov[2][2]),
        sqvar(cov[1][0]),
        sqvar(cov[0][2]),
        sqvar(cov[2][1]),
        if age.is_finite() { age } else { 0.0 },
        ratio
    );
    if let Some(v) = vel {
        let f = |x: f32| if x.is_finite() { x } else { 0.0 };
        let _ = write!(out, " {:10.5} {:10.5} {:10.5}", f(v[1]), f(v[0]), f(v[2]));
    }
    out.push('\n');
}

/// The `.pos` text of a track.
#[must_use]
pub fn track_pos(res: &PpkResult, track: &[EpochSolution], h: &PosHeader) -> String {
    let with: Vec<&EpochSolution> = track.iter().filter(|e| e.has_position()).collect();
    let mut s = header(
        res,
        h,
        track.first().map(|e| e.t_tag),
        track.last().map(|e| e.t_tag),
    );
    for e in with {
        line(
            &mut s,
            e.t,
            h.time_decimals,
            e.llh,
            e.status.q(),
            e.n_sats,
            &e.cov_enu,
            e.age,
            e.ratio,
            h.velocity.then_some(e.vel_enu),
        );
    }
    s
}

/// The `.pos` text of the result's track ([`PpkResult::epochs`]).
#[must_use]
pub fn epochs_pos(res: &PpkResult, h: &PosHeader) -> String {
    track_pos(res, &res.epochs, h)
}

/// The `.pos` text of the events that have a position: the same layout, one line per event, at the event's
/// time. The position is the point wanted (antenna plus the event's offset); `Q` is the worst status of the
/// epochs interpolated between; `ns` and `ratio` are the smallest among those epochs.
#[must_use]
pub fn events_pos(res: &PpkResult, h: &PosHeader) -> String {
    let mut s = header(
        res,
        h,
        res.epochs.first().map(|e| e.t_tag),
        res.epochs.last().map(|e| e.t_tag),
    );
    for ev in res
        .events
        .iter()
        .filter(|e| e.status != EventStatus::NoSolution)
    {
        let used: Vec<&EpochSolution> = ev
            .epochs_used
            .iter()
            .filter_map(|&i| res.epochs.get(i))
            .collect();
        let ns = used.iter().map(|e| e.n_sats).min().unwrap_or(0);
        let ratio = used.iter().map(|e| e.ratio).fold(f32::INFINITY, f32::min);
        let age = used.iter().map(|e| e.age).fold(0.0, f32::max);
        let c = [
            [ev.sigma_enu[0].powi(2), 0.0, 0.0],
            [0.0, ev.sigma_enu[1].powi(2), 0.0],
            [0.0, 0.0, ev.sigma_enu[2].powi(2)],
        ];
        line(
            &mut s,
            ev.t,
            h.time_decimals,
            ev.llh,
            ev.fix.q(),
            ns,
            &c,
            age,
            if ratio.is_finite() { ratio } else { 0.0 },
            None,
        );
    }
    s
}

fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_owned()
    }
}

/// The column names of [`events_csv`].
pub const EVENTS_CSV_HEADER: &str = "index,label,gps_week,gps_tow_s,gpst,status,fix,q,method,latitude_deg,\
longitude_deg,height_m,x_m,y_m,z_m,antenna_latitude_deg,antenna_longitude_deg,antenna_height_m,sd_e_m,sd_n_m,\
sd_u_m,gap_s,missing_epochs,epochs_used,epochs_q,why";

fn event_csv_line(out: &mut String, ev: &EventResult) {
    let (w, tow) = gps_week_tow(ev.t);
    let a = llh_deg(ev.antenna_xyz);
    let status = match ev.status {
        EventStatus::Ok => "ok",
        EventStatus::Gap => "gap",
        EventStatus::NoSolution => "none",
    };
    let f = |v: f64, d: usize| {
        if v.is_finite() {
            format!("{v:.d$}")
        } else {
            String::new()
        }
    };
    let used: Vec<String> = ev.epochs_used.iter().map(ToString::to_string).collect();
    let qs: Vec<String> = ev.epochs_fix.iter().map(|s| s.q().to_string()).collect();
    let _ = writeln!(
        out,
        "{},{},{},{:.6},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        ev.index + 1,
        csv_field(&ev.label),
        w,
        tow,
        gpst_text(ev.t, 6),
        status,
        ev.fix.name(),
        ev.fix.q(),
        ev.method.name(),
        f(ev.llh[0], 10),
        f(ev.llh[1], 10),
        f(ev.llh[2], 4),
        f(ev.xyz[0], 4),
        f(ev.xyz[1], 4),
        f(ev.xyz[2], 4),
        f(a[0], 10),
        f(a[1], 10),
        f(a[2], 4),
        f(ev.sigma_enu[0], 4),
        f(ev.sigma_enu[1], 4),
        f(ev.sigma_enu[2], 4),
        f(ev.gap_s, 3),
        ev.missing,
        used.join(";"),
        qs.join(";"),
        csv_field(&ev.why)
    );
}

/// One CSV line per event (with or without a position), under [`EVENTS_CSV_HEADER`].
#[must_use]
pub fn events_csv(res: &PpkResult) -> String {
    let mut s = String::from(EVENTS_CSV_HEADER);
    s.push('\n');
    for ev in &res.events {
        event_csv_line(&mut s, ev);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ppk::FixStatus;

    #[test]
    fn time_text_rounds_and_carries() {
        // 2024-02-08 00:00:49 GPST = week 2300, 345649 s
        let t = 2300.0 * 604_800.0 + 345_649.0;
        assert_eq!(gpst_text(t, 3), "2024/02/08 00:00:49.000");
        assert_eq!(gpst_text(t + 10.9996, 3), "2024/02/08 00:01:00.000");
        assert_eq!(gpst_text(t + 0.696_773, 6), "2024/02/08 00:00:49.696773");
        assert_eq!(gpst_text(t, 0), "2024/02/08 00:00:49");
    }

    #[test]
    fn pos_lines_have_the_rtklib_columns() {
        let mut res = PpkResult {
            base_marker: [4_395_470.63, 617_743.11, 4_565_607.21],
            ..PpkResult::default()
        };
        let mut e = EpochSolution::none(2300.0 * 604_800.0 + 345_649.2);
        e.status = FixStatus::Fixed;
        e.set_xyz([4_395_520.0, 617_800.0, 4_565_650.0]);
        e.cov_enu = [
            [1.6e-5, -1e-6, 0.0],
            [-1e-6, 2.5e-5, 4e-6],
            [0.0, 4e-6, 1e-4],
        ];
        e.n_sats = 31;
        e.ratio = 42.31;
        e.age = 0.0;
        res.epochs = vec![e, EpochSolution::none(0.0)];
        let h = PosHeader::new(&PpkOptions::default());
        let text = epochs_pos(&res, &h);
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('%')).collect();
        assert_eq!(body.len(), 1);
        let f: Vec<&str> = body[0].split_whitespace().collect();
        assert_eq!(f.len(), 15);
        assert_eq!(f[0], "2024/02/08");
        assert_eq!(f[1], "00:00:49.200");
        assert_eq!(f[5], "1");
        assert_eq!(f[6], "31");
        assert_eq!(f[7], "0.0050"); // sdn
        assert_eq!(f[8], "0.0040"); // sde
        assert_eq!(f[9], "0.0100"); // sdu
        assert_eq!(f[10], "-0.0010"); // sdne
        assert_eq!(f[12], "0.0020"); // sdun
        assert_eq!(f[14], "42.3");
        // the fixed-width layout of RTKLIB: time (23), then 14.9, 14.9, 10.4 ...
        assert_eq!(&body[0][23..24], " ");
        assert_eq!(body[0][24..38].trim().len(), 12);
        assert!(text.contains("% ref pos   :"));
        assert!(
            text.lines()
                .any(|l| l.starts_with("%  GPST") && l.ends_with("ratio"))
        );
    }
}
