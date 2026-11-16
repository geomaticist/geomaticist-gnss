//! `gnss ppk` end to end with files: simulated base and rover observations written as RINEX, the simulated
//! orbits as an SP3 file, the events as a `Timestamp.MRK` file; the command reads them, solves, and writes
//! the `.pos` and CSV outputs.

#![allow(
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::many_single_char_names,
    clippy::too_many_lines
)]

use std::fmt::Write as _;

use clap::Parser;
use geomaticist_gnss_cli::{Cli, run};
use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::time::{TICKS, from_ticks, gps_week_tow};
use geomaticist_gnss_engine::ppk::enu_at;
use geomaticist_gnss_engine::ppk::synth::{SynthOrbit, SynthReceiver, enu_point, simulate};
use geomaticist_gnss_rinex::station::{StationHeader, station_to_rinex3};

const T0: f64 = 2300.0 * 604_800.0 + 345_600.0;
const BASE: [f64; 3] = [4_395_470.63, 617_743.11, 4_565_607.21];

fn path(t: f64) -> [f64; 3] {
    let s = t - T0;
    enu_point(
        BASE,
        [120.0 * (0.1 * s).sin() + 40.0, 60.0 * (0.2 * s).sin(), 90.0],
    )
}

/// The constellation as SP3 text: positions every 5 minutes around the session, clocks in microseconds.
fn sp3(orbit: &SynthOrbit) -> String {
    let mut s = String::from("#cP2024  2  8  0  0  0.00000000      25 ORBIT IGS20 FIT  SIM\n");
    for k in -12..=12 {
        let t = T0 + 300.0 * f64::from(k);
        let c = from_ticks((t * TICKS as f64).round() as i64);
        let _ = writeln!(
            s,
            "*  {:4} {:2} {:2} {:2} {:2} {:11.8}",
            c.year,
            c.month,
            c.day,
            c.hour,
            c.minute,
            c.sec_ticks as f64 / TICKS as f64
        );
        for sat in &orbit.sats {
            let p = sat.position(t - orbit.t0);
            let clk = sat.clk + sat.drift * (t - orbit.t0);
            let _ = writeln!(
                s,
                "P{}{:14.6}{:14.6}{:14.6}{:14.6}",
                sat.sat,
                p[0] / 1e3,
                p[1] / 1e3,
                p[2] / 1e3,
                clk * 1e6
            );
        }
    }
    s.push_str("EOF\n");
    s
}

#[test]
fn files_in_positions_out() {
    let dir = std::env::temp_dir().join(format!("gnss-ppk-e2e-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let orbit = SynthOrbit::constellation(
        &[GnssSystem::Gps, GnssSystem::Galileo, GnssSystem::BeiDou],
        6,
        4,
        T0,
    );
    let epochs: Vec<f64> = (0..=200).map(|k| T0 + 0.2 * f64::from(k)).collect();
    let base = simulate(&orbit, &SynthReceiver::new("BASE", 3), &|_| BASE, &epochs);
    let rover = simulate(&orbit, &SynthReceiver::new("ROVR", 4), &path, &epochs);
    let header = StationHeader {
        program: "test".into(),
        ..StationHeader::default()
    };
    let p = |name: &str| dir.join(name);
    std::fs::write(p("base.obs"), station_to_rinex3(&base, &header)).unwrap();
    std::fs::write(p("rover.obs"), station_to_rinex3(&rover, &header)).unwrap();
    std::fs::write(p("orbit.sp3"), sp3(&orbit)).unwrap();
    // three photos, the last one after the data
    let (week, tow) = gps_week_tow(T0);
    let mut mrk = String::new();
    for (k, dt) in [10.07, 20.33, 60.0].iter().enumerate() {
        let _ = writeln!(
            mrk,
            "{}\t{:.6}\t[{week}]\t  -120,N\t    45,E\t   193,V\t46.0,Lat\t8.0,Lon\t500.0,Ellh\t1.0, 1.0, 2.0\t16,Q",
            k + 1,
            tow + dt
        );
    }
    std::fs::write(p("flight_Timestamp.MRK"), mrk).unwrap();
    let xyz = format!("{},{},{}", BASE[0], BASE[1], BASE[2]);
    let s = |x: &std::path::Path| x.to_str().unwrap().to_owned();
    let cli = Cli::try_parse_from([
        "gnss".to_owned(),
        "ppk".into(),
        s(&p("base.obs")),
        s(&p("rover.obs")),
        "--base-xyz".into(),
        xyz,
        "--sp3".into(),
        s(&p("orbit.sp3")),
        "--mrk".into(),
        s(&p("flight_Timestamp.MRK")),
        "-o".into(),
        s(&p("out.pos")),
        "--json".into(),
    ])
    .unwrap();
    let v: serde_json::Value = serde_json::from_str(&run(&cli).unwrap()).unwrap();
    assert_eq!(v["summary"]["n_epochs"], 201);
    assert!(
        v["summary"]["fix_percent"].as_f64().unwrap() > 95.0,
        "{}",
        v["summary"]
    );
    assert_eq!(v["events"].as_array().unwrap().len(), 3);
    assert_eq!(v["events"][0]["status"], "ok");
    assert_eq!(v["events"][0]["fix"], "fixed");
    assert_eq!(v["events"][2]["status"], "nosolution");
    assert_eq!(v["written"].as_array().unwrap().len(), 3);
    // the first photo: the camera is the antenna plus the offset; the antenna is on the path
    let ant: Vec<f64> = v["events"][0]["antenna_xyz"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    let cam: Vec<f64> = v["events"][0]["xyz"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect();
    let t = v["events"][0]["t"].as_f64().unwrap();
    let d = enu_at(path(t), [ant[0], ant[1], ant[2]]);
    assert!(
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < 0.03,
        "{d:?}"
    );
    let o = enu_at([ant[0], ant[1], ant[2]], [cam[0], cam[1], cam[2]]);
    assert!(
        (o[0] - 0.045).abs() < 1e-6 && (o[1] + 0.120).abs() < 1e-6 && (o[2] + 0.193).abs() < 1e-6,
        "{o:?}"
    );
    // the files: RTKLIB layout, one line per epoch; two event lines; the CSV has all three events
    let pos = std::fs::read_to_string(p("out.pos")).unwrap();
    let body: Vec<&str> = pos.lines().filter(|l| !l.starts_with('%')).collect();
    assert_eq!(body.len(), 201);
    assert!(body.iter().all(|l| l.split_whitespace().count() == 15));
    assert!(pos.contains("% pos mode  : Kinematic"));
    let ev = std::fs::read_to_string(p("out_events.pos")).unwrap();
    assert_eq!(ev.lines().filter(|l| !l.starts_with('%')).count(), 2);
    let csv = std::fs::read_to_string(p("out_events.csv")).unwrap();
    assert_eq!(csv.lines().count(), 4);
    // the events can also be written into the rover file and read back from it
    let ins = Cli::try_parse_from([
        "gnss".to_owned(),
        "events".into(),
        "insert".into(),
        s(&p("rover.obs")),
        "--mrk".into(),
        s(&p("flight_Timestamp.MRK")),
        "-o".into(),
        s(&p("rover_events.obs")),
        "--json".into(),
    ])
    .unwrap();
    let w: serde_json::Value = serde_json::from_str(&run(&ins).unwrap()).unwrap();
    assert_eq!(
        (w["inserted"].as_u64(), w["after"].as_u64()),
        (Some(2), Some(1))
    );
    let _ = std::fs::remove_dir_all(&dir);
}
