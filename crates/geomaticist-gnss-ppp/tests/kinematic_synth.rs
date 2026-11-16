//! Kinematic precise point positioning end to end on simulated data: a receiver on a known path, no base,
//! observations generated with the engine's own models (`ppk::synth`).

#![allow(
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::many_single_char_names
)]

use std::collections::BTreeMap;

use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::event::Event;
use geomaticist_gnss_engine::ppk::synth::{SynthOrbit, SynthReceiver, enu_point, simulate};
use geomaticist_gnss_engine::ppk::{Direction, EventStatus, FixStatus, enu_at};
use geomaticist_gnss_engine::tropo::TropoInputs;
use geomaticist_gnss_ppp::kinematic::{KinematicInput, KinematicOptions, solve};

const T0: f64 = 2300.0 * 604_800.0 + 345_600.0;
const HOME: [f64; 3] = [4_395_470.63, 617_743.11, 4_565_607.21];

/// A slow figure of eight, 120 m up.
fn path(t: f64) -> [f64; 3] {
    let s = t - T0;
    enu_point(
        HOME,
        [
            300.0 * (0.004 * s).sin(),
            150.0 * (0.008 * s).sin(),
            120.0 + 5.0 * (0.003 * s).sin(),
        ],
    )
}

fn err(t: f64, xyz: [f64; 3]) -> f64 {
    let d = enu_at(path(t), xyz);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

#[test]
fn a_moving_receiver_converges_and_events_follow() {
    let orbit = SynthOrbit::constellation(&[GnssSystem::Gps, GnssSystem::Galileo], 6, 4, T0);
    let mut rx = SynthReceiver::new("ROVR", 7);
    rx.signals.insert('E', vec!["L1C".into(), "L5Q".into()]);
    rx.signals.remove(&'C');
    rx.clock = 2.0e-4;
    // two hours at 30 s
    let epochs: Vec<f64> = (0..=240).map(|k| T0 + 30.0 * f64::from(k)).collect();
    let st = simulate(&orbit, &rx, &path, &epochs);
    let orb = orbit.orbits();
    let fcn = BTreeMap::new();
    let osb = BTreeMap::new();
    let tropo = TropoInputs::default();
    let mut opts = KinematicOptions::default();
    // the simulation has no relativity, wind-up or tides
    opts.ppp.shapiro = false;
    opts.ppp.windup = false;
    opts.ppp.tides = false;
    opts.ppp.pole_tide = false;
    opts.ppp.shadow = false;
    opts.ppp.systems = "GE".into();
    opts.keep_passes = true;
    let events = [Event::at(T0 + 3600.0 + 7.3), Event::at(T0 - 5.0)];
    let mut input = KinematicInput {
        obs: &st,
        orbits: &orb,
        antex: None,
        fcn: &fcn,
        osb: &osb,
        erp: None,
        blq: None,
        tropo: &tropo,
        x0: None,
        events: &events,
        occupations: &[],
    };
    let res = solve(&mut input, &opts, &mut |_| true).unwrap();
    assert_eq!(res.epochs.len(), 241);
    assert!(
        res.epochs.iter().all(|e| e.status == FixStatus::Ppp),
        "{:?}",
        res.summary
    );
    assert_eq!(res.epochs[10].status.q(), 6);
    // the receive time is the tag less the clock offset
    assert!((res.epochs[10].t_tag - res.epochs[10].t - 2.0e-4).abs() < 1e-6);
    // forward: metres or decimetres at the start, converged in the second hour
    let late: Vec<f64> = res.forward[150..].iter().map(|e| err(e.t, e.xyz)).collect();
    let worst_late = late.iter().copied().fold(0.0, f64::max);
    assert!(worst_late < 0.25, "forward, second hour: {worst_late} m");
    // merged: the whole track within a few decimetres, its middle better
    let all: Vec<f64> = res.epochs.iter().map(|e| err(e.t, e.xyz)).collect();
    let worst = all.iter().copied().fold(0.0, f64::max);
    assert!(worst < 0.6, "merged: {worst} m");
    let mut mid = all[60..180].to_vec();
    mid.sort_by(f64::total_cmp);
    assert!(
        mid[mid.len() / 2] < 0.12,
        "merged, middle: median {} m",
        mid[mid.len() / 2]
    );
    // the σ reported covers the error
    for e in &res.epochs[60..180] {
        let s = e.sigma_enu();
        let s3 = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
        assert!(err(e.t, e.xyz) < 4.0 * s3, "{} vs {s3}", err(e.t, e.xyz));
    }
    // events: the same interface as the differential solution
    assert_eq!(res.events[0].status, EventStatus::Ok);
    assert_eq!(res.events[0].fix, FixStatus::Ppp);
    assert!(err(res.events[0].t, res.events[0].antenna_xyz) < 0.3);
    assert_eq!(res.events[1].status, EventStatus::NoSolution);
    // one direction only
    opts.direction = Direction::Forward;
    opts.keep_passes = false;
    let f = solve(&mut input, &opts, &mut |_| true).unwrap();
    assert_eq!(f.epochs.len(), 241);
    assert_eq!(f.forward.len(), 0);
}
