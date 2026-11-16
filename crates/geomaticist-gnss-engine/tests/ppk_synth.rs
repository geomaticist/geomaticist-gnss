//! Kinematic processing end to end on simulated data: a rover on a known trajectory against a static base,
//! observations generated with the engine's own models (`ppk::synth`). The truth is known exactly, so the
//! tests state errors in millimetres.

#![allow(
    clippy::unwrap_used,
    clippy::many_single_char_names,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::field_reassign_with_default,
    clippy::too_many_lines
)]

use std::collections::BTreeMap;

use geomaticist_gnss_core::event::Event;
use geomaticist_gnss_core::{GnssSystem, Sat, StationObs};
use geomaticist_gnss_engine::ppk::synth::{
    SynthOrbit, SynthReceiver, SynthSlip, enu_point, simulate,
};
use geomaticist_gnss_engine::ppk::{
    self, BaseStation, Direction, EpochSolution, EventInterp, EventStatus, FixStatus, PpkInput,
    PpkOptions, PpkResult, enu_at,
};
use geomaticist_gnss_engine::tropo::TropoInputs;

const T0: f64 = 2300.0 * 604_800.0 + 345_600.0;
const BASE: [f64; 3] = [4_395_470.63, 617_743.11, 4_565_607.21];

/// A drone-like path around the base: a figure of eight 100 m up, up to 17 m/s and 2.5 m/s².
fn path(t: f64) -> [f64; 3] {
    let s = t - T0;
    enu_point(
        BASE,
        [
            150.0 * (0.08 * s).sin() + 60.0,
            80.0 * (0.16 * s).sin() - 40.0,
            100.0 + 3.0 * (0.05 * s).sin(),
        ],
    )
}

struct Scene {
    orbit: SynthOrbit,
    base: StationObs,
    rover: StationObs,
}

fn scene(seconds: f64, rate_hz: f64, rover: &SynthReceiver, base_rate_hz: f64) -> Scene {
    let orbit = SynthOrbit::constellation(
        &[GnssSystem::Gps, GnssSystem::Galileo, GnssSystem::BeiDou],
        6,
        4,
        T0,
    );
    let grid = |hz: f64| -> Vec<f64> {
        (0..=(seconds * hz).round() as usize)
            .map(|k| T0 + k as f64 / hz)
            .collect()
    };
    let mut b = SynthReceiver::new("BASE", 11);
    b.clock = 3.0e-7;
    b.clock_drift = 2.0e-9;
    let base = simulate(&orbit, &b, &|_| BASE, &grid(base_rate_hz));
    let rover = simulate(&orbit, rover, &path, &grid(rate_hz));
    Scene { orbit, base, rover }
}

fn run(sc: &Scene, opts: &PpkOptions, events: &[Event]) -> PpkResult {
    let orb = sc.orbit.orbits();
    let fcn: BTreeMap<Sat, i32> = BTreeMap::new();
    let tropo = TropoInputs::default();
    let mut input = PpkInput {
        base: &sc.base,
        rover: &sc.rover,
        orbits: &orb,
        fcn: &fcn,
        base_station: BaseStation {
            marker_xyz: BASE,
            delta_hen: None,
        },
        antex: None,
        tropo: &tropo,
        events,
        occupations: &[],
    };
    ppk::solve(&mut input, opts, &mut |_| true).unwrap()
}

/// 3-D error of an epoch against the path at its receive time, m.
fn err(e: &EpochSolution) -> f64 {
    let d = enu_at(path(e.t), e.xyz);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn rover_rx(seed: u64) -> SynthReceiver {
    let mut r = SynthReceiver::new("ROVR", seed);
    // a receiver clock 0.3 ms off and drifting: the positions belong to the receive time, not the tag
    r.clock = 3.0e-4;
    r.clock_drift = 1.0e-7;
    r
}

#[test]
fn fixes_a_flight_to_millimetres() {
    let sc = scene(90.0, 5.0, &rover_rx(5), 5.0);
    let mut o = PpkOptions::default();
    o.keep_passes = true;
    let res = run(&sc, &o, &[]);
    let s = &res.summary;
    assert_eq!(s.n_epochs, 451);
    assert!(s.fix_percent > 99.0, "{s:?}");
    assert!(
        s.fix_percent_forward > 95.0 && s.fix_percent_backward > 95.0,
        "{s:?}"
    );
    assert_eq!(s.fb_conflicts, 0);
    assert!(s.fb_max < 0.01, "{s:?}");
    // 2 mm of phase noise at the zenith: a few millimetres in position, the height weakest
    let fixed: Vec<f64> = res
        .epochs
        .iter()
        .filter(|e| e.status == FixStatus::Fixed)
        .map(err)
        .collect();
    let rms = (fixed.iter().map(|e| e * e).sum::<f64>() / fixed.len() as f64).sqrt();
    let worst = fixed.iter().copied().fold(0.0, f64::max);
    assert!(rms < 0.008, "rms of the fixed errors {rms} m");
    assert!(worst < 0.03, "worst fixed error {worst} m");
    assert!(s.slips < 60, "{s:?}");
    // the receive time is the tag less the receiver clock offset (0.3 ms and drifting)
    let e = &res.epochs[100];
    let clock = 3.0e-4 + 1.0e-7 * (e.t_tag - T0);
    assert!((e.t_tag - e.t - clock).abs() < 1e-6, "{}", e.t_tag - e.t);
    // realistic σ: no fixed epoch errs by more than 4 σ (3-D)
    for e in res.epochs.iter().filter(|e| e.status == FixStatus::Fixed) {
        let sg = e.sigma_enu();
        let s3 = (sg[0] * sg[0] + sg[1] * sg[1] + sg[2] * sg[2]).sqrt();
        assert!(err(e) < 4.0 * s3, "{} vs σ {s3}", err(e));
    }
    // the Doppler velocity agrees with the path's
    let v = res.epochs[200].vel_enu;
    let (a, b) = (
        path(res.epochs[200].t - 0.05),
        path(res.epochs[200].t + 0.05),
    );
    let tv = enu_at(a, b);
    for i in 0..3 {
        assert!((f64::from(v[i]) - tv[i] / 0.1).abs() < 0.1, "{v:?} {tv:?}");
    }
    // the first forward epochs are not fixed (the integers must repeat for a second first)
    assert_ne!(res.forward[0].status, FixStatus::Fixed);
    assert_eq!(res.forward[300].status, FixStatus::Fixed);
    assert_eq!(res.backward[0].status, FixStatus::Fixed);
}

#[test]
fn recovers_from_slips_and_a_gap() {
    let mut r = rover_rx(8);
    let g = |prn: u8| Sat::new(GnssSystem::Gps, prn);
    // which GPS satellites are in view is not known in advance: slip several
    for prn in 1..=24u8 {
        if prn % 3 == 0 {
            // an unflagged one-cycle slip on L1, and a flagged large one on L2 later
            r.slips.push(SynthSlip {
                sat: g(prn),
                band: 1,
                t: T0 + 20.0,
                cycles: 1.0,
                flagged: false,
            });
            r.slips.push(SynthSlip {
                sat: g(prn),
                band: 2,
                t: T0 + 40.0,
                cycles: -37.0,
                flagged: true,
            });
        }
    }
    // an unflagged half-cycle slip on every Galileo satellite at once
    for prn in 1..=24u8 {
        r.slips.push(SynthSlip {
            sat: Sat::new(GnssSystem::Galileo, prn),
            band: 1,
            t: T0 + 30.0,
            cycles: 0.5,
            flagged: false,
        });
    }
    // two epochs lost
    r.missing = vec![T0 + 50.0, T0 + 50.2];
    let sc = scene(80.0, 5.0, &r, 5.0);
    let res = run(&sc, &PpkOptions::default(), &[]);
    let s = &res.summary;
    assert!(s.slips > 40, "{s:?}");
    assert!(s.fix_percent > 97.0, "{s:?}");
    let worst = res
        .epochs
        .iter()
        .filter(|e| e.status == FixStatus::Fixed)
        .map(err)
        .fold(0.0, f64::max);
    assert!(worst < 0.03, "worst fixed error {worst} m");
    assert_eq!(s.fb_conflicts, 0);
}

/// A receiver under stress: phases dropped for two epochs and back with a jump of whole cycles and the
/// loss-of-lock flag held up for a second, as a hard brake produces on the highest satellites. The jumps
/// are taken out and the ambiguities carry on; nothing has to be found again.
#[test]
fn repairs_whole_cycle_slips_across_a_short_break() {
    let clean = scene(40.0, 5.0, &rover_rx(13), 5.0);
    let base_line = run(&clean, &PpkOptions::default(), &[]);
    let mut r = rover_rx(13);
    for prn in 1..=24u8 {
        for (sys, band) in [
            (GnssSystem::Gps, 1u8),
            (GnssSystem::Gps, 2),
            (GnssSystem::Galileo, 1),
        ] {
            if prn % 2 == 0 {
                let sat = Sat::new(sys, prn);
                r.phase_gaps.push((sat, band, T0 + 20.0, T0 + 20.2));
                r.slips.push(SynthSlip {
                    sat,
                    band,
                    t: T0 + 20.4,
                    cycles: f64::from(prn) - 30.0,
                    flagged: true,
                });
                r.flags.push((sat, band, T0 + 20.4, T0 + 21.4));
            }
        }
    }
    let sc = scene(40.0, 5.0, &r, 5.0);
    let res = run(&sc, &PpkOptions::default(), &[]);
    let s = &res.summary;
    assert!(s.slips_repaired >= 8, "{s:?}");

    // about half are repaired here: the simulated noise grows as 1 / sin(elevation), and the low
    // satellites, whose phases are too noisy to tell a whole cycle by, start again
    let restarted = s.slips - base_line.summary.slips;
    assert!(
        s.slips_repaired + 3 >= restarted,
        "repaired {} restarted {restarted}",
        s.slips_repaired
    );
    assert!(s.fix_percent > 99.0, "{s:?}");
    assert_eq!(s.fb_conflicts, 0);
    let worst = res
        .epochs
        .iter()
        .filter(|e| e.status == FixStatus::Fixed)
        .map(err)
        .fold(0.0, f64::max);
    assert!(worst < 0.03, "worst fixed error {worst} m");
    // the epochs right after the break are fixed, with all the satellites back
    let k = res
        .epochs
        .iter()
        .position(|e| (e.t_tag - (T0 + 20.4)).abs() < 1e-3)
        .unwrap();
    assert_eq!(res.epochs[k].status, FixStatus::Fixed);
    assert!(res.epochs[k].n_sats + 1 >= res.epochs[k - 5].n_sats);
}

#[test]
fn a_doppler_of_the_other_sign_is_recognised() {
    let mut r = rover_rx(17);
    r.doppler_inverted = true;
    let sc = scene(20.0, 5.0, &r, 5.0);
    let res = run(&sc, &PpkOptions::default(), &[]);
    assert!(
        res.warnings.iter().any(|w| w.contains("Doppler")),
        "{:?}",
        res.warnings
    );
    let e = &res.epochs[50];
    let tv = enu_at(path(e.t - 0.05), path(e.t + 0.05));
    for i in 0..3 {
        assert!(
            (f64::from(e.vel_enu[i]) - tv[i] / 0.1).abs() < 0.1,
            "{:?} {tv:?}",
            e.vel_enu
        );
    }
    // and the Hermite events use those velocities
    let mut o = PpkOptions::default();
    o.event_interp = EventInterp::Hermite2;
    let ev = [Event::at(T0 + 10.07)];
    let h = run(&sc, &o, &ev);
    assert_eq!(h.events[0].method, EventInterp::Hermite2);
    let d = enu_at(path(h.events[0].t), h.events[0].antenna_xyz);
    assert!(
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() < 0.015,
        "{d:?}"
    );
}

#[test]
fn a_slower_base_is_interpolated() {
    // base at 1 Hz, rover at 5 Hz
    let sc = scene(60.0, 5.0, &rover_rx(21), 1.0);
    let res = run(&sc, &PpkOptions::default(), &[]);
    assert!(res.summary.fix_percent > 99.0, "{:?}", res.summary);
    let worst = res
        .epochs
        .iter()
        .filter(|e| e.status == FixStatus::Fixed)
        .map(err)
        .fold(0.0, f64::max);
    assert!(worst < 0.03, "worst fixed error {worst} m");
    assert!(res.epochs[3].age > 0.39 && res.epochs[3].age < 0.41);
}

#[test]
fn events_on_a_fixed_track() {
    let sc = scene(60.0, 5.0, &rover_rx(33), 5.0);
    let events: Vec<Event> = (0..25)
        .map(|k| Event::at(T0 + 3.0 + 2.137 * f64::from(k)))
        .chain([Event::at(T0 - 1.0), Event::at(T0 + 61.0)])
        .collect();
    // the limits are the noise of the fixed positions (a few millimetres, the height most) plus, for the
    // linear form, the curvature of the path between epochs
    for (method, limit) in [(EventInterp::Linear2, 0.035), (EventInterp::Cubic4, 0.02)] {
        let mut o = PpkOptions::default();
        o.event_interp = method;
        let res = run(&sc, &o, &events);
        assert_eq!(res.events.len(), 27);
        assert_eq!(res.summary.events_failed, 2);
        for ev in &res.events[..25] {
            assert_eq!(ev.status, EventStatus::Ok);
            assert_eq!(ev.fix, FixStatus::Fixed);
            let d = enu_at(path(ev.t), ev.antenna_xyz);
            let e = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!(e < limit, "{method:?}: {e} m at {}", ev.t - T0);
        }
        assert!(res.events[25].why.contains("before"));
        assert!(res.events[26].why.contains("after"));
    }
}

/// One epoch per second, and one every ten seconds: the limits that are given in seconds (the longest
/// break an ambiguity survives, the confirmation of an integer) must not turn every epoch into a break.
#[test]
fn low_rates() {
    for (hz, seconds, least) in [(1.0, 120.0, 95.0), (0.1, 1200.0, 90.0)] {
        let sc = scene(seconds, hz, &rover_rx(41), hz);
        let res = run(&sc, &PpkOptions::default(), &[]);
        let s = &res.summary;
        assert!(s.fix_percent > least, "{hz} Hz: {s:?}");
        // new ambiguities are first appearances and risings, not one per epoch
        assert!(s.slips < 80, "{hz} Hz: {s:?}");
        let worst = res
            .epochs
            .iter()
            .filter(|e| e.status == FixStatus::Fixed)
            .map(err)
            .fold(0.0, f64::max);
        assert!(worst < 0.04, "{hz} Hz: worst fixed error {worst} m");
    }
}

#[test]
fn directions_and_float_only() {
    let sc = scene(40.0, 2.0, &rover_rx(3), 2.0);
    let mut o = PpkOptions::default();
    o.direction = Direction::Forward;
    let f = run(&sc, &o, &[]);
    assert!(f.summary.fix_percent > 90.0, "{:?}", f.summary);
    o.direction = Direction::Backward;
    let b = run(&sc, &o, &[]);
    assert!(b.summary.fix_percent > 90.0, "{:?}", b.summary);
    assert_eq!(b.epochs[0].status, FixStatus::Fixed);
    o.ar = ppk::ArMode::Off;
    let fl = run(&sc, &o, &[]);
    assert_eq!(fl.summary.n_fixed, 0);
    // a float solution converges backward: the first epochs are within a few decimetres
    assert!(err(&fl.epochs[0]) < 0.5, "{}", err(&fl.epochs[0]));
}

#[test]
fn cancel_and_bad_input() {
    let sc = scene(10.0, 1.0, &rover_rx(3), 1.0);
    let orb = sc.orbit.orbits();
    let fcn = BTreeMap::new();
    let tropo = TropoInputs::default();
    let mut input = PpkInput {
        base: &sc.base,
        rover: &sc.rover,
        orbits: &orb,
        fcn: &fcn,
        base_station: BaseStation {
            marker_xyz: BASE,
            delta_hen: None,
        },
        antex: None,
        tropo: &tropo,
        events: &[],
        occupations: &[],
    };
    let r = ppk::solve(&mut input, &PpkOptions::default(), &mut |_| false);
    assert_eq!(r.unwrap_err(), ppk::PpkError::Cancelled);
    input.base_station.marker_xyz = [1.0, 2.0, 3.0];
    let r = ppk::solve(&mut input, &PpkOptions::default(), &mut |_| true);
    assert!(matches!(r, Err(ppk::PpkError::Input(_))));
}
