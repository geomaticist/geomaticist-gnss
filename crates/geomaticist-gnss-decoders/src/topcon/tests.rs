//! Synthetic logs: built with the module's own encoders ([`frame::push_message`], [`obs::encode_rd`],
//! [`obs::encode_rp`]) and byte writers for the navigation messages, then read back.

use super::*;
use geomaticist_gnss_core::GnssSystem;
use geomaticist_gnss_core::consts::CLIGHT;
use obs::{SatObs, Signal, carrier, encode_rd, encode_rp, sat_usi};

struct Dir(std::path::PathBuf);

impl Dir {
    fn new(tag: &str) -> Self {
        let d = std::env::temp_dir().join(format!("topcon_mod_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("temp dir");
        Dir(d)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// GPS week of the synthetic log (starts on Sunday 2026-03-15).
const WEEK: i64 = 2410;
const TOW0: u32 = 100_000_000;

fn sat(sys: GnssSystem, prn: u8, fcn: Option<i8>, sigs: &[&[u8; 2]], k: i64) -> SatObs {
    let base: i64 = 1_100_000_000 + i64::from(prn) * 1_000_000 + k * 37;
    let f0 = carrier(sys, sigs[0][0] - b'0', fcn).expect("f0");
    let dop0 = 1500.0 + f64::from(prn) * 10.0;
    let signals = sigs
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let f = carrier(sys, s[0] - b'0', fcn).expect("f");
            let i = i as i64;
            Signal {
                sig: **s,
                range: base + i * 150 - 40,
                phase: base * obs::PHASE_PER_RANGE + 12_345 + i * 100_003 + k * 7,
                doppler: (dop0 * f / f0 * 1000.0).round() as i32,
                snr16: (40 + i as u16) * 16 + 3,
                lock: u16::try_from(10 + k).expect("lock"),
                settled: true,
                half_cycle_known: true,
            }
        })
        .map(|mut g| {
            if g.sig == sigs[0][..] {
                g.range = base;
            }
            g
        })
        .collect();
    SatObs {
        sys,
        prn,
        usi: sat_usi(sys, prn, fcn).expect("usi"),
        fcn,
        info: 0,
        signals,
    }
}

fn epoch(k: i64) -> Epoch {
    Epoch {
        counter: u16::try_from(k).expect("counter"),
        week10: (WEEK % 1024) as u16,
        tow_ms: TOW0 + u32::try_from(k).expect("k") * 1000,
        clock: 123_456 - i32::try_from(k).expect("k") * 1000,
        header_rest: [0; 14],
        sats: vec![
            sat(GnssSystem::Gps, 5, None, &[b"1C", b"1W", b"2W"], k),
            sat(GnssSystem::Glonass, 3, Some(-2), &[b"1C", b"2P"], k),
            sat(GnssSystem::Galileo, 11, None, &[b"1X", b"5X"], k),
            sat(GnssSystem::BeiDou, 30, None, &[b"2I", b"7D"], k),
        ],
    }
}

fn put(b: &mut [u8], o: usize, v: &[u8]) {
    b[o..o + v.len()].copy_from_slice(v);
}

fn ge_body() -> Vec<u8> {
    let mut b = vec![0u8; 122];
    b[0] = 5;
    put(&mut b, 1, &100_000u32.to_le_bytes());
    b[5] = 1;
    put(&mut b, 6, &77i16.to_le_bytes());
    put(&mut b, 8, &100_800i32.to_le_bytes());
    put(&mut b, 14, &((WEEK % 1024) as i16).to_le_bytes());
    put(&mut b, 16, &(-1.0e-8f32).to_le_bytes());
    put(&mut b, 24, &1.0e-12f32.to_le_bytes());
    put(&mut b, 28, &1.0e-4f32.to_le_bytes());
    put(&mut b, 32, &100_800i32.to_le_bytes());
    put(&mut b, 36, &77i16.to_le_bytes());
    put(&mut b, 38, &5153.6f64.to_le_bytes());
    put(&mut b, 46, &0.01f64.to_le_bytes());
    put(&mut b, 54, &0.5f64.to_le_bytes());
    put(&mut b, 70, &0.3f64.to_le_bytes());
    b
}

fn ne_body() -> Vec<u8> {
    let mut b = vec![0u8; 94];
    b[0] = 3;
    b[1] = (-2i8) as u8;
    // tb 06:45 Moscow = 03:45 UTC; tk 06:30
    put(&mut b, 4, &23_400i32.to_le_bytes());
    put(&mut b, 8, &24_300i32.to_le_bytes());
    b[12] = 24;
    put(&mut b, 15, &12_345.678f64.to_le_bytes());
    put(&mut b, 23, &(-9_876.5f64).to_le_bytes());
    put(&mut b, 31, &19_000.25f64.to_le_bytes());
    put(&mut b, 39, &1.5f32.to_le_bytes());
    put(&mut b, 71, &(-2.5e-5f32).to_le_bytes());
    b
}

fn le_body() -> Vec<u8> {
    let mut b = vec![0u8; 129];
    b[0] = 11;
    put(&mut b, 1, &100_000u32.to_le_bytes());
    put(&mut b, 9, &2.0e-9f32.to_le_bytes());
    put(&mut b, 17, &99_600i32.to_le_bytes());
    put(&mut b, 21, &(-3.0e-4f64).to_le_bytes());
    put(&mut b, 37, &55i16.to_le_bytes());
    put(&mut b, 59, &5440.6f64.to_le_bytes());
    put(&mut b, 123, &99_600i32.to_le_bytes());
    b[127] = 107;
    b[128] = 4;
    b
}

fn de_body() -> Vec<u8> {
    let mut b = vec![0u8; 155];
    b[0] = 30;
    put(&mut b, 1, &99_990u32.to_le_bytes());
    put(&mut b, 5, &((WEEK - 1356) as u16).to_le_bytes());
    put(&mut b, 9, &3.5f32.to_le_bytes());
    put(&mut b, 14, &97_200i32.to_le_bytes());
    put(&mut b, 18, &1.0e-5f32.to_le_bytes());
    b[30] = 9;
    put(&mut b, 31, &97_200i32.to_le_bytes());
    put(&mut b, 35, &5282.6f64.to_le_bytes());
    put(&mut b, 123, &1.0e-8f32.to_le_bytes());
    b
}

fn uo_body() -> Vec<u8> {
    let mut b = vec![0u8; 23];
    put(&mut b, 8, &(-2.0e-15f32).to_le_bytes());
    put(&mut b, 12, &319_488u32.to_le_bytes());
    put(&mut b, 16, &((WEEK % 1024) as u16).to_le_bytes());
    b[18] = 18;
    b[19] = 7;
    put(&mut b, 20, &((WEEK % 1024) as u16 + 3).to_le_bytes());
    b[22] = 18;
    b
}

fn io_body() -> Vec<u8> {
    let mut b = vec![0u8; 38];
    put(&mut b, 6, &1.5e-8f32.to_le_bytes());
    put(&mut b, 22, &98_304.0f32.to_le_bytes());
    b
}

/// A whole synthetic log of `n` epochs; `rp`: with the signal lists.
fn log(n: i64, rp: bool) -> Vec<u8> {
    let mut f = Vec::new();
    frame::push_message(
        &mut f,
        *b"JP",
        b"RLOGF JPS Test VR Receiver Log File                ",
        false,
    );
    f.extend_from_slice(b"\r\n");
    frame::push_message(&mut f, *b"PM", b"rcv/ver/main=\"9.9+1\",", true);
    frame::push_message(&mut f, *b"PM", b"rcv/id=\"RID1\",", true);
    frame::push_message(&mut f, *b"PM", b"rcv/sn=\"SN-1\",", true);
    frame::push_message(
        &mut f,
        *b"PM",
        b"/par={ref={{{{W84,N1d00m00.0s,E2d00m00.0s,+10.0000}}},{off,1,{off,0.5}},",
        true,
    );
    frame::push_message(
        &mut f,
        *b"PM",
        b"  any,{{0.0000,0.0000,0.0000},{0.0000,0.0000,1.5000},{0.0060,0.0016,-0.0059},TEST_ANT        NONE,,0},{\"x\",\"y\"}}}",
        true,
    );
    let mut rd = Vec::new();
    rd.extend_from_slice(&2026u16.to_le_bytes());
    rd.extend_from_slice(&[3, 16, 0]);
    frame::push_message(&mut f, *b"RD", &rd, false);
    for k in 0..n {
        let e = epoch(k);
        frame::push_message(&mut f, *b"~~", &[0, 0, 0, 0], false);
        frame::push_message(&mut f, *b"rD", &encode_rd(&e).expect("rD"), false);
        if rp {
            frame::push_message(&mut f, *b"rP", &encode_rp(&e).expect("rP"), false);
        }
        let mut pv = Vec::new();
        for v in [4_395_000.0f64, 618_000.0, 4_566_000.0] {
            pv.extend_from_slice(&(v + k as f64).to_le_bytes());
        }
        pv.extend_from_slice(&[0; 20]);
        pv.push(1);
        frame::push_message(&mut f, *b"PV", &pv, false);
        f.extend_from_slice(b"\r\n");
        if k == 1 {
            frame::push_message(&mut f, *b"UO", &uo_body(), false);
            frame::push_message(&mut f, *b"IO", &io_body(), false);
            for _ in 0..2 {
                frame::push_message(&mut f, *b"GE", &ge_body(), false);
            }
            frame::push_message(&mut f, *b"NE", &ne_body(), false);
            frame::push_message(&mut f, *b"LE", &le_body(), false);
            frame::push_message(&mut f, *b"DE", &de_body(), false);
        }
    }
    f
}

fn read(data: &[u8]) -> (Scan, Vec<(i64, Epoch)>) {
    let mut v = Vec::new();
    let sc = scan(data, 0, &mut |t, e| v.push((t, e.clone())));
    (sc, v)
}

#[test]
fn epochs_round_trip_with_rp() {
    let data = log(5, true);
    assert!(recognises(&data));
    let (sc, v) = read(&data);
    assert!(sc.has_rp);
    assert_eq!(sc.skipped, 0);
    assert_eq!(v.len(), 5);
    for (k, (t, e)) in v.iter().enumerate() {
        let want = epoch(k as i64);
        assert_eq!(e, &want);
        assert_eq!(
            *t,
            (WEEK * 604_800_000 + i64::from(want.tow_ms)) * (TICKS / 1000)
        );
    }
    assert_eq!(sc.dropped, Dropped::default());
}

#[test]
fn signals_inferred_without_rp() {
    let (sc, v) = read(&log(5, false));
    assert!(!sc.has_rp);
    assert_eq!(v.len(), 5);
    for (k, (_, e)) in v.iter().enumerate() {
        assert_eq!(e, &epoch(k as i64));
    }
}

#[test]
fn header_facts_and_navigation() {
    let (sc, _) = read(&log(3, true));
    let i = &sc.info;
    assert_eq!(i.model, "Test VR");
    assert_eq!(i.receiver_id, "RID1");
    assert_eq!(i.firmware, "9.9+1");
    assert_eq!(i.serial, "SN-1");
    assert_eq!(i.antenna, "TEST_ANT        NONE");
    assert_eq!(i.antenna_offset, Some([1.5, 0.0, 0.0]));
    assert_eq!(i.date, Some((2026, 3, 16)));
    assert_eq!(i.leap, Some(18));
    let a = sc.approx.expect("approx");
    assert!((a[0] - 4_395_001.0).abs() < 1e-9);
    // the GE logged twice is one record
    assert_eq!(sc.nav.len(), 4);
    let g = sc.nav.iter().find(|e| e.rec.sat == "G05").expect("G05");
    assert_eq!(
        g.rec.epoch.micros_since_gps0(),
        (WEEK * 604_800 + 100_800) * 1_000_000
    );
    assert_eq!(g.rec.vals[0], Some(f64::from(1.0e-4f32)));
    assert_eq!(g.rec.vals[10], Some(5153.6));
    assert_eq!(g.rec.vals[11], Some(100_800.0));
    assert_eq!(g.rec.vals[21], Some(WEEK as f64));
    assert_eq!(g.rec.vals[28], Some(6.0));
    let r = sc.nav.iter().find(|e| e.rec.sat == "R03").expect("R03");
    assert_eq!((r.rec.epoch.hour, r.rec.epoch.minute), (3, 45));
    assert_eq!(r.rec.vals[0], Some(f64::from(2.5e-5f32)));
    assert_eq!(r.rec.vals[3], Some(12_345.678));
    assert_eq!(r.rec.vals[10], Some(-2.0));
    let e = sc.nav.iter().find(|e| e.rec.sat == "E11").expect("E11");
    assert_eq!(e.rec.vals[20], Some(258.0));
    assert_eq!(e.rec.vals[10], Some(5440.6));
    assert_eq!(e.rec.vals[21], Some(WEEK as f64));
    let c = sc.nav.iter().find(|e| e.rec.sat == "C30").expect("C30");
    assert_eq!(c.rec.vals[21], Some((WEEK - 1356) as f64));
    assert_eq!(c.rec.vals[25], Some(f64::from(3.5f32) * 1e-9));
    assert_eq!(
        c.rec.epoch.micros_since_gps0(),
        (WEEK * 604_800 + 97_200) * 1_000_000
    );
    let (ga, gb) = sc.gps_iono.expect("iono");
    assert_eq!((ga[0], gb[0]), (f64::from(1.5e-8f32), 98_304.0));
    assert!(sc.bds_iono.is_some());
    let u = sc.gps_utc.expect("utc");
    assert_eq!((u.wnt, u.wnlsf, u.dtlsf), (WEEK, WEEK + 3, 18));
}

#[test]
fn damage_loses_one_epoch() {
    let mut data = log(5, true);
    let third = Messages::new(&data)
        .filter(|m| &m.id == b"rD")
        .nth(2)
        .expect("rD")
        .offset;
    data[third + 30] ^= 0x10;
    let (sc, v) = read(&data);
    assert_eq!(v.len(), 4);
    assert!(sc.skipped > 0);
    assert_eq!(sc.gaps, 1);
    // the rP that lost its rD is not used for anything
    assert_eq!(v[2].1, epoch(3));
}

#[test]
fn convert_writes_readable_rinex() {
    let dir = Dir::new("convert");
    let raw = dir.0.join("syn1.tps");
    std::fs::write(&raw, log(5, true)).expect("write");
    assert!(can_convert(&raw));
    assert!(crate::makers::Maker::Topcon.can_decode(&raw));
    let mut msgs = Vec::new();
    let c = convert_with(
        &raw,
        &dir.0,
        &mut |m| msgs.push(m.to_owned()),
        "20260101 000000 UTC",
        None,
    )
    .expect("convert");
    assert!(c.warnings.is_empty(), "{:?}", c.warnings);
    assert_eq!(c.obs.file_name().expect("name"), "syn1.26o");
    assert_eq!(c.nav.len(), 1);
    let text = std::fs::read_to_string(&c.obs).expect("obs");
    assert!(text.contains("RID1                TPS Test VR         9.9+1"));
    assert!(text.contains("TEST_ANT        NONE"));
    assert!(text.contains("  1 R03 -2"));
    assert!(text.contains("G L1W  0.25000"));
    assert!(!text.contains('\r'));
    let st = geomaticist_gnss_rinex::obs::read_obs_kinds(&c.obs, "", None, None, "GRECJ", "CLSD")
        .expect("read back");
    assert_eq!(st.t.len(), 5);
    assert_eq!(st.delta_hen, [1.5, 0.0, 0.0]);
    let want = epoch(2);
    for s in &want.sats {
        let sat = geomaticist_gnss_core::system::Sat::new(s.sys, s.prn);
        let got = st.sats.get(&sat).expect("satellite");
        for g in &s.signals {
            let t = String::from_utf8_lossy(&g.sig).into_owned();
            let val = |k: &str| {
                let code =
                    geomaticist_gnss_core::signal::ObsCode::new(format!("{k}{t}").as_bytes())
                        .expect("code");
                got.get(&code).expect("signal").values[2]
            };
            assert!((val("C") - g.range_m()).abs() < 5e-4, "{} C{t}", sat_id(s));
            let f = carrier(s.sys, g.sig[0] - b'0', s.fcn).expect("f");
            let shift = rinex::phase_shift(s.sys, g.sig).unwrap_or(0.0);
            let l = g.phase as f64 * obs::PHASE_UNIT * f / CLIGHT - shift;
            assert!((val("L") - l).abs() < 6e-4, "{} L{t}", sat_id(s));
            assert!((val("D") - g.doppler_hz()).abs() < 5e-4);
            assert!((val("S") - g.snr()).abs() < 1e-3);
        }
    }
    let nav = geomaticist_gnss_rinex::nav::read_nav(&c.nav);
    assert!(nav.bad.is_empty(), "{:?}", nav.bad);
    assert_eq!(nav.kep.len(), 3);
    assert_eq!(nav.glo.len(), 1);
    // to another folder through the crate's entry point
    let out2 = dir.0.join("out2");
    std::fs::create_dir_all(&out2).expect("dir");
    let c2 = crate::convert(&raw, &out2, &mut |_| {}).expect("convert");
    assert!(c2.obs.exists());
}

fn sat_id(s: &SatObs) -> String {
    format!("{}{:02}", s.sys.letter(), s.prn)
}

#[test]
fn not_a_log() {
    assert!(!recognises(b"hello world, not GREIS"));
    let mut f = Vec::new();
    frame::push_message(&mut f, *b"JP", b"RLOGF JPS X Receiver Log File", false);
    frame::push_message(&mut f, *b"PM", b"rcv/id=\"A\",", true);
    // framed, but without measurements
    assert!(!recognises(&f));
}
