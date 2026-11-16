use super::*;

const WN: u16 = 2300;

/// A block body (the bytes after `WNc`) written at the reference's block offsets.
#[derive(Default)]
struct Body(Vec<u8>);

impl Body {
    fn put(&mut self, off: usize, bytes: &[u8]) -> &mut Self {
        let o = off - 14;
        if self.0.len() < o + bytes.len() {
            self.0.resize(o + bytes.len(), 0);
        }
        self.0[o..o + bytes.len()].copy_from_slice(bytes);
        self
    }
    fn f4(&mut self, off: usize, v: f32) -> &mut Self {
        self.put(off, &v.to_le_bytes())
    }
    fn f8(&mut self, off: usize, v: f64) -> &mut Self {
        self.put(off, &v.to_le_bytes())
    }
    fn u2(&mut self, off: usize, v: u16) -> &mut Self {
        self.put(off, &v.to_le_bytes())
    }
    fn u4(&mut self, off: usize, v: u32) -> &mut Self {
        self.put(off, &v.to_le_bytes())
    }
    fn len(&mut self, n: usize) -> &mut Self {
        self.0.resize(n - 14, 0);
        self
    }
}

/// A MeasEpochChannelType1 sub-block (20 bytes).
#[allow(clippy::too_many_arguments)]
fn t1(
    typ: u8,
    svid: u8,
    code_mm: u64,
    carr: i32,
    dop: i32,
    cn0: u8,
    lock: u16,
    info: u8,
    n2: u8,
) -> Vec<u8> {
    let mut s = vec![0u8; 20];
    s[1] = typ;
    s[2] = svid;
    s[3] = (code_mm >> 32) as u8 & 15;
    s[4..8].copy_from_slice(&(code_mm as u32).to_le_bytes());
    s[8..12].copy_from_slice(&dop.to_le_bytes());
    s[12..14].copy_from_slice(&(carr as u16).to_le_bytes());
    s[14] = (carr >> 16) as u8;
    s[15] = cn0;
    s[16..18].copy_from_slice(&lock.to_le_bytes());
    s[18] = info;
    s[19] = n2;
    s
}

/// A MeasEpochChannelType2 sub-block (12 bytes); offsets in mm, cycles/1000 and 1e-4 Hz.
fn t2(typ: u8, lock: u8, cn0: u8, code_off: i32, carr: i32, dop_off: i32, info: u8) -> Vec<u8> {
    let mut q = vec![0u8; 12];
    q[0] = typ;
    q[1] = lock;
    q[2] = cn0;
    q[3] = ((code_off >> 16) as u8 & 7) | (((dop_off >> 16) as u8 & 31) << 3);
    q[4] = (carr >> 16) as u8;
    q[5] = info;
    q[6..8].copy_from_slice(&(code_off as u16).to_le_bytes());
    q[8..10].copy_from_slice(&(carr as u16).to_le_bytes());
    q[10..12].copy_from_slice(&(dop_off as u16).to_le_bytes());
    q
}

fn meas(tow_ms: u32, flags: u8, subs: &[Vec<u8>], n1: u8) -> Vec<u8> {
    let mut body = vec![n1, 20, 12, flags, 0, 0];
    for s in subs {
        body.extend_from_slice(s);
    }
    block(MEAS_EPOCH, 1, tow_ms, WN, &body)
}

fn find<'a>(r: &'a Sbf, sat: Sat, code: &str) -> &'a geomaticist_gnss_core::obs::Signal {
    r.obs.sats[&sat]
        .iter()
        .find(|(c, _)| c.to_string() == code)
        .map_or_else(|| panic!("{sat:?} {code}"), |(_, s)| s)
}

#[test]
fn crc_check_value() {
    // CRC-16/XMODEM check value
    assert_eq!(crc16(b"123456789"), 0x31C3);
}

#[test]
fn meas_epoch_signals_and_slips() {
    let pr = 21_000_000_123u64; // mm
    let tow = 345_600_000u32;
    let first = meas(
        tow,
        0,
        &[
            // G05 L1CA with an L2C sub-block
            t1(0, 5, pr, 1_234_567, -12_345_678, 180, 500, 0, 1),
            t2(3, 200, 160, -70_000, -2_000_000, -5, 0),
            // R03 L1CA, frequency number +1 in ObsInfo
            t1(8, 40, 20_000_000_000, 1000, 50_000, 170, 100, 9 << 3, 0),
            // E11 E1
            t1(17, 81, 24_000_000_000, 0, 0, 176, 100, 0, 0),
            // C20 B1I
            t1(28, 160, 22_000_000_000, 0, 0, 160, 100, 0, 0),
            // J01 L1C through the extended signal number (31 + ObsInfo)
            t1(31, 181, 37_000_000_000, 0, 0, 160, 100, 0, 0),
            // J01 L6: no RINEX 3 code here
            t1(27, 181, 37_000_000_000, 0, 0, 160, 100, 0, 0),
            // a second antenna: left out
            t1(1 << 5, 6, pr, 0, 0, 160, 100, 0, 0),
        ],
        8,
    );
    let second = meas(
        tow + 1000,
        0,
        &[t1(0, 5, pr + 200, 1_234_567, 0, 180, 3, 4, 0)],
        1,
    );
    let d = [b"junk".to_vec(), first, second].concat();
    let r = decode(&d, "TEST");
    let t0 = f64::from(WN) * WEEK + 345_600.0;
    assert_eq!(r.obs.t, vec![t0, t0 + 1.0]);
    let g5 = Sat::new(GnssSystem::Gps, 5);
    let f1 = 1575.42e6;
    let f2 = 1227.60e6;
    let p1 = 21_000_000.123;
    assert_eq!(find(&r, g5, "C1C").values[0], p1);
    let l1 = find(&r, g5, "L1C");
    assert!((l1.values[0] - (p1 * f1 / CLIGHT + 1234.567)).abs() < 1e-6);
    assert_eq!(find(&r, g5, "D1C").values[0], -1234.5678);
    assert_eq!(find(&r, g5, "S1C").values[0], 55.0);
    // the lock time went down and the half-cycle flag is set
    assert_eq!(l1.lli.as_ref().unwrap()[0], 0);
    assert_eq!(l1.lli.as_ref().unwrap()[1], 3);
    // L2C from the Type2 offsets
    let p2 = p1 - 70.0;
    assert!((find(&r, g5, "C2L").values[0] - p2).abs() < 1e-9);
    let l2 = find(&r, g5, "L2L").values[0];
    assert!((l2 - (p2 * f2 / CLIGHT - 2000.0)).abs() < 1e-6);
    let d2 = find(&r, g5, "D2L").values[0];
    assert!((d2 - (-1234.5678 * f2 / f1 - 0.0005)).abs() < 1e-9);
    assert_eq!(find(&r, g5, "S2L").values[0], 50.0);
    // GLONASS: slot 3, channel +1
    let r3 = Sat::new(GnssSystem::Glonass, 3);
    assert_eq!(r.glonass_fcn[&r3], 1);
    let fr = 1602.0e6 + 0.5625e6;
    let lr = find(&r, r3, "L1C").values[0];
    assert!((lr - (20_000_000.0 * fr / CLIGHT + 1.0)).abs() < 1e-6);
    assert!(find(&r, Sat::new(GnssSystem::Galileo, 11), "C1C").values[0] > 0.0);
    assert!(find(&r, Sat::new(GnssSystem::BeiDou, 20), "C2I").values[0] > 0.0);
    let j1 = Sat::new(GnssSystem::Qzss, 1);
    assert!(find(&r, j1, "C1L").values[0] > 0.0);
    assert!(!r.obs.sats.contains_key(&Sat::new(GnssSystem::Gps, 6)));
    assert!(
        r.warnings
            .iter()
            .any(|w| w.starts_with("1 SBF measurements"))
    );
    assert_eq!(r.blocks[&MEAS_EPOCH], 2);
}

#[test]
fn scrambled_epochs_are_left_out() {
    let d = meas(
        1000,
        0x80,
        &[t1(0, 5, 21_000_000_000, 0, 0, 180, 5, 0, 0)],
        1,
    );
    let r = decode(&d, "X");
    assert_eq!(r.obs.t.len(), 0);
    assert!(r.warnings.iter().any(|w| w.contains("scrambled")));
}

#[test]
fn bad_crc_is_counted_and_skipped() {
    let good = meas(1000, 0, &[t1(0, 5, 21_000_000_000, 0, 0, 180, 5, 0, 0)], 1);
    let mut bad = meas(2000, 0, &[t1(0, 5, 21_000_000_000, 0, 0, 180, 5, 0, 0)], 1);
    bad[20] ^= 0xFF;
    let d = [good.clone(), bad, good.clone()].concat();
    assert!(recognises(&d));
    let r = decode(&d, "X");
    assert_eq!(r.bad_crc, 1);
    assert_eq!(r.obs.t.len(), 1);
    assert!(
        r.warnings
            .iter()
            .any(|w| w.starts_with("1 SBF blocks failed"))
    );
    assert!(!recognises(&good));
}

/// The Kepler part shared by GPSNav, QZSNav and BDSNav (offsets 36 .. 136).
fn kepler(b: &mut Body) {
    b.f4(36, 0.0)
        .f4(40, 1e-12)
        .f4(44, 1e-4)
        .f4(48, 10.0)
        .f4(52, 1e-9)
        .f8(56, 0.25)
        .f4(64, 1e-6)
        .f8(68, 0.01)
        .f4(76, 2e-6)
        .f8(80, 5153.7)
        .u4(88, 352_800)
        .f4(92, 3e-8)
        .f8(96, -0.5)
        .f4(104, 4e-8)
        .f8(108, 0.3)
        .f4(116, 200.0)
        .f8(120, 0.75)
        .f4(128, -2.5e-9)
        .f4(132, 1e-10);
}

#[test]
fn gps_and_qzss_nav() {
    let mut b = Body::default();
    b.put(14, &[5]);
    kepler(&mut b);
    b.u2(16, WN % 1024)
        .put(18, &[1, 0, 0, 0])
        .u2(22, 77)
        .put(24, &[77, 77, 0])
        .f4(28, -1e-8)
        .u4(32, 352_800)
        .u2(136, WN % 1024)
        .u2(138, WN % 1024);
    let g = block(GPS_NAV, 0, 345_000_000, WN, &b.0);
    let mut q = Body(b.0.clone());
    q.put(14, &[183]).put(26, &[1]);
    let j = block(QZS_NAV, 0, 345_000_000, WN, &q.0);
    let r = decode(&[g.clone(), g, j].concat(), "X");
    assert_eq!(r.nav.len(), 2);
    let e = r.nav.iter().find(|e| e.rec.sat == "G05").unwrap();
    assert_eq!(e.rec.epoch, gps_time(i64::from(WN), 352_800));
    assert_eq!(e.rec.vals.len(), 29);
    assert_eq!(e.rec.vals[0], Some(f64::from(1e-4f32)));
    assert_eq!(e.rec.vals[6], Some(0.25 * GPS_PI));
    assert_eq!(e.rec.vals[21], Some(f64::from(WN)));
    assert_eq!(e.rec.vals[26], Some(77.0));
    assert_eq!(e.rec.vals[27], Some(345_000.0));
    assert_eq!(e.rec.vals[28], Some(4.0));
    let j = r.nav.iter().find(|e| e.rec.sat == "J03").unwrap();
    assert_eq!(j.rec.vals[28], Some(1.0));
}

#[test]
fn glonass_nav() {
    let mut b = Body::default();
    b.put(14, &[40, 9])
        .f8(16, 10_000.0)
        .f8(24, -20_000.0)
        .f8(32, 5_000.0)
        .f4(40, 1.5)
        .f4(44, -2.5)
        .f4(48, 0.5)
        .f4(52, 1e-9)
        .f4(64, 1e-12)
        .f4(68, 2e-5)
        .u4(76, 345_618)
        .u2(80, WN % 1024)
        .put(84, &[1, 0])
        .u2(86, 1000)
        .len(96);
    // ReceiverTime: GPS - UTC = 18 s
    let mut t = Body::default();
    t.put(14, &[26, 1, 2, 3, 4, 5, 18, 0]);
    let d = [
        block(GLO_NAV, 1, 345_000_000, WN, &b.0),
        block(RX_TIME, 0, 345_000_000, WN, &t.0),
    ]
    .concat();
    let r = decode(&d, "X");
    assert_eq!(r.nav.len(), 1);
    let e = &r.nav[0];
    assert_eq!(e.rec.sat, "R03");
    assert_eq!(e.rec.epoch, gps_time(i64::from(WN), 345_600));
    assert_eq!(e.rec.vals[0], Some(-f64::from(2e-5f32)));
    // received 345 000 s GPS = 344 982 s UTC, in the frame that starts at 344 970 s
    assert_eq!(e.rec.vals[2], Some(344_970.0));
    assert_eq!(e.rec.vals[3], Some(10_000.0));
    assert_eq!(e.rec.vals[10], Some(1.0));
    assert_eq!(e.rec.vals[14], Some(1.0));
    assert_eq!(r.nav_header.leap, Some(18));
}

#[test]
fn galileo_nav() {
    let mut b = Body::default();
    b.put(14, &[81, 2])
        .f8(16, 5440.6)
        .f8(24, 0.1)
        .f8(32, 0.0002)
        .f8(40, 0.31)
        .f8(48, -0.2)
        .f8(56, 0.9)
        .u4(100, 352_800)
        .u4(104, 352_800)
        .f8(116, -3e-4)
        .u2(124, WN % 4096)
        .u2(126, WN % 4096)
        .u2(128, 91)
        // E1-B and E5b valid; E5b signal-health 1
        .u2(130, 0x0001 | 0x0010 | (1 << 6))
        .put(133, &[255, 107])
        .f4(136, 1e-9)
        .f4(140, F_DNU as f32)
        .len(152);
    let r = decode(&block(GAL_NAV, 0, 345_000_000, WN, &b.0), "X");
    let e = &r.nav[0];
    assert_eq!(e.rec.sat, "E11");
    assert_eq!(e.rec.vals.len(), 28);
    assert_eq!(e.rec.vals[0], Some(-3e-4));
    assert_eq!(e.rec.vals[3], Some(91.0));
    assert_eq!(e.rec.vals[20], Some(517.0));
    assert_eq!(e.rec.vals[21], Some(f64::from(WN)));
    assert!((e.rec.vals[23].unwrap() - (2.0 + 7.0 * 0.16)).abs() < 1e-12);
    assert_eq!(e.rec.vals[24], Some(f64::from(1u16 << 7)));
    assert_eq!(e.rec.vals[26], None);
    assert_eq!(
        gal_health(0x0100 | (1 << 9) | (2 << 10)),
        (1 << 3) | (2 << 4)
    );
}

#[test]
fn beidou_nav_and_iono() {
    let mut b = Body::default();
    b.put(14, &[160]);
    kepler(&mut b);
    let bweek = i64::from(WN) - BDS_WEEK0;
    b.u2(16, bweek as u16)
        .put(18, &[2, 0, 5, 6])
        .f4(24, 1e-9)
        .f4(28, F_DNU as f32)
        .u4(32, 352_800)
        .u2(136, bweek as u16)
        .u2(138, bweek as u16);
    let mut ion = Body::default();
    ion.put(14, &[160]).f4(16, 1e-8).f4(32, 9e4).len(48);
    let d = [
        block(BDS_NAV, 0, 345_014_000, WN, &b.0),
        block(BDS_ION, 0, 345_014_000, WN, &ion.0),
    ]
    .concat();
    let r = decode(&d, "X");
    let e = &r.nav[0];
    assert_eq!(e.rec.sat, "C20");
    assert_eq!(e.rec.epoch, gps_time(i64::from(WN), 352_800));
    assert_eq!(e.rec.vals[3], Some(6.0));
    assert_eq!(e.rec.vals[21], Some(bweek as f64));
    assert_eq!(e.rec.vals[23], Some(4.85));
    assert_eq!(e.rec.vals[26], None);
    // reception 345 014 s GPST is 345 000 s BDT
    assert_eq!(e.rec.vals[27], Some(345_000.0));
    assert_eq!(e.rec.vals[28], Some(5.0));
    let (a, g) = r.nav_header.bds_iono.unwrap();
    assert_eq!(a[0], f64::from(1e-8f32));
    assert_eq!(g[0], f64::from(9e4f32));
}

fn setup_block() -> Vec<u8> {
    let mut b = Body::default();
    let put_text = |b: &mut Body, off: usize, s: &str| {
        b.put(off, s.as_bytes());
    };
    put_text(&mut b, 16, "SITE");
    put_text(&mut b, 156, "1234567");
    put_text(&mut b, 176, "ENGINE");
    put_text(&mut b, 196, "4.14.4");
    put_text(&mut b, 236, "SEPCHOKE_B3E6   SPKE");
    put_text(&mut b, 328, "mosaic-X5");
    b.f4(256, 1.5).f4(260, 0.0).f4(264, 0.0).len(368);
    block(RX_SETUP, 3, 1000, WN, &b.0)
}

#[test]
fn setup_and_position() {
    let mut c = Body::default();
    c.put(14, &[1, 0])
        .f8(16, 4_000_000.0)
        .f8(24, 1_000_000.0)
        .f8(32, 4_800_000.0)
        .len(88);
    let mut none = Body::default();
    none.put(14, &[0, 1]).len(88);
    let d = [
        setup_block(),
        block(PVT_CART, 2, 1000, WN, &c.0),
        block(PVT_CART, 2, 2000, WN, &none.0),
        meas(1000, 0, &[t1(0, 5, 21_000_000_000, 0, 0, 180, 5, 0, 0)], 1),
    ]
    .concat();
    let r = decode(&d, "FILE");
    assert_eq!(r.obs.name, "SITE");
    assert_eq!(r.obs.antenna, "SEPCHOKE_B3E6   SPKE");
    assert_eq!(r.obs.receiver, "SEPT MOSAIC-X5");
    assert_eq!(r.obs.delta_hen, [1.5, 0.0, 0.0]);
    assert_eq!(r.obs.approx, [4_000_000.0, 1_000_000.0, 4_800_000.0]);
    let s = r.setup.unwrap();
    assert_eq!(s.rx_serial, "1234567");
    assert_eq!(s.rx_version, "4.14.4");
}

#[test]
fn geodetic_position() {
    let p = geodetic_to_ecef(0.0, 0.0, 0.0);
    assert_eq!(p, [6_378_137.0, 0.0, 0.0]);
    let q = geodetic_to_ecef(std::f64::consts::FRAC_PI_2, 0.0, 0.0);
    assert!((q[2] - 6_356_752.314_245).abs() < 1e-3);
}

#[test]
fn convert_writes_obs_and_nav() {
    let dir = std::env::temp_dir().join(format!("sbf-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut b = Body::default();
    b.put(14, &[5]);
    kepler(&mut b);
    b.u2(16, WN % 1024)
        .u4(32, 352_800)
        .u2(136, WN % 1024)
        .u2(138, WN % 1024);
    let d = [
        setup_block(),
        meas(
            345_600_000,
            0,
            &[t1(0, 5, 21_000_000_000, 0, 0, 180, 5, 0, 0)],
            1,
        ),
        block(GPS_NAV, 0, 345_000_000, WN, &b.0),
    ]
    .concat();
    let raw = dir.join("site1.sbf");
    std::fs::write(&raw, d).unwrap();
    let c = convert(&raw, &dir, &mut |_| {}).unwrap();
    let obs = std::fs::read_to_string(&c.obs).unwrap();
    assert!(obs.contains("SEPT MOSAIC-X5"));
    assert!(obs.contains("SEPCHOKE_B3E6   SPKE"));
    assert_eq!(c.nav.len(), 1);
    let nav = std::fs::read_to_string(&c.nav[0]).unwrap();
    assert!(nav.contains("G05 2024"), "{nav}");
    assert_eq!(c.info["receiver_serial"], "1234567");
    let _ = std::fs::remove_dir_all(&dir);
}
