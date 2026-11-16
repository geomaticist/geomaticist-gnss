//! Simulated observations of a receiver on a known trajectory, built with the engine's own models, for
//! tests and examples: with the truth known to the millimetre, the kinematic solution and the event
//! interpolation can be checked end to end without any recorded data.
//!
//! [`SynthOrbit`] is a constellation on circular orbits (an [`ExternalOrbit`], so the solver uses it like
//! any other orbit source). [`simulate`] writes the code, phase, Doppler and signal strength a receiver would
//! record: geometric range with the Earth-rotation term, light time, satellite and receiver clocks, the
//! troposphere of the standard atmosphere, integer ambiguities, optional cycle slips, and noise from a
//! seeded generator (the same seed gives the same numbers on every machine).

use std::collections::BTreeMap;

use geomaticist_gnss_core::consts::{CLIGHT, OMGE};
use geomaticist_gnss_core::obs::{SatSignals, Signal};
use geomaticist_gnss_core::signal::{code_key, freq};
use geomaticist_gnss_core::{GnssSystem, ObsCode, Sat, StationObs};

use super::prepare::slant_tropo;
use crate::models::{GM_EARTH, azel, ecef_to_llh, enu_matrix, geometry};
use crate::orbits::{ExternalOrbit, OrbitSource, Orbits};
use crate::tropo::TropoInputs;

/// One simulated satellite on a circular orbit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynthSat {
    /// The satellite.
    pub sat: Sat,
    /// Orbit radius, m.
    pub radius: f64,
    /// Inclination, rad.
    pub inc: f64,
    /// Right ascension of the ascending node at t = 0, rad.
    pub raan: f64,
    /// Argument of latitude at t = 0, rad.
    pub u0: f64,
    /// Clock offset, s.
    pub clk: f64,
    /// Clock drift, s/s.
    pub drift: f64,
}

impl SynthSat {
    /// ECEF position at GPS seconds `t` (m).
    #[must_use]
    pub fn position(&self, t: f64) -> [f64; 3] {
        // seconds from a reference near the data keep the angles small and exact
        let n = (GM_EARTH / self.radius.powi(3)).sqrt();
        let u = self.u0 + n * t;
        let (xo, yo) = (self.radius * u.cos(), self.radius * u.sin());
        let om = self.raan - OMGE * t;
        let (ci, si) = (self.inc.cos(), self.inc.sin());
        [
            xo * om.cos() - yo * ci * om.sin(),
            xo * om.sin() + yo * ci * om.cos(),
            yo * si,
        ]
    }
}

/// A simulated constellation; times are taken relative to `t0` so that the orbit angles stay small.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthOrbit {
    /// The satellites.
    pub sats: Vec<SynthSat>,
    /// Reference time, GPS seconds.
    pub t0: f64,
}

impl SynthOrbit {
    /// `planes × per_plane` satellites for each system given, on inclined circular orbits (a Walker
    /// pattern), with small clock offsets.
    #[must_use]
    pub fn constellation(systems: &[GnssSystem], planes: u8, per_plane: u8, t0: f64) -> Self {
        let mut sats = Vec::new();
        for (si, sys) in systems.iter().enumerate() {
            let (radius, inc) = match sys {
                GnssSystem::Galileo => (29_600e3, 56f64.to_radians()),
                GnssSystem::BeiDou => (27_900e3, 55f64.to_radians()),
                GnssSystem::Glonass => (25_500e3, 64.8f64.to_radians()),
                _ => (26_560e3, 55f64.to_radians()),
            };
            for p in 0..planes {
                for k in 0..per_plane {
                    let prn = p * per_plane + k + 1;
                    let fp = f64::from(p) / f64::from(planes);
                    let fk = f64::from(k) / f64::from(per_plane);
                    sats.push(SynthSat {
                        sat: Sat::new(*sys, prn),
                        radius,
                        inc,
                        raan: std::f64::consts::TAU * fp + 0.4 * si as f64,
                        u0: std::f64::consts::TAU * (fk + 0.37 * fp) + 0.9 * si as f64,
                        clk: 1e-5 * (f64::from(prn) - 12.0) + 3e-6 * si as f64,
                        drift: 1e-11 * (f64::from(prn % 5) - 2.0),
                    });
                }
            }
        }
        Self { sats, t0 }
    }

    /// The constellation as the solver's orbit source.
    #[must_use]
    pub fn orbits(&self) -> Orbits {
        Orbits::new(vec![OrbitSource::External(Box::new(self.clone()))])
    }

    fn find(&self, sat: Sat) -> Option<&SynthSat> {
        self.sats.iter().find(|s| s.sat == sat)
    }
}

impl ExternalOrbit for SynthOrbit {
    fn has(&self, sat: Sat) -> bool {
        self.find(sat).is_some()
    }

    fn center_of_mass(&self) -> bool {
        false
    }

    fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        match self.find(sat) {
            Some(s) => (
                t.iter().map(|&x| s.position(x - self.t0)).collect(),
                t.iter().map(|&x| s.clk + s.drift * (x - self.t0)).collect(),
                vec![true; t.len()],
            ),
            None => (
                vec![[0.0; 3]; t.len()],
                vec![0.0; t.len()],
                vec![false; t.len()],
            ),
        }
    }
}

/// A small seeded generator of normal deviates (xorshift64* and Box-Muller).
#[derive(Debug, Clone)]
pub struct Noise(u64);

impl Noise {
    /// A generator from a seed.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }

    /// A uniform deviate in (0, 1).
    pub fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let v = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11;
        (v as f64 + 0.5) / (1u64 << 53) as f64
    }

    /// A standard normal deviate.
    pub fn normal(&mut self) -> f64 {
        let (a, b) = (self.uniform(), self.uniform());
        (-2.0 * a.ln()).sqrt() * (std::f64::consts::TAU * b).cos()
    }
}

/// A cycle slip to put into the simulated phases.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SynthSlip {
    /// The satellite.
    pub sat: Sat,
    /// RINEX band number.
    pub band: u8,
    /// From this epoch tag on (GPS seconds) ...
    pub t: f64,
    /// ... the phase is this many cycles higher.
    pub cycles: f64,
    /// Whether the receiver flags it (loss-of-lock indicator).
    pub flagged: bool,
}

/// What a simulated receiver records.
#[derive(Debug, Clone, PartialEq)]
pub struct SynthReceiver {
    /// Station name.
    pub name: String,
    /// RINEX 3 phase codes per system letter, e.g. `'G' -> ["L1C", "L2W"]` (code, Doppler and signal
    /// strength are written with the same band and attribute).
    pub signals: BTreeMap<char, Vec<String>>,
    /// Receiver clock offset at the first epoch, s.
    pub clock: f64,
    /// Receiver clock drift, s/s.
    pub clock_drift: f64,
    /// Phase noise at the zenith, m.
    pub phase_noise: f64,
    /// Code noise at the zenith, m.
    pub code_noise: f64,
    /// Doppler noise, Hz.
    pub doppler_noise: f64,
    /// Elevation below which nothing is recorded, degrees.
    pub cutoff: f64,
    /// Seed of the noise and of the ambiguities.
    pub seed: u64,
    /// Cycle slips.
    pub slips: Vec<SynthSlip>,
    /// Epoch tags at which the receiver records nothing (missing epochs).
    pub missing: Vec<f64>,
    /// Spans in which a phase carries the loss-of-lock flag although it is whole (receivers keep the flag
    /// up for a while): satellite, band, from, to (epoch tags).
    pub flags: Vec<(Sat, u8, f64, f64)>,
    /// Spans in which a phase is not recorded at all (code and Doppler are): satellite, band, from, to.
    pub phase_gaps: Vec<(Sat, u8, f64, f64)>,
    /// The Doppler is written with the opposite sign of the RINEX definition, as some receivers do.
    pub doppler_inverted: bool,
    /// Write Doppler.
    pub doppler: bool,
}

impl SynthReceiver {
    /// A dual-frequency GPS, single-frequency Galileo and BeiDou receiver with a few millimetres of phase
    /// noise.
    #[must_use]
    pub fn new(name: &str, seed: u64) -> Self {
        let mut signals = BTreeMap::new();
        signals.insert('G', vec!["L1C".to_owned(), "L2W".to_owned()]);
        signals.insert('E', vec!["L1C".to_owned()]);
        signals.insert('C', vec!["L2I".to_owned()]);
        Self {
            name: name.to_owned(),
            signals,
            clock: 0.0,
            clock_drift: 0.0,
            phase_noise: 0.002,
            code_noise: 0.25,
            doppler_noise: 0.05,
            cutoff: 5.0,
            seed,
            slips: Vec::new(),
            missing: Vec::new(),
            flags: Vec::new(),
            phase_gaps: Vec::new(),
            doppler_inverted: false,
            doppler: true,
        }
    }
}

/// The noise-free pseudorange (m) of `sat` at the receive time `t` for a receiver at `x`: range with
/// light time and Earth rotation, less the satellite clock, plus the troposphere. `None` below `cutoff`.
fn range_of(
    orbit: &SynthOrbit,
    s: &SynthSat,
    t: f64,
    x: [f64; 3],
    cutoff: f64,
) -> Option<(f64, f64)> {
    let mut tau = 0.075;
    let mut rho = 0.0;
    let mut e = [0.0; 3];
    for _ in 0..4 {
        let sp = s.position(t - tau - orbit.t0);
        let g = geometry(sp, x);
        rho = g.0;
        e = g.1;
        tau = rho / CLIGHT;
    }
    let (lat, lon, h) = ecef_to_llh(x);
    let (_, el, _) = azel(e, &enu_matrix(lat, lon));
    if el < cutoff {
        return None;
    }
    let trop = slant_tropo(
        "niell",
        0.5,
        &TropoInputs::default(),
        t,
        (lat, lon, h),
        &[el],
    )[0];
    let dts = s.clk + s.drift * (t - tau - orbit.t0);
    Some((rho - CLIGHT * dts + trop, el))
}

/// Simulates the observations of `rx` at the epoch tags `epochs` for an antenna whose ECEF position at GPS
/// time `t` is `position(t)`.
#[must_use]
pub fn simulate(
    orbit: &SynthOrbit,
    rx: &SynthReceiver,
    position: &dyn Fn(f64) -> [f64; 3],
    epochs: &[f64],
) -> StationObs {
    let tags: Vec<f64> = epochs
        .iter()
        .copied()
        .filter(|t| !rx.missing.iter().any(|m| (m - t).abs() < 1e-6))
        .collect();
    let n = tags.len();
    let mut st = StationObs::new(rx.name.clone(), format!("{}.obs", rx.name));
    st.t.clone_from(&tags);
    st.power_fail = vec![false; n];
    let t_first = tags.first().copied().unwrap_or(0.0);
    st.approx = position(t_first);
    let mut rng = Noise::new(rx.seed);
    let cutoff = rx.cutoff.to_radians();
    let mut order = 0u32;
    for s in &orbit.sats {
        let Some(codes) = rx.signals.get(&s.sat.sys.letter()) else {
            continue;
        };
        struct Sig {
            band: u8,
            lam: f64,
            amb: f64,
            l: Vec<f64>,
            lli: Vec<i8>,
            c: Vec<f64>,
            d: Vec<f64>,
            snr: Vec<f64>,
            code: String,
        }
        let mut bands: Vec<Sig> = Vec::new();
        for code in codes {
            let key = code_key(code.as_bytes(), s.sat.sys.letter(), 3, 3.04);
            let Ok(f) = freq(s.sat, key.band, Some(0)) else {
                continue;
            };
            bands.push(Sig {
                band: key.band,
                lam: CLIGHT / f,
                amb: (rng.uniform() * 2.0e6 - 1.0e6).round(),
                l: vec![f64::NAN; n],
                lli: vec![0; n],
                c: vec![f64::NAN; n],
                d: vec![f64::NAN; n],
                snr: vec![f64::NAN; n],
                code: code.clone(),
            });
        }
        let mut any = false;
        for (i, &tag) in tags.iter().enumerate() {
            let dt = rx.clock + rx.clock_drift * (tag - t_first);
            let t = tag - dt;
            let x = position(t);
            let Some((p, el)) = range_of(orbit, s, t, x, cutoff) else {
                continue;
            };
            any = true;
            // range rate over ±50 ms, with the receiver clock drift
            const H: f64 = 0.05;
            let rate = match (
                range_of(orbit, s, t - H, position(t - H), -1.0),
                range_of(orbit, s, t + H, position(t + H), -1.0),
            ) {
                (Some(a), Some(b)) => (b.0 - a.0) / (2.0 * H) + CLIGHT * rx.clock_drift,
                _ => f64::NAN,
            };
            let w = 1.0 / el.sin().max(0.1);
            for b in &mut bands {
                let slip: f64 = rx
                    .slips
                    .iter()
                    .filter(|sl| sl.sat == s.sat && sl.band == b.band && tag >= sl.t - 1e-6)
                    .map(|sl| sl.cycles)
                    .sum();
                let flagged = rx.slips.iter().any(|sl| {
                    sl.sat == s.sat && sl.band == b.band && sl.flagged && (tag - sl.t).abs() < 1e-6
                });
                let range = p + CLIGHT * dt;
                b.c[i] = range + rx.code_noise * w * rng.normal();
                b.l[i] = (range + rx.phase_noise * w * rng.normal()) / b.lam + b.amb + slip;
                let within = |v: &[(Sat, u8, f64, f64)]| {
                    v.iter().any(|f| {
                        f.0 == s.sat && f.1 == b.band && tag >= f.2 - 1e-6 && tag <= f.3 + 1e-6
                    })
                };
                b.lli[i] = i8::from(flagged || within(&rx.flags));
                if within(&rx.phase_gaps) {
                    b.l[i] = f64::NAN;
                }
                let sign = if rx.doppler_inverted { 1.0 } else { -1.0 };
                b.d[i] = sign * rate / b.lam + rx.doppler_noise * rng.normal();
                b.snr[i] = 30.0 + 20.0 * el.sin();
            }
        }
        if !any {
            continue;
        }
        let mut sigs = SatSignals::new();
        for b in bands {
            let tail = &b.code[1..];
            let kinds: [(char, Vec<f64>, Option<Vec<i8>>); 4] = [
                ('C', b.c, None),
                ('L', b.l, Some(b.lli)),
                ('D', b.d, None),
                ('S', b.snr, None),
            ];
            for (kind, values, lli) in kinds {
                if kind == 'D' && !rx.doppler {
                    continue;
                }
                let name = format!("{kind}{tail}");
                let Ok(code) = name.parse::<ObsCode>() else {
                    continue;
                };
                sigs.insert(
                    code,
                    Signal {
                        key: code_key(name.as_bytes(), s.sat.sys.letter(), 3, 3.04),
                        values,
                        lli,
                        first_seen: order,
                    },
                );
                order += 1;
            }
        }
        st.sats.insert(s.sat, sigs);
    }
    st
}

/// ECEF of a point given as east / north / up metres from `origin`.
#[must_use]
pub fn enu_point(origin: [f64; 3], enu: [f64; 3]) -> [f64; 3] {
    let (lat, lon, _) = ecef_to_llh(origin);
    let d = crate::linalg::mat3t_vec(&enu_matrix(lat, lon), enu);
    [origin[0] + d[0], origin[1] + d[1], origin[2] + d[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbit_and_noise_are_sane() {
        let orb = SynthOrbit::constellation(&[GnssSystem::Gps, GnssSystem::Galileo], 6, 4, 1.4e9);
        assert_eq!(orb.sats.len(), 48);
        let s = orb.sats[0];
        let (a, b) = (s.position(0.0), s.position(1.0));
        let r = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        assert!((r - 26_560e3).abs() < 1.0);
        let v = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt();
        assert!(v > 2000.0 && v < 4500.0, "{v}");
        let mut n = Noise::new(7);
        let x: Vec<f64> = (0..20_000).map(|_| n.normal()).collect();
        let m = x.iter().sum::<f64>() / x.len() as f64;
        let sd = (x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / x.len() as f64).sqrt();
        assert!(m.abs() < 0.03 && (sd - 1.0).abs() < 0.03, "{m} {sd}");
    }

    #[test]
    fn a_static_receiver_sees_enough_satellites() {
        let t0 = 1.4e9;
        let orb = SynthOrbit::constellation(
            &[GnssSystem::Gps, GnssSystem::Galileo, GnssSystem::BeiDou],
            6,
            4,
            t0,
        );
        let x = [4_395_470.6, 617_743.1, 4_565_607.2];
        let st = simulate(
            &orb,
            &SynthReceiver::new("BASE", 1),
            &|_| x,
            &[t0, t0 + 1.0],
        );
        assert!(st.sats.len() >= 15, "{}", st.sats.len());
        let g = st
            .sats
            .iter()
            .find(|(s, _)| s.sys == GnssSystem::Gps)
            .unwrap();
        assert!(g.1.contains_key(&"D2W".parse::<ObsCode>().unwrap()));
    }
}
