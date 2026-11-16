//! Broadcast orbits and clocks (GPS, Galileo, BeiDou, QZSS Keplerian; GLONASS by RK4) and the composite
//! `Orbits` over broadcast and precise sources.
//!
//! Records are selected by signed age. The precise part is `geomaticist_gnss_products::Precise`.
//!
//! The Kepler iteration runs over the whole time array with the shared convergence test, and the GLONASS RK4 keeps its 60 s steps from tb to each time. The RK4's chain
//! of whole steps from one record is the same for every time it serves, so it is computed once per
//! (satellite, record, direction) and cached; only the last, partial step is per time, so the orbit
//! numbers stay identical.

use std::collections::HashMap;
use std::sync::Mutex;

use geomaticist_gnss_core::consts::CLIGHT;
use geomaticist_gnss_core::{GnssSystem, KeplerEph, NavData, Sat};
use geomaticist_gnss_products::Precise;

/// Gravitational constant per system.
#[must_use]
pub fn mu(sys: GnssSystem) -> Option<f64> {
    match sys {
        GnssSystem::Gps | GnssSystem::Qzss => Some(3.986_005_0e14),
        GnssSystem::Galileo | GnssSystem::BeiDou => Some(3.986_004_418e14),
        _ => None,
    }
}

/// Earth rotation rate per system.
#[must_use]
pub fn omg(sys: GnssSystem) -> f64 {
    match sys {
        GnssSystem::BeiDou => 7.292_115e-5,
        _ => 7.292_115_146_7e-5,
    }
}

/// How old a used ephemeris may be (RTKLIB's MAXDTOE per system).
#[must_use]
pub fn max_age(sys: GnssSystem) -> f64 {
    match sys {
        GnssSystem::Galileo => 14400.0,
        GnssSystem::BeiDou => 21600.0,
        GnssSystem::Glonass => 1800.0,
        _ => 7200.0,
    }
}

/// How far ahead of the time its toe may lie (B3).
#[must_use]
pub fn fwd_age(sys: GnssSystem) -> f64 {
    match sys {
        GnssSystem::Galileo | GnssSystem::BeiDou => 600.0,
        GnssSystem::Glonass => 1800.0,
        _ => 7200.0,
    }
}

const MU_GLO: f64 = 3.986_004_4e14;
const J2_GLO: f64 = 1.082_625_7e-3;
const RE_GLO: f64 = 6_378_136.0;
const OMGE_GLO: f64 = 7.292_115e-5;
/// GLONASS RK4 step (s).
pub const GLO_STEP: f64 = 60.0;
/// GPS / QZSS / BeiDou URA index 15.
pub const URA_MAX: f64 = 6144.0;
/// Galileo SISA NAPA bound.
pub const SISA_MAX: f64 = 6.0;

/// Whether the broadcast accuracy field says the ephemeris is unusable.
#[must_use]
pub fn accuracy_unusable(sys: GnssSystem, acc: f64) -> bool {
    match sys {
        GnssSystem::Galileo => acc < 0.0 || acc > SISA_MAX,
        GnssSystem::Gps | GnssSystem::Qzss | GnssSystem::BeiDou => acc > URA_MAX,
        _ => false,
    }
}

/// The result of `_select` for one time.
#[derive(Debug, Clone, Copy)]
struct Pick {
    idx: usize,
    ok: bool,
    age: f64,
}

/// The record for each time: the usable one nearest in toe among those whose signed age `t - toe` lies
/// within `[-FWD_AGE, MAX_AGE]`; where none, `ok` false and the nearest usable record indexed.
fn select(sys: GnssSystem, toe: &[f64], good: &[bool], t: &[f64]) -> Vec<Pick> {
    let cand: Vec<usize> = (0..toe.len()).filter(|&i| good[i]).collect();
    if cand.is_empty() {
        return t
            .iter()
            .map(|_| Pick {
                idx: 0,
                ok: false,
                age: f64::NAN,
            })
            .collect();
    }
    let (mx, fw) = (max_age(sys), fwd_age(sys));
    t.iter()
        .map(|&tt| {
            let mut best_d = f64::INFINITY;
            let mut best_j = 0usize;
            let mut best_abs = f64::INFINITY;
            let mut best_abs_j = 0usize;
            for (j, &ci) in cand.iter().enumerate() {
                let a = tt - toe[ci];
                let d = if a <= mx && a >= -fw {
                    a.abs()
                } else {
                    f64::INFINITY
                };
                if d < best_d {
                    best_d = d;
                    best_j = j;
                }
                if a.abs() < best_abs {
                    best_abs = a.abs();
                    best_abs_j = j;
                }
            }
            let ok = best_d.is_finite();
            let j = if ok { best_j } else { best_abs_j };
            Pick {
                idx: cand[j],
                ok,
                age: tt - toe[cand[j]],
            }
        })
        .collect()
}

/// GLONASS state derivative.
#[inline]
fn deq(x: &[f64; 6], acc: &[f64; 3]) -> [f64; 6] {
    let r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    let r3 = r2 * r2.sqrt();
    let a = 1.5 * J2_GLO * MU_GLO * RE_GLO * RE_GLO / r2 / r3;
    let b = 5.0 * x[2] * x[2] / r2;
    let c = -MU_GLO / r3 - a * (1.0 - b);
    let w2 = OMGE_GLO * OMGE_GLO;
    [
        x[3],
        x[4],
        x[5],
        (c + w2) * x[0] + 2.0 * OMGE_GLO * x[4] + acc[0],
        (c + w2) * x[1] - 2.0 * OMGE_GLO * x[3] + acc[1],
        (c - 2.0 * a) * x[2] + acc[2],
    ]
}

/// One RK4 step of `h` seconds.
#[inline]
fn rk4_step(x: &[f64; 6], acc: &[f64; 3], h: f64) -> [f64; 6] {
    let k1 = deq(x, acc);
    let mut y = [0.0; 6];
    for i in 0..6 {
        y[i] = x[i] + k1[i] * h / 2.0;
    }
    let k2 = deq(&y, acc);
    for i in 0..6 {
        y[i] = x[i] + k2[i] * h / 2.0;
    }
    let k3 = deq(&y, acc);
    for i in 0..6 {
        y[i] = x[i] + k3[i] * h;
    }
    let k4 = deq(&y, acc);
    let mut out = [0.0; 6];
    for i in 0..6 {
        out[i] = x[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * h / 6.0;
    }
    out
}

/// The cached chain of whole ±60 s steps from one GLONASS record.
#[derive(Debug, Default, Clone)]
struct GloChain {
    fwd: Vec<[f64; 6]>,
    bwd: Vec<[f64; 6]>,
}

/// Broadcast orbits and clocks from `read_nav` data.
#[derive(Debug)]
pub struct Broadcast {
    /// The navigation data.
    pub nav: NavData,
    glo_cache: Mutex<HashMap<(Sat, usize), GloChain>>,
}

impl Broadcast {
    /// Wraps navigation data.
    #[must_use]
    pub fn new(nav: NavData) -> Self {
        Self {
            nav,
            glo_cache: Mutex::new(HashMap::new()),
        }
    }

    /// Whether the satellite has any record.
    #[must_use]
    pub fn has(&self, sat: Sat) -> bool {
        self.nav.kep.contains_key(&sat) || self.nav.glo.contains_key(&sat)
    }

    /// Healthy records whose accuracy does not say unusable.
    fn usable_kep(sat: Sat, recs: &[KeplerEph]) -> Vec<bool> {
        recs.iter()
            .map(|r| r.svh == 0.0 && !accuracy_unusable(sat.sys, r.acc))
            .collect()
    }

    fn records(&self, sat: Sat) -> Option<(Vec<f64>, Vec<bool>)> {
        if sat.sys == GnssSystem::Glonass {
            let g = self.nav.glo.get(&sat)?;
            Some((
                g.iter().map(|r| r.toe).collect(),
                g.iter().map(|r| r.svh == 0.0).collect(),
            ))
        } else {
            let k = self.nav.kep.get(&sat)?;
            Some((k.iter().map(|r| r.toe).collect(), Self::usable_kep(sat, k)))
        }
    }

    /// The signed age `t - toe` (s) of the ephemeris `position` uses at each time, NaN where none is within
    /// the caps.
    #[must_use]
    pub fn age_used(&self, sat: Sat, t: &[f64]) -> Vec<f64> {
        let Some((toe, good)) = self.records(sat) else {
            return vec![f64::NAN; t.len()];
        };
        select(sat.sys, &toe, &good, t)
            .into_iter()
            .map(|p| if p.ok { p.age } else { f64::NAN })
            .collect()
    }

    /// `|age|` (s) of the usable ephemeris used, infinite where none.
    #[must_use]
    pub fn age(&self, sat: Sat, t: &[f64]) -> Vec<f64> {
        self.age_used(sat, t)
            .into_iter()
            .map(|a| {
                if a.is_finite() {
                    a.abs()
                } else {
                    f64::INFINITY
                }
            })
            .collect()
    }

    /// Seconds from each time to the nearest usable ephemeris, either way (infinite without one).
    #[must_use]
    pub fn nearest(&self, sat: Sat, t: &[f64]) -> Vec<f64> {
        let Some((toe, good)) = self.records(sat) else {
            return vec![f64::INFINITY; t.len()];
        };
        let toes: Vec<f64> = toe
            .iter()
            .zip(&good)
            .filter(|(_, g)| **g)
            .map(|(x, _)| *x)
            .collect();
        if toes.is_empty() {
            return vec![f64::INFINITY; t.len()];
        }
        t.iter()
            .map(|&tt| {
                toes.iter()
                    .map(|x| (tt - x).abs())
                    .fold(f64::INFINITY, f64::min)
            })
            .collect()
    }

    /// Healthy ephemerides within reach of `t0 … t1` left out because their accuracy says unusable:
    /// (how many, how many usable remain, the worst accuracy).
    #[must_use]
    pub fn accuracy_out(&self, sat: Sat, t0: f64, t1: f64) -> (usize, usize, f64) {
        let Some(k) = self.nav.kep.get(&sat) else {
            return (0, 0, 0.0);
        };
        let mx = self.max_age(sat);
        let (mut n_bad, mut n_ok, mut worst) = (0usize, 0usize, 0.0f64);
        for r in k {
            let near = r.toe >= t0 - mx && r.toe <= t1 + mx && r.svh == 0.0;
            if !near {
                continue;
            }
            if accuracy_unusable(sat.sys, r.acc) {
                n_bad += 1;
                worst = worst.max(r.acc.abs());
            } else {
                n_ok += 1;
            }
        }
        (n_bad, n_ok, worst)
    }

    /// `MAX_AGE` of the satellite's system.
    #[must_use]
    pub fn max_age(&self, sat: Sat) -> f64 {
        max_age(sat.sys)
    }

    /// Position (m, ECEF), clock (s) and usability at each time.
    #[must_use]
    pub fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        if sat.sys == GnssSystem::Glonass {
            return self.glo(sat, t);
        }
        if self.nav.kep.contains_key(&sat) {
            return self.kepler(sat, t);
        }
        (
            vec![[0.0; 3]; t.len()],
            vec![0.0; t.len()],
            vec![false; t.len()],
        )
    }

    fn kepler(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        let k = &self.nav.kep[&sat];
        let toe: Vec<f64> = k.iter().map(|r| r.toe).collect();
        let good = Self::usable_kep(sat, k);
        let picks = select(sat.sys, &toe, &good, t);
        let n = t.len();
        let mu = mu(sat.sys).unwrap_or(3.986_005_0e14);
        let omge = omg(sat.sys);
        let mut a_ = vec![0.0; n];
        let mut tk = vec![0.0; n];
        let mut m = vec![0.0; n];
        let mut e_ = vec![0.0; n];
        for i in 0..n {
            let g = &k[picks[i].idx];
            tk[i] = t[i] - g.toe;
            a_[i] = g.sqrta * g.sqrta;
            let nn = (mu / (a_[i] * a_[i] * a_[i])).sqrt() + g.deln;
            m[i] = g.m0 + nn * tk[i];
            e_[i] = m[i];
        }
        // Newton on the whole array with a shared convergence test
        for _ in 0..30 {
            let mut mx = 0.0f64;
            for i in 0..n {
                let g = &k[picks[i].idx];
                let d = (m[i] - e_[i] + g.e * e_[i].sin()) / (1.0 - g.e * e_[i].cos());
                e_[i] += d;
                mx = mx.max(d.abs());
            }
            if mx < 1e-14 {
                break;
            }
        }
        let geo = sat.sys == GnssSystem::BeiDou && (sat.prn <= 5 || sat.prn >= 59);
        let (c5, s5) = ((-5.0f64).to_radians().cos(), (-5.0f64).to_radians().sin());
        let mut xyz = Vec::with_capacity(n);
        let mut clk = Vec::with_capacity(n);
        for i in 0..n {
            let g = &k[picks[i].idx];
            let (sin_e, cos_e) = (e_[i].sin(), e_[i].cos());
            let e = g.e;
            let mut u = ((1.0 - e * e).sqrt() * sin_e).atan2(cos_e - e) + g.omg;
            let mut r = a_[i] * (1.0 - e * cos_e);
            let mut inc = g.i0 + g.idot * tk[i];
            let (s2, c2) = ((2.0 * u).sin(), (2.0 * u).cos());
            u = u + g.cus * s2 + g.cuc * c2;
            r = r + g.crs * s2 + g.crc * c2;
            inc = inc + g.cis * s2 + g.cic * c2;
            let (x, y) = (r * u.cos(), r * u.sin());
            let cosi = inc.cos();
            let p = if geo {
                let o = g.omg0 + g.omgd * tk[i] - omge * g.toes;
                let (sin_o, cos_o) = (o.sin(), o.cos());
                let xg = x * cos_o - y * cosi * sin_o;
                let yg = x * sin_o + y * cosi * cos_o;
                let zg = y * inc.sin();
                let (sino, coso) = ((omge * tk[i]).sin(), (omge * tk[i]).cos());
                [
                    xg * coso + yg * sino * c5 + zg * sino * s5,
                    -xg * sino + yg * coso * c5 + zg * coso * s5,
                    -yg * s5 + zg * c5,
                ]
            } else {
                let o = g.omg0 + (g.omgd - omge) * tk[i] - omge * g.toes;
                let (sin_o, cos_o) = (o.sin(), o.cos());
                [
                    x * cos_o - y * cosi * sin_o,
                    x * sin_o + y * cosi * cos_o,
                    y * inc.sin(),
                ]
            };
            xyz.push(p);
            let dt = t[i] - g.toc;
            let mut c = g.f0 + g.f1 * dt + g.f2 * dt * dt;
            c -= 2.0 * (mu * a_[i]).sqrt() * e * sin_e / (CLIGHT * CLIGHT);
            clk.push(c);
        }
        (xyz, clk, picks.iter().map(|p| p.ok).collect())
    }

    /// The state after `n` whole steps of `sign * 60 s` from record `idx`, from the cache.
    fn glo_chain_state(&self, sat: Sat, idx: usize, n: usize, forward: bool) -> [f64; 6] {
        let g = &self.nav.glo[&sat][idx];
        let mut cache = self
            .glo_cache
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let chain = cache.entry((sat, idx)).or_default();
        let list = if forward {
            &mut chain.fwd
        } else {
            &mut chain.bwd
        };
        if list.is_empty() {
            list.push([g.pos[0], g.pos[1], g.pos[2], g.vel[0], g.vel[1], g.vel[2]]);
        }
        let h = if forward { GLO_STEP } else { -GLO_STEP };
        while list.len() <= n {
            let last = list[list.len() - 1];
            list.push(rk4_step(&last, &g.acc, h));
        }
        list[n]
    }

    fn glo(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        let n = t.len();
        let Some(g) = self.nav.glo.get(&sat) else {
            return (vec![[0.0; 3]; n], vec![0.0; n], vec![false; n]);
        };
        let toe: Vec<f64> = g.iter().map(|r| r.toe).collect();
        let good: Vec<bool> = g.iter().map(|r| r.svh == 0.0).collect();
        let picks = select(GnssSystem::Glonass, &toe, &good, t);
        let mut xyz = Vec::with_capacity(n);
        let mut clk = Vec::with_capacity(n);
        for i in 0..n {
            let p = picks[i];
            let rec = &g[p.idx];
            let dt = t[i] - rec.toe;
            // the whole steps: rem -= h for each, exactly as the vectorised loop subtracts
            let mut rem = dt;
            let mut n_full = 0usize;
            let forward = rem > 0.0;
            let h = if forward { GLO_STEP } else { -GLO_STEP };
            while rem.abs() > 1e-9 && rem.abs() >= GLO_STEP {
                rem -= h;
                n_full += 1;
            }
            let mut x = if n_full == 0 {
                [
                    rec.pos[0], rec.pos[1], rec.pos[2], rec.vel[0], rec.vel[1], rec.vel[2],
                ]
            } else {
                self.glo_chain_state(sat, p.idx, n_full, forward)
            };
            if rem.abs() > 1e-9 {
                let hp = rem.signum() * rem.abs().min(GLO_STEP);
                x = rk4_step(&x, &rec.acc, hp);
            }
            xyz.push([x[0], x[1], x[2]]);
            clk.push(-rec.taun + rec.gamn * dt);
        }
        (xyz, clk, picks.iter().map(|p| p.ok).collect())
    }
}

/// An orbit source supplied by the caller: any object with `has` / `position` /
/// `center_of_mass`, e.g. synthetic constellations in tests, reached through the
/// Python bindings (`geomaticist-gnss-py`). Like a precise orbit it has no broadcast ages.
pub trait ExternalOrbit: Send + Sync + std::fmt::Debug {
    /// Whether the source has the satellite.
    fn has(&self, sat: Sat) -> bool;
    /// Whether positions are the centre of mass.
    fn center_of_mass(&self) -> bool;
    /// Position, clock, ok at the times.
    fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>);
}

/// One orbit source: anything with `has` / `position` (`Broadcast`, `Precise`).
#[derive(Debug)]
pub enum OrbitSource {
    /// Broadcast ephemerides (antenna phase centre).
    Broadcast(Broadcast),
    /// SP3 (+ CLK): centre of mass.
    Precise(Precise),
    /// A caller-supplied source ([`ExternalOrbit`]; the Python bridge's duck-typed orbits).
    External(Box<dyn ExternalOrbit>),
}

impl OrbitSource {
    /// Whether the source has the satellite.
    #[must_use]
    pub fn has(&self, sat: Sat) -> bool {
        match self {
            OrbitSource::Broadcast(b) => b.has(sat),
            OrbitSource::Precise(p) => p.has(sat),
            OrbitSource::External(e) => e.has(sat),
        }
    }

    /// Whether positions are the centre of mass.
    #[must_use]
    pub fn center_of_mass(&self) -> bool {
        match self {
            OrbitSource::Broadcast(_) => false,
            OrbitSource::Precise(_) => true,
            OrbitSource::External(e) => e.center_of_mass(),
        }
    }

    /// Position, clock, ok.
    #[must_use]
    pub fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
        match self {
            OrbitSource::Broadcast(b) => b.position(sat, t),
            OrbitSource::Precise(p) => p.position(sat, t),
            OrbitSource::External(e) => e.position(sat, t),
        }
    }
}

/// Several sources, the first that has a satellite wins.
#[derive(Debug, Default)]
pub struct Orbits {
    /// The sources in priority order.
    pub sources: Vec<OrbitSource>,
}

impl Orbits {
    /// Builds the composite.
    #[must_use]
    pub fn new(sources: Vec<OrbitSource>) -> Self {
        Self { sources }
    }

    /// The source that has the satellite.
    #[must_use]
    pub fn source(&self, sat: Sat) -> Option<&OrbitSource> {
        self.sources.iter().find(|s| s.has(sat))
    }

    /// Whether any source has the satellite.
    #[must_use]
    pub fn has(&self, sat: Sat) -> bool {
        self.source(sat).is_some()
    }

    /// (share of the times without a usable broadcast ephemeris within the caps, how far the nearest one's
    /// toe is there); `(0, 0)` for a precise orbit.
    #[must_use]
    pub fn orbit_gap(&self, sat: Sat, t: &[f64]) -> (f64, f64) {
        match self.source(sat) {
            None => (1.0, f64::INFINITY),
            // a source without `age` (precise, or any other) -> (0, 0)
            Some(OrbitSource::Precise(_) | OrbitSource::External(_)) => (0.0, 0.0),
            Some(OrbitSource::Broadcast(b)) => {
                let a = b.age(sat, t);
                let mx = b.max_age(sat);
                let near = b.nearest(sat, t);
                let bad: Vec<bool> = a.iter().map(|&x| x > mx).collect();
                let n_bad = bad.iter().filter(|&&x| x).count();
                let frac = if a.is_empty() {
                    0.0
                } else {
                    n_bad as f64 / a.len() as f64
                };
                let nearest = if n_bad > 0 {
                    near.iter()
                        .zip(&bad)
                        .filter(|(_, b)| **b)
                        .map(|(x, _)| *x)
                        .fold(f64::INFINITY, f64::min)
                } else {
                    0.0
                };
                (frac, nearest)
            }
        }
    }

    /// `Broadcast::age_used` of the source that has the satellite; `None` for a precise orbit or no source.
    #[must_use]
    pub fn age_used(&self, sat: Sat, t: &[f64]) -> Option<Vec<f64>> {
        match self.source(sat)? {
            OrbitSource::Broadcast(b) => Some(b.age_used(sat, t)),
            OrbitSource::Precise(_) | OrbitSource::External(_) => None,
        }
    }

    /// `Broadcast::accuracy_out` of the source that has the satellite; `(0, 0, 0)` otherwise.
    #[must_use]
    pub fn accuracy_out(&self, sat: Sat, t0: f64, t1: f64) -> (usize, usize, f64) {
        match self.source(sat) {
            Some(OrbitSource::Broadcast(b)) => b.accuracy_out(sat, t0, t1),
            _ => (0, 0, 0.0),
        }
    }

    /// Position, clock, ok and whether the positions are the centre of mass.
    #[must_use]
    pub fn position(&self, sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>, bool) {
        match self.source(sat) {
            None => (
                vec![[0.0; 3]; t.len()],
                vec![0.0; t.len()],
                vec![false; t.len()],
                false,
            ),
            Some(s) => {
                let (x, c, ok) = s.position(sat, t);
                (x, c, ok, s.center_of_mass())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selection_by_signed_age() {
        // Galileo: a toe 20 min ahead is not used when an older one within 4 h exists; a toe 5 min ahead is
        let toe = [0.0, 3600.0 * 3.0, 3600.0 * 3.0 + 1200.0];
        let good = [true, true, true];
        let p = select(GnssSystem::Galileo, &toe, &good, &[3600.0 * 3.0 - 300.0]);
        assert_eq!(p[0].idx, 1);
        assert!(p[0].ok && (p[0].age + 300.0).abs() < 1e-9);
        let p = select(GnssSystem::Galileo, &toe, &good, &[3600.0 * 8.0]);
        assert!(!p[0].ok && p[0].idx == 2);
        let p = select(GnssSystem::Gps, &toe, &good, &[-1800.0]);
        assert!(p[0].ok && p[0].idx == 0);
        assert!(accuracy_unusable(GnssSystem::Galileo, -1.0));
        assert!(!accuracy_unusable(GnssSystem::Gps, 2.0));
    }

    #[derive(Debug)]
    struct Fixed;
    impl ExternalOrbit for Fixed {
        fn has(&self, sat: Sat) -> bool {
            sat.prn == 1
        }
        fn center_of_mass(&self) -> bool {
            false
        }
        fn position(&self, _sat: Sat, t: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>, Vec<bool>) {
            (
                vec![[1.0, 2.0, 3.0]; t.len()],
                vec![0.0; t.len()],
                vec![true; t.len()],
            )
        }
    }

    #[test]
    fn external_source_behaves_like_python_duck_typing() {
        let orb = Orbits::new(vec![OrbitSource::External(Box::new(Fixed))]);
        let g01: Sat = "G01".parse().unwrap();
        let g02: Sat = "G02".parse().unwrap();
        assert!(orb.has(g01) && !orb.has(g02));
        let (x, _c, ok, com) = orb.position(g01, &[0.0, 1.0]);
        assert_eq!(x[1], [1.0, 2.0, 3.0]);
        assert!(ok.iter().all(|&o| o) && !com);
        assert_eq!(orb.orbit_gap(g01, &[0.0]), (0.0, 0.0));
        assert!(orb.age_used(g01, &[0.0]).is_none());
        assert_eq!(orb.accuracy_out(g01, 0.0, 1.0), (0, 0, 0.0));
    }

    #[test]
    fn rk4_chain_matches_step_by_step() {
        let x0 = [1.5e7, 1.0e7, 1.9e7, -1000.0, 2500.0, -800.0];
        let acc = [1e-9, -2e-9, 0.5e-9];
        let mut x = x0;
        for _ in 0..3 {
            x = rk4_step(&x, &acc, 60.0);
        }
        let g = geomaticist_gnss_core::GloEph {
            toe: 0.0,
            taun: 0.0,
            gamn: 0.0,
            pos: [x0[0], x0[1], x0[2]],
            vel: [x0[3], x0[4], x0[5]],
            acc,
            svh: 0.0,
            fcn: 1.0,
        };
        let mut nav = NavData::default();
        let r: Sat = "R03".parse().unwrap();
        nav.glo.insert(r, vec![g]);
        let b = Broadcast::new(nav);
        let (p, _, ok) = b.position(r, &[180.0, 200.0]);
        assert!(ok[0]);
        assert_eq!(p[0], [x[0], x[1], x[2]]);
        let y = rk4_step(&x, &acc, 20.0);
        assert_eq!(p[1], [y[0], y[1], y[2]]);
    }
}
