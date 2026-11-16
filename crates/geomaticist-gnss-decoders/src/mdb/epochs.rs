//! The observation record (type 0x77) of a Leica MDB log: one per epoch, bit-packed, most significant bit
//! first.
//!
//! ```text
//! record body: [u8 0x77][u16 GPS week][u32 ms of week][6 bytes, OPEN][u8 satellites]   then per satellite:
//!   [6 PRN][8 GLONASS channel + 45, 0 for GPS][7 OPEN][3 signals]                       then per signal:
//!   first:  [4 code][5 flags][11 SNR][8 OPEN][8 OPEN][u32 range]       [s20 phase - range][s32 Doppler]
//!   others: [4 code][5 flags][11 SNR][8 OPEN][8 OPEN][s14 range - first][s20 phase - range][s32 Doppler]
//! zero bits to the next byte, then the record's check.
//! ```
//!
//! Range 0.02 m, phase - range 0.0005 m (phase = (range + that) / wavelength), Doppler 2^-10 Hz, SNR 0.05 dB-Hz.

/// Record type of an observation record.
pub const T_EPOCH: u8 = 0x77;
/// Seconds per week.
pub const WEEK_S: i64 = 604_800;
/// Speed of light, m/s.
pub const C_LIGHT: i64 = 299_792_458;
/// The channel field is k + 45 for GLONASS.
pub const GLO_K0: i64 = 45;
/// f1: loss of lock since the last epoch (RINEX LLI 1).
pub const F_LLI: u32 = 0x08;
/// f3: phase valid.
pub const F_PHASE: u32 = 0x02;

/// Signal code -> RINEX band + attribute.
#[must_use]
pub fn rinex_code(system: char, code: u32) -> Option<&'static str> {
    match (system, code) {
        ('G' | 'R', 0) => Some("1C"),
        ('G', 6) => Some("2W"),
        ('G', 8) => Some("2S"),
        ('R', 4) => Some("2P"),
        _ => None,
    }
}

/// One signal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Signal {
    /// The 4-bit signal code.
    pub code: u32,
    /// "1C", "2W" ... or `None` for a code not known yet.
    pub rinex: Option<&'static str>,
    /// The 5 flag bits.
    pub flags: u32,
    /// SNR, 0.05 dB-Hz units.
    pub snr_raw: i64,
    /// Range, 0.02 m units.
    pub range_raw: i64,
    /// Phase - range, 0.0005 m units.
    pub phase_raw: i64,
    /// Doppler, 2^-10 Hz units.
    pub doppler_raw: i64,
}

impl Signal {
    /// Loss of lock.
    #[must_use]
    pub fn lli(&self) -> bool {
        self.flags & F_LLI != 0
    }
    /// Phase valid.
    #[must_use]
    pub fn has_phase(&self) -> bool {
        self.flags & F_PHASE != 0
    }
}

/// One satellite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sat {
    /// 'G' or 'R'.
    pub system: char,
    /// PRN / slot.
    pub prn: i64,
    /// GLONASS frequency channel k.
    pub channel: Option<i64>,
    /// Signals.
    pub signals: Vec<Signal>,
}

impl Sat {
    /// "G05".
    #[must_use]
    pub fn name(&self) -> String {
        format!("{}{:02}", self.system, self.prn)
    }
}

/// One epoch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Epoch {
    /// GPS week.
    pub week: i64,
    /// ms of week.
    pub ms: i64,
    /// Satellites.
    pub sats: Vec<Sat>,
}

impl Epoch {
    /// `week * 604800 + ms / 1000.0`.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn gps_seconds(&self) -> f64 {
        (self.week * WEEK_S) as f64 + self.ms as f64 / 1000.0
    }
}

/// Most-significant-bit-first reader.
struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    n: usize,
}

impl Bits<'_> {
    fn u(&mut self, w: usize) -> Result<u64, String> {
        let p = self.pos + w;
        if p > self.n {
            return Err("observation record ends inside a field".into());
        }
        // fields are at most 32 bits: an 8-byte window always holds one
        let (byte, off) = (self.pos / 8, self.pos % 8);
        let mut acc = 0u64;
        for k in 0..8 {
            acc = (acc << 8) | u64::from(self.data.get(byte + k).copied().unwrap_or(0));
        }
        self.pos = p;
        Ok(if w == 0 { 0 } else { (acc << off) >> (64 - w) })
    }

    #[allow(clippy::cast_possible_wrap)]
    fn s(&mut self, w: usize) -> Result<i64, String> {
        let x = self.u(w)?;
        Ok(if x >> (w - 1) != 0 {
            x as i64 - (1i64 << w)
        } else {
            x as i64
        })
    }
}

/// One 0x77 record body (from the type byte to before the check) -> [`Epoch`].
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
pub fn decode(body: &[u8]) -> Result<Epoch, String> {
    if body.is_empty() || body[0] != T_EPOCH || body.len() < 14 {
        return Err("not an observation record".into());
    }
    let week = i64::from(u16::from_le_bytes([body[1], body[2]]));
    let ms = i64::from(u32::from_le_bytes([body[3], body[4], body[5], body[6]]));
    let nsat = body[13];
    let mut b = Bits {
        data: &body[14..],
        pos: 0,
        n: (body.len() - 14) * 8,
    };
    let mut sats = Vec::with_capacity(usize::from(nsat));
    for _ in 0..nsat {
        let prn = b.u(6)? as i64;
        let ch = b.u(8)? as i64;
        let _o1 = b.u(7)?;
        let nsig = b.u(3)?;
        let system = if ch == 0 { 'G' } else { 'R' };
        let mut sat = Sat {
            system,
            prn,
            channel: (system == 'R').then_some(ch - GLO_K0),
            signals: Vec::new(),
        };
        let mut first = 0i64;
        for j in 0..nsig {
            let code = b.u(4)? as u32;
            let flags = b.u(5)? as u32;
            let snr = b.u(11)? as i64;
            let _o2 = b.u(8)?;
            let _o3 = b.u(8)?;
            let rng = if j == 0 {
                first = b.u(32)? as i64;
                first
            } else {
                first + b.s(14)?
            };
            let ph = b.s(20)?;
            let dop = b.s(32)?;
            sat.signals.push(Signal {
                code,
                rinex: rinex_code(system, code),
                flags,
                snr_raw: snr,
                range_raw: rng,
                phase_raw: ph,
                doppler_raw: dop,
            });
        }
        sats.push(sat);
    }
    if b.n - b.pos >= 8 {
        return Err(format!(
            "observation record has {} bytes left over",
            (b.n - b.pos) / 8
        ));
    }
    Ok(Epoch { week, ms, sats })
}

/// Carrier frequency, Hz (exact integers; `None` for an unknown signal).
#[must_use]
pub fn frequency(system: char, rinex: &str, channel: Option<i64>) -> Option<i64> {
    let band = rinex.as_bytes().first()?;
    match (system, band) {
        ('G', b'1') => Some(1_575_420_000),
        ('G', b'2') => Some(1_227_600_000),
        ('G', b'5') => Some(1_176_450_000),
        ('R', b'1') => Some(1_602_000_000 + 562_500 * channel?),
        ('R', b'2') => Some(1_246_000_000 + 437_500 * channel?),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bits_and_signs() {
        let mut b = Bits {
            data: &[0b1010_0000, 0xFF],
            pos: 0,
            n: 16,
        };
        assert_eq!(b.u(3).unwrap(), 0b101);
        assert_eq!(b.s(5).unwrap(), 0);
        assert_eq!(b.s(4).unwrap(), -1);
        assert!(b.u(5).is_err());
    }

    #[test]
    fn frequencies() {
        assert_eq!(
            frequency('R', "1C", Some(-7)),
            Some(1_602_000_000 - 7 * 562_500)
        );
        assert_eq!(frequency('G', "2W", None), Some(1_227_600_000));
        assert_eq!(frequency('G', "7X", None), None);
    }
}
