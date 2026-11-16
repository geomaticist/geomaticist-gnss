//! Every receiver sample log is recognised, and the supported ones decode to plausible RINEX observations.
//!
//! The samples are RTKLIB's public test logs (BSD-2-Clause), kept on this machine only, in
//! `$GEOM_TEST_DATA_DIR/rcv_samples`. Without them the test prints a SKIPPED line (and appends it to
//! `$GEOM_TEST_DATA_REPORT` when that is set) and passes.

use std::io::Write as _;
use std::path::PathBuf;

use geomaticist_gnss_decoders::sniff::Sniffed;
use geomaticist_gnss_decoders::{RawKind, identify, near_from_name};

const TEST: &str = "rcv_samples_local";

fn skip(why: &str) {
    let line = format!("SKIPPED {TEST}: {why}");
    eprintln!("{line}");
    if let Some(r) = std::env::var_os("GEOM_TEST_DATA_REPORT")
        && let Ok(mut f) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(r)
    {
        let _ = writeln!(f, "{line}");
    }
}

fn dir() -> Option<PathBuf> {
    let Some(root) = std::env::var_os("GEOM_TEST_DATA_DIR") else {
        skip("GEOM_TEST_DATA_DIR is not set");
        return None;
    };
    let d = PathBuf::from(root).join("rcv_samples");
    if d.is_dir() {
        Some(d)
    } else {
        skip("no rcv_samples folder under GEOM_TEST_DATA_DIR");
        None
    }
}

#[test]
fn every_sample_is_recognised_and_the_supported_ones_decode() {
    let Some(dir) = dir() else { return };
    let want: &[(&str, Sniffed)] = &[
        ("GMSD7_20121014.rtcm3", Sniffed::Supported(RawKind::Rtcm3)),
        ("testglo.rtcm3", Sniffed::Supported(RawKind::Rtcm3)),
        ("ubx_20080526.ubx", Sniffed::Supported(RawKind::Ubx)),
        ("oemv_200911218.gps", Sniffed::Supported(RawKind::NovAtel)),
        ("javad_20110115.jps", Sniffed::Supported(RawKind::Greis)),
        ("oem3_20090410.gps", Sniffed::Unsupported("NovAtel OEM3")),
        (
            "cres_20080526.bin",
            Sniffed::Unsupported("Hemisphere (Crescent) binary"),
        ),
        (
            "ss2_20080517.log",
            Sniffed::Unsupported("NovAtel Superstar II"),
        ),
        ("gw10_20110121.sbas", Sniffed::Unsupported("Furuno GW-10")),
        ("testglo.rtcm2", Sniffed::Unsupported("RTCM 2")),
    ];
    let out = std::env::temp_dir().join(format!("geom-rcv-samples-{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let mut problems = Vec::new();
    for (name, kind) in want {
        let p = dir.join(name);
        if !p.exists() {
            problems.push(format!("{name}: missing"));
            continue;
        }
        let got = identify(&p);
        if got != *kind {
            problems.push(format!("{name}: identified as {got:?}, want {kind:?}"));
            continue;
        }
        let res = geomaticist_gnss_decoders::convert(&p, &out, &mut |_| {});
        let Sniffed::Supported(_) = kind else {
            let msg = res.err().map(|e| e.to_string()).unwrap_or_default();
            if !msg.contains("recognised as") {
                problems.push(format!("{name}: message {msg:?}"));
            }
            continue;
        };
        let c = match res {
            Ok(c) => c,
            Err(e) => {
                problems.push(format!("{name}: {e}"));
                continue;
            }
        };
        let st =
            geomaticist_gnss_rinex::obs::read_obs_kinds(&c.obs, "X", None, None, "GRECJS", "CLDS")
                .expect("our RINEX reads back");
        let epochs = st.t.len();
        let mut codes = 0usize;
        let mut sane = 0usize;
        let mut phases = 0usize;
        for sigs in st.sats.values() {
            for (code, s) in sigs {
                if code.to_string().starts_with('L') {
                    phases += s.values.iter().filter(|v| v.is_finite()).count();
                }
                if code.to_string().starts_with('C') {
                    for v in s.values.iter().filter(|v| v.is_finite()) {
                        codes += 1;
                        sane += usize::from((1.8e7..4.6e7).contains(v));
                    }
                }
            }
        }
        let when = near_from_name(&p).map_or(String::new(), |near| {
            let off = (st.t[0] - near) / 86_400.0;
            if off.abs() > 2.0 {
                problems.push(format!(
                    "{name}: first epoch {off:.1} days from the file name's date"
                ));
            }
            format!(", {off:+.2} d from the name's date")
        });
        println!(
            "{name}: {epochs} epochs, {} satellites, {codes} pseudoranges ({sane} in range), {phases} phases{when}; \
             warnings {:?}",
            st.sats.len(),
            c.warnings
        );
        if epochs < 2 || codes == 0 || phases == 0 || sane * 100 < codes * 95 {
            problems.push(format!(
                "{name}: {epochs} epochs, {sane} of {codes} pseudoranges in range"
            ));
        }
    }
    let _ = std::fs::remove_dir_all(&out);
    assert!(problems.is_empty(), "{problems:#?}");
}
