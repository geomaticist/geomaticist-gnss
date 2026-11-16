//! Read speed: `cargo run --release --example read_speed -- FILE...` times the engine reader
//! (`obs::read_obs`, whole file) and the file index (`RinexObs::open`) on each file.

use std::path::Path;
use std::time::Instant;

use geomaticist_gnss_rinex::obs::read_obs;
use geomaticist_gnss_rinex::reader::RinexObs;

fn main() -> Result<(), geomaticist_gnss_rinex::RinexError> {
    for f in std::env::args().skip(1) {
        let p = Path::new(&f);
        let t = Instant::now();
        let st = read_obs(p, "STA", None, None, "GRECJ")?;
        let a = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let obs = RinexObs::open(p)?;
        let b = t.elapsed().as_secs_f64();
        println!(
            "{:<45} read_obs {a:6.3} s ({} epochs, {} sats)  index {b:6.3} s ({} epochs)",
            p.file_name()
                .map_or(String::new(), |n| n.to_string_lossy().into_owned()),
            st.t.len(),
            st.sats.len(),
            obs.epoch_ticks.len()
        );
    }
    Ok(())
}
