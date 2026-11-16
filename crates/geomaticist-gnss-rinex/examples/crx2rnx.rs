//! `crx2rnx IN OUT`: decodes a (gzipped) CRINEX file to RINEX with the crate's own decoder.

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: crx2rnx IN OUT");
        std::process::exit(2);
    }
    let t = std::time::Instant::now();
    let data = geomaticist_gnss_rinex::open::read_plain(std::path::Path::new(&args[1]))?;
    eprintln!(
        "decoded {} bytes in {:.3} s",
        data.len(),
        t.elapsed().as_secs_f64()
    );
    std::fs::write(&args[2], data)?;
    Ok(())
}
