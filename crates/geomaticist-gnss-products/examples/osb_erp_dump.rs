//! Dumps `read_osb` and `read_erp` of the given files as JSON (for cross-checks):
//! `cargo run --example osb_erp_dump -- OSB_FILE ERP_FILE`.

use geomaticist_gnss_products::{read_erp, read_osb};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() != 2 {
        return Err("usage: osb_erp_dump OSB_FILE ERP_FILE".into());
    }
    let osb = read_osb(&[&a[0]])?;
    let erp = read_erp(&[&a[1]])?;
    let osb: serde_json::Map<String, serde_json::Value> = osb
        .into_iter()
        .map(|((s, o), v)| (format!("{s} {o}"), serde_json::json!(v)))
        .collect();
    let out = serde_json::json!({"osb": osb, "erp": [erp.mjd, erp.xp, erp.yp]});
    println!("{out}");
    Ok(())
}
