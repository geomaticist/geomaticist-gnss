//! `gnss`: the command line (see [`geomaticist_gnss_cli`]).

use std::process::ExitCode;

use clap::Parser;

fn main() -> ExitCode {
    let cli = geomaticist_gnss_cli::Cli::parse();
    match geomaticist_gnss_cli::run(&cli) {
        Ok(out) => {
            println!("{out}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
