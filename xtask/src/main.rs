//! `cargo xtask <task>` for the engine repository.
//!
//! - `check-deps`: dependency rules over `cargo metadata`: R1, R2, R8, plus the
//!   licence field and the `geomaticist-gnss-` name prefix of every crate.
//! - `ci`: `fmt --check`, `clippy -D warnings`, `test`, `check-deps`, `cargo deny check`.

mod graph;

use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::Instant;

use graph::{Edges, Graph, Report, families};

/// Engine crates that must never depend on `geomaticist-gnss-decoders` (rule R2): every member but the
/// decoders themselves, xtask and the Python bindings (which expose `decode_raw`).
const R2_EXEMPT: &[&str] = &["geomaticist-gnss-decoders", "xtask", "geomaticist-gnss-py"];

/// The Python bindings crate: the only crate allowed to depend on pyo3, and no crate may depend on it
/// (rule "py").
const PY_BRIDGE: &str = "geomaticist-gnss-py";

/// Members that are tooling, not part of the library.
const TOOLING: &[&str] = &["xtask"];

const LICENSE: &str = "MIT OR Apache-2.0";

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
}

fn check_deps(root: &Path) -> Result<(), String> {
    let g = Graph::load(root)?;
    let mut r = Report::default();

    // R1: nothing from the app repo (no path dependency outside this workspace) and no
    // GDAL / PROJ / network / GUI crate anywhere in an engine crate's build (dev-deps included).
    r.rule("R1");
    for (name, id) in &g.members {
        if TOOLING.contains(&name.as_str()) {
            continue;
        }
        for (dep, via) in g.transitive(id, Edges::All) {
            if families::gdal_proj(&dep) || families::network(&dep) || families::gui(&dep) {
                r.fail("R1", format!("{name} pulls in {dep} ({via})"));
            }
        }
    }
    for p in g.all_packages() {
        let outside = p.source.is_none() && !p.manifest_path.as_std_path().starts_with(&g.root);
        if outside {
            r.fail(
                "R1",
                format!(
                    "{} is a path dependency outside this repository ({})",
                    p.name, p.manifest_path
                ),
            );
        }
    }

    // R2: no processing crate depends on geomaticist-gnss-decoders, directly or not, in any dependency kind.
    r.rule("R2");
    for (name, id) in &g.members {
        if R2_EXEMPT.contains(&name.as_str()) {
            continue;
        }
        if let Some(via) = g
            .transitive(id, Edges::All)
            .get("geomaticist-gnss-decoders")
        {
            r.fail(
                "R2",
                format!("{name} depends on geomaticist-gnss-decoders ({via})"),
            );
        }
    }

    // py: pyo3 only in the bridge, and nothing depends on the bridge.
    r.rule("py");
    for (name, id) in &g.members {
        if name == PY_BRIDGE {
            continue;
        }
        for (dep, via) in g.transitive(id, Edges::All) {
            if dep == "pyo3" || dep.starts_with("pyo3-") || dep == "numpy" {
                r.fail(
                    "py",
                    format!("{name} pulls in {dep} ({via}): pyo3 only in {PY_BRIDGE}"),
                );
            }
            if dep == PY_BRIDGE {
                r.fail(
                    "py",
                    format!("{name} depends on the bridge {PY_BRIDGE} ({via})"),
                );
            }
        }
    }

    // R8: every crate has unit tests.
    r.rule("R8");
    for (name, id) in &g.members {
        if !g.has_tests(id) {
            r.fail("R8", format!("{name} has no #[test]"));
        }
    }

    // Licence: every library crate carries the dual licence.
    r.rule("licence");
    for (name, id) in &g.members {
        let lic = g.package(id).license.as_deref().unwrap_or("");
        if lic != LICENSE {
            r.fail(
                "licence",
                format!("{name} has license {lic:?}, expected {LICENSE:?}"),
            );
        }
    }

    // Names: every library crate carries the engine prefix.
    r.rule("names");
    for name in g.members.keys() {
        if !TOOLING.contains(&name.as_str()) && !name.starts_with("geomaticist-gnss-") {
            r.fail(
                "names",
                format!("{name} lacks the geomaticist-gnss- prefix"),
            );
        }
    }
    r.finish()
}

fn run(root: &Path, label: &str, program: &str, args: &[&str]) -> Result<(), String> {
    println!("==> {label}: {program} {}", args.join(" "));
    let t = Instant::now();
    let status = Command::new(program)
        .args(args)
        .current_dir(root)
        .status()
        .map_err(|e| format!("{label}: cannot run {program}: {e}"))?;
    println!("    {label}: {:.1} s", t.elapsed().as_secs_f64());
    if status.success() {
        Ok(())
    } else {
        Err(format!("{label} failed ({status})"))
    }
}

fn ci(root: &Path) -> Result<(), String> {
    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());
    let t = Instant::now();
    run(root, "fmt", &cargo, &["fmt", "--all", "--check"])?;
    run(
        root,
        "clippy",
        &cargo,
        &[
            "clippy",
            "--workspace",
            "--all-targets",
            "--all-features",
            "--",
            "-D",
            "warnings",
        ],
    )?;
    run(
        root,
        "test",
        &cargo,
        &["test", "--workspace", "--all-features"],
    )?;
    println!("==> check-deps");
    check_deps(root)?;
    run(root, "deny", &cargo, &["deny", "--all-features", "check"])?;
    println!("ci: all green in {:.1} s", t.elapsed().as_secs_f64());
    Ok(())
}

fn main() -> ExitCode {
    let root = workspace_root();
    let task = std::env::args().nth(1).unwrap_or_default();
    let result = match task.as_str() {
        "check-deps" => check_deps(&root),
        "ci" => ci(&root),
        _ => {
            eprintln!("usage: cargo xtask <ci | check-deps>");
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn families() {
        use super::families::{gdal_proj, gui, network, toolkit};
        assert!(gdal_proj("gdal-sys") && network("reqwest") && toolkit("egui_extras"));
        assert!(gui("wgpu") && !toolkit("wgpu") && !gui("geom-render"));
    }

    #[test]
    fn r2_rule_holds_here() {
        super::check_deps(&super::workspace_root()).unwrap();
    }
}
