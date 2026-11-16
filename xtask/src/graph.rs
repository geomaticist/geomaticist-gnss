//! Dependency graph from `cargo metadata`, shared by the rule checks.
//!
//! This file is identical in `geomaticist-gnss/xtask` and `geomaticist-rs/xtask`; change both together.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::{Path, PathBuf};

use cargo_metadata::{DependencyKind, Metadata, MetadataCommand, Package, PackageId};

/// Which edges to follow.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Edges {
    /// Normal and build dependencies (what ships).
    Shipping,
    /// Every kind, dev-dependencies included.
    All,
}

/// The resolved graph of one workspace.
pub struct Graph {
    /// Workspace root (used by the engine repo's R1 check).
    #[allow(dead_code)]
    pub root: PathBuf,
    packages: BTreeMap<PackageId, Package>,
    /// Edges: package -> [(dependency, kinds)].
    edges: BTreeMap<PackageId, Vec<(PackageId, Vec<DependencyKind>)>>,
    /// Workspace members by name.
    pub members: BTreeMap<String, PackageId>,
}

impl Graph {
    /// Runs `cargo metadata` (all features, every target platform) in `dir`.
    pub fn load(dir: &Path) -> Result<Self, String> {
        let meta: Metadata = MetadataCommand::new()
            .current_dir(dir)
            .features(cargo_metadata::CargoOpt::AllFeatures)
            .exec()
            .map_err(|e| format!("cargo metadata: {e}"))?;
        let resolve = meta
            .resolve
            .as_ref()
            .ok_or("cargo metadata: no resolve graph")?;
        let packages: BTreeMap<PackageId, Package> = meta
            .packages
            .iter()
            .map(|p| (p.id.clone(), p.clone()))
            .collect();
        let edges = resolve
            .nodes
            .iter()
            .map(|n| {
                let deps = n
                    .deps
                    .iter()
                    .map(|d| (d.pkg.clone(), d.dep_kinds.iter().map(|k| k.kind).collect()))
                    .collect();
                (n.id.clone(), deps)
            })
            .collect();
        let members = meta
            .workspace_members
            .iter()
            .map(|id| (packages[id].name.to_string(), id.clone()))
            .collect();
        Ok(Self {
            root: meta.workspace_root.into_std_path_buf(),
            packages,
            edges,
            members,
        })
    }

    /// Package by id.
    pub fn package(&self, id: &PackageId) -> &Package {
        &self.packages[id]
    }

    /// Every package in the graph (members and dependencies).
    pub fn all_packages(&self) -> impl Iterator<Item = &Package> {
        self.packages.values()
    }

    fn follows(kinds: &[DependencyKind], edges: Edges) -> bool {
        edges == Edges::All || kinds.iter().any(|k| *k != DependencyKind::Development)
    }

    /// Direct dependency names of a package (used by the app repo's layer rules).
    #[allow(dead_code)]
    pub fn direct(&self, id: &PackageId, edges: Edges) -> BTreeSet<String> {
        self.edges
            .get(id)
            .into_iter()
            .flatten()
            .filter(|(_, k)| Self::follows(k, edges))
            .map(|(d, _)| self.packages[d].name.to_string())
            .collect()
    }

    /// Transitive dependencies of a package: name -> the path that reaches it (for messages).
    /// With `Edges::All` only the first hop may be a dev-dependency (dev-deps of deps never build).
    pub fn transitive(&self, id: &PackageId, edges: Edges) -> BTreeMap<String, String> {
        let mut seen: BTreeMap<String, String> = BTreeMap::new();
        let mut queue = VecDeque::new();
        let start = self.packages[id].name.to_string();
        for (d, k) in self.edges.get(id).into_iter().flatten() {
            if Self::follows(k, edges) {
                queue.push_back((d.clone(), format!("{start} -> {}", self.packages[d].name)));
            }
        }
        while let Some((pid, via)) = queue.pop_front() {
            let name = self.packages[&pid].name.to_string();
            if seen.contains_key(&name) {
                continue;
            }
            for (d, k) in self.edges.get(&pid).into_iter().flatten() {
                if Self::follows(k, Edges::Shipping) {
                    queue.push_back((d.clone(), format!("{via} -> {}", self.packages[d].name)));
                }
            }
            seen.insert(name, via);
        }
        seen
    }

    /// True when a member has at least one `#[test]` (or `#[cfg(test)]`) under `src/` or `tests/`.
    pub fn has_tests(&self, id: &PackageId) -> bool {
        let dir = self.packages[id]
            .manifest_path
            .parent()
            .map_or_default(|p| p.as_std_path().to_path_buf());
        ["src", "tests"]
            .iter()
            .any(|sub| tree_contains(&dir.join(sub), "#[test]"))
    }
}

fn tree_contains(dir: &Path, needle: &str) -> bool {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return false;
    };
    rd.flatten().any(|e| {
        let p = e.path();
        if p.is_dir() {
            tree_contains(&p, needle)
        } else {
            p.extension().is_some_and(|x| x == "rs")
                && std::fs::read_to_string(&p).is_ok_and(|s| s.contains(needle))
        }
    })
}

/// Crate-name families the rules ban.
pub mod families {
    /// GDAL / PROJ bindings (C libraries).
    pub fn gdal_proj(name: &str) -> bool {
        matches!(name, "gdal" | "gdal-sys" | "gdal-src" | "proj" | "proj-sys")
    }
    /// Network stacks and async runtimes.
    pub fn network(name: &str) -> bool {
        matches!(
            name,
            "reqwest"
                | "tokio"
                | "hyper"
                | "ureq"
                | "curl"
                | "curl-sys"
                | "isahc"
                | "surf"
                | "async-std"
                | "attohttpc"
        )
    }
    /// GUI toolkits (the toolkit itself, not the renderer).
    pub fn toolkit(name: &str) -> bool {
        name == "egui"
            || name == "eframe"
            || name.starts_with("egui_")
            || name.starts_with("egui-")
            || name == "slint"
            || name.starts_with("slint-")
            || name.starts_with("i-slint-")
            || name == "iced"
            || name.starts_with("iced_")
            || name == "tauri"
            || name.starts_with("gpui")
            || name == "gtk"
            || name == "gtk4"
    }
    /// GPU, windowing and any toolkit: what an engine crate never pulls in.
    pub fn gui(name: &str) -> bool {
        toolkit(name) || matches!(name, "wgpu" | "winit" | "glow" | "glutin" | "softbuffer")
    }
}

/// Collected rule violations.
#[derive(Default)]
pub struct Report {
    lines: Vec<String>,
    checked: Vec<&'static str>,
}

impl Report {
    /// Records that a rule was checked.
    pub fn rule(&mut self, name: &'static str) {
        self.checked.push(name);
    }
    /// Records a violation.
    pub fn fail(&mut self, rule: &str, msg: impl AsRef<str>) {
        self.lines.push(format!("{rule}: {}", msg.as_ref()));
    }
    /// Prints the result; `Err` when anything failed.
    pub fn finish(self) -> Result<(), String> {
        if self.lines.is_empty() {
            println!("check-deps: OK ({})", self.checked.join(", "));
            Ok(())
        } else {
            for l in &self.lines {
                eprintln!("  {l}");
            }
            Err(format!("check-deps: {} violation(s)", self.lines.len()))
        }
    }
}
