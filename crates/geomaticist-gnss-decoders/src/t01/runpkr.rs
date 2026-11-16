//! Runs Trimble's runpkr00 to unpack a packed raw log (.T01, .T00, .T02, .R00) into its .DAT and .TGD files,
//! which [`super::records`] reads.
//!
//! It uses the user's own installed copy, located on the machine; runpkr00 is never bundled.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crate::DecodeError;

/// The packed formats runpkr00 unpacks (lower case, without the dot).
pub const RAW_EXTS: &[&str] = &["t01", "t00", "t02", "r00"];

fn exe_name() -> &'static str {
    if cfg!(windows) {
        "runpkr00.exe"
    } else {
        "runpkr00"
    }
}

fn on_path() -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|d| d.join(exe_name()))
        .find(|p| p.is_file())
}

/// runpkr00: a configured path, `$RUNPKR00`, a copy beside this program (or in its `bin/`), then `PATH`.
#[must_use]
pub fn find_runpkr00(configured: Option<&Path>) -> Option<PathBuf> {
    let mut cands: Vec<PathBuf> = Vec::new();
    if let Some(c) = configured {
        cands.push(c.to_path_buf());
    }
    if let Some(v) = std::env::var_os("RUNPKR00").filter(|v| !v.is_empty()) {
        cands.push(PathBuf::from(v));
    }
    if let Some(dir) = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(Path::to_path_buf))
    {
        cands.push(dir.join(exe_name()));
        cands.push(dir.join("bin").join(exe_name()));
    }
    if let Some(p) = on_path() {
        cands.push(p);
    }
    cands.into_iter().find(|c| c.is_file())
}

struct TempDir(PathBuf);

impl TempDir {
    fn new() -> Result<Self, DecodeError> {
        let base = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.subsec_nanos());
        for k in 0..100u32 {
            let p = base.join(format!("t01dec_{}_{nanos}_{k}", std::process::id()));
            if std::fs::create_dir(&p).is_ok() {
                return Ok(Self(p));
            }
        }
        Err(DecodeError::Io(format!(
            "couldn't make a temporary folder in {}",
            base.display()
        )))
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn run(
    exe: &Path,
    args: &[&std::ffi::OsStr],
    cwd: &Path,
    timeout: Duration,
) -> Result<(Option<i32>, String), DecodeError> {
    let mut child = Command::new(exe)
        .args(args)
        .current_dir(cwd)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| {
            DecodeError::Runpkr(format!(
                "runpkr00 didn't start ({e}). On Windows it needs Microsoft's Visual C++ 2008 SP1 (x86) runtime."
            ))
        })?;
    // drain both pipes so a chatty child never blocks
    let drain = |p: Option<Box<dyn std::io::Read + Send>>| {
        std::thread::spawn(move || {
            let mut buf = Vec::new();
            if let Some(mut p) = p {
                let _ = p.read_to_end(&mut buf);
            }
            buf
        })
    };
    let out_t = drain(
        child
            .stdout
            .take()
            .map(|p| Box::new(p) as Box<dyn std::io::Read + Send>),
    );
    let err_t = drain(
        child
            .stderr
            .take()
            .map(|p| Box::new(p) as Box<dyn std::io::Read + Send>),
    );
    let start = Instant::now();
    let status = loop {
        match child.try_wait() {
            Ok(Some(s)) => break s,
            Ok(None) if start.elapsed() > timeout => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(DecodeError::Runpkr(format!(
                    "runpkr00 timed out after {} s",
                    timeout.as_secs()
                )));
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(20)),
            Err(e) => return Err(DecodeError::Runpkr(format!("runpkr00 failed ({e})"))),
        }
    };
    let so = out_t.join().unwrap_or_default();
    let se = err_t.join().unwrap_or_default();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&so),
        String::from_utf8_lossy(&se)
    );
    Ok((status.code(), text))
}

/// `(.DAT bytes, .TGD bytes)` of one raw log. runpkr00 runs on a copy in a temporary folder (it writes its
/// outputs beside the name it is given), once with `-d` and once with `-g -d`.
pub fn unpack(raw: &Path, exe: Option<&Path>) -> Result<(Vec<u8>, Vec<u8>), DecodeError> {
    let exe = match exe {
        Some(e) => e.to_path_buf(),
        None => find_runpkr00(None)
            .ok_or_else(|| DecodeError::Runpkr("Trimble's runpkr00 isn't available".into()))?,
    };
    let name = raw
        .file_name()
        .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
    let td = TempDir::new()?;
    let ext = raw
        .extension()
        .map_or_else(String::new, |e| e.to_string_lossy().to_uppercase());
    let src = td.0.join(if ext.is_empty() {
        "in".to_owned()
    } else {
        format!("in.{ext}")
    });
    std::fs::copy(raw, &src).map_err(|e| DecodeError::io(raw, &e))?;
    let mut outs = Vec::new();
    for (tag, flags) in [("d", &["-d"][..]), ("g", &["-g", "-d"][..])] {
        let base = td.0.join(tag);
        let mut args: Vec<&std::ffi::OsStr> = flags.iter().map(std::ffi::OsStr::new).collect();
        args.push(src.as_os_str());
        args.push(base.as_os_str());
        let (code, text) = run(&exe, &args, &td.0, Duration::from_secs(600))?;
        let found: Vec<PathBuf> = ["dat", "tgd"]
            .iter()
            .map(|x| base.with_extension(x))
            .filter(|p| p.exists())
            .collect();
        let Some(last) = found.last() else {
            let t = text.trim();
            let tail: String = {
                let chars: Vec<char> = t.chars().collect();
                chars[chars.len().saturating_sub(300)..].iter().collect()
            };
            let msg = if cfg!(windows) && matches!(code, Some(-1_073_741_515 | 14001)) {
                "it couldn't start: install Microsoft's Visual C++ 2008 SP1 (x86) runtime"
                    .to_owned()
            } else if tail.is_empty() {
                format!("exit {}", code.unwrap_or(-1))
            } else {
                tail
            };
            return Err(DecodeError::Runpkr(format!(
                "runpkr00 wrote nothing for {name} ({msg})"
            )));
        };
        outs.push(std::fs::read(last).map_err(|e| DecodeError::io(last, &e))?);
    }
    let tgd = outs.pop().unwrap_or_default();
    let dat = outs.pop().unwrap_or_default();
    Ok((dat, tgd))
}
