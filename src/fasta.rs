//! Small hand-written FASTA reader.
//!
//! Genomes here are at most ~10 Mb, so the whole file is read into memory and
//! every sequence is uppercased. When a `.fai` index is missing it is created
//! with `samtools faidx`; the index is otherwise unused, but the plan requires
//! the assemblies to be indexed.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};

/// FASTA sequences plus their original file order.
pub struct Fasta {
    order: Vec<String>,
    seqs: HashMap<String, Vec<u8>>,
    /// Full header text after `>` per contig (name + description), for tags
    /// like `circular=true`.
    headers: HashMap<String, String>,
}

impl Fasta {
    /// Load a FASTA file into memory, uppercasing all bases.
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading FASTA {}", path.display()))?;
        let mut order: Vec<String> = Vec::new();
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
        let mut headers: HashMap<String, String> = HashMap::new();
        let mut cur: Option<String> = None;
        let mut buf: Vec<u8> = Vec::new();
        for line in text.lines() {
            let line = line.trim_end();
            if let Some(rest) = line.strip_prefix('>') {
                if let Some(name) = cur.take() {
                    seqs.insert(name, std::mem::take(&mut buf));
                }
                let name = rest.split_whitespace().next().unwrap_or("").to_string();
                if name.is_empty() {
                    bail!("FASTA {}: empty sequence name", path.display());
                }
                if seqs.contains_key(&name) {
                    bail!("FASTA {}: duplicate sequence name {name:?}", path.display());
                }
                order.push(name.clone());
                headers.insert(name.clone(), rest.trim().to_string());
                cur = Some(name);
            } else {
                for b in line.bytes() {
                    if !b.is_ascii_whitespace() {
                        buf.push(b.to_ascii_uppercase());
                    }
                }
            }
        }
        if let Some(name) = cur.take() {
            seqs.insert(name, buf);
        }
        if order.is_empty() {
            bail!("FASTA {}: no sequences", path.display());
        }
        Ok(Self { order, seqs, headers })
    }

    /// Contig names in file order.
    pub fn order(&self) -> &[String] {
        &self.order
    }

    /// Full header text (after `>`) for a contig, or `None` if absent.
    pub fn header_of(&self, name: &str) -> Option<&str> {
        self.headers.get(name).map(|s| s.as_str())
    }

    /// Sequence bases for a contig, or `None` if absent.
    pub fn seq(&self, name: &str) -> Option<&[u8]> {
        self.seqs.get(name).map(|v| v.as_slice())
    }

    /// Length of a contig, or `None` if absent.
    pub fn len_of(&self, name: &str) -> Option<u64> {
        self.seqs.get(name).map(|v| v.len() as u64)
    }
}

/// Ensure `<path>.fai` exists, creating it with `samtools faidx` if needed.
pub fn ensure_fai(path: &Path) -> Result<()> {
    let mut fai = path.as_os_str().to_owned();
    fai.push(".fai");
    let fai = PathBuf::from(fai);
    if fai.exists() {
        return Ok(());
    }
    let out = Command::new("samtools")
        .arg("faidx")
        .arg(path)
        .output()
        .with_context(|| format!("running samtools faidx on {}", path.display()))?;
    if !out.status.success() {
        bail!(
            "samtools faidx {} failed: {}",
            path.display(),
            String::from_utf8_lossy(&out.stderr).trim()
        );
    }
    if !fai.exists() {
        bail!("samtools faidx {} did not create {}", path.display(), fai.display());
    }
    Ok(())
}
