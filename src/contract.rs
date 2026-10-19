//! Input-contract check for `cramhook translate` (Plan section 2).
//!
//! The input CRAM must have been produced by the documented minimap2 command
//! shape. This module parses the SAM header and verifies:
//!
//! * `@HD` says `SO:coordinate`;
//! * exactly one minimap2 `@PG` line, whose `CL` uses only the allowed options
//!   and one of the two supported presets, and includes both `--secondary=no`
//!   and `-Y`;
//! * any other `@PG` line is a `samtools sort` invocation, or a `cramhook lift`
//!   `@PG` (a lifted CRAM used as input; its `DS` `remap=` command must use the
//!   same preset) followed by its own `samtools sort`;
//! * `@SQ` names and lengths equal the old FASTA's `.fai` (a wrong reference is
//!   a hard error, not a contract fallback).
//!
//! Everything is returned as a reason string on failure so the caller can log
//! it, defer every read and still write outputs (unless `--strict`).

use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::diff::Preset;

/// The parsed contract when valid.
#[derive(Debug, Clone, PartialEq)]
pub struct Contract {
    pub preset: Preset,
    /// The input is a `cramhook lift` chain (it carries at least one cramhook
    /// `lift` `@PG`). The caller must additionally require its own
    /// eligibility/safety options to be at their defaults before trusting an
    /// inherited `XL:A:T` record (see `TranslateConfig::eligibility_is_default`).
    pub chained: bool,
}

/// `lift` CLI options that do **not** affect which records are translated or
/// how: they take one value.
const LIFT_NEUTRAL_VALUE_OPTS: &[&str] = &[
    "--in", "--old-ref", "--new-ref", "--out", "--threads", "--stats", "--min-input-mb", "--preset",
];

/// `lift` CLI flags that do **not** affect which records are translated or how.
const LIFT_NEUTRAL_FLAG_OPTS: &[&str] = &["--strict", "--keep-tmp", "--xl-tag", "--no-xl-tag"];

/// `lift` options that **do** affect eligibility or safety. Any of these in an
/// earlier lift's `CL` means an inherited `XL:A:T` record's stripped `s2` (or the
/// clips an earlier lift screened) cannot be trusted, so the contract fails and
/// the input is fully remapped.
const LIFT_GATED_OPTS: &[&str] = &[
    "--min-mapq",
    "--max-clip",
    "--margin",
    "--max-de",
    "--t2-win",
    "--t2-step",
    "--t2-flank",
    "--t2-junction-flank",
    "--clip-k",
    "--clip-min-hits",
    "--clip-min-len",
    "--clip-flank",
    "--unsafe-no-clip-screen",
    "--clip-subtract-old",
    "--max-s2-ratio",
    "--max-minus-indel",
    "--max-zdrop",
];

/// `lift` options every recorded `CL` must contain (the required paths).
const LIFT_REQUIRED_OPTS: &[&str] = &["--in", "--old-ref", "--new-ref", "--out"];

/// Split a shell word list as written by `translate::shell_quote` / `shell_join`:
/// plain words and single-quoted words, with an embedded `'` written as the
/// `'\''` idiom. Any double quote, any backslash outside that idiom, or an
/// unbalanced quote returns `None` (an unparseable `CL`, so the chain is
/// rejected and the input remapped).
pub fn shell_split(s: &str) -> Option<Vec<String>> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut started = false;
    let mut in_single = false;
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if in_single {
            if c == '\'' {
                in_single = false;
            } else {
                cur.push(c);
            }
            continue;
        }
        match c {
            '\'' => {
                in_single = true;
                started = true;
            }
            '\\' => {
                // Only the `'\''` idiom escapes a quote; any other backslash is
                // not a form our writer produces.
                if chars.peek() == Some(&'\'') {
                    chars.next();
                    cur.push('\'');
                    started = true;
                } else {
                    return None;
                }
            }
            '"' => return None,
            c if c.is_whitespace() => {
                if started {
                    out.push(std::mem::take(&mut cur));
                    started = false;
                }
            }
            c => {
                cur.push(c);
                started = true;
            }
        }
    }
    if in_single {
        return None;
    }
    if started {
        out.push(cur);
    }
    Some(out)
}

/// Validate a value of a neutral `lift` option, using the CLI's grammar.
fn validate_neutral_value(name: &str, v: &str) -> Result<(), String> {
    if v.is_empty() {
        return Err(format!("option {name} has an empty value"));
    }
    if v.starts_with("--") {
        return Err(format!("option {name} has no value (found {v:?})"));
    }
    match name {
        "--threads" => {
            v.parse::<usize>()
                .map_err(|_| format!("option --threads value {v:?} is not a non-negative integer"))?;
        }
        "--min-input-mb" => {
            let x: f64 = v.parse().map_err(|_| format!("option --min-input-mb value {v:?} is not a number"))?;
            if !x.is_finite() || x < 0.0 {
                return Err(format!("option --min-input-mb value {v:?} must be finite and >= 0"));
            }
        }
        "--preset" => {
            if v != "map-hifi" && v != "lr:hq" {
                return Err(format!("option --preset value {v:?} is not map-hifi or lr:hq"));
            }
        }
        // Paths: non-empty and not another option (a leading `-` would be a flag).
        "--in" | "--old-ref" | "--new-ref" | "--out" | "--stats" if v.starts_with('-') => {
            return Err(format!("option {name} value {v:?} must not start with '-'"));
        }
        _ => {}
    }
    Ok(())
}

/// Validate a previous `cramhook lift` `CL`: only neutral options, each with its
/// one value (flags take none; `--opt value` and `--opt=value` both accepted), no
/// option repeated, every value well-formed, and the required paths present. Any
/// unknown token, a missing/duplicated/invalid value, or an eligibility/safety
/// option is a precise error (the input is then remapped).
fn validate_previous_lift_cl(toks: &[String]) -> Result<(), String> {
    let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut i = 2usize; // toks[0] = binary, toks[1] = "lift"
    while i < toks.len() {
        let t = toks[i].as_str();
        // `--opt=value` is one token.
        let (name, inline) = match t.split_once('=') {
            Some((n, v)) if n.starts_with("--") => (n, Some(v)),
            _ => (t, None),
        };
        if LIFT_GATED_OPTS.contains(&name) {
            return Err(format!(
                "earlier lift used {name}, which affects eligibility/safety; chained input rejected"
            ));
        }
        if LIFT_NEUTRAL_VALUE_OPTS.contains(&name) {
            let value = match inline {
                Some(v) => v,
                None => match toks.get(i + 1) {
                    Some(v) => v.as_str(),
                    None => return Err(format!("option {name} has no value")),
                },
            };
            validate_neutral_value(name, value)?;
            if !seen.insert(name) {
                return Err(format!("option {name} appears more than once"));
            }
            i += if inline.is_some() { 1 } else { 2 };
        } else if inline.is_none() && LIFT_NEUTRAL_FLAG_OPTS.contains(&t) {
            if !seen.insert(t) {
                return Err(format!("option {t} appears more than once"));
            }
            i += 1;
        } else if t.starts_with('-') {
            return Err(format!("unknown option {t} in an earlier lift CL"));
        } else {
            return Err(format!("unexpected argument {t:?} in an earlier lift CL"));
        }
    }
    for req in LIFT_REQUIRED_OPTS {
        if !seen.contains(req) {
            return Err(format!("an earlier lift CL is missing required option {req}"));
        }
    }
    Ok(())
}

/// One parsed header line (`@XX` + `TAG:VALUE` fields).
#[derive(Debug, Clone)]
pub struct HeaderLine {
    pub rec_type: String,
    pub fields: Vec<(String, String)>,
}

impl HeaderLine {
    pub fn field(&self, tag: &str) -> Option<&str> {
        self.fields.iter().find(|(t, _)| t == tag).map(|(_, v)| v.as_str())
    }
}

/// Parse raw SAM header text into lines. Blank lines are skipped.
pub fn parse_header(text: &str) -> Vec<HeaderLine> {
    let mut out = Vec::new();
    for line in text.split('\n') {
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let Some(rec) = parts.next() else { continue };
        if !rec.starts_with('@') {
            continue;
        }
        let mut fields = Vec::new();
        for p in parts {
            if let Some((tag, val)) = p.split_once(':') {
                fields.push((tag.to_string(), val.to_string()));
            }
        }
        out.push(HeaderLine { rec_type: rec[1..].to_string(), fields });
    }
    out
}

/// Parse a minimap2 `@PG` `CL` field, returning the preset.
///
/// Grammar (whitespace-separated): `minimap2`, then the preset in any of the
/// supported spellings, `-t N`/`-tN`, `--secondary=no`, `-Y`, then the old
/// reference path and at least one read file. No other option is allowed.
pub fn parse_minimap2_cl(cl: &str) -> Result<Preset, String> {
    let toks: Vec<&str> = cl.split_ascii_whitespace().collect();
    if toks.is_empty() {
        return Err("empty CL".into());
    }
    let bin = basename(toks[0]);
    if bin != "minimap2" {
        return Err(format!("first CL token is not minimap2: {:?}", toks[0]));
    }

    let mut preset: Option<Preset> = None;
    let mut saw_secondary = false;
    let mut saw_y = false;
    let mut saw_positional = false;
    let mut reads_after_ref = 0usize;
    let mut i = 1usize;
    while i < toks.len() {
        let t = toks[i];
        if saw_positional {
            // After the reference path only read files are allowed; a lone `-`
            // (stdin) is fine, any other option token is a violation.
            if t != "-" && t.starts_with('-') {
                return Err(format!("unexpected option after the reference path: {t:?}"));
            }
            reads_after_ref += 1;
            i += 1;
            continue;
        }
        if t == "-" || !t.starts_with('-') {
            // The first positional is the reference path.
            saw_positional = true;
            reads_after_ref = 0;
            i += 1;
            continue;
        }

        match t {
            "-ax" | "-x" => {
                let Some(p) = toks.get(i + 1) else {
                    return Err(format!("{t} without a preset"));
                };
                if preset.is_some() {
                    return Err("repeated preset".into());
                }
                preset = Some(parse_preset(p)?);
                i += 2;
            }
            "-a" => {
                i += 1;
            }
            "-Y" => {
                saw_y = true;
                i += 1;
            }
            "--secondary=no" => {
                saw_secondary = true;
                i += 1;
            }
            "--secondary=yes" => {
                return Err("--secondary=yes is not allowed (requires --secondary=no)".into());
            }
            "-t" => {
                let Some(n) = toks.get(i + 1) else {
                    return Err("-t without a value".into());
                };
                if n.parse::<u32>().is_err() {
                    return Err(format!("-t has a non-numeric value {n:?}"));
                }
                i += 2;
            }
            _ if t.starts_with("-x") && t.len() > 2 => {
                if preset.is_some() {
                    return Err("repeated preset".into());
                }
                preset = Some(parse_preset(&t[2..])?);
                i += 1;
            }
            _ if t.starts_with("-t") && t.len() > 2 && t[2..].chars().all(|c| c.is_ascii_digit()) => {
                i += 1;
            }
            // Reject combined short-flag bundles (e.g. `-aYx`); only the exact
            // `-ax`, `-xP` and `-tN` spellings are parsed above.
            _ if t.starts_with("-a") && t.len() > 2 => {
                return Err(format!("combined short options are not allowed: {t:?}"));
            }
            other => {
                return Err(format!("unexpected minimap2 option {other:?}"));
            }
        }
    }

    let Some(preset) = preset else {
        return Err("missing preset (-ax/-x)".into());
    };
    if !saw_secondary {
        return Err("missing --secondary=no".into());
    }
    if !saw_y {
        return Err("missing -Y".into());
    }
    if !saw_positional {
        return Err("missing reference path".into());
    }
    if reads_after_ref < 1 {
        return Err("no read files after the reference path".into());
    }
    Ok(preset)
}

fn parse_preset(p: &str) -> Result<Preset, String> {
    match p {
        "map-hifi" => Ok(Preset::MapHifi),
        "lr:hq" => Ok(Preset::LrHq),
        other => Err(format!("unsupported preset {other:?}")),
    }
}

fn basename(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

/// True when `l` looks like a minimap2 `@PG` (via `PN` or an `ID` prefix).
fn is_minimap2(l: &HeaderLine) -> bool {
    l.field("PN") == Some("minimap2") || l.field("ID").is_some_and(|id| id.starts_with("minimap2"))
}

/// True when `l` looks like a cramhook `@PG` (via `PN` or an `ID` prefix).
fn is_cramhook(l: &HeaderLine) -> bool {
    l.field("PN") == Some("cramhook") || l.field("ID").is_some_and(|id| id.starts_with("cramhook"))
}

/// True when `l` looks like a lift `@PG` written before the rename (`cramlift`
/// was the tool's old name): `PN` or an `ID` prefix of `cramlift`.
fn is_old_name_lift(l: &HeaderLine) -> bool {
    l.field("PN") == Some("cramlift") || l.field("ID").is_some_and(|id| id.starts_with("cramlift"))
}

/// No compatibility: fail when the input carries a lift `@PG` written by the
/// tool before the rename.
pub fn reject_old_name_lift(header_text: &str) -> Result<(), String> {
    if parse_header(header_text).iter().filter(|l| l.rec_type == "PG").any(is_old_name_lift) {
        return Err(
            "input was lifted by cramlift (the tool's old name); remap the reads to old.fa and rerun"
                .into(),
        );
    }
    Ok(())
}

/// True when `l` is a `samtools sort` invocation (`samtools view`/`index`/... are
/// not).
fn is_samtools_sort(l: &HeaderLine) -> bool {
    let pn = l.field("PN").unwrap_or("");
    let cl = l.field("CL").unwrap_or("");
    let toks: Vec<&str> = cl.split_ascii_whitespace().collect();
    pn == "samtools" && toks.len() >= 2 && basename(toks[0]) == "samtools" && toks[1] == "sort"
}

/// Check one `cramhook` `@PG` produced by `cramhook lift`: the `CL` must be a
/// `lift` invocation using only neutral options (see the allowlist), and the
/// `DS` `remap=` minimap2 command must be valid and use the same preset as the
/// original mapping.
fn check_cramhook_lift(l: &HeaderLine, preset: Preset) -> Result<(), String> {
    let cl = l.field("CL").unwrap_or("");
    let toks = shell_split(cl).ok_or_else(|| format!("unbalanced quoting in cramhook lift CL {:.60?}", cl))?;
    if toks.len() < 2 || basename(&toks[0]) != "cramhook" || toks[1] != "lift" {
        return Err(format!(
            "cramhook @PG is not a `lift` invocation (CL starts {:.40?})",
            cl
        ));
    }
    validate_previous_lift_cl(&toks)?;
    let ds = l.field("DS").unwrap_or("");
    // A chained `@PG` must carry the argv-encoding marker: an unmarked one is a
    // legacy `.join(" ")` CL whose apostrophes can hide a gated option.
    if !ds.split_ascii_whitespace().any(|t| t == crate::translate::ARGV_ENCODING_MARKER) {
        return Err(
            "lift @PG from an older cramhook (no argv encoding marker); not chaining".into(),
        );
    }
    let Some((_, remap)) = ds.split_once("remap=") else {
        return Err("cramhook lift @PG has no DS remap= command".into());
    };
    let remap = remap.trim();
    if remap.is_empty() {
        return Err("cramhook lift @PG has an empty DS remap= command".into());
    }
    let p = parse_minimap2_cl(remap)
        .map_err(|e| format!("cramhook lift DS remap command is invalid: {e}"))?;
    if p != preset {
        return Err(format!(
            "cramhook lift DS remap preset {} != input preset {}",
            p.minimap2_name(),
            preset.minimap2_name()
        ));
    }
    Ok(())
}

/// Check the minimap2 `@PG` `CL` and every other `@PG` line, returning the
/// preset or a precise violation reason.
pub fn check_contract(header_text: &str) -> Result<Contract, String> {
    let lines = parse_header(header_text);

    let mut hd_ok = false;
    for l in &lines {
        if l.rec_type == "HD" {
            if !l.fields.iter().any(|(t, _)| t == "SO") {
                return Err("@HD is missing SO".into());
            }
            match l.field("SO") {
                Some("coordinate") => hd_ok = true,
                Some(other) => return Err(format!("@HD SO:{other} is not coordinate")),
                None => return Err("@HD is missing SO".into()),
            }
        }
    }
    if !hd_ok {
        return Err("@HD SO:coordinate is required".into());
    }

    // The single minimap2 `@PG` must be the first `@PG`; the produced chain
    // always starts at the original mapping.
    let pgs: Vec<&HeaderLine> = lines.iter().filter(|l| l.rec_type == "PG").collect();
    reject_old_name_lift(header_text)?;

    let mm_positions: Vec<usize> =
        pgs.iter().enumerate().filter(|(_, l)| is_minimap2(l)).map(|(i, _)| i).collect();
    if mm_positions.is_empty() {
        return Err("no minimap2 @PG line".into());
    }
    if mm_positions.len() > 1 {
        return Err(format!("{} minimap2 @PG lines (expected exactly one)", mm_positions.len()));
    }
    if mm_positions[0] != 0 {
        return Err("the minimap2 @PG must be the first @PG line".into());
    }

    let mm = pgs[0];
    let cl = mm.field("CL").ok_or_else(|| "minimap2 @PG has no CL".to_string())?;
    let preset = parse_minimap2_cl(cl)?;

    // After minimap2: `S* (L S+)*`. Every `@PG` is either a `samtools sort` or a
    // `cramhook lift` `@PG` (a lifted input), and each lift must be followed by
    // its own `samtools sort`. Anything else is rejected with a precise reason.
    let mut chained = false;
    let mut i = 1usize;
    while i < pgs.len() && is_samtools_sort(pgs[i]) {
        i += 1;
    }
    while i < pgs.len() {
        let l = pgs[i];
        if !is_cramhook(l) {
            let pn = l.field("PN").unwrap_or("");
            let cl = l.field("CL").unwrap_or("");
            return Err(format!(
                "unexpected @PG line (PN:{pn:?}, CL starts with {:.40?})",
                cl
            ));
        }
        check_cramhook_lift(l, preset)?;
        chained = true;
        i += 1;
        let mut sorts = 0usize;
        while i < pgs.len() && is_samtools_sort(pgs[i]) {
            sorts += 1;
            i += 1;
        }
        if sorts == 0 {
            return Err("cramhook lift @PG is not followed by a samtools sort @PG".into());
        }
    }

    Ok(Contract { preset, chained })
}

/// Read a `.fai` index into `(name, length)` in file order.
pub fn read_fai(path: &Path) -> Result<Vec<(String, u64)>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading fai {}", path.display()))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 {
            bail!("{}:{}: expected at least 2 columns", path.display(), i + 1);
        }
        let len: u64 = f[1]
            .parse()
            .with_context(|| format!("{}:{}: bad length {:?}", path.display(), i + 1, f[1]))?;
        out.push((f[0].to_string(), len));
    }
    Ok(out)
}

/// Compare header `@SQ` names/lengths with the old FASTA `.fai`. A mismatch is
/// a hard error (wrong reference).
pub fn check_sq_against_fai(header_text: &str, fai: &[(String, u64)]) -> Result<()> {
    let lines = parse_header(header_text);
    let mut sq: Vec<(String, u64)> = Vec::new();
    for l in &lines {
        if l.rec_type != "SQ" {
            continue;
        }
        let name = l.field("SN").ok_or_else(|| anyhow::anyhow!("@SQ without SN"))?;
        let len: u64 = l
            .field("LN")
            .ok_or_else(|| anyhow::anyhow!("@SQ {name} without LN"))?
            .parse()
            .with_context(|| format!("@SQ {name} has a non-numeric LN"))?;
        sq.push((name.to_string(), len));
    }
    if sq != fai {
        bail!(
            "@SQ does not match the old reference .fai (wrong reference)\n  header: {sq:?}\n  fai:    {fai:?}"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cl_valid() -> String {
        "minimap2 -ax map-hifi -t 16 --secondary=no -Y old.fa reads.fastq.gz".to_string()
    }

    #[test]
    fn cl_accepts_spellings() {
        assert_eq!(parse_minimap2_cl(&cl_valid()), Ok(Preset::MapHifi));
        assert_eq!(
            parse_minimap2_cl("minimap2 -a -x lr:hq -t8 --secondary=no -Y o.fa r.fq"),
            Ok(Preset::LrHq)
        );
        assert_eq!(
            parse_minimap2_cl("minimap2 -x map-hifi -a --secondary=no -Y o.fa r.fq"),
            Ok(Preset::MapHifi)
        );
        assert_eq!(
            parse_minimap2_cl("minimap2 -xmap-hifi --secondary=no -Y o.fa r.fq"),
            Ok(Preset::MapHifi)
        );
        assert_eq!(
            parse_minimap2_cl("/usr/bin/minimap2 -ax lr:hq --secondary=no -Y o.fa r.fq"),
            Ok(Preset::LrHq)
        );
        // repeated -t is allowed
        assert_eq!(
            parse_minimap2_cl("minimap2 -ax map-hifi -t 4 -t8 --secondary=no -Y o.fa r.fq"),
            Ok(Preset::MapHifi)
        );
    }

    #[test]
    fn cl_rejects_options_after_positionals() {
        // any option token after the reference path is a violation
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no -Y o.fa r.fq -k 11").is_err());
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no -Y o.fa r.fq --MD").is_err());
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no -Y o.fa r.fq -y").is_err());
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no -Y o.fa r.fq r2.fq -Y").is_err());
    }

    #[test]
    fn cl_rejects_repeated_preset() {
        // preset in two spellings
        assert!(
            parse_minimap2_cl("minimap2 -ax map-hifi -x lr:hq --secondary=no -Y o.fa r.fq").is_err()
        );
        assert!(
            parse_minimap2_cl("minimap2 -xmap-hifi -ax lr:hq --secondary=no -Y o.fa r.fq").is_err()
        );
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi -ax map-hifi --secondary=no -Y o.fa r.fq").is_err());
    }

    #[test]
    fn cl_rejects_combined_short_bundles() {
        assert!(parse_minimap2_cl("minimap2 -aYx map-hifi --secondary=no -Y o.fa r.fq").is_err());
        assert!(parse_minimap2_cl("minimap2 -at8 map-hifi --secondary=no -Y o.fa r.fq").is_err());
    }

    #[test]
    fn cl_rejects_violations() {
        // --secondary=yes
        assert!(
            parse_minimap2_cl("minimap2 -ax map-hifi --secondary=yes -Y o.fa r.fq").is_err()
        );
        // missing -Y
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no o.fa r.fq").is_err());
        // missing --secondary=no
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi -Y o.fa r.fq").is_err());
        // -k 15
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi -k 15 --secondary=no -Y o.fa r.fq").is_err());
        // --MD / --cs
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --MD --secondary=no -Y o.fa r.fq").is_err());
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --cs --secondary=no -Y o.fa r.fq").is_err());
        // map-ont preset
        assert!(parse_minimap2_cl("minimap2 -ax map-ont --secondary=no -Y o.fa r.fq").is_err());
        // no read files
        assert!(parse_minimap2_cl("minimap2 -ax map-hifi --secondary=no -Y o.fa").is_err());
        // wrong binary
        assert!(parse_minimap2_cl("bwa -ax map-hifi --secondary=no -Y o.fa r.fq").is_err());
    }

    fn header(mm_cl: &str, extra: &[&str]) -> String {
        let mut s = String::from("@HD\tVN:1.6\tSO:coordinate\n");
        s.push_str(&format!(
            "@PG\tID:minimap2\tPN:minimap2\tVN:2.31\tCL:{mm_cl}\n"
        ));
        for e in extra {
            s.push_str(e);
            s.push('\n');
        }
        s
    }

    #[test]
    fn contract_accepts_sort_only() {
        let h = header(
            &cl_valid(),
            &["@PG\tID:samtools\tPN:samtools\tPP:minimap2\tCL:samtools sort -O CRAM -o old.cram"],
        );
        assert_eq!(
            check_contract(&h),
            Ok(Contract { preset: Preset::MapHifi, chained: false })
        );
    }

    #[test]
    fn contract_rejects_bad_extra_pg() {
        let h = header(
            &cl_valid(),
            &["@PG\tID:samtools.1\tPN:samtools\tCL:/x/samtools view -H old.cram"],
        );
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_two_minimap2() {
        let h = header(
            &cl_valid(),
            &["@PG\tID:minimap2.1\tPN:minimap2\tCL:minimap2 -ax map-hifi --secondary=no -Y o.fa r.fq"],
        );
        assert!(check_contract(&h).is_err());
    }

    /// A `cramhook lift` `@PG` as `cramhook` writes it (M4): a neutral,
    /// all-default CL plus `cl_extra` (which the caller may set to a gated or
    /// malformed option).
    fn cramhook_lift_pg(remap: &str, cl_extra: &str) -> String {
        format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:/x/cramhook lift --in old.cram --old-ref old.fa --new-ref new.fa --out lift.cram --threads 16 --stats s.json{cl_extra}\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap} {m}",
            m = crate::translate::ARGV_ENCODING_MARKER
        )
    }

    fn sort_pg(id: &str) -> String {
        format!("@PG\tID:{id}\tPN:samtools\tVN:1.24\tCL:samtools sort -O CRAM -o out.cram -")
    }

    /// A `cramhook lift` `@PG` with a raw `CL` (for malformed-CL tests).
    fn cramhook_lift_pg_raw(cl: &str, remap: &str) -> String {
        format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:{cl}\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap} {m}",
            m = crate::translate::ARGV_ENCODING_MARKER
        )
    }

    /// A legacy `cramhook lift` `@PG` (older build): `.join(" ")` CL, no argv marker.
    fn legacy_lift_pg_raw(cl: &str, remap: &str) -> String {
        format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:{cl}\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap}"
        )
    }

    fn check_raw_cl(cl: &str) -> Result<Contract, String> {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let pg = cramhook_lift_pg_raw(cl, remap);
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        check_contract(&h)
    }

    fn chain(remap: &str, extras: &[&str]) -> String {
        let mut pgs: Vec<String> = vec![sort_pg("samtools")];
        for (i, e) in extras.iter().enumerate() {
            pgs.push(cramhook_lift_pg(remap, e));
            pgs.push(sort_pg(&format!("samtools.{}", i + 1)));
        }
        let refs: Vec<&str> = pgs.iter().map(String::as_str).collect();
        header(&cl_valid(), &refs)
    }

    #[test]
    fn contract_accepts_chained_lift() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        assert_eq!(
            check_contract(&chain(remap, &[""])),
            Ok(Contract { preset: Preset::MapHifi, chained: true })
        );
    }

    #[test]
    fn contract_accepts_chained_lift_twice_with_neutral_options() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let h = chain(
            remap,
            &[
                "",
                " --min-input-mb=0 --keep-tmp --no-xl-tag --preset map-hifi --strict",
            ],
        );
        assert_eq!(check_contract(&h), Ok(Contract { preset: Preset::MapHifi, chained: true }));
    }

    #[test]
    fn contract_rejects_non_default_previous_lift_options() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        for opt in [
            "--max-s2-ratio 1",
            "--max-s2-ratio=0.25",
            "--clip-min-hits 100000",
            "--unsafe-no-clip-screen",
            "--clip-subtract-old",
            "--min-mapq 0",
            "--max-de 1",
            "--max-clip 0",
            "--margin 1",
            "--t2-win 1",
            "--t2-step 1",
            "--t2-flank 1",
            "--clip-k 1",
            "--clip-min-len 1",
            "--clip-flank 1",
            "--max-minus-indel 0",
            "--max-zdrop 0",
        ] {
            let h = chain(remap, &[&format!(" {opt}")]);
            assert!(check_contract(&h).is_err(), "must reject earlier lift option {opt}");
        }
    }

    #[test]
    fn contract_rejects_unknown_or_malformed_previous_lift_cl() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        for extra in [
            " --frobnicate 1",          // unknown option
            " --clip-k",                // gated option with no value
            " --threads 4 --threads 8", // duplicated value
            " trailing",                // unexpected positional
            " -x lr:hq",                // short option
        ] {
            let h = chain(remap, &[extra]);
            assert!(check_contract(&h).is_err(), "must reject malformed CL extra {extra:?}");
        }
    }

    #[test]
    fn contract_rejects_legacy_lift_without_argv_marker() {
        // A legacy `.join(" ")` CL: apostrophes in the two paths read as quoting
        // and swallow `--max-s2-ratio 1 --stats `, so the gated option is
        // invisible to the option scan. The missing argv marker is the guard.
        let cl = "/x/cramhook lift --in in.cram --old-ref old.fa --new-ref new.fa --out legacy-out' --max-s2-ratio 1 --stats 'legacy-stats.json";
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        // Sanity: the splitter really does hide the gated option in that CL.
        let toks = shell_split(cl).expect("balanced");
        assert!(
            !toks.iter().any(|t| t == "--max-s2-ratio"),
            "the regression CL must hide --max-s2-ratio"
        );
        let pg = legacy_lift_pg_raw(cl, remap);
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        let e = check_contract(&h).unwrap_err();
        assert!(e.contains("older cramhook"), "{e}");
        assert!(e.contains("argv encoding marker"), "{e}");

        // A plain legacy CL (no hidden option) is rejected the same way.
        let pg2 = legacy_lift_pg_raw(
            "/x/cramhook lift --in in.cram --old-ref old.fa --new-ref new.fa --out out.cram",
            remap,
        );
        let h2 = header(&cl_valid(), &[&sort_pg("samtools"), &pg2, &sort_pg("samtools.1")]);
        assert!(check_contract(&h2).is_err(), "a legacy lift @PG must not chain");
    }

    #[test]
    fn contract_rejects_bad_neutral_option_values() {
        // `--threads` garbage / a missing value / a negative; `--min-input-mb`
        // NaN or negative; `--stats` swallowing the next flag; `--out=` (empty);
        // `--preset` not a known value; a path starting with `-`.
        let base = "/x/cramhook lift --in in.cram --old-ref old.fa --new-ref new.fa --out out.cram";
        for cl in [
            format!("{base} --threads garbage"),
            format!("{base} --threads --unsafe-no-clip-screen"),
            format!("{base} --threads -1"),
            format!("{base} --threads="),
            format!("{base} --stats --stats"),
            format!("{base} --min-input-mb NaN"),
            format!("{base} --min-input-mb -1"),
            format!("{base} --min-input-mb --stats x"),
            format!("{base} --out="),
            format!("{base} --in -x"),
            format!("{base} --preset map-ont"),
            format!("{base} --preset"),
        ] {
            assert!(check_raw_cl(&cl).is_err(), "must reject CL {cl:?}");
        }
        // The `--opt=value` neutral spelling is still accepted.
        let ok = format!("{base} --threads=8 --min-input-mb=0 --preset=lr:hq");
        assert_eq!(
            check_raw_cl(&ok),
            Ok(Contract { preset: Preset::MapHifi, chained: true })
        );
    }

    #[test]
    fn contract_rejects_unquoted_old_cl_with_spaced_path() {
        // An older build recorded an unquoted CL; a path containing
        // " --max-s2-ratio 0.25" then splits into tokens and is rejected.
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let pg = format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:/x/cramhook lift --in old.cram --old-ref old.fa --new-ref new.fa --out file --max-s2-ratio 0.25 --threads 16\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap}"
        );
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err(), "unquoted spaced path must be rejected");
    }

    #[test]
    fn contract_accepts_quoted_paths_with_spaces() {
        // `shell_join` quotes such a path; the contract then sees one value token.
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let cl = crate::translate::shell_join(&[
            "/x/cramhook".into(),
            "lift".into(),
            "--in".into(),
            "file --max-s2-ratio 0.25".into(),
            "--old-ref".into(),
            "old.fa".into(),
            "--new-ref".into(),
            "new.fa".into(),
            "--out".into(),
            "out file.cram".into(),
            "--threads".into(),
            "16".into(),
        ]);
        let pg = format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:{cl}\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap} {m}",
            m = crate::translate::ARGV_ENCODING_MARKER
        );
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        assert_eq!(check_contract(&h), Ok(Contract { preset: Preset::MapHifi, chained: true }));
    }

    #[test]
    fn contract_rejects_unbalanced_quoting() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let pg = format!(
            "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:/x/cramhook lift --in 'oops --old-ref o --new-ref n --out x\tPP:samtools\tDS:remap={remap}"
        );
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn shell_quote_and_split_round_trip() {
        let cases: Vec<Vec<String>> = vec![
            vec!["cramhook".into(), "lift".into(), "--in".into(), "old.cram".into()],
            vec!["/a b/cramhook".into(), "lift".into(), "--in".into(), "old file.cram".into()],
            vec![
                "cramhook".into(),
                "lift".into(),
                "--out".into(),
                "file --max-s2-ratio 0.25".into(),
            ],
            // Embedded single quotes (the `'\''` idiom) and other specials.
            vec!["weird 'quote'".into(), "a'b'c".into(), "back\\slash".into(), "$HOME".into()],
            vec![String::new()],
        ];
        for argv in cases {
            let joined = crate::translate::shell_join(&argv);
            assert_eq!(shell_split(&joined).unwrap(), argv, "round trip {joined:?}");
        }
        // Ordinary command lines are unchanged (no added quoting).
        assert_eq!(
            crate::translate::shell_join(&[
                "cramhook".into(),
                "lift".into(),
                "--in".into(),
                "old.cram".into()
            ]),
            "cramhook lift --in old.cram"
        );
        // The `'\''` idiom round-trips a value that is only a quote.
        assert_eq!(shell_split(&crate::translate::shell_join(&["'".into()])).unwrap(), vec!["'"]);
        // Only plain and single-quoted words are accepted: a double quote, a
        // backslash outside the `'\''` idiom, and an unbalanced quote are not.
        for bad in [
            "\"--thre\\ads\"", // double quote + stray backslash
            "foo\\",           // trailing backslash
            "back\\slash",     // backslash outside the idiom
            "a\"b",            // double quote in a plain word
            "'unterminated",   // unbalanced single quote
            "x' y",            // quote then never closed
        ] {
            assert!(shell_split(bad).is_none(), "must reject {bad:?}");
        }
    }

    #[test]
    fn contract_rejects_old_name_lift_pg() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let pg = format!(
            "@PG\tID:cramlift\tPN:cramlift\tVN:0.1.0\tCL:/x/cramlift lift --in old.cram --old-ref old.fa --new-ref new.fa --out lift.cram --threads 16 --stats s.json\tPP:samtools\tDS:old-ref=old.fa new-ref=new.fa remap={remap} {m}",
            m = crate::translate::ARGV_ENCODING_MARKER
        );
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &pg, &sort_pg("samtools.1")]);
        let e = check_contract(&h).unwrap_err();
        assert!(e.contains("old name"), "{e}");
    }

    #[test]
    fn contract_rejects_remap_preset_mismatch() {
        // minimap2 map-hifi but the cramhook remap used lr:hq.
        let remap = "minimap2 -ax lr:hq -t 16 --secondary=no -Y new.fa -";
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &cramhook_lift_pg(remap, ""), &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_unknown_cramhook_subcommand() {
        let pg = "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:/x/cramhook qc --in old.cram\tDS:old-ref=old.fa new-ref=new.fa remap=minimap2 -ax map-hifi --secondary=no -Y new.fa -";
        let h = header(&cl_valid(), &[&sort_pg("samtools"), pg, &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_cramhook_without_remap() {
        let pg = "@PG\tID:cramhook\tPN:cramhook\tVN:0.1.0\tCL:/x/cramhook lift --in old.cram\tDS:old-ref=old.fa new-ref=new.fa";
        let h = header(&cl_valid(), &[&sort_pg("samtools"), pg, &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_cramhook_remap_bad_options() {
        // The remap command must pass the same option rules (needs --secondary=no/-Y).
        let remap = "minimap2 -ax map-hifi -t 16 -Y new.fa -";
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &cramhook_lift_pg(remap, ""), &sort_pg("samtools.1")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_lift_not_followed_by_sort() {
        let remap = "minimap2 -ax map-hifi -t 16 --secondary=no -Y new.fa -";
        let h = header(&cl_valid(), &[&sort_pg("samtools"), &cramhook_lift_pg(remap, "")]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_non_sort_non_lift_program() {
        // A `bwa` @PG is neither a sort nor a cramhook lift.
        let h = header(&cl_valid(), &[&sort_pg("samtools"), "@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref.fa r.fq"]);
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_minimap2_not_first() {
        let h = format!(
            "@HD\tVN:1.6\tSO:coordinate\n{}\n@PG\tID:minimap2\tPN:minimap2\tVN:2.31\tCL:{}\n",
            sort_pg("samtools"),
            cl_valid()
        );
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_missing_pg() {
        let h = "@HD\tVN:1.6\tSO:coordinate\n".to_string();
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn contract_rejects_bad_hd() {
        let h = format!(
            "@HD\tVN:1.6\tSO:unsorted\n@PG\tID:minimap2\tPN:minimap2\tCL:{}\n",
            cl_valid()
        );
        assert!(check_contract(&h).is_err());
    }

    #[test]
    fn sq_check() {
        let h = "@SQ\tSN:chr\tLN:300000\n@SQ\tSN:p1\tLN:45000\n";
        let fai = vec![("chr".to_string(), 300000u64), ("p1".to_string(), 45000u64)];
        assert!(check_sq_against_fai(h, &fai).is_ok());
        let bad = vec![("chr".to_string(), 300001u64)];
        assert!(check_sq_against_fai(h, &bad).is_err());
    }
}
