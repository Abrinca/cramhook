//! Minimal PAF parsing and the short-form `cs` walker.
//!
//! Only what the diff needs: query/target names and intervals, strand, MAPQ,
//! the `tp:A:` tag (primary/secondary) and the `cs:Z:` tag. Anything unexpected
//! returns an error with context instead of panicking. Coordinates, cs bases and
//! cs consumption are validated so a malformed PAF cannot silently corrupt the
//! diff.

use anyhow::{Context, Result, bail};

/// Relative orientation of the query to the target in a PAF record.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(c: char) -> Result<Self> {
        match c {
            '+' => Ok(Strand::Plus),
            '-' => Ok(Strand::Minus),
            other => bail!("invalid PAF strand {other:?}"),
        }
    }

    pub fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }
}

/// One PAF alignment line.
#[derive(Clone, Debug)]
pub struct PafRecord {
    pub qname: String,
    pub qlen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub strand: Strand,
    pub tname: String,
    pub tlen: u64,
    pub tstart: u64,
    pub tend: u64,
    pub mapq: u32,
    /// `tp:A:P` (primary) vs. anything else (secondary/supplementary/N).
    pub primary: bool,
    /// Short-form `cs` string (without the `cs:Z:` prefix).
    pub cs: String,
}

impl PafRecord {
    /// Length of the aligned query interval.
    pub fn q_span(&self) -> u64 {
        self.qend - self.qstart
    }

    /// Length of the aligned target interval.
    pub fn t_span(&self) -> u64 {
        self.tend - self.tstart
    }
}

/// Parse a whole PAF file body. Empty lines are skipped.
pub fn parse_paf(text: &str) -> Result<Vec<PafRecord>> {
    let mut out = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let ctx = || format!("PAF line {}", lineno + 1);
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            bail!("{}: expected at least 12 columns, got {}", ctx(), f.len());
        }
        let parse_u64 = |i: usize| -> Result<u64> {
            f[i].parse::<u64>().with_context(|| format!("{}: column {i}", ctx()))
        };
        let qname = f[0].to_string();
        let qlen = parse_u64(1)?;
        let qstart = parse_u64(2)?;
        let qend = parse_u64(3)?;
        let strand = Strand::from_char(
            f[4].chars().next().with_context(|| format!("{}: empty strand", ctx()))?,
        )?;
        let tname = f[5].to_string();
        let tlen = parse_u64(6)?;
        let tstart = parse_u64(7)?;
        let tend = parse_u64(8)?;
        let mapq = f[11].parse::<u32>().with_context(|| format!("{}: column 11 (mapq)", ctx()))?;

        // Coordinate sanity.
        if qstart > qend || qend > qlen {
            bail!("{}: query interval {qstart}..{qend} out of range 0..{qlen}", ctx());
        }
        if tstart > tend || tend > tlen {
            bail!("{}: target interval {tstart}..{tend} out of range 0..{tlen}", ctx());
        }

        let mut primary = false;
        let mut cs: Option<String> = None;
        for tag in &f[12..] {
            if let Some(v) = tag.strip_prefix("tp:A:") {
                primary = v == "P";
            } else if let Some(v) = tag.strip_prefix("cs:Z:") {
                cs = Some(v.to_string());
            }
        }
        let cs = cs.with_context(|| format!("{}: missing cs:Z: tag", ctx()))?;

        // The cs must account for exactly the PAF spans.
        let (cq, ct) = walk_cs(&cs, |_, _, _| {})
            .with_context(|| format!("{}: invalid cs", ctx()))?;
        if cq != qend - qstart {
            bail!(
                "{}: cs consumes {cq} query bases but the PAF span is {}",
                ctx(),
                qend - qstart
            );
        }
        if ct != tend - tstart {
            bail!(
                "{}: cs consumes {ct} target bases but the PAF span is {}",
                ctx(),
                tend - tstart
            );
        }

        out.push(PafRecord {
            qname,
            qlen,
            qstart,
            qend,
            strand,
            tname,
            tlen,
            tstart,
            tend,
            mapq,
            primary,
            cs,
        });
    }
    Ok(out)
}

/// One step of a short-form `cs` walk.
///
/// `Match`/`Sub`/`Ins`/`Del` all carry the number of consumed bases. See
/// [`walk_cs`] for the cursor conventions.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CsStep {
    /// `:N` — N target bases and N query bases match.
    Match(u64),
    /// `*xy` — one target base and one query base differ.
    Sub,
    /// `+seq` — query-only bases (insertion to the target).
    Ins(u64),
    /// `-seq` — target-only bases (deletion from the query).
    Del(u64),
}

/// Walk a short-form `cs` string in target-forward order.
///
/// The callback is invoked for every op *before* the cursors advance, and
/// receives the target cursor `t` (bases consumed on the target so far) and the
/// query offset `qoff` (bases consumed on the query so far). Both start at 0.
///
/// For both strands the walk is in the forward direction of the *target*. The
/// `cs` encodes the query reverse-complemented for `-` alignments, so the caller
/// converts `qoff` to an actual query coordinate depending on the strand:
///
/// * `+`: query coordinate = `qstart + qoff + i`
/// * `-`: query coordinate = `qend - 1 - (qoff + i)`
///
/// All `cs` operations are validated; malformed input yields an error. Returns
/// `(query_consumed, target_consumed)` so callers can check consistency with the
/// PAF spans.
pub fn walk_cs<F>(cs: &str, mut f: F) -> Result<(u64, u64)>
where
    F: FnMut(CsStep, u64, u64),
{
    if cs.is_empty() {
        bail!("cs: empty");
    }
    let b = cs.as_bytes();
    let mut i = 0usize;
    let mut t: u64 = 0;
    let mut q: u64 = 0;
    while i < b.len() {
        match b[i] {
            b':' => {
                i += 1;
                let s = i;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                if s == i {
                    bail!("cs: ':' without a length");
                }
                let n: u64 = std::str::from_utf8(&b[s..i])
                    .unwrap()
                    .parse()
                    .context("cs: invalid match length")?;
                f(CsStep::Match(n), t, q);
                t = t.checked_add(n).context("cs: target cursor overflow")?;
                q = q.checked_add(n).context("cs: query cursor overflow")?;
            }
            b'*' => {
                if i + 2 >= b.len() {
                    bail!("cs: truncated substitution");
                }
                if !is_base(b[i + 1]) || !is_base(b[i + 2]) {
                    bail!(
                        "cs: invalid substitution bases {:?}",
                        String::from_utf8_lossy(&b[i + 1..i + 3])
                    );
                }
                f(CsStep::Sub, t, q);
                i += 3;
                t = t.checked_add(1).context("cs: target cursor overflow")?;
                q = q.checked_add(1).context("cs: query cursor overflow")?;
            }
            b'+' => {
                i += 1;
                let s = i;
                while i < b.len() && !is_op(b[i]) {
                    if !is_base(b[i]) {
                        bail!("cs: invalid '+' base {:?}", b[i] as char);
                    }
                    i += 1;
                }
                let n = (i - s) as u64;
                if n == 0 {
                    bail!("cs: '+' without bases");
                }
                f(CsStep::Ins(n), t, q);
                q = q.checked_add(n).context("cs: query cursor overflow")?;
            }
            b'-' => {
                i += 1;
                let s = i;
                while i < b.len() && !is_op(b[i]) {
                    if !is_base(b[i]) {
                        bail!("cs: invalid '-' base {:?}", b[i] as char);
                    }
                    i += 1;
                }
                let n = (i - s) as u64;
                if n == 0 {
                    bail!("cs: '-' without bases");
                }
                f(CsStep::Del(n), t, q);
                t = t.checked_add(n).context("cs: target cursor overflow")?;
            }
            other => bail!("cs: unexpected operation byte {:?}", other as char),
        }
    }
    Ok((q, t))
}

#[inline]
fn is_op(b: u8) -> bool {
    matches!(b, b':' | b'*' | b'+' | b'-')
}

/// `acgtn`, case-insensitive (minimap2 may emit lowercase soft-masked bases).
#[inline]
fn is_base(b: u8) -> bool {
    matches!(b.to_ascii_lowercase(), b'a' | b'c' | b'g' | b't' | b'n')
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(cs: &str) -> Vec<(CsStep, u64, u64)> {
        let mut v = Vec::new();
        walk_cs(cs, |s, t, q| v.push((s, t, q))).unwrap();
        v
    }

    #[test]
    fn plus_example() {
        // 10 match, 1 sub, 3 ins, 5 match, 2 del, 100 match
        let v = collect(":10*ac+ggg:5-tt:100");
        assert_eq!(
            v,
            vec![
                (CsStep::Match(10), 0, 0),
                (CsStep::Sub, 10, 10),
                (CsStep::Ins(3), 11, 11),
                (CsStep::Match(5), 11, 14),
                (CsStep::Del(2), 16, 19),
                (CsStep::Match(100), 18, 19),
            ]
        );
    }

    #[test]
    fn minus_query_intervals() {
        // Query qstart=0, qend=20, strand '-'. First :8 covers old [12,20),
        // the sub sits at old 11, and the trailing :11 covers old [0,11).
        // (cs is written in target-forward order with the query revcomp'd.)
        let qstart = 0u64;
        let qend = 20u64;
        let mut ivs: Vec<(u64, u64)> = Vec::new();
        let mut sub_pos = None;
        walk_cs(":8*ag:11", |step, _t, qoff| match step {
            CsStep::Match(n) => {
                let lo = qend - qoff - n;
                let hi = qend - qoff;
                assert!(qstart <= lo && hi <= qend);
                ivs.push((lo, hi));
            }
            CsStep::Sub => {
                sub_pos = Some(qend - 1 - qoff);
            }
            _ => unreachable!(),
        })
        .unwrap();
        assert_eq!(ivs, vec![(12, 20), (0, 11)]);
        assert_eq!(sub_pos, Some(11));
    }

    #[test]
    fn plus_query_intervals() {
        let qstart = 5u64;
        let mut ivs = Vec::new();
        walk_cs(":8*ag:11", |step, _t, qoff| {
            if let CsStep::Match(n) = step {
                ivs.push((qstart + qoff, qstart + qoff + n));
            }
        })
        .unwrap();
        assert_eq!(ivs, vec![(5, 13), (14, 25)]);
    }

    #[test]
    fn walk_returns_consumed_totals() {
        assert_eq!(walk_cs(":10*ac+ggg:5-tt:100", |_, _, _| {}).unwrap(), (119, 118));
    }

    #[test]
    fn accepts_lowercase_and_n() {
        assert!(walk_cs(":10*ac+acgtn:5-ttNN:100", |_, _, _| {}).is_ok());
    }

    #[test]
    fn rejects_empty_cs() {
        assert!(walk_cs("", |_, _, _| {}).is_err());
    }

    #[test]
    fn rejects_unknown_op() {
        for cs in [":10X", "=10", ":10~ag", ":10 5", "10", ":10:" ] {
            assert!(walk_cs(cs, |_, _, _| {}).is_err(), "should reject {cs:?}");
        }
    }

    #[test]
    fn rejects_invalid_bases() {
        for cs in [":10*xy", ":10*1a", "+aZ", "-x", ":5+z1:5", ":5-1:5", ":5+ :5"] {
            assert!(walk_cs(cs, |_, _, _| {}).is_err(), "should reject {cs:?}");
        }
    }

    #[test]
    fn rejects_truncated_and_empty_ops() {
        assert!(walk_cs(":10*", |_, _, _| {}).is_err());
        assert!(walk_cs(":10*a", |_, _, _| {}).is_err());
        assert!(walk_cs(":10+", |_, _, _| {}).is_err());
        assert!(walk_cs("+", |_, _, _| {}).is_err());
        assert!(walk_cs("-", |_, _, _| {}).is_err());
    }

    // --- parse_paf validation -------------------------------------------------

    /// Build one PAF line. Columns 9/10 (residue matches / alignment block
    /// length) are alignment properties, not derivable from the query interval
    /// alone, so the caller passes `block_len` explicitly. This helper is also
    /// used to build deliberately malformed lines (e.g. `qstart > qend`);
    /// deriving those columns from such an interval previously underflowed.
    fn paf_line(
        qlen: u64,
        qstart: u64,
        qend: u64,
        strand: char,
        tlen: u64,
        tstart: u64,
        tend: u64,
        block_len: u64,
        cs: &str,
    ) -> String {
        format!(
            "q\t{qlen}\t{qstart}\t{qend}\t{strand}\tt\t{tlen}\t{tstart}\t{tend}\t{block_len}\t{block_len}\t60\ttp:A:P\tcs:Z:{cs}"
        )
    }

    #[test]
    fn parse_paf_accepts_valid() {
        let line = paf_line(100, 10, 30, '+', 200, 5, 25, 20, ":20");
        let recs = parse_paf(&line).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].q_span(), 20);
        assert_eq!(recs[0].t_span(), 20);
    }

    #[test]
    fn parse_paf_rejects_query_range() {
        // qstart > qend (block length is irrelevant to this validation)
        assert!(parse_paf(&paf_line(100, 30, 10, '+', 200, 0, 10, 0, ":10")).is_err());
        // qend > qlen
        assert!(parse_paf(&paf_line(100, 0, 120, '+', 200, 0, 120, 120, ":120")).is_err());
    }

    #[test]
    fn parse_paf_rejects_target_range() {
        assert!(parse_paf(&paf_line(100, 0, 10, '+', 200, 30, 10, 10, ":10")).is_err());
        assert!(parse_paf(&paf_line(100, 0, 120, '+', 200, 0, 130, 120, ":120")).is_err());
    }

    #[test]
    fn parse_paf_rejects_consumption_mismatch() {
        // span 20 but cs consumes 19 query/target
        assert!(parse_paf(&paf_line(100, 10, 30, '+', 200, 5, 25, 20, ":19")).is_err());
        // span 20, cs consumes 20 target but only 15 query
        assert!(parse_paf(&paf_line(100, 10, 30, '+', 200, 5, 25, 20, ":15+aaaaa")).is_err());
        // span 20, cs consumes 20 query but only 15 target
        assert!(parse_paf(&paf_line(100, 10, 30, '+', 200, 5, 25, 20, ":15-aaaaa")).is_err());
    }

    #[test]
    fn parse_paf_rejects_bad_cs() {
        assert!(parse_paf(&paf_line(100, 0, 10, '+', 100, 0, 10, 10, "")).is_err());
        assert!(parse_paf(&paf_line(100, 0, 10, '+', 100, 0, 10, 10, ":10X")).is_err());
    }
}
