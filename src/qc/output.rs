//! `qc` writers: VCF, BED, JSON summary, bedGraph and candidates TSV.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use anyhow::{Context, Result};
use serde::Serialize;

use super::{
    Alerts, Candidate, ClassFit, ContigCtx, Finding, Findings, NA_REASONS, QcParams, Region, Tier, VarType,
    Xlink, collect_alerts, collect_findings_counts, finding_rows, out_path,
};
use super::circularity::{CircCounts, CircResult, ContigCircularity};
use super::link::{HapRow, LinkageStats};
use crate::diff::Preset;
use crate::scan::ScanResult;
use crate::version::versioned;

pub(crate) struct Outputs<'a> {
    pub input: &'a Path,
    pub ref_path: &'a Path,
    pub prefix: &'a Path,
    pub preset: Preset,
    pub params: &'a QcParams,
    pub indel_min_af: f64,
    pub lowcov_abs: u32,
    /// Contig names in FASTA file order (the JSON `contigs` key order).
    pub fasta_order: &'a [String],
    pub scan: &'a ScanResult,
    pub refs: &'a [&'a [u8]],
    pub ctxs: &'a [ContigCtx],
    pub names: &'a [&'a str],
    pub cands: &'a [Candidate],
    pub regions: &'a [Region],
    pub findings: &'a [Finding],
    pub fits: &'a [ClassFit],
    /// The `linkage` JSON block.
    pub linkage: &'a LinkageStats,
    /// `P.haplotypes.tsv` rows.
    pub haplotypes: &'a [HapRow],
    /// Per-contig circularity evidence and warnings.
    pub circ: &'a CircResult,
    /// Per-contig count of `depth_step` peaks dropped as origin artefacts.
    pub depth_step_origin_dropped: &'a [u64],
    /// Cross-contig clip links (`clip_to=`/`link=` notes, `P.json` `xlinks`).
    pub xlinks: &'a [Xlink],
}

/// Counts per type for one tier.
#[derive(Clone, Debug, Default, Serialize)]
pub struct TypeCounts {
    #[serde(rename = "SNP")]
    pub snp: u64,
    #[serde(rename = "INS")]
    pub ins: u64,
    #[serde(rename = "DEL")]
    pub del: u64,
    pub total: u64,
}

/// A map serialized in insertion order (`serde_json` without `preserve_order`
/// sorts `Value` maps, so the contig order is kept here explicitly).
#[derive(Clone, Debug)]
pub struct OrderedMap<V>(pub Vec<(String, V)>);

impl<V: Serialize> Serialize for OrderedMap<V> {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.collect_map(self.0.iter().map(|(k, v)| (k, v)))
    }
}

/// One per-contig warning.
#[derive(Clone, Debug, Serialize)]
pub struct ContigWarning {
    pub code: String,
    /// Which thresholds triggered: `abs`/`rel` (union) or the NA reasons.
    pub reason: String,
    pub bp: u64,
    pub fraction: f64,
    pub message: String,
}

/// One top-level warning (per code, summed over contigs).
#[derive(Clone, Debug, Serialize)]
pub struct Warning {
    pub code: String,
    pub contigs: Vec<String>,
    pub bp: u64,
    pub fraction: f64,
    pub message: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct ContigSummary {
    pub length_bp: u64,
    pub circular: bool,
    /// Which source decided `circular`: `verdict`, `header` or `flag`.
    pub circular_by: &'static str,
    /// Sum of `M`/`=`/`X` bases on this contig.
    pub aligned_bases: u64,
    pub depth_mean: f64,
    pub depth_mode: u32,
    /// The mode ignores zero-depth positions (baseline of the covered part).
    pub depth_mode_excludes_zero: bool,
    /// Share of positions with depth > 0.
    pub covered_fraction: f64,
    /// depth mode ÷ depth mode of the longest contig.
    pub copy_number: Option<f64>,
    pub assessable_bp: u64,
    pub low_coverage_bp: u64,
    pub low_coverage_fraction: f64,
    pub alerts: Alerts,
    pub findings: Findings,
    pub warnings: Vec<ContigWarning>,
    /// Clipped records whose clip lies at a contig end (origin crossings on
    /// circular contigs; noted for step 2, no signal yet).
    pub edge_clip_records: u64,
    /// `depth_step` peaks dropped as origin artefacts on a wrapping contig.
    pub depth_step_origin_dropped: u64,
    /// Why the depth bands were skipped on this contig (`low_coverage` when the
    /// contig depth mode is below the absolute low-coverage threshold).
    pub depth_skipped: Option<&'static str>,
    /// Circularity evidence (verdict, junction quality, end reads).
    pub circularity: ContigCircularity,
}

#[derive(Clone, Debug, Serialize)]
pub struct RecordsSummary {
    /// Primary mapped records.
    pub primary: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub duplicate: u64,
    pub qcfail: u64,
    pub unmapped: u64,
    /// Counted records with a hard clip.
    pub hard_clipped: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ReadsSummary {
    /// Primary mapped + unmapped records (each read once).
    pub total: u64,
    pub mapped: u64,
    pub unmapped: u64,
    pub unmapped_fraction: f64,
    pub bases_total: u64,
    pub bases_mapped: u64,
    pub bases_unmapped: u64,
    /// Mean aligned reference span of a primary mapped record (bp); the depth
    /// taper near linear contig ends extends this far.
    pub mean_aligned_span: f64,
    pub records: RecordsSummary,
}

#[derive(Clone, Debug, Serialize)]
pub struct AssemblySummary {
    pub length_bp: u64,
    pub contigs: usize,
    pub assessable_bp: u64,
    pub assessable_fraction: f64,
    pub low_coverage_bp: u64,
    pub low_coverage_fraction: f64,
    /// Assessable bases that are not low-coverage by `abs`.
    pub confident_bp: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Timing {
    pub scan: f64,
    pub qc: f64,
    /// Wall time of the circularity pass (its own second scan of the contig ends).
    pub circularity: f64,
    pub total: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct Summary {
    pub analysis_type: &'static str,
    pub version: String,
    pub preset: String,
    pub inputs: BTreeMap<String, String>,
    pub indel_min_af: f64,
    pub lowcov_abs: u32,
    pub reads: ReadsSummary,
    pub assembly: AssemblySummary,
    /// Tiered calls: VCF sites + tiered BED regions, per tier.
    pub alerts: Alerts,
    /// Findings: tiered calls grouped into one problem location each.
    pub findings: Findings,
    /// Cross-contig clip links (the `P.json` `xlinks` count).
    pub xlinks: u64,
    /// Per code, summed over contigs.
    pub warnings: Vec<Warning>,
    /// −10·log10((red + orange site calls + 1) / (confident bases + 1)).
    pub qc_qv: f64,
    /// Contigs keyed by name, in assembly order.
    pub contigs: OrderedMap<ContigSummary>,
    /// Candidates (pre-filter passed) per type.
    pub candidates: TypeCounts,
    pub background: Vec<ClassFit>,
    pub params: QcParams,
    pub threads: usize,
    pub chunks: usize,
    pub timing_s: Timing,
    /// Linkage second pass: regions, pairs, groups and unlinked sites.
    pub linkage: LinkageStats,
    /// Per-verdict contig counts (circular / linear / unclear).
    pub circularity: CircCounts,
}

/// Low-coverage phrase from the reason bits: `below 20x` (abs),
/// `below 0.25x of contig depth` (rel), or both.
fn lowcov_phrase(bits: u8, low_abs: u32, low_rel: f64) -> String {
    match (bits & 1 != 0, bits & 2 != 0) {
        (true, true) => format!("below {low_abs}x or {low_rel}x of contig depth"),
        (true, false) => format!("below {low_abs}x"),
        (false, true) => format!("below {low_rel}x of contig depth"),
        (false, false) => String::new(),
    }
}

/// Parenthetical note after the low-coverage phrase, derived from the reason
/// bits: `abs` low coverage marks calls `LowCov`; `rel`-only low coverage is a
/// `LOWCOV` annotation and the calls themselves stay `PASS`.
fn lowcov_calls_note(bits: u8) -> &'static str {
    if bits & 1 != 0 {
        "calls marked LowCov"
    } else {
        "LOWCOV annotation; calls remain PASS"
    }
}

/// `qc_qv` = −10·log10((flagged + 1) / (confident + 1)), clamped to ≥ 0. The
/// pseudo-counts keep it finite when there are no confident bases; the clamp
/// guards against more flagged calls than confident bases (error probability
/// above 1).
fn qc_qv(flagged: u64, confident: u64) -> f64 {
    (-10.0 * ((flagged as f64 + 1.0) / (confident as f64 + 1.0)).log10()).max(0.0)
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.3}"),
        None => "NA".to_string(),
    }
}

/// Format a depth value for a `link=` note: no decimals when the value is a
/// whole count, one decimal otherwise.
fn fmt_depth(v: f64) -> String {
    if v.fract() == 0.0 {
        format!("{}", v as u64)
    } else {
        format!("{v:.1}")
    }
}

fn create(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?))
}

fn qual(ln_p: f64) -> f64 {
    (-10.0 * ln_p / std::f64::consts::LN_10).clamp(0.0, 999.0)
}

/// Remove duplicate entries from one finding's cross-contig link notes, keeping
/// the first occurrence (and order) of each. Two `clip_hotspot`s that link to
/// each other would otherwise write the same two notes on both findings.
pub(crate) fn dedup_notes(v: &mut Vec<String>) {
    let mut seen = std::collections::HashSet::new();
    v.retain(|n| seen.insert(n.clone()));
}

impl Outputs<'_> {
    pub fn write_all(&self, scan_s: f64, qc_s: f64, circ_s: f64, t0: Instant) -> Result<Summary> {
        crate::translate::ensure_parent(self.prefix)?;
        self.write_vcf()?;
        self.write_bed()?;
        self.write_findings()?;
        self.write_bedgraph()?;
        self.write_haplotypes()?;
        if self.params.candidates {
            self.write_candidates()?;
        }
        let mut s = self.summary(scan_s, qc_s, circ_s);
        s.timing_s.total = t0.elapsed().as_secs_f64();
        let p = out_path(self.prefix, ".json");
        let mut w = create(&p)?;
        serde_json::to_writer_pretty(&mut w, &versioned(&s))?;
        writeln!(w)?;
        w.flush()?;
        Ok(s)
    }

    fn write_vcf(&self) -> Result<()> {
        let p = out_path(self.prefix, ".vcf");
        let mut w = create(&p)?;
        writeln!(w, "##fileformat=VCFv4.2")?;
        writeln!(w, "##format_version={}", crate::version::FORMAT_VERSION)?;
        writeln!(w, "##source=cramhook qc {}", env!("CARGO_PKG_VERSION"))?;
        writeln!(w, "##reference={}", self.ref_path.display())?;
        writeln!(w, "##cramhook_preset={}", self.preset.minimap2_name())?;
        for (i, n) in self.names.iter().enumerate() {
            writeln!(w, "##contig=<ID={},length={}>", n, self.refs[i].len())?;
        }
        writeln!(w, "##FILTER=<ID=LowCov,Description=\"DP < {}: too few reads to trust the call\">", self.lowcov_abs)?;
        writeln!(w, "##INFO=<ID=TIER,Number=1,Type=String,Description=\"What the reads show: RED, ORANGE, YELLOW or GRAY (homopolymer indel, below yellow)\">")?;
        writeln!(w, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"SNP, INS or DEL\">")?;
        writeln!(w, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alt reads / DP\">")?;
        writeln!(w, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads with a base or a deletion at the position (unlike bcftools, deletions are counted); secondary, duplicate and QC-fail records excluded, supplementary included\">")?;
        writeln!(w, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Ref and alt read counts (indel ref = DP minus all indel alleles at the anchor)\">")?;
        writeln!(w, "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Ref and alt read counts, forward strand\">")?;
        writeln!(w, "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Ref and alt read counts, reverse strand\">")?;
        writeln!(w, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Ref-forward, ref-reverse, alt-forward, alt-reverse reads\">")?;
        writeln!(w, "##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Homopolymer run length (SNP: run containing the position; indel: run starting after the anchor)\">")?;
        writeln!(w, "##INFO=<ID=STR,Number=1,Type=Integer,Description=\"Tandem-repeat span (period 2-4, >= 3 copies) at the position / after the anchor; 0 if none\">")?;
        writeln!(w, "##INFO=<ID=CLASS,Number=1,Type=String,Description=\"Background class used for the test\">")?;
        writeln!(w, "##INFO=<ID=PV,Number=1,Type=Float,Description=\"Beta-binomial upper-tail p-value against the class background\">")?;
        writeln!(w, "##INFO=<ID=QV,Number=1,Type=Float,Description=\"Benjamini-Hochberg q-value (m = assessable positions of the class)\">")?;
        writeln!(w, "##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Share of MAPQ 0 reads in DP\">")?;
        writeln!(w, "##INFO=<ID=MIXED,Number=0,Type=Flag,Description=\"Alt fraction in the mixed band\">")?;
        writeln!(w, "##INFO=<ID=HPGRAY,Number=0,Type=Flag,Description=\"Gray: a minority/strand-biased yellow INS/DEL in a homopolymer >= --gray-hp-min (TIER=GRAY; not counted in red/orange/yellow)\">")?;
        writeln!(w, "##INFO=<ID=STRAND,Number=0,Type=Flag,Description=\"Strand-balance gate failed\">")?;
        writeln!(w, "##INFO=<ID=LOWCOV,Number=0,Type=Flag,Description=\"Low coverage: DP < {} (abs) or DP < {} x the contig depth mode (rel); only abs sets FILTER=LowCov\">", self.lowcov_abs, self.params.lowcov_rel)?;
        writeln!(w, "##INFO=<ID=LINKED,Number=1,Type=Integer,Description=\"Haplotype group id from the linkage second pass\">")?;
        writeln!(w, "##INFO=<ID=UNLINKED,Number=0,Type=Flag,Description=\"Tested against >= 2 partner sites within the linkage window and linked to none\">")?;
        if self.params.ins_group_min_len > 0 {
            writeln!(w, "##INFO=<ID=GROUPED,Number=1,Type=Integer,Description=\"Exact insertion alleles merged into this candidate (only when > 1)\">")?;
        }
        writeln!(w, "##INFO=<ID=ORIGIN,Number=0,Type=Flag,Description=\"Indel before the first reference base (right-padded: anchor base is the last REF/ALT base)\">")?;
        writeln!(w, "##INFO=<ID=REASON,Number=1,Type=String,Description=\"Tier rule that fired\">")?;
        writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (haploid)\">")?;
        writeln!(w, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Ref and alt read counts\">")?;
        writeln!(w, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth as INFO/DP\">")?;
        writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")?;
        for c in self.cands.iter().filter(|c| c.tier.tier != Tier::None) {
            // FILTER is confidence only: low `abs` depth is not trusted.
            let filter = if c.lowcov & 1 != 0 { "LowCov" } else { "PASS" };
            let (r, a) = (c.reff + c.refr, c.alt());
            let pv = c.ln_p.exp();
            let mut info = format!(
                "TIER={};TYPE={};AF={:.4};DP={};AD={},{};ADF={},{};ADR={},{};DP4={},{},{},{};HP={};STR={};CLASS={};PV={:.3e};QV={:.3e};MQ0F={:.3}",
                c.tier.tier.name().to_ascii_uppercase(),
                c.vtype.name(),
                c.af(),
                c.dp,
                r,
                a,
                c.reff,
                c.altf,
                c.refr,
                c.altr,
                c.reff,
                c.refr,
                c.altf,
                c.altr,
                c.hp,
                c.str_span,
                super::class_name(c.class),
                pv,
                c.q,
                c.mq0f
            );
            if c.tier.mixed {
                info.push_str(";MIXED");
            }
            if c.tier.tier == Tier::Gray {
                info.push_str(";HPGRAY");
            }
            if c.tier.strand {
                info.push_str(";STRAND");
            }
            if c.origin {
                info.push_str(";ORIGIN");
            }
            if c.lowcov != 0 {
                info.push_str(";LOWCOV");
            }
            if c.grouped > 1 {
                info.push_str(&format!(";GROUPED={}", c.grouped));
            }
            if let Some(g) = c.linked_group {
                info.push_str(&format!(";LINKED={g}"));
            }
            if c.unlinked {
                info.push_str(";UNLINKED");
            }
            info.push_str(";REASON=");
            info.push_str(c.tier.reason);
            writeln!(
                w,
                "{}\t{}\t.\t{}\t{}\t{:.1}\t{}\t{}\tGT:AD:DP\t1:{},{}:{}",
                self.names[c.tid],
                c.pos0 + 1,
                String::from_utf8_lossy(&c.ref_allele),
                String::from_utf8_lossy(&c.alt_allele),
                qual(c.ln_p),
                filter,
                info,
                r,
                a,
                c.dp
            )?;
        }
        w.flush()?;
        Ok(())
    }

    fn write_bed(&self) -> Result<()> {
        let p = out_path(self.prefix, ".bed");
        let mut w = create(&p)?;
        writeln!(w, "#contig\tstart\tend\tlabel\tscore\ttier\treason\tdepth_ratio\thighmapq_share\textra")?;
        for (i, r) in self.regions.iter().enumerate() {
            // Cross-contig clip links annotate the region's `extra` with
            // `clip_to`/`clip_to_reads`; a two-sided hotspot can carry one pair
            // per side, joined with `|` (parallel lists).
            let mut extra = r.extra.clone();
            let mut tos: Vec<String> = Vec::new();
            let mut reads: Vec<String> = Vec::new();
            for l in self.xlinks.iter().filter(|l| l.region == i) {
                tos.push(format!("{}:{}-{}", self.names[l.tid2], l.start2, l.end2));
                reads.push(l.reads.to_string());
            }
            if !tos.is_empty() {
                extra.push_str(&format!(";clip_to={};clip_to_reads={}", tos.join("|"), reads.join("|")));
            }
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.names[r.tid],
                r.start,
                r.end,
                r.label,
                r.score,
                r.tier,
                r.reason,
                fmt_opt(r.depth_ratio),
                fmt_opt(r.hiq_share),
                if extra.is_empty() { "." } else { &extra }
            )?;
        }
        w.flush()?;
        Ok(())
    }

    /// The `link=finding_<k>@<contig>:<s>-<e>,reads=<n>,f=<f>,partner_depth=<d>`
    /// note a finding that contains one end of a cross-contig link gets,
    /// pointing at the other end (`.@...` when no finding contains it).
    /// `partner_depth` is the depth at that other end as seen from this finding.
    fn xlink_note(
        other: Option<usize>,
        contig: &str,
        s: u64,
        e: u64,
        reads: u64,
        f: f64,
        partner_depth: f64,
    ) -> String {
        let k = other.map(|i| format!("finding_{}", i + 1)).unwrap_or_else(|| ".".to_string());
        format!(
            "link={k}@{contig}:{s}-{e},reads={reads},f={f:.4},partner_depth={}",
            fmt_depth(partner_depth)
        )
    }

    /// Per-finding `link=` notes from the cross-contig links (the hotspot side
    /// points at the partner and the partner finding at the hotspot). Findings
    /// touched by several links get them joined with `|`; identical notes (two
    /// hotspots linking to each other) are deduplicated, first occurrence first.
    fn xlink_notes(&self) -> Vec<Vec<String>> {
        let mut notes: Vec<Vec<String>> = vec![Vec::new(); self.findings.len()];
        for l in self.xlinks {
            if let Some(f1) = l.finding1 {
                // Hotspot finding: the other end is the partner breakpoint.
                notes[f1].push(Self::xlink_note(
                    l.finding2,
                    self.names[l.tid2],
                    l.start2,
                    l.end2,
                    l.reads,
                    l.f,
                    l.partner_depth,
                ));
            }
            if let Some(f2) = l.finding2 {
                // Partner finding: the other end is the hotspot side, whose depth
                // is the denominator of `f`.
                notes[f2].push(Self::xlink_note(
                    l.finding1,
                    self.names[l.tid1],
                    l.start1,
                    l.end1,
                    l.reads,
                    l.f,
                    l.hotspot_depth,
                ));
            }
        }
        for n in &mut notes {
            dedup_notes(n);
        }
        notes
    }

    fn write_findings(&self) -> Result<()> {
        let p = out_path(self.prefix, ".findings.bed");
        let mut w = create(&p)?;
        writeln!(w, "#contig\tstart\tend\tlabel\tscore\ttier\tlikely\textra")?;
        let notes = self.xlink_notes();
        for (n, f) in self.findings.iter().enumerate() {
            let members: Vec<String> =
                f.members.iter().map(|(l, c)| format!("{l}:{c}")).collect();
            let mut extra = format!(
                "evidence={};members={};lowcov={}",
                f.evidence_list(),
                members.join(","),
                f.lowcov as u8
            );
            // Cross-contig clip links: `;link=<other finding>@<contig>:<s>-<e>,...`.
            if !notes[n].is_empty() {
                extra.push_str(";");
                extra.push_str(&notes[n].join("|"));
            }
            // A finding that wraps the origin is written like a wrapped region
            // (see `depth_segment`/`snp_cluster`): two rows, both carrying the
            // same `finding_<n>` id (and tier/likely/extra), so every reader can
            // treat the id as the finding.
            let len = self.refs[f.tid].len() as u64;
            for (s, e) in finding_rows(f.start, f.end, len) {
                writeln!(
                    w,
                    "{}\t{}\t{}\tfinding_{}\t{}\t{}\t{}\t{}",
                    self.names[f.tid],
                    s,
                    e,
                    n + 1,
                    f.score,
                    f.tier,
                    f.likely,
                    extra
                )?;
            }
        }
        w.flush()?;
        Ok(())
    }

    fn write_bedgraph(&self) -> Result<()> {
        let p = out_path(self.prefix, ".depth.bedgraph");
        let mut w = create(&p)?;
        for (tid, c) in self.scan.contigs.iter().enumerate() {
            let n = c.len;
            let mut i = 0;
            while i < n {
                let s = i;
                let v = c.aligned_dp(s);
                while i < n && c.aligned_dp(i) == v {
                    i += 1;
                }
                writeln!(w, "{}\t{}\t{}\t{}", self.names[tid], s, i, v)?;
            }
        }
        w.flush()?;
        Ok(())
    }

    /// `P.haplotypes.tsv`: one row per (group, read) the linkage pass found.
    fn write_haplotypes(&self) -> Result<()> {
        let p = out_path(self.prefix, ".haplotypes.tsv");
        let mut w = create(&p)?;
        writeln!(w, "group\tcontig\tspan\tread\tstrand\talleles\tsa")?;
        for h in self.haplotypes {
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}", h.group, h.contig, h.span, h.read, h.strand, h.alleles, h.sa)?;
        }
        w.flush()?;
        Ok(())
    }

    fn write_candidates(&self) -> Result<()> {        let p = out_path(self.prefix, ".candidates.tsv");
        let mut w = create(&p)?;
        // The `grouped` column (exact alleles merged, 1 when not grouped) is
        // only written when INS grouping is on, so `--ins-group-min-len 0`
        // restores the exact-sequence schema byte for byte.
        let grouped_col = self.params.ins_group_min_len > 0;
        let mut head = String::from(
            "contig\tpos0\ttype\tlen\tref\talt\tdp\tdpf\tdpr\talt_reads\talt_fwd\talt_rev\tref_fwd\tref_rev\taf\thp\tstr\tmq0f\tclass\tpv\tqv\twilson_lb\tstrand_ok\tsignificant\ttier\treason\tflags\tlowcov",
        );
        if grouped_col {
            head.push_str("\tgrouped");
        }
        writeln!(w, "{head}")?;
        for c in self.cands {
            let mut flags = Vec::new();
            if c.tier.mixed {
                flags.push("MIXED");
            }
            if c.tier.tier == Tier::Gray {
                flags.push("HPGRAY");
            }
            if c.tier.strand {
                flags.push("STRAND");
            }
            if c.origin {
                flags.push("ORIGIN");
            }
            write!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.3}\t{}\t{:.4e}\t{:.4e}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.names[c.tid],
                c.pos0,
                c.vtype.name(),
                c.len,
                String::from_utf8_lossy(&c.ref_allele),
                String::from_utf8_lossy(&c.alt_allele),
                c.dp,
                c.dpf,
                c.dpr,
                c.alt(),
                c.altf,
                c.altr,
                c.reff,
                c.refr,
                c.af(),
                c.hp,
                c.str_span,
                c.mq0f,
                super::class_name(c.class),
                c.ln_p.exp(),
                c.q,
                c.wilson,
                c.strand_ok as u8,
                c.significant as u8,
                c.tier.tier.name(),
                c.tier.reason,
                if flags.is_empty() { ".".to_string() } else { flags.join(",") },
                { let r = super::lowcov_reason(c.lowcov); if r.is_empty() { ".".to_string() } else { r } }
            )?;
            if grouped_col {
                write!(w, "\t{}", c.grouped)?;
            }
            writeln!(w)?;
        }
        w.flush()?;
        Ok(())
    }

    fn summary(&self, scan_s: f64, qc_s: f64, circ_s: f64) -> Summary {
        let rc = &self.scan.records;
        let reads_total = rc.primary_mapped + rc.unmapped;
        let reads = ReadsSummary {
            total: reads_total,
            mapped: rc.primary_mapped,
            unmapped: rc.unmapped,
            unmapped_fraction: if reads_total == 0 { 0.0 } else { rc.unmapped as f64 / reads_total as f64 },
            bases_total: rc.mapped_bp + rc.unmapped_bp,
            bases_mapped: rc.mapped_bp,
            bases_unmapped: rc.unmapped_bp,
            mean_aligned_span: if rc.primary_mapped > 0 {
                rc.aligned_span as f64 / rc.primary_mapped as f64
            } else {
                0.0
            },
            records: RecordsSummary {
                primary: rc.primary_mapped,
                supplementary: rc.supplementary,
                secondary: rc.secondary,
                duplicate: rc.duplicate,
                qcfail: rc.qcfail,
                unmapped: rc.unmapped,
                hard_clipped: rc.hard_clipped_records,
            },
        };
        let longest = (0..self.ctxs.len()).max_by_key(|&i| (self.refs[i].len(), std::cmp::Reverse(i)));
        let base_mode = longest.map(|i| self.ctxs[i].mode).unwrap_or(0);
        let genome: u64 = self.refs.iter().map(|r| r.len() as u64).sum();
        let (per_alerts, alerts) = collect_alerts(self.ctxs.len(), self.cands, self.regions);
        let (per_findings, findings) = collect_findings_counts(self.ctxs.len(), self.findings);

        let mut contigs = Vec::new();
        let mut assessable = 0u64;
        let mut confident = 0u64;
        let mut low_total = 0u64;
        // (code, contig, bp, reason, message, reason bits, lowcov red calls)
        // per contig warning, in contig order.
        let mut cwarns: Vec<(String, String, u64, String, String, u8, u64)> = Vec::new();
        for (i, ctx) in self.ctxs.iter().enumerate() {
            let len = self.refs[i].len() as u64;
            let name = self.names[i].to_string();
            let a = ctx.na.iter().filter(|&&x| x == 0).count() as u64;
            let low_bp = ctx.lowcov.iter().filter(|&&x| x != 0).count() as u64;
            let abs_bp = ctx.lowcov.iter().filter(|&&x| x & 1 != 0).count() as u64;
            let na_bp = ctx.na.iter().filter(|&&x| x != 0).count() as u64;
            assessable += a;
            confident += a - abs_bp;
            low_total += low_bp;
            let too_much = |bp: u64| bp >= 1000 || (len > 0 && bp as f64 / len as f64 >= 0.05);
            let frac = |bp: u64| if len == 0 { 0.0 } else { bp as f64 / len as f64 };

            let mut warnings = Vec::new();
            if too_much(low_bp) {
                let mut bits = 0u8;
                for &x in &ctx.lowcov {
                    bits |= x;
                }
                let reason = super::lowcov_reason(bits);
                let fraction = frac(low_bp);
                let red = per_alerts[i].counts.lowcov.red;
                let mut message = format!(
                    "{name}: {:.0} % of bases {} ({})",
                    100.0 * fraction,
                    lowcov_phrase(bits, self.lowcov_abs, self.params.lowcov_rel),
                    lowcov_calls_note(bits)
                );
                if red > 0 {
                    message.push_str(&format!("; {red} low-coverage red calls"));
                }
                cwarns.push(("low_coverage".to_string(), name.clone(), low_bp, reason.clone(), message.clone(), bits, red));
                warnings.push(ContigWarning { code: "low_coverage".to_string(), reason, bp: low_bp, fraction, message });
            }
            if too_much(na_bp) {
                let reasons: Vec<&str> =
                    (1..NA_REASONS.len()).filter(|&r| ctx.na.iter().any(|&x| x as usize == r)).map(|r| NA_REASONS[r]).collect();
                let reason = reasons.join(",");
                let fraction = frac(na_bp);
                let message = format!("{name}: {:.0} % of bases not assessable", 100.0 * fraction);
                cwarns.push(("not_assessable".to_string(), name.clone(), na_bp, reason.clone(), message.clone(), 0, 0));
                warnings.push(ContigWarning { code: "not_assessable".to_string(), reason, bp: na_bp, fraction, message });
            }

            contigs.push((name, ContigSummary {
                length_bp: len,
                circular: ctx.circular,
                circular_by: ctx.circular_by,
                aligned_bases: self.scan.contigs[i].aligned_bases,
                depth_mean: (ctx.mean * 1000.0).round() / 1000.0,
                depth_mode: ctx.mode,
                depth_mode_excludes_zero: true,
                covered_fraction: if ctx.dp.is_empty() {
                    0.0
                } else {
                    ctx.dp.iter().filter(|&&d| d > 0).count() as f64 / ctx.dp.len() as f64
                },
                copy_number: if base_mode > 0 {
                    Some((ctx.mode as f64 / base_mode as f64 * 1000.0).round() / 1000.0)
                } else {
                    None
                },
                assessable_bp: a,
                low_coverage_bp: low_bp,
                low_coverage_fraction: frac(low_bp),
                alerts: per_alerts[i].clone(),
                findings: per_findings[i].clone(),
                warnings,
                edge_clip_records: self.scan.contigs[i].edge_clips,
                depth_step_origin_dropped: self.depth_step_origin_dropped.get(i).copied().unwrap_or(0),
                depth_skipped: ctx.depth_skipped,
                circularity: self.circ.per_contig[i].clone(),
            }));
        }

        // The JSON `contigs` follow FASTA order, not BAM @SQ order.
        let rank: std::collections::HashMap<&str, usize> =
            self.fasta_order.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        contigs.sort_by_key(|(n, _)| rank.get(n.as_str()).copied().unwrap_or(usize::MAX));

        // One top-level warning per code, summed over contigs (assembly order).
        let mut warnings = Vec::new();
        for code in ["low_coverage", "not_assessable"] {
            let mut names = Vec::new();
            let mut bp = 0u64;
            let mut bits = 0u8;
            let mut red = 0u64;
            for (c, name, b, _, _, bs, r) in &cwarns {
                if c == code {
                    names.push(name.clone());
                    bp += b;
                    bits |= bs;
                    red += r;
                }
            }
            if names.is_empty() {
                continue;
            }
            let fraction = if genome == 0 { 0.0 } else { bp as f64 / genome as f64 };
            let message = if code == "low_coverage" {
                let mut m = format!(
                    "{}: {:.0} % of the assembly {} ({})",
                    names.join(", "),
                    100.0 * fraction,
                    lowcov_phrase(bits, self.lowcov_abs, self.params.lowcov_rel),
                    lowcov_calls_note(bits)
                );
                if red > 0 {
                    m.push_str(&format!("; {red} low-coverage red calls"));
                }
                m
            } else {
                format!("{}: {:.0} % of the assembly not assessable", names.join(", "), 100.0 * fraction)
            };
            warnings.push(Warning { code: code.to_string(), contigs: names, bp, fraction, message });
        }

        // Circularity warnings (header/verdict mismatch, junction overlaps/gaps,
        // unclear contigs).
        for w in &self.circ.warnings {
            let fraction = if genome == 0 { 0.0 } else { w.bp as f64 / genome as f64 };
            warnings.push(Warning {
                code: w.code.to_string(),
                contigs: w.contigs.clone(),
                bp: w.bp,
                fraction,
                message: w.message.clone(),
            });
        }

        let mut cand_counts = TypeCounts::default();
        for c in self.cands {
            match c.vtype {
                VarType::Snp => cand_counts.snp += 1,
                VarType::Ins => cand_counts.ins += 1,
                VarType::Del => cand_counts.del += 1,
            }
            cand_counts.total += 1;
        }
        // qc_qv: red + orange site calls that are trusted (FILTER PASS) over
        // confident bases (assessable, not low by abs).
        let site_ro = self
            .cands
            .iter()
            .filter(|c| matches!(c.tier.tier, Tier::Red | Tier::Orange) && c.lowcov & 1 == 0)
            .count() as u64;
        let qc_qv = qc_qv(site_ro, confident);
        let mut inputs = BTreeMap::new();
        inputs.insert("in".to_string(), self.input.display().to_string());
        inputs.insert("ref".to_string(), self.ref_path.display().to_string());
        Summary {
            analysis_type: "cramhook qc",
            version: env!("CARGO_PKG_VERSION").to_string(),
            preset: self.preset.minimap2_name().to_string(),
            inputs,
            indel_min_af: self.indel_min_af,
            lowcov_abs: self.lowcov_abs,
            reads,
            assembly: AssemblySummary {
                length_bp: genome,
                contigs: self.ctxs.len(),
                assessable_bp: assessable,
                assessable_fraction: if genome == 0 { 0.0 } else { assessable as f64 / genome as f64 },
                low_coverage_bp: low_total,
                low_coverage_fraction: if genome == 0 { 0.0 } else { low_total as f64 / genome as f64 },
                confident_bp: confident,
            },
            alerts,
            findings,
            xlinks: self.xlinks.len() as u64,
            warnings,
            qc_qv,
            contigs: OrderedMap(contigs),
            candidates: cand_counts,
            background: self.fits.to_vec(),
            params: self.params.clone(),
            threads: self.params.threads,
            chunks: self.scan.chunks,
            timing_s: Timing { scan: scan_s, qc: qc_s, circularity: circ_s, total: 0.0 },
            linkage: self.linkage.clone(),
            circularity: self.circ.counts,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Several PASS red/orange calls can share one confident base, so the raw
    /// error rate can exceed 1 and the QV go negative; it must clamp to 0.
    #[test]
    fn qc_qv_clamped_to_zero_and_finite() {
        // 174 flagged calls over 100 confident bases: error rate 1.75 → raw QV −2.387.
        let q = qc_qv(174, 100);
        assert_eq!(q, 0.0, "{q}");
        assert!(q.is_finite());
        // More calls per base stays clamped, never negative or non-finite.
        assert_eq!(qc_qv(1_000, 1), 0.0);
        // Zero confident bases stays finite and non-negative.
        let z = qc_qv(0, 0);
        assert!(z.is_finite() && z >= 0.0, "{z}");
        // A clean assembly still yields a positive QV.
        let clean = qc_qv(0, 999);
        assert!(clean > 0.0 && clean.is_finite(), "{clean}");
    }

    /// Two hotspots that link to each other write the same two notes; the
    /// per-finding notes are deduplicated, first occurrence and order kept
    /// (reviewer issue 3).
    #[test]
    fn dedup_notes_keeps_first_occurrence() {
        let mut v = vec!["a".to_string(), "b".to_string(), "a".to_string(), "c".to_string(), "b".to_string()];
        dedup_notes(&mut v);
        assert_eq!(v, vec!["a".to_string(), "b".to_string(), "c".to_string()]);
        let mut empty: Vec<String> = Vec::new();
        dedup_notes(&mut empty);
        assert!(empty.is_empty());
        let mut one = vec!["x".to_string(), "x".to_string(), "x".to_string()];
        dedup_notes(&mut one);
        assert_eq!(one, vec!["x".to_string()]);
    }

    #[test]
    fn contigs_serialize_in_order() {
        let m = OrderedMap(vec![
            ("b".to_string(), 1u32),
            ("a".to_string(), 2u32),
            ("c".to_string(), 3u32),
        ]);
        let s = serde_json::to_string(&m).unwrap();
        let ib = s.find("\"b\"").unwrap();
        let ia = s.find("\"a\"").unwrap();
        let ic = s.find("\"c\"").unwrap();
        assert!(ib < ia && ia < ic, "{s}");
    }
}
