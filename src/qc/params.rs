//! `qc` tunable parameters, param enums and validators.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use anyhow::{Result, bail};

use crate::diff::Preset;

/// Tunable `qc` parameters (starting values; see the spec).
#[derive(Clone, Debug, serde::Serialize)]
pub struct QcParams {
    /// Read preset; `None` = from the minimap2 `@PG` line.
    #[serde(skip)]
    pub preset: Option<Preset>,
    #[serde(skip)]
    pub threads: usize,
    /// Contigs treated as linear (from `--linear`); all others wrap iff the
    /// circularity verdict is `circular`/`unclear` or the FASTA header says
    /// `circular=true`.
    pub linear: Vec<String>,
    /// Write `P.candidates.tsv`.
    pub candidates: bool,
    /// Not assessable below this depth.
    pub min_depth: u32,
    /// Not assessable above this MAPQ0 share.
    pub max_mq0_frac: f64,
    /// Candidate pre-filter: minimum alt reads.
    pub min_alt: u32,
    /// Candidate pre-filter: minimum SNP alt fraction.
    pub snp_min_af: f64,
    /// Candidate pre-filter: minimum indel alt fraction (`None` = preset default).
    pub indel_min_af: Option<f64>,
    /// Low coverage: assessable positions below this absolute depth are
    /// low-coverage by `abs` (`None` = preset default).
    pub lowcov_abs: Option<u32>,
    /// Low coverage: assessable positions below this fraction of the contig
    /// depth mode are low-coverage by `rel`.
    pub lowcov_rel: f64,
    /// BED low-coverage regions: merge runs separated by gaps below this.
    pub lowcov_merge: u32,
    /// BED low-coverage regions: drop regions shorter than this.
    pub lowcov_min_len: u32,
    /// Background fit uses positions with alt fraction below this.
    pub fit_max_af: f64,
    /// Indel classes with fewer fit sites fall back to the pooled INS/DEL fit.
    pub min_fit_sites: u64,
    /// BH q-value threshold.
    pub fdr: f64,
    /// Strand balance: each strand needs at least this depth.
    pub strand_min_dp: u32,
    /// Strand balance: weaker-strand AF >= ratio * stronger-strand AF.
    pub strand_ratio: f64,
    /// Red: minimum alt reads.
    pub red_min_alt: u32,
    /// Red: Wilson 95 % lower bound above this.
    pub red_min_wilson: f64,
    /// Orange (MIXED) band lower bound (upper bound = 0.5; a site with AF > 0.5
    /// that misses the red bar is `weak_majority` yellow).
    pub mixed_min_af: f64,
    /// Strand gate also applies to yellow minority sites (AF < mixed_min_af):
    /// failing it gives tier none (`strand_bias`). `false` restores the old
    /// rule (yellow with a STRAND flag).
    pub yellow_strand_gate: bool,
    /// Gray (below yellow): a minority/strand-biased yellow call, of type INS
    /// or DEL, whose homopolymer run (INFO `HP`) is at least `gray_hp_min` is
    /// reported as `GRAY` (`HPGRAY`). Red/orange calls, SNPs and `weak_majority`
    /// (AF > 0.5) calls are never gray. Gray calls stay in the VCF but are not
    /// finding members and are not counted in the red/orange/yellow alert or
    /// tier counts.
    pub gray_homopolymer: bool,
    /// Gray: minimum homopolymer run length (INFO `HP`, the indel's forward run
    /// after the anchor) for a minority/strand-biased yellow INS/DEL.
    pub gray_hp_min: u32,
    pub large_merge_gap: u32,
    pub large_min: u32,
    pub large_piece_min: u32,
    pub large_min_support: u32,
    pub large_min_ratio: f64,
    /// INS events within this distance are clustered.
    pub large_ins_window: u32,
    pub cluster_window: u32,
    pub cluster_min_sites: usize,
    pub cluster_min_af: f64,
    pub cluster_min_alt: u32,
    /// Cluster site refinement A: the SNP also needs a beta-binomial upper-tail
    /// p-value <= this (nominal, no FDR); default 1e-3, 1 = off (old rule).
    pub cluster_site_p: f64,
    /// Cluster refinement B: keep a window only if Poisson P(N >= k | λ·W) <=
    /// this / (assessable bp / W) (Bonferroni), λ = sites / assessable bp per
    /// contig; `None` = off.
    pub cluster_density_p: Option<f64>,
    /// Cluster site refinement C: strand requirement for cluster sites;
    /// default `ratio`, `off` = old rule.
    pub cluster_strand: ClusterStrand,
    /// Group similar INS alleles: only insertions at least this long take part;
    /// 0 disables grouping (exact-sequence behaviour).
    pub ins_group_min_len: u32,
    /// Group similar INS alleles: length tolerance as a fraction of the
    /// representative length (at least 2 bp).
    pub ins_group_len_frac: f64,
    /// Group similar INS alleles: maximum edit distance to the representative
    /// sequence as a fraction of its length.
    pub ins_group_max_diff: f64,
    /// Clip hotspots: minimum clip length to count (bp).
    pub clip_min_len: u32,
    /// Clip hotspots: single-linkage merge gap between clip positions (bp).
    pub clip_merge: u32,
    /// Clip hotspots: minimum cluster reads for a candidate.
    pub clip_min_reads: u32,
    /// Clip hotspots: candidate (low bar) ratio floor (cluster count /
    /// aligned-side depth); a cluster below it is not reported.
    pub clip_yellow: f64,
    /// Clip hotspots: high significance bar (`s = -log10 p`); at or above it a
    /// call is orange (f <= 0.5) or red (f > 0.5), below it yellow.
    pub clip_high_s: f64,
    /// Clip hotspots: distance within which a left and a right cluster pair (bp).
    pub clip_pair_dist: u32,
    /// Clip hotspots: skip clusters within this of a contig end (bp).
    pub clip_edge: u32,
    /// Depth bands: window length in bp.
    pub depth_window: u32,
    /// Depth bands: `low` band ratio ceiling.
    pub depth_low: f64,
    /// Depth bands: `half` band ratio floor.
    pub depth_half_min: f64,
    /// Depth bands: `half` band ratio ceiling.
    pub depth_half_max: f64,
    /// Depth bands: `half` band high-MAPQ share ceiling.
    pub depth_half_max_mapq: f64,
    /// Depth bands: `high` band ratio floor.
    pub depth_high: f64,
    /// Depth bands: `high` band high-MAPQ share floor.
    pub depth_high_min_mapq: f64,
    /// Depth bands: gap (in windows) allowed inside a merged region.
    pub depth_gap_windows: usize,
    /// Depth bands: drop regions shorter than this (bp).
    pub depth_min_len: u64,
    /// Depth bands: local-baseline span (bp); the baseline is the median of the
    /// window means within ±this of the window (whole contig when shorter than
    /// 2× the span).
    pub depth_baseline_span: u32,
    /// Depth steps: contrast window `w` in bp (one side).
    pub step_window: u32,
    /// Depth steps: local median/MAD neighbourhood half-width (bp).
    pub step_span: u32,
    /// Depth steps: non-maximum suppression half-width (bp).
    pub step_nms: u32,
    /// Depth steps: yellow `|z|` threshold (low bar).
    pub step_yellow: f64,
    /// Depth steps: red `|z|` threshold (high bar).
    pub step_red: f64,
    /// Findings: single-linkage merge gap between member intervals (bp).
    pub finding_merge: u32,
    /// Findings: a site call joins a finding that has a region-type member only
    /// within this distance (bp) of one of those regions; region-region and
    /// site-site joins still use `finding_merge`.
    pub finding_site_join: u32,
    /// Findings: link the two ends of a `split_junction` cluster and of a
    /// `depth_segment` (bp apart at most, the shorter way round on a circular
    /// contig) so their grouped calls pool into one finding; 0 disables linking
    /// (today's behaviour).
    pub finding_link_max: u32,
    /// Linkage: sub-threshold candidate sites down to this alt fraction.
    pub link_min_af: f64,
    /// Linkage: group candidate sites within this distance (bp).
    pub link_window: u32,
    /// Linkage: minimum alt/alt reads for a linked pair.
    pub link_min_reads: u32,
    /// Linkage: minimum concordance a/(a+b) and a/(a+c).
    pub link_min_conc: f64,
    /// Linkage: maximum one-sided Fisher exact p for a linked pair.
    pub link_p: f64,
    /// Linkage: maximum one-sided Fisher exact p for each strand table (the
    /// strand check only guards against strand-specific errors; significance
    /// comes from the pooled test).
    pub link_strand_p: f64,
    /// Linkage conflict: minimum reads carrying the alt at X or Y (a+b+c).
    pub link_conflict_min_reads: u32,
    /// Linkage conflict: maximum alt/alt share a/(a+b+c) below which the pair
    /// is a decisive incompatibility.
    pub link_conflict_max_share: f64,
    /// Linkage: what to do with tested-but-unlinked minority sites.
    pub unlinked_action: UnlinkedAction,
    /// Circularity: an alignment within this many bp of a contig end reaches it.
    pub circ_end_tol: u32,
    /// Circularity: minimum spanning reads for a `circular` verdict.
    pub circ_min_reads: u32,
    /// Circularity: minimum spanning / min(end, start) reads for `circular`.
    pub circ_min_frac: f64,
    /// Junction reads: minimum |A - B| between the two breakpoints (bp).
    pub junction_min_dist: u32,
    /// Junction reads: maximum read-coordinate gap/overlap between the two
    /// segments (bp); extra unaligned sequence is an insertion, not a junction.
    pub junction_max_gap: u32,
    /// Junction reads: cluster breakpoints within this of each other (bp).
    pub junction_tol: u32,
    /// Junction reads: skip a junction with either breakpoint within this of a
    /// contig end, on a circular contig (bp; the wrap is a special case).
    pub junction_edge_excl: u32,
    /// Junction reads: minimum distinct reads for a candidate.
    pub junction_min_reads: u32,
    /// Junction reads: reporting gate (support / mean breakpoint depth); every
    /// reported junction is orange (f <= 0.5) or red (f > 0.5).
    pub junction_yellow: f64,
    /// Cross-contig clip links: maximum read-coordinate gap *or overlap*
    /// between the two alignments of one read that join across contigs (bp).
    /// Larger than `junction_max_gap`: a repeat-mediated join legitimately
    /// shares a few hundred bp of read between the two alignments.
    pub xlink_max_gap: u32,
    #[serde(skip)]
    pub chunk_size: Option<u64>,
}

impl Default for QcParams {
    fn default() -> Self {
        QcParams {
            preset: None,
            threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8),
            linear: Vec::new(),
            candidates: false,
            min_depth: 5,
            max_mq0_frac: 0.5,
            min_alt: 2,
            snp_min_af: 0.08,
            indel_min_af: None,
            lowcov_abs: None,
            lowcov_rel: 0.25,
            lowcov_merge: 50,
            lowcov_min_len: 100,
            fit_max_af: 0.3,
            min_fit_sites: 1000,
            fdr: 0.01,
            strand_min_dp: 3,
            strand_ratio: 1.0 / 3.0,
            red_min_alt: 5,
            red_min_wilson: 0.50,
            mixed_min_af: 0.40,
            yellow_strand_gate: true,
            gray_homopolymer: true,
            gray_hp_min: 5,
            large_merge_gap: 50,
            large_min: 50,
            large_piece_min: 10,
            large_min_support: 3,
            large_min_ratio: 0.3,
            large_ins_window: 100,
            cluster_window: 500,
            cluster_min_sites: 5,
            cluster_min_af: 0.10,
            cluster_min_alt: 2,
            cluster_site_p: 1e-3,
            cluster_density_p: None,
            cluster_strand: ClusterStrand::Ratio,
            ins_group_min_len: 10,
            ins_group_len_frac: 0.10,
            ins_group_max_diff: 0.20,
            clip_min_len: 50,
            clip_merge: 50,
            clip_min_reads: 3,
            clip_yellow: 0.15,
            clip_high_s: 100.0,
            clip_pair_dist: 50,
            clip_edge: 100,
            depth_window: 500,
            depth_low: 0.25,
            depth_half_min: 0.35,
            depth_half_max: 0.65,
            depth_half_max_mapq: 0.25,
            depth_high: 1.75,
            depth_high_min_mapq: 0.75,
            depth_gap_windows: 1,
            depth_min_len: 1000,
            depth_baseline_span: 50_000,
            step_window: 100,
            step_span: 25_000,
            step_nms: 500,
            step_yellow: 8.0,
            step_red: 20.0,
            finding_merge: 1000,
            finding_site_join: 50,
            finding_link_max: 100_000,
            link_min_af: 0.05,
            link_window: 5000,
            link_min_reads: 3,
            link_min_conc: 0.8,
            link_p: 1e-6,
            link_strand_p: 0.01,
            link_conflict_min_reads: 6,
            link_conflict_max_share: 0.3,
            unlinked_action: UnlinkedAction::Drop,
            circ_end_tol: 50,
            circ_min_reads: 3,
            circ_min_frac: 0.3,
            junction_min_dist: 1000,
            junction_max_gap: 200,
            junction_tol: 200,
            junction_edge_excl: 2000,
            junction_min_reads: 5,
            junction_yellow: 0.02,
            xlink_max_gap: 1000,
            chunk_size: None,
        }
    }
}

/// Strand requirement for minority-SNP cluster sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum ClusterStrand {
    /// No strand requirement.
    Off,
    /// Alt reads on both strands.
    Both,
    /// The site strand gate (each strand depth >= strand_min_dp, weaker-strand
    /// AF >= strand_ratio × stronger-strand AF).
    Ratio,
}

/// What the linkage pass does with a tested-but-unlinked minority site.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum UnlinkedAction {
    /// Keep the site's tier and add INFO `UNLINKED`.
    Flag,
    /// Drop the site to tier none (reason `unlinked`; kept in candidates.tsv).
    Drop,
}

/// Default indel pre-filter alt fraction per preset.
pub fn default_indel_min_af(p: Preset) -> f64 {
    match p {
        Preset::MapHifi => 0.08,
        Preset::LrHq => 0.15,
    }
}

/// Default low-coverage absolute depth threshold per preset.
pub fn default_lowcov_abs(p: Preset) -> u32 {
    match p {
        Preset::MapHifi => 10,
        Preset::LrHq => 20,
    }
}

/// `--lowcov-rel` must be a finite, non-negative number.
pub fn validate_lowcov_rel(v: f64) -> Result<()> {
    if !v.is_finite() || v < 0.0 {
        bail!("--lowcov-rel must be a finite, non-negative number, got {v}");
    }
    Ok(())
}

/// `--ins-group-len-frac` and `--ins-group-max-diff` must be finite,
/// non-negative fractions.
pub fn validate_ins_group_frac(v: f64) -> Result<()> {
    if !v.is_finite() || v < 0.0 {
        bail!("--ins-group-len-frac/--ins-group-max-diff must be finite, non-negative numbers, got {v}");
    }
    Ok(())
}

/// Variant type; the order is the VCF sort order at one position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum VarType {
    Snp,
    Del,
    Ins,
}

impl VarType {
    pub fn name(self) -> &'static str {
        match self {
            VarType::Snp => "SNP",
            VarType::Del => "DEL",
            VarType::Ins => "INS",
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lowcov_rel_validation() {
        assert!(validate_lowcov_rel(0.0).is_ok());
        assert!(validate_lowcov_rel(0.25).is_ok());
        assert!(validate_lowcov_rel(-0.1).is_err());
        assert!(validate_lowcov_rel(f64::NAN).is_err());
        assert!(validate_lowcov_rel(f64::INFINITY).is_err());
    }
    #[test]
    fn ins_group_frac_validation() {
        assert!(validate_ins_group_frac(0.0).is_ok());
        assert!(validate_ins_group_frac(0.2).is_ok());
        assert!(validate_ins_group_frac(-0.1).is_err());
        assert!(validate_ins_group_frac(f64::NAN).is_err());
        assert!(validate_ins_group_frac(f64::INFINITY).is_err());
    }
}
