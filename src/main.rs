//! cramhook command line interface.

use std::path::PathBuf;

use anyhow::Result;
use clap::{Args, Parser, Subcommand};

use cramhook::{Preset, note};

// The static musl build swaps musl's slow allocator for mimalloc; glibc builds are unchanged.
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// `--max-de`: must be a finite, non-negative number.
fn parse_max_de(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("invalid number {s:?}"))?;
    if !v.is_finite() || v < 0.0 {
        return Err("must be finite and >= 0".to_string());
    }
    Ok(v)
}

/// `--max-s2-ratio`: must be a finite, non-negative number.
fn parse_max_s2_ratio(s: &str) -> Result<f64, String> {
    parse_max_de(s)
}

/// `--min-input-mb`: must be a finite, non-negative number.
fn parse_min_input_mb(s: &str) -> Result<f64, String> {
    parse_max_de(s)
}

/// Default `--threads`: the number of available cores (8 if that query fails).
fn default_threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8)
}

/// `--max-zdrop`: must be a non-negative integer.
fn parse_max_zdrop(s: &str) -> Result<i64, String> {
    let v: i64 = s.parse().map_err(|_| format!("invalid integer {s:?}"))?;
    if v < 0 {
        return Err("must be >= 0".to_string());
    }
    Ok(v)
}

#[derive(Parser)]
#[command(name = "cramhook", version, about = "Lift reads from an old assembly to a new one")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Full pipeline: translate clean reads, remap the rest, and write a sorted CRAM.
    Lift(LiftArgs),
    /// Compare two assemblies; write identity blocks, taint and novel regions.
    Diff(DiffArgs),
    /// Translate old-alignment records to the new assembly (debug/intermediate).
    Translate(TranslateArgs),
    /// Assembly QC: flag sites and regions where the reads disagree with the assembly.
    Qc(QcArgs),
    /// Extract the reads that disagree with the assembly (one FASTQ + reasons table).
    Fish(FishArgs),
}

#[derive(Args)]
struct QcArgs {
    /// Input BAM/CRAM (coordinate-sorted, indexed) mapped to the assembly.
    #[arg(long)]
    r#in: PathBuf,
    /// Assembly FASTA (the alignment reference).
    #[arg(long)]
    r#ref: PathBuf,
    /// Output prefix; writes `<prefix>.{vcf,bed,json,depth.bedgraph,gc.bedgraph}` (+ `.candidates.tsv`).
    #[arg(long)]
    out_prefix: PathBuf,
    /// Read preset (thresholds); default: from the minimap2 @PG line.
    #[arg(long)]
    preset: Option<Preset>,
    /// Threads (default: all available cores).
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// Contigs treated as linear (comma-separated), overriding the circularity
    /// verdict and a FASTA `circular=true` header. All other contigs wrap iff
    /// the verdict is `circular`/`unclear` or the header says `circular=true`.
    #[arg(long, value_delimiter = ',')]
    linear: Vec<String>,
    /// Also write `<prefix>.candidates.tsv` (every pre-filtered candidate with all statistics).
    #[arg(long)]
    candidates: bool,
    /// Positions below this depth are not assessable.
    #[arg(long, default_value_t = 5)]
    min_depth: u32,
    /// Positions above this MAPQ0 share are not assessable.
    #[arg(long, default_value_t = 0.5)]
    max_mq0_frac: f64,
    /// Candidate pre-filter: minimum alt reads.
    #[arg(long, default_value_t = 2)]
    min_alt: u32,
    /// Candidate pre-filter: minimum SNP alt fraction.
    #[arg(long, default_value_t = 0.08)]
    snp_min_af: f64,
    /// Candidate pre-filter: minimum indel alt fraction (default: map-hifi 0.08, lr:hq 0.15).
    #[arg(long)]
    indel_min_af: Option<f64>,
    /// Low coverage: assessable positions below this absolute depth (default: map-hifi 10, lr:hq 20).
    #[arg(long)]
    lowcov_abs: Option<u32>,
    /// Low coverage: assessable positions below this fraction of the contig depth mode.
    #[arg(long, default_value_t = 0.25)]
    lowcov_rel: f64,
    /// BED low-coverage regions: merge runs separated by gaps below this many bp.
    #[arg(long, default_value_t = 50)]
    lowcov_merge: u32,
    /// BED low-coverage regions: drop regions shorter than this many bp.
    #[arg(long, default_value_t = 100)]
    lowcov_min_len: u32,
    /// Benjamini-Hochberg q-value threshold for significance.
    #[arg(long, default_value_t = 0.01)]
    fdr: f64,
    /// Strand balance: weaker-strand AF must be >= this fraction of the stronger-strand AF.
    #[arg(long, default_value_t = 1.0 / 3.0)]
    strand_ratio: f64,
    /// Red: minimum alt reads.
    #[arg(long, default_value_t = 5)]
    red_min_alt: u32,
    /// Do not gate yellow minority sites (AF < 0.40) on strand balance: keep
    /// them yellow with a STRAND flag (old rule, for measurement).
    #[arg(long)]
    no_yellow_strand_gate: bool,
    /// Gray: minimum homopolymer run length (INFO `HP`) for a minority/
    /// strand-biased yellow INS/DEL to be reported as GRAY instead.
    #[arg(long, default_value_t = 5)]
    gray_hp_min: u32,
    /// Do not mark homopolymer yellow indels gray: keep them yellow.
    #[arg(long)]
    no_gray_homopolymer: bool,
    /// Large indels: merge same-type indels of one read separated by at most this many bp.
    #[arg(long, default_value_t = 50)]
    large_merge_gap: u32,
    /// Large indels: minimum merged size in bp.
    #[arg(long, default_value_t = 50)]
    large_min: u32,
    /// Large indels: only indels of at least this size take part in the merge.
    #[arg(long, default_value_t = 10)]
    large_piece_min: u32,
    /// Large indels: minimum supporting reads.
    #[arg(long, default_value_t = 3)]
    large_min_support: u32,
    /// Large indels: minimum support / depth ratio.
    #[arg(long, default_value_t = 0.3)]
    large_min_ratio: f64,
    /// SNP clusters: window length in bp.
    #[arg(long, default_value_t = 500)]
    cluster_window: u32,
    /// SNP clusters: minimum sites per window.
    #[arg(long, default_value_t = 5)]
    cluster_min_sites: usize,
    /// SNP clusters, refinement A: a site also needs a beta-binomial p <= this
    /// (nominal, no FDR); 1 = off.
    #[arg(long, default_value_t = 1e-3)]
    cluster_site_p: f64,
    /// SNP clusters, refinement B: Poisson window density test at this
    /// family-wise level (Bonferroni over windows); default off.
    #[arg(long)]
    cluster_density_p: Option<f64>,
    /// SNP clusters, refinement C: strand requirement for cluster sites
    /// (`off` with `--cluster-site-p 1` restores the old rule).
    #[arg(long, value_enum, default_value_t = cramhook::qc::ClusterStrand::Ratio)]
    cluster_strand: cramhook::qc::ClusterStrand,
    /// Group similar insertions of the same site: only insertions at least this
    /// long take part; 0 restores exact-sequence counting.
    #[arg(long, default_value_t = 10)]
    ins_group_min_len: u32,
    /// Group similar insertions: length tolerance as a fraction of the
    /// representative length (at least 2 bp).
    #[arg(long, default_value_t = 0.10)]
    ins_group_len_frac: f64,
    /// Group similar insertions: maximum edit distance to the representative
    /// sequence as a fraction of its length.
    #[arg(long, default_value_t = 0.20)]
    ins_group_max_diff: f64,
    /// Clip hotspots: only clips of at least this length count.
    #[arg(long, default_value_t = 50)]
    clip_min_len: u32,
    /// Clip hotspots: merge clip positions within this gap.
    #[arg(long, default_value_t = 50)]
    clip_merge: u32,
    /// Clip hotspots: minimum cluster reads.
    #[arg(long, default_value_t = 3)]
    clip_min_reads: u32,
    /// Clip hotspots: candidate (low bar) ratio floor (cluster count /
    /// aligned-side depth); a cluster below it is not reported.
    #[arg(long, default_value_t = 0.15)]
    clip_yellow: f64,
    /// Clip hotspots: high significance bar (`s = -log10 p`); at or above it a
    /// call is orange (f <= 0.5) or red (f > 0.5), below it yellow.
    #[arg(long, default_value_t = 100.0)]
    clip_high_s: f64,
    /// Clip hotspots: pair left/right clusters within this distance.
    #[arg(long, default_value_t = 50)]
    clip_pair_dist: u32,
    /// Clip hotspots: skip clusters within this of a contig end.
    #[arg(long, default_value_t = 100)]
    clip_edge: u32,
    /// GC track: window length in bp; 0 disables `<prefix>.gc.bedgraph`.
    #[arg(long, default_value_t = 100)]
    gc_window: u32,
    /// Depth bands: window length in bp.
    #[arg(long, default_value_t = 500)]
    depth_window: u32,
    /// Depth bands: `low` band ratio ceiling (red).
    #[arg(long, default_value_t = 0.25)]
    depth_low: f64,
    /// Depth bands: `half` band ratio floor.
    #[arg(long, default_value_t = 0.35)]
    depth_half_min: f64,
    /// Depth bands: `half` band ratio ceiling.
    #[arg(long, default_value_t = 0.65)]
    depth_half_max: f64,
    /// Depth bands: `half` band high-MAPQ share ceiling.
    #[arg(long, default_value_t = 0.25)]
    depth_half_max_mapq: f64,
    /// Depth bands: `high` band ratio floor.
    #[arg(long, default_value_t = 1.75)]
    depth_high: f64,
    /// Depth bands: `high` band high-MAPQ share floor.
    #[arg(long, default_value_t = 0.75)]
    depth_high_min_mapq: f64,
    /// Depth bands: gap (windows of another or no band) allowed inside a region.
    #[arg(long, default_value_t = 1)]
    depth_gap_windows: usize,
    /// Depth bands: drop regions shorter than this many bp.
    #[arg(long, default_value_t = 1000)]
    depth_min_len: u64,
    /// Depth bands: local-baseline span in bp (median of window means within
    /// ±this; whole contig when shorter than 2× the span).
    #[arg(long, default_value_t = 50_000)]
    depth_baseline_span: u32,
    /// Depth steps: contrast window `w` in bp (one side).
    #[arg(long, default_value_t = 100)]
    step_window: u32,
    /// Depth steps: local median/MAD neighbourhood half-width in bp.
    #[arg(long, default_value_t = 25_000)]
    step_span: u32,
    /// Depth steps: non-maximum suppression half-width in bp.
    #[arg(long, default_value_t = 500)]
    step_nms: u32,
    /// Depth steps: yellow `|z|` threshold (low bar).
    #[arg(long, default_value_t = 8.0)]
    step_yellow: f64,
    /// Depth steps: red `|z|` threshold (high bar).
    #[arg(long, default_value_t = 20.0)]
    step_red: f64,
    /// Findings: merge member intervals within this many bp into one finding.
    #[arg(long, default_value_t = 1000)]
    finding_merge: u32,
    /// Findings: a site call joins a finding with a region-type member only
    /// within this many bp of one of those regions.
    #[arg(long, default_value_t = 50)]
    finding_site_join: u32,
    /// Findings: link the two ends of a `split_junction` cluster and of a
    /// `depth_segment` within this many bp (shorter way round on a circular
    /// contig), pooling their calls into one finding; 0 disables linking.
    #[arg(long, default_value_t = 100_000)]
    finding_link_max: u32,
    /// Linkage: sub-threshold candidate SNP sites down to this alt fraction.
    #[arg(long, default_value_t = 0.05)]
    link_min_af: f64,
    /// Linkage: group candidate sites within this many bp.
    #[arg(long, default_value_t = 5000)]
    link_window: u32,
    /// Linkage: minimum alt/alt reads for a linked pair.
    #[arg(long, default_value_t = 3)]
    link_min_reads: u32,
    /// Linkage: minimum concordance a/(a+b) and a/(a+c) for a linked pair.
    #[arg(long, default_value_t = 0.8)]
    link_min_conc: f64,
    /// Linkage: maximum one-sided Fisher exact p for a linked pair.
    #[arg(long, default_value_t = 1e-6)]
    link_p: f64,
    /// Linkage: maximum one-sided Fisher exact p for each strand table
    /// (guards strand-specific errors; significance comes from --link-p).
    #[arg(long, default_value_t = 0.01)]
    link_strand_p: f64,
    /// Linkage conflict: minimum reads carrying the alt at X or Y (a+b+c).
    #[arg(long, default_value_t = 6)]
    link_conflict_min_reads: u32,
    /// Linkage conflict: maximum alt/alt share a/(a+b+c) below which a pair is
    /// a decisive incompatibility (two haplotypes).
    #[arg(long, default_value_t = 0.3)]
    link_conflict_max_share: f64,
    /// Linkage: what to do with tested-but-unlinked minority sites (flag|drop).
    #[arg(long, value_enum, default_value_t = cramhook::qc::UnlinkedAction::Drop)]
    unlinked_action: cramhook::qc::UnlinkedAction,
    /// Circularity: an alignment within this many bp of a contig end reaches it.
    #[arg(long, default_value_t = 50)]
    circ_end_tol: u32,
    /// Circularity: minimum spanning reads for a `circular` verdict.
    #[arg(long, default_value_t = 3)]
    circ_min_reads: u32,
    /// Circularity: minimum spanning / min(end, start) reads for `circular`.
    #[arg(long, default_value_t = 0.3)]
    circ_min_frac: f64,
    /// Junction reads: minimum |A - B| between the two breakpoints.
    #[arg(long, default_value_t = 1000)]
    junction_min_dist: u32,
    /// Junction reads: maximum read-coordinate gap/overlap between the two
    /// segments; extra unaligned sequence is an insertion, not a junction.
    #[arg(long, default_value_t = 200)]
    junction_max_gap: u32,
    /// Junction reads: cluster breakpoints within this of each other (ONT
    /// breakpoints jitter; 200 recovers the analysis' support without merging
    /// distinct junctions).
    #[arg(long, default_value_t = 200)]
    junction_tol: u32,
    /// Junction reads: skip a junction with either breakpoint within this of a
    /// contig end, on a circular contig (the wrap is the special case of both).
    #[arg(long, default_value_t = 2000)]
    junction_edge_excl: u32,
    /// Junction reads: minimum distinct reads for a candidate.
    #[arg(long, default_value_t = 5)]
    junction_min_reads: u32,
    /// Junction reads: reporting gate (support / mean breakpoint depth); every
    /// reported junction is orange (f <= 0.5) or red (f > 0.5).
    #[arg(long, default_value_t = 0.02)]
    junction_yellow: f64,
    /// Cross-contig clip links: maximum read-coordinate gap or overlap between
    /// the two alignments of one read that join across contigs.
    #[arg(long, default_value_t = 1000)]
    xlink_max_gap: u32,
    /// Debug: scan chunk size in bp.
    #[arg(long, hide = true)]
    chunk_size: Option<u64>,
    /// Write machine-readable progress to stderr (one JSON object per line).
    #[arg(long)]
    progress: bool,
}

/// `fish --scope`.
#[derive(Clone, Copy, clap::ValueEnum)]
enum ScopeArg {
    Findings,
    Genome,
}

/// `fish --min-tier`.
#[derive(Clone, Copy, clap::ValueEnum)]
enum TierArg {
    Yellow,
    Orange,
    Red,
}

/// `fish --bait-records`.
#[derive(Clone, Copy, clap::ValueEnum)]
enum BaitRecordsArg {
    Any,
    Primary,
}

#[derive(Args)]
struct FishArgs {
    /// Input BAM/CRAM (coordinate-sorted, indexed) mapped to the assembly.
    #[arg(long)]
    r#in: PathBuf,
    /// Assembly FASTA (the alignment reference).
    #[arg(long)]
    r#ref: PathBuf,
    /// Output prefix; writes `<prefix>.{fastq.gz,fish.tsv,fish.json}`.
    #[arg(long)]
    out_prefix: PathBuf,
    /// Threads (default: all available cores).
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// Linear contigs (comma-separated); all others are circular.
    #[arg(long, value_delimiter = ',')]
    linear: Vec<String>,
    /// `P.findings.bed`; required by `--locus` and `--scope findings`, and the
    /// default scope becomes `findings` when it is given.
    #[arg(long)]
    findings: Option<PathBuf>,
    /// `P.haplotypes.tsv`; linked selection requires it (there is no positive
    /// `--linked` flag).
    #[arg(long)]
    haplotypes: Option<PathBuf>,
    /// Bait FASTA(s) (e.g. a phage or plasmid contig). Repeatable; a
    /// multi-record FASTA is fine. Adds the `bait` category: reads that a full
    /// remap to assembly + bait places on a bait contig. Writes
    /// `<prefix>.bait_ref.fa` and `<prefix>.bait.cram` (+ index).
    #[arg(long)]
    bait: Vec<PathBuf>,
    /// Tiny-input fallback threshold in MB (10^6 bytes) for the bait `lift`
    /// session. Default: the per-preset value `cramhook lift` uses (lr:hq 20,
    /// map-hifi 3); 0 disables the fallback so reads are translated too.
    #[arg(long, value_parser = parse_min_input_mb)]
    bait_min_input_mb: Option<f64>,
    /// Which records select a read for `bait`: any counted record on a bait
    /// contig, primary or supplementary (default), or only a read whose primary
    /// record is on a bait contig. `primary` drops reads that reach the bait
    /// only through a supplementary alignment (its primary lies elsewhere).
    #[arg(long, value_enum, default_value_t = BaitRecordsArg::Any)]
    bait_records: BaitRecordsArg,
    /// Scope for clip/large-indel selection: findings (default with --findings) or genome.
    #[arg(long, value_enum)]
    scope: Option<ScopeArg>,
    /// Minimum finding tier for the findings scope, --locus and linked selection.
    #[arg(long, value_enum, default_value_t = TierArg::Yellow)]
    min_tier: TierArg,
    /// bp around each finding for the findings scope / finding column.
    #[arg(long, default_value_t = 1000)]
    finding_pad: u64,
    /// Do not select unmapped reads.
    #[arg(long)]
    no_unmapped: bool,
    /// Do not select terminal clips.
    #[arg(long)]
    no_clips: bool,
    /// Do not select large indels.
    #[arg(long)]
    no_indels: bool,
    /// Do not select haplotype-group reads.
    #[arg(long)]
    no_linked: bool,
    /// Also select any record overlapping a finding ± --locus-flank.
    #[arg(long)]
    locus: bool,
    /// bp added around a finding for `--locus`.
    #[arg(long, default_value_t = 5000)]
    locus_flank: u64,
    /// Terminal clip length threshold.
    #[arg(long, default_value_t = 50)]
    clip_min_len: u32,
    /// Minimum merged large-indel size.
    #[arg(long, default_value_t = 50)]
    indel_min_len: u32,
    /// Merge same-type indels of one read separated by at most this many bp.
    #[arg(long, default_value_t = 50)]
    indel_merge_gap: u32,
    /// Only indels of at least this size take part in the merge.
    #[arg(long, default_value_t = 10)]
    indel_piece_min: u32,
    /// Write machine-readable progress to stderr (one JSON object per line).
    #[arg(long)]
    progress: bool,
}

#[derive(Args)]
struct LiftArgs {
    /// Input CRAM mapped to the old assembly.
    #[arg(long)]
    r#in: PathBuf,
    /// Old assembly FASTA (the CRAM reference).
    #[arg(long)]
    old_ref: PathBuf,
    /// New assembly FASTA.
    #[arg(long)]
    new_ref: PathBuf,
    /// Output CRAM (coordinate-sorted, indexed).
    #[arg(long)]
    out: PathBuf,
    /// Threads (default: all available cores).
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// Write stats JSON to this path.
    #[arg(long)]
    stats: Option<PathBuf>,
    /// Tiny-input fallback threshold in MB (10^6 bytes): if the input CRAM is
    /// smaller than this, remap every read instead of lifting. Overrides the
    /// per-preset default (lr:hq 20, map-hifi 3); 0 disables the fallback.
    #[arg(long, value_parser = parse_min_input_mb)]
    min_input_mb: Option<f64>,
    /// Exit with an error on a contract violation instead of deferring all reads.
    #[arg(long)]
    strict: bool,
    /// Keep temporary files (deferred FASTQ copy and diff outputs) for debugging.
    #[arg(long)]
    keep_tmp: bool,
    /// Override the read preset when the input contract is violated and it
    /// cannot be recovered from the minimap2 @PG line.
    #[arg(long)]
    preset: Option<Preset>,
    /// Add an `XL:A:T`/`XL:A:R` tag to every output record (default; kept for
    /// compatibility).
    #[arg(long)]
    xl_tag: bool,
    /// Disable the `XL:A:T`/`XL:A:R` tag (default: on).
    #[arg(long)]
    no_xl_tag: bool,
    /// Minimum MAPQ for translation.
    #[arg(long, default_value_t = 60)]
    min_mapq: u8,
    /// Maximum soft-clip (total for single, outer for split) in bp; 0 = unlimited.
    #[arg(long, default_value_t = 200)]
    max_clip: u64,
    /// Window margin in bp.
    #[arg(long, default_value_t = 500)]
    margin: u64,
    /// Divergence (`de`) threshold; overrides the preset default (0.01 hifi, 0.03 ont).
    #[arg(long, value_parser = parse_max_de)]
    max_de: Option<f64>,
    /// T2 projection window length in bp.
    #[arg(long, default_value_t = 10_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_win: u64,
    /// T2 projection window step in bp.
    #[arg(long, default_value_t = 2_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_step: u64,
    /// T2 projection flank added around novel intervals in bp.
    #[arg(long, default_value_t = 10_000)]
    t2_flank: u64,
    /// T2 junction taint: half-width in bp of the interval added at every new
    /// adjacency (deletion/join); 0 disables.
    #[arg(long, default_value_t = 50)]
    t2_junction_flank: u64,
    /// Clip screen: k-mer length, 1..=32.
    #[arg(long, default_value_t = 17, value_parser = clap::value_parser!(u32).range(1..=32))]
    clip_k: u32,
    /// Clip screen: minimum k-mers shared with novel regions to defer a clip.
    #[arg(long, default_value_t = 2)]
    clip_min_hits: u32,
    /// Clip screen: minimum soft-clip length in bp to screen.
    #[arg(long, default_value_t = 20)]
    clip_min_len: u64,
    /// Clip screen: extra bp added around novel intervals / edit sites (0 = none).
    #[arg(long, default_value_t = 0)]
    clip_flank: u64,
    /// Disable the fail-safe clip screen (chimeric reads may lose supplementary records).
    #[arg(long)]
    unsafe_no_clip_screen: bool,
    /// Debug: clip screen Variant B' (drop k-mers that also occur in old.fa).
    #[arg(long, hide = true)]
    clip_subtract_old: bool,
    /// Defer reads whose competing chain score s2 is >= this fraction of s1 (MAPQ instability).
    #[arg(long, default_value_t = 0.5, value_parser = parse_max_s2_ratio)]
    max_s2_ratio: f64,
    /// Defer '-' block translations whose CIGAR has an indel of at least this many bp.
    #[arg(long, default_value_t = 50)]
    max_minus_indel: u64,
    /// Defer '-' block translations whose old-alignment path Z-drop statistic
    /// exceeds this (minimap2's -z; 0 is a legitimate strict value).
    #[arg(long, default_value_t = 400, value_parser = parse_max_zdrop)]
    max_zdrop: i64,
    /// Write machine-readable progress to stderr (one JSON object per line).
    #[arg(long)]
    progress: bool,
}

#[derive(Args)]
struct TranslateArgs {
    /// Input CRAM mapped to the old assembly.
    #[arg(long)]
    r#in: PathBuf,
    /// Old assembly FASTA (the CRAM reference).
    #[arg(long)]
    old_ref: PathBuf,
    /// New assembly FASTA.
    #[arg(long)]
    new_ref: PathBuf,
    /// Output prefix; writes `<prefix>.{translated.bam,deferred.fq,deferred.tsv,stats.json}`.
    #[arg(long)]
    out_prefix: PathBuf,
    /// Threads for minimap2 and CRAM decompression (default: all available cores).
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// Exit with an error on a contract violation instead of deferring all reads.
    #[arg(long)]
    strict: bool,
    /// Add an `XL:A:T` tag to every translated record (default; kept for
    /// compatibility).
    #[arg(long)]
    xl_tag: bool,
    /// Disable the `XL:A:T` tag (default: on).
    #[arg(long)]
    no_xl_tag: bool,
    /// Minimum MAPQ for translation.
    #[arg(long, default_value_t = 60)]
    min_mapq: u8,
    /// Maximum soft-clip (total for single, outer for split) in bp; 0 = unlimited.
    #[arg(long, default_value_t = 200)]
    max_clip: u64,
    /// Window margin in bp.
    #[arg(long, default_value_t = 500)]
    margin: u64,
    /// Divergence (`de`) threshold; overrides the preset default (0.01 hifi, 0.03 ont).
    #[arg(long, value_parser = parse_max_de)]
    max_de: Option<f64>,
    /// T2 projection window length in bp.
    #[arg(long, default_value_t = 10_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_win: u64,
    /// T2 projection window step in bp.
    #[arg(long, default_value_t = 2_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_step: u64,
    /// T2 projection flank added around novel intervals in bp.
    #[arg(long, default_value_t = 10_000)]
    t2_flank: u64,
    /// T2 junction taint: half-width in bp of the interval added at every new
    /// adjacency (deletion/join); 0 disables.
    #[arg(long, default_value_t = 50)]
    t2_junction_flank: u64,
    /// Clip screen: k-mer length, 1..=32.
    #[arg(long, default_value_t = 17, value_parser = clap::value_parser!(u32).range(1..=32))]
    clip_k: u32,
    /// Clip screen: minimum k-mers shared with novel regions to defer a clip.
    #[arg(long, default_value_t = 2)]
    clip_min_hits: u32,
    /// Clip screen: minimum soft-clip length in bp to screen.
    #[arg(long, default_value_t = 20)]
    clip_min_len: u64,
    /// Clip screen: extra bp added around novel intervals / edit sites (0 = none).
    #[arg(long, default_value_t = 0)]
    clip_flank: u64,
    /// Disable the fail-safe clip screen (chimeric reads may lose supplementary records).
    #[arg(long)]
    unsafe_no_clip_screen: bool,
    /// Debug: clip screen Variant B' (drop k-mers that also occur in old.fa).
    #[arg(long, hide = true)]
    clip_subtract_old: bool,
    /// Defer reads whose competing chain score s2 is >= this fraction of s1 (MAPQ instability).
    #[arg(long, default_value_t = 0.5, value_parser = parse_max_s2_ratio)]
    max_s2_ratio: f64,
    /// Defer '-' block translations whose CIGAR has an indel of at least this many bp.
    #[arg(long, default_value_t = 50)]
    max_minus_indel: u64,
    /// Defer '-' block translations whose old-alignment path Z-drop statistic
    /// exceeds this (minimap2's -z; 0 is a legitimate strict value).
    #[arg(long, default_value_t = 400, value_parser = parse_max_zdrop)]
    max_zdrop: i64,
}

#[derive(Args)]
struct DiffArgs {
    /// Old assembly FASTA.
    #[arg(long)]
    old_ref: PathBuf,
    /// New assembly FASTA.
    #[arg(long)]
    new_ref: PathBuf,
    /// Read preset used for the T2 probe.
    #[arg(long)]
    preset: Preset,
    /// Output prefix; writes `<prefix>.{diff.paf,t2.paf,blocks.bed,taint.bed,novel.bed,summary.json}`.
    #[arg(long)]
    out_prefix: PathBuf,
    /// Threads for minimap2 (default: all available cores).
    #[arg(long, default_value_t = default_threads())]
    threads: usize,
    /// T2 projection window length in bp.
    #[arg(long, default_value_t = 10_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_win: u64,
    /// T2 projection window step in bp.
    #[arg(long, default_value_t = 2_000, value_parser = clap::value_parser!(u64).range(1..))]
    t2_step: u64,
    /// T2 projection flank added around novel intervals in bp.
    #[arg(long, default_value_t = 10_000)]
    t2_flank: u64,
    /// T2 junction taint: half-width in bp of the interval added at every new
    /// adjacency (deletion/join); 0 disables.
    #[arg(long, default_value_t = 50)]
    t2_junction_flank: u64,
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Lift(a) => run_lift_cmd(a),
        Cmd::Diff(a) => run_diff_cmd(a),
        Cmd::Translate(a) => run_translate_cmd(a),
        Cmd::Qc(a) => run_qc_cmd(a),
        Cmd::Fish(a) => run_fish_cmd(a),
    };
    match r {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            // Same text and exit code as returning the `Result` from `main`,
            // but a closed stderr cannot panic.
            note!("Error: {e:?}");
            std::process::ExitCode::FAILURE
        }
    }
}

/// `--progress`: a stderr reporter for `cmd`, or a disabled one.
fn progress(on: bool, cmd: &'static str) -> cramhook::Progress {
    if on { cramhook::Progress::stderr(cmd) } else { cramhook::Progress::off() }
}

fn run_lift_cmd(a: LiftArgs) -> Result<()> {
    let tcfg = cramhook::TranslateConfig {
        strict: a.strict,
        xl_tag: !a.no_xl_tag,
        min_mapq: a.min_mapq,
        max_clip: a.max_clip,
        margin: a.margin,
        max_de: a.max_de,
        t2_win: a.t2_win,
        t2_step: a.t2_step,
        t2_flank: a.t2_flank,
        t2_junction_flank: a.t2_junction_flank,
        clip_k: a.clip_k,
        clip_min_hits: a.clip_min_hits,
        clip_min_len: a.clip_min_len,
        clip_flank: a.clip_flank,
        unsafe_no_clip_screen: a.unsafe_no_clip_screen,
        clip_subtract_old: a.clip_subtract_old,
        max_s2_ratio: a.max_s2_ratio,
        max_minus_indel: a.max_minus_indel,
        max_zdrop: a.max_zdrop,
        // Rule 7d scoring is set from the resolved read preset at runtime.
        scoring: Default::default(),
        // Chain state, set from the input contract by `lift`/`translate`.
        inherited_xl_ok: false,
    };
    let opts = cramhook::LiftOptions { strict: a.strict, keep_tmp: a.keep_tmp, min_input_mb: a.min_input_mb };
    let stats = cramhook::lift_with_progress(
        &a.r#in,
        &a.old_ref,
        &a.new_ref,
        &a.out,
        a.stats.as_deref(),
        a.threads,
        &tcfg,
        &opts,
        a.preset,
        &progress(a.progress, "lift"),
    )?;
    note!("contract ok: {}", stats.contract_ok);
    if let Some(r) = &stats.contract_reason {
        note!("contract reason: {r}");
    }
    note!("preset: {}", stats.preset);
    note!(
        "records in: primary {}, supplementary {}, secondary {}, unmapped {}",
        stats.records_in.primary,
        stats.records_in.supplementary,
        stats.records_in.secondary,
        stats.records_in.unmapped
    );
    note!(
        "translated: {} reads ({} records); remapped: {} reads; out primary: {}",
        stats.translated_reads,
        stats.translated_records,
        stats.remapped_reads,
        stats.out_primary_records
    );
    note!(
        "origin-crossing: {} translated, {} deferred",
        stats.origin_crossing_translated, stats.origin_crossing_deferred
    );
    for (reason, n) in &stats.deferred_by_reason {
        note!("  deferred {reason}: {n}");
    }
    let t = &stats.timings_s;
    note!(
        "timings: diff {:.2}s, pass {:.2}s, drain {:.2}s, sort+index {:.2}s, index {:.2}s, minimap2 {:.2}s, children_cpu {:.1}s, total {:.2}s",
        t.diff, t.translate_pass, t.writer_drain, t.sort_index, t.minimap2_index, t.minimap2_total, t.children_cpu, t.total
    );
    Ok(())
}

fn run_translate_cmd(a: TranslateArgs) -> Result<()> {
    let cfg = cramhook::TranslateConfig {
        strict: a.strict,
        xl_tag: !a.no_xl_tag,
        min_mapq: a.min_mapq,
        max_clip: a.max_clip,
        margin: a.margin,
        max_de: a.max_de,
        t2_win: a.t2_win,
        t2_step: a.t2_step,
        t2_flank: a.t2_flank,
        t2_junction_flank: a.t2_junction_flank,
        clip_k: a.clip_k,
        clip_min_hits: a.clip_min_hits,
        clip_min_len: a.clip_min_len,
        clip_flank: a.clip_flank,
        unsafe_no_clip_screen: a.unsafe_no_clip_screen,
        clip_subtract_old: a.clip_subtract_old,
        max_s2_ratio: a.max_s2_ratio,
        max_minus_indel: a.max_minus_indel,
        max_zdrop: a.max_zdrop,
        // Rule 7d scoring is set from the resolved read preset at runtime.
        scoring: Default::default(),
        // Chain state, set from the input contract by `lift`/`translate`.
        inherited_xl_ok: false,
    };
    let stats =
        cramhook::translate(&a.r#in, &a.old_ref, &a.new_ref, &a.out_prefix, &cfg, a.threads)?;
    note!("contract ok: {}", stats.contract_ok);
    if let Some(r) = &stats.contract_reason {
        note!("contract reason: {r}");
    }
    note!("preset: {}", stats.preset);
    note!(
        "records in: primary {}, supplementary {}, secondary {}, unmapped {}",
        stats.records_in.primary,
        stats.records_in.supplementary,
        stats.records_in.secondary,
        stats.records_in.unmapped
    );
    note!(
        "translated: {} reads ({} records); deferred: {} reads",
        stats.translated_reads, stats.translated_records, stats.deferred_reads
    );
    note!(
        "origin-crossing: {} translated, {} deferred",
        stats.origin_crossing_translated, stats.origin_crossing_deferred
    );
    if stats.qual_missing > 0 {
        note!("deferred reads with missing QUAL: {}", stats.qual_missing);
    }
    for (reason, n) in &stats.deferred_by_reason {
        note!("  deferred {reason}: {n}");
    }
    Ok(())
}

fn run_qc_cmd(a: QcArgs) -> Result<()> {
    let params = cramhook::QcParams {
        preset: a.preset,
        threads: a.threads,
        linear: a.linear,
        candidates: a.candidates,
        min_depth: a.min_depth,
        max_mq0_frac: a.max_mq0_frac,
        min_alt: a.min_alt,
        snp_min_af: a.snp_min_af,
        indel_min_af: a.indel_min_af,
        lowcov_abs: a.lowcov_abs,
        lowcov_rel: a.lowcov_rel,
        lowcov_merge: a.lowcov_merge,
        lowcov_min_len: a.lowcov_min_len,
        fdr: a.fdr,
        strand_ratio: a.strand_ratio,
        red_min_alt: a.red_min_alt,
        yellow_strand_gate: !a.no_yellow_strand_gate,
        gray_homopolymer: !a.no_gray_homopolymer,
        gray_hp_min: a.gray_hp_min,
        large_merge_gap: a.large_merge_gap,
        large_min: a.large_min,
        large_piece_min: a.large_piece_min,
        large_min_support: a.large_min_support,
        large_min_ratio: a.large_min_ratio,
        cluster_window: a.cluster_window,
        cluster_min_sites: a.cluster_min_sites,
        cluster_site_p: a.cluster_site_p,
        cluster_density_p: a.cluster_density_p,
        cluster_strand: a.cluster_strand,
        ins_group_min_len: a.ins_group_min_len,
        ins_group_len_frac: a.ins_group_len_frac,
        ins_group_max_diff: a.ins_group_max_diff,
        clip_min_len: a.clip_min_len,
        clip_merge: a.clip_merge,
        clip_min_reads: a.clip_min_reads,
        clip_yellow: a.clip_yellow,
        clip_high_s: a.clip_high_s,
        clip_pair_dist: a.clip_pair_dist,
        clip_edge: a.clip_edge,
        gc_window: a.gc_window,
        depth_window: a.depth_window,
        depth_low: a.depth_low,
        depth_half_min: a.depth_half_min,
        depth_half_max: a.depth_half_max,
        depth_half_max_mapq: a.depth_half_max_mapq,
        depth_high: a.depth_high,
        depth_high_min_mapq: a.depth_high_min_mapq,
        depth_gap_windows: a.depth_gap_windows,
        depth_min_len: a.depth_min_len,
        depth_baseline_span: a.depth_baseline_span,
        step_window: a.step_window,
        step_span: a.step_span,
        step_nms: a.step_nms,
        step_yellow: a.step_yellow,
        step_red: a.step_red,
        finding_merge: a.finding_merge,
        finding_site_join: a.finding_site_join,
        finding_link_max: a.finding_link_max,
        link_min_af: a.link_min_af,
        link_window: a.link_window,
        link_min_reads: a.link_min_reads,
        link_min_conc: a.link_min_conc,
        link_p: a.link_p,
        link_strand_p: a.link_strand_p,
        link_conflict_min_reads: a.link_conflict_min_reads,
        link_conflict_max_share: a.link_conflict_max_share,
        unlinked_action: a.unlinked_action,
        circ_end_tol: a.circ_end_tol,
        circ_min_reads: a.circ_min_reads,
        circ_min_frac: a.circ_min_frac,
        junction_min_dist: a.junction_min_dist,
        junction_max_gap: a.junction_max_gap,
        junction_tol: a.junction_tol,
        junction_edge_excl: a.junction_edge_excl,
        junction_min_reads: a.junction_min_reads,
        junction_yellow: a.junction_yellow,
        xlink_max_gap: a.xlink_max_gap,
        chunk_size: a.chunk_size,
        progress: progress(a.progress, "qc"),
        ..Default::default()
    };
    let s = cramhook::run_qc(&a.r#in, &a.r#ref, &a.out_prefix, &params)?;
    note!("preset: {}", s.preset);
    note!(
        "records: {} primary mapped, {} supplementary, {} unmapped ({:.2}% of reads)",
        s.reads.records.primary, s.reads.records.supplementary, s.reads.unmapped,
        100.0 * s.reads.unmapped_fraction
    );
    note!(
        "assembly: {} bp assessable ({:.4}), {} bp confident, {} bp low coverage",
        s.assembly.assessable_bp, s.assembly.assessable_fraction, s.assembly.confident_bp, s.assembly.low_coverage_bp
    );
    note!(
        "alerts: {} red, {} orange, {} yellow, {} gray (PASS); lowcov: {} red, {} orange, {} yellow, {} gray",
        s.alerts.counts.red, s.alerts.counts.orange, s.alerts.counts.yellow, s.alerts.counts.gray,
        s.alerts.counts.lowcov.red, s.alerts.counts.lowcov.orange, s.alerts.counts.lowcov.yellow,
        s.alerts.counts.lowcov.gray
    );
    for w in &s.warnings {
        note!("warning {}: {} bp ({:.4}) on {:?}", w.code, w.bp, w.fraction, w.contigs);
    }
    note!("qc_qv: {:.2}", s.qc_qv);
    note!("timings: scan {:.2}s, qc {:.2}s, total {:.2}s", s.timing_s.scan, s.timing_s.qc, s.timing_s.total);
    Ok(())
}

fn run_fish_cmd(a: FishArgs) -> Result<()> {
    let params = cramhook::FishParams {
        threads: a.threads,
        linear: a.linear,
        findings: a.findings,
        haplotypes: a.haplotypes,
        bait: a.bait,
        bait_min_input_mb: a.bait_min_input_mb,
        bait_records: match a.bait_records {
            BaitRecordsArg::Any => cramhook::BaitRecords::Any,
            BaitRecordsArg::Primary => cramhook::BaitRecords::Primary,
        },
        scope: a.scope.map(|s| match s {
            ScopeArg::Findings => cramhook::FishScope::Findings,
            ScopeArg::Genome => cramhook::FishScope::Genome,
        }),
        min_tier: match a.min_tier {
            TierArg::Yellow => cramhook::FishTier::Yellow,
            TierArg::Orange => cramhook::FishTier::Orange,
            TierArg::Red => cramhook::FishTier::Red,
        },
        finding_pad: a.finding_pad,
        no_unmapped: a.no_unmapped,
        no_clips: a.no_clips,
        no_indels: a.no_indels,
        no_linked: a.no_linked,
        locus: a.locus,
        locus_flank: a.locus_flank,
        clip_min_len: a.clip_min_len,
        indel_min_len: a.indel_min_len,
        indel_merge_gap: a.indel_merge_gap,
        indel_piece_min: a.indel_piece_min,
        progress: progress(a.progress, "fish"),
        ..Default::default()
    };
    let s = cramhook::run_fish(&a.r#in, &a.r#ref, &a.out_prefix, &params)?;
    note!("scope: {}", s.params.scope);
    note!(
        "selected: {} reads, {} bp",
        s.reads.reads, s.reads.bp
    );
    for (reason, c) in &s.reads.by_reason {
        note!("  {reason}: {} reads, {} bp", c.reads, c.bp);
    }
    for n in &s.notes {
        note!("note: {n}");
    }
    if let Some(b) = &s.bait {
        let fb = b.fallback.as_deref().map(|f| format!(", fallback {f}")).unwrap_or_default();
        note!(
            "bait: {} reads, {} bp ({} MAPQ 0); lift: {} translated, {} remapped{}",
            b.reads, b.bp, b.mapq0_reads, b.translated_reads, b.remapped_reads, fb
        );
    }
    note!(
        "records: {} primary, {} supplementary, {} secondary, {} duplicate, {} qcfail, {} unmapped",
        s.records.primary,
        s.records.supplementary,
        s.records.secondary,
        s.records.duplicate,
        s.records.qcfail,
        s.records.unmapped
    );
    note!("timings: pass {:.2}s, total {:.2}s", s.timing_s.pass, s.timing_s.total);
    Ok(())
}

fn run_diff_cmd(a: DiffArgs) -> Result<()> {
    let params = cramhook::DiffParams {
        win: a.t2_win,
        step: a.t2_step,
        flank: a.t2_flank,
        junction_flank: a.t2_junction_flank,
    };
    let run = cramhook::run_diff_with(
        &a.old_ref,
        &a.new_ref,
        a.preset,
        a.threads,
        &a.out_prefix,
        params,
    )?;
    let s = &run.summary;
    note!("primary alignments: {}", s.n_primary_alignments);
    note!("blocks:  {} ({} bp)", s.blocks.count, s.blocks.bp);
    note!("taint T1: {} ({} bp)", s.taint.t1.count, s.taint.t1.bp);
    note!("taint T2: {} ({} bp)", s.taint.t2.count, s.taint.t2.bp);
    note!("novel:   {} ({} bp)", s.novel.count, s.novel.bp);
    note!("windows: {}", s.windows);
    note!(
        "timings: diff {:.2}s, t1 {:.2}s, novel {:.2}s, t2 {:.2}s, total {:.2}s",
        s.timings_s.diff, s.timings_s.t1, s.timings_s.novel, s.timings_s.t2, s.timings_s.total
    );
    Ok(())
}
