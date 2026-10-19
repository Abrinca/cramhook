//! `cramhook qc`: per-site signals from the `scan` counters.
//!
//! Counters → candidates (pre-filter) → beta-binomial background per class
//! (fitted genome-wide by moments) → p-value, Benjamini–Hochberg q-value →
//! strand balance → yellow / orange / red tiers. Plus minority-SNP clusters,
//! large indels and not-assessable regions. Outputs: VCF (sites), BED
//! (regions), JSON (summary), bedGraph (depth) and optionally a candidates
//! TSV. All thresholds are starting values (see `QcParams`).

pub(crate) mod alerts;
pub(crate) mod calling;
pub(crate) mod circularity;
pub(crate) mod findings;
pub(crate) mod junctions;
pub(crate) mod likely;
pub(crate) mod link;
pub(crate) mod output;
pub(crate) mod params;
pub(crate) mod regions;
pub(crate) mod sig;
pub(crate) mod sites;
pub(crate) mod stats;
pub(crate) mod tiers;

#[cfg(test)]
pub(crate) mod test_util;

pub(crate) use alerts::*;
pub(crate) use calling::*;
pub(crate) use findings::*;
pub(crate) use likely::*;
pub use params::*;
pub(crate) use regions::*;
pub(crate) use sig::*;
pub(crate) use sites::*;
pub(crate) use tiers::*;

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Result, bail};

use crate::fasta::Fasta;
use crate::scan::{self, ScanParams};
use stats::{BbFit, MomentSums};

/// `--progress` phases of `qc`, in order.
pub const QC_PHASES: [&str; 5] = ["scan", "circularity", "call", "linkage", "output"];

/// Run `cramhook qc`. Returns the JSON summary that was written.
pub fn run_qc(input: &Path, ref_path: &Path, out_prefix: &Path, params: &QcParams) -> Result<output::Summary> {
    params.progress.plan(&QC_PHASES);
    params.progress.finish(run_qc_inner(input, ref_path, out_prefix, params))
}

fn run_qc_inner(input: &Path, ref_path: &Path, out_prefix: &Path, params: &QcParams) -> Result<output::Summary> {
    let t0 = Instant::now();
    let progress = &params.progress;
    validate_lowcov_rel(params.lowcov_rel)?;
    validate_ins_group_frac(params.ins_group_len_frac)?;
    validate_ins_group_frac(params.ins_group_max_diff)?;
    let preset = match params.preset {
        Some(p) => p,
        None => {
            let h = scan::read_header_text(input, ref_path)?;
            match crate::lift::best_effort_preset(&h) {
                Some(p) => p,
                None => bail!(
                    "cannot determine the read preset from the minimap2 @PG line of {}; pass --preset",
                    input.display()
                ),
            }
        }
    };
    let indel_min_af = params.indel_min_af.unwrap_or_else(|| default_indel_min_af(preset));
    let lowcov_abs = params.lowcov_abs.unwrap_or_else(|| default_lowcov_abs(preset));
    let fasta = Fasta::load(ref_path)?;
    for l in &params.linear {
        if fasta.seq(l).is_none() {
            bail!("--linear: contig {l:?} is not in {}", ref_path.display());
        }
    }
    let sp = ScanParams {
        threads: params.threads.max(1),
        large_merge_gap: params.large_merge_gap,
        large_min: params.large_min,
        large_piece_min: params.large_piece_min,
        clip_min_len: params.clip_min_len,
        chunk_size: params.chunk_size,
        progress: progress.clone(),
    };
    progress.start("scan");
    let scan = scan::scan(input, ref_path, &fasta, &sp)?;
    progress.end();
    let scan_s = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();

    // Circularity pass (its own scan of the contig ends; also derives the
    // split-read junction events). It needs only `scan`/`fasta`, so it runs
    // before the per-contig contexts: the wrap decision (verdict/header/flag)
    // then drives every later wrap (depth bands, depth step, clusters and the
    // junction edge exclusion), so no input to this pass depends on the wrap
    // flag itself.
    progress.start("circularity");
    let circ = circularity::run_circularity(input, ref_path, &scan, &fasta, params)?;
    progress.end();
    progress.start("call");
    let circ_s = circ.timing_s;

    let linear: HashSet<&str> = params.linear.iter().map(|s| s.as_str()).collect();
    let refs: Vec<&[u8]> = scan.contigs.iter().map(|c| fasta.seq(&c.name).unwrap()).collect();
    let ctxs: Vec<ContigCtx> = scan
        .contigs
        .iter()
        .zip(&refs)
        .enumerate()
        .map(|(tid, (c, r))| {
            let (wrap, by) = circularity::wrap_decision(&circ.per_contig[tid], linear.contains(c.name.as_str()));
            let mut ctx = contig_ctx(c, r, wrap, params, lowcov_abs, params.lowcov_rel);
            ctx.circular_by = by;
            ctx
        })
        .collect();

    // Indel slices per contig (indels are sorted by tid).
    let mut ranges = vec![(0usize, 0usize); scan.contigs.len()];
    {
        let mut i = 0;
        while i < scan.indels.len() {
            let t = scan.indels[i].tid as usize;
            let s = i;
            while i < scan.indels.len() && scan.indels[i].tid as usize == t {
                i += 1;
            }
            ranges[t] = (s, i);
        }
    }

    // Grouped INS alleles per contig (also used by the background fit).
    let groups: Vec<Vec<InsGroup>> = scan
        .contigs
        .iter()
        .enumerate()
        .map(|(tid, _)| {
            let (a, b) = ranges[tid];
            ins_groups(&scan.indels[a..b], params)
        })
        .collect();

    // Background fit.
    let mut sums = [MomentSums::default(); NCLASS];
    let mut m = [0u64; NCLASS];
    for (tid, c) in scan.contigs.iter().enumerate() {
        let (a, b) = ranges[tid];
        accumulate_background(c, refs[tid], &ctxs[tid], &scan.indels[a..b], &groups[tid], params, &mut sums, &mut m);
    }
    let mut pooled = [MomentSums::default(); 2];
    for b in 0..NB {
        pooled[0].add(&sums[1 + b]);
        pooled[1].add(&sums[1 + NB + b]);
    }
    let pooled_fit = [stats::fit_bb(&pooled[0]), stats::fit_bb(&pooled[1])];
    let mut fits: Vec<ClassFit> = (0..NCLASS)
        .map(|c| {
            let own = stats::fit_bb(&sums[c]);
            let fallback = c > 0 && sums[c].sites < params.min_fit_sites;
            let fit = if fallback { pooled_fit[if c <= NB { 0 } else { 1 }] } else { own };
            ClassFit { class: class_name(c), fit, fit_sites: sums[c].sites, m: m[c], fallback, candidates: 0 }
        })
        .collect();

    // Candidates.
    let mut cands: Vec<Candidate> = Vec::new();
    let mut cluster_sites: Vec<Vec<u32>> = vec![Vec::new(); scan.contigs.len()];
    for (tid, c) in scan.contigs.iter().enumerate() {
        let (a, b) = ranges[tid];
        collect_candidates(
            tid,
            c,
            refs[tid],
            &ctxs[tid],
            &scan.indels[a..b],
            &groups[tid],
            params,
            indel_min_af,
            lowcov_abs,
            &fits[0].fit,
            &mut cands,
            &mut cluster_sites[tid],
        );
    }

    // p-values, BH per class, tiers.
    let mut by_class: Vec<Vec<usize>> = vec![Vec::new(); NCLASS];
    for (i, cd) in cands.iter_mut().enumerate() {
        let f = &fits[cd.class].fit;
        cd.ln_p = stats::bb_ln_sf(cd.alt(), cd.dp, f.alpha, f.beta);
        by_class[cd.class].push(i);
    }
    for (cl, idx) in by_class.iter().enumerate() {
        fits[cl].candidates = idx.len() as u64;
        let pv: Vec<f64> = idx.iter().map(|&i| cands[i].ln_p.exp()).collect();
        let q = stats::bh(&pv, m[cl] as usize);
        for (k, &i) in idx.iter().enumerate() {
            cands[i].q = q[k];
        }
    }
    for cd in cands.iter_mut() {
        cd.significant = cd.q <= params.fdr;
        cd.tier = tier_of(
            &TierIn { af: cd.af(), alt: cd.alt(), wilson: cd.wilson, strand_ok: cd.strand_ok, significant: cd.significant },
            params,
        );
        // Gray: a minority/strand-biased yellow INS/DEL (`minority` or
        // `strand_bias`) in a homopolymer of at least `gray_hp_min` (INFO `HP`).
        // A `weak_majority` yellow (AF > 0.5 that missed the red bar) stays
        // yellow: a majority HP indel may mean the consensus itself is wrong.
        // Red/orange calls and SNPs are never gray. The reason is kept;
        // `HPGRAY` marks the call. Applied before linkage/findings so gray is
        // not a finding member and the linkage pass (SNPs only) is untouched.
        if params.gray_homopolymer
            && cd.tier.tier == Tier::Yellow
            && matches!(cd.tier.reason, "minority" | "strand_bias")
            && matches!(cd.vtype, VarType::Ins | VarType::Del)
            && cd.hp >= params.gray_hp_min
        {
            cd.tier.tier = Tier::Gray;
        }
    }
    cands.sort_by(|x, y| {
        (x.tid, x.pos0, x.origin, x.vtype, x.len, &x.alt_allele)
            .cmp(&(y.tid, y.pos0, y.origin, y.vtype, y.len, &y.alt_allele))
    });

    // Regions (the circularity pass above supplied the junction events and the
    // per-contig wrap decision).
    let mut regions: Vec<Region> = Vec::new();
    // Sorted junction-event p1 positions per contig: the junction score counts
    // background events in a window by binary search instead of rescanning.
    let junc_starts: Vec<Vec<u32>> = circ
        .junctions
        .iter()
        .map(|evs| {
            let mut v: Vec<u32> = evs.iter().map(|e| e.p1 as u32).collect();
            v.sort_unstable();
            v
        })
        .collect();
    for (tid, ctx) in ctxs.iter().enumerate() {
        not_assessable_regions(tid, &scan.contigs[tid], ctx, &mut regions);
        low_coverage_regions(tid, &scan.contigs[tid], ctx, params, &mut regions);
        snp_cluster_regions(tid, &scan.contigs[tid], ctx, &cluster_sites[tid], params, &mut regions);
        clip_hotspot_regions(tid, &scan.contigs[tid], ctx, params, lowcov_abs, &mut regions);
        junction_regions(
            tid,
            &scan.contigs[tid],
            ctx,
            params,
            &circ.junctions[tid],
            &junc_starts[tid],
            &scan.contigs[tid].name,
            &mut regions,
        );
    }
    large_indel_regions(&scan, &ctxs, params, lowcov_abs, &mut regions);
    let mut depth = Vec::new();
    depth_regions(&scan, &ctxs, params, scan.mean_aligned_span, &mut depth);
    let depth_step_origin_dropped = depth_step_regions(&ctxs, params, lowcov_abs, &mut depth);
    regions.extend(depth);
    regions.sort_by(|a, b| (a.tid, a.start, a.end, a.label).cmp(&(b.tid, b.start, b.end, b.label)));

    // Linkage second pass (step 2d): mark LINKED/UNLINKED sites and emit
    // `linked_haplotype` regions. Runs after the depth bands so a trusted
    // `depth_high` region can make a group red.
    let names: Vec<&str> = scan.contigs.iter().map(|c| c.name.as_str()).collect();
    progress.end();
    progress.start("linkage");
    let link = link::run_linkage(
        input, ref_path, &scan, &ctxs, &refs, &mut cands, &regions, &names, params, lowcov_abs,
    )?;
    progress.end();
    progress.start("output");
    regions.extend(link.regions);
    regions.sort_by(|a, b| (a.tid, a.start, a.end, a.label).cmp(&(b.tid, b.start, b.end, b.label)));

    // Findings: group the tiered calls into one problem location each, then
    // pool the groups the event links connect (`--finding-link-max`).
    let lens: Vec<u64> = refs.iter().map(|r| r.len() as u64).collect();
    let circular: Vec<bool> = ctxs.iter().map(|c| c.circular).collect();
    let findings = findings_linked(
        scan.contigs.len(),
        &cands,
        &regions,
        params.finding_merge,
        params.finding_site_join,
        params.finding_link_max,
        &lens,
        &circular,
    );

    // Cross-contig clip links between `clip_hotspot` regions and cluster
    // partners elsewhere, resolved against the findings just built. Links only
    // annotate regions/findings; tiers and `likely` are unchanged.
    let xlinks = collect_xlinks(&circ.xlinks, &regions, &ctxs, &lens, &findings, params.clip_merge);

    let qc_s = (t1.elapsed().as_secs_f64() - circ_s).max(0.0);
    let outs = output::Outputs {
        input,
        ref_path,
        prefix: out_prefix,
        preset,
        params,
        indel_min_af,
        lowcov_abs,
        fasta_order: fasta.order(),
        scan: &scan,
        refs: &refs,
        ctxs: &ctxs,
        names: &names,
        cands: &cands,
        regions: &regions,
        findings: &findings,
        fits: &fits,
        linkage: &link.stats,
        haplotypes: &link.haplotypes,
        circ: &circ,
        xlinks: &xlinks,
        depth_step_origin_dropped: &depth_step_origin_dropped,
    };
    let summary = outs.write_all(scan_s, qc_s, circ_s, t0)?;
    progress.end();
    Ok(summary)
}

/// Helper for callers: output path for a suffix.
pub(crate) fn out_path(prefix: &Path, suffix: &str) -> PathBuf {
    crate::translate::suffix(prefix, suffix)
}
