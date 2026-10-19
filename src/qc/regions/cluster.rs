//! `qc` not-assessable and minority-SNP cluster regions.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use crate::qc::*;
use crate::scan::*;

pub(crate) fn not_assessable_regions(tid: usize, c: &ContigCounts, ctx: &ContigCtx, out: &mut Vec<Region>) {
    let mut i = 0;
    while i < ctx.na.len() {
        let r = ctx.na[i];
        if r == 0 {
            i += 1;
            continue;
        }
        let s = i;
        while i < ctx.na.len() && ctx.na[i] == r {
            i += 1;
        }
        let (dr, hq) = region_depth(c, ctx, s as u64, i as u64);
        out.push(Region {
            tid,
            start: s as u64,
            end: i as u64,
            label: "not_assessable",
            score: 0,
            tier: "none",
            reason: NA_REASONS[r as usize].to_string(),
            depth_ratio: dr,
            hiq_share: hq,
            extra: String::new(),
            lowcov: false,
            link: None,
            clip_left: None,
            clip_right: None,
        });
    }
}

/// SNP-cluster regions from the qualifying SNP sites (`cluster_sites`, SNP
/// only): gray calls are INS/DEL, so they are never cluster sites and need no
/// exclusion here.
pub(crate) fn snp_cluster_regions(
    tid: usize,
    c: &ContigCounts,
    ctx: &ContigCtx,
    sites: &[u32],
    p: &QcParams,
    out: &mut Vec<Region>,
) {
    let len = c.len as u64;
    let assessable = ctx.na.iter().filter(|&&x| x == 0).count() as u64;
    let kmin = match p.cluster_density_p {
        Some(alpha) => {
            density_min_sites(sites.len() as u64, assessable, len, ctx.circular, p.cluster_window as u64, alpha)
                .max(p.cluster_min_sites)
        }
        None => p.cluster_min_sites,
    };
    for (s, e, n) in cluster_windows(sites, len, ctx.circular, p.cluster_window as u64, kmin) {
        let (dr, hq) = region_depth(c, ctx, s, e);
        // A cluster window is "highly significant" (its gate); tier by the
        // depth fraction f = (r - 1) / r, so red only when f > 0.5 (r > 2).
        let tier = match dr {
            Some(r) if (r - 1.0) / r > 0.5 => "red",
            _ => "orange",
        };
        let wraps = e > len;
        // `kmin` only when the density test (B) is on, so the all-off setting
        // reproduces the original BED.
        let kmin_s = if p.cluster_density_p.is_some() { format!(";kmin={kmin}") } else { String::new() };
        let extra =
            format!("sites={n};window={}{kmin_s}{}", p.cluster_window, if wraps { ";wraps_origin=1" } else { "" });
        let pieces: Vec<(u64, u64)> = if wraps { vec![(s, len), (0, e - len)] } else { vec![(s, e)] };
        for (a, b) in pieces {
            out.push(Region {
                tid,
                start: a,
                end: b,
                label: "snp_cluster",
                score: n.min(1000) as u32,
                tier,
                reason: "minority_snp_cluster".to_string(),
                depth_ratio: dr,
                hiq_share: hq,
                extra: extra.clone(),
                lowcov: false,
                link: None,
                clip_left: None,
                clip_right: None,
            });
        }
    }
}

