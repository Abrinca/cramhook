//! `qc` finding grouping and linking.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

use super::*;

/// One finding: the tiered calls of one problem location grouped together.
#[derive(Clone, Debug)]
pub(crate) struct Finding {
    pub tid: usize,
    pub start: u64,
    pub end: u64,
    /// Highest member tier (`red`/`orange`/`yellow`).
    pub tier: &'static str,
    /// All members are LowCov / `lowcov=1`.
    pub lowcov: bool,
    /// Number of members.
    pub score: u32,
    pub likely: &'static str,
    pub evidence: Evidence,
    /// Member label counts in [`ALERT_TYPES`] order (zero counts dropped).
    pub members: Vec<(&'static str, u32)>,
}

impl Finding {
    /// Comma list of the evidence flags present, in a fixed order (the
    /// `Evidence` field order): sites, then long indels, clips, depth, cluster,
    /// linked.
    pub(crate) fn evidence_list(&self) -> String {
        let e = &self.evidence;
        let mut v: Vec<&str> = Vec::new();
        if e.site_major {
            v.push("site_major");
        }
        if e.site_mixed {
            v.push("site_mixed");
        }
        if e.site_minor {
            v.push("site_minor");
        }
        if e.long_ins {
            v.push("long_ins");
        }
        if e.long_del {
            v.push("long_del");
        }
        if e.clips {
            v.push("clips");
        }
        if e.clips_both {
            v.push("clips_both");
        }
        if e.depth_low {
            v.push("depth_low");
        }
        if e.depth_half {
            v.push("depth_half");
        }
        if e.depth_high {
            v.push("depth_high");
        }
        if e.depth_step {
            v.push("depth_step");
        }
        if e.cluster {
            v.push("cluster");
        }
        if e.linked {
            v.push("linked");
        }
        v.join(",")
    }
}

/// One tiered call before grouping: its interval plus the data the evidence
/// flags and member counts are derived from.
#[derive(Clone, Debug)]
enum Member {
    Site { vtype: VarType, af: f64, len: u32, lowcov: bool },
    Region { label: &'static str, reason: String, lowcov: bool, depth_high: bool },
}

#[derive(Clone, Debug)]
struct MemberIv {
    tid: usize,
    start: u64,
    end: u64,
    tier: Tier,
    m: Member,
}

impl MemberIv {
    fn lowcov(&self) -> bool {
        match &self.m {
            Member::Site { lowcov, .. } | Member::Region { lowcov, .. } => *lowcov,
        }
    }}

/// Group the tiered calls into findings on one contig. Region-type members
/// (tiered BED regions) join each other when their intervals are within `merge`
/// bp, and site calls (VCF records) join each other the same way. A site joins
/// a finding that has a region-type member only when it lies within `site_join`
/// bp of one of those regions; a site farther than that stays out even if it is
/// within `merge` of another site, so an unrelated background region no longer
/// absorbs a nearby small error. A site call's interval is its REF span.
/// Untiered annotations are not members. Output is sorted by contig order then
/// start. This is the base grouping with linking disabled (see
/// [`findings_linked`]); it is only used by tests.
#[cfg(test)]
pub(crate) fn collect_findings(
    ncontigs: usize,
    cands: &[Candidate],
    regions: &[Region],
    merge: u32,
    site_join: u32,
) -> Vec<Finding> {
    findings_linked(ncontigs, cands, regions, merge, site_join, 0, &[], &[])
}

/// Base grouping of the tiered calls: one vector of members per base group,
/// region-type groups first, then the site-only groups, in the order findings
/// are emitted (before the final sort). `merge`/`site_join` are already widened.
fn base_groups(
    ncontigs: usize,
    cands: &[Candidate],
    regions: &[Region],
    merge: u64,
    site_join: u64,
) -> Vec<Vec<MemberIv>> {
    let mut sites: Vec<MemberIv> = Vec::new();
    for c in cands {
        // Untiered calls are not members; gray calls (homopolymer-ambiguous
        // indels, below yellow) never create or join a finding either.
        if matches!(c.tier.tier, Tier::None | Tier::Gray) {
            continue;
        }
        let len = (c.ref_allele.len() as u64).max(1);
        sites.push(MemberIv {
            tid: c.tid,
            start: c.pos0 as u64,
            end: c.pos0 as u64 + len,
            tier: c.tier.tier,
            m: Member::Site { vtype: c.vtype, af: c.af(), len: c.len, lowcov: c.lowcov & 1 != 0 },
        });
    }
    let mut rms: Vec<MemberIv> = Vec::new();
    for r in regions {
        // `depth_segment` (the paired interval) is output only: never a member.
        if r.label == "depth_segment" {
            continue;
        }
        let tier = match r.tier {
            "red" => Tier::Red,
            "orange" => Tier::Orange,
            "yellow" => Tier::Yellow,
            _ => continue,
        };
        rms.push(MemberIv {
            tid: r.tid,
            start: r.start,
            end: r.end,
            tier,
            m: Member::Region {
                label: r.label,
                reason: r.reason.clone(),
                lowcov: r.lowcov,
                depth_high: r.extra.contains("depth_high=1"),
            },
        });
    }
    debug_assert!(sites.iter().chain(rms.iter()).all(|m| m.tid < ncontigs));
    let by_pos = |a: &MemberIv, b: &MemberIv| (a.tid, a.start, a.end).cmp(&(b.tid, b.start, b.end));
    sites.sort_by(by_pos);
    rms.sort_by(by_pos);

    // Region-type members form groups by single linkage within `merge`.
    let mut groups: Vec<Vec<MemberIv>> = Vec::new();
    let mut i = 0;
    while i < rms.len() {
        let tid = rms[i].tid;
        let mut end = rms[i].end;
        let mut j = i + 1;
        while j < rms.len() && rms[j].tid == tid && rms[j].start <= end + merge {
            end = end.max(rms[j].end);
            j += 1;
        }
        groups.push(rms[i..j].to_vec());
        i = j;
    }

    // Every region member with its group index, sorted by (tid, start, end), so
    // a site is matched against the actual region intervals (not the union span
    // of a group, which would cover the gaps between its members).
    let mut r_ivs: Vec<(usize, u64, u64, usize)> = Vec::new();
    for (gi, members) in groups.iter().enumerate() {
        for m in members {
            r_ivs.push((m.tid, m.start, m.end, gi));
        }
    }
    r_ivs.sort_by(|a, b| (a.0, a.1, a.2).cmp(&(b.0, b.1, b.2)));

    // A site joins the first region group it is within `site_join` of; sites
    // that reach none are grouped among themselves below.
    let mut free_sites: Vec<MemberIv> = Vec::new();
    for s in sites {
        let mut attached = false;
        for &(tid, rs, re, gi) in &r_ivs {
            if tid < s.tid {
                continue;
            }
            if tid > s.tid || rs > s.end + site_join {
                break;
            }
            if re + site_join >= s.start {
                groups[gi].push(s.clone());
                attached = true;
                break;
            }
        }
        if !attached {
            free_sites.push(s);
        }
    }

    // Remaining sites join each other within `merge`.
    let mut i = 0;
    while i < free_sites.len() {
        let tid = free_sites[i].tid;
        let mut end = free_sites[i].end;
        let mut j = i + 1;
        while j < free_sites.len() && free_sites[j].tid == tid && free_sites[j].start <= end + merge {
            end = end.max(free_sites[j].end);
            j += 1;
        }
        groups.push(free_sites[i..j].to_vec());
        i = j;
    }
    groups
}

/// Disjoint-set forest over group and link nodes (smallest root wins).
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind { parent: (0..n).collect() }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra.max(rb)] = ra.min(rb);
        }
    }
}

/// The arc(s) a link covers on a contig, as up to two non-wrapping `[lo, hi)`
/// pieces. On a circular contig the shorter way round is taken; a span crossing
/// the origin becomes two pieces. Linear contigs (and zero-length ones) give a
/// single `[min, max]` piece.
fn arc_pieces(tid: usize, a: u64, b: u64, lens: &[u64], circular: &[bool]) -> Vec<(u64, u64)> {
    let (lo, hi) = (a.min(b), a.max(b));
    let len = lens.get(tid).copied().unwrap_or(0);
    if !circular.get(tid).copied().unwrap_or(false) || len == 0 {
        return vec![(lo, hi)];
    }
    let fwd = hi - lo;
    let rev = len.saturating_sub(fwd);
    if fwd <= rev { vec![(lo, hi)] } else { vec![(hi, len), (0, lo)] }
}

/// Minimal arc covering every member interval on a circular contig, written in
/// the `start < len < end` convention (the complement of the largest gap; a
/// tie prefers the wrap gap, so a non-wrapping set keeps `[min, max]`).
fn covering_arc(ms: &[MemberIv], len: u64) -> Option<(u64, u64)> {
    if len == 0 || ms.is_empty() {
        return None;
    }
    let mut ivs: Vec<(u64, u64)> = ms.iter().map(|m| (m.start, m.end)).collect();
    ivs.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for (s, e) in ivs {
        match merged.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => merged.push((s, e)),
        }
    }
    if merged.len() == 1 {
        return Some(merged[0]);
    }
    let k = merged.len();
    let mut best_gap = 0u64;
    let mut best_idx = 0usize;
    for i in 0..k - 1 {
        let gap = merged[i + 1].0.saturating_sub(merged[i].1);
        if gap > best_gap {
            best_gap = gap;
            best_idx = i;
        }
    }
    let wrap_gap = len.saturating_sub(merged[k - 1].1).saturating_add(merged[0].0);
    if wrap_gap >= best_gap {
        return Some((merged[0].0, merged[k - 1].1.min(len)));
    }
    let start = merged[best_idx + 1].0;
    Some((start, start + len.saturating_sub(best_gap)))
}

/// Group the tiered calls into findings, then pool the groups connected by the
/// event links of `regions`: a `split_junction` cluster's two breakpoints and a
/// `depth_segment`'s two ends. A link joins every base group whose interval
/// overlaps its span (any overlap, option b), so every tiered call inside the
/// span is absorbed; overlapping link spans chain transitively. `link_max == 0`
/// disables linking and reproduces the base grouping exactly. `lens` and
/// `circular` (indexed by contig id) drive the circular shorter-way distance
/// and the wrapped span of a finding.
#[allow(clippy::too_many_arguments)]
pub(crate) fn findings_linked(
    ncontigs: usize,
    cands: &[Candidate],
    regions: &[Region],
    merge: u32,
    site_join: u32,
    link_max: u32,
    lens: &[u64],
    circular: &[bool],
) -> Vec<Finding> {
    let groups = base_groups(ncontigs, cands, regions, merge as u64, site_join as u64);
    let finish = |groups: &[Vec<MemberIv>]| -> Vec<Finding> {
        let mut out: Vec<Finding> = groups.iter().map(|g| build_finding(g[0].tid, g)).collect();
        out.sort_by(|a, b| (a.tid, a.start, a.end).cmp(&(b.tid, b.start, b.end)));
        out
    };
    if link_max == 0 {
        return finish(&groups);
    }

    // Event links, deduplicated (a `split_junction` pair and a wrapping
    // `depth_segment` each appear once per region).
    let mut links: Vec<(usize, u64, u64)> = Vec::new();
    for r in regions {
        if let Some((a, b)) = r.link {
            links.push((r.tid, a.min(b), a.max(b)));
        }
    }
    links.sort_unstable();
    links.dedup();

    // Keep only links within `link_max` (shorter way round on circular contigs).
    let link_max = link_max as u64;
    let dist = |tid: usize, lo: u64, hi: u64| -> u64 {
        let len = lens.get(tid).copied().unwrap_or(0);
        let fwd = hi - lo;
        if circular.get(tid).copied().unwrap_or(false) && len > 0 { fwd.min(len - fwd) } else { fwd }
    };
    links.retain(|&(tid, lo, hi)| dist(tid, lo, hi) <= link_max);
    if links.is_empty() {
        return finish(&groups);
    }

    let mut by_tid: Vec<Vec<usize>> = vec![Vec::new(); ncontigs];
    let mut spans: Vec<(u64, u64)> = Vec::with_capacity(groups.len());
    for (gi, g) in groups.iter().enumerate() {
        by_tid[g[0].tid].push(gi);
        spans.push((g.iter().map(|m| m.start).min().unwrap(), g.iter().map(|m| m.end).max().unwrap()));
    }
    let mut uf = UnionFind::new(groups.len() + links.len());
    let link_pieces: Vec<Vec<(u64, u64)>> =
        links.iter().map(|&(tid, lo, hi)| arc_pieces(tid, lo, hi, lens, circular)).collect();
    for (li, &(tid, _, _)) in links.iter().enumerate() {
        let node = groups.len() + li;
        for &(plo, phi) in &link_pieces[li] {
            for &gi in &by_tid[tid] {
                let (gs, ge) = spans[gi];
                if gs <= phi && ge >= plo {
                    uf.union(gi, node);
                }
            }
        }
    }
    // Overlapping link spans chain transitively (same contig).
    for (i, &(ti, _, _)) in links.iter().enumerate() {
        for (j, &(tj, _, _)) in links.iter().enumerate().skip(i + 1) {
            if ti != tj {
                continue;
            }
            let pi = &link_pieces[i];
            let pj = &link_pieces[j];
            if pi.iter().any(|&(a, b)| pj.iter().any(|&(c, d)| a <= d && c <= b)) {
                uf.union(groups.len() + i, groups.len() + j);
            }
        }
    }

    // Components in base-group order (deterministic).
    let mut roots: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for gi in 0..groups.len() {
        let r = uf.find(gi);
        match roots.get(&r) {
            Some(&k) => comps[k].push(gi),
            None => {
                roots.insert(r, comps.len());
                comps.push(vec![gi]);
            }
        }
    }
    let mut out: Vec<Finding> = Vec::with_capacity(comps.len());
    for gis in &comps {
        let tid = groups[gis[0]][0].tid;
        let mut members: Vec<MemberIv> = Vec::new();
        for &gi in gis {
            members.extend(groups[gi].iter().cloned());
        }
        let mut f = build_finding(tid, &members);
        if circular.get(tid).copied().unwrap_or(false)
            && let Some((s, e)) = covering_arc(&members, lens.get(tid).copied().unwrap_or(0))
        {
            f.start = s;
            f.end = e;
        }
        out.push(f);
    }
    out.sort_by(|a, b| (a.tid, a.start, a.end).cmp(&(b.tid, b.start, b.end)));
    out
}

/// Build one finding from its members (`ms` non-empty, one contig, any order).
fn build_finding(tid: usize, ms: &[MemberIv]) -> Finding {
    let start = ms.iter().map(|m| m.start).min().unwrap();
    let end = ms.iter().map(|m| m.end).max().unwrap();
    let mut tier = Tier::Yellow;
    let mut lowcov = true;
    let mut ev = Evidence { only_sites: true, ..Evidence::default() };
    let mut counts: Vec<(&'static str, u32)> = ALERT_TYPES.iter().map(|t| (*t, 0u32)).collect();
    let mut has_clip_left = false;
    let mut has_clip_right = false;
    let mut clips_both = false;
    let mut score = 0u32;
    for m in ms {
        score += 1;
        tier = tier.max(m.tier);
        lowcov &= m.lowcov();
        match &m.m {
            Member::Site { vtype, af, len, .. } => {
                if *af > 0.6 {
                    ev.site_major = true;
                } else if *af >= 0.4 {
                    ev.site_mixed = true;
                } else {
                    ev.site_minor = true;
                }
                match vtype {
                    VarType::Ins if *len >= 50 => ev.long_ins = true,
                    VarType::Del if *len >= 50 => ev.long_del = true,
                    _ => {}
                }
                if let Some(s) = counts.iter_mut().find(|(l, _)| *l == vtype.name()) {
                    s.1 += 1;
                }
            }
            Member::Region { label, reason, depth_high, .. } => {
                ev.only_sites = false;
                match *label {
                    "large_ins" => ev.long_ins = true,
                    "large_del" => ev.long_del = true,
                    "clip_hotspot" => {
                        ev.clips = true;
                        if reason == "clip_both" {
                            clips_both = true;
                        }
                        if reason == "clip_left" {
                            has_clip_left = true;
                        }
                        if reason == "clip_right" {
                            has_clip_right = true;
                        }
                    }
                    "depth_low" => ev.depth_low = true,
                    "depth_half" => ev.depth_half = true,
                    "depth_high" => ev.depth_high = true,
                    "depth_step" => ev.depth_step = true,
                    "snp_cluster" => ev.cluster = true,
                    "linked_haplotype" => {
                        ev.linked = true;
                        // A group red from its own depth ratio / high-MAPQ share
                        // carries the `depth_high` evidence of its finding.
                        if *depth_high {
                            ev.depth_high = true;
                        }
                    }
                    // A split-read junction is clip-style structural evidence
                    // (the `likely` hint is unchanged: `clips` -> `misjoin`).
                    "split_junction" => ev.clips = true,
                    _ => {}
                }
                if let Some(s) = counts.iter_mut().find(|(l, _)| *l == *label) {
                    s.1 += 1;
                }
            }
        }
    }
    if clips_both || (has_clip_left && has_clip_right) {
        ev.clips_both = true;
    }
    counts.retain(|(_, c)| *c > 0);
    let likely = likely_of(&ev);
    Finding { tid, start, end, tier: tier.name(), lowcov, score, likely, evidence: ev, members: counts }
}

/// Findings per tier and `likely` (PASS findings; lowcov ones split out).
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct FindingCounts {
    pub red: u64,
    pub orange: u64,
    pub yellow: u64,
    pub by_likely: std::collections::BTreeMap<String, u64>,
}

impl FindingCounts {
    fn bump(&mut self, tier: &str, likely: &str) {
        match tier {
            "red" => self.red += 1,
            "orange" => self.orange += 1,
            "yellow" => self.yellow += 1,
            _ => {}
        }
        *self.by_likely.entry(likely.to_string()).or_default() += 1;
    }
}

/// Finding counts (per contig or whole run), split like [`Alerts`].
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct Findings {
    #[serde(flatten)]
    pub counts: FindingCounts,
    pub lowcov: FindingCounts,
}

impl Findings {
    fn bump(&mut self, tier: &str, likely: &str, lowcov: bool) {
        if lowcov {
            self.lowcov.bump(tier, likely);
        } else {
            self.counts.bump(tier, likely);
        }
    }
}

/// Finding counts per contig and top-level (the per-contig counts sum to the
/// top-level ones).
pub(crate) fn collect_findings_counts(ncontigs: usize, findings: &[Finding]) -> (Vec<Findings>, Findings) {
    let mut per: Vec<Findings> = (0..ncontigs).map(|_| Findings::default()).collect();
    let mut top = Findings::default();
    for f in findings {
        per[f.tid].bump(f.tier, f.likely, f.lowcov);
        top.bump(f.tier, f.likely, f.lowcov);
    }
    (per, top)
}

/// BED pieces of one finding: a single `[start, end)` row, or two rows
/// `[start, len)` and `[0, end - len)` when the finding wraps the origin
/// (`end > len`, the `start < len < end` convention). Both rows keep the same
/// finding id (the caller's label), like a wrapped `depth_segment`/`snp_cluster`.
pub(crate) fn finding_rows(start: u64, end: u64, len: u64) -> Vec<(u64, u64)> {
    if end > len { vec![(start, len), (0, end - len)] } else { vec![(start, end)] }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::qc::test_util::*;

    /// A site call at `pos0` with the given type, length and alt/dp (AF).
    fn fsite(tid: usize, pos0: u32, vtype: VarType, len: u32, alt: u32, dp: u32, tier: Tier, lowcov: u8) -> Candidate {
        let ref_allele = match vtype {
            VarType::Snp | VarType::Ins => vec![b'A'],
            VarType::Del => vec![b'A'; len as usize + 1],
        };
        Candidate {
            tid,
            pos0,
            origin: false,
            vtype,
            len,
            ref_allele,
            alt_allele: vec![b'C'],
            dp,
            dpf: dp / 2,
            dpr: dp - dp / 2,
            altf: alt,
            altr: 0,
            reff: 0,
            refr: 0,
            hp: 1,
            str_span: 0,
            mq0f: 0.0,
            class: 0,
            grouped: 1,
            ln_p: 0.0,
            q: 1.0,
            wilson: 0.0,
            strand_ok: true,
            significant: true,
            lowcov,
            linked_group: None,
            unlinked: false,
            tier: TierOut { tier, reason: "minority", mixed: false, strand: false },
        }
    }
    /// A tiered region carrying an event link (its two ends).
    fn lreg(
        tid: usize,
        label: &'static str,
        tier: &'static str,
        start: u64,
        end: u64,
        link: Option<(u64, u64)>,
    ) -> Region {
        let mut r = freg(tid, label, tier, start, end, false);
        r.link = link;
        r
    }
    #[test]
    fn finding_grouping_merge_and_tiers() {
        // 999 bp apart (gap 999) joins.
        let a = fsite(0, 0, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let b = fsite(0, 1_000, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let f = collect_findings(1, &[a.clone(), b], &[], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end, f[0].score), (0, 1001, 2));
        // 1001 bp apart does not.
        let c = fsite(0, 1_002, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let f = collect_findings(1, &[a.clone(), c], &[], 1000, 50);
        assert_eq!(f.len(), 2, "{f:?}");
        // Transitive chain (single linkage): 0..900..1800 all join.
        let chain: Vec<Candidate> = (0..3)
            .map(|k| fsite(0, 900 * k, VarType::Snp, 1, 10, 10, Tier::Yellow, 0))
            .collect();
        let f = collect_findings(1, &chain, &[], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end, f[0].score), (0, 1801, 3));
        // Tier = highest member tier; lowcov only when all members are lowcov.
        let red = fsite(0, 0, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let yel_lc = fsite(0, 100, VarType::Snp, 1, 10, 10, Tier::Yellow, 1);
        let f = collect_findings(1, &[red.clone(), yel_lc.clone()], &[], 1000, 50);
        assert_eq!(f.len(), 1);
        assert_eq!((f[0].tier, f[0].lowcov), ("red", false));
        let f = collect_findings(1, &[yel_lc.clone(), fsite(0, 100, VarType::Snp, 1, 10, 10, Tier::Red, 1)], &[], 1000, 50);
        assert_eq!((f[0].tier, f[0].lowcov), ("red", true));
        // A region and a site on different contigs never join.
        let f = collect_findings(2, &[red], &[freg(1, "large_del", "red", 0, 10, false)], 1000, 50);
        assert_eq!(f.len(), 2, "{f:?}");
    }
    #[test]
    fn finding_likely_rules_and_order() {
        let site = |pos0, vtype, len, alt, dp, tier| fsite(0, pos0, vtype, len, alt, dp, tier, 0);
        // Each rule fires on its own evidence, and the order is respected.
        let cases: Vec<(Vec<Candidate>, Vec<Region>, &str)> = vec![
            (vec![], vec![freg(0, "depth_high", "orange", 0, 10, false)], "collapsed_repeat"),
            (vec![], vec![freg(0, "snp_cluster", "orange", 0, 10, false)], "collapsed_repeat"),
            // depth_high beats a cluster member.
            (
                vec![],
                vec![
                    freg(0, "depth_high", "orange", 0, 10, false),
                    freg(0, "snp_cluster", "orange", 0, 10, false),
                ],
                "collapsed_repeat",
            ),
            (vec![], vec![freg(0, "depth_half", "orange", 0, 10, false)], "false_duplication"),
            (vec![], vec![freg(0, "depth_low", "red", 0, 10, false)], "extra_sequence"),
            (vec![], vec![freg(0, "large_del", "red", 0, 10, false)], "extra_sequence"),
            (vec![site(0, VarType::Del, 60, 10, 10, Tier::Red)], vec![], "extra_sequence"),
            (vec![], vec![freg(0, "large_ins", "red", 0, 10, false)], "missing_sequence"),
            (vec![site(0, VarType::Ins, 60, 10, 10, Tier::Red)], vec![], "missing_sequence"),
            // clip_hotspot and split_junction (both folded into `clips`) -> misjoin.
            (vec![], vec![freg(0, "clip_hotspot", "red", 0, 10, false)], "misjoin"),
            (vec![], vec![freg(0, "split_junction", "red", 0, 10, false)], "misjoin"),
            // A `depth_step` as the only evidence -> depth_change; a site call
            // or clip evidence keeps the earlier label.
            (vec![], vec![freg(0, "depth_step", "red", 0, 10, false)], "depth_change"),
            (
                vec![site(0, VarType::Snp, 1, 2, 10, Tier::Yellow)],
                vec![freg(0, "depth_step", "red", 0, 10, false)],
                "unclear",
            ),
            (
                vec![],
                vec![freg(0, "depth_step", "red", 0, 10, false), freg(0, "clip_hotspot", "red", 0, 10, false)],
                "misjoin",
            ),
            // Conflicting direction (both long_ins and long_del) -> unclear,
            // before the long_del/long_ins rules.
            (vec![], vec![freg(0, "large_ins", "red", 0, 10, false), freg(0, "large_del", "red", 0, 10, false)], "unclear"),
            (
                vec![site(0, VarType::Ins, 60, 10, 10, Tier::Red), site(5, VarType::Del, 60, 10, 10, Tier::Red)],
                vec![],
                "unclear",
            ),
            // depth_low is behind the conflict rule too.
            (
                vec![],
                vec![
                    freg(0, "large_ins", "red", 0, 10, false),
                    freg(0, "large_del", "red", 0, 10, false),
                    freg(0, "depth_low", "red", 0, 10, false),
                ],
                "unclear",
            ),
            // Rule 1 still beats the conflict rule.
            (
                vec![],
                vec![
                    freg(0, "depth_high", "orange", 0, 10, false),
                    freg(0, "large_ins", "red", 0, 10, false),
                    freg(0, "large_del", "red", 0, 10, false),
                ],
                "collapsed_repeat",
            ),
            // Only site calls: major -> base_error, mixed -> mixed, minor -> minority.
            (vec![site(0, VarType::Snp, 1, 7, 10, Tier::Red)], vec![], "base_error"),
            (vec![site(0, VarType::Snp, 1, 5, 10, Tier::Orange)], vec![], "mixed"),
            (vec![site(0, VarType::Snp, 1, 2, 10, Tier::Yellow)], vec![], "minority"),
            // A region plus a site is not "only sites": falls through to unclear
            // only when no rule matches; large_del sets long_del, so use the
            // unit-level `likely_of` for the unclear fallback.
        ];
        for (i, (cands, regions, want)) in cases.iter().enumerate() {
            let f = collect_findings(1, cands, regions, 1000, 50);
            assert_eq!(f.len(), 1, "case {i}: {f:?}");
            assert_eq!(f[0].likely, *want, "case {i}: {f:?}");
        }
        // Rule 9 fallback.
        assert_eq!(likely_of(&Evidence::default()), "unclear");
        // AF boundaries: > 0.6 major, 0.4-0.6 mixed, < 0.4 minor.
        let mk = |alt: u32| {
            let f = collect_findings(1, &[fsite(0, 0, VarType::Snp, 1, alt, 100, Tier::Yellow, 0)], &[], 1000, 50);
            f[0].evidence
        };
        assert!(mk(61).site_major && !mk(61).site_mixed);
        assert!(mk(60).site_mixed && !mk(60).site_major);
        assert!(mk(40).site_mixed && !mk(40).site_minor);
        assert!(mk(39).site_minor && !mk(39).site_mixed);
    }
    #[test]
    fn finding_evidence_from_each_label() {
        let ev = |label: &'static str, reason: &str| {
            let mut r = freg(0, label, "red", 0, 10, false);
            r.reason = reason.to_string();
            collect_findings(1, &[], &[r], 1000, 50).remove(0).evidence
        };
        assert!(ev("large_ins", "large_indel").long_ins);
        assert!(ev("large_del", "large_indel").long_del);
        assert!(ev("depth_low", "low").depth_low);
        assert!(ev("depth_half", "half").depth_half);
        assert!(ev("depth_high", "high").depth_high);
        assert!(ev("snp_cluster", "minority_snp_cluster").cluster);
        // A linked group red from its own depth ratio carries depth_high.
        let mut lh = freg(0, "linked_haplotype", "red", 0, 100, false);
        lh.extra = "group=1;sites=2;reads=10;depth_high=1".to_string();
        let f = collect_findings(1, &[], &[lh], 1000, 50);
        assert!(f[0].evidence.linked && f[0].evidence.depth_high, "{:?}", f[0].evidence);
        assert_eq!(f[0].likely, "collapsed_repeat");
        // An orange linked group (no depth evidence) is a subpopulation.
        let mut lh = freg(0, "linked_haplotype", "orange", 0, 100, false);
        lh.extra = "group=1;sites=2;reads=10".to_string();
        let f = collect_findings(1, &[], &[lh], 1000, 50);
        assert!(f[0].evidence.linked && !f[0].evidence.depth_high, "{:?}", f[0].evidence);
        assert_eq!(f[0].likely, "subpopulation");
        let c = ev("clip_hotspot", "clip_both");
        assert!(c.clips && c.clips_both);
        let c = ev("clip_hotspot", "clip_left");
        assert!(c.clips && !c.clips_both);
        // A left and a right hotspot in the same finding set clips_both.
        let mut left = freg(0, "clip_hotspot", "red", 0, 10, false);
        left.reason = "clip_left".to_string();
        let mut right = freg(0, "clip_hotspot", "red", 20, 30, false);
        right.reason = "clip_right".to_string();
        let f = collect_findings(1, &[], &[left, right], 1000, 50);
        assert!(f[0].evidence.clips_both, "{:?}", f[0].evidence);
        // Member counts and the evidence string.
        let f = collect_findings(
            1,
            &[fsite(0, 0, VarType::Snp, 1, 10, 10, Tier::Red, 0), fsite(0, 5, VarType::Del, 60, 10, 10, Tier::Red, 0)],
            &[],
            1000, 50,
        );
        assert_eq!(f[0].members, vec![("SNP", 1), ("DEL", 1)], "{:?}", f[0].members);
        assert_eq!(f[0].evidence_list(), "site_major,long_del");
        assert_eq!(f[0].likely, "extra_sequence");
    }
    #[test]
    fn finding_ref_span_and_merge_boundary() {
        // A DEL's interval is its REF span: a 100 bp deletion at 0 covers
        // [0, 101). Anchor-only would leave a 1099 bp gap; the REF span leaves
        // 999, so the two calls join.
        let del = fsite(0, 0, VarType::Del, 100, 10, 10, Tier::Red, 0);
        let snp = fsite(0, 1_100, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let f = collect_findings(1, &[del, snp], &[], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end), (0, 1_101));
        // Exact merge boundary: starts 1001 bp apart -> gap exactly 1000 -> join.
        let a = fsite(0, 0, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let b = fsite(0, 1_001, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        assert_eq!(collect_findings(1, &[a.clone(), b], &[], 1000, 50).len(), 1);
        // One bp more -> gap 1001 -> split.
        let c = fsite(0, 1_002, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        assert_eq!(collect_findings(1, &[a.clone(), c], &[], 1000, 50).len(), 2);
        // A custom merge of 0 joins only touching/overlapping intervals.
        let d = fsite(0, 1, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        assert_eq!(collect_findings(1, &[a.clone(), d], &[], 0, 50).len(), 1);
        let e = fsite(0, 2, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        assert_eq!(collect_findings(1, &[a, e], &[], 0, 50).len(), 2);
    }
    #[test]
    fn finding_site_region_join_distance() {
        let region = freg(0, "clip_hotspot", "yellow", 1_000, 1_100, false);
        let site = |pos0| fsite(0, pos0, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);

        // A site 49 bp from the region end joins it; 51 bp does not.
        let f = collect_findings(1, &[site(1_149)], &[region.clone()], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!(f[0].score, 2);
        let f = collect_findings(1, &[site(1_151)], &[region.clone()], 1000, 50);
        assert_eq!(f.len(), 2, "{f:?}");

        // 600 bp away is its own finding (the old 1000 bp rule merged it).
        let f = collect_findings(1, &[site(1_700)], &[region.clone()], 1000, 50);
        assert_eq!(f.len(), 2, "{f:?}");
        assert_eq!(f[0].score, 1);

        // Region-type members still join each other within 1000 bp, even when a
        // site-only component sits in between: the regions form one finding and
        // the site (200 bp from each) is its own.
        let far = freg(0, "depth_step", "red", 1_300, 1_350, false);
        let f = collect_findings(1, &[site(1_200)], &[region.clone(), far], 1000, 50);
        let spans: Vec<_> = f.iter().map(|f| (f.start, f.end, f.score)).collect();
        assert_eq!(spans, vec![(1_000, 1_350, 2), (1_200, 1_201, 1)], "{f:?}");

        // Site-site grouping is unchanged.
        assert_eq!(collect_findings(1, &[site(0), site(999)], &[], 1000, 50).len(), 1);
        assert_eq!(collect_findings(1, &[site(0), site(1_002)], &[], 1000, 50).len(), 2);

        // Independent of input order.
        let cands = vec![site(1_149), site(1_700), site(3_000)];
        let regions = vec![region, freg(0, "depth_low", "red", 2_000, 2_100, false)];
        let key = |f: &Finding| {
            (f.tid, f.start, f.end, f.tier, f.lowcov, f.score, f.likely, f.evidence_list(), f.members.clone())
        };
        let a: Vec<_> = collect_findings(1, &cands, &regions, 1000, 50).iter().map(key).collect();
        let mut rc = cands;
        rc.reverse();
        let mut rr = regions;
        rr.reverse();
        let b: Vec<_> = collect_findings(1, &rc, &rr, 1000, 50).iter().map(key).collect();
        assert_eq!(a, b);
    }
    #[test]
    fn finding_site_joins_leftmost_eligible_region_group() {
        // Two region groups far enough apart to stay separate (> merge): a long
        // site call reaches both, and the leftmost group must be chosen.
        let left = freg(0, "clip_hotspot", "red", 0, 10, false);
        let right = freg(0, "depth_step", "red", 2_000, 2_010, false);
        let long_site = fsite(0, 0, VarType::Del, 2_000, 10, 10, Tier::Yellow, 0);
        let key = |f: &Finding| (f.start, f.end, f.score);
        let a: Vec<_> = collect_findings(1, &[long_site.clone()], &[left.clone(), right.clone()], 1000, 50)
            .iter().map(key).collect();
        // The site joins the left group; the right region is its own finding.
        assert_eq!(a, vec![(0, 2_001, 2), (2_000, 2_010, 1)], "{a:?}");
        // Same with shuffled input.
        let b: Vec<_> = collect_findings(1, &[long_site], &[right, left], 1000, 50)
            .iter().map(key).collect();
        assert_eq!(a, b);
    }
    #[test]
    fn finding_attached_site_does_not_pull_free_site_chain() {
        // A region + an attached site are one finding; a free-site chain nearby
        // stays separate instead of being pulled in.
        let region = freg(0, "clip_hotspot", "red", 0, 10, false);
        let attached = fsite(0, 20, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let free1 = fsite(0, 1_000, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let free2 = fsite(0, 2_000, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let spans: Vec<_> = collect_findings(
            1, &[attached, free1, free2], &[region], 1000, 50)
            .iter().map(|f| (f.start, f.end, f.score)).collect();
        assert_eq!(spans, vec![(0, 21, 2), (1_000, 2_001, 2)], "{spans:?}");
    }
    #[test]
    fn finding_mixed_tier_and_lowcov_members() {
        // tier = highest member; lowcov only when every member is lowcov.
        let red_lc = freg(0, "clip_hotspot", "red", 0, 10, true);
        let yellow = fsite(0, 20, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let f = collect_findings(1, &[yellow], &[red_lc], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].tier, f[0].lowcov, f[0].score), ("red", false, 2));
        // All lowcov -> the finding is lowcov.
        let yel_lc_r = freg(0, "depth_step", "yellow", 5_000, 5_010, true);
        let yel_lc_s = fsite(0, 5_020, VarType::Snp, 1, 10, 10, Tier::Yellow, 1);
        let f = collect_findings(1, &[yel_lc_s], &[yel_lc_r], 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].tier, f[0].lowcov, f[0].score), ("yellow", true, 2));
    }
    #[test]
    fn finding_site_join_exact_boundary() {
        // A 50 bp gap joins exactly; 51 bp does not.
        let region = freg(0, "clip_hotspot", "red", 1_000, 1_100, false);
        let at = fsite(0, 1_150, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let over = fsite(0, 1_151, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        assert_eq!(collect_findings(1, &[at], &[region.clone()], 1000, 50).len(), 1);
        assert_eq!(collect_findings(1, &[over], &[region], 1000, 50).len(), 2);
    }
    #[test]
    fn finding_untiered_excluded_and_order_independent() {
        // Untiered annotations (tier `none`) are not members and do not group
        // with a nearby tiered call.
        let ann = vec![
            freg(0, "not_assessable", "none", 0, 100, false),
            freg(0, "low_coverage", "none", 200, 300, false),
        ];
        assert!(collect_findings(1, &[], &ann, 1000, 50).is_empty());
        let site = fsite(0, 150, VarType::Snp, 1, 10, 10, Tier::Red, 0);
        let f = collect_findings(1, &[site], &ann, 1000, 50);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!(f[0].score, 1);
        // Shuffled input gives the same findings.
        let cands = vec![
            fsite(0, 0, VarType::Snp, 1, 10, 10, Tier::Yellow, 0),
            fsite(0, 900, VarType::Del, 60, 10, 10, Tier::Red, 1),
            fsite(0, 1_800, VarType::Ins, 60, 10, 10, Tier::Orange, 0),
        ];
        let key = |f: &Finding| {
            (f.tid, f.start, f.end, f.tier, f.lowcov, f.score, f.likely, f.evidence_list(), f.members.clone())
        };
        let a: Vec<_> = collect_findings(1, &cands, &[], 1000, 50).iter().map(key).collect();
        let mut rev = cands;
        rev.reverse();
        let b: Vec<_> = collect_findings(1, &rev, &[], 1000, 50).iter().map(key).collect();
        assert_eq!(a, b);
        // The three chain (gaps 899) into one finding with conflicting indel
        // evidence (long_ins + long_del) -> unclear; the mixed lowcov does not
        // make the finding lowcov.
        assert_eq!(a.len(), 1, "{a:?}");
        assert_eq!(a[0].6, "unclear");
        assert!(!a[0].4);
    }
    #[test]
    fn finding_indel_length_evidence_boundaries() {
        let ev = |vtype, len| collect_findings(1, &[fsite(0, 0, vtype, len, 10, 10, Tier::Red, 0)], &[], 1000, 50).remove(0);
        // 49 bp is not long, 50 bp is.
        assert!(!ev(VarType::Ins, 49).evidence.long_ins);
        assert!(ev(VarType::Ins, 50).evidence.long_ins);
        assert!(!ev(VarType::Del, 49).evidence.long_del);
        assert!(ev(VarType::Del, 50).evidence.long_del);
        // The hint follows the boundary.
        assert_eq!(ev(VarType::Ins, 49).likely, "base_error");
        assert_eq!(ev(VarType::Ins, 50).likely, "missing_sequence");
        assert_eq!(ev(VarType::Del, 49).likely, "base_error");
        assert_eq!(ev(VarType::Del, 50).likely, "extra_sequence");
    }
    #[test]
    fn finding_counts_split_and_sum() {
        let cands = vec![
            fsite(0, 0, VarType::Snp, 1, 7, 10, Tier::Red, 0),
            fsite(0, 5_000, VarType::Snp, 1, 2, 10, Tier::Yellow, 1),
            fsite(1, 0, VarType::Snp, 1, 5, 10, Tier::Orange, 0),
        ];
        let regions = vec![freg(1, "large_del", "red", 0, 10, true)];
        let findings = collect_findings(2, &cands, &regions, 1000, 50);
        let (per, top) = collect_findings_counts(2, &findings);
        // Recompute from the findings themselves.
        let (mut red, mut orange, mut yellow) = (0u64, 0u64, 0u64);
        let (mut lred, mut lorange, mut lyellow) = (0u64, 0u64, 0u64);
        for f in &findings {
            if f.lowcov {
                if f.tier == "red" {
                    lred += 1
                } else if f.tier == "orange" {
                    lorange += 1
                } else {
                    lyellow += 1
                }
            } else if f.tier == "red" {
                red += 1
            } else if f.tier == "orange" {
                orange += 1
            } else {
                yellow += 1
            }
        }
        assert_eq!((top.counts.red, top.counts.orange, top.counts.yellow), (red, orange, yellow));
        assert_eq!((top.lowcov.red, top.lowcov.orange, top.lowcov.yellow), (lred, lorange, lyellow));
        // Per-contig sums to the top level, and the lowcov finding is split out.
        let sum = |f: fn(&Findings) -> u64| per.iter().map(f).sum::<u64>();
        assert_eq!(sum(|a| a.counts.red + a.counts.orange + a.counts.yellow), red + orange + yellow);
        assert_eq!(sum(|a| a.lowcov.red + a.lowcov.orange + a.lowcov.yellow), lred + lorange + lyellow);
        assert!(lred + lorange + lyellow >= 1, "{top:?}");
        // `by_likely` matches the per-finding hints, split PASS/lowcov.
        let (mut likely, mut likely_low) = (std::collections::BTreeMap::new(), std::collections::BTreeMap::new());
        for f in &findings {
            let m = if f.lowcov { &mut likely_low } else { &mut likely };
            *m.entry(f.likely.to_string()).or_insert(0u64) += 1;
        }
        assert_eq!(top.counts.by_likely, likely);
        assert_eq!(top.lowcov.by_likely, likely_low);
    }
    /// A prophage-like gain: `depth_step` up and down 30 kb apart with a SNP and
    /// a clip hotspot inside; the `depth_segment` links the two ends, so the four
    /// calls pool into one finding of the highest member tier. A small
    /// `--finding-link-max` drops the link and keeps today's separate findings.
    #[test]
    fn finding_link_pools_prophage_gain() {
        let step_up = lreg(0, "depth_step", "red", 9_950, 10_050, None);
        let step_dn = lreg(0, "depth_step", "orange", 39_950, 40_050, None);
        let seg = lreg(0, "depth_segment", "orange", 10_000, 40_000, Some((10_000, 40_000)));
        let clip = lreg(0, "clip_hotspot", "yellow", 37_950, 38_050, None);
        let snp = fsite(0, 12_000, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let far = fsite(0, 200_000, VarType::Snp, 1, 10, 10, Tier::Yellow, 0);
        let regions = vec![step_up, step_dn, seg, clip];
        let cands = vec![snp, far];
        let lens = [1_000_000u64];
        let circ = [false];
        let linked = findings_linked(1, &cands, &regions, 1000, 50, 100_000, &lens, &circ);
        assert_eq!(linked.len(), 2, "{linked:?}"); // pooled gain + far site
        let pf = &linked[0];
        assert_eq!((pf.start, pf.end, pf.tier, pf.score), (9_950, 40_050, "red", 4), "{pf:?}");
        // The 30 kb link is dropped at 10 kb -> today's separate findings.
        let sep = findings_linked(1, &cands, &regions, 1000, 50, 10_000, &lens, &circ);
        assert_eq!(sep.len(), 5, "{sep:?}");
    }
    /// A `split_junction` pair 30 kb apart pools into one finding; 150 kb apart
    /// exceeds the default `--finding-link-max` and stays two.
    #[test]
    fn finding_link_junction_pair_distance() {
        let lens = [1_000_000u64];
        let circ = [false];
        let a = lreg(0, "split_junction", "yellow", 29_800, 30_200, Some((30_000, 60_000)));
        let b = lreg(0, "split_junction", "yellow", 59_800, 60_200, Some((30_000, 60_000)));
        let f = findings_linked(1, &[], &[a, b], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end, f[0].score), (29_800, 60_200, 2));
        let c = lreg(0, "split_junction", "yellow", 59_800, 60_200, Some((60_000, 210_000)));
        let d = lreg(0, "split_junction", "yellow", 209_800, 210_200, Some((60_000, 210_000)));
        let f = findings_linked(1, &[], &[c, d], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 2, "{f:?}");
    }
    /// A link across the origin of a circular contig (the shorter way round) pools
    /// the two ends and writes the finding with `start < len < end`; on a linear
    /// contig the same pair is too far apart.
    #[test]
    fn finding_link_circular_origin_span() {
        let lo = lreg(0, "depth_step", "yellow", 50, 150, None);
        let hi = lreg(0, "depth_step", "yellow", 999_850, 999_950, None);
        let seg = lreg(0, "depth_segment", "yellow", 999_900, 100, Some((999_900, 100)));
        let regions = vec![lo, hi, seg];
        let lens = [1_000_000u64];
        let f = findings_linked(1, &[], &regions, 1000, 50, 100_000, &lens, &[true]);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end), (999_850, 1_000_150), "{f:?}");
        let f = findings_linked(1, &[], &regions, 1000, 50, 100_000, &lens, &[false]);
        assert_eq!(f.len(), 2, "{f:?}");
    }
    /// `--finding-link-max 0` reproduces the base grouping. The expected base
    /// findings are spelled out independently (not via `findings_linked(...,0)`).
    #[test]
    fn finding_link_max_zero_reproduces_base() {
        let a = lreg(0, "split_junction", "yellow", 29_800, 30_200, Some((30_000, 60_000)));
        let b = lreg(0, "split_junction", "yellow", 59_800, 60_200, Some((30_000, 60_000)));
        let snp = fsite(0, 45_000, VarType::Snp, 1, 10, 10, Tier::Orange, 0);
        let regions = vec![a, b];
        let cands = vec![snp];
        let key = |f: &Finding| {
            (f.tid, f.start, f.end, f.tier, f.likely, f.score, f.lowcov, f.members.clone())
        };
        let zero = findings_linked(1, &cands, &regions, 1000, 50, 0, &[1_000_000], &[false]);
        let got: Vec<_> = zero.iter().map(key).collect();
        // Base grouping: the two junction regions and the site stay apart.
        let want: Vec<(usize, u64, u64, &str, &str, u32, bool, Vec<(&str, u32)>)> = vec![
            (0, 29_800, 30_200, "yellow", "misjoin", 1, false, vec![("split_junction", 1)]),
            (0, 45_000, 45_001, "orange", "base_error", 1, false, vec![("SNP", 1)]),
            (0, 59_800, 60_200, "yellow", "misjoin", 1, false, vec![("split_junction", 1)]),
        ];
        assert_eq!(got, want, "{zero:?}");
        // The default does pool the link, so the zero case is not vacuous.
        let linked = findings_linked(1, &cands, &regions, 1000, 50, 100_000, &[1_000_000], &[false]);
        assert_eq!(linked.len(), 1, "{linked:?}");
    }
    /// Two links whose spans overlap chain even when no base group lies in the
    /// overlap (the links must be joined through each other).
    #[test]
    fn finding_link_chains_overlapping_spans_without_member() {
        let lens = [1_000_000u64];
        let circ = [false];
        // link1 [1000,1900], link2 [1850,3000]: overlap [1850,1900], no member.
        let l1 = lreg(0, "split_junction", "yellow", 990, 1_010, Some((1_000, 1_900)));
        let l2 = lreg(0, "split_junction", "yellow", 2_940, 2_960, Some((1_850, 3_000)));
        let f = findings_linked(1, &[], &[l1.clone(), l2], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end, f[0].score), (990, 2_960, 2));
        // Spans that do not overlap stay separate.
        let l3 = lreg(0, "split_junction", "yellow", 2_940, 2_960, Some((2_000, 3_000)));
        let f = findings_linked(1, &[], &[l1, l3], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 2, "{f:?}");
    }
    /// The link distance cap is inclusive: exactly `--finding-link-max` links,
    /// one bp more does not.
    #[test]
    fn finding_link_exact_cap_boundary() {
        let lens = [1_000_000u64];
        let circ = [false];
        let a = lreg(0, "split_junction", "yellow", 0, 20, Some((10, 100_010)));
        let b = lreg(0, "split_junction", "yellow", 100_000, 100_020, Some((10, 100_010)));
        let f = findings_linked(1, &[], &[a, b], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 1, "{f:?}"); // 100_010 - 10 == 100_000
        let c = lreg(0, "split_junction", "yellow", 0, 20, Some((10, 100_011)));
        let d = lreg(0, "split_junction", "yellow", 100_001, 100_021, Some((10, 100_011)));
        let f = findings_linked(1, &[], &[c, d], 1000, 50, 100_000, &lens, &circ);
        assert_eq!(f.len(), 2, "{f:?}"); // one bp over the cap
    }
    /// A contig shorter than `--finding-link-max`: the shorter way round still
    /// links the two ends across the origin.
    #[test]
    fn finding_link_short_circular_contig() {
        let a = lreg(0, "depth_step", "yellow", 50, 150, None);
        let b = lreg(0, "depth_step", "yellow", 4_850, 4_950, None);
        let seg = lreg(0, "depth_segment", "yellow", 4_900, 100, Some((4_900, 100)));
        let f = findings_linked(1, &[], &[a, b, seg], 1000, 50, 100_000, &[5_000], &[true]);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end), (4_850, 5_150), "{f:?}");
    }
    /// An equal largest gap must prefer the non-wrapping arc (100 kb circle,
    /// ends at 25,000 / 75,000 -> [24,900, 75,100)).
    #[test]
    fn finding_link_equal_gap_prefers_non_wrapping() {
        let a = lreg(0, "depth_step", "yellow", 24_900, 25_100, None);
        let b = lreg(0, "depth_step", "yellow", 74_900, 75_100, None);
        let seg = lreg(0, "depth_segment", "yellow", 25_000, 75_000, Some((25_000, 75_000)));
        let f = findings_linked(1, &[], &[a, b, seg], 1000, 50, 100_000, &[100_000], &[true]);
        assert_eq!(f.len(), 1, "{f:?}");
        assert_eq!((f[0].start, f[0].end), (24_900, 75_100), "{f:?}");
    }
    /// A wrapped finding is written as two rows `[start, len)` + `[0, end-len)`
    /// (same id); a non-wrapped one as a single row.
    #[test]
    fn finding_rows_split_at_origin() {
        assert_eq!(finding_rows(100, 200, 1_000), vec![(100, 200)]);
        assert_eq!(finding_rows(9_850, 10_150, 10_000), vec![(9_850, 10_000), (0, 150)]);
    }
}
