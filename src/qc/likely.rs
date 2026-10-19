//! `qc` finding evidence flags and the `likely` hint.
//!
//! Moved verbatim out of `qc/mod.rs` (cleanup-split step 2).

/// Evidence flags of a finding, derived from its members (see the spec).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Evidence {
    pub site_major: bool,
    pub site_mixed: bool,
    pub site_minor: bool,
    pub long_ins: bool,
    pub long_del: bool,
    pub clips: bool,
    pub clips_both: bool,
    pub depth_low: bool,
    pub depth_half: bool,
    pub depth_high: bool,
    /// A `depth_step` region member (`depth_segment` is output only, never a
    /// member, so it never sets this).
    pub depth_step: bool,
    pub cluster: bool,
    /// A `linked_haplotype` region member.
    pub linked: bool,
    /// Every member is a VCF site call (no BED region member).
    pub only_sites: bool,
}

/// The `likely` hint from the evidence flags (first match wins).
pub(crate) fn likely_of(e: &Evidence) -> &'static str {
    if e.depth_high || e.cluster {
        "collapsed_repeat"
    } else if e.depth_half {
        "false_duplication"
    } else if e.long_ins && e.long_del {
        // Conflicting direction: the assembly both lacks and has sequence the
        // reads show. No single explanation fits, so do not guess.
        "unclear"
    } else if e.depth_low || e.long_del {
        "extra_sequence"
    } else if e.long_ins {
        "missing_sequence"
    } else if e.clips {
        "misjoin"
    } else if e.linked {
        // A coherent second haplotype at normal depth: linked without any depth
        // or structural evidence.
        "subpopulation"
    } else if e.only_sites && e.site_major {
        "base_error"
    } else if e.only_sites && e.site_mixed {
        "mixed"
    } else if e.only_sites && e.site_minor {
        "minority"
    } else if e.depth_step && !e.site_major && !e.site_mixed && !e.site_minor {
        // A depth step as the only evidence: a change in covered sequence, not
        // a site or structural signal. A depth step together with site calls
        // (not `only_sites`) stays `unclear`.
        "depth_change"
    } else {
        "unclear"
    }
}

