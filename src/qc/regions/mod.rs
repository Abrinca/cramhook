//! `qc` region-signal detectors.
//!
//! Each submodule emits one family of BED regions; moved verbatim out of
//! `qc/mod.rs` (cleanup-split step 2).

mod clip;
mod cluster;
mod depth;
mod indel;
mod junction;
mod step;
mod xlink;

pub(crate) use clip::*;
pub(crate) use cluster::*;
pub(crate) use depth::*;
pub(crate) use indel::*;
pub(crate) use junction::*;
pub(crate) use step::*;
pub(crate) use xlink::*;
