//! cramhook — lift reads from an old assembly to a new one.
//!
//! Exposes the assembly diff (M2), the input-contract check and the streaming
//! translator (`translate`, M3) as library functions so the later `lift`
//! milestone can reuse them.

pub(crate) mod bamio;
pub mod coords;
pub mod cramcount;
pub mod contract;
pub mod diff;
pub mod fasta;
pub mod fish;
pub mod lift;
pub mod paf;
pub mod progress;
pub mod qc;
pub mod scan;
pub mod translate;
pub mod version;

pub use coords::{NewIv, NewPos, OldIv, OldPos};
pub use diff::{
    Block, Diff, DiffParams, DiffRun, Preset, Summary, compute_diff, compute_diff_with, run_diff,
    run_diff_with, write_outputs,
};
pub use fasta::Fasta;
pub use paf::{PafRecord, Strand};
pub use translate::{Fallback, Params, Reason, Stats, TranslateConfig, translate};
pub use lift::{LiftOptions, lift, lift_with_progress};
pub use progress::Progress;
pub use fish::{
    BaitContig, BaitRecords, BaitSummary, CapSummary, FishParams, FishScope, FishSummary, FishTier, run_fish,
};
pub use qc::{QcParams, run_qc};
