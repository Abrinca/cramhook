# Decisions

One entry per non-obvious decision: what, why, alternatives.

## M1

- **Own read simulator (`scripts/simulate.py`) instead of badread/pbsim3.** Deterministic, encodes
  true origin in read names, circular-aware, no new dependencies. Alternative: badread/pbsim3 via
  pixi (more realistic error profiles); revisit if synthetic results disagree with real data.
- **Feasibility probe added to M1 (`scripts/probe.py`).** Polishing-type edits every few kb could
  leave almost no read windows inside one identity block, which would trigger the abort criterion (7.1).
  The probe measures an optimistic upper bound (no taint, overlapping blocks allowed) on real data
  before M2–M4 are built.
- **Real-data `.fai`:** if the FASTA directory behind a symlink is not writable, the FASTA is copied
  into the case dir. This keeps the source datasets untouched.
- **Baseline noise reference:** `new.fa` + a 50 kb random contig `dummy_noise` (seed 12345).
  Reads on `dummy_noise` count as primary mismatches.
- **Origin-crossing split reads are in v1** (user decision after M1). In ONT FAM20838, 45–62 % of
  reads are split, mostly across the start/end of small circular plasmids. Such a read is translated
  as a unit if the contig junction is preserved (Plan 4.2). Alternative: keep remapping all split
  reads, which leaves ONT cases below the 50 % savings bar.
- **Baseline noise is zero** with the dummy-contig definition, so 7.2 effectively requires exact
  agreement. Kept as is; revisit with evidence only if M4 shows unavoidable diffs. User (after
  M2): a non-zero allowance is acceptable if necessary. Exact agreement is nice to have, not a
  hard requirement.
- **Coverage:** 50×, one seed per scenario/model (as planned). Add seeds if mismatch counts are too
  small to interpret.

- **Equivalent CIGARs count as agreement (user decision during M3).** On reversed contigs,
  translated reads match the oracle in contig/pos/strand/MAPQ, but indel placement inside repeats
  can differ: S3 has ~7 % (HiFi) and ~60 % (ONT) of reads on the reversed plasmid affected. Left-aligning
  like minimap2's `mm_fix_cigar` fixes most of them, but not all, because ksw2 on the reverse strand is not a mirror.
  Alternatives rejected: defer all '-' reads (loses the speedup on reversed contigs), or defer
  ambiguous '-' reads (≈ all ONT '-' reads). Plan 7.2 amended; compare.py reports "equivalent"
  separately.

- **Equivalence redefined (user, after M4):** same contig/pos/end/strand/MAPQ/read length/AS.
  NM and de are now reported only: on FAM14197 (both contigs reversed), 13 co-optimal alignments
  differ in NM by 1.
- **Depth counts deletions (user, after M4):** on FAM14197, depth without deletions was 99.915 %
  identical (max 2); with `samtools depth -J` it is 100.000 % (max 0). The differences were deletions
  shifted within homopolymers by equivalent alignments.

## M2

- **Novel = new positions not covered exactly once by identity-block new intervals.** This
  deliberately includes SNP and indel positions, unaligned new sequence and multiply covered
  target (an expanded/collapsed duplication). Alternative: only truly unaligned new sequence;
  rejected because an alignment ends at every SNP, so those positions are not translated anyway.
- **T2 taint is projected base-by-base, not whole-hit.** Novel intervals are extended by FLANK,
  cut into length-`WIN` windows and mapped with the read preset keeping all hits; then each hit's
  `cs` is walked and only old bases aligned to *novel* window positions are tainted. Rationale:
  a window's flanks map back to their own origin, so tainting whole hits would taint roughly
  ±WIN around every polishing edit. Only bases similar to novel sequence can make a read
  ambiguous; T2 errors also fail unsafe, so the probe is sensitive (`-N 1000000 -p 0.000001`).
- **T2 window constants (initial; tune on S6/S6b): WIN = 10 000, STEP = 2 000, FLANK = 10 000.**
  S6b (diverged IS copy, ~90 % identity) reaches ≥ 97.8 % coverage of every old IS copy with
  these values.
- **T1 is the union of two sweeps:** primary query-interval coverage (catches unaligned contigs
  and overlapping alignment ends) and identity-block coverage (catches sequence inside an
  alignment's deletion gap). Needed because asm5 often represents a collapsed duplication as a
  single alignment with a deletion gap (S10a), where interval coverage alone would report no
  taint. Alternative: interval coverage only, which leaves the deleted copy untranslated by rule 7
  but invisible in `taint.bed`; rejected because the test requires taint over the duplicated
  segment.
- **Coordinates:** old target cursor in the asm5 PAF (`-cx asm5 --cs new.fa old.fa`, query = old)
  is `tstart + t`, where `t` is the cs target offset; `Block::new_start` and the T2 projection
  both add `tstart`. For `-` alignments the cs is written target-forward with the query
  reverse-complemented, so a run at query offset `qoff` covers old `[qend-qoff-n, qend-qoff)`.

## M3

- **`MD` handling: strip it, do not disable `CRAM_OPT_DECODE_MD`.** The CRAM produced by
  `samtools sort` stores neither `MD` nor `NM` (htslib drops both as recomputable), so htslib
  regenerates both on decode. Setting `hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0)` removes `MD`
  *and* `NM` (verified empirically: `MD=0, NM=0`), which breaks the `SA` rewrite and the
  "MAPQ and tags unchanged" contract. Keeping the default decode and calling
  `Record::remove_aux(b"MD")` instead keeps `NM` (equal to minimap2's value, since both count
  the same mismatches+indels against the same reference). This is the documented fallback:
  the `@PG` contract forbids `--MD`/`--cs`, so any `MD` is a decode artifact. A record carrying
  any tag outside `NM ms AS nn tp cm s1 s2 de rl SA` (after stripping `MD`) is deferred with
  reason `tags`.
- **`rust-htslib` 1.0.2 was not on crates.io at first.** It existed only in `ref/rust-htslib`
  (same version), so `Cargo.toml` kept `rust-htslib = "=1.0.2"` and patched `[patch.crates-io]
  rust-htslib = { path = "ref/rust-htslib" }`. `ref/` stays untouched. The `libdeflate` feature
  builds fine; no feature was dropped. **Resolved 2026-10-23:** 1.0.2 is published on crates.io
  (`checksum fc6fe95f…`); the `[patch]` was dropped and the lock now points at the registry. The
  published crate lacks the unreleased post-1.0.2 `Header::to_hashmap -> Result` change that the
  vendored checkout carried; cramhook does not call `to_hashmap`. The requirement is relaxed to
  `1.0.2` (caret) so `cargo update` can pick up 1.0.x fixes; `Cargo.lock` still pins 1.0.2.
  `ref/htslib` (C source) is kept for reference.
- **`SA` rewrite reproduces minimap2's simplification, not the record CIGAR.** minimap2 writes
  each `SA` entry as `rname,pos1,strand,<clip5>S<M><I|D><clip3>S,mapq,NM;` from the *internal*
  `mm_reg1_t` (`ref/minimap2/format.c`): `l_M = min(query_span, ref_span)`, a single `I` or `D`
  for the difference, terminal soft clips, and `NM` equal to the other record's own `NM` tag.
  We recompute it from the other *translated* record's CIGAR (query/ref spans, leading/trailing
  soft clips, flag, MAPQ, `NM`); internal indels are intentionally collapsed exactly as
  minimap2 does. On S1 the rewritten `SA` is byte-identical to the original; on the real ONT
  case all 2949 rewritten supplementary records match the oracle.
- **`SA` tag order preserved with `update_aux`.** We replace in place rather than
  `remove_aux` + `push_aux`, so `SA` keeps its position in the aux block and S1 records stay
  byte-identical apart from `XL`.
- **Contract restart.** A mid-stream contract violation (a secondary record `0x100` or a
  hard-clipped primary) aborts the translate pass; the output files are reopened (truncated)
  and the CRAM is re-streamed once in all-defer mode. With `--strict` the process exits with an
  error instead. Header-level violations are detected before streaming and skip the diff.
- **Missing `QUAL` fill.** A record whose qualities are all `0xff` is treated as missing; its
  deferred FASTQ entry gets `!` for every base and is counted in `stats.qual_missing`.
- **Deferred FASTQ `QUAL` is Phred+33.** htslib stores Phred scores (0..93); the FASTQ writer
  adds 33 (or `!` when missing), and reverse-complements/reverses SEQ/QUAL for `0x10` records.
- **Junction reads are deferred unless the contig junction is preserved (`junction`).** The
  eligibility window is clipped to `[0, contig_len)`, which hid the fact that a record's soft
  clips continue across a circular contig end. Now, if the *unclipped* window
  (`[ref_start − left_clip − margin, ref_end + right_clip + margin)`) crosses a contig boundary,
  translation is only allowed when the block reaches that old boundary and maps it to the
  corresponding new boundary (mirrored for `-`). The test is factored into
  `junction_preserved()` and shared with the origin-crossing path. This fixes the S2 rotation
  and S12 mismatches (both now 100 % CIGAR agreement) and the junction reads in the PacBio case.
- **`-` translations left-align every indel (`left_align_cigar`).** minimap2 left-aligns gaps in
  reference-forward orientation, so reversing an old `-` alignment puts each indel at the right
  end of its equivalent range. After reversing we shift each deletion left while
  `ref[base before the gap] == ref[last base of the gap]` and each insertion while
  `read[base before] == read[last inserted base]`, consuming the preceding match run and growing
  the following one (never crossing another indel or a clip) — the exact `mm_fix_cigar` pass-1
  algorithm (`ref/minimap2/align.c`). Critically, shifting decrements **both** query and reference
  offsets (the preceding match shrinks in both), which an initial implementation missed and which
  left S3 at 88 %; fixing it raised S3/hifi from 88.4 % → 99.19 % and S3/ont 88.8 % → 92.9 %.
  `+` translations are not touched (already byte-identical by construction); `NM`/`AS`/`ms` are
  unaffected (spans do not change) and the invariant check is unaffected.
- **Oracle CIGAR agreement: S2/S12 and all other scenarios 100 %; S3 remains 99.19 % (hifi) /
  92.88 % (ont) — reported, assertion kept.** After the two fixes only S3 (reverse-complemented
  contig) is below the 99.9 % bar, and only in CIGAR: `contig`, `pos`, `strand` and `mapq` agree
  100 %. minimap2's `mm_fix_cigar` is a single greedy pass over the raw ksw2 alignment and is not
  idempotent (re-running it shifts ~22 % of the oracle CIGARs further), so the oracle is not the
  global leftmost placement. Our translated CIGAR is a single greedy pass over the *reversed old*
  alignment; since the banded ksw2 alignment on `revcomp(ref)` is not the exact mirror of the
  alignment on `ref` in long repeat regions, the two greedy normalisations can choose different
  equivalent placements. `mm_fix_cigar` pass 2 (adjacent `I`/`D` collapse) changes nothing, and no
  post-hoc normalisation can recover the oracle without the raw new-reference alignment. On
  S3/hifi 133/143 p1 reads already match the oracle (0/143 before left-aligning), and the
  remaining 10 have identical contig/pos/strand/MAPQ and identical op-letter CIGAR sequences with
  shifted indel lengths. Per the task the 99.9 % assertion is not weakened; the test prints the
  table and every mismatching read.
- **Contract hardening.** `parse_minimap2_cl` now validates the *whole* command line: any option
  token after the first positional (reference path) is a violation, a repeated preset in any
  spelling (`-ax P` plus `-x Q`) is a violation, and combined short-flag bundles other than the
  exact `-ax` are rejected. Repeated `-t` is still allowed. Regression tests cover
  `... ref reads -k 11`, `... ref reads --MD`, `... ref reads -y`, `-ax map-hifi -x lr:hq`,
  `-xmap-hifi -ax lr:hq`, and `-aYx`.
- **SA handling: presence tracked separately; fields validated strictly.** `RecInfo` records
  whether the `SA` tag is present and whether it parsed (`sa_present`, `sa_ok`). `parse_sa`
  validates every field (non-empty rname, positive pos, strand exactly `+`/`-`, valid non-empty
  CIGAR, mapq 0–255, NM ≥ 0, entries terminated by `;` with no empty entry). An absent tag is a
  normal single record; an empty, malformed or wrongly typed tag defers the read with reason
  `split`. For an origin-crossing candidate the `SA` CIGAR of each record must equal
  `sa_cigar(other)`, else the read is deferred with the read-level reason `split_sa_inconsistent`.
- **`de` validation.** A mapped record's `de` must be finite and ≥ 0; otherwise the read is
  deferred with reason `de` (covers NaN, ±inf and negative values).
- **Memory: decide split reads from the first record.** Non-origin split reads are decided on
  their first record (SA ≠ 1 entry → `split_not_two_records`, different contig →
  `split_diff_contigs`, cannot reach a boundary → `split_not_origin`, bad SA → `split`); a
  non-candidate primary is deferred immediately and its name remembered in `seen`, so later
  records are dropped. Only true origin candidates are buffered (`pending`), and a
  supplementary-first non-candidate stores **only** its name + reason (`pending_reason`), never
  its records; when the primary arrives its FASTQ entry is written with the stored reason. If a
  buffered candidate's mate turns out to be a primary without `SA`, the read is deferred as
  `split_sa_inconsistent` from that primary and the buffered records are dropped (ordering
  regression tests in `tests/translate_synthetic.rs`). At EOF only names are collected, never the
  groups. `// ponytail:` marks the `seen` name-set ceiling (a few million reads). Every
  `BufWriter` (FASTQ, TSV, JSON) is flushed explicitly with error propagation.
- **Oracle equivalence (user decision after M3).** "Equivalent" = same contig, pos, strand and
  MAPQ, same reference span and query length, and `AS` and `NM` **present on both records** with
  equal values. `AS`/`NM` come from minimap2's own scoring; ours are carried over unchanged, which
  is valid because the reversed (`-`) alignment has the same score; a missing tag on either side is
  a mismatch (`kind_requires_present_as_and_nm`). The oracle test reports identical / equivalent /
  mismatch per scenario and asserts (identical + equivalent) ≥ 99.9 %, plus equivalent == 0 on `+`
  block translations. A canonical left-align-to-fixpoint check was tried and dropped:
  `mm_fix_cigar`-style left-align is a non-confluent, non-idempotent greedy pass, so S3 has several
  distinct left-aligned fixpoints and the residual was reported as mismatch even though every S3
  read has equal `NM` (and `AS`) — i.e. they are placements of the same-scoring alignment. The new
  definition classifies exactly those as equivalent (S3/hifi 10, S3/ont 131), so all scenarios now
  pass. A sanity check asserts every translated record's `AS` equals the `AS` it had in `old.cram`.
  Real PacBio and ONT cases remain 100 %.

## M4

- **Streaming I/O design (speed first).** At t=0 we spawn the remap `minimap2 ... new.fa -` and
  `samtools sort -@ N -m 1G -O CRAM --reference new.fa --write-index -o OUT##idx##OUT.crai -` (both
  read from a pipe; minimap2's index build overlaps the diff). The translate pass writes translated
  records as **uncompressed BAM** through an htslib `Writer` opened on `/proc/self/fd/<sort stdin>`
  (`CompressionLevel::Uncompressed`) and deferred reads as FASTQ into minimap2's stdin; the CRAM
  reader uses a shared htslib `ThreadPool`. A reader thread drains minimap2's SAM stdout (htslib on
  `/proc/self/fd/<fd>`) into the same writer behind a `Mutex`; it asserts minimap2's `@SQ` order
  equals `new.fa`'s `.fai` order. The reader thread **owns minimap2's stdout `ChildStdout`**, so the
  pipe read end closes on every exit; an error is published to a shared slot that the sink checks
  each record, and the pass stops instead of deadlocking. On pass end minimap2's stdin is closed,
  the thread is joined, then the writer and sort's stdin are dropped, then both children are waited
  on (exit codes checked). This took the PacBio FAM20843 pass from M3's ~5.9 s to
  **0.92 s**, and the whole lift to **3.0 s** vs 23.8 s full remap (7.8×).
- **Single failure path + staged output.** `lift` writes to `OUT.tmp.cram` (`+ .tmp.cram.crai`) and
  promotes only after all checks pass, so a failed run never destroys a previous good output. The
  `Children` holder is created before spawning and kills+reaps in `Drop` unless explicitly
  finished, so every `?`/setup error is covered; the shared `fail` helper additionally joins the
  reader thread, removes the staged files, and appends child exit statuses to the error (e.g.
  `minimap2 exited with status 43`). Promotion moves any existing `OUT`/`OUT.crai` to `.bak`, renames
  the index before the CRAM, and restores the backups on any failure. Covered by tests using
  `CRAMHOOK_MINIMAP2` / `CRAMHOOK_SAMTOOLS` wrappers, a wrong-`@SQ` wrapper (must not hang), injected
  exit 43/42, and a promotion failure (target index path is a directory).
- **Remap preset on fallback.** A header-level contract violation becomes an all-defer fallback with
  the real reason in `stats.contract_reason` (`contract_ok=false`); the remap preset is recovered
  from the minimap2 `@PG` when parseable, else from a new `--preset` override, else hard error. The
  same reason/preset handling applies to a mid-stream violation after the restart.
- **Thread split policy.** `--threads N` is split (CRAM decode, remap minimap2, sort, diff
  minimap2) = `(N/2, N, N/2, N/2)` by default. Measured on FAM20843 / FAM22000 at N=16,
  `(8,16,8,8)` beat `(4,16,4,8)` (default at first), `(16,16,16,16)` and `(2,16,2,4)` (which
  starved sort by ~1 s). Overridable for experiments with `CRAMHOOK_THREAD_SPLIT=d,m,s,diff`.
- **Final check without re-reading the output.** Translated primary count (`stats.translated_reads`)
  + remapped primary count (counted by the reader thread) must equal the input primary count, and
  the remapped primary name set must equal the deferred name set (translated names are already
  unique). This avoids a second CRAM scan.
- **`CRAM_OPT_DECODE_MD` stays default (1).** Decoding with `decode_md=0` removes `MD` *and* `NM`
  (verified). `NM` is required by the `SA` rewrite for split reads and by the equivalence check, so
  the reader cannot use `decode_md=0`; the decoded tags of the output are identical to the oracle's
  (both `MD`/`NM` regenerated on decode, `AS`/`ms`/`de`/... carried through). No time is lost:
  the given floors show `decode_md=0` is not faster here (0.61 s vs 0.45 s decode).
- **Temp files.** Normally none: the intermediate BAM is a pipe and the deferred FASTQ goes straight
  to minimap2's stdin. `--keep-tmp` writes `<out>.tmp/deferred.{fq,tsv}` (a tee of the FASTQ) and
  the diff outputs (`<out>.tmp/diff.*`).
- **Stats.** `--stats` JSON now also records per-phase wall times (`diff`, `translate_pass`,
  `writer_drain`, `sort_index`, `minimap2_total`, `total`) and reaped-children CPU time
  (`cutime+cstime` from `/proc/self/stat`). Both child stderrs are captured to a temp file so the
  console stays quiet; that text is surfaced only on child failure, and `minimap2_index` is parsed
  from the minimap2 `loaded/built the index` line (0 if that line is absent or unparseable).
- **FAM14197 agreement (99.98 %).** Raw CIGAR agreement is 98.32 %; with the M3 equivalence
  (same contig/pos/strand/MAPQ/span + equal `AS`/`NM`/`de`) it is 99.98 %, leaving 13 genuine
  mismatches (10 `-` strand, 3 `+` strand; `AS` equal, `NM` ±1). FAM20843 and FAM22000 are 100 %.

## M5

- **`no_block` vs `junction`.** A record whose unclipped window crosses a contig end is
  `no_block` when no identity block covers the (clipped) window, and `junction` only when a block
  covers it but the boundary junction is not preserved. Previously both were reported as
  `junction`/`split_junction`; `check_record` now returns `NoBlock` unconditionally for an
  uncovered window, so the split classifier maps it to `split_no_block` via `split_reason`.
- **Tuning knobs (CLI, defaults unchanged).** `lift`/`translate` accept `--max-de` (overrides the
  preset `de` threshold: 0.01 hifi / 0.03 ont), `--t2-win 10000`, `--t2-step 2000`,
  `--t2-flank 10000`; `diff` accepts the three `--t2-*` flags. `--min-mapq`, `--max-clip`,
  `--margin` already existed; `--max-clip 0` now means unlimited (documented in `--help`). The
  effective values are recorded in `stats.json`/`summary.json` under `params` (`max_de` is the
  resolved threshold).
- **Fail-unsafe clip screen (`clip_hit`, Variant B).** Chimeric reads (chromosome fragment + a
  30–3000 bp tail from the new plasmid p3) were translated without the p3 supplementary record
  because the tail's soft clip is ≤ 200 bp (inside the `--max-clip` limit) so the read looked clean.
  Raising `--max-clip` only makes more of both presets unsafe, so the 200 limit stays; instead every
  record about to be translated is screened. An in-process k-mer screen builds a `HashSet` of
  canonical 2-bit k-mers of `new.fa` that **contain at least one novel base** or **span an internal
  identity-block boundary** (new edit site; window starts `[b−(k−1), b+(k−1))`, contig-end boundaries
  skipped). A soft clip ≥ `--clip-min-len` (20) hits when it shares ≥ `--clip-min-hits` k-mers, and
  the read is deferred with the diff-level `clip_hit` (`split_clip_hit` for origin-crossing reads).
  N-containing k-mers are skipped; canonical encoding makes it strand-invariant. Single records
  screen both free ends; origin-crossing records screen only the **outer** clip (the free end
  reaching the contig boundary) — the junction-facing clip continues into the mate record and is
  screened there. `--clip-flank` (default **0**) adds extra bp around novel intervals / edit sites.
- **Over-deferral finding and Variant B/B′.** The first implementation screened ±1000 bp around every
  novel interval and boundary, including *unchanged* flank sequence. On small plasmids, origin-
  crossing/wrapping reads whose outer clip is simply more of the same plasmid then hit the flank and
  were deferred even though translating them is correct: `ont_FAM22000__cons-to-polished` deferred
  17,812 reads (translated 89.5 % → 73.6 %) with **0** oracle disagreements, `ont_FAM20838` 2,332.
  Variant B counts only novel-containing and boundary-spanning k-mers (flank 0), which removes the
  uninformative flank. Variant B′ (debug `--clip-subtract-old`) additionally drops any k-mer that
  also occurs anywhere in `old.fa`. Table (need.py; `unnec`/`nec` = recovered-without-screen with
  0 / >0 oracle disagreements):

  | case | flank-1000 (old) | B | B′ |
  |---|---|---|---|
  | FAM22000-cons-to-polished | 17,812 (0 nec) | 72 (0 nec) | 63 (0 nec) |
  | FAM20838-cons-to-polished | 2,332 (0 nec) | 6 (0 nec) | 0 (0 nec) |
  | FAM20843-canu | 4 (0 nec) | 4 (0 nec) | 0 (0 nec) |
  | FAM20843-ilesta | 0 | 0 | 0 |

  B is safe on the chimeras (S5c 49/49, S6c 21/21, both models C2b `differ=0`) and on the clip_join
  repro, so it is the default; B′ is also safe on S5c/S6c (`differ=0`) but is kept measurement-only.
  Synthetic translated-count change (default vs `--unsafe-no-clip-screen`) for B is 0 everywhere
  except the chimera scenarios S5c (−61/−80), S6c (−50/−49) and S6b (−24/−17); B′ additionally keeps
  S6b at 0 and reduces S6c to −33/−35. Screen build time is negligible: ≤ 16 ms on the largest
  4.5 Mb case (FAM20843-ilesta), ~6 µs on cons-to-finalized.
- **Chosen k/H (measured).** `k=17`, `H=2`, `--clip-min-len 20`, `--clip-flank 0`. Swept over
  `k∈{11,13,15,17}`, `H∈{1,2,3}`; every setting was a superset of the minimap2 prototype's hit set
  (0 false negatives) with all chimeras caught. k=17 minimises false positives; H=2 is preferred over
  H=3 (identical extras, one more k-mer of margin); any minimap2 exact seed of k≥19 implies ≥3 shared
  17-mers, so H=2 stays a superset.
- **Verified.** S5c/S6c both models with `--xl-tag`: `scripts/compare.py` reports XL=T supplementary
  `differ=0` (S5c/hifi 1225, S5c/ont 1875, S6c/hifi 944, S6c/ont 1518); the same holds for S9/S12.
  The reviewer's edit-site repro (a 4 kb deletion in `B`; `diff.novel` empty) is a regression test
  (`lift_clip_join_edit_site_supplementary_matches_oracle`, `tests/fixtures/clip_join`): before the
  fix 3 XL=T reads differed from the oracle's supplementary, now they are deferred (`clip_hit`) and
  `differ=0`. Alternatives rejected: lowering `--max-clip` (loses legitimate short-clip translations,
  makes HiFi unsafe), and a two-pass minimap2 remap of clips (breaks the single streaming pass and is
  far slower).
- **Bounds and unsafe escape hatch.** `--clip-k` is validated to `1..=32` at the CLI (clap range)
  and in `ClipScreen::build`; `--max-de` must be finite and ≥ 0 (CLI parser + `TranslateConfig::
  validate`); `--t2-win > 0`, `--t2-step > 0`, `--t2-step ≤ --t2-win` are validated at the CLI
  (clap range) and in `DiffParams::validate`/`TranslateConfig::validate` before any work. Disabling
  the screen is an explicit `--unsafe-no-clip-screen` flag (not `--clip-k 0`): it prints a prominent
  warning to stderr and is recorded in `stats.json` `params.unsafe_no_clip_screen`.

## M2 — circular T2 windows (added after M5)

- **All new contigs are treated as circular when building T2 windows.** The old code clipped a
  novel interval's flank-extended range to `[0, len)`, so a novel interval at the contig origin
  produced *no* window spanning the junction (the `ponytail: no circular wrap` comment). A match
  visible only across the origin could then go unprobed and a similar old region stay untainted —
  T2 errors fail unsafe, so this is exactly the direction that must not happen. Now the arc is
  `[iv.start − flank, iv.end + flank)` modulo the contig length, capped at the whole contig
  (covering more only re-probes the same bases); windows may cross the end, and both the probe
  sequence (`circular_seq`) and the novel-base projection (`circ_pos`, i.e. `(new_start + qc) %
  contig_len`) are read modulo. Window length is capped at `contig_len`, so at most one wrap
  occurs and `new_end − new_start` stays exact. Treating a linear contig as circular adds an
  artificial junction, which can only *add* taint, so the probe stays conservative.
  Alternatives: keep clipping (misses the junction — the bug); merge novel intervals that are
  adjacent across the origin before windowing (optional only — overlapping duplicate windows are
  harmless, correctness comes from the modulo projection). Unit tests cover a window crossing the
  end, a novel interval spanning the origin, and a contig shorter than the window.
- **S6w (`scripts/simulate.py`, `tests/diff_synthetic.rs`).** Like S6 (added plasmid p3 with an
  exact IS copy) but the copy spans p3's origin: `p3 = is_seq[half:] + middle + is_seq[:half]`, so
  the circular arc `[L3 − half, half)` is exactly the chromosomal IS. Reads are circular-aware, so
  the copy is only visible as one sequence across the origin. Measured: the pre-fix binary already
  covered **100 %** of every chromosomal IS copy (both presets) because each 600 bp half is
  independently alignable by minimap2 and both halves together cover the whole old copy — so the
  S6-style ≥ 90 % assertion does *not* distinguish fixed from unfixed for an exact split copy. The
  fix is nonetheless required for correctness (a wrapping window is now actually generated and its
  novel bases projected modulo the origin) and for any match whose halves are not separately
  alignable. S6w is kept as the regression fixture for the circular window/projection path.

## M5: MAPQ-instability rule (competitor s2/s1)

- **Rule (added as rule 3b, right after MAPQ).** Defer a record when minimap2's best competing
  chain score `s2` is `>= max_s2_ratio * s1` (`--max-s2-ratio`, default 0.5), or when `s1`/`s2` are
  missing, or `s1 <= 0`. Reason `competitor` (read-level), `split_competitor` for split units.
  Applied per record via `check_record`, so it covers single records and **both** records of an
  origin-crossing split unit (each uses its own `s1`/`s2`); the split classifier maps it through
  `split_reason`. Rationale: on the hold-out set all 6 MAPQ mismatches were reads with a near-equal
  competitor (`s2/s1 >= 0.77`; all other mismatches had `s2/s1 <= 0.32`) whose MAPQ flips when the
  competitor's score shifts slightly in the new assembly (orientation / minimizer-frequency
  effects). Recorded in `stats.json` `params.max_s2_ratio`.
- **Measured (lift --threads 16 --xl-tag, before = baseline3 / val_base, after = this rule).**
  Tuning set `data/real` (24 cases): translated fraction 0.8462 → 0.8377; T field differences
  mapq 0→0, pos 4→4, cigar 22195→22010; T mismatch total 13→13; T supplementary differ 4→4;
  min depth identical 0.9994→0.9994; the rule deferred 34,442 reads. Hold-out `data/val/real`
  (24 cases): translated fraction 0.8351 → 0.8240; mapq 6→**0**, pos 4→4, cigar 360443→357591;
  T mismatch total 29→**23**; T supplementary differ 6→6; min depth identical 0.9832→0.9832; the
  rule deferred 49,232 reads. No new mismatch case appeared; the six hold-out MAPQ mismatches are
  gone.

## M5: minus-indel rule (7b)

- **Rule (7b, after the window/block/taint rule).** When a record would be translated through a
  `-` block and its CIGAR contains an insertion or deletion `>= max_minus_indel`
  (`--max-minus-indel`, default 50; recorded in `params`), defer it. Diff-level reason `minus_indel`
  (`split_minus_indel` for split units); applied per record via `check_record`, so both records of a
  split unit are checked, each against its own block. `+` blocks are unaffected. Rationale: the
  orientation study showed `-` block differences are minimap2's own non-mirror direction effects, and
  the structural ones (position/supplementary differences that hurt depth) all carry a large indel
  in the old CIGAR (the 8 position mismatches had max indels 96 / 963 / 1453 (×3) / 19260 bp).
- **Measured (lift --threads 16 --xl-tag; baseline → s2rule → minus-indel).** Tuning `data/real`
  (24 cases): translated fraction 0.8462 → 0.8377 → 0.8373; T field diffs pos 4 → 4 → **0**,
  cigar 22195 → 22010 → 21817, mapq 0 → 0 → 0; T mismatch total 13 → 13 → **9**; T supplementary
  differ 4 → 4 → **0**; min depth identical 0.9994 → 0.9994 → **1.0000**; `minus_indel` deferred 644
  reads. Hold-out `data/val/real` (24 cases): translated fraction 0.8351 → 0.8240 → 0.8221; pos
  4 → 4 → **0**, cigar 360443 → 357591 → 353685, mapq 6 → **0** → 0; T mismatch total 29 → 23 → **3**;
  T supplementary differ 6 → 6 → **0**; min depth identical 0.9832 → 0.9832 → **1.0000**;
  `minus_indel` deferred 4,551 reads. The four cases below 99.9 % depth (FAM29046-pflye_05 98.32 %,
  FAM20866 hifiasm/myloasm 99.29 %, FAM20866-raven 99.70 %) all reach **100 %**. Remaining
  non-equivalent CIGAR differences (3 per tuning case canu/lja/pbipa; 1 each in FAM14197-minipolish_05
  and FAM20866 hifiasm/myloasm) are the direction effects the rule intentionally does not chase.

## M6 scope (user, 2026-10-09)

- **Goal:** ONT speed on larger datasets. Tiny datasets don't need a speedup: add an automatic
  full-remap fallback below a simple size threshold, chosen from measurements and logged in stats.
- **Use case:** a cleaned assembly → a slightly different cleaned assembly from the same reads;
  circular contigs are often reoriented to dnaA. The general target is "very similar assemblies
  from the same raw reads". Test proxy = raw assembler output → consensus (slightly harder).
- **Data:** `$DATA_ROOT/ont` (11 ONT samples; a separate, incomplete copy of the
  tree is not used). Read sets > 300 MB downsampled to ~100 MB
  with `seqtk sample` (fixed seed). Known bias: at ~100 MB, fixed costs weigh more than at full
  size. Dev on all samples, no hold-out (the user accepts the overfitting risk).
- **Order:** baseline with a deferral breakdown by read bases → stop and report → build only
  the levers the user approves.
- **PacBio:** must stay correct (all gates); speed may change and is reported.
- **Compute reuse:** old/oracle CRAMs are already cached by `build_crams.py`. Full-remap
  benchmark timings are cached per case, keyed by command, tool versions, inputs and threads.
- **Thread split (measured, not changed):** `samtools sort -@ n` was 0.08–0.13 s slower on 3 of
  4 cases (equal on FAM22000) and used up to +870 MB RSS. Giving the diff minimap2 n threads was
  within noise (worse on ONT). `split_threads` stays (n/2, n, n/2, n/2). Data:
  `scratch/worker-threads/out/`.

## M6 refocus (user, 2026-10-09, after the step-1 baseline)

Replaces the "M6 scope" choices above about data and target:
- **Target:** speed on really large datasets only. Small and downsampled cases are no longer
  the optimisation target; the tiny-data fallback is deprioritised.
- **Set:** `data/m6large/real`, i.e. the 6 ONT samples with reads ≥ ~500 MB at **full size**
  (FAM20853, 20856, 20862, 20866, 22000, 29080), **raw primary → polished**, with the same
  3 primaries as `data/m6`; 18 cases. `data/m6` (step-1 baseline) stays as a reference.
- Unchanged: correctness gates, PacBio stays correct, baseline → stop → user picks levers.

## M6 decisions (user, 2026-10-09)

- **s2/s1 limit stays 0.5.** The FAM20856 MAPQ 60→59 read sits at 0.498. Tightening to 0.49
  (+0.09 % deferred) would fit one read; rejected as overfitting. The other three FAM20856
  residuals (incl. one 1 bp pos shift, depth −1 at one position) are minimap2 direction effects,
  reproduced on revcomp(new.fa) (`scratch/m6-mismatch/Report.md`).
- **CRAM output at `-l 1` by default** (manager's call; the user delegated it). It halves the encode
  tail for ~1.3 % larger files; records are identical. Speedups vs full remap are reported with
  the oracle at the samtools default level, and this is noted in the results.
- **Deferral-reduction levers** (small-edit patching, generalised split units) are parked in
  `ideas.md`, not built.

## M7 (2026-10-10)

- **Use-case set `data/v2/real`** built by `populate_m6.py --kind m7` (24 cases, full-size reads).
  Rotation cases: every circular contig rotated by `Random(7).randint(0.1L, 0.9L)` in file order,
  then the largest contig (the chromosome) reverse-complemented, as *dnaA* reorientation does in
  practice. PacBio: the 6 largest read files in `$DATA_ROOT/pacbio` from samples not used in
  `data/real`/`data/val/real` (`data-2025-12-03` reads are not available locally).
- **Mapping reuse across roles:** an existing `oracle` mapping may serve as another case's `old`
  mapping (and vice versa) if reference and reads content, preset and minimap2 version match and
  the CRAM `@SQ` validates. 21 of 48 mappings were reused.
- **FAM20862 gate failures are reported, not fixed** (M7 measures only). The decision on a rule
  (e.g. score-drop ≥ 400) or on the gate definition is the user's.
- **M8 skipped (user, 2026-10-10):** the gate was not met (0.76 % of input bases recoverable,
  ~2.4 % estimated lift-time saving). Small-edit patching stays parked in `ideas.md`.

## M7b (user, 2026-10-10)

- **Rule 7d (`zdrop`, `--max-zdrop` 400).** On `-` blocks, defer a record whose old alignment's
  ksw2 z-drop statistic (path score with the preset's scoring; mismatches from the decoded `MD`;
  `drop = Smax − S − e2·|Δdiag|`, as `ksw_apply_zdrop`) exceeds minimap2's z-drop (400 for both
  presets). Threshold from minimap2, not fitted. Missing/malformed `MD` on a `-` block defers.
  `N` mismatches are scored −b, not −sc_ambi (conservative). Computed event-wise in `to_info`
  (no per-base loop, no allocation): +1–4 % translate time on `-`-heavy cases.
- **Supplementary gate ignores MAPQ < 20 supplementaries** on both sides (`compare.py`
  `SUPP_MIN_MAPQ`); raw differ and ignored counts stay in the output.
- **Accepted as orientation effects** (as in M5/M6): 1–2 bp start shifts and co-optimal primary
  ties on `-` blocks (FAM20862 `202edc37`, `7bd2e7c9`).

## M9 (2026-10-10)

- **Tiny-input fallback (manager's call, from M6/M7 measurements).** A full remap costs ~0.09 s
  per MB of input CRAM (ONT `lr:hq`) and ~0.6 s/MB (PacBio `map-hifi`); lift has ~1.1–1.5 s fixed
  cost, and cases with a full remap < ~2 s gained 0.4–1.1×. Threshold: input CRAM < 20 MB
  (`lr:hq`) / 3 MB (`map-hifi`), i.e. ≈ 1.8 s estimated full remap → skip the diff and remap all
  reads (reason `tiny_input`, class `fallback`, `stats.fallback`). `--min-input-mb` overrides,
  0 disables. Record-level contract checks still run; a contract violation wins. Measured on the
  16.8 MB ONT FAM20847 ilesta case: 1.08 s (full remap 1.42 s, normal lift 3.26 s), output
  identical to the oracle. Synthetic acceptance passes `--min-input-mb 0`.
- **`--threads` default = available cores** (fallback 8); recorded in `params.threads`.
- **`paf` test helper:** computed `qend − qstart` on a deliberately invalid interval; the block
  length is now passed explicitly. Plain debug `cargo test` passes.
- **`pixi run real-smoke`** (`scripts/real_smoke.py`): 5 real cases covering the use case,
  rotation + `-` chromosome (rule 7d), PacBio, split-heavy ONT, PacBio both contigs reversed;
  correctness gates only, ~3 min; skips cases whose local data is missing, fails on missing results.

## M11 (2026-10-11)

- **Translated records are edited in place** (no clone, no `Record::set` on the hot path). On
  `-` blocks CIGAR/SEQ/QUAL are rewritten in the existing buffer when the CIGAR op count is
  unchanged (the usual case), else `set` is used. The `-` invariant compares against a per-pass
  reverse complement of each new contig (≈ one genome of memory). `bin` is not recomputed after
  the position change (neither was it before); harmless because the output is CRAM via
  `samtools sort`, which does not store `bin`. Output records identical to M10 on 4 real cases.
- **No worker pool:** the non-parallelisable part of the pass (decode handoff, `to_info`,
  classification, defer) is ~2.3 of 3.7 s, and at the default `-K` minimap2's drain is the
  critical path; the gain would not justify the risk to output identity.
- **`-K` stays at minimap2's default:** re-tested on the faster pass, best `-K 20M` + 12 remap
  threads gives −7.3 % median (gate 10 %), and tie order becomes timing-dependent.

## M10b (user, 2026-10-11)

- **`s2` removed from translated records** (the only tag below 99.9 % on identical alignments;
  stale competitor score on the old assembly, not recomputable). Removed together with the
  decode `MD` in one aux pass. Rule 3b still reads it from the input. Other tags kept.
- **`XL` on by default**; `--no-xl-tag` disables, `--xl-tag` is accepted (no-op).
- **Duplicate aux tags defer the read (`tags`)**: duplicated tags are invalid SAM and minimap2
  never writes them; the strip also removes every `MD`/`s2` occurrence. 0 such records in all
  local inputs.

## Plan-v3 refinement (user, 2026-10-11)

After reading CRAQ, Inspector, Flagger, Clair3, NucFreq, bcftools (`scratch/manager/v3refs/`):
- **Supplementary alignments counted** in depth, base counts and clips; secondary, duplicate and
  QC-fail skipped; MAPQ per record.
- **`fish` requires `-Y`** (full `SEQ` in every record), like lift; `qc` accepts hard clips.
- **Minority alleles:** beta-binomial (spread fitted per class) and strand balance; strand
  balance alone if the beta-binomial is slow.
- **Site tiers:** yellow = significant minority < 40 %; orange = 40–60 %; red = > 60 % with
  lower confidence bound > 50 % and a minimum read count.
- **Large indels** in the first milestone; `SA` evidence later, only as a flag.
- **`qc` single-pass**; region pass only if measurement shows a need.
- **`qc_qv`** in the JSON, named and documented as a summary of our own calls.
- **Error-spiking generator** in scope (edit a finalized assembly, remap real reads). `qc` flags
  positions; it does not classify structural variants.

## qc step 1 (2026-10-12)

- **Strand gate on yellow minority sites (user):** a significant site with alt < 40 % that fails
  strand balance gets no tier (reason `strand_bias`, still in `candidates.tsv`). Measured on the
  raw → polished sets: ~85 % of the ~29 000 ONT site false positives were strand-biased minority
  SNPs. Sites ≥ 40 % failing strand balance stay yellow with `STRAND`, so clear errors that fail
  the gate are still flagged. `--no-yellow-strand-gate` restores the flag-only rule.
- **Truth by read support (standard eval view):** raw → polished differences are not all majority
  errors (mixed plasmid populations, repeat copies paired differently, low-depth polishing).
  `qc_eval.py` classes each truth event by the highest candidate alt fraction near it
  (`clear` > 0.6, `mixed` 0.4–0.6, `minority`, `no_support`) or `repeat_unaligned` near
  unaligned/multi regions, and reports recall per class.
- **Depth mode excludes zero-depth positions (manager):** a mostly uncovered contig gets the depth
  of its covered part as baseline; JSON has `depth_mode_excludes_zero` and `covered_fraction`.
- **Clusters merge the qualifying 500 bp windows** (as specified), not the site spans; wider
  regions lower the depth ratio, so fewer clusters reach red (52 → 39 over all sets).
- **Indel counting:** a read counts a left-normalised indel only if the normalised anchor lies in
  a depth-contributing segment of that read (support and DP from the same reads); same-read
  pieces normalising to one key are merged; long insertions are counted by exact sequence;
  origin (position 0) indels use the VCF right-padded form with an `ORIGIN` flag.

## qc clusters (2026-10-12)

- **Cluster sites need significance and strand balance by default (manager, from measurement):**
  a minority-SNP cluster site needs, besides ≥ 2 alt reads and AF ≥ 0.10, a nominal beta-binomial
  p ≤ 1e-3 against the SNP background (`--cluster-site-p`) and the site strand gate
  (`--cluster-strand ratio`). Over the 66 raw → polished cases, cluster false positives went
  from 3,260 to 25 on ONT. The 12× FAM20843 cases went from ~637 clusters to 0, and FAM20862
  (~400×) from 148 to 1, so its clusters were strand-biased ONT errors, not a mixed population.
  PacBio repeat/unaligned recall is unchanged (94 %). ONT mixed-class red+orange recall dropped
  36 → 21 %; 53 of the 61 truth events no longer recalled at any tier are in the 12× FAM20843
  cases, where the noise clusters had covered much of the genome and hit truth by chance.
- **Window density test (`--cluster-density-p`) kept but off:** it cut PacBio repeat/unaligned
  red+orange recall from 94 to 51 %.
- `--cluster-site-p 1 --cluster-strand off` restores the step-1 rule.

## qc low coverage, FILTER and JSON (user, 2026-10-12)

- **Low coverage is its own category, not a cap.** A site with DP < `--lowcov-abs` (lr:hq 20,
  map-hifi 10) or < `--lowcov-rel` (0.25) x the contig depth mode is marked INFO `LOWCOV`; the
  absolute case sets `FILTER=LowCov`. The **tier is never capped**: the statistics already
  account for depth. Rationale (user suggestion, manager's recommendation): the tier says what
  the reads show, FILTER says how far to trust it; nothing is rewritten after the fact, and
  `bcftools view -f PASS` gives the confident subset. Example: ~70 sites at DP 11 with 90–100 %
  of reads disagreeing come out `TIER=RED`, `FILTER=LowCov` (real error, few reads), not yellow.
- **Tier moved from FILTER to INFO `TIER`**; FILTER is reserved for confidence (`PASS`/`LowCov`),
  which is VCF convention and works with `bcftools view -f`.
- **JSON alerts count PASS calls**; low-coverage calls are reported in a parallel `lowcov` block
  (and per type). `warnings` carry `low_coverage` per contig (≥ 1000 bp or ≥ 5 %), with the
  reason (abs/rel) and the number of low-coverage red calls. `low_coverage` BED regions are
  annotations (untiered), like `not_assessable`.
- **Clusters** exclude low-coverage sites; **large-indel** BED regions carry `lowcov=1` when
  their mean depth is below the absolute threshold.
- Consequence measured on `real/ont`: about half of that group's `all`-tier recall is in
  low-coverage calls (PASS recall 48.7 % → 19.9 %); every other set is unaffected (no LowCov
  calls). The 12x FAM20843 genome is 92 % low coverage, so its alerts are 0 red / 1 orange
  PASS plus 0 red / 4 orange / 63 yellow lowcov.

## qc insertion grouping (user, 2026-10-12)

- **Similar insertions at the same anchor are counted together** (length ≥ 10 bp, length
  within max(2, 10 %), banded edit distance ≤ 20 % of the representative's length; greedy by
  support). Long ONT insertions carry read errors, so exact-sequence counting split their
  support: a 474 bp insertion in FAM20866 raven_03 went from 476/837 (orange) to 542/837 (red,
  `GROUPED=57`). LINS red recall: m6large 45 → 82 %, val ONT 62 → 88 %, real ONT 0 → 20 %;
  false positives unchanged.
- **Same anchor only (manager, from measurement):** a ±5 bp window gave identical recall but
  pooled alleles from different anchors, so support and depth no longer came from the same reads
  (AF > 1, overlapping REF/ALT). Same-anchor alleles come from distinct reads that all cover the
  anchor. `--ins-group-min-len 0` restores exact-sequence counting.

## qc clip hotspots (step 2b-1, 2026-10-12)

- **Clip hotspot signal (CRAQ idea):** clips ≥ 50 bp (all terminal S/H ops summed per read end;
  edge clips at contig ends excluded as origin crossings or read-end jitter) are chained per side
  within 50 bp; ratio = cluster count ÷ mean depth over 10 bp on the aligned side; ≥ 3 reads;
  red > 0.75, orange ≥ 0.40, yellow ≥ 0.15; a left and a right cluster within 50 bp form one
  `clip_both` region. BED label `clip_hotspot`, `lowcov=1` below the absolute depth threshold.
- Measured on the spiked assemblies (6 samples x 2 replicates): segment swaps 0 → 100 %,
  inversions 8–25 → 100 %, 20 kb deletions 8 → 100 %, 10 kb duplications 50 → 100 % (any tier;
  red 92–100 % for swaps, inversions and 20 kb deletions). Cost: +176 raw -> polished FPs
  (all `clip_hotspot`, of 6,029) and 6 new control FPs (all clip ratio < 0.4, i.e. yellow).
- **Insertion grouping speed (exact):** one high-depth anchor (602 alleles of ~7 kb) took 24 s;
  the banded edit distance is now bit-parallel (Myers with Ukkonen banding, Edlib-style):
  0.42 s, output byte-identical.

## qc depth bands (step 2b-2, 2026-10-12)

- **Depth bands, no HMM (Plan-v3, Flagger idea):** 500 bp windows; `depth_low` (≤ 0.25×, red),
  `depth_half` (0.35–0.65× with high-MAPQ share ≤ 0.25, orange), `depth_high` (≥ 1.75× with share
  ≥ 0.75, orange; red when a trusted `snp_cluster` or `clip_hotspot` lies within 500 bp); regions
  ≥ 1 kb, one gap window allowed; linear contig ends skipped over the mean aligned span; contigs
  whose depth mode is below the low-coverage threshold get no depth signal.
- **Local baseline (manager, from measurement):** the ratio is the window mean ÷ the median of
  window means within ±50 kb (wrapping on circular contigs), not the contig mode. Growing
  bacteria have more genome copies near the replication origin: the unmodified FAM29080 ONT
  assembly runs from ~130× at the terminus to ~377× near the origin, which flagged 337 kb as
  `depth_high`. With the local baseline the window-ratio SD in the 6 controls fell from 0.27 to
  0.078 and control depth regions from 18 to 1. The contig mode still gives copy numbers.
- **`depth_high` threshold stays 1.75×:** 1.4× added control and raw -> polished FPs (0 → 2,
  18 → 25) but covered only one more collapse cluster; collapses in the raw assemblies are short
  (0.5–2.5 kb) and partial (~1.5×), and their SNP clusters already flag them orange.
- Spiked recall: 10 kb duplications red 17 → 58 % (red+orange 100 %); the depth signal adds little
  where clip hotspots already fire.

## qc bridging reads (step 2c-1, user + manager, 2026-10-12)

- **Two layers (user):** detection ("are the reads consistent with the assembly here?") sets the
  tier and does not care about the cause; explanation (a later `likely=` hint) never changes the
  tier and is kept only if it is accurate on the planted errors.
- **Bridging reads, the type-agnostic continuity check:** a read bridges p if one unbroken
  aligned segment (split at I/D/N ≥ 50 bp) covers p ± 500 bp. f = bridged ÷ DP, normalised by
  the median of f within ±50 kb (expected f < 1 because reads ending near p cannot bridge).
  Regions where f_rel ≤ 0.75; red ≤ 0.25, orange ≤ 0.50; ≥ 3 non-bridging reads; contig ends
  (L + 100 bp) skipped, also on circular contigs (origin-crossing reads are split records).
- Why: per-symptom signals split the evidence of one break between them (at a 5 kb deletion
  some reads clip, others carry a long insertion), leaving each signal orange. Counting reads
  that do *not* confirm continuity counts each read once, whatever its symptom.
- Measured: every planted event kind and size is now 100 % red (5 kb deletions were 17 %,
  10 kb duplications 58 %). Control f_rel SD 0.03; control `bridge_break` 0.47 per Mb (mostly
  yellow, concentrated on short PacBio contigs). Raw -> polished: +119 FPs (+15 red), and 68 of
  244 calls lie in unaligned/multi truth regions (real misjoins of the raw assemblies).
- Known gap: breaks at the origin of circular contigs are not checked (listed as open).

## qc scope decisions (user, 2026-10-13)

- **No origin check for circular contigs.** The use case runs dnaapler, which rotates each
  assembly to start at dnaA. An assembler may get the region around the original start wrong,
  but there is no reason to expect errors at the start of dnaA, so the origin is not checked.
- **No "two signals → red" combination rule.** It was meant to turn split evidence into red; the
  bridging-read signal already makes every planted error red, so the rule would mostly add red
  false alarms.

## qc findings and the `likely` hint (step 2c-2, 2026-10-13)

- **Findings:** all tiered calls (VCF sites and tiered BED regions) within 1 kb of each other form
  one finding (`P.findings.bed`, JSON `findings`); tier = the highest member tier. Existing VCF and
  BED are unchanged. On the planted errors one event gives exactly one finding, except inversions
  (2, one per breakpoint) and swaps (4).
- **`likely` hint from a fixed rule table** (collapsed_repeat, false_duplication, unclear on
  conflicting long insertion + long deletion, extra_sequence, missing_sequence, misjoin,
  base_error, mixed, minority). It never changes a tier.
- Measured on the planted errors (12 spiked runs): the hint is right for 100 % of deletions,
  random insertions, duplications, inversions and swaps, 99–100 % of SNPs and small indels, and
  10 of 12 collapses (the other 2 say `unclear`). Accepted answers: a deletion may be
  `missing_sequence` (reads span it) or `misjoin` (≥ 20 kb, reads cannot span); a short collapse
  may be `missing_sequence`, because long reads spanning it show the missing copy.
- `collapsed_repeat` (2× depth) only fires on real raw assemblies, where the repeat is longer than
  the reads; the planted collapses are 1.0–1.4 kb.

## Comparison tools in separate pixi environments (user, 2026-10-13)

- `testing-clair3` (bioconda clair3 2.0.3, models shipped in the package) and `testing-sniffles2`
  (sniffles 2.x), both with `no-default-feature = true` so their dependencies stay out of the
  default environment. Used only for comparisons (`scratch/cmp-tools/`), never by cramhook.

## Read-end pile-ups stay flagged (user, 2026-10-13)

- FAM29030 PacBio scf3 (circular, 45 kb): read starts concentrate near position 2,436 and depth
  before it drops to 0.27×; reads confirm the sequence everywhere. `bridge_break` and `depth_low`
  flag it red. No fix: the pattern is unusual and worth a look, so flagging it is intended even
  though it is not an assembly error.

## Findings: site calls join regions only within 50 bp; cross-mapped calls in the planted eval (2026-10-13)

- **`--finding-site-join` 50 bp:** region-type calls still group within 1 kb, but a site call joins
  a region's finding only within 50 bp of one of its regions. Before, a real small error within 1 kb
  of an unrelated background signal took that signal's hint (2 of 720 planted small errors got
  `misjoin`); now 1 of 720 (a clip hotspot 25 bp away, genuinely ambiguous). Large events keep one
  finding each (two for inversions, four for swaps).
- **Cross-mapped calls counted separately (planted eval only):** after a planted deletion, reads of
  the deleted sequence map to similar sequence elsewhere and create real-looking calls there. The
  scorer attributes a call to this only with spike-specific evidence: no matching call in the
  control at the projected original position, and at least half of the reads supporting the
  specific signal come from the deleted interval. False alarms away from planted errors over the
  12 spiked runs: 316 → 205 after removing 111 cross-mapped calls (9.3 → 6.0 per Mb).

## qc linkage of minority alleles (step 2d, 2026-10-13)

- **Sparse observations, no regions (redesign after review):** the second pass collects, per counted
  record covering >= 1 candidate SNP site, `(read, site, state)`; read names are interned per contig.
  Observations are consolidated per read name across primary and supplementary records (one state per
  read/site, `other` on disagreement; strand from the primary record or the first record), then pair
  counts are accumulated **per read**. Pairs are counted genome-wide per contig (no regions, no caps),
  so every pair within `--link-window` is counted exactly once for any window size. The earlier
  sub-region splitting (<= 250 kb / 25,000 sites) is removed: it split one haplotype at the overlap
  boundary, dropped cross-boundary pairs and produced a site carrying both `LINKED` and `UNLINKED`.
- **Linkage test:** pooled `a >= --link-min-reads` (3), concordance `a/(a+b)` and `a/(a+c) >=
  --link-min-conc` (0.8), phasing `c/(c+d) <= 0.1` and `b/(a+b) <= 0.2`, one-sided Fisher exact p <=
  `--link-p` (1e-6); the forward and reverse tables each need `a >= 2` and Fisher p <=
  `--link-strand-p` (0.01, a conventional level). The strand threshold is looser on purpose: the
  strand check only guards against strand-specific errors, significance comes from the pooled test.
  Sub-threshold sites down to `--link-min-af` (0.05) are included so linkage can rescue a weak site.
- **Groups with a conflict guard:** linked pairs are processed strongest first (Fisher p, then
  concordance, then positions); union-find merges two components only when no *decisive
  incompatibility* joins them. A pair conflicts when `a+b+c >= --link-conflict-min-reads` (6) and
  `a/(a+b+c) < --link-conflict-max-share` (0.3): enough reads carry the alt at X or Y and the two
  minority alleles mostly sit on different reads (two haplotypes). Underpowered but compatible pairs
  (e.g. a=3, b=c=0) are neutral and do not block a merge. Group reads (alt at >= 2 group sites, ref at
  none) go to `P.haplotypes.tsv`.
- **Unlinked classification once per contig, after grouping:** a candidate minority site (AF < 0.40)
  with >= 2 partners each sharing >= 10 reads and in no group gets `UNLINKED` under
  `--unlinked-action flag` or is dropped to tier none under `drop` (the default, user 2026-10-13:
  halves ONT yellow false alarms with no red/orange recall loss; see below).
- **"Belongs elsewhere":** the 1 kb window over the sorted SA targets holding the most distinct group
  reads (ties -> leftmost), reported only when it holds >= 30 % of the group reads (a read counts once
  even with duplicate SA entries); a shared clip side over >= 30 % is reported as `clipped=`.
- **Copy-number evidence:** a group red from its own depth ratio (>= 1.75) and high-MAPQ share
  (>= 0.75) carries `depth_high=1` in its extra, so its finding gets evidence `depth_high` and
  `likely=collapsed_repeat`; an orange group with no depth/structural evidence is `subpopulation`.
- **Tiers:** `linked_haplotype` regions are orange; red when a trusted `depth_high` region overlaps or
  the group's own depth ratio / high-MAPQ share cross the thresholds.
- **Measured (new design, `--link-strand-p` 0.01):** planted errors: recall unchanged (100 % SNP/DEL,
  99 % INS red); spiked `linked_haplotype` regions 36, 26 with `elsewhere=`. Controls: `flag` 59
  findings (3 red / 3 orange / 53 yellow; 2 linked groups: an orange `subpopulation` and a red
  `collapsed_repeat`), `drop` 43 (3/3/37). Raw -> polished (66 cases, s11 vs s10b): red+orange recall
  up (m6large `repeat_unaligned` 81 -> 86 % red, val pacbio 94.5 -> 100 % red+orange); 164 groups,
  4631 linked sites, 102 119 unlinked, 61 of 74 tiered `snp_cluster` regions coincide; yellow FPs
  unchanged under `flag` and roughly halved under `drop` (m6large 1174 -> 540, real ont 1952 -> 1612).
  ONT noise set: group A 188/274 = 68.6 % linked vs B 316/944 = 33.5 %. Timing (16 threads;
  linkage-pass time, qc total, max RSS): FAM20838 minipolish_04 0.32 s / 0.80 s / 521 MB,
  FAM20862 flye_04 1.02 s / 1.33 s / 838 MB, FAM29080 minipolish_01 1.61 s / 1.83 s / 589 MB
  (all below 1.5 GB).
- **Recommendation:** `drop` (halves ONT yellow FPs with no red/orange recall loss); keep `flag`
  while the `fish` consumer still wants the unlinked sites listed.

## Unlinked minority sites are dropped by default (user, 2026-10-13)

- `--unlinked-action drop` is the default: a minority site (AF < 0.40) tested against ≥ 2 partners
  and linked to none gets no tier (reason `unlinked`), stays in `candidates.tsv`, and leaves the
  VCF. Yellow false alarms over the raw -> polished cases: m6large ONT 1,174 -> 540, real ONT
  1,952 -> 1,612, val ONT 274 -> 205, val PacBio 946 -> 847; unmodified controls 59 -> 43 findings.
  No red/orange recall loss; any-tier recall of minority-class truth events drops (m6large ONT
  40 -> 33 %). `flag` keeps them yellow with INFO `UNLINKED`.
- **The trade-off, counted** over the 66 raw -> polished cases: `drop` removes 1,142 false alarms
  (all yellow, 4,346 -> 3,204) and loses 4 truth events (all m6large ONT, all "minority" class:
  < 40 % of reads disagree); 0 planted errors lost. About 285 false alarms per real signal.
- **Cost by construction:** "unlinked" only means "no partner", not "noise". A real minority
  variant that stands alone (a single-SNP subpopulation) among background-noise candidate sites
  is dropped with the noise (reason `unlinked`, still in `candidates.tsv`). The synthetic test
  `qc_planted_events_and_thread_identity` checks this behaviour. Use `--unlinked-action flag`
  to keep such sites.

## qc circularity evidence per contig (user request, 2026-10-14)

- **JSON only** (not VCF): per contig `circularity` with a verdict (`circular` / `linear` /
  `unclear`), the FASTA header flag (`circular=true|false`, or null), reads spanning the
  end→start junction, reads ending at each end and clipped there, and the junction quality
  (`clean`, `overlap_<n>`: the contig's ends repeat each other; `gap_<n>`: reads carry extra
  sequence at the junction). Warnings for header/verdict mismatches, junction overlaps and gaps.
- A read spans the junction if two of its alignments on the contig are consecutive in its
  same-strand alignment chain, the first reaches the contig end (±50 bp), the second starts at
  position 0, in read order. Alignments on other contigs do not break adjacency: a read going
  end → other contig → start is the "sequence missing at the junction, present elsewhere" case
  (e.g. raw FAM29080 flye contig_8 `gap_1855`, whose missing 1.85 kb sits on contig_3).
- Measured on 334 contigs (15 polished, 66 raw, 6 rotated cases): no contig with
  `circular=true` in its header gets another verdict; 3 raw metaMDBG contigs with
  `circular=false` are spanned by reads. 14 junction defects, all in raw assemblies (e.g. canu
  tig00000013 ends overlap by 5,137 bp). All rotated and polished junctions are `clean`. Sniffles2
  agrees on all 15 main replicons; it reports nothing for small high-copy contigs that read
  evidence shows as circular (two checked: they match circular assembler unitigs exactly).
- Cost: one more pass over each contig, ~0.13 s on the largest case (~6 % of qc time).

## `fish` step 1: categories, own per-read checks, no qc read table (user + manager, 2026-10-14)

- **qc writes no per-read table.** A staged qc evidence pass (`P.reads.tsv`) was dropped: fish
  runs the same per-read checks itself (shared helpers in `src/scan.rs`: `terminal_clips`,
  `raw_indels`, `merge_large`), with its own thresholds, so qc stays unchanged (outputs
  byte-identical to 7c3b69d). fish reads qc's `P.findings.bed` and `P.haplotypes.tsv`.
- **One on/off flag per category** (user): `unmapped`, `clip`, `large_ins`/`large_del`, `linked`
  on by default; `--locus` (all reads overlapping a finding ± flank) off. No `no_bridge` category
  (it duplicated clip/indel). Collapsed repeats are not a design driver (autocycler fixes them).
- **Scope:** `findings` (events within a finding of tier >= `--min-tier` ± `--finding-pad` 1 kb;
  default when `--findings` is given) or `genome`. Unmapped reads are always genome-wide.
  A clip at the end of a circular contig is an origin crossing and does not count.
- **Outputs:** `P.fastq.gz` (each read once, original orientation, from the primary record),
  `P.fish.tsv` (one row per read, reason and finding; an event counts for every finding it falls
  in), `P.fish.json`. Two passes: selection, then sequences for the selected reads only. Hard clips
  are an error (`-Y` required); the record count is checked against the index totals.
- **Eval truth classes** (`scripts/fish_eval.py`, user-approved): for planted deletions and
  collapses a control read is *must-find* if its aligned bases cover >= min(100 bp, event length)
  of the removed sequence (union over its records), *neutral* for 1 bp up to that, *wrong* at 0.
  The earlier single 100 bp rule made 50 bp deletions unscorable (93 % of the "15 % wrong").
- **Measured** (tag `f3`; 12 spiked runs, 6 controls; 16 threads), findings scope: recall of
  must-find reads 99.6 % for deletions (19,032 / 19,108; 3,354 of them only as unmapped) and
  99.6 % for collapses; selected reads at an event with no overlap with it 0.3 % / 0.6 %. Genome
  scope: 100.0 % / 99.9 %, but on the unmodified controls it selects 5.5 % of ONT reads (6.8 %
  of bp) vs 0.3 % in the findings scope. Speed: ~1.3–1.5 s at 16 threads, ~400 MB.

## `fish --bait` via lift (user + manager, 2026-10-14)

- **Mechanism (user: "lift-like logic"):** `--bait` writes `P.bait_ref.fa` = assembly + bait
  contigs (name clash = error) and runs the lift session in-process (old = `--ref`, new =
  `P.bait_ref.fa`, `cramhook lift` defaults) into `P.bait.cram`, kept for IGV. Reads with a
  counted record on a bait contig are selected (reason `bait`); detail gives MAPQ and minimap2
  `s1`/`s2`. lift remaps exactly the reads a full remap could place differently (T2 taint of
  regions similar to the bait, clipped/split/unmapped reads), so no separate candidate filter or
  score comparison is needed. Reads whose original primary is duplicate/QC-fail are skipped, as in
  the other categories.
- **MAPQ 0 bait reads are kept:** interior reads of an integrated element (prophage) fit the
  chromosome copy and the bait equally; minimap2 picks one at random, so about half land on the
  bait. For reassembly both copies are the same sequence.
- **Measured (ONT demo, `data/v2` polished → pbfinal, phage missing from the old assembly;
  truth = reads on the phage contig in `oracle.cram`):** FAM20853 (40.9 kb phage) default fish
  selects 0 % of the 732 truth reads (686 map cleanly inside the prophage), `--locus` 13 % at
  17,730 reads, `--bait` 70 % at 918 reads; FAM29080 (20.3 kb) 21 % / 82 % at 12,727 / 74 % at 401.
  The < 100 % bait recall is the random prophage/phage tie, also present in the oracle. myloasm on
  the `--bait` reads recovers each phage completely at 100 % identity, same as on the truth reads.
  Cost: ~7 s and 2.5 GB at 16 threads (memory is lift's samtools sort and minimap2; fish ~275 MB).
- **qc on the demo:** ONT FAM20853 has no flag at the prophage (free phage ~1.1x, ~5 % clipped
  reads at each end); ONT FAM29080 red; PacBio HiFi on the raw autocycler assembly flags both
  (free phage 1.7x / 1.4x). Not tuned (user decides).
- lift now captures minimap2/samtools stderr in private temp files (shown only on failure) and
  parses the minimap2 index time from it (see M4).

## qc split-read junction signal (user + manager, 2026-10-14)

- **Why:** ONT FAM20853 carries a prophage whose excised circle is ~5 % of reads; qc had no flag
  (clip hotspot needs ratio >= 0.15; here 15 clips / ~330 depth = 0.045) and Sniffles2 2.8.1 drops
  it too (DUP VAF 0.048 < its 0.05 mosaic minimum). Reads whose consecutive alignments join two
  distant positions are a specific signal (`scratch/phage-sens/report.md`).
- **Signal (`split_junction`, BED only; user: no VCF record when the BED has it):** consecutive
  alignments of one read on the same contig, strictly increasing in read coordinates, read gap
  <= 200 bp, reference distance >= 1 kb. Canonical key (pos1, side1, pos2, side2); reverse-
  complement traversals merge; display class `tandem` / `deletion` / `inverted_l` / `inverted_r`.
  Junctions within 2 kb of either end of a circular contig are excluded (origin joins belong to
  the circularity block). Clusters within 200 bp; yellow at >= 5 reads and >= 2 % of depth;
  orange/red at the clip thresholds (0.40/0.75). Two regions per cluster (one per breakpoint,
  `partner=`). Counts as clip evidence for `likely` (no new label). Same pass as circularity.
- **Measured:** all four phage cases flagged, one `tandem` cluster each (ONT FAM20853 15 reads /
  4.7 %, ONT FAM29080 65, PacBio FAM20853 160, PacBio FAM29080 68). Planted errors (sp10 vs
  sp7drop): recall identical; away false alarms 309 -> 348; controls 43 -> 45 findings.
  Raw -> polished (s14 vs s12head, 66 cases): recall identical; +86 false alarms, all
  `split_junction` (72 yellow / 10 orange / 4 red), none within 2 kb of a contig end. Background on
  10 clean genomes: 2 clusters. Time unchanged (1.14 s), memory +2 %.

## qc depth steps (user + manager, 2026-10-15)

- **Why:** depth bands compare 500 bp windows with a ±50 kb median against fixed ratios (1.75x
  high). At the PacBio FAM20853 prophage depth jumps 255x -> 400-540x at both ends, but only the
  middle 8.5 kb crossed 1.75x; PacBio FAM29080 (1.33x, on a 2.35x replication slope) was missed.
  Variants tried in a numpy prototype (`scratch/depth-tune/`): hysteresis, lower threshold +
  minimum length, origin-normalised penalised segmentation (CNV-BAC style; needs the origin,
  user rejected), step detection. The user's idea won: a locally standardised step contrast.
- **Signal (`depth_step`, BED only):** Δ = mean depth over the next 100 bp minus the previous
  100 bp (prefix sums on qc's depth array); z = (Δ - median) / max(1.4826 MAD, sqrt(depth/100)),
  median/MAD over ±25 kb on a 1 kb grid (the Poisson floor binds at most positions; without it
  flat stretches give z in the thousands). Local maxima of |z|, suppression ±500 bp, 1 bp
  refinement. Circular contigs wrap; linear ends count. Tiers: yellow |z| >= 8; |z| >= 20 is red
  when f = |Δ| / max(depth sides) > 0.5, else orange (the user's tier rule, applied to this signal
  first). No origin, no background subtraction, no minimum relative step.
- **`depth_segment`:** an up step paired with a down step (no stronger step between) gives an
  interval, tier of the weaker end, ratio inside/flanks. Output only: not a findings member and not
  counted in alerts (a 40 kb interval would merge every call inside it). Depth bands stay.
- **Measured:** prophage ends exact (0 bp) for PacBio FAM20853 (z ±60), PacBio FAM29080 (±30),
  ONT FAM29080 (±29), all orange (f 0.21-0.41), one segment each; ONT FAM20853 z +5/-6 (below
  yellow; the junction signal flags it). Inspection of 40 background peaks: 21 real structure,
  13 technical (7 contig starts with clipped reads, 3 GC, 1 read-end pile-up, 1 homopolymer,
  1 MAPQ), 6 unclear. Planted errors (sp11 vs sp10): recall unchanged; away false alarms
  R/O/Y 11/81/256 -> 17/120/312. Raw -> polished (66 cases): no existing signal changes tier;
  depth_step false alarms 0.30/0.61/0.46 per Mb (R/O/Y). qc 1.22 s on FAM20838.
- **Pairing (revised, user 2026-10-15):** steps sorted by position are grouped into same-sign runs,
  each represented by its strongest step; a gain joins an up run to the following down run (loss
  mirrored), greedy, each run once; segments whose inside/flank ratio contradicts the direction
  are dropped. The first rule ("no stronger opposite step between") produced 1.2 Mb segments across
  the origin (pbipa_04); plain adjacency cut the PacBio FAM20853 prophage at an inner step.
- **Wrapping is automatic (user):** a contig wraps if its circularity verdict is `circular` or
  `unclear` or the FASTA header says `circular=true`; `--linear` still overrides. JSON
  `circular_by`. No contig in the eval data changes (no `linear` verdicts).
- Eval scripts ignore `depth_segment` rows (output only); rescored baselines sp11/s15b accordingly.

## qc significance scores; `bridge_break` removed (user + manager, 2026-10-15)

- **Step 2 of the tier rule (`scratch/manager/qc-tiers-step-spec.md`):** clip_hotspot, large_ins /
  large_del and split_junction rows carry `s` (-log10 Poisson upper tail of the read count vs the
  signal's local background in ±25 kb, excluding the call; λ floor 0.5) and `f` (disagreeing
  reads / depth). Tiers unchanged until step 3. Calibration on the 6 controls (most control calls
  are real structure): high bars clip s >= 100 (target ~0.1 per Mb on contigs >= 500 kb, not
  tuned to single cases); large indels and junctions keep their existing gate.
- **`bridge_break` removed (user):** it counted reads that do not pass ±500 bp, mixing broken reads
  (already counted by clip / junction / large-indel signals) with plain read starts and ends
  (counted by `depth_step`); FAM14197 scf2:8,697, an unclipped read-start pile-up, became a false
  "disagreement". Measured: 0 planted events only it caught; raw -> polished -3 coincidental SNP
  matches, -119 false alarms, -15 bridge-only findings. Planted hint accuracy unchanged; the misjoin
  hint now comes from clip evidence; 55 raw hints changed (accepted). 26 planted breaks lost red
  (clip f ~0.49: each side holds half the reads) -> step 3 uses (left + right) / depth for
  two-sided hotspots.

## Linked findings (user, 2026-10-15)

- **One problematic region = one finding (user).** Two ends of one event are linked: the two
  regions of a `split_junction` cluster and the two ends of a `depth_segment` (explicit link, not
  parsed text). Every tiered call overlapping a link span is absorbed (user's option b); groups
  are joined by union-find; tier = highest member, `likely` from the pooled evidence, no new
  labels. Links only join ends <= `--finding-link-max` (100 kb, shorter way round on circles);
  overlapping spans chain. `--finding-link-max 0` reproduces the old findings byte for byte.
- Wrapped findings are written as two BED rows with the same id (as other wrapped regions); fish
  and the eval scripts count unique ids.
- **Measured:** planted recall identical; findings per inversion 2 -> 1, collapse 1.5 -> 1, swap
  4 -> 3.75; raw -> polished findings 1,588 -> 1,479, none > 100 kb. Absorbed small planted errors
  take the region's hint (SNP hint 100 -> 96 %, INS 99 -> 97 %, DEL 99 -> 98 %); accepted by the
  user: such small errors usually share the cause of the link (e.g. phage sequence and a SNP
  subpopulation).
- Headers: `topology=circular|linear` is parsed next to `circular=`; a conflict counts as unknown
  with a warning.

## One tier rule for every signal (user, 2026-10-15)

- **Rule (user):** red = highly significant and > 50 % of reads disagree; orange = highly
  significant and <= 50 %; yellow = reported but below the high bar. A finding keeps its highest
  member tier; no "two signals -> red". Spec: `scratch/manager/qc-tiers-step-spec.md`.
- **Per signal:** clip_hotspot s >= 100 (`--clip-high-s`), f = clips / depth, for two-sided hotspots
  (left + right) / depth; large indels, split junctions, snp clusters, linked haplotypes, depth
  bands: the existing gate is the high bar, tier by f; depth_step: |z| >= 20. Sites: red = Wilson
  LB > 0.5 and >= 5 alt reads (was AF > 0.60); a strand-balanced majority site that misses red is
  yellow `weak_majority` (was orange), reported regardless of q; mixed 0.40-0.50 stays orange;
  AF 0.50-0.60 without a confident bound is yellow (user). Removed: `--clip-orange`, `--clip-red`,
  `--cluster-red-depth`, the depth_high red upgrade, the linked-haplotype depth_high trigger.
- **Measured (sp23/sp24 vs sp22; s26 vs s23):** no planted event loses red or orange; red gains
  (del 5 kb 2 -> 12 events); 14 of the 26 breaks red only via bridge_break are red again, 12 stay
  orange (large-indel member f 0.31-0.43). Hint accuracy unchanged. Raw -> polished: any-tier
  recall and total false alarms unchanged; 34 truth events orange -> yellow (16 weak_majority,
  18 AF 0.50-0.60), none red -> lower. Phage ends orange (f < 0.5) except PacBio FAM29080 red.
- **depth_segment pairing, third revision:** alternating runs allow two matchings; take the one
  pairing more runs, then the smaller total length. Greedy from the leftmost run paired the origin
  run on FAM20862 / FAM22000 flye contigs into 2.4 Mb "losses" and lost a 2.4 kb 60x gain.

## Depth track shows aligned bases (user, 2026-10-15)

- `P.depth.bedgraph` counts reads with an aligned base (deletions excluded), as IGV's coverage and
  `samtools depth -a` (identical to samtools on real data, max difference 0). qc's own DP (base or
  deletion) and all signals are unchanged; a large deletion is reported by `large_del`, not by a
  depth dip. Reason: the file is for viewing, and a flat track over a 450 bp deletion that IGV shows
  as a dip confused review.

## Gray tier for homopolymer indels (user, 2026-10-15)

- **New tier `gray` below yellow:** a minority or strand-biased yellow INS/DEL whose reference
  homopolymer context (INFO `HP`) is >= `--gray-hp-min` (5) becomes gray. It stays in the VCF
  (`TIER=GRAY`, `HPGRAY`, normal FILTER), is counted separately in the JSON, never creates or joins
  a finding and is excluded from `qc_qv`; the eval scripts treat it as no call. Default on for ONT
  and HiFi; `--no-gray-homopolymer` disables it.
- **Why (user):** qc should not judge homopolymer length noise; a polisher handles it. Gray instead
  of hiding keeps the calls visible and countable. Red/orange HP indels (the consensus itself
  disagrees with most reads) and `weak_majority` HP indels stay as they are; SNPs are untouched.
- **Measured:** raw -> polished (66 cases) 421 ONT and 103 PacBio calls gray; yellow findings
  ONT 878 -> 528, PacBio 142 -> 51; red/orange recall and false alarms unchanged; 29 truth HP
  indels (hp 5-13) were found only by a now-gray minority call. Planted recall unchanged.

## `likely=depth_change` (user, 2026-10-16)

- **Rule (user):** a finding whose sole evidence is a `depth_step` region becomes `depth_change`
  (a change in covered sequence, not a site or structural signal). The new `Evidence.depth_step`
  flag is set when a finding has a `depth_step` member (`depth_segment` is output only, never a
  member, so it never sets it) and is listed in `evidence_list()` after `depth_high`. The branch
  sits after the `only_sites` rules and before the final `unclear`, so a `depth_step` only becomes
  the label when no site call (`site_major`/`site_mixed`/`site_minor`) and no earlier evidence
  (`clips` -> `misjoin`, depth evidence, ...) is present. Every earlier branch is unchanged.
- **Why (user):** after the `bridge_break` removal the depth-step-only findings fell through to
  `unclear`. In `s28` there were 132 `unclear`: 130 were `depth_step`-only structural candidates
  (no site flags) and 2 also held a `site_minor`; 21 of the population came from the `bridge_break`
  removal (`scratch/likely-unclear/report.md`). 0 of the 130 overlap a truth event (the raw->polished
  truth set has no structural kinds), so there is no accuracy to lose, and the user chose a neutral
  `depth_change` over the report's proposed `misjoin`.
- **Measured:** `s29` vs `s28` (66 cases): exactly 130 findings move `unclear` -> `depth_change`
  (red 21, orange 45, yellow 64); the 2 `site_minor`+`depth_step` findings stay `unclear`. The only
  other change is the new `depth_step` evidence flag. `qc.bed`/`qc.vcf` byte-identical; recall and
  false alarms identical. `sp27` vs `sp26`: 33 findings move `unclear` -> `depth_change` (orange 6,
  yellow 27) and all 11 control `unclear` become `depth_change`, leaving the control `by_likely`
  distribution at `depth_change` 11, `misjoin` 4, `minority` 2 plus one each of `base_error`,
  `collapsed_repeat`, `extra_sequence`, `mixed`, `subpopulation`; `per_event.tsv` and the planted
  recall/hint tables are unchanged.

## Cross-contig clip links (user, 2026-10-16)

- **Rule (user):** when the clipped part of a `clip_hotspot`'s reads aligns consistently to one place
  on **another contig**, note the link on both sides and draw it in IGV. Nothing is merged: tiers and
  `likely` are unchanged and finding ids are unchanged. The circularity pass collects cross-contig
  events (consecutive alignments of one read on different contigs, in query-interval order, forward
  progress in read coordinates, read-coordinate gap *or overlap* <= `--xlink-max-gap`, default
  1000 bp); the `qc` pass matches events whose endpoint falls in a `clip_hotspot` region
  (± `--clip-merge`, on the side the clip is on: a left clip joins an aligned segment to its right),
  clusters the partner breakpoints (same contig and side, within 200 bp) and keeps the largest
  cluster when it has >= 5 reads and >= 50 % of the hotspot's clipped reads on that side. These are
  starting values, not tuned. `f = reads / depth at the hotspot side` (not a mean of both ends).
- **Why:** the motivating case (pacbio_FAM29030 control) is a repeat-mediated join: `scf6:1805` has
  31/81 reads with a long left clip and each also aligns to `scf2:30574-35546` with ~810 bp of
  read-coordinate **overlap** between the two alignments. The same-contig `contig_junctions` rule
  (`--junction-max-gap 200`) treats that overlap as unaligned sequence and never sees the pair, and
  the two ends already sit in findings unlinked (`finding_6` scf6, `finding_3` scf2). The clip-driven
  cross-contig rule recovers it without touching `--junction-max-gap`, and the 50 % side gate keeps
  the repeat-overlap noise out (`scratch/interjunction/report.md`; option 3 there).
- **Measured on the gate values before implementation** (6 controls, `--xlink-max-gap 1000`,
  `scratch/interjunction/gate_check.py`): FAM29030's `scf6:1805` and `scf1:1105025` left-clip
  hotspots each match one partner cluster at 100 % of their side (31/31, 37/37); the `scf2:35539`
  right-clip hotspot matches the same two loci but only 37+31 = 68 of its 532 clipped reads (12.8 %),
  so it fails the 50 % gate and only appears as the mirrored note on `finding_3`. The other five
  controls have 0 links. The data support the starting values.
- **Outputs:** `P.bed` `clip_hotspot` extra gains `;clip_to=<contig>:<start>-<end>;clip_to_reads=<n>`
  (0-based half-open partner cluster `[min, max + 1)`, one pair per side, `|`-joined);
  `P.findings.bed` extra gains `;link=finding_<k>@<contig>:<s>-<e>,reads=<n>,f=<f>,partner_depth=<d>`
  on both ends (`k` is `.` when no finding contains the other end; several links are `|`-joined;
  identical notes on one finding are deduplicated, first occurrence first); `P.json` gains `xlinks`.
- **No BEDPE (user, 2026-10-16):** a `P.links.bedpe` file was written for an igv.js `interact` track
  and then dropped after the IGV trial: igv.js draws an inter-contig pair as a single-end block per
  end, not an arc between contigs, so the file added nothing the notes did not. The partner-end depth
  moved from the BEDPE `partner_depth`/`f_partner` columns into the `partner_depth=` field of the
  `link=` note (there is no separate `f_partner`).
- **Coordinates and depth:** the link's `f` is `reads / depth` at the **matched side's own
  count-weighted median** (kept from `clip_candidates` as `Region.clip_left`/`clip_right`), so it
  equals that side's ratio; the note's `partner_depth` is the depth at the **other** end as seen from
  that finding: the raw depth at the partner breakpoint on the hotspot finding, and the denominator of
  `f` (the hotspot-side depth) on the partner finding's mirrored note. It is written without decimals
  when it is a whole count and with one decimal otherwise. A partner
  breakpoint at the contig end (`len`) is anchored on the last base (`len - 1`, as in
  `junction_regions`/`emit_large`), the BEDPE interval stays inside the contig, and the finding is
  looked up with the same clamped coordinate. `finding_covers` wraps only on a **circular** contig,
  so a terminal breakpoint never matches an unrelated origin finding on a linear contig.
- **Measured:** `s31` vs `s29` (66 cases): `qc.vcf` byte-identical; `qc.bed`/`qc.findings.bed`
  identical after removing the `clip_to=`/`link=` fields; recall and false alarms identical. 56
  links across 17 cases (0 elsewhere); `sp29` vs `sp27`: reports identical except the tag, 6 links
  (all FAM29030: control 2, rep1 2, rep2 2). Time/memory delta within run-to-run noise (< 1 % wall,
  < 3 % RSS on a 648 MB ONT case). The reviewer fixes changed only link fields vs `s30`: same 56
  links, 13 `f` values corrected to the matched side's ratio (e.g. `ont_FAM20838__myloasm_04`
  0.0084 -> 0.6548, `ont_FAM29080__metamdbg_03` 21.15 -> 0.94), 9 partner intervals clamped at a
  contig end, 5 links re-attributed to the finding that really contains the partner breakpoint, and
  the `partner_depth`/`f_partner` columns added; `P.json` `xlinks` unchanged.
  See `scratch/manager/qc-depth-change-xlink-spec.md` part B.

## fish eval: "wrong" counts reads outside the matched finding (user, 2026-10-16)

- **Why:** linked findings join a planted event's finding with its far linked end (20-93 kb). fish
  attributes clip/indel reads anywhere in the finding, and the eval counted reads that miss the
  planted interval as wrong: 0.35 % (`f3`) -> 4.87 % (`fs1`), recall unchanged. With
  `--finding-link-max 0` it is 0.30 % again (`scratch/fish-wrong/report.md`). Those reads sit at
  members of the same problem region, so the metric, not fish, was changed (user).
- **Rule:** a finding-attributed read is `wrong` only if its aligned footprint touches neither the
  planted interval nor any row of a matched finding (both rows of a wrapped finding), without
  fish's 1 kb `--finding-pad` (with it the test is tautological). The old count stays as
  `wrong_near_event`. Rescored: `fs1` default 4.87 % -> 0.19 %, `f3` 0.35 % -> 0.22 %; recall
  columns unchanged.

## Depth steps at the origin (user, 2026-10-16)

- **Why:** a read crossing the arbitrarily chosen origin of a circular contig is represented by the
  aligner as a primary plus a supplementary alignment, soft-clipped at both ends, so the number of
  reads landing just inside each end differs (median `end_reads / start_reads` = 1.34 over the 57
  `s28` origin peaks). The wrapped step contrast `Δ` therefore peaks exactly at the origin. The qc
  wrap arithmetic is correct (`delta = depth_right − depth_left` across the origin, e.g. FAM20853
  contig_1 `274.0 − 460.4 = −186.4`); a non-wrapping computation would have no left window at
  `p = 0`. It is an aligner/representation artefact at the origin, not a coverage change
  (`scratch/contig-start-steps/report.md`).
- **Rule (option B, origin-exact only):** on a contig that wraps (`ctx.circular`: verdict `circular`
  or `unclear`, or a `circular=true` header), a `depth_step` whose 100 bp emitted interval touches
  the origin (`start == 0` or `end == len`) is not emitted. The report's definition reproduces the
  `s28` figures exactly on `s31`/`s29`: 114 calls = 57 peaks (R/O/Y 18/19/20 peaks); all 57 peaks
  are within ±50 bp of the origin, and there are 10 more peaks within ±500 bp (67 peaks / 124
  calls). The literal "refined position `== 0 (mod len)`" matches only 49 of the 57 peaks, so the
  emitted-interval test is used. Linear contigs and non-wrapping contigs are unchanged.
- **Suppression:** the peak is dropped *after* non-maximum suppression (not before), so a dropped
  peak cannot unmask a neighbour it suppressed and every other `depth_step` call stays
  byte-identical. `depth_segment`s pair from the remaining peaks (a dropped origin step no longer
  anchors a segment). This is a deliberate trade-off, not two free guarantees: an origin artefact
  can still suppress a genuine step within `--step-nms` (±500 bp). Measured with a throwaway
  drop-*before*-NMS build on the same inputs (`s32db` / `sp30db`): **1** extra `depth_step` call
  appears in the 66 raw→polished cases, **0** in the 18 spike runs. The one call is
  `real/ont_FAM20838__minipolish_04-to-cons` `utg000003c` pos 334 (BED 284-384), yellow, `z=16.78`,
  334 bp from the origin, suppressed by the origin peak at 0; it overlaps no truth event (the contig
  has none) and qc scores it `in_region` (inside an unaligned/multi region), so recall and FP are
  unchanged (`in_region` +1). Dropping before NMS would therefore recover this one artefact-adjacent
  call at the cost of the byte-identical guarantee; not shipped.
- **JSON:** per-contig `depth_step_origin_dropped` (0 for contigs without one).
- **Measured, `s31` -> `s32` (raw -> polished, 66 cases):** 114 origin `depth_step` rows removed
  (R/O/Y 36/38/40; `depth_step_origin_dropped` total 57); no other `depth_step` row changes.
  Findings (grouped by id, ignoring renumbering and link notes): 1037 -> 969 logical findings; 82
  whose only evidence was a `depth_step` disappear; **21** mixed-evidence findings change, of which
  **8 split into 2-3 findings** (e.g. `val ont_FAM29080__metamdbg_03-to-cons` `finding_13`, a
  wrapped finding, becomes 3, and `val ont_FAM14197__ilesta_01-to-cons` `finding_5` sheds its far
  INS into a new `base_error` finding) and the other 13 change span/members only; **35** new
  findings have no old counterpart; 0 changes are link-note-only. No non-`depth_step` member is
  lost (per-contig non-`depth_step` member totals are identical). `depth_segment` rows 22 removed /
  4 added (17 cases). `.vcf` byte-identical. Recall unchanged (every `recall_*` column and
  `exact_total` 757 identical); false alarms down by the removed calls: FP R/O/Y 277/1504/3199 ->
  259/1466/3169 (total 4980 -> 4894), in-region 1296 -> 1268.
- **Measured, `sp29` -> `sp30` (18 runs):** 18 rows removed (O/Y 6/12; 9 peaks = 18 calls, all
  within ±50 bp; 12 peaks / 21 calls within ±500 bp), 1017 -> 999 logical findings with the 18
  `depth_step`-only ones gone, **0** mixed findings changed, control findings 22 -> 16, away FP
  335 -> 323 (old FP/Mb 9.87 -> 9.52). Recall unchanged; `.vcf` byte-identical.

## Code cleanup (user, 2026-10-16)

- **Layout:** `src/qc/mod.rs` (7.0k lines) split into `params`, `tiers`, `sites`, `sig`, `alerts`,
  `likely`, `findings`, `calling`, `regions/{cluster,clip,junction,indel,depth,step,xlink}` and
  `test_util`; `mod.rs` keeps `run_qc` (~320 lines). `circularity.rs` keeps the verdicts; split-read
  junction and cross-contig event extraction moved to `junctions.rs`. One commit per module, pure
  moves (verified by a code-line multiset diff and byte-identical outputs).
- **Shared helpers:** `src/bamio.rs` holds `is_cram`, `open_indexed` (CRAM required-field mask per
  caller, unchanged), `keep`, `op_parts`, `OP_*`, terminal clips (u32 and a u64 variant for the
  circularity pass) and `aux_int`. Not merged because they differ: the three median helpers, and
  `base_bits` (case-insensitive) vs `base_idx` (uppercase only, lowercase is NA in qc).
- **Not done (user):** flattening the `main.rs` CLI -> `QcParams` mapping (it would change `--help`).
  No dead code was found (clippy clean; no `bridge_break` remnants).
- **Safety net:** `scratch/cleanup/compare.py` against qc `s33`, planted `sp31`, fish `fb0`
  (findings + genome scope) and phage bait `pbb0`: identical after every step (all inputs CRAM);
  two runs of the same binary are identical too.

## T2 junction taint (user, 2026-10-17)

- **The bug.** T2 taint is projected from *novel* new positions (not covered exactly once by
  identity blocks). An edit that makes existing sequence more similar to a read **without adding a
  base** therefore leaves no novel interval: a deletion that joins a copy an old contig held
  fragmented (the repro), a contig join, a strand flip or a circular rotation. A read on the
  unchanged duplicate copy is translated with its stale old-assembly MAPQ while a full remap gives
  MAPQ 0. Rule 3b does not help (the old `s2/s1` was 0.40) and the clip screen only looks at soft
  clips. Spec `scratch/manager/lift-junction-taint-spec.md`; two regressions fail before the fix
  (`lift_junction.rs`: `deletion_join` lift `A:30000 MAPQ 60` vs oracle `MAPQ 0`; `contig_join`
  lift `A:30000 MAPQ 60` vs oracle `X:1000 MAPQ 0`).
- **Design.** At every junction where adjacent new positions come from non-adjacent old positions
  (`compute_junctions`, `src/diff.rs`) a `[j − J, j + J)` interval (`J = --t2-junction-flank`,
  default 50; `0` disables) is added to a **T2-only** set, *not* to `diff.novel`, so the clip
  screen is unchanged. Windows are built over novel ∪ junctions through the same `interval_windows`
  path (extended by `FLANK` 10 kb, `WIN`/`STEP` as today); `project_t2_hits` tests membership in
  novel ∪ junctions, novel winning when both (so `taint_t2` is unchanged). Provenance is kept
  apart (`taint_junction`); a record overlapping only junction taint defers as the new reason
  `taint_junction`/`split_taint_junction`, so `stats.json` `deferred_by_reason` reports the count
  directly. `summary.json` gains `junctions` (merged count/bp), `junctions_detected`,
  `junction_rotation_wraps`, `t2_junction_flank` and `taint.t2_junction`; `diff` writes
  `junctions.bed` and `T2J` rows in `taint.bed`. Junctions are found from **block-end events**
  (every new position where a block ends; the successor is a block that starts there), robust to
  nested blocks (a position strictly inside another block is `novel` and rejected first), and include
  the last-block→first-block **seam** of a new contig (all new contigs are circular, M2). Block starts
  and ends are indexed by position, so the pass is `O(B log B)`, and `flank == 0` skips it entirely. A
  junction interval is normalised modulo the contig length (a left/right wrap is split; `2J >= len`
  gives `[0, len)`). The wrap count is `junction_rotation_wraps` (informational).
- **Self-hit exclusion (reviewer follow-up, 2026-10-17).** The circular-header skip was **wrong**:
  minimap2 indexes contigs linearly, so a competitor copy that spans the old origin of a circular
  contig is split in the old index and becomes contiguous after a rotation (repro `origin_competitor`:
  old s1 9976 / s2 4984, oracle MAPQ 0, the skip translated it at MAPQ 60). Wraps are now junctions
  like any other. To control the cost, a junction window's **self-hit** is ignored: an old base
  aligned to new position `p` counts only when it is not the base that `p`'s own identity block maps
  to (same old contig, position and strand). Ignoring the self-hit is sound by the
  **single-identity-block eligibility-window invariant**: rule 7 translates a read only when its whole
  window lies inside one identity block, and an exact-match block maps the read identically in old and
  new, so the self-hit carries no placement ambiguity. A read whose window would cross the junction
  fails rule 7 and defers as `no_block`; the self-hit cannot hide a competitor, because any competitor
  is a *non-self* hit of the same window. Verified end-to-end by `lift_deletion_source_safety` (the
  source-spanning read defers via `no_block`) and the two `origin_competitor` fixtures.
- **Reviewer regressions.** `lift_junction.rs` grows to 10 tests: `origin_competitor` (plus and
  `_minus`), `nested_junction` (P1b: nested blocks `[0,20000)`/`[1000,2000)`/`[20000,40000)` hid the
  `[0,20000)->[20000,40000)` adjacency from a windows(2) scan; lift MAPQ 60 vs oracle 0) and
  `diff_seam_junction_found` (P2a: a new contig that is a 30 kb arc of a 40 kb old contig; the seam
  at the new origin was never examined). Unit tests: `compute_junctions_wrap_nested_and_seam`,
  `push_junction_normalizes` (left underflow clamped to 0 instead of wrapping to `[960,1000)`; a
  large `J` overflowed the contig), `is_self_hit_isolates_contig_pos_and_strand` (same contig/pos but
  the opposite strand is not a self-hit; ±1 bp shifted hits are not), `junction_self_hit_is_not_tainted`,
  and `compute_junctions_scales_linearly` (80 000 blocks with every boundary a junction: `O(B log B)`,
  <0.1 s release, 0.5 s debug). Reviewer perf assets, `--t2-junction-flank 0 --t2-flank 0`: 2 Mb
  total 0.84 → **0.33 s**, 8 Mb 9.40 → **1.24 s** (the remainder is the asm5 diff, not junction
  detection).
- **Measured (scripts/real_eval.py, 16 threads, `--xl-tag`; `m7b2` = pre-fix, `jt1` = junction
  taint, `jt2` = circular skip, `jt3` = this fix).**
  - `data/real` (24): translated 0.8373 → 0.8368; T mismatch total 9 → 9; field diffs
    contig/pos/strand/mapq 0; cigar 21 818 → 21 384; supp differ 0 → 0; depth 1.0000;
    junction-taint reads 611 (jt1) → 883 (jt3).
  - `data/val/real` (24): translated 0.8221 → **0.8069** (jt1 0.8070); T mismatch 3 → 3; cigar
    ~353 674 → 352 624; supp 0; depth 1.0000; junction-taint reads 35 610 (jt1) → 36 069 (jt3).
  - `data/v2/real` (22 common): translated 0.8947 → **0.8833** (jt1 0.8850, jt2 0.8874); T mismatch
    18 → 18; pos 2, mapq 1; cigar ~662 915 → 644 829; supp 0; depth 1.0000; junction-taint reads
    36 089 (jt1) → 42 356 (jt3). The two cases with a stale `new.fa.fai` (data drift) are excluded.
  - **No correctness regression:** T mismatch totals, field-diff patterns, supplementary differ and
    min depth are unchanged from `jt1`; the extra deferrals are remapped reads.
  - **Cost (accepted for correctness).** The self-hit exclusion removes only the exact source hit (12
    of the 95 junction-taint intervals on FAM29046 scf1); on real, repeat-rich assemblies a 10 kb
    junction window has many **non-self** hits (short repeat/secondary hits under `-N 1e6 -p 1e-6`)
    that cover the `±J` junction positions, and at ~500x they defer many reads. Removing the unsound
    circular skip therefore restores correctness at a real cost on the wrap-rich sets: hold-out
    `data/val/real` translated **0.822 → 0.807** and `data/v2/real` **0.895 → 0.883** (junction-taint
    reads `pacbio_FAM29046` 18 321 → 18 308, `pacbio_FAM29030` 5 512 → 5 498,
    `ont_FAM20856__polished-to-rotated` 7 652 → 7 637). **Accepted** for correctness; T mismatch
    totals, field diffs, supplementary differ and min depth are unchanged.
  - **A fixed spanning-length gate is not sound.** Counterexample (`scratch/lift-junction/rereview/
    gate-example`, reviewer): a 2 kb read whose competitor is split into 800/600/600 bp pieces, old
    `s2/s1 = 789/1977 = 0.40` (below `--max-s2-ratio` so the old copy is not deferred); after the joins
    a full remap gives MAPQ 0. A gate that deferred only when a competitor chain covers, say, >= L bp
    would drop the 800 bp piece and translate the read at MAPQ 60. A sound reduction needs a
    *read-specific* bound on the omitted chain scores (the sum of the pieces a join could splice),
    which is a possible future item.
- **Timing (16 threads, machine idle; `flank 0` → `flank 50`).** Largest cases: ont_FAM22000
  minipolish_05 total 10.87 → 11.13 s, T2 0.78 → 0.85 s; ont_FAM20866 polished-to-rotated
  9.83 → 9.81 s, T2 0.00 → 0.14 s; ont_FAM20866 myloasm_03 9.72 → 9.79 s, T2 0.07 → 0.19 s. Total
  lift wall is within run-to-run noise; RSS unchanged. Aggregate T2 over each set: real 7.1 → 9.2 s,
  val 14.4 → 22.1 s, v2 1.0 → 2.4 s (T2 is a small absolute part of lift). The added cost is ~10
  windows per junction; merging the junction intervals keeps the window count below the naive
  `10 × junctions`.
- **Scenario regressions (branch `lift-scenarios`, 2026-10-17).** S13-S21, S13b/S14b/S17b and the
  random SR1-SR5 (3-5 random edits) were added by `lift-scenarios` with an xfail list for the
  junction class (S18 hifi/ont, S20 ont, S21 hifi). Merged and the xfail list + its ignored tracker
  removed: `translate_synthetic` (`--include-ignored`) now has **every** scenario × model at
  1.00000 agreement / 0 mismatches, including the four former xfails (S18 hifi/ont, S20 ont, S21
  hifi) and SR1-SR5 (all 1.00000, 0 mismatch). `scripts/acceptance.py` on S13-S21 + SR1-SR5 (34
  runs): C1/C2/C2b/C3 pass, 0 gate failures.
- **Scenario truth fix (P2c).** `scripts/simulate.py` derived the S18/S21 deletion offsets and the
  repeat truth from hard-coded positions (258000/282000) and the un-shifted base repeats, and S20
  put the old `RcopyA` on the new contig `chr_B`. Now the truth comes from the same offsets as the
  construction (`a = Rlen//3`, `b = 2*Rlen//3`; separators at `P+a` = 256666 and `P+b+sep` = 283333),
  `_insert_shift_repeats` shifts chr repeats by `Rlen + 2*sep_len` for the fragmentation/join, S21
  records its polishing edits, and `_R_pieces_records` keeps the old copy on old `chr`. S21's
  polishing edits now use the established `old_contig`/`old_start` interval schema (they previously
  used `contig`/`pos`, which `acceptance.py`'s `edit_intervals` silently skipped, and their positions
  were post-translocation intermediates); `_back_project_blocks` maps each position back through the
  translocation and deletion blocks to original old `chr`, and the translocation's own old interval is
  back-projected too. `acceptance.py` now checks that every planted edit is recognised
  (`unrecognised_edits`, non-zero -> exit 1); S1-S21 and SR1-SR5 all pass. Regenerated with
  `pixi run sim --only S18,S20,S21`: `old.fa`, `new.fa` and both `reads.fastq.gz` are
  **byte-identical** (checked by md5 against the pre-fix generator); only `truth.json` changed
  (S18 rRNA3 old 270000 → 330000, new 290000; S21 polishing edits back-projected). Acceptance C6 for
  the affected runs (before → after the truth fix): S18 hifi 73/73 → 72/73, S18 ont 28/28 → 28/28,
  S21 hifi 370/486 → 312/486, S21 ont 356/450 → 322/450 unnecessary diff deferrals; S20 0/10 and
  0/11 both. The residual C6 "unnecessary" are junction deferrals on the intact copy, whose window
  does not overlap the planted edit (C6's definition cannot see a competitor created elsewhere).
- **Perf (reviewer follow-up).** `compute_junctions` was `O(B²)` (every block end scanned all
  blocks); it now indexes block starts/ends by position, rejects `novel` candidates before any scan,
  and returns immediately on `flank == 0`. Reviewer assets with `--t2-junction-flank 0 --t2-flank 0`:
  2 Mb total 0.84 → 0.33 s, 8 Mb 9.40 → 1.24 s. `compute_junctions_scales_linearly` (80k blocks)
  guards it.
- **Verified.** `cargo test` and `cargo test --release` green (276 lib; **10** lift_junction, 19
  lift_synthetic, 11 diff_synthetic, 7 translate_synthetic + 1 ignored random; `--include-ignored`
  also runs SR1-SR5). Safety net with a fresh tag `jt4`: qc `s33` 462 files, spike `sp31` 126, fish
  `fb0` 108, bait `pbb0` 10 — **all byte-identical**.

## lift accepts lifted input (user, 2026-10-16)

User report: `cramhook fish --bait` (and calling `cramhook lift` twice) on a CRAM produced by
`cramhook lift` falls back to remapping every read. On `scratch/wfcmp/FAM20853/run1/B` the bait
`lift` session had `contract_ok=false`, `translated_reads=0`, `remapped_reads=118936`, and the
bait step took `19.99 s` (vs ~7 s for the same session on a minimap2 CRAM).

**Phase 1 findings (real case `run1/B/lift.cram`, pysam + `samtools view -H`).**

- **What lift reads.** Translation reads, per record: MAPQ, flags, CIGAR, POS, contig, SEQ length;
  and the aux tags `de` (float, divergence gate), `NM` (SA rewrite), `s1`/`s2` (rule 3b competitor
  ratio), `SA` (split classification) and the decode-generated `MD` (rule 7d z-drop). The
  unexpected-tag allow-list was `NM ms AS nn tp cm s1 s2 de rl SA`; a record with any other tag is
  deferred with reason `tags`.
- **Stored header** has four `@PG`: `minimap2` (valid `lr:hq` CL), `samtools sort`, `cramhook`
  (`PN:cramhook`, `CL` is the lifted argv, `DS:old-ref=… new-ref=… remap=minimap2 -ax lr:hq …
  --secondary=no -Y new.fa -`) and a second `samtools sort`. The old contract rejected the
  `cramhook` line: `unexpected @PG line (PN:"cramhook", CL starts with …)`. (A `samtools view -H`
  `@PG` is only added when viewing, is not stored, and was and stays rejected.)
- **Tags on lifted records.** All records carry `XL` (`A:T` translated / `A:R` remapped), `AS`,
  `MD` (decode), `NM`, `cm`, `de`, `ms`, `nn`, `rl`, `s1`, `tp` (and `SA` / minimap2's `zd` on
  supplementary). The **translated** (`XL:A:T`) records have **no `s2`** (M10b strips the stale
  old-assembly competitor score); the **remapped** (`XL:A:R`) records keep minimap2's fresh `s2`.
- **Two obstacles beyond the `@PG` chain.**
  1. `XL` was not in the allow-list, so every record would be deferred `tags`.
  2. `s2` is absent on translated records, so rule 3b (`competitor` when `s1`/`s2` are missing)
     would defer them. The chain could therefore never translate a translated record.
  The `de`/`s1` on translated records are the old-assembly values; translation preserves the
  alignment, so they still describe the translated record (the alignment is unchanged), but a
  stripped `s2` cannot be re-derived — that is exactly why rule 3b needs a decision for a chain.

**Decisions (fail-safe allowlist, manager's redesign 2026-10-17).**

- **Contract accepts a chain only when it is all-default.** Exactly one minimap2 `@PG` (first);
  the remaining `@PG`s must match `S* (L S+)*` (`S` = `samtools sort`, `L` = `cramhook lift`), each
  `L` carrying a valid `lift` `CL` **and** a `DS` `remap=` minimap2 command that parses under the
  existing rules and uses the same preset. `Contract::chained` is set when at least one `L` is
  present.
- **Option allowlist for a chain.** A `cramhook lift` `@PG`'s `CL` is shell-split and every token
  must be a neutral option taking one value or none, no option may repeat, and the four required
  paths must be present. Any token outside the allowlist (a gated option, an unknown option, an
  unexpected argument) or a missing/duplicated value fails the contract, so the input is fully
  remapped. The current lift must additionally use default eligibility/safety options
  (`TranslateConfig::eligibility_is_default`, checked in `lift`/`translate`), otherwise the contract
  fails the same way.

  | `lift` option | class | why |
  |---|---|---|
  | `--in`, `--old-ref`, `--new-ref`, `--out` | neutral | paths only |
  | `--threads` | neutral | resources |
  | `--stats` | neutral | output only |
  | `--min-input-mb` | neutral | perf fallback; its worst case is remapping *more* |
  | `--preset` | neutral | fallback for an invalid input only; ignored on a valid chain |
  | `--strict` | neutral | error handling; a successful chain behaves the same |
  | `--keep-tmp` | neutral | debug temp files |
  | `--xl-tag`, `--no-xl-tag` | neutral | output tagging only |
  | `--min-mapq`, `--max-clip`, `--margin`, `--max-de` | gated | read eligibility |
  | `--t2-win`, `--t2-step`, `--t2-flank`, `--t2-junction-flank` | gated | read eligibility (T2 taint) |
  | `--max-s2-ratio`, `--max-minus-indel`, `--max-zdrop` | gated | read eligibility (rules 3b/7b/7d) |
  | `--clip-k`, `--clip-min-hits`, `--clip-min-len`, `--clip-flank`, `--clip-subtract-old` | gated | clip-screen safety |
  | `--unsafe-no-clip-screen` | gated | disables the clip screen |

- **Rule 3b trusts an inherited `XL:A:T` only for an all-default chain.** An input record carrying
  `XL:A:T` had its `s2` stripped by an earlier lift and already passed rule 3b there; since every
  lift in the chain (and the current lift) uses the default `--max-s2-ratio`, that is evidence it
  passes here. `TranslateConfig::inherited_xl_ok` is set to `Contract::chained` after that check, and
  a missing `s2` on an `XL:A:T` record is trusted only when it is set (with `s1` present and
  positive); a present close competitor still defers. `XL:A:R` records keep fresh `s2` and are
  checked normally. Only `s2` needs this: `de`, MAPQ, clips, `s1`, `SA` and the `MD`-based z-drop are
  re-evaluated against the new revision, so they cannot be inherited.
- **Composition argument.** With every lift all-default, each lift screened its own edits with the
  same default clip screen and eligibility thresholds, and `s2` was checked at the default
  `--max-s2-ratio` (0.5). An inherited translated record is therefore as safe as one translated by
  this lift; the only unrecoverable loss is the stripped `s2`, which was checked at the same
  threshold.
- **Unambiguous argv recording and strict splitting.** The `@PG` `CL` shell-quotes each argv element
  (only when needed, so ordinary command lines are unchanged; `translate::shell_quote`/`shell_join`:
  plain words, else single-quoted with each `'` as the `'\''` idiom). `contract::shell_split` accepts
  exactly those forms: any double quote, any backslash outside the `'\''` idiom, or an unbalanced
  quote returns `None` → the chain is rejected and the input remapped. A path containing
  ` --max-s2-ratio 0.25` is one token and cannot smuggle an option; an older unquoted `CL` with such
  a path yields extra tokens and fails the allowlist.
- **argv-encoding marker.** Even the strict splitter can be fooled by a *legacy* `.join(" ")` `CL`
  where apostrophes in two paths pair up and swallow a gated option between them
  (`--out legacy-out' --max-s2-ratio 1 --stats 'legacy-stats.json` reads as one `--out` value). The
  fix does not try to parse such CLs: every `cramhook lift` `@PG` written by this build carries the
  `DS` token `argv=sh1` (`translate::ARGV_ENCODING_MARKER`), and the contract accepts a chained lift
  `@PG` **only** if it carries it. An unmarked lift `@PG` fails with "lift @PG from an older cramhook
  (no argv encoding marker)" and the input is fully remapped. Consequence: a CRAM lifted by a build
  before this change (e.g. `scratch/wfcmp/FAM20853/run1/B/lift.cram`) is not chained; the workflow
  must be rerun.
- **Neutral-option values are validated** with the CLI grammar: `--threads` a non-negative integer,
  `--min-input-mb` finite and ≥ 0, `--preset` `map-hifi`/`lr:hq`, paths non-empty and not starting
  with `-`; a value token starting with `--` means the value is missing; duplicates are rejected.
- **Unique `@PG` id per lift.** `build_output_header` picks the first unused id (`cramhook`,
  `cramhook.1`, …) and sets `PP` to the input's last `@PG` id; the contract matches cramhook `@PG`s
  on `PN`.
- **`XL` is accepted on input and never inherited.** `XL` is in the allow-list; `strip_translated_tags`
  removes an inherited `XL` unconditionally and `add_xl` pushes the fresh `XL:A:T` only when tagging
  is on, so `--no-xl-tag` leaves no inherited `XL:A:R`.
- **M10b unchanged.** `s2` stays stripped from translated records.

**P1 investigation (2026-10-17): found, then fixed on master by the T2 junction taint.** A competitor
improved by **deletions only** was not caught by an all-default chain, but it was not a chaining gap.
The reviewer's repro (`/tmp/opencode/lift-chain-rereview/default_comp.py`): one read maps perfectly to
contig A and to a fragmented copy on B; r0→r1 and r1→r2 each delete a separator, so B's competitor
becomes contiguous. `map0` (r0): `s1=9977 s2=3982` (ratio 0.40); `map1` (r1): `s2=5989` (0.60);
`map2` (r2): `s2=9977` (1.0, full-remap MAPQ 0). Before the fix the chain translated lift 1 (rule 3b
passes at 0.40, `s2` stripped) and lift 2 kept MAPQ 60, and a **direct lift from `map0` straight to
r2** (both deletions in one diff) did too — same input `s2=3982`, ratio 0.40 — confirming a single-lift
gap, not chaining. The reason: **T2 taint** is projected from *novel* intervals (sequence in new but
not old); a deletion adds no base, so the read's window (on the unchanged A) was never tainted, and
the **clip screen** only inspects soft clips (this read is fully aligned `10000M`; a deletion makes no
novel k-mer either). Rule 3b was the only guard and read the old-assembly `s2`, where the competitor
was still fragmented. Master's "T2 junction taint (user, 2026-10-17)" projects taint from the new
adjacencies too, so with the rebase the repro changes: lift 1 now defers the read as
`taint_junction` → remapped, lift 2 remaps → **MAPQ 0**, equal to the full remap; the direct A→C also
defers `taint_junction` → MAPQ 0. Chain and direct still agree exactly (S1 A→B→C vs A→C: translated
1011 / deferred 229, identical sets).

**Files and tests.** `src/contract.rs` (`shell_split`, the option allowlist,
`validate_previous_lift_cl`/`validate_neutral_value`, `check_cramhook_lift`, `Contract::chained`),
`src/translate.rs` (`shell_quote`/`shell_join`/`recorded_argv`, `eligibility_is_default`,
`inherited_xl_ok`, rule 3b, `strip_translated_tags` removing `XL`, unique `@PG` id), `src/lift.rs`,
`src/main.rs`, `tests/lift_synthetic.rs`, `tests/lift_roundtrip.rs`, `README.md` and this section.
Unit: `contract.rs` accepts a
chain (twice) and rejects a preset mismatch, a non-`lift` subcommand, a missing/bad `remap=`, a lift
without a sort, a foreign program, a non-first minimap2, every gated option, unknown/malformed CLs,
bad neutral-option values (`--threads garbage`, `--threads --unsafe-no-clip-screen`, `--stats --stats`,
`--out=`, `--preset map-ont`, paths starting `-`, …), an unquoted spaced path, double quotes, stray
backslashes and unbalanced quotes, a legacy lift `@PG` without the `argv=sh1` marker (whose
`.join(" ")` CL hides a gated option behind paired apostrophes), plus a `shell_quote`/`shell_split`
round trip; `translate.rs`
covers the gated `XL:A:T` rule-3b path, `eligibility_is_default` for every option (incl.
`--t2-junction-flank`), `XL` in `scan_aux`, and `strip`+`add_xl` for a single and a split record
(tagging on/off). Integration (`tests/lift_synthetic.rs`, all through the real CLI):
`lift_chain_three_lifts_matches_full_remap`
(A→B→C→C on `S1/hifi` with a synthetic revision C; contract, the `argv=sh1` marker, unique `@PG` ids,
and **both** the second and third lifts compared to a direct full remap onto C by primary
contig/pos/strand/cigar/mapq and supplementary sets with MAPQ ≥ 20);
`lift_chain_requires_all_default_options` (relaxed-lift-1 and
non-default-current chains must fail the contract and remap, default→default translates);
`lift_chain_rejects_legacy_lift_pg_without_marker` (a marker-stripped first-lift output fails the
contract and remaps to zero translated);
`lift_chain_clip_join_default_equals_full_remap` (clip_join: default→default equals the direct full
remap including the junction45/50/55 supplementaries; relaxed-screen → contract fails and remaps);
`lift_chain_clip_join_old_to_old_then_new_defers_clip_hits` (old→old translates the three junction
reads; old→new defers them `clip_hit` and restores their supplementaries equal to the full remap).
`tests/lift_roundtrip.rs` runs A→B→A and A→B→A→B over S1–S21 × hifi/ont (random SR1–SR5 in the
ignored `roundtrip_random_scenarios`), each step compared to its full remap with the
`translate_synthetic` equivalence criteria. `pixi run cargo test --release` and `pixi run cargo test`
pass (297 lib + 11 diff + 30 fish + 10 lift_junction + 24 lift_synthetic + 1 lift_roundtrip + 17
circularity + 4 linkage + 26 qc + 7 translate; `--include-ignored` adds the two random runners).

**Round-trip results (S1–S21 × hifi/ont, CLI, release).** Every A→B, B→A and A→B→A→B step is
**100 % identical+equivalent** and **0 mismatch** on primary and supplementary signatures against its
full-remap oracle (`oracle.cram` for →B; the original `old.cram` for →A — verified once that a fresh
`minimap2 -ax … -Y old.fa reads` reproduces `old.cram` 100 % identical). A→B→A reproduces the original
A mapping exactly, and A→B→A→B reproduces the B oracle including the inherited-`XL`/stripped-`s2`
step 3. Translated fraction (drift) is flat across the chain for 20 of 26 scenarios; the six where
the B→A step moves (`S6b`, `S10a`, `S10b`, `S16`, `S20`, `S21`) do so only for the middle B→A step
(the →A diff is the reverse edit class) and recover exactly on step 3. The
lowest fractions are S11 (0.18 hifi / 0.37 ont) and S21 (0.49 / 0.59), where most reads are
correctly deferred by the diff/s2 rules. Equivalent (CIGAR-different, placement-equal) pairs, from
the reverse-complement scenarios: S3 (10 hifi / 145 ont), S13 (4 / 55), S19 (69 / 732); NM drifted on
14 and `de` on 10 equivalent pairs (metrics, not a gate). Runtime: 50 s release / 89 s debug (one
thread per model); the random SR1–SR5 round-trip is ~10 s release.

**Chain measurements (S1/hifi, CLI lifts).** Second lift's translated primaries: default→default
**1231**; relaxed (`--max-s2-ratio 1`)→default **0** (`contract ok: false`); default→`--max-s2-ratio 0.4`
**0** (`contract ok: false`). The three-lift A→B→C→C chain: B→C translated **1011**, C→C **1231**;
both match the direct full remap onto C exactly; `@PG` ids `cramhook`, `cramhook.1`, `cramhook.2`.
clip_join: default→default `contract ok: true` (translated 0 / remapped 10) and equals the oracle
incl. junction45/50/55 supplementaries; first lift `--clip-min-hits 100000` → second lift
`contract ok: false` (`earlier lift used --clip-min-hits …`) and remaps.

**Measured on `scratch/wfcmp/FAM20853/run1/B` (ONT `lr:hq`, 16 threads).**

- `cramhook fish --bait` on `run1/B/lift.cram`: before (`run1/B`, old binary) `contract_ok=false`,
  0 translated / 118936 remapped, bait step `19.99 s`; after `contract_ok=true`, **111565
  translated / 7371 remapped**, bait step **5.47 s** (`6.80 s` wall). The selected read set is
  identical: 912 distinct reads, same `reason`+`class` per read (755 `bait`, 104 `clip`, 48
  `unmapped`, 35 `linked`, 14 `large_ins`, 6 `large_del`).
- `scripts/compare.py` (the lift equivalence tool) on that output vs the full-remap oracle
  `run1/A/bait.cram` (reads mapped to `new.fa`+phage): primary `n=118936`, **100.0000 %
  identical, 0 mismatch** (translated and remapped blocks both 100 % on contig/pos/strand/cigar/
  mapq), supplementary 100 % equal, depth 100 % identical (`max_abs_diff=0`).
- `cramhook lift` `lift.cram` -> `new.fa`+phage directly (`--threads 16`): `contract_ok=true`,
  111565 translated / 7371 remapped, `5.43 s` total; `compare.py` vs the same oracle again 100 %
  identical, 0 mismatch.
- Cleanup safety net (fresh tag `lc2`, `scratch/cleanup/compare.py`): qc `s33` vs `lc2` **462/462
  files identical**, spike `sp31` **126/126**, fish `fb0` **108/108**, phage bait `pbb0` **10/10**
  — no output that does not involve a lifted input changed.
- Note: the lifted CRAM's supplementary records carry minimap2's `zd` tag, which is not in the
  allow-list; such records always take the split-classification path, so `zd` never reaches the
  unexpected-tag check and does not need handling (verified on `run1/B`).

## `fish --bait-records` (user, 2026-10-17)

- **Why (user):** a `bait` read can reach the bait *only* through a supplementary alignment (its
  primary record lies elsewhere). In the FAM20853 case study these are 9 chimeric reads: they are
  the entire difference between B' and the standard workflow A, and they move myloasm off the clean
  **circular 40,865 bp** contig (A) onto a 62.7 kb non-circular one. Both selections are wanted, so
  both are kept (user decision).
- **Flag:** `--bait-records any|primary`, default `any` = today's behaviour (a counted record on a
  bait contig, primary or supplementary). `primary` selects a read for `bait` only when its
  **primary** record is on a bait contig; supplementary-only reads are dropped. The reasons-TSV
  `detail` already carries `record=primary|supplementary`, so the TSV needed no change.
- **JSON:** the `bait` block of `P.fish.json` always carries `mode` (`any`/`primary`) and
  `supplementary_only_reads` — the number of reads whose only counted bait record is supplementary,
  reported in both modes (a diagnostic). An absent flag and an explicit `--bait-records any` are
  equivalent (same selection and same JSON apart from timing). Against the previous baseline
  `pbb0`, the only change is these two keys in the `bait` block (no selection change); the new
  reference tag is `pbb1` (see `scratch/cleanup/baseline.md`). The no-bait fish tree `fb0` stays
  **108/108** identical.
- **Measured (ONT FAM20853 `run3/B/lift.cram`, 16 threads, same args as B'/B″):**
  - **B″** = `--bait-records primary --no-clips --no-indels --no-linked`: **794 reads** (746 `bait`
    + 48 `unmapped`), 505 truth reads, recall 69.0 %; myloasm → **1 contig, 40,865 bp, circular**.
    The read set is **set-equal to standard A's 794 reads** (symmetric difference 0) and the
    reassembly reproduces A's `u728ctg` 40,865 bp circular contig.
  - `supplementary_only_reads = 9` in both modes; under `primary` those 9 leave the `bait` reason
    (755 → 746).
  - **B‴** = `--bait-records primary` with the default categories: **912 reads** (`bait` 746, `clip`
    104, `large_ins` 14, `large_del` 6, `linked` 35, `unmapped` 48) = exactly B's 912-read set (the 9
    supplementary-only reads are still selected by `clip`/`linked`), recall 72.1 %; myloasm → 3
    contigs, best 70,622 bp, non-circular — identical to B.
  - Net: `primary` reproduces A when the flanking categories are off (B″) and leaves B unchanged
    when they are on (B‴).
- **Alternative rejected:** a separate `--bait-primary-only` boolean; the two-way enum reads better
  and mirrors `--scope`.

## Round trips on the real data sets (user, 2026-10-17)

- **What:** every case of `data/real` (24), `data/val/real` (24) and `data/v2/real` (22; the two
  with a stale `new.fa.fai` excluded) lifted A→B, then the lifted CRAM B→A, then A→B again (master
  8c57697, 16 threads, 210 lifts, all `contract_ok`). Each step compared with its full remap
  (`compare.py`): `oracle.cram` for →B; `old.cram` for →A, a valid full-remap oracle in all 70
  cases (same minimap2 2.31-r1302, preset, `--secondary=no -Y`, identical `@SQ`).
  Scripts and per-case results: `scratch/real-roundtrip/`.
- **A→B and A→B→A→B vs `oracle.cram`:** identical to the single-lift baselines (jt3/jt4) case by
  case: T mismatch 9 / 3 / 18 (real / val / v2), the known `-` block co-optimal CIGARs and one
  competitor-rule MAPQ (60 vs 55, FAM20840); supplementary differ 0; depth identical.
- **A→B→A→B equals A→B read for read in all 70 cases** (0 equivalent, 0 mismatch): no drift.
- **A→B→A vs `old.cram`:** no contig, strand or MAPQ difference in 7.64 M reads; CIGAR mismatches
  0 / 2 / 3 (plus 161 / 8,972 / 18,302 equivalent CIGARs), all reads translated through a `-`
  block with the same span and score as the oracle. Cause: `-` translations left-align indels
  (`left_align_cigar`, see M2), which is deterministic and composes to the identity on the next
  step, hence the zero drift above. No new mismatch class.
- **Translated fraction is not chain quality:** the route changes between steps (up to +31,727
  translated reads at B→A on v2) while the output is identical; report chains by agreement, not
  by translated fraction.
- Time: ≈ 40 min for all three sets (lifts 1,630 s, comparisons overlapped).

## split_junction near contig ends: origin-only rule tried, not adopted (user, 2026-10-18)

- **Problem:** `split_junction` drops every junction with a breakpoint within
  `--junction-edge-excl` (2 kb) of an end of a wrapping contig. Measured
  (`scratch/interjunction/report.md`): on circular contigs < 200 kb this blanks a median 18–20 % of
  the length; on raw -> polished it suppresses 21 passing clusters that are not the end -> start join.
- **Tried (user):** (1) exclude only the origin join (tandem key, breakpoints within E of opposite
  ends). It reported all 21 clusters, but 18 of them are origin crossings whose continuation was
  placed on a copy of the contig start (+38 false alarms on s34). (2) A read-continuation test: drop
  a near-end cluster when the reads' continuations fit across the origin at least as well as at the
  partner (prototype `scratch/origin-cont/report.md`: clean separation on the existing data; qc
  implementation with a seeded banded aligner, targeted SEQ fetch, whole next segment up to 20 kb).
  Final state (s37 / sp35 / fb4): recall unchanged, 16 of the 18 dropped, +6 scored orange false
  alarms vs s33, +0.2-0.54 s on cases with candidates.
- **Not adopted:** three review rounds found real defects each time (a geometry shortcut that hid a
  genuine 2.2 kb deletion across the origin and an inversion; hard-clipped SEQ; cluster membership;
  seed range; weak-hit votes); the last round still had two medium findings (a double orientation
  conversion that makes the verdict depend on which alignment is primary; origin anchoring that
  ignores an end overlap, so the verdict depends on which repeat copy the aligner chose). The gain
  was 3 genuine near-end junctions; the user's fallback criteria (no remaining high/medium findings,
  false alarms not above s33) were not met. Master's 2 kb exclusion stays; the limitation is
  documented in `docs/validation/qc-validation.md` §8.
- **Parked:** the full change (`scratch/origin-cont/parked/origin-only-and-continuation.patch`
  against `BASE`, plus `origin.rs`), the prototype, and the review reproductions
  (`scratch/origin-cont/review2/`). Picking it up needs the two open findings fixed and, ideally, a
  truth set with misjoins at repeat copies of the contig start (the remaining blind spot: a copy
  longer than the tested continuation).
