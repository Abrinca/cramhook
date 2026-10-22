# cramlift — Plan

> The tool was renamed from cramlift to cramhook on 2026-10-19.

> Historical design document (Plan v1, 2026-10-07). Superseded by `docs/project/plans/Plan-v3.md`;
> the decisions and measured results are in `docs/project/decisions.md`.

## 1. Goal

Bacterial assemblies change in small ways between iterations: rotated contigs, an added or
removed plasmid, a few medium-sized indels, polishing edits. Each time, all reads are mapped
again from scratch. `cramlift` takes an existing CRAM (reads mapped to the old assembly) and
produces a CRAM for the new assembly that is **equivalent to a fresh minimap2 run**, but faster:

1. Compare old and new assembly with minimap2.
2. **Translate** reads that align cleanly inside regions that are unchanged and unambiguous:
   only contig, position and strand change; the alignment itself stays byte-identical.
3. **Remap** all other reads with the user's normal minimap2 command, then merge and sort.

Correctness comes first. When in doubt, a read is remapped. The tool must never produce a
worse result than a full remap; at worst it is not faster.

Downstream consumers of the output: IGV, coverage estimation, Clair3 (SNVs/indels),
Sniffles2 (SVs). They define what "equivalent" means (section 7).

## 2. Input contract (v1)

The input CRAM was produced by exactly this command shape:

```
minimap2 -ax $PRESET -t N --secondary=no -Y old.fa reads.fastq.gz \
  | samtools sort -@ M -O CRAM --reference old.fa -o old.cram
```

- `$PRESET` is `map-hifi` (PacBio HiFi) or `lr:hq` (accurate ONT).
- Coordinate-sorted, single-end long reads, unmapped reads included.
- The minimap2 `@PG` line in the header is parsed to verify this. Allowed options: the
  preset, `-t`, `--secondary=no`, `-Y`. Anything else (or a missing `@PG`) violates the contract.
- Expected tags: `NM ms AS nn tp cm s1 s2 de rl SA`. A record with any other tag
  (e.g. `MD`, `cs`, `MM`) is not translated; its read is remapped.

**Contract violation → fallback:** remap all reads (extracted from the CRAM) and log the reason.
With `--strict`, exit with an error instead.

Old and new assemblies: plain FASTA; create `.fai` if missing. CRAM decoding/encoding always
uses an explicit reference (never network lookup).

## 3. Output

- CRAM, coordinate-sorted, indexed, written with `new.fa` as reference.
- Header: `@HD SO:coordinate`; `@SQ` for every contig of `new.fa` in `.fai` order (including
  contigs with no reads, e.g. a new plasmid); original `@PG` lines kept, plus a `cramlift`
  `@PG` line with the full command line, old/new reference paths and the remap command used.
- Each primary read appears exactly once (mapped or unmapped), like minimap2 output.
- Per-record tag `XL:A:T` (translated) or `XL:A:R` (remapped). Option `--xl-tag` to add these tags, off by default.
- `--stats stats.json`: counts per outcome and per deferral reason, timings of each step.

## 4. Algorithm

### 4.1 Assembly diff

Align old (query) to new (target):

```
minimap2 -cx asm5 --cs --secondary=no -t N new.fa old.fa > diff.paf
```

Primary alignments only. Errors in the diff **fail safe**: a missed or broken alignment
only means fewer identity blocks, so more reads are remapped (slower, never wrong), and the
invariant check (4.3) catches incorrect blocks. Do not make the diff repeat-sensitive
(e.g. `-N 1000000 -p 0.000001`): that would taint every unchanged repeat copy and remap
reads that could be translated. Repeats that exist in both assemblies need no taint; reads
inside them have low MAPQ and are remapped by the MAPQ rule.

Large indels may appear either as long gaps inside one alignment or as separate alignments;
the parser must handle both (both end a block).

From the cs strings derive, in **old coordinates**:

- **Identity blocks:** maximal runs of exact matches (`:n` operations) of primary alignments.
  Each block: `old_contig, old_start, old_end → new_contig, new_start, strand`.
  Every SNP or indel ends a block. A rotation yields two blocks per contig, a reverse-complement
  a block on the `-` strand, a rename just a different `new_contig`.
- **Taint regions** (old coordinates) where translation is not safe even if sequence is unchanged:
  - T1: old sequence covered by no primary alignment (deleted/unaligned) or by more than one
    (overlapping alignment ends).
  - T2: old regions similar to **novel new sequence**. Novel = new-assembly intervals not
    covered exactly once by primary alignment targets (new plasmid, insertion, the extra copy
    of an expanded duplication). This catches e.g. a new plasmid that shares an IS element
    with the chromosome: reads at the chromosomal copy may now be ambiguous.

    **T2 errors fail unsafe** (a missed similar region means a translated read keeps a stale
    MAPQ, and nothing downstream notices), so this step must be at least as sensitive as the
    read mapper: a copy of ~90 % identity can already lower minimap2's MAPQ. Extract novel
    intervals plus flanks, cut them into overlapping read-length windows, and map them with the
    read preset, keeping all hits:

    ```
    minimap2 -cx $PRESET -N 1000000 -p 0.000001 -t N old.fa novel_windows.fa
    ```

    Mark all hit intervals (on `old.fa`) as tainted. Window size, step, flank size and any
    minimum hit length/identity are constants to be tuned on S6/S6b (section 6).

  Regions similar to *deleted* old sequence need no extra taint: reads there had low MAPQ in
  the old mapping and are remapped by the MAPQ rule.

### 4.2 Eligibility (per record, one streaming pass)

A record is **translated** only if all hold:

1. primary (not supplementary), mapped;
2. no `SA` tag, unless the read is an **origin-crossing split read** (below). Other split reads
   are always remapped as a whole: all records of a split read carry `SA`, so supplementary
   records are simply dropped and the read is remapped via its primary;
3. MAPQ ≥ `min_mapq` (default 60);
3b. MAPQ instability from a near-equal competitor: minimap2's competing chain score `s2` must be
   `< max_s2_ratio` × `s1` (default `max_s2_ratio` 0.5); a missing/invalid `s1`/`s2` or `s1` ≤ 0
   also defers. A close competitor can flip the MAPQ when its score shifts in the new assembly.
   Read-level reason `competitor` (`split_competitor` for split units); each record uses its own
   `s1`/`s2`;
4. divergence `de` ≤ preset-specific threshold (initial: 0.01 for map-hifi, 0.03 for lr:hq; tune);
5. total soft-clip ≤ `max_clip` (default 200 bp; a long clip may align to novel sequence elsewhere),
   AND no soft clip of ≥ 20 bp shares k-mers with new sequence (novel intervals ± 1 kb). Added in M5:
   in S5c/S6c (chimeras chromosome + new plasmid), lr:hq creates supplementaries for 50–200 bp tails
   on new sequence, so the length limit alone is unsafe. A hit is a diff-level deferral (`clip_hit`);
6. only expected tags;
7. the **window** `[ref_start − left_clip − margin, ref_end + right_clip + margin)` (default margin
   500 bp), clipped to the contig boundaries, lies inside a single identity block and does not
   overlap any taint region.
7b. a record translated through a `-` block must not contain an insertion or deletion ≥
   `max_minus_indel` (default 50 bp): a large indel on `-` is minimap2's own non-mirror direction
   effect at a structural change and can move the position or drop a supplementary. `+` blocks are
   unaffected. Diff-level reason `minus_indel` (`split_minus_indel` for split units).

Every primary record that fails a rule (and every unmapped read) is **deferred**, with the
first failing rule recorded as reason for the stats. Rules 1–6 (and unmapped) are
**read-level** reasons: they depend only on the read and would defer it whatever the assemblies
are. Rule 7 is the only **diff-level** reason. Stats report the two classes separately.

**Origin-crossing split reads** (added after M1: 45–62 % of reads in some ONT cases are split, most
of them crossing the start/end of small circular contigs). Such a read is translated as a unit if
all of these hold, otherwise it is deferred:

- exactly two records (primary + one supplementary) on the same old contig, and the `SA` tags
  are consistent with the two records;
- one record (A) reaches the contig end: `ref_end + junction_clip ≥ contig_len`, where
  junction_clip is its soft-clip on the side facing the end. The other record (B) reaches the
  contig start: `ref_start − junction_clip ≤ 0`;
- each record passes rules 3–7 individually, with rule 5 applied only to the outer
  (non-junction) clip. Rule 3 uses each record's own MAPQ;
- the contig junction is preserved: A's block reaches the old contig end, B's block starts at
  old 0, both map to the same new contig with the same strand, and old end/old 0 map to
  new end/new 0 (or new 0/new end for `-`).

Both records are translated by 4.3 and their `SA` tags are rewritten from the translated records
(contig, pos, strand, CIGAR, MAPQ, NM). Records of a split read are buffered by name until both
are seen (they are at opposite ends of the contig in the sorted input). Reasons for failure are
reported in the stats like any other deferral.

### 4.3 Translation

- New contig and position from the block. For `-` strand blocks: position computed from the
  old alignment end, CIGAR reversed, SEQ reverse-complemented, QUAL reversed, flag 0x10 toggled.
- MAPQ and all expected tags unchanged (they are strand-independent or depend only on the read).
- **Invariant check (always on):** the reference bases under the alignment in `old.fa` must equal
  the bases at the new location in `new.fa` (reverse-complemented for `-`). A violation is a bug
  in the diff: abort with a clear error, never silently continue.

### 4.4 Remap and merge

- Deferred reads are written to a temporary FASTQ from their primary records (reverse-complement
  SEQ and reverse QUAL if 0x10 is set). Supplementary records are never used as source.
- Remap with the same preset and flags: `minimap2 -ax $PRESET -t N --secondary=no -Y new.fa deferred.fq`.
- Translated and remapped records are merged and sorted. v1: stream both into
  `samtools sort -O CRAM --reference new.fa`, then `samtools index`. Optimize later if profiling
  shows the need.
- Final check: number of primary records out == number of primary records in (unmapped included),
  and every read name appears exactly once as primary.

## 5. Implementation notes

- Rust, single binary. CRAM I/O via `rust-htslib` (pinned). minimap2 and samtools are called as
  external programs. `clap` for the CLI, `anyhow` for errors.
- Use distinct types for old vs. new coordinates (`OldPos`, `NewPos`, intervals likewise) so
  mixing them is a compile error. 0-based half-open intervals everywhere internally.
- Blocks and taint: sorted vectors per contig + binary search (few contigs; no interval trees needed).
- Temporary files in a temp dir; `--keep-tmp` keeps them for debugging.
- Python (pysam, numpy) for test data generation and comparison scripts in `scripts/`.
- Consult `ref/` sources and the pinned crate versions for APIs and semantics; do not code
  crate APIs from memory. SAM/CRAM semantics: `ref/hts-specs`.

CLI sketch:

```
cramlift lift --in old.cram --old-ref old.fa --new-ref new.fa --out new.cram \
              [--threads N] [--stats stats.json] [--strict] [--keep-tmp] [--xl-tag]
cramlift diff --old-ref old.fa --new-ref new.fa --out-prefix PREFIX
              # debug: writes blocks.bed and taint.bed (old coordinates) + summary
```

## 6. Test data

### 6.1 Synthetic scenarios (`scripts/simulate.py`, deterministic seeds)

Base genome: random chromosome (~300 kb) and 1–2 plasmids (20–60 kb), with **planted repeats**
(several ~1.2 kb IS-like copies, a few ~5 kb rRNA-operon-like copies), because MAPQ effects only
show up with repeats. Contigs marked circular in the FASTA header.

Each scenario produces `old.fa`, `new.fa`, `truth.json` (list of edits with coordinates):

| ID  | Change (old → new) |
|-----|--------------------|
| S1  | identical (sanity: everything translated except non-origin-crossing split reads) |
| S2  | rotation of the circular chromosome |
| S3  | reverse-complement of a contig |
| S4  | contig renamed |
| S5  | plasmid added (unique sequence) |
| S6  | plasmid added that carries a copy of a chromosomal IS element |
| S6b | as S6, but the plasmid's IS copy is diverged (~90 % identity); tests T2 sensitivity |
| S7  | plasmid removed |
| S8  | insertion of 5–20 kb |
| S9  | deletion of 5–20 kb |
| S10 | collapsed duplication (two copies in old, one in new) and the reverse |
| S11 | 20–50 SNPs and small indels (polishing-like) |
| S12 | combination: rotation + plasmid added + 3 medium indels (the typical real case) |

Reads are simulated from `new.fa` (the "truer" assembly), circular-aware (reads may cross contig
ends), ~50× coverage, two error/length models: HiFi-like (~0.2–0.5 % errors, ~15 kb) and
ONT-hq-like (~1–2 %, homopolymer-biased indels, broader length distribution). Add a few
chimeric reads and a few random junk reads. Read names encode true origin.

For each scenario and preset: `old.cram` via the standard command, and the **oracle**
`oracle.cram` = standard command on `new.fa`. When creating `oracle.cram`, save how it was created 
(precise command) and the time.

### 6.2 Real data

`data/real/<case>/` (see `data/real/manifest.tsv`): `old.fa`, `new.fa`, `reads.fastq.gz` symlinks,
preset in the manifest. Cases are `cons-to-polished` / `cons-to-finalized` and
`<primary>-to-cons`. `old.cram` and `oracle.cram` are generated by a pixi task.

### 6.3 Comparison (`scripts/compare.py`)

Compares `cramlift` output with the oracle and writes JSON + a short text summary:

- per `XL` category: agreement of primary contig/position/strand, CIGAR, MAPQ; list of mismatches;
- supplementary alignments: same set per read (for remapped reads);
- per-base depth (primary, MAPQ ≥ 20, deletions counted) difference along the genome;
- invariants: read counts, each read exactly once as primary.

Also measure the **baseline noise**: compare the oracle with a second minimap2 run where
the reference differs trivially (e.g. an added unrelated dummy contig). minimap2 depends
slightly on the whole reference (minimizer frequency cutoff), so some disagreement is inherent.

## 7. Acceptance criteria (initial; revise with evidence, never silently)

On all synthetic scenarios, both presets:

1. zero invariant-check violations; read-count invariants hold;
2. translated reads: ≥ 99.9 % identical to the oracle in contig, position, strand, CIGAR and MAPQ,
   and not worse than the baseline noise; every mismatch is listed and explained. A differing CIGAR
   counts as agreement and is reported separately as "equivalent" if contig, pos, end, strand,
   MAPQ, read length and `AS` are all equal, i.e. an equally scored alignment of the same span
   under minimap2's own scoring. `NM` and `de` are reported for equivalent pairs but are not part
   of the definition. Amended after M3/M4, user decisions: on '-' strand translations, minimap2's
   indel placement is not mirror-symmetric. Left-alignment alone does not reconcile all cases, because
   S3 contains co-optimal alignments with different I/D order. Translations on '+' (unflipped)
   blocks must be byte-identical to the oracle (bug guard);
3. per-base depth with deletions counted (as `samtools depth -J`; primary, MAPQ ≥ 20): identical
   at ≥ 99.9 % of positions, max absolute difference ≤ 2. Amended after M4: without deletions,
   equivalent alignments that shift a deletion within a homopolymer cause cosmetic differences;
4. contract violation (e.g. CRAM made with `--secondary=yes`) → output equals a full remap;
5. S1 (identical assemblies): the diff yields exactly one identity block per contig spanning
   the whole contig and no taint; **no read is deferred for a diff-level reason**; every
   translated record is identical to its input record apart from the `XL` tag. Read-level
   deferrals (unmapped, non-origin-crossing split reads, MAPQ < 60 in repeats, `de`/clip limits) are expected and
   reported per reason.
6. All scenarios (reported metric, not a gate): fraction of diff-level deferrals whose window
   does not overlap a planted edit or a T2 region justified by `truth.json`. These are reads
   remapped unnecessarily; the number measures how efficient the diff is, not its correctness.

On real data: same agreement metrics, plus wall-clock time vs. full remap
(minimap2 + sort) measured with `hyperfine`. Speed target (not a hard gate): ≥ 2–3× for
`cons-to-finalized` and typical `primary-to-cons` cases.

### 7.1 Abortion criterion

If it is clear that this approach does not save significant compute time (50%), abort this project.

## 8. Milestones

Work one milestone at a time. After each: commit, then **stop and report** (what was done,
test results, deviations, open questions).

- **M0 — Workspace setup.** Done.
- **M1 — Test data + oracle.** Done (results in `docs/project/results.md`). `simulate.py` for S1–S12, both read models; pixi tasks to build
  `old.cram`/`oracle.cram` for synthetic and real cases; `compare.py` (validated by comparing the
  oracle with itself and with the baseline-noise run). Report the baseline noise.
- **M2 — Diff.** `cramlift diff`: PAF/cs parsing, identity blocks, taint T1/T2, BED output.
  Unit tests against `truth.json`: blocks end exactly at planted edits (an ambiguous indel may
  be placed at any equivalent position; tests accept exactly those positions); S6/S6b taint covers the
  chromosomal IS copy; S2/S3/S4 map correctly.
- **M3 — Contract + translate.** Header/`@PG` parsing and contract check, eligibility rules,
  translation with invariant check, origin-crossing split reads, deferral reasons. Output: translated records + list of
  deferred reads (no remap yet). Unit tests for strand flips and coordinate math.
- **M4 — Full pipeline.** Remap, merge, sort, index, header, `XL` tags, stats JSON, fallback.
  Meet the synthetic acceptance criteria.
- **M5 — Real data.** Run all real cases, compare with oracles, benchmark, tune thresholds per
  preset (`de`, clip, margin, taint parameters) using deferral-reason stats. Record results in
  `docs/project/results.md`.
- **M6 — ONT speed on larger datasets** (scoped with the user 2026-10-09; see
  `docs/project/decisions.md` "M6 scope"). Target use case: very similar assemblies built from the same
  raw reads; test proxy = raw assembler output → consensus (slightly harder). Goal: maximal
  speedup on non-tiny ONT datasets; tiny datasets get an automatic full-remap fallback.
  1. Build the M6 ONT set: all 11 samples of `$DATA_ROOT/ont`,
     cons→polished + 3 raw→cons each; read sets > 300 MB downsampled to ~100 MB
     (`seqtk sample`, fixed seed). Dev on all, no hold-out. Measure a baseline with
     deferrals broken down by read **bases**. Cache full-remap timings (no repeated oracle
     minimap2 runs). **Stop and report** before building rules.
  2. Build the levers the user approves from that breakdown (candidates: split-read unit rule,
     finer T2 taint, no_block).
  3. Tiny-data auto-fallback: simple size threshold, chosen from measurements, logged in stats.
  Constraints: all correctness gates unchanged; PacBio M5 cases + synthetic acceptance must
  still pass (PacBio speed may change; report it). Downstream checks were done in
  `scratch/study-ambiguous-healing` (E5); the samtools-sort/diff thread split was measured and
  left unchanged (2026-10-09).

## 9. Out of scope (v1)

Patching alignments across edits (CIGAR surgery), secondary alignments, inputs from other
aligners or other minimap2 options, paired-end reads, transforming `MD`/`cs`/`MM` tags,
lifting from older than the directly preceding assembly. BAM input may work via htslib but
is not a goal.

## 10. Working rules for the agent

- Do not work in global dirs like `/var/tmp` or `/tmp`. Work in `./scratch/{agent-slug}/`.
- Correctness over speed: if a rule is unclear, defer the read.
- Never weaken an acceptance criterion or edit a test just to make it pass; report instead.
- Keep `docs/project/decisions.md`: one short entry per non-obvious decision (what, why, alternatives).
- Do not modify anything in `ref/`. Do not commit `data/`.
- Pin crate versions; when an API is unclear, read the pinned source/docs.
- Ask before adding dependencies outside pixi/Cargo or changing the input contract.
