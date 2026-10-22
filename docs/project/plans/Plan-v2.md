# cramlift — Plan v2 (proposal, 2026-10-09)

> The tool was renamed from cramlift to cramhook on 2026-10-19.

> Historical design document (Plan v2, 2026-10-09). Completed and superseded by
> `docs/project/plans/Plan-v3.md`; the decisions and measured results are in `docs/project/decisions.md`.

## Where v1 stands

- Correct on 84 real cases: depth matches the full remap everywhere except ±1 at one position in
  two cases, and supplementary agreement is 100 %. The only residuals are documented `-` block
  direction effects and one MAPQ near-tie.
- Speed vs a full remap: full-size ONT raw → polished median 2.8×; PacBio 3.9× (M5 build).
- Lift time on a typical full-size ONT case (FAM20853 myloasm_03: 6.0 s vs 24 s full remap):
  diff 0.6 s, translate pass 1.3 s, remap of deferred reads ~2.2 s, sort + CRAM encode ~1.8 s.

## Is another round worth it?

**Yes, but a small one.** The remaining speed levers are either modest or risky:

- Remapping fewer reads (small-edit patching + general split units, `ideas.md`) would save
  ~9 % of lift time on the current test set. It costs a new subsystem that edits CIGARs, which
  Plan §9 excludes from v1.
- Fixed costs (diff, CRAM decode/encode) are near their floor. Encoding can't start before
  minimap2 finishes, and the thread/memory settings are exhausted (measured).

The bigger gap is evidence, not code. We optimized on raw assembler output → polished, a
proxy that is deliberately harder than the real use case (cleaned assembly → slightly different
cleaned assembly, often rotated to *dnaA*). Patching's value depends on that use case: there the
diff consists mostly of polishing edits, so `no_block` may dominate, and the gain could be far
larger than 9 %. So v2 measures the real use case first and builds patching only if the numbers
justify it.

## M7 — real-use-case validation (measure only)

1. **Use-case set** (`data/v2/`, full size; reuse existing CRAMs as in M6):
   - **cons → polished**, full size, for the 6 large ONT samples (same reads; polishing edits
     only).
   - **Rotation:** polished → polished with each circular contig rotated to a different start
     and one contig reverse-complemented (synthetic *dnaA* reorientation of a real assembly).
   - **Different cleaned assembly:** ONT reads; old = ONT-polished consensus, new = PacBio-finalized
     consensus of the same isolate (`$DATA_ROOT/pacbio/<S>-p1-1.1`). These are independent cleaned
     assemblies of one genome.
   - **PacBio:** cons → finalized for the ~6 largest PacBio samples by read size (not used yet).
2. **Benchmark the current build** on this set, on `data/m6large`, and on the M5 PacBio sets.
   The M6 speed-ups have not been measured on PacBio yet. Report speedups, correctness, and
   deferral by bases.
3. **Stop and report.**

## M8 — conditional: small-edit patching

Gate: in the M7 use-case set, at least 10 % of input bases are deferred as `no_block` because of
small edits, and the M7 analysis estimates at least 20 % less lift time. Otherwise skip M8.

- Needs a user decision to lift the Plan §9 restriction (CIGAR surgery).
- Scope: SNPs and indels ≤ 5 bp inside one diff alignment, with no other rule failing.
- The new CIGAR, `AS` and `NM` must equal what minimap2 would produce. Where that can't be
  guaranteed (e.g. score near a tie), the read is deferred.
- New `XL:A:P` category; equivalence measured separately; same correctness gates as v1.
- Generalised split units (`ideas.md` candidate 2) are added only if M7 shows a meaningful share.

## M9 — hardening (independent of M8) — done

- **Tiny-dataset fallback:** if the input is so small that a full remap is cheaper, remap
  everything. Use a simple threshold, chosen from the M6 small-data measurements, and log it in
  the stats.
- **`--threads` default:** use the available cores instead of 8.
- **Fix the debug-mode overflow** in the `paf` test helper, so plain `cargo test` passes.
- **Real-data regression set:** a small, fast subset of `data/m6large` + PacBio as a one-command
  check (`pixi run real-smoke`) for future changes.

## M10 — remap overlap (measure first; user, 2026-10-10) — done, gate not met

minimap2 reads input in batches of `-K` bases (default 500 M) and maps a batch only when it is
full or at EOF. Deferred input is typically 50–150 M bases, so the remap does not overlap the
streaming pass; it runs afterwards ("writer drain", e.g. 3.1 of 10.9 s).

1. Sweep `-K` (e.g. 10M/25M/50M/500M) via the `CRAMLIFT_MINIMAP2` wrapper on ~4 full-size cases
   (ONT + PacBio, low and high deferral). Records must be identical to the default (check vs
   oracle and vs the default lift output).
2. If ≥ 10 % less lift time: make it the default (and record it in `@PG`/stats).
3. Profile the streaming pass (4.2 s vs ~0.9 s decode floor): CPU in the pass thread or sort
   backpressure? Optional: feed diff-independent deferrals (rules 1–6) to minimap2 during the diff.
4. Stop and report.

## M11 — faster streaming pass (user, 2026-10-11; before M10b) — done

M10 profile: the pass thread is CPU-bound (4.1 s vs 1.1 s decode floor, FAM20866 rotated):
record rebuild (`Record::set`/`bam_set1`) 34 %, `translate_one` 29 %, memmove 10 %.

1. Cut per-record work: edit translated records in place (tid/pos/flag/MD removal/aux) instead of
   rebuilding; on `-` blocks rebuild once from the 4-bit sequence; avoid clones/copies.
2. If the pass is still the critical path: translate in a small worker pool with ordered output.
3. Then re-test `-K 20M–50M` with 12–14 remap threads; adopt only with ≥ 10 % median gain, and
   only with a fixed tie order (byte-reproducible output) or an explicit user decision.
4. Gates: records identical to the current build (md5 of `samtools view --no-PG`, order-
   insensitive and, where order is deterministic, ordered); tests, accept, real-smoke; bench on
   the M10 cases + `data/v2`. Commit, stop and report.

(M10b follows M11.)

## M10b — tag reliability (user, 2026-10-10) — done

1. Measure, on translated reads (`XL:A:T`) of all 90 real cases vs the full remap, how often each
   tag differs: `NM ms AS nn tp cm s1 s2 de rl SA` and MAPQ (split by `+`/`-` block). Read-only
   on existing eval outputs where possible.
2. Decide per tag: keep (≥ 99.9 % equal or only documented orientation effects) or strip from
   translated reads (meaningful drift, cannot be recomputed without minimap2's index). Decide
   whether `XL` becomes the default.
3. README: a short section on what `XL` means, how reliable each tag is on translated reads, and
   which tags are removed.
4. Implement, verify (tests, accept, real-smoke), commit, stop and report.

## Not planned

- Finer T2 taint, more sort/diff threads, pre-sorted merge: measured as not worth it
  (`ideas.md`).
- Lifting across multiple revisions, BAM/paired-end input, other aligners (Plan §9).

## Working rules

Unchanged from `docs/project/plans/Plan-v1.md` §10:
- One milestone at a time; commit, then stop and report.
- Never weaken a correctness gate silently.
- Scratch only under `./scratch/`; never write outside the repo.
- Benchmarks run serialized under `scratch/m6.lock`.
