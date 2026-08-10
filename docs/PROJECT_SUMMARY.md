# Restoration-at-scale simulation framework
## Model development & preliminary parameterization

*Status: Stage 1 (coalescent ancestry, msprime) parameterization complete at prototype scale. Stages 2-4 designed but not yet implemented.*

---

## 1. Project aim

Investigate how seed-collection design for ecological restoration — how many sites to visit, how many maternal trees to sample per site — affects downstream genetic diversity outcomes, at a scale (1,000s-10,000s of seedlings) beyond what prior small-*n* empirical studies have addressed. Approach: population simulation, calibrated where possible against real *Melaleuca quinquenervia* data.

## 2. Pipeline structure (4 stages)

| Stage | Purpose | Tooling |
|---|---|---|
| 1. Coalescent ancestry | Simulate a large, weak-IBD source metapopulation. Treated as ground truth once built — sampled repeatedly for different collection designs without re-simulating. | msprime |
| 2. Collection design sampling | Factorial site x mothers-per-site grid; genotypes for sampled seedlings drawn from the Stage 1 ground truth via mutation overlay / resampling. | msprime |
| 3. Forward restoration simulation | A few generations of the founded population, outcrossing, random planting. Tree-sequence handoff from Stage 1 founders into a forward simulation. | SLiM + pyslim |
| 4. Diversity outcome analysis | Rarefaction-style accumulation curves and a practitioner-facing decision table, comparing restored-population diversity to wild background diversity. | R / Python |

## 3. Key design decisions locked in

| Decision | Choice | Rationale |
|---|---|---|
| Metapopulation geometry | Linear (non-circular) stepping-stone | Matches *M. quinquenervia*'s roughly linear coastal NSW distribution (Guo et al. 2026, Fig. 1a), rather than an arbitrary island model |
| Genome architecture | 11 chromosomes | Real Myrtaceae / *M. quinquenervia* karyotype, confirmed (not guessed) via Guo et al. 2026 |
| Per-chromosome length | **Not yet decided** | Decoupled from migration calibration (see Methodological Lessons, below) — to be set by ROH-resolution needs for the next phase |
| Mating system | Outcrossing, no selfing | Eucalypts and Melaleucas are predominantly outcrossing |
| Planting pattern | Random (no family clustering) | Realistic at restoration scale — clustering can't be controlled the way it can in small plantings |
| Pollen contamination from wild remnant | Deferred, but solved in principle | Reuse Stage 1 ground-truth haplotypes as migrant gamete donors in a SLiM `reproduction()` callback — no need for a second simulation engine or for forward-simulating the (potentially enormous) neighbouring wild population |
| Mother-tree relatedness | Start with unrelated draws | Matches the real collection protocol (trees spaced \u226520 m apart); kinship-imposed comparison deferred to a later iteration |
| Founder hand-off | pyslim recapitation pattern, run in reverse | msprime generates founder genotypes; SLiM receives them as the starting tree sequence for the forward phase (the usual recapitation idea, just pointed the other direction) |
| Deliverable format | Both rarefaction curves AND a decision table | Nested within total collection size N \u2208 {1024, 4096, 16384} |

**Flagged but not yet started:** extending the framework beyond neutral diversity capture to also track capture of the myrtle-rust resistance loci identified in Guo et al. 2026 (two regions on Chromosome 8, one on Chromosome 4) — i.e. does a collection design optimised for neutral diversity under-sample a rare but ecologically critical resistance haplotype? Noted as a "beyond Hoban et al." extension specific to this system, not pursued yet.

## 4. Real-world calibration anchors

All from Guo et al. 2026 (*Molecular Ecology* 35:e70413), *M. quinquenervia* / myrtle rust genomic prediction study:

| Quantity | Value | Used for |
|---|---|---|
| Global Fst (NSW range) | 0.05 | Stage 1 migration-rate calibration target |
| Effective population size | >10\u2075 for ~450,000 generations, then recent contraction | Order-of-magnitude anchor for metapopulation Ne (recent post-contraction value not yet extracted — see paper's Fig. S4 if exact recency matters) |
| Chromosome count | 11 | Genome architecture |
| LD decay | R\u00b2 falls to half-max at 711 bp | Sanity check for recombination-rate realism (not yet used to calibrate `r` directly) |
| Mutation rate (their SMC++ analysis) | 1\u00d710\u207b\u2078 / bp / gen | Adopted directly for our simulations |
| Real collection protocol | 16 sites, 6\u201322 trees/site (mean ~12), trees \u226520 m apart | Empirical bound on "mothers per site" realism in the Stage 2 design grid |

## 5. Stage 1 (msprime) parameterization — progress

### deme_Ne: resolved analytically, not by simulate-then-infer

Running a full demographic-inference method (e.g. SMC++) on our own simulated data to "recover" Ne would be circular — it would just validate our own arithmetic, not provide an independent check. Instead, deme_Ne was derived from a closed-form relationship between subdivision, Fst, and effective size (Whitlock & Barton 1997-style result for structured populations): a subdivided population maintains a *higher* total effective size than an unstructured population of the same total census size, approximately

```
Ne_total ≈ N_sum / (1 - Fst)        where N_sum = K x deme_Ne
```

Rearranged: `deme_Ne = Ne_target x (1 - Fst) / K`. With Ne_target = 1x10^5 (the long-term value from Guo et al. 2026), Fst = 0.05, K = 16:

```
deme_Ne = 100,000 x 0.95 / 16 ≈ 5,938
```

This lands almost exactly on the 6,000 placeholder that had been carried along since the K=8 prototype — reassuring, but now arrived at deliberately rather than by inheritance. [SPECULATIVE: the Whitlock & Barton-style formula is derived under island-model-like assumptions; applicability to a 16-deme linear chain specifically is approximate, not exact.]

A free corroborating check, requiring no extra simulation: branch-mode diversity (Section 5b, no mutations needed) gives the whole-metapopulation pairwise diversity directly in coalescent time units. For a haploid-convention model, expected pairwise branch length ≈ 2 x Ne_global. The 4.5Mb production-architecture run gave Ht ≈ 205,000 generations → implied Ne_global ≈ 102,500 — close to the 100,000 target, from a statistic we were already computing for Fst anyway.

**Open: long-term Ne (>10^5 sustained for ~450,000 generations) vs. the more recent post-contraction value** — the paper reports putative recent contractions without giving the exact recent figure in the main text (see their Fig. S4 if this distinction matters for the question being asked).

### Calibration history (superseded values struck through)

| Parameter | Value | Status |
|---|---|---|
| K (deme count) | 16 | Locked |
| deme_Ne | 6,000 (≈5,938 from theory) | Locked — see derivation above |
| Genome architecture | 11 chromosomes x 4.5 Mb | **Locked** — chosen for ROH resolution; confirmed architecture-independent of m, so no need to recalibrate when this was decided |
| Migration rate (m) | ~~0.0024~~ → **0.0042** | **Revised** — see below |
| Mutation rate | 1x10^-8/bp/gen | Matches Guo et al. 2026 |
| Recombination rate | **5x10^-7/bp/gen** | **Locked** — calibrated against the empirical LD-decay curve (R\u00b2 half-max at 711bp); see derivation below. Supersedes the 1x10^-8 placeholder (50x higher) |

### Important correction: the Fst estimator was still biased, and the migration rate has changed

Refactoring onto tskit's native branch-mode diversity statistic (Section 5b) — done specifically to remove the need for the hand-rolled, small-sample-corrected Nei's Gst — surfaced a real problem rather than just tidying the code: the old estimator (`nei_gst_corrected`) only corrected the within-deme (Hs) term for small-sample bias, leaving the pooled (Ht) term uncorrected. That residual bias was smaller than the one we'd already fixed, but not zero, and it was large enough to matter: recomputing Fst at the old "calibrated" m=0.0024 using the unbiased branch-mode statistic gives **Fst ≈ 0.086**, not 0.05.

Re-calibrating against the trustworthy statistic:

| m | mean Fst (branch-mode) | std |
|---|---|---|
| 0.0038 | 0.0566 | ±0.0026 |
| 0.0040 | 0.0533 | ±0.0023 |
| 0.0042 | 0.0490 | ±0.0027 |
| 0.0044 | 0.0470 | ±0.0027 |
| 0.0046 | 0.0480 | ±0.0014 |

**Revised calibrated value: m ≈ 0.0042.** Confirmed at the full production architecture (11 x 4.5Mb, branch-mode, no mutations simulated): per-chromosome Fst ranged 0.043–0.053, mean 0.0482 — right on target.

### 5b. tskit-native statistics (replaces the hand-rolled estimator)

`branch_gst(ts, K, n_per_deme)` in `stage1_utils.py` computes Fst directly from tree topology via `ts.diversity(..., mode="branch")` — no mutation simulation needed at all for this statistic, no small-sample correction needed (branch-mode diversity is an unbiased pairwise-distance estimator by construction), and it's a peer-reviewed implementation rather than a hand-rolled one. `nei_gst_corrected` is retained in the module for reference/comparison but is superseded.

### 5c. Recombination rate, calibrated against the empirical LD-decay curve

Target: Guo et al. 2026 report R\u00b2 falling to half its near-zero-distance value at 711bp separation.

**Closed-form starting point** (Sved 1971: E[r\u00b2] \u2248 1/(1+4Nc), c = r x distance): solving for the half-decay point gives `r = 1/(4 x deme_Ne x 711) \u2248 5.86x10^-8`. Used deme_Ne (not the global metapopulation Ne) since short-range LD decay is a within-deme process, distinct from the population-structure floor that flattens the real decay curve at long range (that floor is the Fst/migration effect, already locked in separately — folding it into this calibration would double-count structure).

**This theoretical value didn't hold up empirically.** Simulating single-deme genotypes (no migration, isolating the within-deme process) and directly measuring R\u00b2 by physical distance bin, the theoretical r=5.86x10^-8 produced decay far slower than the real curve — even at 1500-2200bp, simulated R\u00b2 was still ~65% of its near-zero value, nowhere near the target ~50% by 711bp. Likely contributors: the idealized formula assumes E[r\u00b2]\u21921 at d=0, but the actual sample-based estimator (finite n, finite sites) never gets there and plateaus above 0 at long range too — the same kind of estimator-reality gap that bit the Fst calibration.

**Empirical sweep, same method as the Fst recalibration:**

| r | near-zero R\u00b2 | 711-1000bp R\u00b2 | ratio (target 0.5) |
|---|---|---|---|
| 1x10^-7 | 0.190 | 0.157 | 0.829 |
| 3x10^-7 | 0.190 | 0.107 | 0.564 |
| 5x10^-7 | 0.181 | 0.088 | 0.488 |
| 1x10^-6 | 0.176 | 0.072 | 0.410 |

**Calibrated value: r \u2248 5x10^-7/bp/gen** — confirmed with 10 replicates: near-zero R\u00b2=0.181, 711-1000bp bin R\u00b2=0.089 (ratio 0.493, essentially exact). About 8.5x higher than the closed-form estimate, and 50x higher than the original 1x10^-8 placeholder.

Caveat: this was calibrated on a single, unstructured deme deliberately, to isolate the short-range decay process from the structure-driven floor. The full structured metapopulation's LD curve will show both effects layered together (fast within-deme decay near r\u22485x10^-7, flattening to a Fst-driven floor at long range) — consistent with, but not yet directly checked against, the shape of the real curve in Guo et al. 2026's Figure 2.

### 5d. Production-scale tractability: the real cost, and a bug caught by integration testing

Once the locked migration rate (m=0.0042) and locked recombination rate (r=5x10^-7) were combined at full production scale (K=16, 11x4.5Mb), a single chromosome became computationally intractable on a single core — neither parameter alone was expensive to simulate; together, the combination of high recombination and 16-deme migration is. Scaling tests at shorter lengths (10kb-300kb) showed superlinear cost growth (~L^1.6-1.9), extrapolating to roughly tens of minutes to over an hour per full-length chromosome on a single core.

**First attempted fix: merge all 16 demes into one ancestral population at T=30,000 generations**, capping how far back migration-driven complexity needs tracking. This gave an apparent ~7x speedup (300kb: 73s -> 10.5s) with Fst unchanged (0.047 vs 0.050) — looked like a clean win.

**It wasn't quite right.** Running it through `stage1_validate_groundtruth.py`'s implied-Ne diagnostic (Section 5b's free Ht/2 byproduct) immediately caught a problem: implied Ne_global = 30,395, nowhere near the 100,000 target. The bug: the ancestral population had been sized at `deme_Ne` (6,000) rather than the metapopulation's true global effective size (~101,053, the same Whitlock & Barton-style figure from Section 5). Beyond the merge point, lineages were drifting at the wrong (much smaller) Ne, truncating total diversity — Fst (a *relative* differentiation measure) was unaffected, but absolute diversity level was wrong, which matters for every downstream capture-rate comparison.

**Corrected: ancestral population sized at 101,053, not 6,000.** Confirmed: implied Ne_global = 101,027 (target 100,000), Fst = 0.0485 (target 0.05) — both correct now. But the real speed benefit of the merge, once sized correctly, is much more modest: ~2x at moderate lengths (300kb: 73s unmerged -> ~33s merged-correctly), not the ~7x first seen. **Tried and ruled out as a further fix:** msprime's `smc_prime` approximate ancestry model gave no measurable improvement over the default (`hudson`) model at this combination of parameters; varying the merge depth (5,000 vs 10,000 vs 30,000 generations) also made no real difference. The dominant cost is fundamentally the deep coalescent time needed to reach Ne~100,000 under r=5x10^-7 — not really avoidable, since that depth of history is what the real species actually has.

**Confirmed on real hardware: 6,886.1s (114.8 minutes, ~1.9 hours) for one full-length chromosome** — right at the upper-middle of the sandbox-extrapolated 30–120 minute bracket, so the extrapolation held up reasonably well despite being projected over more than an order of magnitude of chromosome length.

**Practical upshot:** since the 11 chromosomes are independent, running all of them simultaneously brings total wall-clock for the whole ground-truth build down to ~115 minutes (not 11x that, ~21 hours, which is what serial execution would cost). On a 100-thread machine there's also substantial headroom left over — up to ~9 independent replicate "worlds" (99 of 100 threads) could run in that same ~115-minute window, which is the natural way to spend the spare capacity if replicate-world robustness checks are wanted (see NEXT_PHASE_PLANNING.md Section 2).

**The methodological point worth keeping:** this bug was invisible at the level of any individual calibration (Fst calibration didn't need the merge at all; the merge was only introduced when assembling the full production build) — it only surfaced once the validation script checked a *combination* of locked parameters together. Integration-level diagnostics caught something none of the component-level calibrations could have.

### 5e. Current locked production parameters (cross-check against the scripts)

| Parameter | Value |
|---|---|
| K | 16 |
| deme_Ne | 6,000 |
| m | 0.0042 |
| r | 5x10^-7/bp/gen |
| mu | 1x10^-8/bp/gen |
| Chromosomes | 11 x 4.5Mb |
| ancestral_merge_time | 30,000 generations |
| ancestral_Ne | 101,053 (NOT deme_Ne — see 5d) |
| n_per_deme (ground truth) | 12 |

## 6. Methodological lessons learned

1. **Naive Fst/Gst from small per-deme samples is upward-biased.** With only a handful of haploid samples per deme, apparent Fst plateaued around 0.12\u20130.15 regardless of migration rate — pure sampling noise in the within-deme heterozygosity estimate, not real structure. A Nei (1978)-style small-sample correction (n/(n-1) on Hs) was necessary before migration rate had any visible effect at all. Relevant again later: the same bias will inflate apparent differentiation between real collection sites at low mothers-per-site.

2. **Stepping-stone chain length matters a lot for *global* Fst.** A K=16 chain needed roughly 2.4x the local (neighbour-to-neighbour) migration rate of a K=8 chain to hold the same chain-wide Fst, because indirect gene flow across more steps is weaker even when each step's rate is high. This is a real modelling choice (it sets how many independently-addressable demes are available downstream), not a computational nuisance.

3. **Chromosome count/length doesn't change the *expected* Fst, only its precision.** Porting the K=16 calibration from one 300kb locus to 11 independent (unlinked) chromosomes left the calibrated `m` essentially unchanged — exactly as population-genetic theory predicts, since Fst is a property of drift/migration/Ne, not of how many linked sites you average over. What changed was precision: pooling sites across 11 independent ancestries cut the standard error by roughly 5\u201310x relative to one locus, for similar per-chromosome compute cost. **Practical takeaway: to tighten an Fst estimate cheaply, add independent loci rather than lengthen one.**

4. **A partially-corrected bias estimator is still a biased estimator.** The original small-sample correction fixed the within-deme (Hs) term but left the pooled (Ht) term uncorrected. That gap was smaller than the error it fixed, but not negligible — it shifted the calibrated migration rate by nearly 2x once checked against an unbiased estimator (tskit's branch-mode diversity). **Practical takeaway: prefer a properly unbiased, ideally peer-reviewed estimator over a partial hand-correction, even when the partial correction "seems to work" — it can still be quietly wrong by an amount that matters.**

5. **Closed-form theory is a starting bracket, not a substitute for empirical calibration.** Sved's LD-decay formula gave a clean closed-form estimate for the recombination rate (r\u22485.86x10^-8), but it assumes an idealized estimator that reaches R\u00b2=1 at zero distance — real (and simulated) sample-based R\u00b2 never gets there, and plateaus above zero at long range too. The empirically-fit value (r\u22485x10^-7) differed by 8.5x. Same pattern as the Fst recalibration: trust the simulated/empirical fit over the textbook formula, but use the formula to know where to start looking.

6. **Component-level calibration doesn't catch integration-level bugs.** Every parameter (m, r, deme_Ne) was individually calibrated and individually correct — but combining them at full production scale required a speed fix (the ancestral merge), and the first version of that fix silently broke total diversity (wrong ancestral Ne) in a way none of the component calibrations could have caught, since none of them used the merge. Only `stage1_validate_groundtruth.py`, checking the actual combined production output, surfaced it. **Practical takeaway: always validate the assembled, full-scale artifact, not just its individually-calibrated pieces.**

## 7. Open items

- [x] Deliberately choose deme_Ne — resolved analytically via the Whitlock & Barton-style subdivision/Fst/Ne relationship (Section 5)
- [x] Decide per-chromosome length — locked at 11 x 4.5Mb, chosen for ROH resolution
- [x] Refactor Fst estimator onto tskit-native statistics — done; also caught and corrected a residual bias in the old estimator
- [x] Calibrate recombination rate against the empirical LD-decay curve — done, r\u22485x10^-7/bp/gen (Section 5c)
- [x] Make the production build tractable at full scale — done, with a real bug caught and fixed along the way (Section 5d)
- [x] Build a diagnostics/validation script for the assembled production output — done (`stage1_validate_groundtruth.py`)
- [ ] **Confirm actual full-length (4.5Mb) per-chromosome runtime on the real server** — sandbox extrapolation only, uncertain over more than an order of magnitude; test before committing to a full batch
- [ ] Implement the relatedness-imposed mother-tree sampling comparison (currently deferred)
- [ ] Implement the pollen-contamination module (conceptually solved, not built)
- [ ] Consider extending the framework to track myrtle-rust resistance locus capture, not just neutral diversity
- [ ] Decide whether to target the long-term Ne (>10\u2075) or the more recent post-contraction value (see paper's Fig. S4)
- [x] Confirm the LD-decay shape holds up on the full structured 11x4.5Mb metapopulation — checked, WARN-level deviation, understood and accepted (see Section 5f)

## 5f. Stage 1 ground truth: VALIDATED (world00, full 4.5Mb x 11 production build)

Real `stage1_validate_groundtruth.py` output against the completed production build:

| Check | Result | Target | Status |
|---|---|---|---|
| Fst | 0.0500 | 0.05 | PASS |
| Implied global Ne | 100,945 | 100,000 | PASS |
| LD half-decay ratio at 711bp | 0.580 | 0.5 | WARN (16% deviation) |

Fst and Ne landed almost exactly on target. The LD WARN is explained, not a bug: the recombination calibration (r=5x10^-7, Section 5c) was fit on a fully migration-isolated population, but no deme in the actual locked metapopulation (m=0.0042) is migration-isolated — some of "deme0"'s sampled individuals carry admixed ancestry, and admixture-driven LD (Wahlund effect) decays much more slowly with distance than pure recombination-driven LD, flattening the curve at longer range than the isolated-population calibration predicted. This is a smaller-scale version of the same mechanism behind the long-range plateau in the real paper's LD curve (Guo et al. 2026, Figure 2) — structure bleeding into a same-deme sample because no real deme is migration-free either.

**Decision: accepted as-is for now.** Fst and Ne are the parameters that matter most for the Stage 2/3 diversity-capture questions this framework is built around; the LD figure being 16% off, for an explained reason, doesn't block proceeding. Revisit only if a downstream analysis turns out to be sensitive to the exact short-range LD decay rate specifically.

**Stage 1 is now considered complete.** Next: Stage 2 (collection-design sampling) implementation.
- [ ] Begin Stage 2 (collection-design sampling) implementation

## 8. Scripts

See accompanying files. Two tiers: calibration history (how each parameter was derived) and production (what to actually run on the server).

**Production — run these on the server:**
- `stage1_build_groundtruth.py` — **the** production script. Every locked parameter listed explicitly at the top of the file. Builds the 11-chromosome ground truth via `ProcessPoolExecutor`, saves `.trees` files to `groundtruth/`.
- `stage1_validate_groundtruth.py` — diagnostics script. Run this immediately after a build. Checks Fst (target 0.05), implied global Ne (target ~1e5, free byproduct), and LD half-decay distance (target ~711bp) against the actual produced output, with PASS/WARN/FAIL bands.

**Calibration history — how each locked value was derived, kept for reference/reproducibility:**
- `stage1_utils.py` — shared demography builder, both the legacy (superseded) and current tskit-native Fst estimators, and the chromosome-build helper used by both calibration and production scripts
- `stage1_calibrate_single_locus.py` — original K=8 vs K=16 single-locus calibration sweep (historical — uses the now-superseded estimator)
- `stage1_calibrate_11chrom.py` — architecture port to 11 chromosomes (historical — same superseded estimator)
- `stage1_recalibrate_tskit_native.py` — documents the Fst bias discovery and the revised calibration (m \u2248 0.0042) using the unbiased branch-mode statistic
- `stage1_calibrate_recombination.py` — LD-decay measurement and recombination-rate calibration (r \u2248 5x10^-7); also supplies `pairwise_r2_by_distance`/`bin_pairs`, reused by the validation script

All require only `msprime`, `tskit` (installed alongside msprime), and `numpy`. See RUNBOOK.md for server setup and execution steps.

## 9. Stage 2: collection-design sampling — in progress

**Status: core sampling code written and tested; the actual design grid (specific site counts / mothers-per-site / total N values to run) deliberately not yet decided** — per explicit direction, priority here was robust, general, parameterized code over locking in specific numbers.

### What exists

`stage2_sampling.py`:
- `CollectionDesign` — a dataclass specifying `n_sites`, `mothers_per_site` (single int or per-site list), `site_selection` strategy, `seed`, `relatedness_mode`.
- `validate_design_feasibility()` — checked *before* any sampling happens, not discovered mid-failure. Raises specific, actionable `ValueError`s for: too many sites requested (>K), too many mothers requested at any site (> ground truth's `n_per_deme`), and unimplemented `relatedness_mode` values.
- `select_sites()` — pluggable site-selection strategies: `"random"` (uniform subset), `"even_spacing"` (evenly spaced across the linear chain), `"contiguous"` (random contiguous block). Which of these is most realistic for actual collection patterns is an open question, not resolved here — the point is the code supports comparing them later without rewrites.
- `select_mothers()` — without-replacement draw from the available per-deme individuals.
- `sample_design()` — orchestrates the above: loads the relevant ground-truth chromosome tree sequences for a given world, validates, selects, `ts.simplify()`s down to the sampled individuals, overlays mutations (the ground truth itself is stored mutation-free, ancestry only — see Section 5b's reasoning), returns simplified+mutated tree sequences per chromosome plus full selection metadata.
- `save_sampled_design()` — writes per-chromosome `.trees` files plus a `manifest.json` recording the exact design parameters, selected demes, and selected individuals (full provenance, not just a seed to re-derive it from).

### The capacity ceiling, stated plainly

The ground truth was built with `n_per_deme=12` (sized for Fst/diversity calibration, not for being the operational sampling pool). **12 is therefore the maximum mothers-per-site achievable from the current ground truth, for any design.** This is enforced in code (`validate_design_feasibility` raises immediately, with a message pointing at `stage1_build_groundtruth.py` as the fix) rather than silently capping or under-sampling. If/when a real design grid is decided and it needs more than 12 mothers from a single site, the ground truth needs rebuilding with a larger `n_per_deme` first — that's a Stage 1 action, not something Stage 2's sampling code should paper over.

### Tested

- A normal feasible design end-to-end (sample -> simplify -> mutate -> save), with sample counts checked against expectation.
- The over-capacity rejection (64 mothers/site against a 12-cap ground truth) — confirmed it raises rather than silently returning fewer mothers than requested.
- All three site-selection strategies, including their differing outputs on the same seed.
- Unequal per-site mother counts via a list.
- Mismatched-length mothers list, and too-many-sites — both raise cleanly.
- Reproducibility: identical seed -> identical site/mother selection; different seed -> different selection.

See `stage2_demo.py` for a runnable walkthrough of all of the above (illustrative only, not a real design grid).

### Open for later

- The actual design grid (which site counts, mothers-per-site, and total N values to run) — deliberately deferred.
- Which `site_selection` strategy is most realistic — deferred; the code supports comparing them once that's worth doing.
- Outcome/diversity statistics on the sampled designs (rarefaction-style accumulation curves, the eventual deliverable) — not yet built; `sample_design()`'s output (simplified, mutated tree sequences) is the right input for that, once written.
- Increasing the ground truth's `n_per_deme` if/when a design grid needs more than 12 mothers/site at any single site.
