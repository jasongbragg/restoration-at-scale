# Next phase: scaling genome size and moving to a local server

*Companion to PROJECT_SUMMARY.md. This is forward-looking — none of this has been implemented or tested yet, it's a plan to work from.*

---

## 1. The genome-size decision: locked at 11x4.5Mb

Decided: **11 chromosomes x 4.5Mb**, for the sake of ROH resolution (a 22x2.25Mb split would cap detectable ROH length lower, which matters given the explicit interest in reusing the ROH-nestedness pipeline). Parallelising by chromosome across threads still gives the ~10x-class speedup either way — the choice between 11 and 22 pieces doesn't change that, since both are embarrassingly parallel at the chromosome level.

Per the architecture-independence finding (PROJECT_SUMMARY.md, Section 6), the migration rate doesn't need separate re-calibration for this genome size specifically — m≈0.0042 (the corrected value, see below) was confirmed directly at the full 11x4.5Mb scale.

## 2. Efficient execution on a 100-thread server — and the real parallelism ceiling — and the real parallelism ceiling

**msprime itself is single-threaded per call** — the C library doesn't parallelise internally. The workload here is embarrassingly parallel instead: many independent jobs, each cheap on its own. Parallelise *across* jobs, not within one simulation call.

**Can this actually use 100 threads? Not at Stage 1 alone.** The Stage 1 ground truth only needs to be built *once* — there's no benefit to computing the same chromosome's ancestry twice — so the natural job list for "build the ground truth" has exactly 11 entries (one per chromosome), capping parallelism at 11 regardless of how many cores are available. Two ways to actually use more of a 100-thread machine right now, rather than waiting for Stage 2/3:

- **Multiple independent replicate "ground truth worlds"** (different random seeds), if there's value in checking results aren't an artefact of one particular random ancestral history — e.g. 11 chromosomes x 10 replicate worlds = 110 jobs, which uses the full 100 threads with room to spare. This has genuine scientific value beyond just exercising the harness, not just a way to burn cycles.
- **Calibration/exploration sweeps** — every (parameter value x replicate) combination during a calibration search is independent, and there are usually many more of these than 11 (the migration-rate re-calibration this session alone ran 5 values x 10 replicates = 50 jobs, serially, that could have been simultaneous).

**The real, natural 100-thread workload is Stage 2/3** — every (collection design x replicate) combination is an independent SLiM run, and a realistic version of the N-grid (4 site counts x several N totals x several replicates) easily produces hundreds of independent jobs. Stage 1's 11-way ceiling is a feature of *that specific stage's* one-time-build nature, not a limitation of the harness.

`stage1_parallel_groundtruth.py` (in the scripts folder) is written generically over a job list — `[(chrom, replicate_world) for ...]` — specifically so the same harness code will fan out properly once Stage 2/3 produce hundreds of jobs, rather than needing to be re-engineered later. It defaults to the 11-job (single world) case for now, with a flag to scale up to multiple replicate worlds.

**Practical notes regardless of job count:**
- Use `concurrent.futures.ProcessPoolExecutor` (or shell-level `xargs -P` / GNU `parallel` / a job array) to fan work out.
- Check memory headroom before scaling concurrency up — run one job alone first, note peak RSS, then set concurrency = min(n_jobs, total_RAM / per_job_RAM).
- Each chromosome is saved as its own `.trees` file immediately on completion, so a long batch can be interrupted/resumed without re-doing finished work.

## 3. Efficiency upgrade: move off `genotype_matrix()` + hand-rolled Gst

The calibration prototype materialised a dense genotype matrix and used a hand-corrected Nei's Gst estimator — fine at prototype scale, but worth replacing before scaling up sample sizes and genome size:

- **Use tskit's native statistics** (`ts.Fst()`, `ts.diversity()`, `ts.divergence()`, etc.) instead. In `mode="branch"`, these operate directly on the tree topology and don't require simulating mutations or materialising a genotype matrix at all for many neutral statistics — faster, and they sidestep the small-sample bias issue we had to hand-correct for, since they're proper peer-reviewed estimators.
- They support windowed computation natively, useful for per-chromosome or per-region breakdowns without manual array slicing.
- This is worth doing *before* scaling up rather than after — hand-rolled stats get slower and more bug-prone as data volume grows, and we already found one bias bug in the prototype version.

## 4. Storage strategy

**Use tree sequences (`.trees`, tskit's binary format) as the primary storage unit.** Far more compact than VCF for this kind of data, and they support efficient downstream resampling (`ts.simplify()`) and mutation overlay without re-running the coalescent simulation.

Suggested layout:
```
groundtruth/
  chr01.trees ... chr11.trees       # one-time ancestry build, read-only after creation
designs/
  sites16_mothers64_rep03.trees     # sampled/derived tree sequences per design x replicate
manifest.csv                        # one row per file: all sim parameters + key summary stats
```

A flat manifest (CSV or Parquet) recording every file's parameters (K, m, deme_Ne, chromosome lengths, r, mu, seeds) and headline stats is probably enough — no need for a database system, and it stays consistent with your existing R/Python-centric workflow conventions rather than introducing a new tool.

Tree sequence files are already a fairly compact binary encoding; further compression (zstd/gzip) is a nice-to-have if disk is genuinely tight, not a priority.

## 5. Containerization — worth it?

**msprime/tskit installation is genuinely trivial via pip** on essentially any Linux server, so a container isn't *needed* for that piece alone — your instinct there is right.

**Where it earns its keep is the whole toolchain, once SLiM enters the picture.** SLiM is a compiled C++ binary (not pip-installable), and pinning msprime + tskit + pyslim + SLiM + R consistently across jupiter, this new local server, and rainforest/rcn01 is exactly the kind of cross-machine drift a container prevents — and you've already got Docker workflow experience from the RepeatModeler/RepeatMasker side of things, so this isn't a new skill, just the same pattern applied here.

- No runtime performance cost for CPU-bound simulation work inside a container — the benefit is reproducibility/portability, not speed.
- If this might ever move back onto the HPC, consider building as Singularity/Apptainer (no root daemon, natively supported by most HPC schedulers) — or build a Docker image and convert it (`singularity build my_image.sif docker://<tag>`) to get both.
- **Net call: worth doing once Stage 3 (SLiM) starts, not urgent for the Stage 1 msprime-only work happening right now.** A starter Dockerfile sketch is included alongside this document (`Dockerfile.example`) — untested end-to-end in this environment, treat it as a draft to adapt rather than a finished artifact.

## 6. Suggested concrete next actions, in order

1. ~~Deliberately choose deme_Ne~~ — done, resolved analytically (PROJECT_SUMMARY.md, Section 5)
2. ~~Decide 11x4.5Mb vs 22x2.25Mb~~ — done, locked at 11x4.5Mb
3. ~~Refactor Fst/diversity statistics onto tskit-native methods~~ — done, and it caught a real bias in the old estimator; migration rate revised to m\u22480.0042
4. ~~Stand up the parallel execution harness~~ — done (`stage1_parallel_groundtruth.py`); ceiling at Stage 1 is 11 jobs (one per chromosome) unless running replicate worlds, see Section 2
5. Build the container once SLiM enters the pipeline (not before) — still pending, by design
6. **New, carried forward from this session:** decide whether to target the long-term Ne (>10\u2075) or the post-contraction recent value for deme_Ne, if the distinction matters for the question being asked
7. **New:** calibrate the recombination rate against the empirical LD-decay curve (R\u00b2 half-max at 711bp) if exactness matters before moving to Stage 2
