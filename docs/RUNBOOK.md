# RUNBOOK: running Stage 1 on the server

## 1. Environment setup

```bash
pip install msprime tskit numpy --break-system-packages   # or use a venv
```

That's it for Stage 1 -- no SLiM, no container needed yet (containerization is deferred until Stage 3, see NEXT_PHASE_PLANNING.md).

Copy the whole `restoration_sim_project/` folder to the server. All scripts are self-contained beyond the three packages above.

## 2. Confirmed timing (real server measurement)

**One full-length (4.5Mb) chromosome: 6,886.1s = 114.8 minutes = ~1.9 hours**, measured directly on the server. This landed right at the upper-middle of the sandbox-extrapolated 30–120 minute bracket — the extrapolation held up despite projecting well beyond directly-tested chromosome lengths.

Since the 11 chromosomes are independent, running them all simultaneously brings the whole ground-truth build down to ~115 minutes wall-clock (not ~21 hours, which is what serial execution would cost). **Before launching all 11 at once, check the memory footprint of that single completed run** (peak RSS, e.g. from `top`/`htop` while it was running, or rerun with `/usr/bin/time -v`) — 11 simultaneous jobs need roughly 11x that single job's memory, and that number isn't yet in hand. If headroom is tight, reduce concurrency in `stage1_build_groundtruth.py` rather than letting it default to 11 at once.

With a 100-thread machine, there's also substantial spare capacity in that same ~115-minute window — up to ~9 independent replicate "worlds" (99 of 100 threads) could run alongside the main build, if replicate-world robustness checks are worth having now rather than later (see NEXT_PHASE_PLANNING.md Section 2).

## 3. Run the production build

```bash
nohup python3 stage1_build_groundtruth.py > build.log 2>&1 &
# or, if tmux/screen is available, prefer that over nohup so you can reattach:
tmux new -d -s groundtruth 'python3 stage1_build_groundtruth.py'
```

`REPLICATE_WORLDS` defaults to 10 -- it processes **one world at a time** (11 chromosomes in parallel within a world, ~115 min, then moves to the next world), not all worlds simultaneously. At 500GB RAM this distinction barely matters for safety, but it does mean predictable, steady progress over however many days it takes rather than one huge simultaneous batch. It's also resumable: rerunning the script skips any `world{N}_chr{C}.trees` that already exists, so it's safe to interrupt and restart (e.g. after a disconnect, if not already running under nohup/tmux).

Check progress any time with:
```bash
ls groundtruth/ | wc -l        # how many chromosome files exist so far
tail -f build.log              # if running under nohup
```

Output: `groundtruth/world00_chr00.trees` ... `world09_chr10.trees` (110 files total if all 10 worlds complete).

## 4. Validate the build

```bash
python3 stage1_validate_groundtruth.py
```

Checks the actual produced output against all three locked targets (Fst=0.05, implied Ne~1e5, LD half-decay~711bp) with PASS/WARN/FAIL bands. **Don't skip this** -- it's what caught the ancestral-Ne bug during development (Section 5d of PROJECT_SUMMARY.md); there's no guarantee a different server/msprime version won't surface something else.

If anything reads FAIL: don't proceed to Stage 2 until it's understood. WARN is worth a look but not necessarily a blocker -- the bands are heuristic (\u00b110%/\u00b125%), not hard pass/fail thresholds from theory.

## 5. What "next step" means concretely from here

Once the ground truth validates cleanly:

1. **Start Stage 2 (collection-design sampling)** -- this is genuinely new code, not yet written. For each (site count x mothers-per-site x total N) design point, sample the corresponding individuals from the ground-truth tree sequences and generate their genotypes (mutation overlay + `ts.simplify()` down to the sampled individuals). This is where the 100 threads actually get used heavily -- see NEXT_PHASE_PLANNING.md Section 2 on the parallelism ceiling.
2. **Re-run the LD-decay check on the full structured metapopulation**, not just the isolated single-deme test used for calibration (flagged as open in PROJECT_SUMMARY.md Section 7) -- cheap to do alongside Stage 2's first sampling runs, since you'll already be pulling genotypes from the ground truth.
3. **Hold off on Stage 3 (SLiM) and the container** until Stage 2 has produced sampled genotypes to hand off -- no point setting up SLiM before there's anything for it to receive.

## 6. If full-length chromosomes turn out to be too slow even on the server

Fallback options, roughly in order of preference:
- Reduce `CHROM_LENGTH` in `stage1_build_groundtruth.py` (e.g. to 2Mb) -- shrinks ROH-detection ceiling but is a known, already-discussed tradeoff (see NEXT_PHASE_PLANNING.md Section 1's original 11x4.5Mb vs 22x2.25Mb discussion).
- Reduce `n_per_deme` for the ground truth itself (this doesn't need to be large -- it's the *source* population, not a collection design; Stage 2 resamples from it).
- Revisit `ancestral_merge_time` -- shallower values (5,000-10,000 generations) gave no real speedup in sandbox testing, but that was at small scale; worth a quick re-check at full length if the deeper merge proves to be the bottleneck specifically.
