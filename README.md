# Restoration at scale: simulation framework

Simulation pipeline for investigating how seed-collection design
(sites visited × maternal trees sampled per site) affects downstream
genetic diversity outcomes in ecological restoration plantings — at
scales (1,000s–10,000s of seedlings) beyond what prior small-*n*
empirical studies have addressed.

**Biological context:** calibrated against *Melaleuca quinquenervia*
(wetland paperbark, NSW Australia), using population-genetic parameters
from Guo et al. (2026) *Molecular Ecology* 35:e70413.

---

## Four-stage pipeline

| Stage | What it does | Status |
|---|---|---|
| 1 | Coalescent ancestry (msprime) — build a weak-IBD metapopulation ground truth | Complete, validated |
| 2 | Collection-design sampling — draw sites × mothers from the ground truth | Code written, design grid pending |
| 3 | Forward restoration simulation (SLiM) — a few generations post-planting | Not yet started |
| 4 | Diversity outcome analysis — accumulation curves and decision table | Not yet started |

## Quick start

```bash
pip install msprime tskit numpy

# Build the Stage 1 ground truth (runs for days, one world at a time):
nohup python3 stage1/build_groundtruth.py > build.log 2>&1 &

# Validate a completed world:
python3 stage1/validate_groundtruth.py

# Sample a collection design (once ground truth exists):
python3 stage2/demo.py
```

## Parameters

All locked simulation parameters live in **`params.py`** at the repo root.
Change values there and nowhere else — per-script constants caused real
bugs during development (see `docs/PROJECT_SUMMARY.md` Section 5).

## Data

`data/` holds the simulation outputs (tree sequences, sampled designs).
It is **not committed to git** (files are large). Back it up separately
(S3, NAS, or similar). The scripts are deterministic given the parameters
in `params.py` and their random seeds, so outputs can always be
reproduced from code alone.

## Documentation

- `docs/PROJECT_SUMMARY.md` — full calibration history, parameter derivations, design decisions, methodological lessons
- `docs/NEXT_PHASE_PLANNING.md` — forward planning notes (genome size, parallelism, containers)
- `docs/RUNBOOK.md` — concrete server execution steps

## Working with Claude

This repository is the shared workspace. Claude can read files directly
from GitHub (via the raw URL or API) when given the repo URL. Workflow:
1. Make code changes locally, push to `main`
2. Share the repo URL or specific file URLs with Claude
3. Claude fetches the current code, proposes changes as files for download
4. Download, review, add to repo, push
