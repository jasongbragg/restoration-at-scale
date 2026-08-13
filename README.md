# Restoration at scale: simulation framework

Simulation pipeline for examining the genetic diversity of restoration plantings established using different collection strategies. A large population is simulated with a coalescent model, and mother trees and seedlings are 'sampled' from it for large plantings (e.g. 1e3, 1e4 seedlings). The goal is to understand the diversity that is obtained using different collection strategies, particularly where large plantings are being contemplated. 

The background population is modelled loosely on *Melaleuca quinquenervia*
calibrated using paramaters (Ne, Fst, LD decay) from a study by Guo et al. (2026) *Molecular Ecology* 35:e70413.

---

## The project will ultimately have five stages

| Stage | What it does | Status |
|---|---|---|
| 1 | Coalescent simulation (msprime) of a population across a landscape, with weak-IBD, and rapid LD decay | Complete, validated |
| 2 | Collection-design sampling — draw sites × mothers from the ground truth | Code written, design grid pending |
| 3 | Forward restoration simulation (SLiM) — a few generations post-planting | Not yet started |
| 4 | Diversity outcome analysis — accumulation curves and decision table | Not yet started |
| 5 | Response to selection, including positive selection, purging | Not yet started |


## Parameters

All locked simulation parameters live in **`params.py`** at the repo root.
This includes chromosome length and number. Chromosomes are small for computational tractability. 

## Data

`data/` holds the simulation outputs (tree sequences, sampled designs). 
The scripts are deterministic given the parameters in `params.py` and 
their random seeds, so outputs can always be reproduced from code alone.

