"""
stage2/run_design_grid.py

Batch runner for the Stage 2 fixed-N design grid sweep.

Iterates over every (design, world, rng_seed) combination defined in
DESIGN_TABLE below. Each row fixes total planted N while varying the
allocation across sites, mothers per site, and seeds per mother --
isolating the pure allocation question from the total-effort question.

DESIGN TABLE
------------
Two total N values (4096 and 8192), 9 allocation combinations each,
2 spatial strategies (even_spacing and random), 3 RNG seeds per design,
10 worlds = 1,080 total jobs.

Columns in the output CSV are derived dynamically from the first result
returned by design_stats_dict(). Key groups:
  Provenance : world_idx, n_sites, mothers_per_site, seeds_per_mother,
               n_mothers, total_seeds, site_selection, pollen_pool, seed
  Seed stats : seed_pi_mean, seed_He_mean, seed_n_seg_sites_total,
               seed_roh_mean_n, seed_roh_fraction_genome,
               seed_pi_fraction_of_background, ...
  Maternal   : mat_pi_mean, mat_He_mean, mat_n_seg_sites_total,
               mat_pi_fraction_of_background, ...

NOTE ON N
---------
  n_mothers  = n_sites x mothers_per_site  (trees sampled / genotyped)
  total_seeds = n_mothers x seeds_per_mother  (seedlings planted)
Both are fixed across designs within the same N group.

RESUMABILITY
------------
Safe to interrupt and restart -- completed rows are tracked in the CSV
and skipped on rerun.

Usage:
    nohup python3 stage2/run_design_grid.py > grid.log 2>&1 &
    tail -f grid.log

Output:
    data/results/design_grid_results.csv
"""

import csv
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params
from stage2.sampling import CollectionDesign, sample_design
from stage2.diversity import compute_background, design_stats_dict

# ---------------------------------------------------------------------------
# DESIGN TABLE -- the fixed-N combinations to run
# Each tuple: (sites, mothers_per_site, seeds_per_mother)
# total planted = sites x mothers x seeds
# ---------------------------------------------------------------------------

DESIGNS_4096 = [
    (4,  8,  128),
    (4,  16, 64),
    (4,  32, 32),
    (8,  8,  64),
    (8,  16, 32),
    (8,  32, 16),
    (16, 8,  32),
    (16, 16, 16),
    (16, 32, 8),
]

DESIGNS_8192 = [
    (4,  8,  256),
    (4,  16, 128),
    (4,  32, 64),
    (8,  8,  128),
    (8,  16, 64),
    (8,  32, 32),
    (16, 8,  64),
    (16, 16, 32),
    (16, 32, 16),
]

# ---------------------------------------------------------------------------
# SWEEP PARAMETERS -- edit these to change the sweep
# ---------------------------------------------------------------------------

SELECTION_STRATEGIES = ["even_spacing", "random"]
POLLEN_POOL = "metapopulation"
RNG_SEEDS = [1, 2, 3]

# Which worlds to include. None = all complete worlds found automatically.
WORLD_OVERRIDE = None

# ---------------------------------------------------------------------------
# PATHS AND PARALLELISM
# ---------------------------------------------------------------------------

GROUNDTRUTH_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")
RESULTS_DIR     = os.path.join(os.path.dirname(__file__), "..", "data", "results")
RESULTS_CSV     = os.path.join(RESULTS_DIR, "design_grid_results.csv")
BACKGROUND_DIR  = os.path.join(os.path.dirname(__file__), "..", "data")

N_WORKERS = min(os.cpu_count() or 1, 88)

# ---------------------------------------------------------------------------
# Grid enumeration
# ---------------------------------------------------------------------------

def enumerate_designs():
    """Expand the design table across spatial strategies and RNG seeds."""
    designs = []
    for sites, mothers, seeds in DESIGNS_4096 + DESIGNS_8192:
        if sites > params.K or mothers > params.N_PER_DEME:
            print(f"  SKIP: ({sites}, {mothers}, {seeds}) exceeds ground truth capacity")
            continue
        for strategy in SELECTION_STRATEGIES:
            for rng_seed in RNG_SEEDS:
                designs.append(CollectionDesign(
                    n_sites=sites,
                    mothers_per_site=mothers,
                    seeds_per_mother=seeds,
                    site_selection=strategy,
                    pollen_pool=POLLEN_POOL,
                    seed=rng_seed,
                ))
    return designs


def available_worlds():
    if WORLD_OVERRIDE is not None:
        return WORLD_OVERRIDE
    worlds = []
    for w in range(params.REPLICATE_WORLDS):
        if all(os.path.exists(
                os.path.join(GROUNDTRUTH_DIR, f"world{w:02d}_chr{c:02d}.trees"))
               for c in range(params.N_CHROM)):
            worlds.append(w)
    return worlds


# ---------------------------------------------------------------------------
# Completion tracking
# ---------------------------------------------------------------------------

def _design_key(d, world_idx):
    return (world_idx, d.n_sites,
            d.mothers_per_site if isinstance(d.mothers_per_site, int) else "unequal",
            d.seeds_per_mother, d.site_selection, d.seed)


def load_completed(csv_path):
    done = set()
    if not os.path.exists(csv_path):
        return done
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                done.add((
                    int(row["world_idx"]),
                    int(row["n_sites"]),
                    int(row["mothers_per_site"]) if str(row["mothers_per_site"]).isdigit()
                        else row["mothers_per_site"],
                    int(row["seeds_per_mother"]),
                    row["site_selection"],
                    int(row["seed"]),
                ))
            except (KeyError, ValueError):
                pass
    return done


# ---------------------------------------------------------------------------
# CSV writing -- dynamic columns from first result
# ---------------------------------------------------------------------------

_COLUMNS = None


def write_rows(rows, csv_path, write_header):
    global _COLUMNS
    if not rows:
        return
    if _COLUMNS is None:
        prov = ["world_idx", "n_sites", "mothers_per_site", "seeds_per_mother",
                "n_mothers", "total_seeds", "site_selection", "pollen_pool", "seed"]
        rest = [k for k in rows[0].keys() if k not in prov]
        _COLUMNS = prov + rest
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_COLUMNS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in _COLUMNS})


# ---------------------------------------------------------------------------
# Worker function
# ---------------------------------------------------------------------------

def _run_one(args):
    design, world_idx, background_path = args
    bg = pickle.load(open(background_path, "rb"))
    result = sample_design(
        design, world_idx=world_idx,
        groundtruth_dir=GROUNDTRUTH_DIR,
        K=params.K, n_per_deme=params.N_PER_DEME,
        n_chrom=params.N_CHROM, mu=params.MU,
    )
    return design_stats_dict(result, background=bg)


def _background_path(world_idx):
    return os.path.join(BACKGROUND_DIR, f"background_world{world_idx:02d}.pkl")


def ensure_background(world_idx):
    path = _background_path(world_idx)
    if os.path.exists(path):
        print(f"  background world {world_idx:02d}: loaded from cache")
        return path
    print(f"  background world {world_idx:02d}: computing...")
    t0 = time.time()
    bg = compute_background(GROUNDTRUTH_DIR, world_idx=world_idx)
    pickle.dump(bg, open(path, "wb"))
    print(f"  background world {world_idx:02d}: done in {(time.time()-t0)/60:.1f} min "
          f"(pi={bg.pi_mean:.3e}, He={bg.He_mean:.4f})")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs(RESULTS_DIR, exist_ok=True)

    worlds  = available_worlds()
    designs = enumerate_designs()

    if not worlds:
        sys.exit(f"No complete worlds found in {GROUNDTRUTH_DIR}.")

    n4 = len(DESIGNS_4096) * len(SELECTION_STRATEGIES) * len(RNG_SEEDS)
    n8 = len(DESIGNS_8192) * len(SELECTION_STRATEGIES) * len(RNG_SEEDS)
    print(f"Design table: {len(DESIGNS_4096)} N=4096 combos + "
          f"{len(DESIGNS_8192)} N=8192 combos")
    print(f"× {len(SELECTION_STRATEGIES)} strategies × {len(RNG_SEEDS)} RNG seeds "
          f"= {len(designs)} designs")
    print(f"× {len(worlds)} world(s) = {len(designs)*len(worlds)} total jobs")
    print(f"Workers: {N_WORKERS}")
    print(f"Worlds available: {worlds}\n")

    completed   = load_completed(RESULTS_CSV)
    first_write = not os.path.exists(RESULTS_CSV)
    total_done  = len(completed)
    total_jobs  = len(worlds) * len(designs)
    print(f"Already completed: {total_done} jobs (will be skipped)\n")

    for world_idx in worlds:
        bg_path = ensure_background(world_idx)

        todo = [d for d in designs
                if _design_key(d, world_idx) not in completed]

        if not todo:
            print(f"world {world_idx:02d}: all {len(designs)} jobs done, skipping")
            continue

        print(f"world {world_idx:02d}: {len(todo)}/{len(designs)} jobs to run "
              f"({N_WORKERS} workers)...")
        t0 = time.time()
        new_rows = []

        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(_run_one, (d, world_idx, bg_path)): d
                       for d in todo}
            for fut in as_completed(futures):
                d = futures[fut]
                try:
                    row = fut.result()
                    new_rows.append(row)
                    total_done += 1
                    if len(new_rows) % 10 == 0:
                        write_rows(new_rows, RESULTS_CSV, write_header=first_write)
                        first_write = False
                        new_rows = []
                except Exception as e:
                    print(f"  ERROR world{world_idx:02d} "
                          f"s{d.n_sites}×m{d.mothers_per_site}×spm{d.seeds_per_mother}: "
                          f"{e}")

        if new_rows:
            write_rows(new_rows, RESULTS_CSV, write_header=first_write)
            first_write = False

        elapsed = time.time() - t0
        print(f"world {world_idx:02d}: done in {elapsed:.1f}s "
              f"({elapsed/max(len(todo),1):.1f}s/job avg) -- "
              f"total progress {total_done}/{total_jobs}")

    print(f"\nAll done. Results: {RESULTS_CSV}")
    print("Load in R: df <- read.csv('data/results/design_grid_results.csv')")


