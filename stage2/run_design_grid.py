"""
stage2/run_design_grid.py

Batch runner for the Stage 2 collection-design grid sweep.

Iterates over every (design, world, seed) combination, runs each through
stage2/sampling.sample_design() and stage2/diversity.design_stats_dict(),
and writes results to a flat CSV. Each row is one fully-evaluated
combination -- ready to load directly in R for the accumulation curves
and decision table.

DESIGN GRID
-----------
The grid is defined by SITES_OPTIONS, MOTHERS_OPTIONS, and
SELECTION_STRATEGIES below. Edit those lists to change the sweep.

NOTE ON N: the total N here is the number of *maternal lines* sampled
(genotyped mother trees), not the number of seedlings planted. Each
maternal line can contribute many seeds to a planting. The ground truth
ceiling is params.K × params.N_PER_DEME = 16 × 32 = 512 maternal lines.

PARALLELISM
-----------
Processes one world at a time, all designs for that world in parallel
via ProcessPoolExecutor. This caps peak memory at n_workers simultaneous
jobs, each loading ~77MB of .trees files (11 chromosomes × ~7MB each).
Adjust N_WORKERS if disk I/O becomes a bottleneck.

RESUMABILITY
------------
Completed (design, world, seed) combinations are logged to the CSV as
they finish. On restart, completed rows are skipped. Safe to interrupt
and resume.

Usage:
    python3 stage2/run_design_grid.py

Output:
    data/results/design_grid_results.csv
"""

import csv
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params
from stage2.sampling import CollectionDesign, sample_design
from stage2.diversity import compute_background, design_stats_dict

# ---------------------------------------------------------------------------
# GRID DEFINITION -- edit these to change the sweep
# ---------------------------------------------------------------------------

SITES_OPTIONS = [2, 4, 8, 16]          # n_sites values
MOTHERS_OPTIONS = [2, 4, 8, 16, 32]    # mothers_per_site values
SELECTION_STRATEGIES = ["even_spacing", "random"]

# Multiple seeds per (design, world) give sampling variance independent
# of the ancestral history variance (which comes from multiple worlds).
SEEDS_PER_DESIGN = [1, 2, 3]

# Which worlds to include. Defaults to all available in data/groundtruth/.
# Override here if you want to limit the sweep, e.g. WORLD_OVERRIDE = [0, 1]
WORLD_OVERRIDE = None

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

GROUNDTRUTH_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "results")
RESULTS_CSV = os.path.join(RESULTS_DIR, "design_grid_results.csv")
BACKGROUND_CACHE_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

# cap at available CPUs; each job loads 11 .trees files so disk I/O can
# be a bottleneck at very high concurrency -- reduce if needed
N_WORKERS = min(os.cpu_count() or 1, 44)

# ---------------------------------------------------------------------------
# Grid enumeration and feasibility
# ---------------------------------------------------------------------------

def enumerate_designs():
    """All (CollectionDesign, seed) pairs within the ground-truth capacity."""
    designs = []
    skipped = 0
    for sites, mothers, strategy, seed in product(
        SITES_OPTIONS, MOTHERS_OPTIONS, SELECTION_STRATEGIES, SEEDS_PER_DESIGN
    ):
        if sites > params.K:
            skipped += 1
            continue
        if mothers > params.N_PER_DEME:
            skipped += 1
            continue
        designs.append(CollectionDesign(
            n_sites=sites, mothers_per_site=mothers,
            site_selection=strategy, seed=seed,
        ))
    if skipped:
        print(f"  ({skipped} combinations skipped: exceed K={params.K} or "
              f"N_PER_DEME={params.N_PER_DEME})")
    return designs


def available_worlds():
    """Worlds with all N_CHROM chromosomes present in groundtruth_dir."""
    if WORLD_OVERRIDE is not None:
        return WORLD_OVERRIDE
    worlds = []
    for w in range(params.REPLICATE_WORLDS):
        paths = [
            os.path.join(GROUNDTRUTH_DIR, f"world{w:02d}_chr{c:02d}.trees")
            for c in range(params.N_CHROM)
        ]
        if all(os.path.exists(p) for p in paths):
            worlds.append(w)
    return worlds


# ---------------------------------------------------------------------------
# CSV handling
# ---------------------------------------------------------------------------

# Column order in the output CSV
COLUMNS = [
    "world_idx", "n_sites", "mothers_per_site", "site_selection", "seed",
    "total_n",
    # diversity metrics
    "pi_mean", "He_mean", "n_seg_sites_total", "seg_sites_per_bp",
    # background comparison
    "pi_fraction_of_background", "He_fraction_of_background",
    "seg_sites_fraction_of_background",
    # ROH
    "roh_min_length_bp", "roh_mean_n", "roh_mean_total_length",
    "roh_mean_longest", "roh_fraction_genome",
    # misc
    "n_haploid_samples", "n_diploid_individuals",
]


def load_completed(csv_path):
    """Return a set of (world_idx, n_sites, mothers_per_site, site_selection, seed)
    tuples already present in the CSV."""
    done = set()
    if not os.path.exists(csv_path):
        return done
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            done.add((
                int(row["world_idx"]),
                int(row["n_sites"]),
                int(row["mothers_per_site"])
                  if row["mothers_per_site"].isdigit() else row["mothers_per_site"],
                row["site_selection"],
                int(row["seed"]),
            ))
    return done


def write_rows(rows, csv_path, write_header):
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in COLUMNS})


# ---------------------------------------------------------------------------
# Worker function (runs in a subprocess)
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
    stats = design_stats_dict(result, background=bg)
    return stats


def _background_path(world_idx):
    return os.path.join(BACKGROUND_CACHE_DIR, f"background_world{world_idx:02d}.pkl")


def ensure_background(world_idx):
    """Load cached background or compute and cache it."""
    path = _background_path(world_idx)
    if os.path.exists(path):
        print(f"  background world {world_idx:02d}: loaded from cache")
        return path
    print(f"  background world {world_idx:02d}: computing (this takes a few minutes)...")
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

    worlds = available_worlds()
    designs = enumerate_designs()

    if not worlds:
        sys.exit(f"No complete worlds found in {GROUNDTRUTH_DIR}. "
                 f"Run stage1/build_groundtruth.py first.")

    print(f"Worlds available: {worlds}")
    print(f"Designs in grid: {len(designs)}")
    print(f"Seeds per design: {SEEDS_PER_DESIGN}")
    print(f"Total jobs: {len(worlds) * len(designs)}")
    print(f"Workers: {N_WORKERS}\n")

    completed = load_completed(RESULTS_CSV)
    print(f"Already completed: {len(completed)} jobs (will be skipped)\n")

    first_write = not os.path.exists(RESULTS_CSV)
    total_done = len(completed)
    total_jobs = len(worlds) * len(designs)

    for world_idx in worlds:
        bg_path = ensure_background(world_idx)

        # filter to jobs not yet done for this world
        todo = [
            d for d in designs
            if (world_idx,
                d.n_sites,
                d.mothers_per_site if isinstance(d.mothers_per_site, int)
                  else "unequal",
                d.site_selection,
                d.seed) not in completed
        ]

        if not todo:
            print(f"world {world_idx:02d}: all {len(designs)} jobs already done, skipping")
            continue

        print(f"world {world_idx:02d}: {len(todo)}/{len(designs)} jobs to run "
              f"({N_WORKERS} workers)...")
        t0 = time.time()
        args_list = [(d, world_idx, bg_path) for d in todo]
        new_rows = []

        with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
            futures = {pool.submit(_run_one, a): a for a in args_list}
            for fut in as_completed(futures):
                try:
                    row = fut.result()
                    new_rows.append(row)
                    total_done += 1
                    if len(new_rows) % 10 == 0:
                        # flush to disk every 10 rows so progress is saved
                        write_rows(new_rows, RESULTS_CSV, write_header=first_write)
                        first_write = False
                        new_rows = []
                except Exception as e:
                    d, w, _ = futures[fut]
                    print(f"  ERROR world{w:02d} "
                          f"s{d.n_sites}m{d.mothers_per_site}: {e}")

        # flush any remaining rows
        if new_rows:
            write_rows(new_rows, RESULTS_CSV, write_header=first_write)
            first_write = False

        elapsed = time.time() - t0
        print(f"world {world_idx:02d}: done in {elapsed:.1f}s "
              f"({elapsed/len(todo):.1f}s/job avg). "
              f"Total: {total_done}/{total_jobs}")

    print(f"\nAll done. Results at: {RESULTS_CSV}")
    print("Next: load in R with read.csv() and build accumulation curves.")
