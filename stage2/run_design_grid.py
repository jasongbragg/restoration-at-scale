"""
stage2/run_design_grid.py

Batch runner for the Stage 2 collection-design grid sweep.

Iterates over every (design, world, rng_seed) combination, runs each through
stage2/sampling.sample_design() and stage2/diversity.design_stats_dict(),
and writes results to a flat CSV. Each row is one fully-evaluated combination
-- one row per (design × world × rng_seed) -- ready to load in R for the
accumulation curves and decision table.

OUTPUT COLUMNS
--------------
The CSV column list is derived dynamically from the first result returned by
design_stats_dict() -- not hardcoded here. This keeps the runner in sync
automatically when diversity metrics are added or renamed, without needing
a matching edit in two places.

Key column groups:
  Provenance:      world_idx, n_sites, mothers_per_site, seeds_per_mother,
                   n_mothers, total_seeds, site_selection, pollen_pool, seed
  Seed stats:      seed_pi_mean, seed_He_mean, seed_n_seg_sites_total, ...
                   seed_roh_mean_n, seed_roh_fraction_genome, ...
                   seed_pi_fraction_of_background, ...
  Maternal stats:  mat_pi_mean, mat_He_mean, mat_n_seg_sites_total, ...
                   mat_pi_fraction_of_background, ...
                   (no mat_roh_* -- not meaningful for haploid maternal cohort)

NOTE ON N
---------
  n_mothers    = n_sites × mothers_per_site  (maternal lines sampled / genotyped)
  total_seeds  = n_mothers × seeds_per_mother (individuals planted in restoration)

RESUMABILITY
------------
Completed (world, n_sites, mothers_per_site, seeds_per_mother,
           site_selection, rng_seed) tuples are tracked in the CSV.
Safe to interrupt and restart -- completed rows are skipped.

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
from itertools import product

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params
from stage2.sampling import CollectionDesign, sample_design
from stage2.diversity import compute_background, design_stats_dict

# ---------------------------------------------------------------------------
# GRID DEFINITION -- edit these to change the sweep
# ---------------------------------------------------------------------------

SITES_OPTIONS      = [2, 4, 8, 16]           # n_sites values
MOTHERS_OPTIONS    = [2, 4, 8, 16, 32]        # mothers_per_site values
SELECTION_STRATEGIES = ["even_spacing", "random"]
SEEDS_PER_MOTHER   = 5                        # seeds retained per maternal line;
                                               # total planted = mothers × this
POLLEN_POOL        = "metapopulation"          # pollen drawn from all K demes
                                               # (realistic for wind-pollinated Melaleuca)

# Multiple RNG seeds per (design, world) give variance from site/mother
# selection, independent of the world-to-world ancestral history variance.
RNG_SEEDS = [1, 2, 3]

# Which worlds to process. None = all complete worlds found automatically.
# Set e.g. WORLD_OVERRIDE = [0] to run a quick test on world 00 only.
WORLD_OVERRIDE = None

# ---------------------------------------------------------------------------
# PATHS AND PARALLELISM
# ---------------------------------------------------------------------------

GROUNDTRUTH_DIR    = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")
RESULTS_DIR        = os.path.join(os.path.dirname(__file__), "..", "data", "results")
RESULTS_CSV        = os.path.join(RESULTS_DIR, "design_grid_results.csv")
BACKGROUND_DIR     = os.path.join(os.path.dirname(__file__), "..", "data")

# Each job loads 11 .trees files (~7 MB each = ~77 MB) and runs
# sampling + diversity. Disk I/O is the practical limit at high concurrency.
# Tune down if you see disk saturation (e.g. set to 44).
N_WORKERS = min(os.cpu_count() or 1, 88)

# ---------------------------------------------------------------------------
# Grid enumeration
# ---------------------------------------------------------------------------

def enumerate_designs():
    """All (CollectionDesign) instances within the ground-truth capacity."""
    designs = []
    skipped = 0
    for sites, mothers, strategy, rng_seed in product(
        SITES_OPTIONS, MOTHERS_OPTIONS, SELECTION_STRATEGIES, RNG_SEEDS
    ):
        if sites > params.K or mothers > params.N_PER_DEME:
            skipped += 1
            continue
        designs.append(CollectionDesign(
            n_sites=sites,
            mothers_per_site=mothers,
            seeds_per_mother=SEEDS_PER_MOTHER,
            site_selection=strategy,
            pollen_pool=POLLEN_POOL,
            seed=rng_seed,
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
        if all(os.path.exists(
                os.path.join(GROUNDTRUTH_DIR, f"world{w:02d}_chr{c:02d}.trees"))
               for c in range(params.N_CHROM)):
            worlds.append(w)
    return worlds


# ---------------------------------------------------------------------------
# Completion tracking
# ---------------------------------------------------------------------------

def _completion_key(row):
    """Tuple used to identify a completed job in the CSV."""
    return (
        int(row["world_idx"]),
        int(row["n_sites"]),
        int(row["mothers_per_site"]) if str(row["mothers_per_site"]).isdigit()
            else row["mothers_per_site"],
        int(row["seeds_per_mother"]),
        row["site_selection"],
        int(row["seed"]),
    )


def _design_key(d, world_idx):
    return (
        world_idx,
        d.n_sites,
        d.mothers_per_site if isinstance(d.mothers_per_site, int) else "unequal",
        d.seeds_per_mother,
        d.site_selection,
        d.seed,
    )


def load_completed(csv_path):
    done = set()
    if not os.path.exists(csv_path):
        return done
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                done.add(_completion_key(row))
            except (KeyError, ValueError):
                pass   # malformed row -- skip
    return done


# ---------------------------------------------------------------------------
# CSV writing -- dynamic columns from first result
# ---------------------------------------------------------------------------

_COLUMNS = None   # set from first result; None until then


def write_rows(rows, csv_path, write_header):
    global _COLUMNS
    if not rows:
        return
    if _COLUMNS is None:
        # derive column order from the first result; provenance columns first
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
        sys.exit(f"No complete worlds found in {GROUNDTRUTH_DIR}. "
                 f"Run stage1/build_groundtruth.py first.")

    print(f"Grid: {len(designs)} designs × {len(worlds)} world(s) "
          f"× 1 (seeds_per_mother={SEEDS_PER_MOTHER})")
    print(f"Worlds: {worlds}")
    print(f"Total jobs: {len(designs) * len(worlds)}  |  Workers: {N_WORKERS}")
    print(f"Columns: derived dynamically from first result\n")

    completed  = load_completed(RESULTS_CSV)
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
                          f"s{d.n_sites}m{d.mothers_per_site}: {e}")

        if new_rows:
            write_rows(new_rows, RESULTS_CSV, write_header=first_write)
            first_write = False

        elapsed = time.time() - t0
        print(f"world {world_idx:02d}: done in {elapsed:.1f}s "
              f"({elapsed/max(len(todo),1):.1f}s/job avg) -- "
              f"total {total_done}/{total_jobs}")

    print(f"\nAll done. Results at: {RESULTS_CSV}")
    print("Load in R: df <- read.csv('data/results/design_grid_results.csv')")
