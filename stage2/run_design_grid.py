"""
stage2/run_design_grid.py

Batch runner for the Stage 2 fixed-N design grid sweep.

MEMORY DESIGN
-------------
The naive approach (sample_design() builds all 11 chromosome geno_seeds
at once) would use 12-24 GB per job at full scale -- far too much with
88 workers. The worker here processes ONE CHROMOSOME AT A TIME, accumulates
statistics, then frees the array before loading the next chromosome.
Peak memory per job = one chromosome's geno_seeds:
  N=4096 designs: ~1.1 GB/job  → 88 workers = ~97 GB total
  N=8192 designs: ~2.2 GB/job  → 88 workers = ~194 GB total
Both are comfortably within the 500 GB server limit.

DESIGN TABLE
------------
Two total-N values (4096 and 8192), 9 allocation combinations each,
2 spatial strategies, 3 RNG seeds, 10 worlds = 1,080 total jobs.

RESUMABILITY
------------
Safe to interrupt and restart -- completed rows are tracked in the CSV.

Usage:
    nohup python3 stage2/run_design_grid.py > grid.log 2>&1 &
    tail -f grid.log
"""

import csv
import gc
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import tskit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params

from stage2.sampling import (
    CollectionDesign, validate_design_feasibility,
    select_sites, select_diploid_mothers, _build_chrom_geno,
)
from stage2.diversity import (
    compute_background, compare_to_background, DiversityResult,
    _stats_from_geno, _roh_from_diploid_geno, _aggregate_roh,
    _haplotype_stats_from_geno,
)

# ---------------------------------------------------------------------------
# DESIGN TABLE
# ---------------------------------------------------------------------------

DESIGNS_1024 = [
    (4,  8,  32), (4,  16, 16),  (4,  32, 8),
    (8,  8,  16),  (8,  16, 8),  (8,  32, 4),
    (16, 8,  8),  (16, 16, 4),  (16, 32, 2),
]  

DESIGNS_2048 = [
    (4,  8,  64), (4,  16, 32),  (4,  32, 16),
    (8,  8,  32),  (8,  16, 16),  (8,  32, 8),
    (16, 8,  16),  (16, 16, 8),  (16, 32, 4),
]

DESIGNS_4096 = [
    (4,  8,  128), (4,  16, 64),  (4,  32, 32),
    (8,  8,  64),  (8,  16, 32),  (8,  32, 16),
    (16, 8,  32),  (16, 16, 16),  (16, 32, 8),
]
DESIGNS_8192 = [
    (4,  8,  256), (4,  16, 128), (4,  32, 64),
    (8,  8,  128), (8,  16, 64),  (8,  32, 32),
    (16, 8,  64),  (16, 16, 32),  (16, 32, 16),
]

SELECTION_STRATEGIES = ["even_spacing", "random"]
POLLEN_POOL        = "local"   # paternal gametes from same deme as mother -- preserves Fst
RNG_SEEDS          = [1, 2, 3]
WORLD_OVERRIDE     = None

ROH_MIN_LENGTH     = 100_000   # bp
HAP_WINDOW_SIZE    = 100_000   # bp

# ---------------------------------------------------------------------------
# PATHS AND PARALLELISM
# ---------------------------------------------------------------------------

GROUNDTRUTH_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")
RESULTS_DIR     = os.path.join(os.path.dirname(__file__), "..", "data", "results")
RESULTS_CSV     = os.path.join(RESULTS_DIR, "design_grid_results.csv")
BACKGROUND_DIR  = os.path.join(os.path.dirname(__file__), "..", "data")

# Safe with streaming: peak per job = 1-2 GB regardless of design size
N_WORKERS = min(os.cpu_count() or 1, 88)

# ---------------------------------------------------------------------------
# Grid enumeration
# ---------------------------------------------------------------------------

def enumerate_designs():
    designs = []
    for sites, mothers, seeds in DESIGNS_1024 + DESIGNS_2048 + DESIGNS_4096 + DESIGNS_8192:
        if sites > params.K or mothers > params.N_PER_DEME // 2:
            print(f"  SKIP ({sites},{mothers},{seeds}): exceeds ground truth capacity")
            continue
        for strategy in SELECTION_STRATEGIES:
            for rng_seed in RNG_SEEDS:
                designs.append(CollectionDesign(
                    n_sites=sites, mothers_per_site=mothers, seeds_per_mother=seeds,
                    site_selection=strategy, pollen_pool=POLLEN_POOL, seed=rng_seed,
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
                done.add((int(row["world_idx"]), int(row["n_sites"]),
                           int(row["mothers_per_site"]) if str(row["mothers_per_site"]).isdigit()
                               else row["mothers_per_site"],
                           int(row["seeds_per_mother"]),
                           row["site_selection"], int(row["seed"])))
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
# Streaming worker: one chromosome at a time, never holds all 11 at once
# ---------------------------------------------------------------------------

def _run_one(args):
    """Process chromosomes sequentially, accumulating statistics as we go.

    This keeps peak memory at one chromosome's geno_seeds array (1-2 GB)
    rather than all eleven (12-24 GB), making it safe to run 88 simultaneously
    on a 500 GB server regardless of design size.
    """
    import msprime as _msprime   # local import: each subprocess gets its own

    design, world_idx, background_path = args
    bg = pickle.load(open(background_path, "rb"))

    validate_design_feasibility(design, K=params.K, n_per_deme=params.N_PER_DEME)

    n_per_deme        = params.N_PER_DEME
    n_diploid_per_deme = n_per_deme // 2
    rng = np.random.default_rng(design.seed)

    site_demes = select_sites(params.K, design.n_sites, design.site_selection, rng)
    mothers_by_site = [
        select_diploid_mothers(n_diploid_per_deme, n_mothers, rng)
        for n_mothers in design.mothers_list()
    ]
    n_seeds   = sum(len(m) for m in mothers_by_site) * design.seeds_per_mother
    n_mothers = sum(len(m) for m in mothers_by_site)

    # per-chromosome accumulators
    seed_pi_list, seed_He_list, seed_roh_per_chrom, seed_hap_list = [], [], [], []
    mat_pi_list,  mat_He_list,  mat_roh_per_chrom,  mat_hap_list  = [], [], [], []
    seed_afs_total = mat_afs_total = None
    seed_n_seg = mat_n_seg = 0
    total_bp = 0.0

    for c in range(params.N_CHROM):
        path = os.path.join(GROUNDTRUTH_DIR,
                            f"world{world_idx:02d}_chr{c:02d}.trees")
        ts = tskit.load(path)
        ts = _msprime.sim_mutations(
            ts, rate=params.MU, random_seed=design.seed + 1000 * c + 1)

        gd = _build_chrom_geno(
            ts, site_demes, mothers_by_site, n_per_deme,
            design.seeds_per_mother, design.pollen_pool, params.R, rng,
        )

        chrom_len  = gd["seq_length"]
        positions  = gd["positions"]
        total_bp  += chrom_len

        if gd["geno_seeds"].shape[0] > 0:
            # seeds statistics
            g = gd["geno_seeds"]
            st = _stats_from_geno(g, chrom_len)
            seed_pi_list.append(st["pi"])
            seed_He_list.append(st["He"])
            seed_n_seg  += st["n_seg"]
            seed_afs_total = (st["afs"] if seed_afs_total is None
                               else seed_afs_total + st["afs"])
            seed_roh_per_chrom.append(
                _roh_from_diploid_geno(g, positions, chrom_len,
                                        ROH_MIN_LENGTH, n_seeds))
            seed_hap_list.append(_haplotype_stats_from_geno(
                g, positions, chrom_len, HAP_WINDOW_SIZE,
                n_seeds, use_pairs=True))

            # maternal statistics (same geno_maternal, different columns)
            gm = gd["geno_maternal"]
            stm = _stats_from_geno(gm, chrom_len)
            mat_pi_list.append(stm["pi"])
            mat_He_list.append(stm["He"])
            mat_n_seg  += stm["n_seg"]
            mat_afs_total = (stm["afs"] if mat_afs_total is None
                              else mat_afs_total + stm["afs"])
            mat_roh_per_chrom.append(
                _roh_from_diploid_geno(gm, positions, chrom_len,
                                        ROH_MIN_LENGTH, n_mothers))
            mat_hap_list.append(_haplotype_stats_from_geno(
                gm, positions, chrom_len, HAP_WINDOW_SIZE,
                n_mothers, use_pairs=False))

        # *** free memory before next chromosome ***
        del ts, gd
        gc.collect()

    # aggregate seeds
    seed_roh_stats = _aggregate_roh(seed_roh_per_chrom, total_bp)
    seed_result = DiversityResult(
        n_individuals=n_seeds, n_chromosomes=params.N_CHROM, total_bp=total_bp,
        pi_mean=float(np.mean(seed_pi_list)) if seed_pi_list else 0.0,
        pi_per_chrom=seed_pi_list,
        n_seg_sites_total=seed_n_seg,
        seg_sites_per_bp=seed_n_seg / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(seed_He_list)) if seed_He_list else 0.0,
        He_per_chrom=seed_He_list,
        afs=seed_afs_total.tolist() if seed_afs_total is not None else [],
        roh_min_length_bp=ROH_MIN_LENGTH, **seed_roh_stats,
    )
    seed_result = compare_to_background(seed_result, bg)

    # aggregate maternal
    mat_roh_stats = _aggregate_roh(mat_roh_per_chrom, total_bp)
    mat_result = DiversityResult(
        n_individuals=n_mothers, n_chromosomes=params.N_CHROM, total_bp=total_bp,
        pi_mean=float(np.mean(mat_pi_list)) if mat_pi_list else 0.0,
        pi_per_chrom=mat_pi_list,
        n_seg_sites_total=mat_n_seg,
        seg_sites_per_bp=mat_n_seg / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(mat_He_list)) if mat_He_list else 0.0,
        He_per_chrom=mat_He_list,
        afs=mat_afs_total.tolist() if mat_afs_total is not None else [],
        roh_min_length_bp=ROH_MIN_LENGTH, **mat_roh_stats,
    )
    mat_result = compare_to_background(mat_result, bg)

    # haplotype diversity
    def _avg(lst, key):
        vals = [s[key] for s in lst if key in s]
        return float(np.mean(vals)) if vals else 0.0

    hap_dict = {}
    for prefix, lst in [("hap_seed_", seed_hap_list), ("hap_mat_", mat_hap_list)]:
        for k in ("n_distinct_mean", "h_mean", "n_eff_hap_mean", "frac_windows_novar"):
            hap_dict[f"{prefix}{k}"] = _avg(lst, k)
    hap_dict["hap_window_size_bp"] = HAP_WINDOW_SIZE

    out = ({"seed_" + k: v for k, v in seed_result.to_dict().items()} |
           {"mat_"  + k: v for k, v in mat_result.to_dict().items()}  |
           hap_dict)
    out |= {
        "world_idx":        world_idx,
        "n_sites":          design.n_sites,
        "mothers_per_site": design.mothers_per_site
                            if isinstance(design.mothers_per_site, int)
                            else "unequal",
        "seeds_per_mother": design.seeds_per_mother,
        "n_mothers":        n_mothers,
        "total_seeds":      design.total_n(),
        "site_selection":   design.site_selection,
        "pollen_pool":      design.pollen_pool,
        "seed":             design.seed,
    }
    return out


# ---------------------------------------------------------------------------
# Background
# ---------------------------------------------------------------------------

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

    print(f"Design table: {len(DESIGNS_1024)} N=1024 combos + {len(DESIGNS_2048)} N=2048 combos + {len(DESIGNS_4096)} N=4096 combos + "
          f"{len(DESIGNS_8192)} N=8192 combos")
    print(f"× {len(SELECTION_STRATEGIES)} strategies × {len(RNG_SEEDS)} RNG seeds "
          f"= {len(designs)} designs")
    print(f"× {len(worlds)} world(s) = {len(designs)*len(worlds)} total jobs")
    print(f"Workers: {N_WORKERS}  |  Memory model: streaming (1 chrom at a time)")
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
            print(f"world {world_idx:02d}: all {len(designs)} done, skipping")
            continue

        print(f"world {world_idx:02d}: {len(todo)}/{len(designs)} jobs to run...")
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
                          f"s{d.n_sites}×m{d.mothers_per_site}×"
                          f"spm{d.seeds_per_mother}: {e}")

        if new_rows:
            write_rows(new_rows, RESULTS_CSV, write_header=first_write)
            first_write = False

        elapsed = time.time() - t0
        print(f"world {world_idx:02d}: done in {elapsed:.1f}s "
              f"({elapsed/max(len(todo),1):.1f}s/job avg) -- "
              f"total {total_done}/{total_jobs}")

    print(f"\nAll done. Results: {RESULTS_CSV}")
    print("Load in R: df <- read.csv('data/results/design_grid_results.csv')")

