"""
stage1/validate_groundtruth.py

Checks a completed ground truth (built by stage1/build_groundtruth.py)
against all locked calibration targets. Run this after every build --
it caught a real bug (wrong ancestral Ne) during development that no
component-level calibration test could have surfaced.

Checks:
  1. Fst           -- target 0.05 (branch-mode, no mutations required)
  2. Implied Ne    -- target ~1e5 (free byproduct of the Fst calculation)
  3. LD half-decay -- target ~711bp (requires mutations overlaid on top)

Usage: python3 stage1/validate_groundtruth.py [world_idx]
  world_idx defaults to 0.
"""

import glob
import os
import sys

import msprime
import numpy as np
import tskit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params
from stage1.utils import branch_gst
from stage1.calibration.calibrate_recombination import (
    pairwise_r2_by_distance, bin_pairs
)

DATADIR = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")
TARGETS = {"fst": 0.05, "ne_global": 1e5, "ld_half_decay_bp": 711}


def status(value, target, rel_tol_pass=0.10, rel_tol_warn=0.25):
    rel_dev = abs(value - target) / target
    if rel_dev <= rel_tol_pass:
        return "PASS"
    elif rel_dev <= rel_tol_warn:
        return "WARN"
    return "FAIL"


def check_fst_and_ne(chrom_paths):
    Ht_list, fst_list = [], []
    all_samples = list(range(params.K * params.N_PER_DEME))
    sample_sets = [
        list(range(i * params.N_PER_DEME, (i + 1) * params.N_PER_DEME))
        for i in range(params.K)
    ]
    for p in chrom_paths:
        ts = tskit.load(p)
        Ht = ts.diversity(sample_sets=[all_samples], mode="branch")[0]
        Hs = ts.diversity(sample_sets=sample_sets, mode="branch").mean()
        Ht_list.append(Ht)
        fst_list.append((Ht - Hs) / Ht)
    return float(np.mean(fst_list)), float(np.mean(Ht_list)) / 2, fst_list


def check_ld_decay(chrom_path,
                    bin_edges=(0, 200, 450, 711, 1000, 1500, 2200, 3000)):
    """LD decay within deme0 only -- isolates the within-deme process that
    r was calibrated against. Note: deme0 in the real metapopulation still
    receives gene flow (m=0.0042) from its neighbour, so the measured
    half-decay distance will be slightly longer than the isolated-deme
    calibration target due to admixture-driven LD. This is expected and
    understood -- see docs/PROJECT_SUMMARY.md Section 5f."""
    ts = tskit.load(chrom_path)
    ts = msprime.sim_mutations(ts, rate=params.MU, random_seed=12345)
    geno = ts.genotype_matrix()
    positions = np.array([s.position for s in ts.sites()])

    geno_d0 = geno[:, :params.N_PER_DEME]
    counts = geno_d0.sum(axis=1)
    poly = (counts > 0) & (counts < params.N_PER_DEME)
    n_dropped = (~poly).sum()
    if n_dropped > 0:
        print(f"  (dropped {n_dropped}/{len(poly)} sites monomorphic within deme0)")
    geno_d0 = geno_d0[poly]
    positions = positions[poly]

    pairs = pairwise_r2_by_distance(geno_d0, positions, max_dist=bin_edges[-1])
    return bin_pairs(pairs, bin_edges)


if __name__ == "__main__":
    world_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    pattern = os.path.join(DATADIR, f"world{world_idx:02d}_chr*.trees")
    chrom_paths = sorted(glob.glob(pattern))

    if not chrom_paths:
        sys.exit(f"No tree sequences found at {pattern} -- "
                 f"run stage1/build_groundtruth.py first.")
    print(f"Found {len(chrom_paths)} chromosome(s) for world {world_idx:02d} "
          f"(n_per_deme={params.N_PER_DEME})\n")

    print("=" * 60)
    print("CHECK 1: Fst (target 0.05)")
    print("=" * 60)
    fst_mean, ne_implied, fst_list = check_fst_and_ne(chrom_paths)
    s = status(fst_mean, TARGETS["fst"])
    print(f"[{s}] mean Fst = {fst_mean:.4f} (target {TARGETS['fst']})")
    print(f"      per-chromosome: [{min(fst_list):.3f}, {max(fst_list):.3f}]")

    print("\n" + "=" * 60)
    print("CHECK 2: implied global Ne (target ~1e5)")
    print("=" * 60)
    s = status(ne_implied, TARGETS["ne_global"])
    print(f"[{s}] implied Ne_global = {ne_implied:,.0f} (target {TARGETS['ne_global']:,.0f})")

    print("\n" + "=" * 60)
    print("CHECK 3: LD half-decay distance, deme0 (target ~711bp)")
    print("=" * 60)
    rows = check_ld_decay(chrom_paths[0])
    near_zero = rows[0][3]
    for lo, hi, n, r2 in rows:
        marker = " <-- 711bp" if lo <= 711 < hi else ""
        print(f"  {lo:>5}-{hi:<5} bp | n={n:6d} | mean R2={r2:.4f}{marker}")
    target_bin = next((row for row in rows if row[0] <= 711 < row[1]), None)
    if target_bin and not np.isnan(target_bin[3]):
        ratio = target_bin[3] / near_zero
        s = status(ratio, 0.5, rel_tol_pass=0.15, rel_tol_warn=0.35)
        print(f"[{s}] ratio at 711bp bin = {ratio:.3f} (target 0.5)")
        if s == "WARN":
            print("  NOTE: WARN here is expected and understood -- admixture-driven")
            print("  LD from neighbouring demes adds a slow-decaying floor to the")
            print("  within-deme curve. See docs/PROJECT_SUMMARY.md Section 5f.")

    print(f"\nPASS/WARN bands: ±10% / ±25% relative deviation from target.")
