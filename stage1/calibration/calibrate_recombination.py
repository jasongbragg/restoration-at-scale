"""
stage1_calibrate_recombination.py

Calibrates recombination rate (r) against the empirical LD-decay curve
from Guo et al. 2026: R^2 falls to half its near-zero-distance value at
711 bp separation.

Approach: single-deme simulations (deme_Ne=6000, no migration) isolate
the within-deme recombination/drift process that drives SHORT-RANGE LD
decay, deliberately avoiding the population-structure floor that flattens
the real decay curve at long range (that floor is a Fst/migration effect,
already locked in separately -- conflating the two would double-count
structure). SNP density at this Ne/mu is low (~1 per 4167bp), so pairs
are pooled across several replicate chromosomes to get usable counts per
distance bin.
"""

import msprime
import numpy as np


def simulate_single_deme(r, L=2_000_000, Ne=6000, mu=1e-8, n=40, seed=1):
    ts = msprime.sim_ancestry(
        samples=n, population_size=Ne, sequence_length=L,
        recombination_rate=r, random_seed=seed, ploidy=1,
    )
    ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed + 1)
    geno = ts.genotype_matrix()
    positions = np.array([s.position for s in ts.sites()])
    return geno, positions


def pairwise_r2_by_distance(geno, positions, max_dist=3000):
    """All site pairs within max_dist, as (distance, r2) tuples."""
    n_sites = geno.shape[0]
    if n_sites < 2:
        return []
    corr = np.corrcoef(geno)
    r2 = corr ** 2
    pairs = []
    for i in range(n_sites):
        dists = np.abs(positions[i + 1:] - positions[i])
        close = np.where(dists <= max_dist)[0]
        for k in close:
            j = i + 1 + k
            pairs.append((dists[k], r2[i, j]))
    return pairs


def bin_pairs(pairs, bin_edges):
    """Mean r2 and pair count within each [edge[k], edge[k+1]) bin."""
    pairs = np.array(pairs)
    out = []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (pairs[:, 0] >= lo) & (pairs[:, 0] < hi)
        n = mask.sum()
        # nanmean, not mean: a single undefined (zero-variance-site) pair
        # would otherwise poison the whole bin's average via propagation
        mean_r2 = np.nanmean(pairs[mask, 1]) if n > 0 else float("nan")
        out.append((lo, hi, n, mean_r2))
    return out


def ld_decay_curve(r, n_replicates=8, L=2_000_000, Ne=6000, mu=1e-8, n=40,
                    bin_edges=(0, 200, 450, 711, 1000, 1500, 2200, 3000)):
    all_pairs = []
    for rep in range(n_replicates):
        geno, positions = simulate_single_deme(r, L=L, Ne=Ne, mu=mu, n=n,
                                                 seed=1000 * rep + 1)
        all_pairs.extend(pairwise_r2_by_distance(geno, positions,
                                                   max_dist=bin_edges[-1]))
    return bin_pairs(all_pairs, bin_edges)


if __name__ == "__main__":
    r_theory = 5.86e-8
    print(f"Testing theoretical estimate: r = {r_theory:.3e}")
    print(f"{'bin (bp)':>14} | {'n_pairs':>7} | {'mean r2':>8}")
    print("-" * 38)
    rows = ld_decay_curve(r_theory)
    for lo, hi, n_pairs, mean_r2 in rows:
        print(f"{lo:>6}-{hi:<6} | {n_pairs:7d} | {mean_r2:8.4f}")

    near_zero = rows[0][3]
    print("-" * 38)
    print(f"Near-zero-distance r2 (reference 'max'): {near_zero:.4f}")
    print(f"Half-max target: {near_zero/2:.4f}")
    print("Looking for which bin straddles 711bp and whether its mean r2")
    print("is close to half the near-zero reference value.")
