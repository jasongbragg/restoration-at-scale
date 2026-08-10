"""
stage1/utils.py

Shared utilities for Stage 1 (coalescent ancestry) simulation and
calibration. Imports locked parameters from params.py at the repo root
rather than defining its own constants.

For calibration history (how each parameter was derived), see the scripts
in stage1/calibration/ and docs/PROJECT_SUMMARY.md Section 5.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import msprime
import numpy as np
import params


def build_demography(K=None, deme_Ne=None, m=None,
                     ancestral_merge_time=None, ancestral_Ne=None):
    """Linear (non-circular) stepping-stone chain of K demes.

    All arguments default to the locked values in params.py.

    ancestral_merge_time/ancestral_Ne: if both are set, all K demes are
    merged into a single ancestral population at this depth.

    IMPORTANT: ancestral_Ne must be the structured metapopulation's true
    global effective size (params.ANCESTRAL_NE ~101,053), NOT params.DEME_NE.
    Using DEME_NE as the ancestral size was a real bug caught by the
    integration-level validation script -- it silently truncated total
    diversity while leaving Fst unchanged. See docs/PROJECT_SUMMARY.md
    Section 5d.
    """
    K = K if K is not None else params.K
    deme_Ne = deme_Ne if deme_Ne is not None else params.DEME_NE
    m = m if m is not None else params.M
    if ancestral_merge_time is None:
        ancestral_merge_time = params.ANCESTRAL_MERGE_TIME
    if ancestral_Ne is None:
        ancestral_Ne = params.ANCESTRAL_NE

    demography = msprime.Demography()
    for i in range(K):
        demography.add_population(name=f"deme{i}", initial_size=deme_Ne)
    for i in range(K - 1):
        demography.set_migration_rate(source=f"deme{i}", dest=f"deme{i+1}", rate=m)
        demography.set_migration_rate(source=f"deme{i+1}", dest=f"deme{i}", rate=m)
    if ancestral_merge_time is not None:
        demography.add_population(name="ancestral", initial_size=ancestral_Ne)
        demography.add_population_split(
            time=ancestral_merge_time,
            derived=[f"deme{i}" for i in range(K)],
            ancestral="ancestral",
        )
        demography.sort_events()
    return demography


def branch_gst(ts, K=None, n_per_deme=None):
    """Unbiased Gst-style Fst from a tree sequence using tskit's branch-mode
    diversity. No mutations required; no small-sample bias of the kind that
    plagued the earlier hand-rolled Nei's Gst estimator (see
    docs/PROJECT_SUMMARY.md Section 6, lesson 4).
    """
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME
    sample_sets = [list(range(i * n_per_deme, (i + 1) * n_per_deme)) for i in range(K)]
    all_samples = list(range(K * n_per_deme))
    Ht = ts.diversity(sample_sets=[all_samples], mode="branch")[0]
    Hs = ts.diversity(sample_sets=sample_sets, mode="branch").mean()
    return float((Ht - Hs) / Ht)


def build_chromosome_ts(chrom_idx, seed=1, K=None, deme_Ne=None, m=None,
                         L=None, r=None, n_per_deme=None,
                         ancestral_merge_time=None, ancestral_Ne=None):
    """Build one chromosome's ancestry tree sequence.

    All arguments default to the locked values in params.py.
    This is the natural unit of work for the parallel harness in
    stage1/build_groundtruth.py -- one call per chromosome, entirely
    independent of all others.
    """
    K = K if K is not None else params.K
    deme_Ne = deme_Ne if deme_Ne is not None else params.DEME_NE
    m = m if m is not None else params.M
    L = L if L is not None else params.CHROM_LENGTH
    r = r if r is not None else params.R
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME

    chrom_seed = seed + 1000 * chrom_idx
    demography = build_demography(K=K, deme_Ne=deme_Ne, m=m,
                                   ancestral_merge_time=ancestral_merge_time,
                                   ancestral_Ne=ancestral_Ne)
    samples = {f"deme{i}": n_per_deme for i in range(K)}
    ts = msprime.sim_ancestry(
        samples=samples, demography=demography, sequence_length=L,
        recombination_rate=r, random_seed=chrom_seed, ploidy=1,
    )
    return ts


# ---------------------------------------------------------------------------
# Legacy estimator -- retained for reference only, superseded by branch_gst
# ---------------------------------------------------------------------------

def nei_gst_corrected(geno, K, n_per_deme):
    """
    SUPERSEDED by branch_gst(). Retained only so the calibration history
    scripts in stage1/calibration/ still work for reference.

    This estimator only corrected the within-deme (Hs) term for small-sample
    bias, leaving the pooled (Ht) term uncorrected -- a residual bias that
    caused m=0.0024 to be calibrated against the wrong Fst value (~0.086,
    not 0.05). See docs/PROJECT_SUMMARY.md Section 6, lessons 1 and 4.
    Do NOT use this for any new analysis.
    """
    n_sites = geno.shape[0]
    geno3 = geno.reshape(n_sites, K, n_per_deme)
    p_k = geno3.mean(axis=2)
    p_total = geno.mean(axis=1)
    Ht = 2 * p_total * (1 - p_total)
    correction = n_per_deme / (n_per_deme - 1)
    Hs = (correction * 2 * p_k * (1 - p_k)).mean(axis=1)
    valid = Ht > 0
    fst_per_site = (Ht[valid] - Hs[valid]) / Ht[valid]
    return float(fst_per_site.mean()), int(valid.sum())
