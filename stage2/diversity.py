"""
stage2/diversity.py

Diversity outcome statistics for Stage 2 collection designs.

Takes the output of stage2/sampling.sample_design() and computes
genetic diversity metrics for the generated seed cohort.

TWO MODES
---------
Seeds mode (sample_result contains 'seed_structure'):
    The tree sequence has been simplified to the unique haplotype set
    involved in the seeds. seed_structure[(i)] = (maternal_pos, paternal_pos)
    gives each seed's two haplotype positions within that ordering.
    Statistics are computed from the per-seed diploid genotype matrix.

Background mode (no seed_structure, e.g. compute_background()):
    The tree sequence has all K*N_PER_DEME haploid samples. Statistics
    use tskit native site-mode functions directly.

METRICS
-------
pi          : nucleotide diversity per bp (Nei 1987 unbiased estimator
              from genotype matrix: mean of n/(n-1)*2pq across sites)
He          : expected heterozygosity (mean 2pq across seg sites)
n_seg_sites : number of segregating sites across all chromosomes
AFS         : folded allele frequency spectrum (summed across chromosomes)
ROH         : runs of homozygosity per diploid seed (or per synthetic
              diploid pair in background mode). min_length configurable.
"""

import copy
import glob
import os
import sys
from dataclasses import dataclass, asdict
from typing import Optional

import msprime
import numpy as np
import tskit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class DiversityResult:
    n_individuals: int          # seeds (or haploid samples / 2 for background)
    n_chromosomes: int
    total_bp: float

    pi_mean: float              # per-bp, mean across chromosomes
    pi_per_chrom: list

    n_seg_sites_total: int
    seg_sites_per_bp: float

    He_mean: float
    He_per_chrom: list

    afs: list                   # folded, summed across chromosomes

    roh_min_length_bp: int
    roh_mean_n: float
    roh_mean_total_length: float
    roh_mean_longest: float
    roh_fraction_genome: float

    pi_fraction_of_background: Optional[float] = None
    He_fraction_of_background: Optional[float] = None
    seg_sites_fraction_of_background: Optional[float] = None

    def to_dict(self):
        d = asdict(self)
        d.pop("pi_per_chrom")
        d.pop("He_per_chrom")
        d.pop("afs")
        return d


# ---------------------------------------------------------------------------
# Genotype-matrix statistics (used for seeds mode)
# ---------------------------------------------------------------------------

def _stats_from_geno(geno: np.ndarray, seq_length: float):
    """
    Compute pi, He, n_seg_sites, and AFS from a (n_sites, n_haploid) geno
    matrix (0/1 values). n_haploid = 2 * n_seeds for diploid seeds.

    pi uses the Nei (1987) unbiased estimator: n/(n-1) * 2pq per site.
    Returns dict with keys: pi, He, n_seg, afs_counts.
    """
    if geno.shape[0] == 0:
        n_hap = geno.shape[1]
        return dict(pi=0.0, He=0.0, n_seg=0,
                    afs=np.zeros(n_hap // 2 + 1))

    n_hap = geno.shape[1]
    counts = geno.sum(axis=1)           # derived allele count per site
    p = counts / n_hap                  # derived allele frequency
    seg = (p > 0) & (p < 1)

    He_vals = 2 * p[seg] * (1 - p[seg])
    He = float(He_vals.mean()) if seg.sum() > 0 else 0.0

    # Nei unbiased pi: correction factor n/(n-1)
    correction = n_hap / (n_hap - 1) if n_hap > 1 else 1.0
    pi_per_site = correction * 2 * p[seg] * (1 - p[seg])
    pi = float(pi_per_site.mean() / seq_length) if seg.sum() > 0 else 0.0

    # folded AFS: bincount of min(count, n-count)
    folded_counts = np.minimum(counts[seg].astype(int), n_hap - counts[seg].astype(int))
    afs = np.bincount(folded_counts, minlength=n_hap // 2 + 1).astype(float)

    return dict(pi=pi, He=He, n_seg=int(seg.sum()), afs=afs)


# ---------------------------------------------------------------------------
# Background-mode statistics (tskit native, no seed_structure)
# ---------------------------------------------------------------------------

def _pi_from_ts(ts: tskit.TreeSequence) -> float:
    """Per-bp nucleotide diversity (site mode). ts.diversity with
    span_normalise=True already returns per-bp; do NOT divide again."""
    if ts.num_mutations == 0:
        return 0.0
    n = ts.num_samples
    return float(ts.diversity(sample_sets=[list(range(n))], mode="site")[0])


def _seg_sites_from_ts(ts: tskit.TreeSequence) -> int:
    """Number of segregating sites: ts.segregating_sites returns per-bp,
    multiply by sequence_length to get raw count."""
    if ts.num_mutations == 0:
        return 0
    return int(round(ts.segregating_sites(mode="site") * ts.sequence_length))


def _He_from_geno(geno: np.ndarray) -> float:
    if geno.shape[0] == 0:
        return 0.0
    n = geno.shape[1]
    p = geno.sum(axis=1) / n
    seg = (p > 0) & (p < 1)
    return float((2 * p[seg] * (1 - p[seg])).mean()) if seg.sum() > 0 else 0.0


def _afs_from_ts(ts: tskit.TreeSequence, n_samples: int) -> np.ndarray:
    if ts.num_mutations == 0:
        return np.zeros(n_samples // 2 + 1)
    afs = ts.allele_frequency_spectrum(
        sample_sets=[list(range(n_samples))],
        mode="site", polarised=False, span_normalise=False,
    )
    return np.array(afs[: n_samples // 2 + 1])


# ---------------------------------------------------------------------------
# ROH (works for both modes -- takes the per-individual genotype columns)
# ---------------------------------------------------------------------------

def _roh_from_pairs(geno: np.ndarray, positions: np.ndarray,
                     pairs: list, chrom_length: float, min_length: int):
    """
    Compute ROH for each individual defined by pairs of column indices.
    pairs: list of (col_a, col_b) giving the two haplotypes of each diploid.
    Returns a list of ROH length lists, one per individual.
    """
    roh_by_ind = []
    for col_a, col_b in pairs:
        h1 = geno[:, col_a]
        h2 = geno[:, col_b]
        het_pos = positions[h1 != h2]
        boundaries = np.concatenate([[0.0], het_pos, [chrom_length]])
        roh = [
            float(boundaries[k + 1] - boundaries[k])
            for k in range(len(boundaries) - 1)
            if boundaries[k + 1] - boundaries[k] >= min_length
        ]
        roh_by_ind.append(roh)
    return roh_by_ind


def _aggregate_roh(roh_per_chrom: list, total_bp: float) -> dict:
    if not roh_per_chrom:
        return dict(roh_mean_n=0.0, roh_mean_total_length=0.0,
                    roh_mean_longest=0.0, roh_fraction_genome=0.0)
    n_ind = len(roh_per_chrom[0])
    all_roh = [[] for _ in range(n_ind)]
    for chrom_roh in roh_per_chrom:
        for i, ind_roh in enumerate(chrom_roh):
            all_roh[i].extend(ind_roh)
    totals = [sum(r) for r in all_roh]
    counts = [len(r) for r in all_roh]
    longest = [max(r) if r else 0.0 for r in all_roh]
    return dict(
        roh_mean_n=float(np.mean(counts)),
        roh_mean_total_length=float(np.mean(totals)),
        roh_mean_longest=float(np.mean(longest)),
        roh_fraction_genome=float(np.mean(totals)) / total_bp if total_bp > 0 else 0.0,
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_diversity(sample_result: dict,
                       roh_min_length: int = 100_000) -> DiversityResult:
    """
    Compute diversity statistics from one sample_design() result.

    Handles both seeds mode (seed_structure present) and background mode.
    """
    chrom_ts = sample_result["chrom_ts"]
    seed_structure = sample_result.get("seed_structure", None)
    chroms = sorted(chrom_ts.keys())

    if not chroms:
        raise ValueError("sample_result['chrom_ts'] is empty")

    seeds_mode = seed_structure is not None
    n_individuals = (len(seed_structure) if seeds_mode
                     else chrom_ts[chroms[0]].num_samples // 2)

    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0
    roh_per_chrom = []

    for c in chroms:
        ts = chrom_ts[c]
        chrom_len = ts.sequence_length
        total_bp += chrom_len

        if seeds_mode:
            # --- seeds mode: build per-seed genotype matrix ---
            if ts.num_mutations == 0:
                pi_per_chrom.append(0.0)
                He_per_chrom.append(0.0)
                continue

            geno_all = ts.genotype_matrix()  # (n_sites, n_unique_samples)
            positions = np.array([v.position for v in ts.variants()])

            # build seed genotype matrix: interleaved maternal/paternal columns
            seed_cols = [col for pair in seed_structure for col in pair]
            geno_seeds = geno_all[:, seed_cols]  # (n_sites, 2*n_seeds)

            st = _stats_from_geno(geno_seeds, chrom_len)
            pi_per_chrom.append(st["pi"])
            He_per_chrom.append(st["He"])
            n_seg_total += st["n_seg"]
            afs_total = st["afs"] if afs_total is None else afs_total + st["afs"]

            # ROH: use seed_structure pairs directly
            if geno_seeds.shape[0] > 0:
                roh_by_ind = _roh_from_pairs(
                    geno_seeds, positions, seed_structure, chrom_len, roh_min_length
                )
                roh_per_chrom.append(roh_by_ind)

        else:
            # --- background mode: use tskit native statistics ---
            pi_per_chrom.append(_pi_from_ts(ts))
            n_seg_total += _seg_sites_from_ts(ts)

            if ts.num_mutations > 0:
                geno = ts.genotype_matrix()
                positions = np.array([v.position for v in ts.variants()])
                He_per_chrom.append(_He_from_geno(geno))
                n_samp = ts.num_samples
                afs_c = _afs_from_ts(ts, n_samp)
                afs_total = afs_c if afs_total is None else afs_total + afs_c
                # consecutive-pair ROH for background
                pairs = [(2 * i, 2 * i + 1) for i in range(n_samp // 2)]
                roh_by_ind = _roh_from_pairs(
                    geno, positions, pairs, chrom_len, roh_min_length
                )
                roh_per_chrom.append(roh_by_ind)
            else:
                He_per_chrom.append(0.0)

    n_afs = (2 * n_individuals if seeds_mode
             else chrom_ts[chroms[0]].num_samples)
    if afs_total is None:
        afs_total = np.zeros(n_afs // 2 + 1)

    roh_stats = _aggregate_roh(roh_per_chrom, total_bp)

    return DiversityResult(
        n_individuals=n_individuals,
        n_chromosomes=len(chroms),
        total_bp=total_bp,
        pi_mean=float(np.mean(pi_per_chrom)),
        pi_per_chrom=pi_per_chrom,
        n_seg_sites_total=n_seg_total,
        seg_sites_per_bp=n_seg_total / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(He_per_chrom)) if He_per_chrom else 0.0,
        He_per_chrom=He_per_chrom,
        afs=afs_total.tolist(),
        roh_min_length_bp=roh_min_length,
        roh_mean_n=roh_stats["roh_mean_n"],
        roh_mean_total_length=roh_stats["roh_mean_total_length"],
        roh_mean_longest=roh_stats["roh_mean_longest"],
        roh_fraction_genome=roh_stats["roh_fraction_genome"],
    )


def compute_background(groundtruth_dir: str, world_idx: int = 0,
                        roh_min_length: int = 100_000) -> DiversityResult:
    """
    Compute diversity from the full ground truth (wild-population baseline).
    Mutations are overlaid before computing site-mode statistics.
    Cache this result -- it is expensive and reused for every design comparison.
    """
    pattern = os.path.join(groundtruth_dir, f"world{world_idx:02d}_chr*.trees")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(
            f"No ground truth files at {pattern}. "
            f"Run stage1/build_groundtruth.py first."
        )
    chrom_ts = {}
    for p in paths:
        c = int(os.path.basename(p).split("_chr")[1].split(".")[0])
        ts = tskit.load(p)
        ts = msprime.sim_mutations(ts, rate=params.MU,
                                    random_seed=world_idx * 1000 + c + 1)
        chrom_ts[c] = ts
    # background has no seed_structure
    return compute_diversity({"chrom_ts": chrom_ts}, roh_min_length=roh_min_length)


def compare_to_background(result: DiversityResult,
                            background: DiversityResult) -> DiversityResult:
    r = copy.copy(result)
    r.pi_fraction_of_background = (
        result.pi_mean / background.pi_mean if background.pi_mean > 0 else None)
    r.He_fraction_of_background = (
        result.He_mean / background.He_mean if background.He_mean > 0 else None)
    r.seg_sites_fraction_of_background = (
        result.seg_sites_per_bp / background.seg_sites_per_bp
        if background.seg_sites_per_bp > 0 else None)
    return r



def compute_maternal_diversity(sample_result: dict) -> DiversityResult:
    """
    Compute diversity statistics on the maternal lines (haploid cohort).

    Uses the same simplified tree sequences already in sample_result --
    no extra file loading. Maternal haplotypes sit at the 'maternal_positions'
    within the simplified ts sample ordering.

    These are the genotypes you would have if you genotyped the sampled
    mother trees directly -- data typically collected in the field before
    seedlings are propagated.

    Statistics are haploid (one chromosome per mother tree): pi, He,
    n_seg_sites, AFS. ROH is not meaningful for a haploid cohort of
    distinct individuals and is returned as zeros.
    """
    chrom_ts = sample_result["chrom_ts"]
    maternal_positions = sample_result.get("maternal_positions")
    if maternal_positions is None:
        raise ValueError(
            "'maternal_positions' not found in sample_result. "
            "Re-run sample_design() with the current version of sampling.py."
        )

    chroms = sorted(chrom_ts.keys())
    n_mothers = len(maternal_positions)
    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0

    for c in chroms:
        ts = chrom_ts[c]
        chrom_len = ts.sequence_length
        total_bp += chrom_len

        if ts.num_mutations == 0:
            pi_per_chrom.append(0.0)
            He_per_chrom.append(0.0)
            continue

        geno_all = ts.genotype_matrix()
        geno_mat = geno_all[:, maternal_positions]   # (n_sites, n_mothers)

        st = _stats_from_geno(geno_mat, chrom_len)
        pi_per_chrom.append(st["pi"])
        He_per_chrom.append(st["He"])
        n_seg_total += st["n_seg"]
        afs_total = st["afs"] if afs_total is None else afs_total + st["afs"]

    if afs_total is None:
        afs_total = np.zeros(n_mothers // 2 + 1)

    return DiversityResult(
        n_individuals=n_mothers,
        n_chromosomes=len(chroms),
        total_bp=total_bp,
        pi_mean=float(np.mean(pi_per_chrom)),
        pi_per_chrom=pi_per_chrom,
        n_seg_sites_total=n_seg_total,
        seg_sites_per_bp=n_seg_total / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(He_per_chrom)) if He_per_chrom else 0.0,
        He_per_chrom=He_per_chrom,
        afs=afs_total.tolist(),
        roh_min_length_bp=0,
        roh_mean_n=0.0,
        roh_mean_total_length=0.0,
        roh_mean_longest=0.0,
        roh_fraction_genome=0.0,
    )


def design_stats_dict(sample_result: dict,
                       background: Optional[DiversityResult] = None,
                       roh_min_length: int = 100_000) -> dict:
    """
    Compute diversity for BOTH maternal lines and seeds, return as one
    flat dict suitable for writing to CSV.

    Maternal stats are prefixed 'mat_'; seed stats are prefixed 'seed_'.
    Design provenance fields are unprefixed.

    This is the primary output for the batch runner -- one row per
    (design x world x seed) in the results CSV.
    """
    # --- seeds ---
    seed_result = compute_diversity(sample_result, roh_min_length=roh_min_length)
    if background is not None:
        seed_result = compare_to_background(seed_result, background)

    # --- maternal lines ---
    mat_result = compute_maternal_diversity(sample_result)
    if background is not None:
        mat_result = compare_to_background(mat_result, background)

    seed_dict = {"seed_" + k: v for k, v in seed_result.to_dict().items()}
    mat_dict = {"mat_" + k: v for k, v in mat_result.to_dict().items()
                if not k.startswith("roh")}   # roh not meaningful for haploid maternal

    out = {**seed_dict, **mat_dict}

    design = sample_result["design"]
    out["world_idx"] = sample_result["world_idx"]
    out["n_sites"] = design.n_sites
    out["mothers_per_site"] = (design.mothers_per_site
                                if isinstance(design.mothers_per_site, int)
                                else "unequal")
    out["seeds_per_mother"] = design.seeds_per_mother
    out["n_mothers"] = design.n_mothers()
    out["total_seeds"] = design.total_n()
    out["site_selection"] = design.site_selection
    out["pollen_pool"] = design.pollen_pool
    out["seed"] = design.seed

    return out

