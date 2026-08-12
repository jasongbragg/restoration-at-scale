"""
stage2/diversity.py

Diversity outcome statistics for Stage 2 collection designs.

Takes the output of stage2/sampling.sample_design() and computes
genetic diversity metrics for:
  - the SEED COHORT (diploid, with meiotic recombination)
  - the MATERNAL LINES (haploid, both chromosomes of each mother)

INPUT FORMAT
------------
sample_result['chrom_geno'] is the primary data source, a dict keyed by
chromosome index with entries:
  geno_seeds    : (n_sites, 2*n_seeds) int8 -- recombinant diploid seeds
                  cols 2i = maternal gamete, cols 2i+1 = paternal gamete
  geno_maternal : (n_sites, 2*n_mothers) int8 -- maternal haploid data
                  cols 2i, 2i+1 = the two chromosomes of mother i
  positions     : (n_sites,) float64 -- site positions
  seq_length    : float

BACKGROUND MODE
---------------
compute_background() still uses tree sequences (branch-mode diversity,
no mutations needed). This path is separate from the seed/maternal stats.

METRICS
-------
Standard: pi, He, n_seg_sites, AFS, ROH
Haplotype: n_distinct_mean, h_mean (haplotype diversity), n_eff_hap_mean
           per sliding window of window_size bp.

For ROH: the seed cohort is diploid (pairs cols 2i, 2i+1). ROH lengths
reflect the actual recombination history because meiotic recombination
was modelled during seed generation.

For haplotype metrics: because meiosis IS modelled, each seed's maternal
gamete is genuinely recombinant. hap_seed_* metrics reflect true mosaic
diversity. hap_mat_* metrics are on the maternal haploid cohort (both
chromosomes per mother pooled as independent haplotypes).
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
    n_individuals: int
    n_chromosomes: int
    total_bp: float
    pi_mean: float
    pi_per_chrom: list
    n_seg_sites_total: int
    seg_sites_per_bp: float
    He_mean: float
    He_per_chrom: list
    afs: list
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
# Core statistics from genotype matrix
# ---------------------------------------------------------------------------

def _stats_from_geno(geno: np.ndarray, seq_length: float):
    """pi, He, n_seg, AFS from a (n_sites, n_haploid) genotype matrix."""
    if geno.shape[0] == 0:
        n = geno.shape[1]
        return dict(pi=0.0, He=0.0, n_seg=0, afs=np.zeros(n // 2 + 1))
    n = geno.shape[1]
    counts = geno.sum(axis=1)
    p = counts / n
    seg = (p > 0) & (p < 1)
    He = float((2 * p[seg] * (1 - p[seg])).mean()) if seg.sum() > 0 else 0.0
    correction = n / (n - 1) if n > 1 else 1.0
    pi = float((correction * 2 * p[seg] * (1 - p[seg])).mean() / seq_length) \
        if seg.sum() > 0 else 0.0
    folded = np.minimum(counts[seg].astype(int), n - counts[seg].astype(int))
    afs = np.bincount(folded, minlength=n // 2 + 1).astype(float)
    return dict(pi=pi, He=He, n_seg=int(seg.sum()), afs=afs)


# ---------------------------------------------------------------------------
# ROH from genotype matrix (diploid pairs)
# ---------------------------------------------------------------------------

def _roh_from_diploid_geno(geno: np.ndarray, positions: np.ndarray,
                             chrom_length: float, min_length: int,
                             n_individuals: int):
    """
    ROH per diploid individual from genotype matrix.
    For seeds: cols 2i (maternal) and 2i+1 (paternal).
    For mothers: cols 2i and 2i+1 are the two maternal chromosomes.
    """
    roh_by_ind = []
    for i in range(n_individuals):
        h1 = geno[:, 2 * i]
        h2 = geno[:, 2 * i + 1]
        het_pos = positions[h1 != h2]
        boundaries = np.concatenate([[0.0], het_pos, [chrom_length]])
        roh = [float(boundaries[k + 1] - boundaries[k])
               for k in range(len(boundaries) - 1)
               if boundaries[k + 1] - boundaries[k] >= min_length]
        roh_by_ind.append(roh)
    return roh_by_ind


def _aggregate_roh(roh_per_chrom, total_bp):
    if not roh_per_chrom:
        return dict(roh_mean_n=0.0, roh_mean_total_length=0.0,
                    roh_mean_longest=0.0, roh_fraction_genome=0.0)
    n_ind = len(roh_per_chrom[0])
    all_roh = [[] for _ in range(n_ind)]
    for chrom_roh in roh_per_chrom:
        for i, r in enumerate(chrom_roh):
            all_roh[i].extend(r)
    totals  = [sum(r) for r in all_roh]
    counts  = [len(r) for r in all_roh]
    longest = [max(r) if r else 0.0 for r in all_roh]
    return dict(
        roh_mean_n=float(np.mean(counts)),
        roh_mean_total_length=float(np.mean(totals)),
        roh_mean_longest=float(np.mean(longest)),
        roh_fraction_genome=float(np.mean(totals)) / total_bp if total_bp > 0 else 0.0,
    )


# ---------------------------------------------------------------------------
# Haplotype diversity (SNP window)
# ---------------------------------------------------------------------------

def _haplotype_stats_from_geno(geno, positions, seq_length, window_size,
                                 n_individuals, use_pairs=False):
    """
    Haplotype diversity per sliding window.
    If use_pairs=True, each 'individual' is a diploid pair (cols 2i, 2i+1),
    and its haplotype in a window is the concatenation of both allele vectors.
    If use_pairs=False, each column is an independent haploid haplotype.
    """
    n_sites, n_samp = geno.shape
    if n_sites == 0:
        return dict(n_distinct_mean=0.0, h_mean=0.0,
                    n_eff_hap_mean=0.0, frac_windows_novar=1.0)

    n_windows = max(1, int(seq_length // window_size))
    per_window = []
    n_empty = 0

    for w in range(n_windows):
        lo, hi = w * window_size, (w + 1) * window_size
        mask = (positions >= lo) & (positions < hi)
        if mask.sum() == 0:
            n_empty += 1
            per_window.append((1, 0.0, 1.0))
            continue
        g_win = geno[mask]
        if use_pairs:
            haps = [tuple(g_win[:, 2*i]) + tuple(g_win[:, 2*i+1])
                    for i in range(n_individuals)]
        else:
            haps = [tuple(g_win[:, i]) for i in range(n_samp)]
        counts = {}
        for h in haps:
            counts[h] = counts.get(h, 0) + 1
        freqs = np.array(list(counts.values()), dtype=float) / len(haps)
        per_window.append((len(counts),
                            float(1 - (freqs**2).sum()),
                            float(1 / (freqs**2).sum())))

    return dict(
        n_distinct_mean=float(np.mean([s[0] for s in per_window])),
        h_mean=float(np.mean([s[1] for s in per_window])),
        n_eff_hap_mean=float(np.mean([s[2] for s in per_window])),
        frac_windows_novar=n_empty / n_windows,
    )


# ---------------------------------------------------------------------------
# Background mode (tree-sequence based, no chrom_geno)
# ---------------------------------------------------------------------------

def _pi_from_ts(ts):
    if ts.num_mutations == 0:
        return 0.0
    n = ts.num_samples
    return float(ts.diversity(sample_sets=[list(range(n))], mode="site")[0])


def _seg_sites_from_ts(ts):
    if ts.num_mutations == 0:
        return 0
    return int(round(ts.segregating_sites(mode="site") * ts.sequence_length))


def compute_background(groundtruth_dir: str, world_idx: int = 0,
                        roh_min_length: int = 100_000) -> DiversityResult:
    """Wild-population baseline from full ground truth tree sequences."""
    pattern = os.path.join(groundtruth_dir, f"world{world_idx:02d}_chr*.trees")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No ground truth files at {pattern}.")

    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0
    roh_per_chrom = []
    n_per_deme = params.N_PER_DEME
    K = params.K
    n_ind = K * n_per_deme // 2   # diploid individuals in background

    for p in paths:
        ts = tskit.load(p)
        ts = msprime.sim_mutations(
            ts, rate=params.MU,
            random_seed=world_idx * 1000 +
                int(os.path.basename(p).split("_chr")[1].split(".")[0]) + 1)
        chrom_len = ts.sequence_length
        total_bp += chrom_len

        pi_per_chrom.append(_pi_from_ts(ts))
        n_seg_total += _seg_sites_from_ts(ts)

        if ts.num_mutations > 0:
            geno = ts.genotype_matrix()
            n = geno.shape[1]
            p_arr = geno.sum(axis=1) / n
            seg = (p_arr > 0) & (p_arr < 1)
            He_per_chrom.append(float((2*p_arr[seg]*(1-p_arr[seg])).mean())
                                  if seg.sum() > 0 else 0.0)
            folded = np.minimum(geno.sum(axis=1)[seg].astype(int),
                                 n - geno.sum(axis=1)[seg].astype(int))
            afs_c = np.bincount(folded, minlength=n // 2 + 1).astype(float)
            afs_total = afs_c if afs_total is None else afs_total + afs_c
            positions = np.array([v.position for v in ts.variants()])
            roh_by_ind = _roh_from_diploid_geno(
                geno, positions, chrom_len, roh_min_length, n_ind)
            roh_per_chrom.append(roh_by_ind)
        else:
            He_per_chrom.append(0.0)

    n_samp = params.K * params.N_PER_DEME
    if afs_total is None:
        afs_total = np.zeros(n_samp // 2 + 1)
    roh_stats = _aggregate_roh(roh_per_chrom, total_bp)

    return DiversityResult(
        n_individuals=n_ind,
        n_chromosomes=len(paths),
        total_bp=total_bp,
        pi_mean=float(np.mean(pi_per_chrom)),
        pi_per_chrom=pi_per_chrom,
        n_seg_sites_total=n_seg_total,
        seg_sites_per_bp=n_seg_total / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(He_per_chrom)),
        He_per_chrom=He_per_chrom,
        afs=afs_total.tolist(),
        roh_min_length_bp=roh_min_length,
        **roh_stats,
    )


# ---------------------------------------------------------------------------
# Seed and maternal diversity from chrom_geno
# ---------------------------------------------------------------------------

def compute_diversity(sample_result: dict,
                       roh_min_length: int = 100_000) -> DiversityResult:
    """Diversity of the SEED COHORT from chrom_geno (diploid, recombinant)."""
    chrom_geno = sample_result["chrom_geno"]
    chroms = sorted(chrom_geno.keys())
    n_seeds = sample_result["n_seeds"]

    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0
    roh_per_chrom = []

    for c in chroms:
        gd = chrom_geno[c]
        geno = gd["geno_seeds"]        # (n_sites, 2*n_seeds)
        positions = gd["positions"]
        chrom_len = gd["seq_length"]
        total_bp += chrom_len

        if geno.shape[0] == 0:
            pi_per_chrom.append(0.0)
            He_per_chrom.append(0.0)
            continue

        st = _stats_from_geno(geno, chrom_len)
        pi_per_chrom.append(st["pi"])
        He_per_chrom.append(st["He"])
        n_seg_total += st["n_seg"]
        afs_total = st["afs"] if afs_total is None else afs_total + st["afs"]
        roh_per_chrom.append(
            _roh_from_diploid_geno(geno, positions, chrom_len,
                                    roh_min_length, n_seeds))

    if afs_total is None:
        afs_total = np.zeros(2 * n_seeds // 2 + 1)
    roh_stats = _aggregate_roh(roh_per_chrom, total_bp)

    return DiversityResult(
        n_individuals=n_seeds,
        n_chromosomes=len(chroms),
        total_bp=total_bp,
        pi_mean=float(np.mean(pi_per_chrom)),
        pi_per_chrom=pi_per_chrom,
        n_seg_sites_total=n_seg_total,
        seg_sites_per_bp=n_seg_total / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(He_per_chrom)),
        He_per_chrom=He_per_chrom,
        afs=afs_total.tolist(),
        roh_min_length_bp=roh_min_length,
        **roh_stats,
    )


def compute_maternal_diversity(sample_result: dict,
                                roh_min_length: int = 100_000) -> DiversityResult:
    """
    Diversity of the MATERNAL LINES from chrom_geno.
    geno_maternal has both chromosomes of each mother (cols 2i, 2i+1).
    Standard stats (pi, He, seg_sites) treat all 2*n_mothers columns as
    a haploid cohort. ROH is computed per DIPLOID mother (pairs 2i, 2i+1),
    reflecting within-mother heterozygosity -- meaningful since these are
    real diploid individuals with two distinct chromosomes.
    """
    chrom_geno = sample_result["chrom_geno"]
    chroms = sorted(chrom_geno.keys())
    n_mothers = sample_result["design"].n_mothers()

    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0
    roh_per_chrom = []

    for c in chroms:
        gd = chrom_geno[c]
        geno = gd["geno_maternal"]     # (n_sites, 2*n_mothers)
        positions = gd["positions"]
        chrom_len = gd["seq_length"]
        total_bp += chrom_len

        if geno.shape[0] == 0:
            pi_per_chrom.append(0.0)
            He_per_chrom.append(0.0)
            continue

        st = _stats_from_geno(geno, chrom_len)
        pi_per_chrom.append(st["pi"])
        He_per_chrom.append(st["He"])
        n_seg_total += st["n_seg"]
        afs_total = st["afs"] if afs_total is None else afs_total + st["afs"]
        roh_per_chrom.append(
            _roh_from_diploid_geno(geno, positions, chrom_len,
                                    roh_min_length, n_mothers))

    if afs_total is None:
        afs_total = np.zeros(2 * n_mothers // 2 + 1)
    roh_stats = _aggregate_roh(roh_per_chrom, total_bp)

    return DiversityResult(
        n_individuals=n_mothers,
        n_chromosomes=len(chroms),
        total_bp=total_bp,
        pi_mean=float(np.mean(pi_per_chrom)),
        pi_per_chrom=pi_per_chrom,
        n_seg_sites_total=n_seg_total,
        seg_sites_per_bp=n_seg_total / total_bp if total_bp > 0 else 0.0,
        He_mean=float(np.mean(He_per_chrom)),
        He_per_chrom=He_per_chrom,
        afs=afs_total.tolist(),
        roh_min_length_bp=roh_min_length,
        **roh_stats,
    )


def compute_haplotype_diversity(sample_result: dict,
                                 window_size: int = 100_000) -> dict:
    """
    SNP-window haplotype diversity for seed cohort and maternal lines.
    Returns flat dict with hap_seed_* and hap_mat_* prefixed columns.

    For seeds (use_pairs=True): each seed's 'haplotype' in a window is
    its maternal-gamete SNPs concatenated with paternal-gamete SNPs.
    Because meiotic recombination IS modelled, maternal gametes are
    genuinely recombinant mosaics -- hap_seed_* reflects true diversity.

    For maternal lines (use_pairs=False): each of the 2*n_mothers haploid
    columns is treated as an independent haplotype.
    """
    chrom_geno = sample_result["chrom_geno"]
    n_seeds = sample_result["n_seeds"]
    n_mothers = sample_result["design"].n_mothers()
    chroms = sorted(chrom_geno.keys())

    seed_lists, mat_lists = [], []
    for c in chroms:
        gd = chrom_geno[c]
        if gd["geno_seeds"].shape[0] == 0:
            continue
        positions = gd["positions"]
        chrom_len = gd["seq_length"]
        seed_lists.append(_haplotype_stats_from_geno(
            gd["geno_seeds"], positions, chrom_len,
            window_size, n_seeds, use_pairs=True))
        mat_lists.append(_haplotype_stats_from_geno(
            gd["geno_maternal"], positions, chrom_len,
            window_size, n_mothers, use_pairs=False))

    def _avg(lst, key):
        vals = [s[key] for s in lst if key in s]
        return float(np.mean(vals)) if vals else 0.0

    out = {}
    for prefix, lst in [("hap_seed_", seed_lists), ("hap_mat_", mat_lists)]:
        for k in ("n_distinct_mean", "h_mean", "n_eff_hap_mean",
                   "frac_windows_novar"):
            out[f"{prefix}{k}"] = _avg(lst, k)
    out["hap_window_size_bp"] = window_size
    return out


# ---------------------------------------------------------------------------
# Comparison and convenience functions
# ---------------------------------------------------------------------------

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


def design_stats_dict(sample_result: dict,
                       background: Optional[DiversityResult] = None,
                       roh_min_length: int = 100_000,
                       hap_window_size: int = 100_000) -> dict:
    """
    Compute all diversity metrics for one design result and return as a
    single flat dict suitable for writing to CSV.

    Column groups:
      seed_*   : seed cohort (diploid, recombinant)
      mat_*    : maternal lines (diploid mothers)
      hap_seed_*, hap_mat_* : haplotype diversity per SNP window
      Provenance: n_sites, mothers_per_site, seeds_per_mother, etc.
    """
    seed_r = compute_diversity(sample_result, roh_min_length=roh_min_length)
    mat_r  = compute_maternal_diversity(sample_result, roh_min_length=roh_min_length)
    hap_d  = compute_haplotype_diversity(sample_result, window_size=hap_window_size)

    if background is not None:
        seed_r = compare_to_background(seed_r, background)
        mat_r  = compare_to_background(mat_r, background)

    out = ({"seed_" + k: v for k, v in seed_r.to_dict().items()} |
           {"mat_"  + k: v for k, v in mat_r.to_dict().items()} |
           hap_d)

    design = sample_result["design"]
    out |= {
        "world_idx":        sample_result["world_idx"],
        "n_sites":          design.n_sites,
        "mothers_per_site": design.mothers_per_site
                            if isinstance(design.mothers_per_site, int)
                            else "unequal",
        "seeds_per_mother": design.seeds_per_mother,
        "n_mothers":        design.n_mothers(),
        "total_seeds":      design.total_n(),
        "site_selection":   design.site_selection,
        "pollen_pool":      design.pollen_pool,
        "seed":             design.seed,
    }
    return out

