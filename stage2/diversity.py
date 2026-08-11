"""
stage2/diversity.py

Diversity outcome statistics for Stage 2 (collection-design sampling).
Takes the output of stage2/sampling.sample_design() and computes the
genetic diversity metrics that are the primary deliverable of the analysis.

METRICS
-------
The following statistics are computed across all chromosomes in a sampled
design and compared against the wild-population background:

  pi (nucleotide diversity)
      tskit site-mode diversity, averaged across chromosomes and normalised
      to per-bp. Matches what you would compute from real sequencing data.

  He (expected heterozygosity)
      Mean 2pq across all segregating sites across all chromosomes.
      Computed from the allele frequency spectrum rather than from paired
      genotypes, so it doesn't depend on the haploid/diploid pairing.

  n_seg_sites / seg_sites_per_bp
      Total number of segregating sites (summed across chromosomes) and
      the density per base pair.

  Allele frequency spectrum (AFS)
      Folded, summed across all chromosomes. Shape (n_samples//2 + 1,).

  ROH (runs of homozygosity)
      Computed on synthetic diploid individuals formed by pairing
      consecutive haploid samples: (0,1), (2,3), ...
      See PLOIDY NOTE below.

PLOIDY NOTE
-----------
The ground truth was simulated with ploidy=1 (haploid), meaning each
"mother tree" in the collection is represented by one haploid chromosome
per genome-chromosome. For ROH and observed heterozygosity, a diploid
model is needed. This module pairs consecutive haploid samples to form
synthetic diploid individuals:
  individual 0: haploid samples 0 and 1
  individual 1: haploid samples 2 and 3
  ...
This is a reasonable approximation for a randomly-mating outcrossing
population. The effective diploid sample size is half the haploid sample
count. For n_sites=4, mothers_per_site=8 (total N=32 haploid samples),
we get 16 synthetic diploid individuals.

If a design requests an odd number of total haploid samples, the last
sample is silently dropped from ROH/Ho computation (flagged in output).

OUTPUT FORMAT
-------------
compute_diversity() returns a DiversityResult dataclass. Use .to_dict()
to get a flat dict suitable for writing to CSV/JSON, or pass the whole
result to compare_to_background() to get proportional capture metrics.
"""

import os
import sys
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np
import tskit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class DiversityResult:
    """Diversity statistics for one sampled design (all chromosomes pooled)."""

    # --- provenance --------------------------------------------------------
    n_haploid_samples: int          # total samples in the simplified ts
    n_diploid_individuals: int      # = n_haploid_samples // 2 (for ROH)
    n_chromosomes: int
    total_bp: float

    # --- nucleotide diversity (pi) -----------------------------------------
    pi_mean: float                  # mean per-bp pi across chromosomes
    pi_per_chrom: list              # per-chromosome per-bp values

    # --- segregating sites -------------------------------------------------
    n_seg_sites_total: int
    seg_sites_per_bp: float

    # --- expected heterozygosity -------------------------------------------
    He_mean: float                  # mean 2pq across all seg sites
    He_per_chrom: list              # per-chromosome values

    # --- allele frequency spectrum (folded) --------------------------------
    afs: list                       # shape (n_haploid_samples//2 + 1,)

    # --- ROH (runs of homozygosity, synthetic diploid pairs) ---------------
    roh_min_length_bp: int          # minimum ROH length used
    roh_mean_n: float               # mean number of ROH per individual
    roh_mean_total_length: float    # mean total ROH bp per individual
    roh_mean_longest: float         # mean length of longest ROH
    roh_fraction_genome: float      # mean fraction of genome in ROH
    roh_odd_samples_dropped: bool   # True if one haploid sample was dropped

    # --- comparison to background (filled in by compare_to_background) -----
    pi_fraction_of_background: Optional[float] = None
    He_fraction_of_background: Optional[float] = None
    seg_sites_fraction_of_background: Optional[float] = None

    def to_dict(self):
        d = asdict(self)
        # flatten per-chrom lists to summary stats only for CSV friendliness
        d.pop("pi_per_chrom")
        d.pop("He_per_chrom")
        d.pop("afs")
        return d


# ---------------------------------------------------------------------------
# Core statistics from tskit native functions
# ---------------------------------------------------------------------------

def _pi_from_ts(ts: tskit.TreeSequence) -> float:
    """Per-bp nucleotide diversity (site mode) for all samples in ts.

    ts.diversity(span_normalise=True) -- the default -- already returns
    per-bp diversity. Do NOT divide by sequence_length again; the earlier
    version of this function did so and produced values ~4e-10 instead of
    the correct ~2e-3 for Ne=100,000, mu=1e-8.
    """
    if ts.num_mutations == 0:
        return 0.0
    n = ts.num_samples
    pi = ts.diversity(sample_sets=[list(range(n))], mode="site")[0]
    return float(pi)  # already per-bp


def _seg_sites_from_ts(ts: tskit.TreeSequence) -> int:
    """Number of segregating sites (site mode).
    ts.segregating_sites(mode='site') returns a per-bp rate; multiply by
    sequence_length to get the raw count."""
    if ts.num_mutations == 0:
        return 0
    return int(round(ts.segregating_sites(mode="site") * ts.sequence_length))


def _He_from_geno(geno: np.ndarray) -> float:
    """Mean expected heterozygosity (2pq) across all sites in geno matrix.
    geno: (n_sites, n_haploid_samples), values 0/1."""
    if geno.shape[0] == 0:
        return 0.0
    n = geno.shape[1]
    p = geno.sum(axis=1) / n
    seg = (p > 0) & (p < 1)
    if seg.sum() == 0:
        return 0.0
    return float((2 * p[seg] * (1 - p[seg])).mean())


def _afs_from_ts(ts: tskit.TreeSequence, n_samples: int) -> np.ndarray:
    """Folded allele frequency spectrum (site mode), shape (n//2 + 1,)."""
    if ts.num_mutations == 0:
        return np.zeros(n_samples // 2 + 1)
    afs = ts.allele_frequency_spectrum(
        sample_sets=[list(range(n_samples))],
        mode="site", polarised=False, span_normalise=False,
    )
    # tskit returns shape (n_samples + 1,) for folded; slice to (n//2 + 1,)
    return np.array(afs[: n_samples // 2 + 1])


# ---------------------------------------------------------------------------
# ROH on synthetic diploid pairs
# ---------------------------------------------------------------------------

def _roh_from_geno(geno: np.ndarray, positions: np.ndarray,
                    chrom_length: float, min_length: int) -> list:
    """
    For each synthetic diploid individual (paired consecutive haploid
    samples), find ROH >= min_length bp. Returns a list of lists of ROH
    lengths, one inner list per diploid individual.

    geno:       (n_sites, n_haploid_samples)
    positions:  (n_sites,) physical positions
    """
    n_hap = geno.shape[1]
    n_dip = n_hap // 2
    roh_by_ind = []

    for i in range(n_dip):
        h1 = geno[:, 2 * i]
        h2 = geno[:, 2 * i + 1]
        het_pos = positions[h1 != h2]
        # homozygous stretches are gaps between consecutive het sites
        boundaries = np.concatenate([[0.0], het_pos, [chrom_length]])
        roh = [
            float(boundaries[k + 1] - boundaries[k])
            for k in range(len(boundaries) - 1)
            if boundaries[k + 1] - boundaries[k] >= min_length
        ]
        roh_by_ind.append(roh)

    return roh_by_ind


def _aggregate_roh(roh_per_chrom_per_ind: list, total_bp: float) -> dict:
    """
    Aggregate per-chromosome, per-individual ROH lists into summary stats.

    roh_per_chrom_per_ind: list of (per-individual ROH lists), one entry
        per chromosome. Each entry is a list of length n_diploid_individuals
        containing lists of ROH lengths for that individual on that chrom.
    """
    if not roh_per_chrom_per_ind:
        return dict(roh_mean_n=0.0, roh_mean_total_length=0.0,
                    roh_mean_longest=0.0, roh_fraction_genome=0.0)

    n_dip = len(roh_per_chrom_per_ind[0])
    # flatten across chromosomes per individual
    all_roh = [[] for _ in range(n_dip)]
    for chrom_roh in roh_per_chrom_per_ind:
        for ind_idx, ind_roh in enumerate(chrom_roh):
            all_roh[ind_idx].extend(ind_roh)

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
# Main entry points
# ---------------------------------------------------------------------------

def compute_diversity(sample_result: dict,
                       roh_min_length: int = 100_000) -> DiversityResult:
    """
    Compute diversity statistics from one sample_design() result.

    Parameters
    ----------
    sample_result : dict
        Output of stage2.sampling.sample_design(). Must contain key
        'chrom_ts': {chrom_idx: tskit.TreeSequence}.
    roh_min_length : int
        Minimum ROH length in bp (default 100kb).

    Returns
    -------
    DiversityResult
    """
    chrom_ts = sample_result["chrom_ts"]
    chroms = sorted(chrom_ts.keys())

    if not chroms:
        raise ValueError("sample_result['chrom_ts'] is empty")

    n_samples = chrom_ts[chroms[0]].num_samples
    n_dip = n_samples // 2
    odd_drop = (n_samples % 2 != 0)

    pi_per_chrom, He_per_chrom = [], []
    afs_total = None
    n_seg_total = 0
    total_bp = 0.0
    roh_per_chrom = []

    for c in chroms:
        ts = chrom_ts[c]
        chrom_len = ts.sequence_length
        total_bp += chrom_len

        # nucleotide diversity
        pi_per_chrom.append(_pi_from_ts(ts))

        # segregating sites
        n_seg_total += _seg_sites_from_ts(ts)

        # He and ROH need the genotype matrix
        if ts.num_mutations > 0:
            geno = ts.genotype_matrix()
            positions = np.array([v.position for v in ts.variants()])

            He_per_chrom.append(_He_from_geno(geno))

            # ROH on synthetic diploid pairs
            roh_by_ind = _roh_from_geno(geno, positions, chrom_len, roh_min_length)
            roh_per_chrom.append(roh_by_ind)

            # AFS
            afs_c = _afs_from_ts(ts, n_samples)
            afs_total = afs_c if afs_total is None else afs_total + afs_c
        else:
            He_per_chrom.append(0.0)

    if afs_total is None:
        afs_total = np.zeros(n_samples // 2 + 1)

    roh_stats = _aggregate_roh(roh_per_chrom, total_bp)

    return DiversityResult(
        n_haploid_samples=n_samples,
        n_diploid_individuals=n_dip,
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
        roh_odd_samples_dropped=odd_drop,
    )


def compute_background(groundtruth_dir: str, world_idx: int = 0,
                        roh_min_length: int = 100_000) -> DiversityResult:
    """
    Compute diversity statistics from the full ground-truth tree sequences
    (without any sampling/simplification). Used as the wild-population
    baseline for comparison.

    Ground truth .trees files are stored mutation-free (ancestry only);
    mutations are overlaid here at a fixed seed before computing statistics,
    consistent with how sample_design() works. The seed is deterministic
    given world_idx so the background is reproducible.

    This is relatively expensive (K * N_PER_DEME samples per chromosome).
    Cache the result rather than recomputing it for every design comparison:
        import pickle
        bg = compute_background("data/groundtruth", world_idx=0)
        pickle.dump(bg, open("data/background_world00.pkl", "wb"))
    """
    import glob
    import msprime

    pattern = os.path.join(groundtruth_dir,
                            f"world{world_idx:02d}_chr*.trees")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(
            f"No ground truth files found at {pattern}. "
            f"Run stage1/build_groundtruth.py first."
        )

    chrom_ts = {}
    for p in paths:
        c = int(os.path.basename(p).split("_chr")[1].split(".")[0])
        ts = tskit.load(p)
        # overlay mutations -- ground truth is stored ancestry-only;
        # seed is deterministic per (world, chrom) so background is reproducible
        ts = msprime.sim_mutations(ts, rate=params.MU,
                                    random_seed=world_idx * 1000 + c + 1)
        chrom_ts[c] = ts

    return compute_diversity({"chrom_ts": chrom_ts}, roh_min_length=roh_min_length)


def compare_to_background(result: DiversityResult,
                            background: DiversityResult) -> DiversityResult:
    """
    Fill in the *_fraction_of_background fields by dividing the sampled
    design's statistics by the background (wild population) values.
    Returns a new DiversityResult with those fields populated.
    """
    import copy
    r = copy.copy(result)
    r.pi_fraction_of_background = (
        result.pi_mean / background.pi_mean
        if background.pi_mean > 0 else None
    )
    r.He_fraction_of_background = (
        result.He_mean / background.He_mean
        if background.He_mean > 0 else None
    )
    r.seg_sites_fraction_of_background = (
        result.seg_sites_per_bp / background.seg_sites_per_bp
        if background.seg_sites_per_bp > 0 else None
    )
    return r


# ---------------------------------------------------------------------------
# Convenience: run one design and return a flat dict (for CSV/JSON output)
# ---------------------------------------------------------------------------

def design_stats_dict(sample_result: dict, background: Optional[DiversityResult] = None,
                       roh_min_length: int = 100_000) -> dict:
    """
    Convenience wrapper: compute diversity, optionally compare to background,
    merge the design provenance from sample_result, and return a flat dict.
    """
    from dataclasses import asdict

    result = compute_diversity(sample_result, roh_min_length=roh_min_length)
    if background is not None:
        result = compare_to_background(result, background)

    out = result.to_dict()

    # add design provenance fields
    design = sample_result["design"]
    out["world_idx"] = sample_result["world_idx"]
    out["n_sites"] = design.n_sites
    out["mothers_per_site"] = (
        design.mothers_per_site if isinstance(design.mothers_per_site, int)
        else "unequal"
    )
    out["total_n"] = design.total_n()
    out["site_selection"] = design.site_selection
    out["seed"] = design.seed

    return out

