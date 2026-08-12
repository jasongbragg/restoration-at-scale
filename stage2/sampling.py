"""
stage2/sampling.py

Stage 2: draw collection-design samples from the Stage 1 ground truth,
and generate diploid seeds via one round of open pollination with proper
meiotic recombination on the maternal side.

DIPLOID MOTHER TREES
--------------------
The ground truth is simulated with ploidy=1 (haploid), with N_PER_DEME
haploid samples per deme. Consecutive pairs of haploid samples represent
the two chromosomes of one diploid mother tree:
  diploid mother i in deme d:
    haploid sample A = d * N_PER_DEME + 2*i
    haploid sample B = d * N_PER_DEME + 2*i + 1

This means the maximum number of diploid mothers per site is
N_PER_DEME // 2 (currently 32 at N_PER_DEME=64).

MEIOTIC RECOMBINATION
---------------------
For each seed, the maternal gamete is a RECOMBINANT haplotype formed by
crossing over between the two maternal chromosomes (A and B above).
Crossovers are drawn from a Poisson(r * L) distribution and positions
from Uniform(0, L). The recombinant gamete is computed using numpy
searchsorted for vectorized crossover assignment.

The paternal gamete is drawn directly (without additional recombination)
from the pollen pool -- a random haploid sample from the metapopulation.
This is a reasonable approximation given the large pollen pool size and
the random sampling of pollen donors.

OUTPUT FORMAT
-------------
sample_design() returns a 'chrom_geno' dict (keyed by chromosome index)
rather than tree sequences for the seeds. Each entry contains:
  geno_seeds    : (n_sites, 2*n_seeds) int8 -- diploid seed genotypes.
                  Columns 2i = maternal gamete of seed i (recombinant),
                  columns 2i+1 = paternal gamete of seed i.
  geno_maternal : (n_sites, n_mothers) int8 -- maternal line haploid
                  genotypes (both haplotypes concatenated per mother).
                  Cols 2i, 2i+1 are the two haplotypes of mother i.
  positions     : (n_sites,) float64 -- site positions in bp
  seq_length    : float -- chromosome length

The 'chrom_ts' dict (tree sequences of unique parental haplotypes) is
also retained for potential future use (e.g. tree-based statistics or
Stage 3 SLiM integration). These are the un-recombined parental sequences.
"""

import json
import os
import sys
from dataclasses import dataclass, asdict
from typing import Union

import msprime
import numpy as np
import tskit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params

DATADIR = os.path.join(os.path.dirname(__file__), "..", "data")


# ---------------------------------------------------------------------------
# Design specification
# ---------------------------------------------------------------------------

@dataclass
class CollectionDesign:
    """One collection-design sampling regime.

    n_sites         : number of demes to treat as collection sites.
    mothers_per_site: int or per-site list. Maximum = N_PER_DEME // 2
                      (each diploid mother occupies 2 haploid samples).
    seeds_per_mother: seeds generated per maternal line.
    site_selection  : 'random', 'even_spacing', or 'contiguous'.
    pollen_pool     : 'metapopulation' (any deme) or 'local' (same deme).
    seed            : RNG seed for full reproducibility.
    relatedness_mode: only 'unrelated' implemented.
    """
    n_sites: int
    mothers_per_site: Union[int, list]
    seeds_per_mother: int = 1
    site_selection: str = "random"
    pollen_pool: str = "local"   # local = paternal gametes from same deme as mother
                                  # (biologically correct for wind-pollinated species
                                  # with predominantly local pollen dispersal; preserves
                                  # Fst in the planted cohort. metapopulation = draw
                                  # from global frequencies, which homogenises structure)
    seed: int = 1
    relatedness_mode: str = "unrelated"

    def mothers_list(self):
        if isinstance(self.mothers_per_site, int):
            return [self.mothers_per_site] * self.n_sites
        if len(self.mothers_per_site) != self.n_sites:
            raise ValueError(
                f"mothers_per_site has length {len(self.mothers_per_site)} "
                f"but n_sites={self.n_sites}"
            )
        return list(self.mothers_per_site)

    def n_mothers(self):
        return sum(self.mothers_list())

    def total_n(self):
        return self.n_mothers() * self.seeds_per_mother

    def label(self):
        mps = self.mothers_per_site if isinstance(self.mothers_per_site, int) \
            else "unequal"
        return (f"{self.site_selection[:4]}_s{self.n_sites}_m{mps}"
                f"_spm{self.seeds_per_mother}_N{self.total_n()}_seed{self.seed}")


# ---------------------------------------------------------------------------
# Capacity validation
# ---------------------------------------------------------------------------

def validate_design_feasibility(design: CollectionDesign, K=None, n_per_deme=None):
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME
    n_diploid_per_deme = n_per_deme // 2  # each mother needs 2 haploid samples

    if n_per_deme % 2 != 0:
        raise ValueError(f"N_PER_DEME={n_per_deme} must be even (each diploid "
                          f"mother occupies 2 consecutive haploid samples).")
    if design.n_sites > K:
        raise ValueError(
            f"Design requests n_sites={design.n_sites} but K={K}.")
    if design.n_sites < 1:
        raise ValueError(f"n_sites must be >=1.")
    if design.seeds_per_mother < 1:
        raise ValueError(f"seeds_per_mother must be >=1.")
    if design.pollen_pool not in ("metapopulation", "local"):
        raise ValueError(f"pollen_pool must be 'metapopulation' or 'local'.")

    for i, n in enumerate(design.mothers_list()):
        if n > n_diploid_per_deme:
            raise ValueError(
                f"Design requests {n} diploid mothers at site {i}, but the "
                f"ground truth has N_PER_DEME={n_per_deme} haploid samples = "
                f"{n_diploid_per_deme} diploid mothers per deme. Either reduce "
                f"mothers_per_site to <={n_diploid_per_deme}, or increase "
                f"N_PER_DEME in params.py and rebuild the ground truth."
            )
        if n < 1:
            raise ValueError(f"mothers_per_site at site {i} must be >=1.")

    if design.relatedness_mode != "unrelated":
        raise NotImplementedError(
            f"relatedness_mode='{design.relatedness_mode}' not implemented.")


# ---------------------------------------------------------------------------
# Site and mother selection
# ---------------------------------------------------------------------------

def select_sites(K, n_sites, strategy, rng):
    if strategy == "random":
        return sorted(rng.choice(K, size=n_sites, replace=False).tolist())
    elif strategy == "even_spacing":
        idx = np.unique(np.round(np.linspace(0, K - 1, n_sites)).astype(int))
        if len(idx) < n_sites:
            raise ValueError(
                f"even_spacing K={K}, n_sites={n_sites} produced only "
                f"{len(idx)} unique demes.")
        return sorted(idx.tolist())
    elif strategy == "contiguous":
        start = int(rng.integers(0, K - n_sites + 1))
        return list(range(start, start + n_sites))
    else:
        raise ValueError(f"Unknown site_selection '{strategy}'.")


def select_diploid_mothers(n_diploid_available, n_mothers, rng):
    """Select n_mothers diploid mother indices (0-based) without replacement."""
    return sorted(rng.choice(n_diploid_available, size=n_mothers,
                              replace=False).tolist())


def _haploid_pair(deme, mother_idx, n_per_deme):
    """Global haploid sample IDs for diploid mother mother_idx in deme."""
    base = deme * n_per_deme + 2 * mother_idx
    return base, base + 1


# ---------------------------------------------------------------------------
# Meiotic gamete generation
# ---------------------------------------------------------------------------

def _meiotic_gamete(col0, col1, positions, chrom_length, r, rng):
    """
    Generate a recombinant maternal gamete from a diploid mother (col0, col1).

    Crossovers: n ~ Poisson(r * chrom_length), positions ~ Uniform(0, L).
    np.searchsorted assigns each site to the active parental haplotype
    without Python-level loops.

    Returns a 1D int8 array of length len(col0).
    """
    n_co = rng.poisson(r * chrom_length)
    if n_co == 0:
        return col0.copy() if rng.integers(2) == 0 else col1.copy()
    co_pos = np.sort(rng.uniform(0, chrom_length, n_co))
    n_left = np.searchsorted(co_pos, positions, side='right')
    start = rng.integers(2)
    active = (n_left + start) % 2   # 0 -> col0, 1 -> col1
    return np.where(active == 0, col0, col1).astype(np.int8)


# ---------------------------------------------------------------------------
# Build per-chromosome seed and maternal genotype matrices
# ---------------------------------------------------------------------------

def _build_chrom_geno(ts, site_demes, mothers_by_site, n_per_deme,
                       seeds_per_mother, pollen_pool, r, rng):
    """
    Given a full ground-truth chromosome tree sequence, build:
      geno_seeds    : (n_sites, 2*n_seeds) recombinant diploid seed genotypes
      geno_maternal : (n_sites, 2*n_mothers) maternal haplotype genotypes
                      (both haplotypes of each mother, cols 2i and 2i+1)
      positions     : (n_sites,) site positions
      seq_length    : float

    Seeds are generated as: recombinant maternal gamete + paternal gamete.
    Maternal gamete: meiotic recombination between the two paired maternal
      haploid samples (samples 2i and 2i+1 for mother i in a deme).
    Paternal gamete: a random haploid sample from the pollen pool
      (metapopulation = all K*n_per_deme samples, or local = same deme).
      No additional meiosis simulated on the paternal side.
    """
    K_total = ts.num_samples // n_per_deme  # inferred from tree sequence
    chrom_length = ts.sequence_length

    if ts.num_sites == 0:
        # no mutations on this chromosome -- genotypes are all zero
        n_mothers_total = sum(len(m) for m in mothers_by_site)
        n_seeds = n_mothers_total * seeds_per_mother
        return dict(
            geno_seeds=np.zeros((0, 2 * n_seeds), dtype=np.int8),
            geno_maternal=np.zeros((0, 2 * n_mothers_total), dtype=np.int8),
            positions=np.array([]),
            seq_length=chrom_length,
        )

    geno_full = ts.genotype_matrix()   # (n_sites, K_total * n_per_deme)
    positions = np.array([v.position for v in ts.variants()])
    n_sites = geno_full.shape[0]
    n_per_deme_half = n_per_deme // 2  # diploid mothers available per deme

    # Pre-compute allele frequencies for the virtual pollen pool.
    # Rather than drawing paternal gametes from one of the K*n_per_deme fixed
    # ground-truth haplotypes (which limits pollen diversity to 1024 options
    # and makes rare alleles noisy), we draw from the ALLELE FREQUENCY
    # DISTRIBUTION itself: at each site, the paternal allele is 1 with
    # probability = the observed frequency in the full ground truth.
    #
    # This is biologically appropriate for a wind-pollinated species with low
    # correlated paternity (as observed in real M. quinquenervia collections):
    # each seed's pollen donor is effectively an independent random draw from
    # the metapopulation frequency distribution. It gives an effectively
    # infinite pollen pool, removes the N_PER_DEME ceiling on paternal
    # diversity, and makes rare allele tracking in the SFS more accurate
    # (no finite-sample noise from drawing from 1024 fixed haplotypes).
    #
    # Tradeoff: LD in the paternal contribution is not preserved (each site
    # drawn independently). For unrelated random pollen donors from a large
    # population this is a good approximation -- inter-site correlations
    # in pollen average out across many independent donors.
    #
    # pollen_pool="local" uses only haplotypes from the same deme.
    # pollen_pool="metapopulation" uses all K demes (default).
    all_col_idx = np.arange(geno_full.shape[1])
    # metapopulation-wide allele frequencies (used for virtual pool)
    p_meta = geno_full.mean(axis=1)   # (n_sites,) float, per-site frequency

    n_mothers_total = sum(len(m) for m in mothers_by_site)
    n_seeds = n_mothers_total * seeds_per_mother

    geno_seeds    = np.empty((n_sites, 2 * n_seeds), dtype=np.int8)
    geno_maternal = np.empty((n_sites, 2 * n_mothers_total), dtype=np.int8)

    seed_col = 0
    mat_col = 0
    for deme, mother_indices in zip(site_demes, mothers_by_site):
        deme_base = deme * n_per_deme

        # local allele frequencies (used when pollen_pool="local")
        if pollen_pool == "local":
            p_pollen = geno_full[:, deme_base:deme_base + n_per_deme].mean(axis=1)
        else:
            p_pollen = p_meta

        for mother_idx in mother_indices:
            col_a = deme_base + 2 * mother_idx
            col_b = deme_base + 2 * mother_idx + 1

            hap_a = geno_full[:, col_a]
            hap_b = geno_full[:, col_b]

            geno_maternal[:, mat_col]     = hap_a
            geno_maternal[:, mat_col + 1] = hap_b
            mat_col += 2

            # generate seeds: recombinant maternal gamete + virtual paternal gamete
            for _ in range(seeds_per_mother):
                mat_gamete = _meiotic_gamete(hap_a, hap_b, positions,
                                              chrom_length, r, rng)
                # paternal gamete: draw each allele independently from the
                # pollen frequency distribution (virtual infinite pool)
                # Bernoulli draw: site i is allele 1 with prob p_pollen[i]
                pat_gamete = (rng.random(n_sites) < p_pollen).astype(np.int8)

                geno_seeds[:, seed_col]     = mat_gamete
                geno_seeds[:, seed_col + 1] = pat_gamete
                seed_col += 2

    return dict(
        geno_seeds=geno_seeds,
        geno_maternal=geno_maternal,
        positions=positions,
        seq_length=chrom_length,
    )


# ---------------------------------------------------------------------------
# Core sampling entry point
# ---------------------------------------------------------------------------

def sample_design(design: CollectionDesign, world_idx: int = 0,
                   groundtruth_dir: str = None,
                   K: int = None, n_per_deme: int = None,
                   n_chrom: int = None, mu: float = None, r: float = None):
    """
    Draw one collection-design realization from the ground truth, generating
    diploid seeds via open pollination with meiotic recombination.

    Returns
    -------
    dict with:
      'site_demes'     : list of deme indices used as collection sites
      'mothers_by_site': list of diploid mother indices per site (0-based
                         within each deme; maps to haploid samples 2i, 2i+1)
      'n_seeds'        : total seeds generated
      'chrom_geno'     : {chrom_idx: dict with geno_seeds, geno_maternal,
                         positions, seq_length}
      'chrom_ts'       : {chrom_idx: TreeSequence of unique parental haplotypes}
                         (retained for future use; no additional recombination
                          overlaid)
      'design'         : CollectionDesign used
      'world_idx'      : ground-truth world index
    """
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME
    n_chrom = n_chrom if n_chrom is not None else params.N_CHROM
    mu = mu if mu is not None else params.MU
    r = r if r is not None else params.R

    if groundtruth_dir is None:
        groundtruth_dir = os.path.join(DATADIR, "groundtruth")

    validate_design_feasibility(design, K=K, n_per_deme=n_per_deme)

    n_diploid_per_deme = n_per_deme // 2
    rng = np.random.default_rng(design.seed)

    site_demes = select_sites(K, design.n_sites, design.site_selection, rng)
    mothers_by_site = [
        select_diploid_mothers(n_diploid_per_deme, n_mothers, rng)
        for n_mothers in design.mothers_list()
    ]
    n_seeds = sum(len(m) for m in mothers_by_site) * design.seeds_per_mother

    chrom_geno = {}
    chrom_ts   = {}

    for c in range(n_chrom):
        path = os.path.join(groundtruth_dir,
                            f"world{world_idx:02d}_chr{c:02d}.trees")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"Ground truth file not found: {path} -- "
                f"has stage1/build_groundtruth.py finished world {world_idx}?"
            )
        ts = tskit.load(path)
        ts = msprime.sim_mutations(ts, rate=mu,
                                    random_seed=design.seed + 1000 * c + 1)

        chrom_geno[c] = _build_chrom_geno(
            ts, site_demes, mothers_by_site, n_per_deme,
            design.seeds_per_mother, design.pollen_pool, r, rng,
        )
        chrom_ts[c] = ts   # retain full ts for future use

    return {
        "site_demes":      site_demes,
        "mothers_by_site": mothers_by_site,
        "n_seeds":         n_seeds,
        "chrom_geno":      chrom_geno,
        "chrom_ts":        chrom_ts,
        "design":          design,
        "world_idx":       world_idx,
    }


# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------

def save_sampled_design(result: dict, label: str = None,
                         designs_dir: str = None):
    if designs_dir is None:
        designs_dir = os.path.join(DATADIR, "designs")
    if label is None:
        label = f"world{result['world_idx']:02d}_{result['design'].label()}"

    out = os.path.join(designs_dir, label)
    os.makedirs(out, exist_ok=True)

    # save genotype matrices as compressed numpy archives
    for c, gd in result["chrom_geno"].items():
        np.savez_compressed(
            os.path.join(out, f"chr{c:02d}_geno.npz"),
            geno_seeds=gd["geno_seeds"],
            geno_maternal=gd["geno_maternal"],
            positions=gd["positions"],
            seq_length=np.array([gd["seq_length"]]),
        )

    manifest = {
        "design":      asdict(result["design"]),
        "world_idx":   result["world_idx"],
        "site_demes":  result["site_demes"],
        "mothers_by_site": result["mothers_by_site"],
        "n_seeds":     result["n_seeds"],
        "params_snapshot": {
            "K": params.K, "N_PER_DEME": params.N_PER_DEME,
            "N_CHROM": params.N_CHROM, "CHROM_LENGTH": params.CHROM_LENGTH,
            "MU": params.MU, "R": params.R,
        },
    }
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    return out

