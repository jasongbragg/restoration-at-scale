"""
stage2/sampling.py

Stage 2: draw collection-design samples from the Stage 1 ground truth,
and generate diploid seeds via one round of open pollination.

BIOLOGICAL MODEL
----------------
At restoration scale you plant seeds, not cuttings of the sampled mother
trees. Each seed is the product of:
  - A maternal gamete: one haplotype drawn from the mother tree
  - A paternal gamete: one haplotype from a random pollen donor

This module implements that crossing step. sample_design() returns diploid
seeds. Each seed's genotype is the combination of two haploid ancestral
chromosomes from the ground truth tree sequence.

OUTPUT STRUCTURE
----------------
sample_design() returns a result dict containing:
  'chrom_ts'       : {chrom_idx: simplified+mutated TreeSequence}
                     Simplified to the UNIQUE set of haploid samples
                     involved (maternal + paternal, deduplicated).
  'seed_structure' : list of (maternal_pos, paternal_pos) per seed,
                     giving each seed's two haplotype positions within
                     the simplified ts's sample ordering.

The diversity module uses seed_structure to construct per-seed diploid
genotype matrices for statistics. ROH and heterozygosity treat position
pairs (maternal_pos, paternal_pos) as the two alleles of each diploid.

POLLEN POOL
-----------
pollen_pool="metapopulation" (default): paternal gametes drawn from any
  individual in the full K-deme ground truth. Appropriate for
  wind-pollinated species with wide pollen dispersal (Eucalyptus,
  Melaleuca). Selfing is excluded.
pollen_pool="local": paternal gametes drawn from the same deme only.
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
    mothers_per_site: single int (same count at every site) or a list of
                      length n_sites for unequal designs.
    seeds_per_mother: seeds retained per maternal line. Total seeds planted
                      = n_sites * mothers_per_site * seeds_per_mother.
    site_selection  : 'random', 'even_spacing', or 'contiguous'.
    pollen_pool     : 'metapopulation' (any deme) or 'local' (same deme).
    seed            : RNG seed for full reproducibility.
    relatedness_mode: only 'unrelated' implemented (PROJECT_SUMMARY.md).
    """
    n_sites: int
    mothers_per_site: Union[int, list]
    seeds_per_mother: int = 1
    site_selection: str = "random"
    pollen_pool: str = "metapopulation"
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
        """Total number of maternal lines sampled."""
        return sum(self.mothers_list())

    def total_n(self):
        """Total seeds generated = n_mothers * seeds_per_mother.
        This is the number of individuals planted."""
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
    """Raise ValueError with specific, actionable message if infeasible."""
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME

    if design.n_sites > K:
        raise ValueError(
            f"Design requests n_sites={design.n_sites} but ground truth has "
            f"only K={K} demes."
        )
    if design.n_sites < 1:
        raise ValueError(f"n_sites must be >=1, got {design.n_sites}")
    if design.seeds_per_mother < 1:
        raise ValueError(f"seeds_per_mother must be >=1, got {design.seeds_per_mother}")
    if design.pollen_pool not in ("metapopulation", "local"):
        raise ValueError(
            f"pollen_pool must be 'metapopulation' or 'local', "
            f"got '{design.pollen_pool}'"
        )
    for i, n in enumerate(design.mothers_list()):
        if n > n_per_deme:
            raise ValueError(
                f"Design requests {n} mothers at site {i}, but the ground "
                f"truth was built with N_PER_DEME={n_per_deme}. Either reduce "
                f"mothers_per_site to <={n_per_deme}, or increase N_PER_DEME "
                f"in params.py and rebuild the ground truth."
            )
        if n < 1:
            raise ValueError(f"mothers_per_site at site {i} must be >=1, got {n}")
    if design.relatedness_mode != "unrelated":
        raise NotImplementedError(
            f"relatedness_mode='{design.relatedness_mode}' not implemented."
        )


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
                f"{len(idx)} unique demes. Try 'random' instead."
            )
        return sorted(idx.tolist())
    elif strategy == "contiguous":
        start = int(rng.integers(0, K - n_sites + 1))
        return list(range(start, start + n_sites))
    else:
        raise ValueError(
            f"Unknown site_selection '{strategy}'. "
            f"Use 'random', 'even_spacing', or 'contiguous'."
        )


def select_mothers(n_per_deme, n_mothers, rng):
    return sorted(rng.choice(n_per_deme, size=n_mothers, replace=False).tolist())


def _global_sample_ids(site_demes, mothers_by_site, n_per_deme):
    """Map (deme, local_idx) -> global haploid sample ID."""
    ids = []
    for deme, mothers in zip(site_demes, mothers_by_site):
        ids.extend(deme * n_per_deme + m for m in mothers)
    return ids


# ---------------------------------------------------------------------------
# Seed generation
# ---------------------------------------------------------------------------

def _generate_seed_pairs(maternal_ids, K, n_per_deme, seeds_per_mother,
                          pollen_pool, rng):
    """
    For each maternal haploid sample, generate seeds_per_mother diploid seeds
    via open pollination. Returns a flat list:
      [maternal_0, paternal_0, maternal_0, paternal_1, ..., maternal_1, ...]
    (seeds_per_mother consecutive pairs per mother, no self-pollination).
    """
    all_samples = np.arange(K * n_per_deme)
    pairs = []
    for maternal_id in maternal_ids:
        maternal_deme = int(maternal_id // n_per_deme)
        if pollen_pool == "local":
            pool = np.arange(maternal_deme * n_per_deme,
                              (maternal_deme + 1) * n_per_deme)
        else:
            pool = all_samples
        pool = pool[pool != maternal_id]  # no selfing
        for _ in range(seeds_per_mother):
            paternal_id = int(rng.choice(pool))
            pairs.extend([maternal_id, paternal_id])
    return pairs


def _build_seed_structure(seed_pairs):
    """
    Deduplicate seed_pairs for ts.simplify() and build the seed_structure
    mapping back to positions in the simplified tree sequence.

    Returns:
      unique_samples    : list of unique sample IDs (for ts.simplify)
      seed_structure    : list of (maternal_pos, paternal_pos) per seed
      maternal_positions: list of unique maternal positions in the simplified
                          ts, one per distinct mother tree (in collection order)
    """
    seen = {}
    for s in seed_pairs:
        if s not in seen:
            seen[s] = len(seen)
    unique_samples = list(seen.keys())

    seed_structure = [
        (seen[seed_pairs[i]], seen[seed_pairs[i + 1]])
        for i in range(0, len(seed_pairs), 2)
    ]

    # unique maternal positions: first element of each pair, deduplicated,
    # preserving collection order (one entry per distinct mother tree)
    mat_seen = {}
    for mat_pos, _ in seed_structure:
        if mat_pos not in mat_seen:
            mat_seen[mat_pos] = len(mat_seen)
    maternal_positions = list(mat_seen.keys())

    return unique_samples, seed_structure, maternal_positions


# ---------------------------------------------------------------------------
# Core sampling entry point
# ---------------------------------------------------------------------------

def sample_design(design: CollectionDesign, world_idx: int = 0,
                   groundtruth_dir: str = None,
                   K: int = None, n_per_deme: int = None,
                   n_chrom: int = None, mu: float = None):
    """
    Draw one collection-design realization from the ground truth, generating
    diploid seeds via open pollination.

    Returns
    -------
    dict with:
      'site_demes'     : which demes were selected as collection sites
      'mothers_by_site': local mother indices selected within each site
      'n_seeds'        : total seeds generated
      'seed_structure' : list of (maternal_pos, paternal_pos) per seed,
                         indexing into the simplified ts sample ordering.
                         Pass to stage2.diversity.compute_diversity().
      'chrom_ts'       : {chrom_idx: simplified+mutated TreeSequence}
                         Simplified to the unique set of haploid samples
                         involved (no duplicate samples).
      'design'         : CollectionDesign used (provenance)
      'world_idx'      : ground-truth world index
    """
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME
    n_chrom = n_chrom if n_chrom is not None else params.N_CHROM
    mu = mu if mu is not None else params.MU

    if groundtruth_dir is None:
        groundtruth_dir = os.path.join(DATADIR, "groundtruth")

    validate_design_feasibility(design, K=K, n_per_deme=n_per_deme)

    rng = np.random.default_rng(design.seed)

    site_demes = select_sites(K, design.n_sites, design.site_selection, rng)
    mothers_by_site = [
        select_mothers(n_per_deme, n_mothers, rng)
        for n_mothers in design.mothers_list()
    ]
    maternal_ids = _global_sample_ids(site_demes, mothers_by_site, n_per_deme)

    seed_pairs = _generate_seed_pairs(
        maternal_ids, K, n_per_deme,
        design.seeds_per_mother, design.pollen_pool, rng,
    )
    unique_samples, seed_structure, maternal_positions = _build_seed_structure(seed_pairs)
    n_seeds = len(seed_structure)

    chrom_ts = {}
    for c in range(n_chrom):
        path = os.path.join(groundtruth_dir,
                            f"world{world_idx:02d}_chr{c:02d}.trees")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"Ground truth file not found: {path} -- "
                f"has stage1/build_groundtruth.py finished world {world_idx}?"
            )
        ts = tskit.load(path)
        ts = ts.simplify(samples=unique_samples)
        ts = msprime.sim_mutations(ts, rate=mu,
                                    random_seed=design.seed + 1000 * c + 1)
        chrom_ts[c] = ts

    return {
        "site_demes": site_demes,
        "mothers_by_site": mothers_by_site,
        "n_seeds": n_seeds,
        "seed_structure": seed_structure,
        "maternal_positions": maternal_positions,
        "chrom_ts": chrom_ts,
        "design": design,
        "world_idx": world_idx,
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

    for c, ts in result["chrom_ts"].items():
        ts.dump(os.path.join(out, f"chr{c:02d}.trees"))

    manifest = {
        "design": asdict(result["design"]),
        "world_idx": result["world_idx"],
        "site_demes": result["site_demes"],
        "mothers_by_site": result["mothers_by_site"],
        "n_seeds": result["n_seeds"],
        "seed_structure": result["seed_structure"],
        "params_snapshot": {
            "K": params.K,
            "N_PER_DEME": params.N_PER_DEME,
            "N_CHROM": params.N_CHROM,
            "CHROM_LENGTH": params.CHROM_LENGTH,
            "MU": params.MU,
        },
    }
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    return out

