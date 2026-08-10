"""
stage2/sampling.py

Stage 2: draw collection-design samples from the Stage 1 ground truth.

Deliberately parameterized rather than tuned to specific numbers -- the
actual site/mothers/N grid is a separate decision (see
docs/PROJECT_SUMMARY.md Section 9). This module should work for ANY
design within the ground truth's capacity and fail loudly for anything
that exceeds it.

Key constraint: the ground truth's N_PER_DEME (currently 32, set in
params.py) is the hard ceiling on mothers-per-site. If a design needs
more, rebuild the ground truth first.
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


@dataclass
class CollectionDesign:
    """One collection-design sampling regime.

    n_sites: number of demes to treat as collection sites.
    mothers_per_site: single int (same count at every site) or a list of
        length n_sites for unequal designs.
    site_selection: 'random', 'even_spacing', or 'contiguous'.
        Which is most realistic for real collection patterns is an open
        question -- all three are supported so they can be compared later.
    seed: reproducibility. Same design + seed always gives the same
        sites and mothers.
    relatedness_mode: only 'unrelated' (independent draws, matching the
        real >=20m spacing protocol) is implemented. Placeholder for the
        deferred kinship-imposed comparison (PROJECT_SUMMARY.md Sec 3).
    """
    n_sites: int
    mothers_per_site: Union[int, list]
    site_selection: str = "random"
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

    def total_n(self):
        return sum(self.mothers_list())

    def label(self):
        """Short human-readable label for use as a directory name."""
        mps = self.mothers_per_site if isinstance(self.mothers_per_site, int) \
            else "unequal"
        return (f"{self.site_selection[:4]}_s{self.n_sites}_m{mps}"
                f"_N{self.total_n()}_seed{self.seed}")


def validate_design_feasibility(design: CollectionDesign, K=None, n_per_deme=None):
    """Raise ValueError with specific, actionable message if the design
    can't be satisfied from the current ground truth. Called before
    any tree-sequence operations -- don't skip even for 'obviously fine'
    designs, since it's cheap and failing early beats failing cryptically."""
    K = K if K is not None else params.K
    n_per_deme = n_per_deme if n_per_deme is not None else params.N_PER_DEME

    if design.n_sites > K:
        raise ValueError(
            f"Design requests n_sites={design.n_sites} but ground truth has "
            f"only K={K} demes. Reduce n_sites, or rebuild with more demes."
        )
    if design.n_sites < 1:
        raise ValueError(f"n_sites must be >=1, got {design.n_sites}")

    for i, n in enumerate(design.mothers_list()):
        if n > n_per_deme:
            raise ValueError(
                f"Design requests {n} mothers at site {i}, but the ground "
                f"truth was built with N_PER_DEME={n_per_deme} -- that's the "
                f"hard ceiling, not a sampling-code limitation. Either reduce "
                f"mothers_per_site to <={n_per_deme}, or increase N_PER_DEME "
                f"in params.py and rebuild the ground truth."
            )
        if n < 1:
            raise ValueError(f"mothers_per_site at site {i} must be >=1, got {n}")

    if design.relatedness_mode != "unrelated":
        raise NotImplementedError(
            f"relatedness_mode='{design.relatedness_mode}' not implemented. "
            f"Only 'unrelated' exists so far."
        )


def select_sites(K, n_sites, strategy, rng):
    """Return a sorted list of n_sites deme indices from range(K)."""
    if strategy == "random":
        return sorted(rng.choice(K, size=n_sites, replace=False).tolist())
    elif strategy == "even_spacing":
        idx = np.unique(np.round(np.linspace(0, K - 1, n_sites)).astype(int))
        if len(idx) < n_sites:
            raise ValueError(
                f"even_spacing with K={K}, n_sites={n_sites} produced only "
                f"{len(idx)} unique demes. Try 'random' for this combination."
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
    """Without-replacement draw of n_mothers local indices from range(n_per_deme)."""
    return sorted(rng.choice(n_per_deme, size=n_mothers, replace=False).tolist())


def _global_sample_ids(site_demes, mothers_by_site, n_per_deme):
    """Map (deme, local_idx) pairs to global sample IDs as assigned by
    msprime when the ground truth was built with deme-contiguous blocks
    (i.e. deme i owns indices [i*n_per_deme, (i+1)*n_per_deme) )."""
    ids = []
    for deme, mothers in zip(site_demes, mothers_by_site):
        ids.extend(deme * n_per_deme + m for m in mothers)
    return ids


def sample_design(design: CollectionDesign, world_idx: int = 0,
                   groundtruth_dir: str = None,
                   K: int = None, n_per_deme: int = None,
                   n_chrom: int = None, mu: float = None):
    """Draw one realization of a collection design from the ground truth.

    Returns a dict with:
      'site_demes'     : which demes were selected as collection sites
      'mothers_by_site': local mother indices selected within each site
      'chrom_ts'       : {chrom_idx: simplified+mutated TreeSequence}
      'design'         : the CollectionDesign used (provenance)
      'world_idx'      : which ground-truth world this came from

    The ground truth tree sequences are stored without mutations (ancestry
    only). Mutations are overlaid here at the per-chromosome level, with
    a seed derived from the design seed so results are fully reproducible.
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
    sample_ids = _global_sample_ids(site_demes, mothers_by_site, n_per_deme)

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
        ts = ts.simplify(samples=sample_ids)
        ts = msprime.sim_mutations(ts, rate=mu,
                                    random_seed=design.seed + 1000 * c + 1)
        chrom_ts[c] = ts

    return {
        "site_demes": site_demes,
        "mothers_by_site": mothers_by_site,
        "chrom_ts": chrom_ts,
        "design": design,
        "world_idx": world_idx,
    }


def save_sampled_design(result: dict, label: str = None,
                         designs_dir: str = None):
    """Save each chromosome's simplified+mutated tree sequence plus a JSON
    manifest recording full provenance (exact demes and individuals selected,
    not just the seed)."""
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
        "realized_total_n": sum(len(m) for m in result["mothers_by_site"]),
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
