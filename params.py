"""
params.py  --  SINGLE SOURCE OF TRUTH for all locked simulation parameters.

Every script in this repository imports from here rather than defining its
own constants. This exists because scattered per-script constants caused two
real bugs during development:
  1. build_chromosome_ts()'s r default was still 1e-8 (the calibration
     placeholder) when the production harness was first written, so an
     early draft would have silently built the ground truth at the wrong
     recombination rate.
  2. The ancestral population was sized at DEME_NE (6000) rather than
     ANCESTRAL_NE (101,053), which silently truncated total diversity --
     caught only by the integration-level validation script.

Both would have been impossible if there had always been one place where
parameters live. Any future parameter change (e.g. n_per_deme, m, r)
should be made HERE and nowhere else.

Cross-reference: PROJECT_SUMMARY.md Section 5e for the full calibration
derivation of each value.
"""

# ---------------------------------------------------------------------------
# Metapopulation structure
# ---------------------------------------------------------------------------
K = 16                    # number of demes, linear stepping-stone chain
DEME_NE = 6000            # per-deme effective size (theory: K*deme_Ne/(1-Fst)
                          # = global Ne target -- see PROJECT_SUMMARY.md Sec 5)
ANCESTRAL_NE = 101_053    # global Ne at the ancestral merge -- MUST NOT be
                          # DEME_NE (that was a real bug caught by validation)
ANCESTRAL_MERGE_TIME = 30_000  # generations; collapses migration structure
                                # beyond this depth, making full-length
                                # chromosomes tractable (~2x speedup)

# ---------------------------------------------------------------------------
# Calibrated population-genetic parameters
# (all calibrated empirically against Guo et al. 2026, M. quinquenervia NSW)
# ---------------------------------------------------------------------------
M = 0.0042        # migration rate/generation -- calibrated against Fst=0.05
                  # using tskit branch-mode diversity (unbiased estimator)
R = 5e-7          # recombination rate/bp/generation -- calibrated against
                  # LD half-decay at 711 bp (Guo et al. 2026 Fig 2)
MU = 1e-8         # mutation rate/bp/generation -- from Guo et al. 2026

# ---------------------------------------------------------------------------
# Genome architecture
# ---------------------------------------------------------------------------
N_CHROM = 11              # real M. quinquenervia karyotype
CHROM_LENGTH = 4_500_000  # bp per chromosome (chosen for ROH resolution)

# ---------------------------------------------------------------------------
# Ground truth sampling
# ---------------------------------------------------------------------------
N_PER_DEME = 32           # individuals simulated per deme in the ground truth.
                          # This is the hard ceiling on mothers-per-site for
                          # any Stage 2 collection design -- if a design needs
                          # more than this, rebuild the ground truth first.
                          # Was 12 in the initial build (sized for Fst/diversity
                          # calibration only); increased to 32 to support the
                          # full field-realistic design grid.

# ---------------------------------------------------------------------------
# Replicate worlds
# ---------------------------------------------------------------------------
REPLICATE_WORLDS = 10     # independent ground-truth builds (different random
                          # seeds, same parameters) for robustness checking.
                          # One world ~= 115 min on the server at full scale.
