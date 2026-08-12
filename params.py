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
parameters live. Any future parameter change should be made HERE and nowhere
else.

Cross-reference: docs/PROJECT_SUMMARY.md Section 5e for the full calibration
derivation of each value.
"""

# ---------------------------------------------------------------------------
# Metapopulation structure
# ---------------------------------------------------------------------------
K = 16                    # number of demes, linear stepping-stone chain
DEME_NE = 6000            # per-deme effective size
ANCESTRAL_NE = 101_053    # global Ne at the ancestral merge -- MUST NOT be
                          # DEME_NE (that was a real bug caught by validation)
ANCESTRAL_MERGE_TIME = 30_000  # generations

# ---------------------------------------------------------------------------
# Calibrated population-genetic parameters
# ---------------------------------------------------------------------------
M = 0.0042        # migration rate/generation -- calibrated against Fst=0.05
R = 5e-7          # recombination rate/bp/generation -- calibrated against
                  # LD half-decay at 711 bp (Guo et al. 2026 Fig 2)
MU = 1e-8         # mutation rate/bp/generation -- from Guo et al. 2026

# ---------------------------------------------------------------------------
# Genome architecture
# ---------------------------------------------------------------------------
N_CHROM = 11              # real M. quinquenervia karyotype
CHROM_LENGTH = 4_500_000  # bp per chromosome

# ---------------------------------------------------------------------------
# Ground truth sampling
# ---------------------------------------------------------------------------
N_PER_DEME = 64           # HAPLOID samples simulated per deme in the ground
                          # truth. Represents N_PER_DEME/2 = 32 DIPLOID mother
                          # trees per site (each mother = 2 consecutive haploid
                          # samples). Must be even. Increased from 32 to 64 to
                          # support the full design table (max 32 diploid
                          # mothers/site requires 64 haploid samples/deme).

# ---------------------------------------------------------------------------
# Replicate worlds
# ---------------------------------------------------------------------------
REPLICATE_WORLDS = 10
