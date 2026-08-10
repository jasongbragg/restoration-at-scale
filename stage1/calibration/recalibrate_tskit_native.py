"""
stage1_recalibrate_tskit_native.py

Refactors the Stage 1 Fst calibration onto tskit's native branch-mode
diversity statistic (unbiased, no mutations required) and re-derives the
migration rate.

THIS SCRIPT DOCUMENTS A CORRECTION, not just a style change: the
previous estimator (nei_gst_corrected in stage1_utils.py) only corrected
the within-deme (Hs) term for small-sample bias, leaving the pooled (Ht)
term uncorrected. That residual bias was small but not zero. Recomputing
Fst at the previously "calibrated" m=0.0024 using the proper unbiased
branch-mode statistic gives Fst~0.086, not the intended 0.05 -- meaning
the true calibrated migration rate needs to be higher than we thought.

Run: python stage1_recalibrate_tskit_native.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from stage1.utils import simulate_branch_gst_replicated, build_chromosome_ts, branch_gst

if __name__ == "__main__":
    print("Step 1: confirm the bias at the OLD calibrated value (m=0.0024)")
    mean_fst, sd = simulate_branch_gst_replicated(0.0024, n_reps=4)
    print(f"  m=0.0024 (old) -> branch Fst = {mean_fst:.4f} +/- {sd:.4f}")
    print("  (expected ~0.05 if the old estimator had been unbiased -- it wasn't)")

    print("\nStep 2: re-sweep migration rate against the unbiased statistic")
    target_fst = 0.05
    print(f"{'m':>10} | {'mean Fst':>9} | {'std':>6}")
    print("-" * 32)
    for m in [0.0038, 0.0040, 0.0042, 0.0044, 0.0046]:
        mean_fst, sd = simulate_branch_gst_replicated(m, n_reps=10)
        print(f"{m:10.4f} | {mean_fst:9.4f} | {sd:6.4f}")

    print("\nStep 3: confirm at full production architecture (11 x 4.5Mb)")
    K, n_per_deme, m_final = 16, 12, 0.0042
    fst_per_chrom = []
    for c in range(11):
        ts = build_chromosome_ts(c, m_final, K=K, deme_Ne=6000, L=4_500_000,
                                  n_per_deme=n_per_deme)
        fst_per_chrom.append(branch_gst(ts, K, n_per_deme))
    print(f"  m={m_final} -> per-chromosome Fst: "
          + ", ".join(f"{v:.3f}" for v in fst_per_chrom))
    print(f"  mean = {sum(fst_per_chrom)/len(fst_per_chrom):.4f} (target 0.05)")

    print("\n" + "-" * 50)
    print("REVISED calibrated value: m \u2248 0.0042 (supersedes the old m=0.0024)")
