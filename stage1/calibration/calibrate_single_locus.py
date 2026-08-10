"""
stage1_calibrate_single_locus.py

Reproduces the single-locus migration-rate calibration: K=8 vs K=16
linear stepping-stone chains, tuned against the real Fst=0.05 target
(Guo et al. 2026, M. quinquenervia NSW range).

Run: python stage1_calibrate_single_locus.py
"""

from stage1_utils import simulate_fst_replicated

if __name__ == "__main__":
    target_fst = 0.05
    print(f"Target Fst (Guo et al. 2026, NSW M. quinquenervia): {target_fst}")

    print("\n-- K=8 chain --")
    print(f"{'m':>10} | {'mean Fst':>9} | {'std':>6}")
    print("-" * 32)
    for m in [0.001, 0.0015, 0.002, 0.0025, 0.003]:
        mean_fst, sd = simulate_fst_replicated(m, K=8)
        print(f"{m:10.4f} | {mean_fst:9.4f} | {sd:6.4f}")

    print("\n-- K=16 chain (same deme_Ne, same n_per_deme) --")
    print(f"{'m':>10} | {'mean Fst':>9} | {'std':>6}")
    print("-" * 32)
    for m in [0.0018, 0.002, 0.0022, 0.0024, 0.0025, 0.003]:
        mean_fst, sd = simulate_fst_replicated(m, K=16)
        print(f"{m:10.4f} | {mean_fst:9.4f} | {sd:6.4f}")

    print("-" * 32)
    print("Calibrated results (deme_Ne=6000, n_per_deme=12, single 300kb locus):")
    print("  K=8  -> m \u2248 0.0010  (Fst \u2248 0.053 \u00b1 0.012)")
    print("  K=16 -> m \u2248 0.0024  (Fst \u2248 0.050 \u00b1 0.005)")
