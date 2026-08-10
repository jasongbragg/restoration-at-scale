"""
stage1/build_groundtruth.py

Production script for building the Stage 1 ground-truth ancestry.
All parameters are imported from params.py at the repo root -- do not
define constants here. If you need to change a parameter, change it in
params.py and it will propagate everywhere automatically.

Confirmed timing on real server hardware: 6,886.1s (~115 min) per
full-length chromosome at the locked parameters. Processes one world at
a time (11 chromosomes in parallel within a world, then moves to the
next world). Designed to run unattended for days:
  - already-completed chromosomes are skipped on restart (resumable)
  - each chromosome is saved immediately on completion

Usage:
    nohup python3 stage1/build_groundtruth.py > build.log 2>&1 &
    # or: tmux new -s groundtruth 'python3 stage1/build_groundtruth.py'

Output: data/groundtruth/world00_chr00.trees ... world{N}_chr10.trees
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import params
from stage1.utils import build_chromosome_ts, branch_gst

OUTDIR = os.path.join(os.path.dirname(__file__), "..", "data", "groundtruth")


def _job(args):
    chrom_idx, world_idx = args
    seed = 1 + 100_000 * world_idx
    ts = build_chromosome_ts(chrom_idx, seed=seed)
    path = os.path.join(OUTDIR, f"world{world_idx:02d}_chr{chrom_idx:02d}.trees")
    ts.dump(path)
    fst = branch_gst(ts)
    return chrom_idx, world_idx, path, fst


def build_one_world(world_idx):
    """Build all N_CHROM chromosomes for one world in parallel.
    Skips chromosomes whose output file already exists."""
    todo = []
    for c in range(params.N_CHROM):
        path = os.path.join(OUTDIR, f"world{world_idx:02d}_chr{c:02d}.trees")
        if not os.path.exists(path):
            todo.append((c, world_idx))

    if not todo:
        print(f"world {world_idx:02d}: already complete, skipping")
        return

    print(f"world {world_idx:02d}: {len(todo)}/{params.N_CHROM} chromosome(s) to build "
          f"(n_per_deme={params.N_PER_DEME})")
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=min(len(todo), os.cpu_count() or 1)) as pool:
        futures = [pool.submit(_job, job) for job in todo]
        for fut in as_completed(futures):
            chrom_idx, w, path, fst = fut.result()
            print(f"  done: world{w:02d} chr{chrom_idx:02d} -> {os.path.basename(path)} "
                  f"(Fst={fst:.4f})")
            results.append(fst)
    elapsed = time.time() - t0
    print(f"world {world_idx:02d}: done in {elapsed:.1f}s ({elapsed/60:.1f} min) -- "
          f"Fst range [{min(results):.3f}, {max(results):.3f}], "
          f"mean={sum(results)/len(results):.4f} (target 0.05)")


if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"Locked parameters from params.py:")
    print(f"  K={params.K}, deme_Ne={params.DEME_NE}, ancestral_Ne={params.ANCESTRAL_NE}")
    print(f"  m={params.M}, r={params.R}, mu={params.MU}")
    print(f"  {params.N_CHROM} x {params.CHROM_LENGTH}bp, n_per_deme={params.N_PER_DEME}")
    print(f"  ancestral_merge_time={params.ANCESTRAL_MERGE_TIME}")
    print(f"Working through {params.REPLICATE_WORLDS} world(s) "
          f"(os.cpu_count()={os.cpu_count()})\n")

    overall_t0 = time.time()
    for world_idx in range(params.REPLICATE_WORLDS):
        build_one_world(world_idx)

    print(f"\nAll done. Total elapsed: {(time.time()-overall_t0)/60:.1f} min.")
    print("Next: run stage1/validate_groundtruth.py")
