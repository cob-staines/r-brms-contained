#!/bin/bash
# ============================================================
# THREADED timing/memory probe: identical to hpc/timing_probe.sbatch
# (baseline, N_THREADS_PER_CHAIN=1) except N_THREADS_PER_CHAIN=2 here.
# Compare this job's "Fit time" line against the baseline's to measure
# whether within-chain reduce_sum threading is actually worth using for
# this model/data - it parallelizes the likelihood sum across the ~54k
# observations within each chain's gradient evaluations, but does NOT
# reduce how many gradient evaluations NUTS needs per iteration (that's
# governed by tree depth), so the real speedup is model/data-dependent
# and worth measuring rather than assuming from general guidance.
#
# Submit with: sbatch hpc/timing_probe_threaded.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 8                    # must equal N_CORES * N_THREADS_PER_CHAIN below
#SBATCH --mem=64G                # generous headroom for measuring, not a final estimate - a 16GB
                                  # local VM already crashed at full production settings
#SBATCH -t 12:00:00              # TODO: unverified guess - bump further if you still hit timeouts
#SBATCH -J size_bd_timing_probe_threaded
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail
mkdir -p logs models

module load singularity 2>/dev/null || true

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "Pulling container image..."
  singularity pull "$IMAGE" docker://ghcr.io/cob-staines/r-brms-contained:latest
fi

# 4 chains x 2 threads/chain = 8 total cores. Same warmup=300/iter=400 as
# the baseline probe - keep these matched for a fair comparison.
export N_CHAINS=4
export N_CORES=4
export N_THREADS_PER_CHAIN=2
export N_WARMUP=300
export N_ITER=400
export RUN_LOO=FALSE   # not meaningful at this small a draw count, and skips loo's own RAM/time cost here

# stdbuf forces line-buffered stdout/stderr so progress flushes to the log
# file live instead of sitting in a block buffer until the process exits -
# lets you `tail -f` a still-running job instead of waiting blind.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_modeling.R
