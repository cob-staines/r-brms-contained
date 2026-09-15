#!/bin/bash
# ============================================================
# Timing probe: runs the REAL model (same data, formula, family,
# adapt_delta) but with a small iter/warmup, to estimate how long
# the full run in run_size_bd_modeling.sbatch will take.
#
# Caveat: warmup (step-size/mass-matrix adaptation) is often slower
# per-iteration than post-warmup sampling, especially with a high
# adapt_delta like 0.98 here. A short probe still includes full
# warmup, so it's a reasonable estimate, but treat the extrapolation
# as a lower bound and pad it (e.g. x1.5) when setting -t on the real job.
#
# Submit with: sbatch hpc/timing_probe.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 1                    # must equal N_CORES * N_THREADS_PER_CHAIN below
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -J size_bd_timing_probe
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

# 1 chain, 100 warmup + 100 sampling iterations - swap in for the full
# 4 chains / 2000 warmup / 4000 iter production run once you have a number.
# Keep N_THREADS_PER_CHAIN matched to whatever run_size_bd_modeling.sbatch
# uses, since threading changes per-iteration wall time.
export N_CHAINS=1
export N_CORES=1
export N_THREADS_PER_CHAIN=1
export N_WARMUP=10
export N_ITER=20

singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" Rscript size_bd_modeling.R
