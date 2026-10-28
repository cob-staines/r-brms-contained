#!/bin/bash
# ============================================================
# Bd load cutoff runs: fits size_bd_modeling.R (m11d) once per Bd load
# cutoff as a Slurm array job, one cutoff per task, all in parallel.
#
# Submit with: sbatch hpc/run_size_bd_cutoffs.sbatch
#
# Array index -> cutoff (response bd_load_<cutoff>_density):
#   0 co0     no cutoff (any positive qPCR result)
#   1 co30
#   2 co100   main analysis
#   3 co300
#   4 co1000
# Rerun a subset by overriding --array on the command line:
#   sbatch --array=2 hpc/run_size_bd_cutoffs.sbatch      # co100 only
#   sbatch --array=0,4 hpc/run_size_bd_cutoffs.sbatch    # co0 and co1000
#
# Required input (gitignored - rsync it up first, see README):
#   DATA_FILE  data/bd_model_env_unified_2026-10-28.RData
# Override with --export=ALL,DATA_FILE=data/<file>.RData
#
# Output: models/unified_model_results_zln_m11d_<cutoff>_<timestamp>.rds
# plus a _with_loo.rds copy (the naming size_bd_unified/02 expects).
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0-4
#SBATCH -c 8                    # must equal N_CORES * N_THREADS_PER_CHAIN below
#SBATCH --mem=48G                # the co100 fit (Sep 2026) completed at these settings; check
                                  # `sacct -j <jobid> --format=JobID,MaxRSS,Elapsed` after this run
                                  # and tighten. co0 has the most nonzero responses (55% vs 51% for
                                  # co100), so it is the largest.
#SBATCH -t 5-00:00:00             # the co100 fit took ~56 h per chain (warmup + sampling) at
                                  # 4000 iter / 2000 warmup / 2 threads per chain
#SBATCH -J size_bd_cutoff
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models

module load singularity 2>/dev/null || true

CUTOFFS=(co0 co30 co100 co300 co1000)
CUTOFF="${CUTOFFS[$SLURM_ARRAY_TASK_ID]}"

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "r-brms.sif not found: run hpc/pull_image.sbatch first" >&2
  exit 1
fi

export DATA_FILE="${DATA_FILE:-data/bd_model_env_unified_2026-10-28.RData}"
export BD_LOAD_VAR="bd_load_${CUTOFF}_density"
export RUN_TAG="m11d_${CUTOFF}"

# N_CORES * N_THREADS_PER_CHAIN must equal -c above.
export N_CHAINS=4
export N_ITER=4000
export N_WARMUP=2000
export N_CORES=4
export N_THREADS_PER_CHAIN=2

echo "Task $SLURM_ARRAY_TASK_ID: cutoff $CUTOFF ($BD_LOAD_VAR), data $DATA_FILE, run tag $RUN_TAG"

# stdbuf forces line-buffered stdout/stderr so progress flushes to the log
# file live - lets you `tail -f` this task's log to see iteration pace.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_modeling.R
