#!/bin/bash
# ============================================================
# Stage 1 contrast models (size_bd_contrast): fits
# size_bd_contrast_modeling.R via brms/cmdstanr as a Slurm array job,
# one job (species x variant) per task.
#
# Array index -> job: the array_index column of data/01_contrast_jobs.csv
# (alphabetical by species then variant, 0-based; printed at the top of
# each task's log). 32 jobs (0-31) as of the 2026-11-24 tiers: 31 species
# plus the a. canorus tier A sensitivity fit (3). If the jobs file changes,
# update the array indices in the groups below.
#
# To submit all four groups at once: bash hpc/submit_size_bd_contrast.sh
#
# Resource groups (sbatch options on the command line override the #SBATCH
# defaults below, which are the small group). Estimates from the local timing
# test in size_bd_contrast/02_contrast_modeling.qmd (4 chains, 300 warmup +
# 100 sampling, scaled to 2000 + 2000; slowest chain; at adapt_delta 0.95):
#   s. flotator      n 0.6k   A x B   est 0.13 h   (1 divergence in the test)
#   r. catesbeiana   n 1.8k   A x B   est 0.51 h
#   r. clamitans     n 2.6k   A x A   est 0.70 h
# Across these, time grows ~n^1.2; other jobs are extrapolated with n^1.2
# (likely) to n^1.6 (upper bound for the time limits). ADAPT_DELTA defaults
# to 0.99 in the R script (~1.5-2x slower than 0.95; the estimates below
# include that); r. sierrae runs at 0.95.
#
#   small (n < 1k, est <= ~0.5 h at 0.99): 23 jobs, 4 h / 8G (the defaults)
#     sbatch hpc/run_size_bd_contrast_modeling.sbatch
#   intermediate (n 1-3.4k, est ~0.5-2 h at 0.99): a. canorus main (2) and
#     sens_tierA (3), p. crucifer (16), r. cascadae (19), r. catesbeiana (20),
#     r. clamitans (21)
#     sbatch --array=2,3,16,19-21 --mem=16G -t 08:00:00 hpc/run_size_bd_contrast_modeling.sbatch
#   large (n 5.3-6.8k, est ~3-6.5 h at 0.99): n. viridescens (15), r. muscosa (23)
#     sbatch --array=15,23 --mem=24G -t 16:00:00 hpc/run_size_bd_contrast_modeling.sbatch
#   sierrae (n 26.6k, A x A, est ~11-29 h at 0.95): r. sierrae (25)
#     sbatch --array=25 --mem=48G -t 2-00:00:00 --export=ALL,ADAPT_DELTA=0.95 hpc/run_size_bd_contrast_modeling.sbatch
#
# Divergences: check the convergence table in 02_contrast_modeling.qmd.
#   - r. sierrae diverging at 0.95: resubmit at 0.99 with more time, e.g.
#       sbatch --array=25 --mem=48G -t 3-00:00:00 --export=ALL,ADAPT_DELTA=0.99 hpc/run_size_bd_contrast_modeling.sbatch
#   - a job still diverging at 0.99 points to the model (e.g. a contrast sd
#     with too few co-sampled levels to inform it), not the sampler: look at
#     the model (tier) before raising adapt_delta further.
# A resubmission overwrites that job's model files in OUT_DIR. Record any
# per-species change in contrast_adapt_delta() in 02_contrast_modeling.qmd.
#
# Memory: not measured. Check actual use and trim with
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
#
# Required inputs (gitignored - rsync them up first, see README):
#   DATA_FILE     data/bd_model_env_unified_2026-11-12.RData (the m11g data file)
#   JOBS_FILE     data/01_contrast_jobs.csv (from size_bd_contrast/intermediates/)
# Override with --export=ALL,DATA_FILE=...,JOBS_FILE=... (also ADAPT_DELTA)
# Output goes to models/size_bd_contrast/ (OUT_DIR). If a task times out
# during LOO, the raw fit is already saved as *_no_loo.rds.
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0,1,4-14,17,18,22,24,26-31   # small group; see above
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=8G                # small group; see above
#SBATCH -t 04:00:00             # small group (est <= ~0.5 h); see above
#SBATCH -J size_bd_contrast
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models/size_bd_contrast

module load singularity 2>/dev/null || true

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "Pulling container image..."
  singularity pull "$IMAGE" docker://ghcr.io/cob-staines/r-brms-contained:latest
fi

# N_CORES * N_THREADS_PER_CHAIN equals -c: threads per chain follow -c, so
# `sbatch -c 8 ...` on the command line runs 2 threads per chain.
export N_CORES=4
export N_THREADS_PER_CHAIN=$(( ${SLURM_CPUS_PER_TASK:-4} / N_CORES ))

# stdbuf forces line-buffered stdout/stderr so progress flushes to the log
# file live - lets you `tail -f` this task's log to see iteration pace.
# SLURM_ARRAY_TASK_ID is passed through to pick this task's job.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_contrast_modeling.R
