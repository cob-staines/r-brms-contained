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
# Submission scripts: hpc/submit_size_bd_contrast.sh (round 1),
# hpc/submit_size_bd_contrast_round2.sh, hpc/submit_size_bd_contrast_round3.sh
# (current: all 32 jobs).
#
# Resource groups (sbatch options on the command line override the #SBATCH
# defaults below, which are the small group). Measured with sacct in rounds 1-2
# (MaxRSS, elapsed; jobs 1326767-1327598). Round 3 runs every job at 2000
# warmup + 4000 sampling (16000 draws), which roughly doubles the draw-related
# memory (r. clamitans 2.5 -> 5.2 GB, r. sierrae 16.5 -> 34 GB from 8000 to
# 16000 draws) and adds ~1.5x run time. Requests below allow ~1.3-2x the
# expected peak WITH loo (RUN_LOO=TRUE), so they hold either way; loo is off by
# default. ADAPT_DELTA 0.99 (r. sierrae 0.95).
#   small (n < 1k): 23 jobs, peak 1.9-2.4 GB, <= 10 min at 8000 draws
#     -> 1 h / 4G (the defaults)
#   intermediate (n 1-3.4k): a. canorus main (2) and sens_tierA (3),
#     p. crucifer (16), r. cascadae (19), r. catesbeiana (20), r. clamitans (21);
#     peak <= 5.2 GB, <= 29 min at 16000 draws (r. clamitans)
#     -> 2 h / 8G
#   large (n 5.3-6.8k): n. viridescens (15; 3.0 GB, 24 min at 8000 draws, now
#     tier D), r. muscosa (23; 7.7 GB, 2.0 h at 16000 draws)
#     -> 5 h / 12G
#   sierrae (n 26.6k): r. sierrae (25); 34 GB with loo (less without it, not yet
#     measured), 5.9 h at 16000 draws
#     -> 16 h / 40G
#
# Divergences: check the convergence table in 02_contrast_modeling.qmd.
#   - r. sierrae diverging at 0.95: resubmit at 0.99 with more time, e.g.
#       sbatch --array=25 --mem=40G -t 1-08:00:00 --export=ALL,ADAPT_DELTA=0.99 hpc/run_size_bd_contrast_modeling.sbatch
#   - a job still diverging at 0.99 points to the model (e.g. a contrast sd
#     with too few co-sampled levels to inform it), not the sampler: look at
#     the model (tier) before raising adapt_delta further.
# A resubmission overwrites that job's model files in OUT_DIR. Record any
# per-species change in contrast_adapt_delta() in 02_contrast_modeling.qmd.
#
# Memory: measured in rounds 1-2 (above). Check a new round with
#   sacct -j <jobid> --format=JobID%22,Elapsed,MaxRSS,ReqMem | grep -E "batch|JobID"
#
# Required inputs (gitignored - rsync them up first, see README):
#   DATA_FILE     data/bd_model_env_unified_2026-11-12.RData (the m11g data file)
#   JOBS_FILE     data/01_contrast_jobs.csv (from size_bd_contrast/intermediates/)
# Override with --export=ALL,DATA_FILE=...,JOBS_FILE=... (also ADAPT_DELTA)
# Output goes to models/size_bd_contrast/ (OUT_DIR). LOO is off by default
# (RUN_LOO=FALSE; not used by the analysis). With RUN_LOO=TRUE the raw fit is
# saved first as *_no_loo.rds, so a task timing out during LOO keeps its fit.
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0,1,4-14,17,18,22,24,26-31   # small group; see above
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=4G                # small group (peak <= 2.4 GB); see above
#SBATCH -t 01:00:00             # small group (<= 10 min at 8000 draws); see above
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
