#!/bin/bash
# ============================================================
# Submit the stage 1 contrast models (size_bd_contrast) in their four
# resource groups. Each group is a separate submission of the same job
# script (run_size_bd_contrast_modeling.sbatch): command-line options
# override that script's #SBATCH defaults, so each group requests only the
# memory and time its jobs need (see the script's header for the estimates).
#
# Run from the repo root (the job script expects data/, logs/, models/
# relative to it):
#   bash hpc/submit_size_bd_contrast.sh
#
# Array indices = array_index in data/01_contrast_jobs.csv (32 jobs,
# 2026-11-24 tiers). If the jobs file changes, update the groups below.
#
# Before submitting: data/bd_model_env_unified_2026-11-12.RData and
# data/01_contrast_jobs.csv must exist, and r-brms.sif should already be
# pulled (hpc/pull_image.sbatch), so the array tasks don't race to pull it.
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_contrast_modeling.sbatch

# fail early if inputs are missing (otherwise every array task fails on its own)
for f in data/bd_model_env_unified_2026-11-12.RData data/01_contrast_jobs.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done

# small: 23 jobs, n < 1k. uses the job script's defaults
# (2 h, 8G, ADAPT_DELTA 0.99 from the R script)
sbatch "$JOB"

# intermediate: a. canorus main (2) and sens_tierA (3), p. crucifer (16),
# r. cascadae (19), r. catesbeiana (20), r. clamitans (21); n 1-3.4k
sbatch --array=2,3,16,19-21 --mem=16G -t 04:00:00 "$JOB"

# large: n. viridescens (15), r. muscosa (23); n 5.3-6.8k
sbatch --array=15,23 --mem=24G -t 08:00:00 "$JOB"

# sierrae: r. sierrae (25), n 26.6k, at adapt_delta 0.95 (see the job script)
sbatch --array=25 --mem=48G -t 1-00:00:00 --export=ALL,ADAPT_DELTA=0.95 "$JOB"

echo
echo "Submitted. Check status with: squeue -u \$USER"
