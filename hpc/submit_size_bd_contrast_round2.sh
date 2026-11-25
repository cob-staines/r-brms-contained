#!/bin/bash
# ============================================================
# Round 2 of the stage 1 contrast models (size_bd_contrast): refits only
# the jobs changed by the model structure rules R1-R3, plus longer runs
# for the three large species that mixed slowly in round 1. The other 10
# jobs keep their round 1 fits. Decisions: size_bd_contrast/02_contrast_modeling.qmd
# ("Fit rounds") and README "Model structure rules".
#
# Run from the repo root:
#   bash hpc/submit_size_bd_contrast_round2.sh
#
# Before submitting: data/01_contrast_jobs.csv must be the NEW jobs file
# (from 01_contrast_tiers.Rmd after R1-R3; it has the sigma_by_population
# column). Array indices are unchanged from round 1. Refits overwrite the
# round 1 model files in models/size_bd_contrast/ - copy those somewhere
# first if you want to keep them for comparison.
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_contrast_modeling.sbatch

for f in data/bd_model_env_unified_2026-11-12.RData data/01_contrast_jobs.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done
if ! head -1 data/01_contrast_jobs.csv | grep -q sigma_by_population; then
  echo "data/01_contrast_jobs.csv is the round 1 file (no structure columns): copy the new one" >&2
  exit 1
fi

# structure changed (R1-R3), small: 17 jobs, n < 1k. job script defaults (2 h, 8G, adapt_delta 0.99)
sbatch --array=0,1,4-6,8-12,17,22,26-29,31 "$JOB"

# structure changed, intermediate: a. canorus main (2; load now "fixed")
sbatch --array=2 --mem=16G -t 04:00:00 "$JOB"

# structure changed, large: n. viridescens (15; site-year and event levels dropped, so faster
# than round 1)
sbatch --array=15 --mem=24G -t 08:00:00 "$JOB"

# same model, more draws (round 1: Rhat 1.010-1.023, no divergences). N_ITER=6000 with the
# default 2000 warmup = 4000 draws per chain (2x), ~1.5x the run time. check round 1 times with
#   sacct -j <round 1 jobid> --format=JobID,Elapsed,MaxRSS
# and adjust -t if needed.
#   r. clamitans (21), r. muscosa (23)
sbatch --array=21 --mem=16G -t 06:00:00 --export=ALL,N_ITER=6000 "$JOB"
sbatch --array=23 --mem=24G -t 12:00:00 --export=ALL,N_ITER=6000 "$JOB"
#   r. sierrae (25), still at adapt_delta 0.95
sbatch --array=25 --mem=48G -t 1-12:00:00 --export=ALL,ADAPT_DELTA=0.95,N_ITER=6000 "$JOB"

echo
echo "Submitted round 2. Check status with: squeue -u \$USER"
