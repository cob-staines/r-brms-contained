#!/bin/bash
# ============================================================
# Round 3 of the stage 1 contrast models (size_bd_contrast): refits ALL 32
# jobs, because the event definition changed for every model (every tier has
# an event-level intercept):
#   - event = survey run (dates at a site with no gap > 5 days; runs > 14
#     days split at their longest gap), replacing 14-day bins of day of year
#   - consolidated contrast levels rule: tier D replaces "fixed";
#     n. viridescens (15) C -> D, l. chiricahuensis (12) prevalence A -> C
#   - every job at N_ITER 6000 (2000 warmup + 4000 sampling per chain, the R
#     script's new default; round 2 used 6000 only for the three slow-mixing
#     large species)
#   - LOO off by default (RUN_LOO=FALSE in the R script; not used)
# Memory and time from sacct in rounds 1-2 (see the job script header):
# requests allow ~1.3-2x the expected round 3 peak, with loo or without.
# Decisions: size_bd_contrast/README.md ("Definitions", "Contrast levels").
#
# Run from the repo root:
#   bash hpc/submit_size_bd_contrast_round3.sh
#
# Before submitting:
#   - pull this repo: size_bd_contrast_modeling.R must have the survey-run
#     add_contrast_groups() and accept tier "D"
#   - data/01_contrast_jobs.csv must be the round 3 file (from
#     01_contrast_tiers.Rmd; it has the d_reason_hu column). array indices
#     are unchanged from rounds 1-2
#   - refits overwrite the round 2 model files in models/size_bd_contrast/:
#     copy those somewhere first if you want them for comparison
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_contrast_modeling.sbatch

for f in data/bd_model_env_unified_2026-11-12.RData data/01_contrast_jobs.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done
if ! head -1 data/01_contrast_jobs.csv | grep -q d_reason_hu; then
  echo "data/01_contrast_jobs.csv is not the round 3 file (no d_reason_hu column): copy the new one" >&2
  exit 1
fi
if ! grep -q 'Sys.getenv("N_ITER", "6000")' size_bd_contrast_modeling.R; then
  echo "size_bd_contrast_modeling.R does not default to N_ITER 6000: git pull first" >&2
  exit 1
fi
if ! grep -q "event_max_gap" size_bd_contrast_modeling.R; then
  echo "size_bd_contrast_modeling.R still has the 14-day bin events: git pull first" >&2
  exit 1
fi

# small: 23 jobs, n < 1k, incl. l. chiricahuensis (12). job script defaults
# (1 h, 4G, adapt_delta 0.99): peak <= 2.4 GB, <= 10 min at 8000 draws
sbatch "$JOB"

# intermediate: a. canorus main (2) and sens_tierA (3), p. crucifer (16),
# r. cascadae (19), r. catesbeiana (20), r. clamitans (21); n 1-3.4k.
# peak <= 5.2 GB and <= 29 min at 16000 draws (r. clamitans)
sbatch --array=2,3,16,19-21 --mem=8G -t 02:00:00 "$JOB"

# large: n. viridescens (15; now tier D, faster than round 2's 24 min),
# r. muscosa (23; 7.7 GB, 2.0 h at 16000 draws; borderline Rhat 1.012 in round 2)
sbatch --array=15,23 --mem=12G -t 05:00:00 "$JOB"

# r. sierrae (25), n 26.6k, at adapt_delta 0.95 (0.99 would roughly double the
# run; no divergences at 0.95 in round 2). 34 GB with loo, 5.9 h at 16000 draws
sbatch --array=25 --mem=40G -t 16:00:00 --export=ALL,ADAPT_DELTA=0.95 "$JOB"

echo
echo "Submitted round 3 (32 jobs). Check status with: squeue -u \$USER"
