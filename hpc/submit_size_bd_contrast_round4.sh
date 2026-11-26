#!/bin/bash
# ============================================================
# Round 4 of the stage 1 contrast models (size_bd_contrast): refits the four
# jobs whose model changes when R3 (sigma by population) is raised from
# >= 3 sites with >= 10 positive swabs to >= 3 sites with >= 20:
#   p. regilla (17), s. albomaculata (28), s. flotator (30), s. sila (31)
# all now sigma ~ 1. in round 3 the site-level sd of sigma was uninformed in
# all four, and s. flotator had 25 divergences in its funnel. no tier changes;
# the other 28 jobs keep their round 3 fits.
# Decisions: size_bd_contrast/README.md ("Nuisance structure", decision log).
#
# Run from the repo root:
#   bash hpc/submit_size_bd_contrast_round4.sh
#
# Before submitting:
#   - data/01_contrast_jobs.csv must be the round 4 file (from
#     01_contrast_tiers.Rmd with min_pos_sigma_pop = 20): checked below on
#     s. flotator's sigma_by_population. array indices are unchanged
#   - no change to size_bd_contrast_modeling.R is needed: it reads
#     sigma_by_population from the jobs file
#   - refits overwrite the round 3 model files of these four jobs in
#     models/size_bd_contrast/: copy them first if you want them for comparison
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_contrast_modeling.sbatch

for f in data/bd_model_env_unified_2026-11-12.RData data/01_contrast_jobs.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done

# s. flotator must now have sigma_by_population FALSE (TRUE in the round 3 file)
# (awk, not R: R lives in the container, not on the login node)
flotator_sigma=$(awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) { gsub(/"/, "", $i); col[$i] = i }; next }
  { gsub(/"/, "") } $col["taxon_capture"] == "silverstoneia_flotator" && $col["variant"] == "main" { print $col["sigma_by_population"] }' \
  data/01_contrast_jobs.csv)
if [ "$flotator_sigma" != "FALSE" ]; then
  echo "data/01_contrast_jobs.csv is not the round 4 file (s. flotator still has sigma by population): copy the new one" >&2
  exit 1
fi

# all small (n 176-600): job script defaults (1 h, 4G, adapt_delta 0.99);
# round 3 peak <= 2.4 GB, <= 10 min for jobs this size at 16000 draws
sbatch --array=17,28,30,31 "$JOB"

echo
echo "Submitted round 4 (4 jobs). Check status with: squeue -u \$USER"
