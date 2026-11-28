#!/bin/bash
# ============================================================
# GAMs of Bd by body size (size_bd_contrast, GAM option 1): submits all 14
# species as four array jobs, one per resource group, so every species runs
# in parallel with a request sized to it. Groups and why: header of
# hpc/run_size_bd_contrast_gam_modeling.sbatch (time from the timing
# test in 06_gam_modeling.qmd; memory provisional until a first sacct).
#
# Run from the repo root:
#   bash hpc/submit_size_bd_contrast_gam.sh
#
# Before submitting:
#   - git pull (this repo), so size_bd_contrast_gam_modeling.R carries the
#     current copy of gam_model_functions.R
#   - copy data/06_gam_model_data.rds and data/06_gam_species.csv from
#     size_bd_contrast/intermediates/ (the same run of 06_gam_modeling.qmd)
#   - fits overwrite models/size_bd_contrast_gam/gam_m4_<species>.rds
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_contrast_gam_modeling.sbatch

for f in data/06_gam_model_data.rds data/06_gam_species.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done

# the array indices below assume these 14 species (alphabetical); stop if the table changed
# (awk, not R: R lives in the container, not on the login node)
n_species=$(awk 'END { print NR - 1 }' data/06_gam_species.csv)
if [ "$n_species" != "14" ]; then
  echo "data/06_gam_species.csv has $n_species species, not 14: update the array indices" >&2
  exit 1
fi

# small (n < 1k): job script defaults (1 h, 4G)
sbatch --array=1,2,3,4,11,13 "$JOB"
# intermediate (n 1-3k)
sbatch --array=0,6,8,9 --mem=8G -t 02:00:00 "$JOB"
# large (n 5-7k)
sbatch --array=5,7,10 --mem=12G -t 05:00:00 "$JOB"
# r. sierrae (n 26.8k)
sbatch --array=12 --mem=40G -t 1-00:00:00 "$JOB"

echo
echo "Submitted 4 array jobs (14 species). Check status with: squeue -u \$USER"
