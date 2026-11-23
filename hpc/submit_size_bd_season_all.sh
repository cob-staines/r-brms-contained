#!/bin/bash
# ============================================================
# Submit the all-species seasonal models in their three resource groups.
# Each group is a separate submission of the same job script
# (run_size_bd_season_all_modeling.sbatch): command-line options override
# that script's #SBATCH defaults, so each group requests only the memory
# and time its species need (see the script's header for the estimates).
#
# Run from the repo root (the job script expects data/, logs/, models/
# relative to it):
#   bash hpc/submit_size_bd_season_all.sh
#
# Array indices = array_index in data/01_season_species_all.csv (34 species,
# 2026-11 screen). If the species file changes, update the groups below.
#
# Before submitting: data/bd_model_env_unified_2026-11-12.RData and
# data/01_season_species_all.csv must exist, and r-brms.sif should already
# be pulled (hpc/pull_image.sbatch), so the array tasks don't race to pull it.
# ============================================================
set -euo pipefail

JOB=hpc/run_size_bd_season_all_modeling.sbatch

# fail early if inputs are missing (otherwise every array task fails on its own)
for f in data/bd_model_env_unified_2026-11-12.RData data/01_season_species_all.csv r-brms.sif; do
  if [ ! -f "$f" ]; then
    echo "missing: $f" >&2
    exit 1
  fi
done

# small: 28 species, n < 1.5k. uses the job script's defaults
# (4 h, 8G, ADAPT_DELTA 0.99 from the R script)
sbatch "$JOB"

# intermediate: n. viridescens (16), r. cascadae (20), r. catesbeiana (21),
# r. clamitans (22), r. muscosa (24); n 1.8-6.8k
sbatch --array=16,20-22,24 --mem=16G -t 12:00:00 "$JOB"

# large: r. sierrae (26), n 26.6k, at adapt_delta 0.95 (see the job script)
sbatch --array=26 --mem=48G -t 2-00:00:00 --export=ALL,ADAPT_DELTA=0.95 "$JOB"

echo
echo "Submitted. Check status with: squeue -u \$USER"
