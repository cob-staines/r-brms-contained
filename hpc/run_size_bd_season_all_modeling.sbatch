#!/bin/bash
# ============================================================
# Per-species seasonal stage x Bd models, all species: fits
# size_bd_season_all_modeling.R via brms/cmdstanr as a Slurm array job,
# one species per task.
#
# Array index -> species: the array_index column of
# data/01_season_species_all.csv (alphabetical, 0-based; printed at the top
# of each task's log). 34 species (0-33) as of the 2026-11 screen. If the
# species file changes, update the array indices in the groups below.
#
# Resource groups (sbatch options on the command line override the #SBATCH
# defaults below, which are the small group). Estimates from the local timing
# test in size_bd_season/02_season_all_species.qmd (4 chains, 300 warmup +
# 100 sampling, scaled to 2000 + 2000; slowest chain):
#   s. flotator      n 0.6k   est 0.07 h   (6 divergences in the test)
#   r. catesbeiana   n 1.8k   est 0.44 h
#   n. viridescens   n 5.3k   est 1.5 h    (10 divergences in the test)
# Time grows faster than n; other species are extrapolated from these.
#
#   small (n < 1.5k, est <= ~0.3 h): 28 species, 4 h / 8G (the defaults)
#     sbatch hpc/run_size_bd_season_all_modeling.sbatch
#   intermediate (n 1.8-6.8k, est ~0.4-3 h): n. viridescens (16), r. cascadae (20),
#     r. catesbeiana (21), r. clamitans (22), r. muscosa (24)
#     sbatch --array=16,20-22,24 --mem=16G -t 12:00:00 hpc/run_size_bd_season_all_modeling.sbatch
#   large (n 26.6k, est ~13-15 h): r. sierrae (26)
#     sbatch --array=26 --mem=48G -t 2-00:00:00 hpc/run_size_bd_season_all_modeling.sbatch
#
# Headroom (~4x or more) covers slower nodes, the warmup extrapolation (the
# least reliable part), and a resubmission at a higher ADAPT_DELTA (roughly
# 1.5-2x slower).
#
# Divergences: fit first at the default ADAPT_DELTA (0.95), check the
# convergence table in 02_season_all_species.qmd, then resubmit only the
# species with divergences at a higher target, keeping each species' group
# resources, e.g.:
#   sbatch --array=32 --export=ALL,ADAPT_DELTA=0.99 hpc/run_size_bd_season_all_modeling.sbatch
#   sbatch --array=16 --mem=16G -t 12:00:00 --export=ALL,ADAPT_DELTA=0.99 hpc/run_size_bd_season_all_modeling.sbatch
# A resubmission overwrites that species' model files in OUT_DIR.
#
# Memory: not measured. Check actual use and trim with
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
#
# Required inputs (gitignored - rsync them up first, see README):
#   DATA_FILE     data/bd_model_env_unified_2026-11-12.RData (the m11g data file)
#   SPECIES_FILE  data/01_season_species_all.csv (from size_bd_season/intermediates/)
# Override with --export=ALL,DATA_FILE=...,SPECIES_FILE=... (also ADAPT_DELTA)
# Output goes to models/size_bd_season_all/ (OUT_DIR). If a task times out
# during LOO, the raw fit is already saved as *_no_loo.rds.
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0-15,17-19,23,25,27-33   # small group; see above
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=8G                # small group; see above
#SBATCH -t 04:00:00             # small group (est <= ~0.3 h); see above
#SBATCH -J size_bd_season_all
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models/size_bd_season_all

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
# SLURM_ARRAY_TASK_ID is passed through to pick this task's species.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_season_all_modeling.R
