#!/bin/bash
# ============================================================
# Per-species seasonal stage x Bd models: fits size_bd_season_modeling.R
# via brms/cmdstanr as a Slurm array job, one species per task, so all
# species run in parallel instead of back to back.
#
# Submit with: sbatch hpc/run_size_bd_season_modeling.sbatch
#
# Array index -> species (alphabetical order of season_species in
# size_bd_season_modeling.R, printed at the top of each task's log):
#   0 anaxyrus_canorus    (n ~1.0k)
#   1 rana_cascadae       (n ~3.4k)
#   2 rana_catesbeiana    (n ~1.8k)
#   3 rana_clamitans      (n ~2.6k)
#   4 rana_muscosa        (n ~6.8k)
#   5 rana_sierrae        (n ~26.6k)
# If season_species changes, update --array below to match.
#
# Rerun a subset of species by overriding --array on the command line:
#   sbatch --array=1,5 hpc/run_size_bd_season_modeling.sbatch
#
# Resource tiers (sbatch options on the command line override the #SBATCH
# defaults below; -c sets threads per chain, see N_THREADS_PER_CHAIN).
# Times come from the local timing test in size_bd_season/02_season_modeling.qmd
# (4 chains, 1 thread each, 300 warmup + 100 sampling, scaled to 2000 + 2000):
#   anaxyrus_canorus  n  1.0k   est ~0.1 h
#   rana_cascadae     n  3.4k   est ~0.2 h
#   rana_sierrae      n 26.6k   est ~6.5 h
# Run time grows faster than n (~n^1.6 between cascadae and sierrae), so the
# untested species are interpolated: catesbeiana / clamitans ~0.2 h, muscosa
# ~0.7 h. Requests allow ~4-8x headroom for a slower cluster node, the
# warmup extrapolation (the least reliable part), per-species posterior
# geometry, and the LOO step after sampling. The default (-t 1 day) covers
# every species, so a plain `sbatch hpc/run_size_bd_season_modeling.sbatch`
# runs all 6. Shorter requests for the small species start sooner when the
# queue is busy (backfill), so optionally split the submission:
#   sbatch --array=0-4 --mem=16G -t 06:00:00   hpc/run_size_bd_season_modeling.sbatch  # n 1-6.8k, est <= ~0.7 h
#   sbatch --array=5   --mem=48G -t 1-00:00:00 hpc/run_size_bd_season_modeling.sbatch  # n 26.6k, est ~6.5 h
# Memory was not measured locally (still a guess; 48G matches the GAM tier for
# r. sierrae). After the first run, check actual use and trim:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
# If a task times out, the raw fit is saved before LOO, so check models/
# before resubmitting.

# Required input (gitignored - rsync it up first, see README):
#   DATA_FILE   data/bd_model_env_unified_2026-11-12.RData (the m11g data file)
# Override with --export=ALL,DATA_FILE=...
# Output goes to models/size_bd_season/ (OUT_DIR).
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0-5
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=48G               # guess sized for r. sierrae; see tiers above
#SBATCH -t 1-00:00:00           # covers r. sierrae (est ~6.5 h); watch progress with `tail -f`
#SBATCH -J size_bd_season
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models/size_bd_season

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
  stdbuf -oL -eL Rscript size_bd_season_modeling.R
