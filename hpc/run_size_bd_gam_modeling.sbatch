#!/bin/bash
# ============================================================
# Per-species size ~ Bd GAMs: fits size_bd_gam_modeling.R via
# brms/cmdstanr as a Slurm array job, one species per task, so
# all species run in parallel instead of back to back.
#
# Submit with: sbatch hpc/run_size_bd_gam_modeling.sbatch
#
# Array index -> species is the alphabetical order of taxon_capture
# in DATA_FILE (printed at the top of each task's log):
#   0 colostethus_panamansis      7 rana_catesbeiana
#   1 lithobates_chiricahuensis   8 rana_clamitans
#   2 lithobates_sphenocephalus   9 rana_muscosa
#   3 lithobates_warszewitschii  10 rana_pipiens
#   4 notophthalmus_viridescens  11 rana_sierrae
#   5 pseudacris_crucifer        12 silverstoneia_flotator
#   6 rana_cascadae
# If the species in the data change, update --array below to match.
#
# Rerun a subset of species by overriding --array on the command line:
#   sbatch --array=6,11 hpc/run_size_bd_gam_modeling.sbatch
#
# Resource tiers (sbatch options on the command line override the #SBATCH
# defaults below; -c sets threads per chain, see N_THREADS_PER_CHAIN):
#   sbatch --array=0,1,2,3,10,12 -c 4 --mem=16G -t 00:30:00  hpc/run_size_bd_gam_modeling.sbatch  # n < 1k
#   sbatch --array=5,7,8         -c 4 --mem=16G -t 03:00:00  hpc/run_size_bd_gam_modeling.sbatch  # n 1-3k
#   sbatch --array=4,6,9         -c 8 --mem=24G -t 3-00:00:00 hpc/run_size_bd_gam_modeling.sbatch # n 5-7k
#   sbatch --array=11            -c 8 --mem=48G -t 3-00:00:00 hpc/run_size_bd_gam_modeling.sbatch # n 27k
#
# Required inputs (both gitignored - rsync them up first, see README):
#   DATA_FILE         data/brms_data_gam_2026-10-28.RData
#   PRIOR_MODEL_FILE  default: newest models/unified_model_results_zln_m11d_co100_*_with_loo.rds
#                     (written by hpc/run_size_bd_cutoffs.sbatch, task 2)
# Override either with --export=ALL,DATA_FILE=...,PRIOR_MODEL_FILE=...
#
# To start automatically once the co100 fit finishes:
#   sbatch --dependency=afterok:<cutoff jobid>_2 hpc/run_size_bd_gam_modeling.sbatch
# Output goes to models/size_bd_gam/ (OUT_DIR).
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=0-12
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=24G                # TODO: guess. Each task loads the ~850MB (compressed) unified
                                  #       model to build the env prior, then frees it before fitting.
                                  #       Check `sacct -j <jobid> --format=JobID,MaxRSS` after a run
                                  #       (rana_sierrae, ~27k rows, is the largest) and adjust.
#SBATCH -t 3-00:00:00             # TODO: guess; watch live progress with `tail -f` rather than
                                  #       trusting this. Cluster max is 30 days.
#SBATCH -J size_bd_gam
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models/size_bd_gam

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
  stdbuf -oL -eL Rscript size_bd_gam_modeling.R
