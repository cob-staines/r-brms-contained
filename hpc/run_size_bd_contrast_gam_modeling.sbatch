#!/bin/bash
# ============================================================
# GAMs of Bd by body size (size_bd_contrast, GAM option 1): fits
# size_bd_contrast_gam_modeling.R via brms/cmdstanr as a Slurm array
# job, one species per task, so all species run in parallel.
# Replaces hpc/run_size_bd_gam_modeling.sbatch (m3), kept until retired.
#
# Submit all four resource groups with:
#   bash hpc/submit_size_bd_contrast_gam.sh
#
# Array index -> species: alphabetical order of taxon_capture in
# SPECIES_FILE (printed at the top of each task's log). GAM n in brackets:
#   0 anaxyrus_canorus          (1969)    7 rana_cascadae       (5447)
#   1 colostethus_panamansis     (783)    8 rana_catesbeiana    (1876)
#   2 lithobates_chiricahuensis  (459)    9 rana_clamitans      (2701)
#   3 lithobates_sphenocephalus  (970)   10 rana_muscosa        (6843)
#   4 lithobates_warszewitschii  (565)   11 rana_pipiens         (584)
#   5 notophthalmus_viridescens (5255)   12 rana_sierrae       (26784)
#   6 pseudacris_crucifer       (1244)   13 silverstoneia_flotator (602)
# If the species table changes, update the indices here and in the
# submit script.
#
# Resource groups (sbatch options on the command line override the #SBATCH
# defaults below, which are the small group).
# Time: from the timing test in 06_gam_modeling.qmd ("Timing test", 2026-11-27):
# short local runs, extrapolated at 8000 draws with run time ~ n^1.2 (likely)
# to n^1.6 (upper bound). The fitted exponent (~2.2) is not used (driven by the
# smallest species; see the chunk). Measured: r. catesbeiana 0.31 h, r. clamitans
# 0.38 h, both faster than their stage 1 models.
# Memory: not yet measured. From the stage 1 fits of the same species (sacct,
# rounds 1-3; header of hpc/run_size_bd_contrast_modeling.sbatch), which had the
# same nuisance intercepts at twice the draws (16000 vs 8000), so generous.
# Replace with sacct after the first run.
#   small (n < 1k): 1, 2, 3, 4, 11, 13    est <= 0.07 h       -> 1 h / 4G (the defaults)
#   intermediate (n 1-3k): 0, 6, 8, 9     est <= 0.38 h       -> 2 h / 8G
#   large (n 5-7k): 5, 7, 10              est 0.85-1.7 h      -> 5 h / 12G
#   sierrae (n 26.8k): 12                 est 6-15 h          -> 24 h / 40G
#     (stage 1: 5.9 h for a more complex model at 16000 draws; a timeout loses
#     the whole fit, so the limit sits well above the 1.6 bound)
#
# Rerun a subset by overriding --array, e.g.
#   sbatch --array=10 --mem=12G -t 05:00:00 hpc/run_size_bd_contrast_gam_modeling.sbatch
# Divergences: printed at the end of each log; check the fit health table in
# 06_gam_modeling.qmd. To resubmit at a higher target:
#   --export=ALL,ADAPT_DELTA=0.99 (about 1.5-2x the run time)
# A resubmission overwrites that species' fit in OUT_DIR.
#
# Memory / time of a finished run:
#   sacct -j <jobid> --format=JobID%22,Elapsed,MaxRSS,ReqMem | grep -E "batch|JobID"
#
# Required inputs (gitignored - rsync them up first, see README), both written
# by size_bd_contrast/06_gam_modeling.qmd into size_bd_contrast/intermediates/:
#   DATA_FILE     data/06_gam_model_data.rds
#   SPECIES_FILE  data/06_gam_species.csv
# Override with --export=ALL,DATA_FILE=...,SPECIES_FILE=...
# Output: models/size_bd_contrast_gam/gam_m4_<species>.rds (OUT_DIR); copy into
# size_bd_contrast/models/ in size_bd_load_analysis.
#
# Run hpc/pull_image.sbatch once first if r-brms.sif doesn't exist yet -
# concurrent array tasks all trying to pull it the first time will race.
# ============================================================
#SBATCH -p grit_nodes
#SBATCH --array=1,2,3,4,11,13   # small group; see above
#SBATCH -c 4                    # multiple of N_CORES (4); -c 8 gives 2 threads per chain
#SBATCH --mem=4G                # small group; provisional, see above
#SBATCH -t 01:00:00             # small group; provisional, see above
#SBATCH -J size_bd_contrast_gam
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs models/size_bd_contrast_gam

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
  stdbuf -oL -eL Rscript size_bd_contrast_gam_modeling.R
