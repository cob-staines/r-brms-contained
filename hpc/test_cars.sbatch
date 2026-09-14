#!/bin/bash
# ============================================================
# Quick sanity-check job: fits brms_cars_test.R (a ~50-row toy
# model) to confirm the container + Slurm setup work end to end
# BEFORE spending core-hours on the real long-running model.
#
# Submit with: sbatch hpc/test_cars.sbatch
# ============================================================
#SBATCH --job-name=r-brms-test
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=REPLACE_ME       # TODO: run `sinfo` and pick a partition
#SBATCH --account=REPLACE_ME         # TODO: run `sacctmgr show associations user=$USER` for your account/allocation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail
mkdir -p logs models

module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "Pulling container image..."
  apptainer pull "$IMAGE" docker://ghcr.io/cob-staines/r-brms-contained:latest
fi

apptainer exec --bind "$PWD:/project" --pwd /project "$IMAGE" Rscript brms_cars_test.R
