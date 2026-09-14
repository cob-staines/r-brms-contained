#!/bin/bash
# ============================================================
# Quick sanity-check job: fits brms_cars_test.R (a ~50-row toy
# model) to confirm the container + Slurm setup work end to end
# BEFORE spending core-hours on the real long-running model.
#
# Submit with: sbatch hpc/test_cars.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -J r-brms-test
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail
mkdir -p logs models

module load singularity 2>/dev/null || true

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "Pulling container image..."
  singularity pull "$IMAGE" docker://ghcr.io/cob-staines/r-brms-contained:latest
fi

singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" Rscript brms_cars_test.R
