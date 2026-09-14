#!/bin/bash
# ============================================================
# Long-running job: fits size_bd_modeling.R via brms/cmdstanr.
# Run hpc/test_cars.sbatch successfully first to validate the
# container + Slurm setup before submitting this.
#
# Submit with: sbatch hpc/run_size_bd_modeling.sbatch
# ============================================================
#SBATCH --job-name=size_bd_model
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=REPLACE_ME       # TODO: run `sinfo` and pick a partition
#SBATCH --account=REPLACE_ME         # TODO: run `sacctmgr show associations user=$USER` for your account/allocation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4            # matches cores=4/chains=4 in size_bd_modeling.R (one core per chain)
#SBATCH --mem=32G                    # TODO: unverified guess for a model with this many random-effect groups.
                                      #       After the first run, check actual usage with `seff <jobid>`
                                      #       or `sacct -j <jobid> --format=MaxRSS` and tune this.
#SBATCH --time=24:00:00              # TODO: unverified guess. Time the first run and tighten this so
                                      #       your job queues faster (shorter requests often schedule sooner).

set -euo pipefail
mkdir -p logs models

module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true

IMAGE="$PWD/r-brms.sif"
if [ ! -f "$IMAGE" ]; then
  echo "Pulling container image..."
  apptainer pull "$IMAGE" docker://ghcr.io/cob-staines/r-brms-contained:latest
fi

apptainer exec --bind "$PWD:/project" --pwd /project "$IMAGE" Rscript size_bd_modeling.R
