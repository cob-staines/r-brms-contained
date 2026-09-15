#!/bin/bash
# ============================================================
# Long-running job: fits size_bd_modeling.R via brms/cmdstanr.
# Run hpc/test_cars.sbatch successfully first to validate the
# container + Slurm setup before submitting this.
#
# Submit with: sbatch hpc/run_size_bd_modeling.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 4                    # must equal N_CORES * N_THREADS_PER_CHAIN below (see brms threading:
                                 # https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html)
#SBATCH --mem=32G                # TODO: unverified guess for a model with this many random-effect groups.
                                  #       After the first run, check actual usage with `seff <jobid>`
                                  #       or `sacct -j <jobid> --format=MaxRSS` and tune this.
#SBATCH -t 5-00:00:00             # TODO: unverified guess. Timing probes on this model showed highly
                                  #       non-linear warmup cost (Stan's windowed adaptation barely
                                  #       engages on a short warmup, so small probes don't extrapolate
                                  #       reliably) - set generously and watch live progress instead
                                  #       (`tail -f logs/<this job>.out`) rather than trusting a probe
                                  #       estimate alone. Cluster max is 30 days.
#SBATCH -J size_bd_model
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

# N_CORES * N_THREADS_PER_CHAIN must equal -c above. Threading (>1) splits
# each chain's likelihood evaluation across threads via reduce_sum - useful
# if you have spare CPUs and want per-chain speedup rather than just more
# parallel chains.
export N_CORES=4
export N_THREADS_PER_CHAIN=1

# stdbuf forces line-buffered stdout/stderr so progress flushes to the log
# file live instead of sitting in a block buffer until the process exits -
# lets you `tail -f` this job's log to see real iteration pace as it runs.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_modeling.R
