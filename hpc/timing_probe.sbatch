#!/bin/bash
# ============================================================
# Timing/memory probe: runs the REAL model (same data, formula, family,
# adapt_delta) but with a small iter/warmup, to gauge how long and how much
# memory the full run in run_size_bd_modeling.sbatch will need.
#
# Runtime caveat: Stan's warmup adaptation is NOT linear in warmup length -
# its default windowed adaptation schedule barely engages for a very short
# warmup (e.g. 10-50 iterations), so tiny probes can finish suspiciously
# fast without doing the real adaptation work, then a slightly bigger one
# can suddenly hit expensive tree-depth blowups as adaptation properly
# kicks in. Don't trust linear extrapolation from a short probe here -
# use the live progress in this job's log (refresh is set in
# size_bd_modeling.R, flushed immediately via stdbuf below) to watch
# actual iteration pace with `tail -f logs/<this job>.out`, and check the
# per-chain warmup/sampling breakdown this script prints at the end.
#
# Memory, unlike runtime, scales close to linearly with chains and kept
# draws, so this probe's memory usage IS a reasonable basis for
# extrapolation. N_CHAINS is set to match production (4) rather than 1 for
# this reason: each chain is a separate OS process (memory isn't shared
# between them), so testing at chains=1 would underestimate real memory
# need by up to 4x. Since chains run in parallel, this costs no extra
# wall-clock time, only more CPU/memory requested. After this job finishes,
# check `sacct -j <jobid> --format=MaxRSS` for the real number, then also
# see how it grows if you bump N_WARMUP/N_ITER here before extrapolating
# to the full 2000/4000.
#
# Submit with: sbatch hpc/timing_probe.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 4                    # must equal N_CORES * N_THREADS_PER_CHAIN below
#SBATCH --mem=64G                # generous headroom for measuring, not a final estimate - a 16GB
                                  # local VM already crashed at full production settings
#SBATCH -t 04:00:00              # TODO: unverified guess - bump further if you still hit timeouts
#SBATCH -J size_bd_timing_probe
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

# 4 chains (matching production) but a small warmup/iter - swap N_WARMUP/
# N_ITER for the full 2000/4000 production run once you have a number.
# Keep N_THREADS_PER_CHAIN matched to whatever run_size_bd_modeling.sbatch
# uses, since threading changes per-iteration wall time and memory.
export N_CHAINS=4
export N_CORES=4
export N_THREADS_PER_CHAIN=1
export N_WARMUP=50
export N_ITER=100

# stdbuf forces line-buffered stdout/stderr so progress flushes to the log
# file live instead of sitting in a block buffer until the process exits -
# lets you `tail -f` a still-running job instead of waiting blind.
singularity exec --bind "$PWD:/project" --pwd /project "$IMAGE" \
  stdbuf -oL -eL Rscript size_bd_modeling.R
