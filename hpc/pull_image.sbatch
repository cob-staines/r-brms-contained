#!/bin/bash
# ============================================================
# One-off job: pulls/converts the published Docker image into a
# local r-brms.sif. Run this once before submitting the other
# hpc/*.sbatch jobs, if you're not confident they'll get enough
# memory to build the SIF themselves (squashfs conversion of a
# multi-GB image is RAM-hungry, separate from the memory the
# actual R model run needs). Once r-brms.sif exists, the other
# scripts skip pulling and just use the cached file.
#
# Submit with: sbatch hpc/pull_image.sbatch
# ============================================================
#SBATCH -p grit_nodes
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH -J pull-r-brms-image
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module load singularity 2>/dev/null || true

singularity pull "$PWD/r-brms.sif" docker://ghcr.io/cob-staines/r-brms-contained:latest
