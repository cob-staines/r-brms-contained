# ============================================================
# Dockerfile: R + brms + cmdstanr + tidyverse
# Base: rocker/tidyverse (Ubuntu 24.04, R + tidyverse pre-installed)
# ============================================================

FROM rocker/tidyverse:4.4.2

# ---- System dependencies --------------------------------
# libicu-dev : required by stringi
# pandoc     : required by bayesplot, loo, rstantools
# libglpk-dev: required by igraph (brms dependency)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libicu-dev \
    pandoc \
    libglpk-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---- R packages -----------------------------------------
# brms depends on rstan, which must be installed as a pre-compiled binary.
# Including stan-dev.r-universe.dev in repos causes R to pick a newer
# dev version of rstan/StanHeaders that has no binary, falls back to source,
# and fails with an Eigen 3.3/3.4 API mismatch. So we install brms and its
# Stan dependencies from PPM only (rocker's default, pre-compiled binaries),
# then install cmdstanr separately from stan-dev (it is not on CRAN/PPM).
RUN Rscript -e " \
    options(Ncpus = parallel::detectCores()); \
    install.packages(c('brms', 'posterior', 'bayesplot', 'loo'), \
                     repos = getOption('repos')); \
    install.packages('cmdstanr', \
                     repos = c('https://stan-dev.r-universe.dev', getOption('repos'))); \
    pkgs <- c('brms', 'cmdstanr', 'posterior', 'bayesplot', 'loo'); \
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; \
    if (length(missing)) stop('Failed to install: ', paste(missing, collapse = ', ')) \
"

# ---- Install CmdStan ------------------------------------
# Install to /opt rather than the default $HOME-relative path: this build
# runs as root (HOME=/root, mode 700), but containers run on HPC via
# Singularity execute as the invoking host user, not root, who can neither
# find CmdStan via cmdstanr's default $HOME-based lookup nor read into
# /root at all. /opt/cmdstan/current is a stable, world-readable path, and
# CMDSTAN points cmdstanr at it regardless of the runtime user's $HOME.
RUN mkdir -p /opt/cmdstan && \
    Rscript -e " \
    cmdstanr::install_cmdstan(dir = '/opt/cmdstan', cores = parallel::detectCores()) \
" && \
    ln -s "$(ls -d /opt/cmdstan/cmdstan-*)" /opt/cmdstan/current && \
    chmod -R a+rX /opt/cmdstan
ENV CMDSTAN=/opt/cmdstan/current

# ---- Set working directory ------------------------------
WORKDIR /project

# ---- Default command ------------------------------------
# Runs an interactive R session; override with your script at runtime
CMD ["R"]
