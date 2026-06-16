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
# cmdstanr requires a separate step to install the CmdStan binaries
RUN Rscript -e " \
    cmdstanr::install_cmdstan(cores = parallel::detectCores()) \
"

# ---- Set working directory ------------------------------
WORKDIR /project

# ---- Default command ------------------------------------
# Runs an interactive R session; override with your script at runtime
CMD ["R"]
