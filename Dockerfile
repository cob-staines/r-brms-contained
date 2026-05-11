# ============================================================
# Dockerfile: R + brms + cmdstanr + tidyverse
# Base: rocker/tidyverse (Debian-based, R + tidyverse pre-installed)
# ============================================================

FROM rocker/tidyverse:4.4.2

# ---- System dependencies --------------------------------
# libicu-dev : required by stringi
# pandoc     : required by bayesplot, loo, rstantools
# libglpk-dev: required by igraph (brms dependency)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libicu-dev \
    pandoc \
    pandoc-citeproc \
    libglpk-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---- R packages -----------------------------------------
# Install from source to avoid binary compatibility issues on Linux
RUN Rscript -e " \
    options(pkgType = 'source', Ncpus = parallel::detectCores()); \
    install.packages(c( \
        'brms', \
        'cmdstanr', \
        'posterior', \
        'bayesplot', \
        'loo' \
    ), repos = c( \
        cmdstanr = 'https://stan-dev.r-universe.dev', \
        CRAN     = 'https://cloud.r-project.org' \
    )) \
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
