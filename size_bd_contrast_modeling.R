# ============================================================
# Run script: stage 1 contrast models, one per species (and per
# sensitivity variant). Script version of
# size_bd_contrast/02_contrast_modeling.qmd (size_bd_load_analysis repo),
# "Model fitting" section. Model description: size_bd_contrast/README.md.
#
# Each job is a hurdle-lognormal model of one species with no weather or
# climate terms. stage_c (-0.5 adult, +0.5 subadult) is the subadult -
# adult contrast; intercepts vary by population, site-year and 14-day
# event, and the contrast varies at the levels the data support (tier):
#   A: population + site-year + event   B: population + site-year
#   C: population                       none (mu only): no stage term
# Tiers are set separately for prevalence (hu) and load (mu).
#
# Jobs (species, variant, tiers) come from JOBS_FILE
# (01_contrast_jobs.csv, written by size_bd_contrast/01_contrast_tiers.Rmd).
# Copy it to data/ alongside the model data.
#
# Job selection (first match wins):
#   JOB=<taxon_capture>_<variant>   fit just that job, e.g. rana_clamitans_main
#   SLURM_ARRAY_TASK_ID=<i>         fit the job with array_index == i
#                                   (0-based) - see
#                                   hpc/run_size_bd_contrast_modeling.sbatch
#   neither                         fit all jobs sequentially
# ============================================================

# load packages. tidyverse (dplyr, lubridate) comes with the rocker/tidyverse base image.
library(tidyverse)
library(brms)
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS stage contrast models (size_bd_contrast, stage 1):\n\n")

# ------------------------------------------------------------
# Run parameters
# ------------------------------------------------------------
# n_cores below is chains run in parallel; total CPUs actually used is
# n_cores * n_threads, so request that many in --cpus-per-task/-c in Slurm.
n_chains <- as.integer(Sys.getenv("N_CHAINS", "4"))
n_iter <- as.integer(Sys.getenv("N_ITER", "4000"))
n_warmup <- as.integer(Sys.getenv("N_WARMUP", "2000"))
n_cores <- as.integer(Sys.getenv("N_CORES", "4"))
n_threads <- as.integer(Sys.getenv("N_THREADS_PER_CHAIN", "1"))

data_file <- Sys.getenv("DATA_FILE", "data/bd_model_env_unified_2026-11-12.RData")
jobs_file <- Sys.getenv("JOBS_FILE", "data/01_contrast_jobs.csv")
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_contrast")
run_loo <- as.logical(Sys.getenv("RUN_LOO", "TRUE"))

# sampler target acceptance rate. a higher target makes the sampler take smaller steps, which
# follow tight posterior curvature (many small groups: events, site-years; contrast sds near 0)
# without divergences, at ~1.5-2x the run time. it changes how the posterior is explored, not
# the model. default 0.99; r. sierrae is submitted with ADAPT_DELTA=0.95 (large and well
# informed; 0.99 would roughly double an ~11-29 h run). keep in sync with contrast_adapt_delta()
# in 02_contrast_modeling.qmd.
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.99"))

# ============================================================
# BEGIN copy of size_bd_contrast/contrast_data_functions.R (constants and add_contrast_groups
# only; the co-sampling counts and tier rules are not needed here). keep in sync.
# ============================================================
response_col = "bd_load_co100_density"

event_days = 14            # event = population x year x 14-day bin of day of year

# stage_c: -0.5 adult, +0.5 subadult. the intercept is then the midpoint of the two stages and
# the stage_c coefficient the full subadult - adult contrast; a group's stage_c sd is how much
# the contrast varies at that level.
# pop_year / pop_year_event: explicit ids (rather than population:year:... in the formula) so
# parameter names stay readable. 14-day bins are cut at the new year.
add_contrast_groups = function(df) {
  df %>%
    mutate(life_stage_simple = factor(life_stage_simple, levels = c("adult", "subadult")),
           stage_c = if_else(life_stage_simple == "subadult", 0.5, -0.5),
           event_bin = (yday(date) - 1) %/% event_days,
           pop_year = paste(population, year, sep = "_"),
           pop_year_event = paste(population, year, event_bin, sep = "_"))
}
# ============================================================
# END copy of contrast_data_functions.R
# ============================================================

# ============================================================
# BEGIN verbatim copy of size_bd_contrast/contrast_model_functions.R. keep in sync.
# ============================================================
# ---- data ------------------------------------------------------------------------------------
# one species' swabs with stage_c (-0.5 adult / +0.5 subadult), pop_year and pop_year_event
# (14-day events). droplevels so random-effect terms have no empty levels.
prep_contrast_data = function(brms_data, species_name) {
  brms_data %>%
    filter(taxon_capture == species_name) %>%
    add_contrast_groups() %>%
    droplevels()
}

# ---- formula ---------------------------------------------------------------------------------
# random-effect terms by tier. intercepts vary at every level in every tier (they absorb site,
# year and event-wide swings in Bd, e.g. outbreaks, whether or not both stages were swabbed);
# contrast slopes (stage_c) only at the levels the data support:
#   A: population + site-year + event      B: population + site-year
#   C: population                          none (mu only): no stage term
contrast_re_terms = function(tier) {
  switch(tier,
         A = "(1 + stage_c | population) + (1 + stage_c | pop_year) + (1 + stage_c | pop_year_event)",
         B = "(1 + stage_c | population) + (1 + stage_c | pop_year) + (1 | pop_year_event)",
         C = "(1 + stage_c | population) + (1 | pop_year) + (1 | pop_year_event)",
         none = "(1 | population) + (1 | pop_year) + (1 | pop_year_event)",
         stop("unknown tier: ", tier))
}

contrast_fixed_terms = function(tier) if (tier == "none") "1" else "1 + stage_c"

# hurdle lognormal: mu = mean log load among positives (load part), hu = probability of a zero
# (prevalence part, logit scale), sigma = residual sd of log load, varying by population.
# tiers are set separately for mu (tier_mu, load counts) and hu (tier_hu, prevalence counts).
contrast_formula = function(tier_hu, tier_mu, response = response_col) {
  bf(as.formula(paste(response, "~", contrast_fixed_terms(tier_mu), "+", contrast_re_terms(tier_mu))),
     as.formula(paste("hu ~", contrast_fixed_terms(tier_hu), "+", contrast_re_terms(tier_hu))),
     sigma ~ 1 + (1 | population))
}

# ---- priors ----------------------------------------------------------------------------------
# groups whose contrast (stage_c) varies, by tier (matches contrast_re_terms)
contrast_slope_groups = function(tier) {
  switch(tier,
         A = c("population", "pop_year", "pop_year_event"),
         B = c("population", "pop_year"),
         C = "population",
         none = character(0),
         stop("unknown tier: ", tier))
}

# contrast sd priors, one per group with a stage_c slope (brms needs the group named for a
# coefficient-level sd prior). dpar = "" is the mu part.
contrast_sd_priors = function(tier, dpar = "") {
  groups = contrast_slope_groups(tier)
  if (length(groups) == 0) return(NULL)
  do.call(c, lapply(groups, function(g) {
    set_prior("exponential(2)", class = "sd", coef = "stage_c", group = g, dpar = dpar)
  }))
}

# fixed effects and intercept sds as in the unified / seasonal models. the contrast sds get
# exponential(2): tighter than the intercept sds (exponential(1)), so a consistent contrast is the
# default and variation in it must be earned from the data. prior median 0.35, 95% below 1.5
# (log / logit units). brms errors on priors that match no parameter, so stage terms only get
# priors when the part has them (tier != "none").
contrast_priors = function(tier_hu, tier_mu) {
  p = c(
    # mu: log load density among positives
    prior(normal(5, 4), class = "Intercept"),
    prior(exponential(1), class = "sd"),
    # hu: logit of non-detection
    prior(normal(0, 2), class = "Intercept", dpar = "hu"),
    prior(normal(0, 2), class = "b", dpar = "hu"),
    prior(exponential(1), class = "sd", dpar = "hu"),
    contrast_sd_priors(tier_hu, dpar = "hu"),
    # sigma: log residual sd of log load
    prior(normal(0, 1), class = "Intercept", dpar = "sigma"),
    prior(exponential(1), class = "sd", dpar = "sigma"),
    # intercept-contrast correlations (every tier of hu has at least one)
    prior(lkj(2), class = "cor")
  )
  if (tier_mu != "none") {
    p = c(p,
          prior(normal(0, 3), class = "b"),
          contrast_sd_priors(tier_mu))
  }
  p
}

# contrast sd prior quantiles, for the prior vs posterior check
contrast_sd_prior_quantile = function(p) qexp(p, rate = 2)
# ============================================================
# END copy of contrast_model_functions.R
# ============================================================

# ------------------------------------------------------------
# Load data, choose jobs
# ------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat(paste0("Loading job table from: ", jobs_file, "\n"))
contrast_jobs <- read.csv(jobs_file, stringsAsFactors = FALSE) %>%
  arrange(array_index) %>%
  mutate(job = paste(taxon_capture, variant, sep = "_"))

cat(paste0("Loading data from: ", data_file, "\n\n"))
load(data_file)

missing_species <- setdiff(contrast_jobs$taxon_capture, unique(as.character(brms_data$taxon_capture)))
if (length(missing_species) > 0) {
  stop(paste0("species not in data: ", paste(missing_species, collapse = ", ")))
}

job_env <- Sys.getenv("JOB", "")
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (job_env != "") {
  if (!job_env %in% contrast_jobs$job) {
    stop(paste0("JOB '", job_env, "' not in jobs file. Available: ", paste(contrast_jobs$job, collapse = ", ")))
  }
  run_jobs <- filter(contrast_jobs, job == job_env)
} else if (array_task_id != "") {
  i <- as.integer(array_task_id)
  if (!i %in% contrast_jobs$array_index) {
    stop(paste0("SLURM_ARRAY_TASK_ID ", array_task_id, " out of range: array_index runs 0-", max(contrast_jobs$array_index)))
  }
  run_jobs <- filter(contrast_jobs, array_index == i)
} else {
  run_jobs <- contrast_jobs
}

cat("All jobs (array index: job, tier_hu x tier_mu):\n")
cat(paste0("  ", contrast_jobs$array_index, ": ", contrast_jobs$job, " (",
           contrast_jobs$tier_hu, " x ", contrast_jobs$tier_mu, ")", collapse = "\n"), "\n\n")
cat(paste0("Fitting: ", paste(run_jobs$job, collapse = ", "), "\n\n"))

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------
model_contrast_job <- function(job_row) {
  cat(paste0("==== Modeling stage contrast for: ", job_row$job,
             " (tier hu ", job_row$tier_hu, ", tier mu ", job_row$tier_mu, ") ====\n\n"))

  model_data <- prep_contrast_data(brms_data, job_row$taxon_capture)

  # group sizes: how many levels each random effect has (the contrast sds at a level are only as
  # well informed as its number of levels, mainly the co-sampled ones)
  cat(paste0("n = ", nrow(model_data),
             " (adult ", sum(model_data$life_stage_simple == "adult"),
             ", subadult ", sum(model_data$life_stage_simple == "subadult"), ")",
             "; populations = ", n_distinct(model_data$population),
             ", site-years = ", n_distinct(model_data$pop_year),
             ", events = ", n_distinct(model_data$pop_year_event),
             "; adapt_delta = ", adapt_delta, "\n\n"))

  model_formula <- contrast_formula(job_row$tier_hu, job_row$tier_mu)
  model_priors <- contrast_priors(job_row$tier_hu, job_row$tier_mu)
  print(model_formula)
  print(model_priors)

  # refresh: print progress regularly, flushed live to the log (see stdbuf
  # in hpc/*.sbatch) so a running job's pace is visible via `tail -f`
  fit_start_time <- Sys.time()
  fit <- brm(
    formula  = model_formula,
    data     = model_data,
    family   = hurdle_lognormal(),
    prior    = model_priors,
    control  = list(adapt_delta = adapt_delta, max_treedepth = 12),
    chains   = n_chains,
    iter     = n_iter,
    warmup   = n_warmup,
    seed     = 128,
    cores    = n_cores,
    threads  = threading(n_threads),
    backend  = "cmdstanr",
    refresh  = max(1, n_iter %/% 50)
  )
  cat(paste0(
    "\nFit time (", n_chains, " chains x ", n_threads, " threads/chain, ",
    n_iter, " iter, ", n_warmup, " warmup): ",
    format(Sys.time() - fit_start_time), "\n\n"
  ))

  # save raw fit FIRST so it's never lost to a loo failure. the with-loo fit is saved under the
  # name 02_contrast_modeling.qmd expects (size_bd_contrast_model_zln_c1_<species>_<variant>.rds).
  model_file <- file.path(out_dir, paste0("size_bd_contrast_model_zln_c1_", job_row$job, ".rds"))
  raw_file <- sub("\\.rds$", "_no_loo.rds", model_file)
  saveRDS(fit, file = raw_file)
  cat(paste0("Model fit (no loo) saved to:\n\n\t", raw_file, "\n\n"))

  if (run_loo) {
    tryCatch({
      cat("Adding LOO...\n")
      fit <- add_criterion(fit, "loo")
      saveRDS(fit, file = model_file)
      cat(paste0("Model fit with loo saved to:\n\n\t", model_file, "\n\n"))
    }, error = function(e) {
      cat("(could not compute/save loo: ", conditionMessage(e), ")\n\n")
    })
  }

  cat(paste0("Done: ", job_row$job, "\n\n"))
}

for (k in seq_len(nrow(run_jobs))) {
  model_contrast_job(run_jobs[k, ])
}
