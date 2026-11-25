# ============================================================
# Run script: stage 1 contrast models, one per species (and per
# sensitivity variant). Script version of
# size_bd_contrast/02_contrast_modeling.qmd (size_bd_load_analysis repo),
# "Model fitting" section. Model description: size_bd_contrast/README.md.
#
# Each job is a hurdle-lognormal model of one species with no weather or
# climate terms. stage_c (-0.5 adult, +0.5 subadult) is the subadult -
# adult contrast; intercepts vary by population, site-year and survey
# event, and the contrast varies at the levels the data support (tier):
#   A: population + site-year + event   B: population + site-year
#   C: population                       D: none (mean contrast only)
#   none (mu only): no stage term
# Tiers are set separately for prevalence (hu) and load (mu), by the
# contrast levels rule (README "Contrast levels"). Nuisance structure
# (README "Nuisance structure"): site-year / event levels a part cannot
# separate are left out (R2); sigma ~ 1 where few sites have positives
# (R3). All read from the job table's structure columns.
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
# 2000 warmup + 4000 sampling per chain = 16000 draws, for every job (round 3). round 2 ran
# most jobs at N_ITER 4000 and only the slow-mixing large species (r. clamitans, r. muscosa,
# r. sierrae) at 6000; one setting for all keeps the fits comparable in Monte Carlo precision.
n_iter <- as.integer(Sys.getenv("N_ITER", "6000"))
n_warmup <- as.integer(Sys.getenv("N_WARMUP", "2000"))
n_cores <- as.integer(Sys.getenv("N_CORES", "4"))
n_threads <- as.integer(Sys.getenv("N_THREADS_PER_CHAIN", "1"))

data_file <- Sys.getenv("DATA_FILE", "data/bd_model_env_unified_2026-11-12.RData")
jobs_file <- Sys.getenv("JOBS_FILE", "data/01_contrast_jobs.csv")
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_contrast")
# loo is off by default: nothing in the contrast analysis uses it (02_contrast_modeling.qmd only
# tabulated Pareto k, and acted on none), and for r. sierrae its pointwise log-likelihood matrix
# (draws x swabs) is a large part of peak memory. set RUN_LOO=TRUE to add it.
run_loo <- as.logical(Sys.getenv("RUN_LOO", "FALSE"))

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

# event = a run of survey dates (of the species) at one site with no gap > event_max_gap days;
# runs longer than event_max_span days are split at their longest gap until each piece fits
event_max_gap = 5
event_max_span = 14

# stage_c: -0.5 adult, +0.5 subadult. the intercept is then the midpoint of the two stages and
# the stage_c coefficient the full subadult - adult contrast; a group's stage_c sd is how much
# the contrast varies at that level.
# pop_year / pop_year_event: explicit ids (rather than population:year:... in the formula) so
# parameter names stay readable.
#
# events follow the survey design rather than the calendar: consecutive-day visits stay together,
# weekly revisits (gaps of 7 days) are separate events, and an intensive resurvey period longer
# than 2 weeks is split where the surveys paused longest. built per species x population from
# the species' own swab dates (each model sees one species). the year is part of the id, so an
# event stays inside its site-year (none crossed the new year in the 2026-11 data).

# split one run of sorted dates (numeric) at its longest gap, recursively, until every piece
# spans <= max_span days. returns piece numbers 1, 2, ...
split_long_run = function(d, max_span) {
  if (length(d) < 2 || max(d) - min(d) + 1 <= max_span) return(rep(1L, length(d)))
  j = which.max(diff(d))
  a = split_long_run(d[1:j], max_span)
  b = split_long_run(d[(j + 1):length(d)], max_span)
  c(a, b + max(a))
}

# event number of each date in a sorted vector of distinct dates (numeric) at one site: a new run
# starts after a gap > max_gap days; long runs are then split by split_long_run
event_ids = function(d, max_gap = event_max_gap, max_span = event_max_span) {
  run = cumsum(c(TRUE, diff(d) > max_gap))
  out = integer(length(d))
  k = 0L
  for (r in unique(run)) {
    i = which(run == r)
    piece = split_long_run(d[i], max_span)
    out[i] = piece + k
    k = k + max(piece)
  }
  out
}

add_contrast_groups = function(df) {
  events = df %>%
    distinct(taxon_capture, population, date) %>%
    arrange(taxon_capture, population, date) %>%
    group_by(taxon_capture, population) %>%
    mutate(event_k = event_ids(as.numeric(date))) %>%
    ungroup()

  df %>%
    left_join(events, by = c("taxon_capture", "population", "date")) %>%
    mutate(life_stage_simple = factor(life_stage_simple, levels = c("adult", "subadult")),
           stage_c = if_else(life_stage_simple == "subadult", 0.5, -0.5),
           pop_year = paste(population, year, sep = "_"),
           pop_year_event = paste(population, year, event_k, sep = "_"))
}
# ============================================================
# END copy of contrast_data_functions.R
# ============================================================

# ============================================================
# BEGIN verbatim copy of size_bd_contrast/contrast_model_functions.R. keep in sync.
# ============================================================
# formula, prior and data-prep builders for the stage 1 contrast models (sourced by
# 02_contrast_modeling.qmd and 03_contrast_analysis.Rmd; copied verbatim into
# r-brms-contained/size_bd_contrast_modeling.R). keep the two in sync. model description:
# README.md ("Stage 1 model", "Model structure rules").
#
# the builders take one job: a row of 01_contrast_jobs.csv (or a list) with
#   tier_hu, tier_mu        "A" / "B" / "C" / "D" (+ "none" for mu; "fixed" = old name for D)
#   hu_has_pop_year, hu_has_event, mu_has_pop_year, mu_has_event
#                           whether that part has the site-year / event level (R2)
#   sigma_by_population     sigma ~ (1 | population) if TRUE, sigma ~ 1 if FALSE (R3)
#
# requires: tidyverse, brms, and contrast_data_functions.R (add_contrast_groups, response_col)

# ---- data ------------------------------------------------------------------------------------
# one species' swabs with stage_c (-0.5 adult / +0.5 subadult), pop_year and pop_year_event
# (survey-run events). droplevels so random-effect terms have no empty levels.
prep_contrast_data = function(brms_data, species_name) {
  brms_data %>%
    filter(taxon_capture == species_name) %>%
    add_contrast_groups() %>%
    droplevels()
}

# ---- structure of one part -------------------------------------------------------------------
# grouping levels present in a part: population always; site-year and event unless dropped (R2)
contrast_levels = function(has_pop_year = TRUE, has_event = TRUE) {
  c("population", if (has_pop_year) "pop_year", if (has_event) "pop_year_event")
}

# levels whose contrast (stage_c) varies, by tier, among the levels present:
#   A: population + site-year + event      B: population + site-year
#   C: population                          D / none: no contrast sd
contrast_slope_groups = function(tier, has_pop_year = TRUE, has_event = TRUE) {
  by_tier = switch(tier,
                   A = c("population", "pop_year", "pop_year_event"),
                   B = c("population", "pop_year"),
                   C = "population",
                   D = character(0),
                   fixed = character(0),   # round 1-2 name for D
                   none = character(0),
                   stop("unknown tier: ", tier))
  intersect(by_tier, contrast_levels(has_pop_year, has_event))
}

# random-effect terms of a part. intercepts vary at every level present (they absorb site, year
# and event-wide swings in Bd, e.g. outbreaks, whether or not both stages were swabbed);
# contrast slopes only at the levels the tier allows.
contrast_re_terms = function(tier, has_pop_year = TRUE, has_event = TRUE) {
  slopes = contrast_slope_groups(tier, has_pop_year, has_event)
  terms = map_chr(contrast_levels(has_pop_year, has_event), function(g) {
    if (g %in% slopes) paste0("(1 + stage_c | ", g, ")") else paste0("(1 | ", g, ")")
  })
  paste(terms, collapse = " + ")
}

# fixed part: a mean contrast unless the part has no stage term ("none")
contrast_fixed_terms = function(tier) if (tier == "none") "1" else "1 + stage_c"

# ---- formula ---------------------------------------------------------------------------------
# hurdle lognormal: mu = mean log load among positives (load part), hu = probability of a zero
# (prevalence part, logit scale), sigma = residual sd of log load among positives (varying by
# population where the data support it, R3). tiers and levels are set separately for mu (load
# counts, positive swabs) and hu (prevalence counts, all swabs).
contrast_formula = function(job, response = response_col) {
  mu_rhs = paste(contrast_fixed_terms(job$tier_mu), "+",
                 contrast_re_terms(job$tier_mu, job$mu_has_pop_year, job$mu_has_event))
  hu_rhs = paste(contrast_fixed_terms(job$tier_hu), "+",
                 contrast_re_terms(job$tier_hu, job$hu_has_pop_year, job$hu_has_event))
  sigma_rhs = if (job$sigma_by_population) "1 + (1 | population)" else "1"

  bf(as.formula(paste(response, "~", mu_rhs)),
     as.formula(paste("hu ~", hu_rhs)),
     as.formula(paste("sigma ~", sigma_rhs)))
}

# ---- priors ----------------------------------------------------------------------------------
# contrast sd priors, one per group with a stage_c slope (brms needs the group named for a
# coefficient-level sd prior). dpar = "" is the mu part.
contrast_sd_priors = function(groups, dpar = "") {
  if (length(groups) == 0) return(NULL)
  do.call(c, lapply(groups, function(g) {
    set_prior("exponential(2)", class = "sd", coef = "stage_c", group = g, dpar = dpar)
  }))
}

# fixed effects and intercept sds as in the unified / seasonal models. the contrast sds get
# exponential(2): tighter than the intercept sds (exponential(1)), so a consistent contrast is the
# default and variation in it must be earned from the data. prior median 0.35, 95% below 1.5
# (log / logit units). brms errors on priors that match no parameter, so each prior is only
# added when its parameter exists (stage terms, contrast sds, correlations, sigma sd).
contrast_priors = function(job) {
  hu_slopes = contrast_slope_groups(job$tier_hu, job$hu_has_pop_year, job$hu_has_event)
  mu_slopes = contrast_slope_groups(job$tier_mu, job$mu_has_pop_year, job$mu_has_event)

  p = c(
    # mu: log load density among positives
    prior(normal(5, 4), class = "Intercept"),
    prior(exponential(1), class = "sd"),
    # hu: logit of non-detection (hu always has a stage term)
    prior(normal(0, 2), class = "Intercept", dpar = "hu"),
    prior(normal(0, 2), class = "b", dpar = "hu"),
    prior(exponential(1), class = "sd", dpar = "hu"),
    contrast_sd_priors(hu_slopes, dpar = "hu"),
    # sigma: log residual sd of log load
    prior(normal(0, 1), class = "Intercept", dpar = "sigma")
  )
  if (job$tier_mu != "none") p = c(p, prior(normal(0, 3), class = "b"))
  p = c(p, contrast_sd_priors(mu_slopes))
  if (job$sigma_by_population) p = c(p, prior(exponential(1), class = "sd", dpar = "sigma"))
  # intercept-contrast correlations exist only where a contrast slope does
  if (length(hu_slopes) + length(mu_slopes) > 0) p = c(p, prior(lkj(2), class = "cor"))
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

# the structure columns (added with rules R1-R3, 2026-11-25) must be present: an older jobs file
# would silently fit the round 1 structure otherwise
structure_cols <- c("tier_hu", "tier_mu", "hu_has_pop_year", "hu_has_event", "mu_has_pop_year",
                    "mu_has_event", "sigma_by_population")
missing_cols <- setdiff(structure_cols, names(contrast_jobs))
if (length(missing_cols) > 0) {
  stop(paste0("jobs file lacks: ", paste(missing_cols, collapse = ", "),
              ". re-run size_bd_contrast/01_contrast_tiers.Rmd and copy the new 01_contrast_jobs.csv"))
}

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

  cat(paste0("levels: hu site-year ", job_row$hu_has_pop_year, ", hu event ", job_row$hu_has_event,
             "; mu site-year ", job_row$mu_has_pop_year, ", mu event ", job_row$mu_has_event,
             "; sigma by population ", job_row$sigma_by_population, "\n\n"))
  model_formula <- contrast_formula(job_row)
  model_priors <- contrast_priors(job_row)
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

  # the fit is saved under the name 02_contrast_modeling.qmd and 03_contrast_analysis.Rmd expect
  # (size_bd_contrast_model_zln_c1_<species>_<variant>.rds). with loo, the raw fit is saved
  # FIRST as *_no_loo.rds so it is never lost to a loo failure, and the with-loo fit then takes
  # the standard name.
  model_file <- file.path(out_dir, paste0("size_bd_contrast_model_zln_c1_", job_row$job, ".rds"))
  if (!run_loo) {
    saveRDS(fit, file = model_file)
    cat(paste0("Model fit saved to:\n\n\t", model_file, "\n\n"))
  } else {
    raw_file <- sub("\\.rds$", "_no_loo.rds", model_file)
    saveRDS(fit, file = raw_file)
    cat(paste0("Model fit (no loo) saved to:\n\n\t", raw_file, "\n\n"))
  }

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
