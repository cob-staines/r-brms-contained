# ============================================================
# Run script: GAM of Bd by body size, one per species
# (size_bd_contrast GAM option 1). Script version of
# size_bd_contrast/06_gam_modeling.qmd (size_bd_load_analysis repo),
# "Fit" section. Model description: size_bd_contrast/README.md
# ("Downstream models" -> "GAM"). Replaces size_bd_gam_modeling.R
# (m3: m11g climate prior), kept until that model is retired.
#
# Each species gets a hurdle-lognormal model with a thin plate SVL
# smooth (plus a metamorph baseline) and the stage 1 nuisance
# intercepts (population, site-year, survey event, as the species'
# stage 1 model has them); no weather or climate terms, no m11g prior.
#
# Inputs, both written by 06_gam_modeling.qmd (copy to data/):
#   DATA_FILE     data/06_gam_model_data.rds   GAM data with event ids
#   SPECIES_FILE  data/06_gam_species.csv      species and their stage 1
#                                              structure flags (R2 / R3)
#
# Species selection (first match wins):
#   SPECIES=<taxon_capture>    fit just that species
#   SLURM_ARRAY_TASK_ID=<i>    fit the i-th species (0-based, sorted
#                              alphabetically) - see
#                              hpc/run_size_bd_contrast_gam_modeling.sbatch
#   neither                    fit all species sequentially
# ============================================================

# load packages. tidyverse comes with the rocker/tidyverse base image.
library(tidyverse)
library(brms)
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS GAMs of Bd by body size (size_bd_contrast, GAM option 1):\n\n")

# ============================================================
# constant from size_bd_contrast/contrast_data_functions.R (the GAM data already carry the
# event ids, so nothing else from it is needed). keep in sync.
# ============================================================
response_col = "bd_load_co100_density"

# ============================================================
# BEGIN verbatim copy of size_bd_contrast/gam_model_functions.R. keep in sync (diff the two).
# ============================================================
# formula and prior builders for the GAM (sourced by 06_gam_modeling.qmd; copied verbatim into
# r-brms-contained/size_bd_contrast_gam_modeling.R). keep the two in sync. model description:
# README.md ("Downstream models" -> "GAM").
#
# the builders take one species row of 06_gam_species.csv (or a list) with the stage 1 structure
# flags: hu_has_pop_year, hu_has_event, mu_has_pop_year, mu_has_event, sigma_by_population.
#
# requires: tidyverse, brms, and response_col (contrast_data_functions.R)

# ---- formula ---------------------------------------------------------------------------------
# per part (mu = log load among positives, hu = non-detection): one SVL smooth + the stage 1
# nuisance intercepts of that part.
#
# size: metamorph baseline (the intercept) + a thin plate smooth over subadults and adults.
# a numeric-by smooth is not centered, so it carries the post-metamorph offset itself; adding
# is_metamorph as well is rank deficient. species without metamorphs get a plain (centered)
# s(svl_mm), and the intercept is then the mean level.
#
# nuisance: (1 | population), plus (1 | pop_year) and (1 | pop_year_event) where stage 1 kept
# them for that part (R2 flags); sigma varies by population where stage 1's does (R3).
gam_formula = function(sp_row, model_data, response = response_col) {
  svl_term = if (any(model_data$life_stage_simple == "metamorph")) {
    "s(svl_mm, by = is_not_metamorph)"
  } else {
    "s(svl_mm)"
  }
  part_rhs = function(has_pop_year, has_event) {
    groups = c("population", if (has_pop_year) "pop_year", if (has_event) "pop_year_event")
    paste(c("1", svl_term, paste0("(1 | ", groups, ")")), collapse = " + ")
  }
  mu_rhs = part_rhs(sp_row$mu_has_pop_year, sp_row$mu_has_event)
  hu_rhs = part_rhs(sp_row$hu_has_pop_year, sp_row$hu_has_event)
  sigma_rhs = if (sp_row$sigma_by_population) "1 + (1 | population)" else "1"

  bf(as.formula(paste(response, "~", mu_rhs)),
     as.formula(paste("hu ~", hu_rhs)),
     as.formula(paste("sigma ~", sigma_rhs)))
}

# ---- priors ----------------------------------------------------------------------------------
# the stage 1 nuisance priors (contrast_priors() in contrast_model_functions.R, without the stage
# and contrast-sd terms): intercepts, intercept sds, sigma. the smooth terms keep brms' defaults
# (flat on its linear part `b` (class bs in the Stan code), student_t(3, 0, 2.5) on its
# wiggliness `sds`), as in the earlier GAM fits.
# note: with a metamorph baseline, the Intercept is the metamorph level; without metamorphs it
# is the mean level. normal(5, 4) / normal(0, 2) are wide for either.
gam_priors = function(sp_row) {
  p = c(
    # mu: log load density among positives
    prior(normal(5, 4), class = "Intercept"),
    prior(exponential(1), class = "sd"),
    # hu: logit of non-detection
    prior(normal(0, 2), class = "Intercept", dpar = "hu"),
    prior(exponential(1), class = "sd", dpar = "hu"),
    # sigma: log residual sd of log load
    prior(normal(0, 1), class = "Intercept", dpar = "sigma")
  )
  if (sp_row$sigma_by_population) p = c(p, prior(exponential(1), class = "sd", dpar = "sigma"))
  p
}

# ---- run settings ----------------------------------------------------------------------------
# shared by the local fits (06) and the cluster script, so both give the same model and draws:
# 4 chains x (2000 warmup + 2000 sampling) = 8000 draws
gam_n_chains = 4
gam_n_iter = 4000
gam_n_warmup = 2000
gam_adapt_delta = 0.95
gam_max_treedepth = 12
gam_seed = 124

# model file name (without .rds) for one species: brm(file = ) adds .rds
gam_model_name = function(sp) paste0("gam_m4_", sp)
# ============================================================
# END copy of size_bd_contrast/gam_model_functions.R
# ============================================================

# ------------------------------------------------------------
# Run parameters
# ------------------------------------------------------------
# defaults are the shared run settings above (same model and draws as a local fit in 06).
# n_cores is chains run in parallel; total CPUs used is n_cores * n_threads, so request that
# many with -c in Slurm.
n_chains <- as.integer(Sys.getenv("N_CHAINS", gam_n_chains))
n_iter <- as.integer(Sys.getenv("N_ITER", gam_n_iter))
n_warmup <- as.integer(Sys.getenv("N_WARMUP", gam_n_warmup))
n_cores <- as.integer(Sys.getenv("N_CORES", n_chains))
n_threads <- as.integer(Sys.getenv("N_THREADS_PER_CHAIN", "1"))
# sampler target acceptance rate: smaller steps for tight curvature (many small groups), at
# ~1.5-2x run time. changes how the posterior is explored, not the model.
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", gam_adapt_delta))

data_file <- Sys.getenv("DATA_FILE", "data/06_gam_model_data.rds")
species_file <- Sys.getenv("SPECIES_FILE", "data/06_gam_species.csv")
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_contrast_gam")

# ------------------------------------------------------------
# Load data, choose species
# ------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat(paste0("Loading data from: ", data_file, "\n"))
gam_data <- readRDS(data_file)
cat(paste0("Loading species table from: ", species_file, "\n\n"))
gam_species_table <- read_csv(species_file, show_col_types = FALSE)

species_list <- sort(gam_species_table$taxon_capture)

# the data and the species table must come from the same run of 06
missing_data <- setdiff(species_list, unique(gam_data$taxon_capture))
if (length(missing_data) > 0) stop("species in SPECIES_FILE but not in DATA_FILE: ", paste(missing_data, collapse = ", "))
if (!all(c("pop_year", "pop_year_event") %in% names(gam_data))) stop("DATA_FILE has no pop_year / pop_year_event: use 06_gam_model_data.rds")

species_env <- Sys.getenv("SPECIES", "")
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (species_env != "") {
  if (!species_env %in% species_list) {
    stop(paste0("SPECIES '", species_env, "' not in the species table. Available: ", paste(species_list, collapse = ", ")))
  }
  run_species <- species_env
} else if (array_task_id != "") {
  i <- as.integer(array_task_id) + 1
  if (i > length(species_list)) {
    stop(paste0("SLURM_ARRAY_TASK_ID ", array_task_id, " out of range: only ", length(species_list), " species (0-", length(species_list) - 1, ")"))
  }
  run_species <- species_list[i]
} else {
  run_species <- species_list
}

cat("All species (array index: taxon_capture):\n")
cat(paste0("  ", seq_along(species_list) - 1, ": ", species_list, collapse = "\n"), "\n\n")
cat(paste0("Fitting: ", paste(run_species, collapse = ", "), "\n\n"))

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------
fit_gam_species <- function(sp) {
  cat(paste0("==== GAM of Bd by SVL for: ", sp, " ====\n\n"))

  sp_row <- filter(gam_species_table, taxon_capture == sp)
  model_data <- droplevels(filter(gam_data, taxon_capture == sp))
  model_formula <- gam_formula(sp_row, model_data)
  model_prior <- gam_priors(sp_row)

  cat(paste0("n = ", nrow(model_data),
             " (metamorphs ", sum(model_data$life_stage_simple == "metamorph"),
             "), events = ", n_distinct(model_data$pop_year_event), "\n\n"))
  print(model_formula)
  print(model_prior)

  # refresh: print progress regularly, flushed live to the log (see stdbuf in the sbatch) so a
  # running job's pace is visible via `tail -f`
  fit_start_time <- Sys.time()
  fit <- brm(
    formula  = model_formula,
    data     = model_data,
    family   = hurdle_lognormal(),
    prior    = model_prior,
    control  = list(adapt_delta = adapt_delta, max_treedepth = gam_max_treedepth),
    chains   = n_chains,
    iter     = n_iter,
    warmup   = n_warmup,
    seed     = gam_seed,
    cores    = n_cores,
    threads  = threading(n_threads),
    backend  = "cmdstanr",
    refresh  = max(1, n_iter %/% 50)
  )
  cat(paste0(
    "\nFit time (", n_chains, " chains x ", n_threads, " threads/chain, ",
    n_iter, " iter, ", n_warmup, " warmup, adapt_delta ", adapt_delta, "): ",
    format(Sys.time() - fit_start_time), "\n"
  ))
  np <- nuts_params(fit)
  cat(paste0("Divergences: ", sum(np$Value[np$Parameter == "divergent__"]), "\n\n"))

  # saved under the name 06 / 07 expect: copy into size_bd_contrast/models/
  model_file <- file.path(out_dir, paste0(gam_model_name(sp), ".rds"))
  saveRDS(fit, file = model_file)
  cat(paste0("Model fit saved to:\n\n\t", model_file, "\n\n"))
}

for (sp in run_species) {
  fit_gam_species(sp)
}
