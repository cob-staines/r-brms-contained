# ============================================================
# Run script: per-species seasonal stage x Bd models, all species
# Script version of size_bd_season/02_season_all_species.qmd
# (size_bd_load_analysis repo), "Model fitting" section.
# Keep data prep, formula, priors and settings in sync with it.
#
# Each species gets the same hurdle-lognormal model: a seasonal smooth
# shared by both life stages plus a shrinkage DIFFERENCE smooth
# (subadult - adult), which collapses to a constant stage difference
# when the data cannot support a changing one. No weather terms; the
# same random-effect structure as the GAM / stage-curve models.
#
# Species, order and seasonal origins come from SPECIES_FILE
# (01_season_species_all.csv, written by size_bd_season/01_season_screen.Rmd).
# Copy it to data/ alongside the model data.
#
# Species selection (first match wins):
#   SPECIES=<taxon_capture>    fit just that species
#   SLURM_ARRAY_TASK_ID=<i>    fit the species with array_index == i
#                              (alphabetical, 0-based) - see
#                              hpc/run_size_bd_season_all_modeling.sbatch
#   neither                    fit all species sequentially
# ============================================================

# load packages
library(brms)
library(cmdstanr)
library(dplyr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS seasonal stage x Bd models (all species):\n\n")

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
species_file <- Sys.getenv("SPECIES_FILE", "data/01_season_species_all.csv")
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_season_all")
run_loo <- as.logical(Sys.getenv("RUN_LOO", "TRUE"))

# sampler target acceptance rate. 0.95 matches 02_season_all_species.qmd; raise it (e.g. 0.99)
# when resubmitting species with divergences: smaller steps follow tight curvature (such as a
# smooth sd shrinking toward 0) better, at the cost of more steps per iteration (slower)
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.95"))

# ------------------------------------------------------------
# Load data, choose species
# ------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat(paste0("Loading species list from: ", species_file, "\n"))
season_species_all <- read.csv(species_file, stringsAsFactors = FALSE)
season_species_all <- season_species_all[order(season_species_all$array_index), ]
season_species <- season_species_all$taxon_capture

cat(paste0("Loading data from: ", data_file, "\n\n"))
load(data_file)

missing_species <- setdiff(season_species, unique(as.character(brms_data$taxon_capture)))
if (length(missing_species) > 0) {
  stop(paste0("species not in data: ", paste(missing_species, collapse = ", ")))
}

species_env <- Sys.getenv("SPECIES", "")
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (species_env != "") {
  if (!species_env %in% season_species) {
    stop(paste0("SPECIES '", species_env, "' not in species file. Available: ", paste(season_species, collapse = ", ")))
  }
  run_species <- species_env
} else if (array_task_id != "") {
  i <- as.integer(array_task_id)
  if (!i %in% season_species_all$array_index) {
    stop(paste0("SLURM_ARRAY_TASK_ID ", array_task_id, " out of range: array_index runs 0-", max(season_species_all$array_index)))
  }
  run_species <- season_species_all$taxon_capture[season_species_all$array_index == i]
} else {
  run_species <- season_species
}

cat("All species (array index: taxon_capture):\n")
cat(paste0("  ", season_species_all$array_index, ": ", season_species, collapse = "\n"), "\n\n")
cat(paste0("Fitting: ", paste(run_species, collapse = ", "), "\n\n"))

# season covariate: months since the species' seasonal origin (middle of its longest sampling
# gap). stage_o: ordered copy of life stage, so by = stage_o builds only the subadult difference
# smooth.
prep_season_all_data <- function(brms_data, species_name, origin_doy) {
  brms_data %>%
    filter(taxon_capture == species_name) %>%
    mutate(life_stage_simple = factor(life_stage_simple, levels = c("adult", "subadult")),
           stage_o = ordered(life_stage_simple, levels = c("adult", "subadult")),
           doy = as.integer(format(date, "%j")),
           season_m = ((doy - origin_doy) %% 365) / 30.44) %>%
    # drop unused levels so random-effect terms have no empty columns
    droplevels()
}

# ------------------------------------------------------------
# Model structure (see 02_season_all_species.qmd for the reasoning)
# ------------------------------------------------------------
season_all_formula <- bf(
  bd_load_co100_density ~ 1 + life_stage_simple +
    s(season_m, bs = "ts", k = 5) +
    s(season_m, by = stage_o, bs = "ts", k = 5) +
    (1 | population:life_stage_simple) +
    (0 + life_stage_simple | population:year),

  hu ~ 1 + life_stage_simple +
    s(season_m, bs = "ts", k = 5) +
    s(season_m, by = stage_o, bs = "ts", k = 5) +
    (1 | population:life_stage_simple) +
    (0 + life_stage_simple | population:year),

  sigma ~ 1 + (1 | population)
)

sds_shared <- 's(season_m, bs = "ts", k = 5)'
sds_diff <- 's(season_m, by = stage_o, bs = "ts", k = 5)'

season_all_priors <- c(
  # mu: log load density among positives
  prior(normal(5, 4), class = "Intercept"),
  prior(normal(0, 3), class = "b"),                                  # subadult offset
  set_prior("student_t(3, 0, 2.5)", class = "sds", coef = sds_shared),
  set_prior("student_t(3, 0, 1)", class = "sds", coef = sds_diff),
  prior(exponential(1), class = "sd"),

  # hu: logit of non-detection
  prior(normal(0, 2), class = "Intercept", dpar = "hu"),
  prior(normal(0, 2), class = "b", dpar = "hu"),                     # subadult offset
  set_prior("student_t(3, 0, 2.5)", class = "sds", coef = sds_shared, dpar = "hu"),
  set_prior("student_t(3, 0, 1)", class = "sds", coef = sds_diff, dpar = "hu"),
  prior(exponential(1), class = "sd", dpar = "hu"),

  # sigma: log residual sd of log load
  prior(normal(0, 1), class = "Intercept", dpar = "sigma"),
  prior(exponential(1), class = "sd", dpar = "sigma")
)

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------
model_season_all_by_species <- function(species_name) {
  cat(paste0("==== Modeling seasonal stage x Bd (all-species form) for: ", species_name, " ====\n\n"))

  origin_doy <- season_species_all$origin_doy[season_species_all$taxon_capture == species_name]
  model_data <- prep_season_all_data(brms_data, species_name, origin_doy)
  cat(paste0("n = ", nrow(model_data), ", seasonal origin day of year = ", origin_doy,
             ", adapt_delta = ", adapt_delta, "\n\n"))
  print(season_all_formula)

  # refresh: print progress regularly, flushed live to the log (see stdbuf
  # in hpc/*.sbatch) so a running job's pace is visible via `tail -f`
  fit_start_time <- Sys.time()
  fit <- brm(
    formula  = season_all_formula,
    data     = model_data,
    family   = hurdle_lognormal(),
    prior    = season_all_priors,
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
  # name the analysis notebooks expect (size_bd_season_all_model_zln_d1_<species>.rds).
  model_file <- file.path(out_dir, paste0("size_bd_season_all_model_zln_d1_", species_name, ".rds"))
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

  cat(paste0("Done: ", species_name, "\n\n"))
}

for (sp in run_species) {
  model_season_all_by_species(sp)
}
