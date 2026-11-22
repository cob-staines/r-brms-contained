# ============================================================
# Run script: per-species seasonal stage x Bd models
# Script version of size_bd_season/02_season_modeling.qmd
# (size_bd_load_analysis repo), "Model fitting" section.
# Keep data prep, formula, priors and settings in sync with it.
#
# Each species gets its own hurdle-lognormal model with a separate
# seasonal smooth per life stage (stage x season), no weather terms,
# and the same random-effect structure as the GAM models.
#
# Species selection (first match wins):
#   SPECIES=<taxon_capture>    fit just that species
#   SLURM_ARRAY_TASK_ID=<i>    fit the i-th species (0-based, alphabetical
#                              order of season_species below) - see
#                              hpc/run_size_bd_season_modeling.sbatch
#   neither                    fit all species sequentially
# ============================================================

# load packages
library(brms)
library(cmdstanr)
library(dplyr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS seasonal stage x Bd models:\n\n")

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
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_season")
run_loo <- as.logical(Sys.getenv("RUN_LOO", "TRUE"))

# species with enough co-sampled weeks of both life stages (screen tier 1, 10 / 6), minus the
# two-window panamanian species. alphabetical: array index 0-5
season_species <- sort(c("anaxyrus_canorus",
                         "rana_cascadae",
                         "rana_catesbeiana",
                         "rana_clamitans",
                         "rana_muscosa",
                         "rana_sierrae"))

# ------------------------------------------------------------
# Load data, choose species
# ------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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
    stop(paste0("SPECIES '", species_env, "' not in season_species. Available: ", paste(season_species, collapse = ", ")))
  }
  run_species <- species_env
} else if (array_task_id != "") {
  i <- as.integer(array_task_id) + 1
  if (i > length(season_species)) {
    stop(paste0("SLURM_ARRAY_TASK_ID ", array_task_id, " out of range: only ", length(season_species), " species (0-", length(season_species) - 1, ")"))
  }
  run_species <- season_species[i]
} else {
  run_species <- season_species
}

cat("All species (array index: taxon_capture):\n")
cat(paste0("  ", seq_along(season_species) - 1, ": ", season_species, collapse = "\n"), "\n\n")
cat(paste0("Fitting: ", paste(run_species, collapse = ", "), "\n\n"))

# season covariate: months from july 1 (day 182). centering and scaling only changes the units
# of the smooth's linear (unpenalized) part, which makes its prior easier to set.
prep_season_data <- function(brms_data, species_name) {
  brms_data %>%
    filter(taxon_capture == species_name) %>%
    mutate(life_stage_simple = factor(life_stage_simple, levels = c("adult", "subadult")),
           doy = as.integer(format(date, "%j")),
           season_m = (doy - 182) / 30.44) %>%
    # drop unused levels so random-effect terms have no empty columns
    droplevels()
}

# ------------------------------------------------------------
# Model structure (see 02_season_modeling.qmd for the reasoning)
# ------------------------------------------------------------
season_formula <- bf(
  bd_load_co100_density ~ 1 + life_stage_simple +
    s(season_m, by = life_stage_simple, k = 5) +
    (1 | population:life_stage_simple) +
    (0 + life_stage_simple | population:year),

  hu ~ 1 + life_stage_simple +
    s(season_m, by = life_stage_simple, k = 5) +
    (1 | population:life_stage_simple) +
    (0 + life_stage_simple | population:year),

  sigma ~ 1 + (1 | population)
)

season_priors <- c(
  # mu: log load density among positives
  prior(normal(5, 4), class = "Intercept"),
  prior(normal(0, 3), class = "b"),                                              # subadult offset
  prior(normal(0, 2), class = "b", coef = "sseason_m:life_stage_simpleadult_1"),
  prior(normal(0, 2), class = "b", coef = "sseason_m:life_stage_simplesubadult_1"),
  prior(student_t(3, 0, 2.5), class = "sds"),
  prior(exponential(1), class = "sd"),

  # hu: logit of non-detection
  prior(normal(0, 2), class = "Intercept", dpar = "hu"),
  prior(normal(0, 2), class = "b", dpar = "hu"),                                 # subadult offset
  prior(normal(0, 2), class = "b", coef = "sseason_m:life_stage_simpleadult_1", dpar = "hu"),
  prior(normal(0, 2), class = "b", coef = "sseason_m:life_stage_simplesubadult_1", dpar = "hu"),
  prior(student_t(3, 0, 2.5), class = "sds", dpar = "hu"),
  prior(exponential(1), class = "sd", dpar = "hu"),

  # sigma: log residual sd of log load
  prior(normal(0, 1), class = "Intercept", dpar = "sigma"),
  prior(exponential(1), class = "sd", dpar = "sigma")
)

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------
model_season_by_species <- function(species_name) {
  cat(paste0("==== Modeling seasonal stage x Bd for: ", species_name, " ====\n\n"))

  model_data <- prep_season_data(brms_data, species_name)
  cat(paste0("n = ", nrow(model_data), "\n\n"))
  print(season_formula)

  # refresh: print progress regularly, flushed live to the log (see stdbuf
  # in hpc/*.sbatch) so a running job's pace is visible via `tail -f`
  fit_start_time <- Sys.time()
  fit <- brm(
    formula  = season_formula,
    data     = model_data,
    family   = hurdle_lognormal(),
    prior    = season_priors,
    control  = list(adapt_delta = 0.95, max_treedepth = 12),
    chains   = n_chains,
    iter     = n_iter,
    warmup   = n_warmup,
    seed     = 127,
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
  # name the analysis notebooks expect (size_bd_season_model_zln_s1_<species>.rds).
  model_file <- file.path(out_dir, paste0("size_bd_season_model_zln_s1_", species_name, ".rds"))
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
  model_season_by_species(sp)
}
