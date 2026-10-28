# ============================================================
# Run script: per-species size (SVL) GAM ~ Bd load models
# Script version of size_bd_gam/02_size_bd_GAM_modeling.qmd
# (size_bd_load_analysis repo), "Model fitting" section.
#
# Each species gets its own hurdle-lognormal GAM, with a
# multivariate normal prior on the env slopes taken from the
# new-taxon predictive distribution of the unified model (m11d).
#
# Species selection (first match wins):
#   SPECIES=<taxon_capture>    fit just that species
#   SLURM_ARRAY_TASK_ID=<i>    fit the i-th species (0-based, sorted
#                              alphabetically) - see
#                              hpc/run_size_bd_gam_modeling.sbatch
#   neither                    fit all species sequentially
# ============================================================

# load packages
library(brms)
library(cmdstanr)
library(dplyr)
library(tibble)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS size ~ Bd GAM models:\n\n")

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

data_file <- Sys.getenv("DATA_FILE", "data/brms_data_gam_2026-10-28.RData")
# env prior source: the newest m11d co100 fit in models/ (written there by
# hpc/run_size_bd_cutoffs.sbatch; the timestamp in the name sorts chronologically),
# unless PRIOR_MODEL_FILE names one explicitly
prior_model_file <- Sys.getenv("PRIOR_MODEL_FILE", "")
if (prior_model_file == "") {
  prior_candidates <- Sys.glob("models/unified_model_results_zln_m11d_co100_*_with_loo.rds")
  if (length(prior_candidates) == 0) {
    stop("no models/unified_model_results_zln_m11d_co100_*_with_loo.rds found: ",
         "fit co100 first (sbatch --array=2 hpc/run_size_bd_cutoffs.sbatch) or set PRIOR_MODEL_FILE")
  }
  prior_model_file <- sort(prior_candidates, decreasing = TRUE)[1]
}
out_dir <- Sys.getenv("OUT_DIR", "models/size_bd_gam")
run_loo <- as.logical(Sys.getenv("RUN_LOO", "TRUE"))

# ------------------------------------------------------------
# Load data, choose species
# ------------------------------------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat(paste0("Loading data from: ", data_file, "\n\n"))
load(data_file)

species_list <- sort(unique(as.character(brms_data$taxon_capture)))

species_env <- Sys.getenv("SPECIES", "")
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (species_env != "") {
  if (!species_env %in% species_list) {
    stop(paste0("SPECIES '", species_env, "' not in data. Available: ", paste(species_list, collapse = ", ")))
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
# New-taxon predictive prior for the env slopes (precip, temp, temp^2)
# ------------------------------------------------------------
# per draw: fixed effect + a new taxon offset drawn from the m11d taxon RE distribution,
# N(0, Sigma) with Sigma from sd_taxon_capture and cor_taxon_capture. this uses no
# species-specific r_taxon_capture, so each species' data is not counted twice
# (once in the prior, again in its GAM). the prior is therefore the same for all species.
env_terms <- c("precip_mm_10_z", "temp_c_30_z", "Itemp_c_30_zE2")

new_taxon_env_draws <- function(model, dpar_prefix = "") {
  draws <- as.data.frame(model)
  terms <- paste0(dpar_prefix, env_terms)

  b <- as.matrix(draws[, paste0("b_", terms)])
  sds <- as.matrix(draws[, paste0("sd_taxon_capture__", terms)])

  offsets <- t(sapply(seq_len(nrow(draws)), function(i) {
    R <- diag(length(terms))
    for (j in 1:(length(terms) - 1)) {
      for (k in (j + 1):length(terms)) {
        R[j, k] <- R[k, j] <- draws[[paste0("cor_taxon_capture__", terms[j], "__", terms[k])]][i]
      }
    }
    Sigma <- diag(sds[i, ]) %*% R %*% diag(sds[i, ])
    MASS::mvrnorm(1, mu = rep(0, length(terms)), Sigma = Sigma)
  }))

  total <- b + offsets
  colnames(total) <- c("precip_linear", "temp_linear", "temp_quad")
  as_tibble(total)
}

cat(paste0("Loading unified model for env prior from: ", prior_model_file, "\n\n"))
unified_model <- readRDS(prior_model_file)

set.seed(124)
new_taxon_hu <- new_taxon_env_draws(unified_model, "hu_")
new_taxon_b <- new_taxon_env_draws(unified_model, "")

# free the unified model - only the prior draws are needed from here on
rm(unified_model)
invisible(gc())

cat("Env prior (new-taxon predictive) summary:\n")
print(bind_rows(hu = summarise(new_taxon_hu, across(everything(), list(mean = mean, sd = sd))),
                mu = summarise(new_taxon_b, across(everything(), list(mean = mean, sd = sd))),
                .id = "dpar"))
cat("\n")

# ------------------------------------------------------------
# Model structure
# ------------------------------------------------------------
# multivariate normal prior on the env slopes [precip, temp, temp^2], from the
# new-taxon predictive draws (keeps their correlations; independent normals do not).
# b and b_hu are ordered [precip_mm_10_z, temp_c_30_z, Itemp_c_30_zE2]: the intercept
# is separate (centered) and smooth fixed parts are in bs, so b[1:3] is the env block.
# added to lprior (not target) so priorsense sees it; the b coefs keep brms' default
# flat prior, so this is their only prior. same prior for all species.
build_env_stanvars <- function(new_taxon_b, new_taxon_hu) {
  stanvar(x = colMeans(new_taxon_b), name = "env_prior_mean_b") +
    stanvar(x = cov(new_taxon_b), name = "env_prior_cov_b") +
    stanvar(x = colMeans(new_taxon_hu), name = "env_prior_mean_hu") +
    stanvar(x = cov(new_taxon_hu), name = "env_prior_cov_hu") +
    stanvar(scode = paste(
      "  lprior += multi_normal_lpdf(b[1:3] | env_prior_mean_b, env_prior_cov_b);",
      "  lprior += multi_normal_lpdf(b_hu[1:3] | env_prior_mean_hu, env_prior_cov_hu);",
      sep = "\n"),
      block = "tparameters", position = "end")
}

env_stanvars <- build_env_stanvars(new_taxon_b, new_taxon_hu)

# guard: the stanvars assume b[1:3] / b_hu[1:3] are the env terms in this order
check_env_order <- function(model_formula, model_data) {
  sdata <- make_standata(model_formula, data = model_data, family = hurdle_lognormal())
  for (x in c("X", "X_hu")) {
    cols <- setdiff(colnames(sdata[[x]]), "Intercept")
    if (!identical(cols[1:3], env_terms)) {
      stop(paste0(x, " columns are not [", paste(env_terms, collapse = ", "), "]: ", paste(cols, collapse = ", ")))
    }
  }
}

# size structure: metamorph baseline (intercept) + SVL smooth over non-metamorphs.
# numeric-by smooths are uncentered, so the smooth's constant is the non-metamorph
# offset; adding is_metamorph as well (as in m2) is rank deficient.
# species without metamorphs (is_not_metamorph all 1) get a plain centered smooth.
# (1 | kg_code_3) dropped: species span only 2-4 KG zones; population REs absorb site climate.
build_spp_formula <- function(model_data) {
  svl_term <- if (any(model_data$is_metamorph == 1)) {
    "s(svl_mm, by = is_not_metamorph)"
  } else {
    "s(svl_mm)"
  }

  rhs <- paste("1 +",
               svl_term,
               "+ precip_mm_10_z + temp_c_30_z + I(temp_c_30_z^2)",
               "+ (1 | population:life_stage_simple)",
               "+ (0 + life_stage_simple | population:year)")

  bf(as.formula(paste("bd_load_co100_density ~", rhs)),
     as.formula(paste("hu ~", rhs)),
     sigma ~ 1 + (1 | population))
}

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------
model_svl_by_species <- function(species_name) {
  cat(paste0("==== Modeling Bd vs SVL for: ", species_name, " ====\n\n"))

  # drop unused life stages (e.g. no metamorphs) so RE terms have no empty columns
  model_data <- droplevels(filter(brms_data, taxon_capture == species_name))
  cat(paste0("n = ", nrow(model_data), "\n\n"))

  m3_formula <- build_spp_formula(model_data)
  print(m3_formula)
  check_env_order(m3_formula, model_data)

  # refresh: print progress regularly, flushed live to the log (see stdbuf
  # in hpc/*.sbatch) so a running job's pace is visible via `tail -f`
  fit_start_time <- Sys.time()
  m3_zln <- brm(
    formula  = m3_formula,
    data     = model_data,
    family   = hurdle_lognormal(),
    stanvars = env_stanvars,
    control  = list(adapt_delta = 0.95, max_treedepth = 12),
    chains   = n_chains,
    iter     = n_iter,
    warmup   = n_warmup,
    seed     = 124,
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

  # save raw fit FIRST so it's never lost to a loo failure (loo has a history
  # of exhausting RAM). the with-loo fit is saved under the name the
  # analysis notebooks expect (size_bd_gam_model_zln_m3_<species>.rds).
  model_file <- file.path(out_dir, paste0("size_bd_gam_model_zln_m3_", species_name, ".rds"))
  raw_file <- sub("\\.rds$", "_no_loo.rds", model_file)
  saveRDS(m3_zln, file = raw_file)
  cat(paste0("Model fit (no loo) saved to:\n\n\t", raw_file, "\n\n"))

  if (run_loo) {
    tryCatch({
      cat("Adding LOO...\n")
      m3_zln <- add_criterion(m3_zln, "loo")
      saveRDS(m3_zln, file = model_file)
      cat(paste0("Model fit with loo saved to:\n\n\t", model_file, "\n\n"))
    }, error = function(e) {
      cat("(could not compute/save loo: ", conditionMessage(e), ")\n\n")
    })
  }

  cat(paste0("Done: ", species_name, "\n\n"))
}

for (sp in run_species) {
  model_svl_by_species(sp)
}
