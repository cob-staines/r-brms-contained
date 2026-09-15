# ============================================================
# Run script: size ~ Bd load model
# Structure and execution are split into sections below: edit
# "Model structure" to change what the model is, edit
# "Run parameters" / the env vars passed into this script to
# change how it's run (locally, as a timing probe, or full-scale
# on the cluster — see hpc/*.sbatch).
# ============================================================

# load packages
library(brms)
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS model:\n\n")

# ------------------------------------------------------------
# Model structure: family, formula, priors
# ------------------------------------------------------------

model_family <- hurdle_lognormal()

model_formula <- bf(
  bd_load_co_density ~ 1 +
    life_stage_simple * taxon_capture +
    (precip_mm_10_z  + temp_c_30_z + I(temp_c_30_z^2)) +
    (0 + precip_mm_10_z + temp_c_30_z + I(temp_c_30_z^2) | taxon_capture) +
    (1 | kg_code_3) +
    (1 | gr(taxon_capture:population:life_stage_simple, by = taxon_capture)) +
    (0 + life_stage_simple | gr(taxon_capture:population:year, by = taxon_capture)),

  hu ~ 1 +
    life_stage_simple * taxon_capture +
    (precip_mm_10_z + temp_c_30_z + I(temp_c_30_z^2)) +
    (0 + precip_mm_10_z + temp_c_30_z + I(temp_c_30_z^2) | taxon_capture) +
    (1 | kg_code_3) +
    (1 | gr(taxon_capture:population:life_stage_simple, by = taxon_capture)) +
    (0 + life_stage_simple | gr(taxon_capture:population:year, by = taxon_capture)),

  sigma ~ 1 + taxon_capture + (1 | taxon_capture:population)
)

model_priors <- c(
  # b
  prior(normal(5, 4), class = "Intercept"),
  prior(normal(0, 3),   class = "b"),

  # hu
  prior(normal(0, 2), class = "Intercept", dpar = "hu"),
  prior(normal(0, 2), class = "b", dpar = "hu"),

  # Sigma
  prior(normal(0, 1), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.3), class = "b", dpar = "sigma"),

  # Random effects
  prior(exponential(1), class = "sd")
)

# ------------------------------------------------------------
# Run parameters
# ------------------------------------------------------------
# Override via env vars, e.g. for a scaled-down timing probe or to enable
# within-chain multithreading (see brms threading vignette:
# https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html).
# n_cores below is chains run in parallel; total CPUs actually used is
# n_cores * n_threads, so request that many in --cpus-per-task/-c in Slurm.
n_chains <- as.integer(Sys.getenv("N_CHAINS", "4"))
n_iter <- as.integer(Sys.getenv("N_ITER", "4000"))
n_warmup <- as.integer(Sys.getenv("N_WARMUP", "2000"))
n_cores <- as.integer(Sys.getenv("N_CORES", "4"))
n_threads <- as.integer(Sys.getenv("N_THREADS_PER_CHAIN", "1"))

# ------------------------------------------------------------
# Execution
# ------------------------------------------------------------

# create output folder
if (!dir.exists("models")) dir.create("models")

# load data
load("data/bd_model_env_unified_2026-09-11.RData")

# fit model
# refresh: print progress every ~50 iterations (min 1), flushed live to the
# log (see stdbuf in hpc/*.sbatch) so a running job's pace is visible via
# `tail -f` rather than only knowable after it finishes or times out.
fit_start_time <- Sys.time()
m11b_zln <- brm(
  formula = model_formula,
  data = brms_data,
  family = model_family,
  prior = model_priors,
  sample_prior = "yes",
  control = list(adapt_delta = 0.98),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  seed = 126,
  cores = n_cores,
  threads = threading(n_threads),
  backend = "cmdstanr",
  refresh = max(1, n_iter %/% 50)
)
cat(paste0(
  "\nFit time (", n_chains, " chains x ", n_threads, " threads/chain, ",
  n_iter, " iter, ", n_warmup, " warmup): ",
  format(Sys.time() - fit_start_time), "\n\n"
))
# warmup vs sampling time per chain, since warmup (adaptation) cost often
# doesn't scale linearly with warmup length - useful for judging whether a
# short probe's timing is representative of the real warmup length
cat("Per-chain warmup/sampling breakdown:\n")
print(m11b_zln$fit$time())
cat("\n")

# save model
model_fit_file = paste0("models/unified_model_results_zln_life_tax_pt2_rgrtaxacorpt2poplsyrkg3_p10t30.rds", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(m11b_zln, file = model_fit_file)

cat(paste0("Done. Model fit saved to:\n\n\t", model_fit_file, "\n\n"))
