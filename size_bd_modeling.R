# ============================================================
# brms test script: Bayesian regression on the 'cars' dataset
# cars: 50 observations of speed (mph) and stopping distance (ft)
# Model: dist ~ speed (simple linear regression, Bayesian)
# ============================================================

# load packages
library(brms)
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning BRMS model:\n\n")

# create output folder
if (!dir.exists("models")) dir.create("models")

# load data
load("data/bd_model_env_unified_2026-04-22.RData")

# set prior
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

# fit model
# isolate scales between taxa
m11c_zln <- brm(
  formula = bf(
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
    
    sigma ~ 1 + taxon_capture + (1 | gr(taxon_capture:population, by = taxon_capture))
  ),
  data = brms_data,
  family = hurdle_lognormal(),
  prior = model_priors,
  sample_prior = "yes",
  control = list(adapt_delta = 0.98),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 126,
  cores = 4,
  backend = "cmdstanr"
)

# save model
model_fit_file = paste0("models/unified_model_results_zln_life_tax_pt2_rgrtaxacorpt2poplsyrkg3tax_p10t30_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(m11c_zln, file = model_fit_file)

cat(paste0("Done. Model fit saved to:\n\n\t", model_fit_file, "\n\n"))
