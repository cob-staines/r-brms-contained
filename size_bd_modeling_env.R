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
if (!dir.exists("models/env")) dir.create("models/env")

# load data
load("data/bd_model_env_unified_2026-05-30.RData")

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

# model fit function
run_zln_model = function(f_formula, file_base) {
  m_zln <- brm(
    formula = f_formula,
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
  model_fit_file = paste0("models/env/", file_base, format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  saveRDS(m_zln, file = model_fit_file)
  
  cat(paste0("Done. Model fit saved to:\n\n\t", model_fit_file, "\n\n"))
}

# fit model
# e_A_B_C_D where
## A is tmax_lag {10, 30}
## B is precip lag {10, 30}
## C is kg_code {1, 2, 3}
## D is boolean to include kg_code as RE in sigma {0, 1}
# 2 * 2 * 3 * 2 = 24 models

# build models
build_model_formula = function(a, b, c, d) {
  
  temp_var = case_when(
    a == 10 ~ temp_c_10_z,
    a == 30 ~ temp_c_30_z,
    .default = NA
  )
  
  precip_var = case_when(
    b == 10 ~ precip_mm_10_z,
    b == 30 ~ precip_mm_30_z,
    .default = NA
  )
  
  kg_var_1 = case_when(
    c == 1 ~ kg_code_1,
    c == 2 ~ kg_code_2,
    c == 3 ~ kg_code_3,
    .default = NA
  )
  
  kg_var_2 = case_when(
    c == 0 ~ NULL,
    c == 1 ~ kg_var_1,
    .default = NA
  )
  
  f_xx = bf(
    bd_load_co_density ~ 1 + 
      life_stage_simple * taxon_capture +
      (precip_var  + temp_var + I(temp_var^2)) +
      (0 + precip_var + temp_var + I(temp_var^2) | taxon_capture) +
      (1 | kg_var_1),
    
    hu ~ 1 + 
      life_stage_simple * taxon_capture +
      (precip_var + temp_var + I(temp_var^2)) +
      (0 + precip_var + temp_var + I(temp_var^2) | taxon_capture) +
      (1 | kg_var_1),
    
    sigma ~ 1 + taxon_capture + (1 | kg_var_2)
  )
  
  file_base = paste0("unified_model_results_zln_env_", a, "_", b, "_", c, "_", d, "_")
  
  return(liste(fxx, file_base))
}

# purr build all model formulas
build_model_formula()

# purr run each model sequentially
run_zln_model()
