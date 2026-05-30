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
if (!dir.exists("models/s2")) dir.create("models/s2")

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
  model_fit_file = paste0("models/s2/", file_base, format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  saveRDS(m_zln, file = model_fit_file)
  
  cat(paste0("Done. Model fit saved to:\n\n\t", model_fit_file, "\n\n"))
}

# fit model
# isolate scales between taxa
f_11e = bf(
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
  
  sigma ~ 1 + taxon_capture + (1 | taxon_capture:population) + (1 | kg_code_3)
)

base_11e = "unified_model_results_zln_life_tax_pt2_rgrtaxacorpt2poplsyrkg3_sigma_taxpopkg_p10t30_s02"

# variations

## A: random effect structure in mu and hu terms
### 1: none
### 2: (1 | taxon_capture:population:life_stage_simple)
### 3: (1 | gr(taxon_capture:population:life_stage_simple, by = taxon_capture))
### 4: ### 3: (1 | gr(taxon_capture:population:life_stage_simple, by = taxon_capture)) + (0 + life_stage_simple | gr(taxon_capture:population:year, by = taxon_capture))
### 5: ### 3: (1 | gr(taxon_capture:population:life_stage_simple, by = taxon_capture)) + (0 + life_stage_simple || gr(taxon_capture:population:year, by = taxon_capture))

## B: random effect strucutre in sigma term
### 1: none
### 2: (1 | population)
### 3: (1 | taxon_capture:population)
### 4: (1 | gr(taxon_capture:population), by = taxon_capture)
