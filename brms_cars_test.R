# ============================================================
# brms test script: Bayesian regression on the 'cars' dataset
# cars: 50 observations of speed (mph) and stopping distance (ft)
# Model: dist ~ speed (simple linear regression, Bayesian)
# ============================================================

# load packages
library(brms)
library(cmdstanr)

cmdstanr::check_cmdstan_toolchain()   # confirms Stan is found

cat("\nRunning test cars BRMS model:\n\n")

# create output folder
if (!dir.exists("models")) dir.create("models")

# fit model
fit_01 <- brm(
  formula = dist ~ speed,
  data    = cars,
  family  = gaussian(),
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 125,
  cores   = 4,
  backend = "cmdstanr"
)

# save model
model_fit_file = paste0("models/fit_01_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(fit_01, file = model_fit_file)

cat(paste0("Done. Model fit saved to:\n\n\t", getwd(), "/", model_fit_file, "\n\n"))
