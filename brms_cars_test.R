# ============================================================
# brms test script: Bayesian regression on the 'cars' dataset
# cars: 50 observations of speed (mph) and stopping distance (ft)
# Model: dist ~ speed (simple linear regression, Bayesian)
# ============================================================

# restore packages from renv
renv::restore()

# create output folder
if (!dir.exists("models")) dir.create("models")

# fit model
fit_01 <- brm(
  formula = dist ~ speed,
  data    = cars,
  family  = gaussian(),
  chains  = 2,
  iter    = 2000,
  warmup  = 1000,
  seed    = 42,
  cores   = 2
)

# compute loo
fit_01 <- add_criterion(fit_01, "loo")

# save model
saveRDS(fit_01, file = "models/fit_01.rds")
