## ============================================================
## sim/00_extract_truth.R
## Fit M_full on the two real campaigns and store posterior medians
## as the ground-truth parameter set for the simulation DGP.
## Also stores the empirical covariance of the standardised
## covariates, so the DGP reproduces the real |r|~0.99 collinearity
## (which the model's independent X* prior does NOT assume — this
## builds a realistic degree of prior misspecification into the DGP).
## ============================================================
rm(list = ls())
library(here)
SIM_DIR  <- here::here("Simulation")
STAN_DIR <- here::here("Stan")
source(here::here("Simulation", "data_prep.R"))

suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)
SEED <- 1234L

sd_all <- load_campaign_stan_data()

sm_full <- stan_model(file.path(STAN_DIR, "eiv_gp.stan"), model_name = "EIV_GP")

fit_full <- function(stan_data, seed_offset) {
  sampling(sm_full, data = stan_data, chains = 4L, iter = 2000L,
           warmup = 1000L, thin = 1L, seed = SEED + seed_offset,
           refresh = 0L, control = list(adapt_delta = 0.95,
                                        max_treedepth = 12))
}

extract_truth <- function(sd_obj, fit, label) {
  post <- rstan::extract(fit)
  K    <- ncol(post$beta)
  Xobs <- sd_obj$stan$X_obs                       # already standardised
  list(
    label   = label,
    S       = sd_obj$stan$S,
    K       = K,
    d       = sd_obj$stan$d,
    tau_y   = sd_obj$stan$tau_y,
    tau_x   = sd_obj$stan$tau_x,                   # S x K
    beta    = apply(post$beta,    2, median),
    mu_x    = apply(post$mu_x,    2, median),
    sigma_x = apply(post$sigma_x, 2, median),
    sigma_u = median(post$sigma_u),
    ell     = median(post$ell),
    # Empirical mean/covariance of the standardised observed covariates:
    # used by the DGP to inject the real cross-covariate collinearity.
    Xbar    = colMeans(Xobs),
    Sigma_X = cov(Xobs),
    cov_names = sd_obj$meta$cov_names
  )
}

fj <- fit_full(sd_all$jan$stan, 0L)
fl <- fit_full(sd_all$jul$stan, 1L)

truth <- list(
  jan = extract_truth(sd_all$jan, fj, "Jan-26"),
  jul = extract_truth(sd_all$jul, fl, "Jul-25")
)

saveRDS(truth, file.path(SIM_DIR, "truth_params.rds"))
for (nm in names(truth)) {
  tr <- truth[[nm]]
  print(round(setNames(tr$beta, tr$cov_names), 3))
}
