## 03_refit_save_samples.R
##
## Refit M_full on Jan-26 and Jul-25 and save the posterior quantities
## required for Fig 3 (PPC) and Fig 4 (GP residual u_s).
##
## Output: Results/mfull_posterior_samples.rds with elements jan, jul.
## Each is a list with depth (bin midpoints in m), y_obs (back-converted
## salinity), y_rep (draws x S, in original salinity units), and u
## (draws x S, posterior of the GP residual on the standardised scale).
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(mvtnorm); library(here)
})
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

SEED <- 1234L
source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jan(), load_jul(), build_bins(), aggregate_to_bins(), vars_phys

build_data <- function(loader) {
  df <- loader()
  bins <- build_bins(df); agg <- aggregate_to_bins(df, bins)
  cov_mean <- paste0(vars_phys, "_mean"); cov_se <- paste0(vars_phys, "_se")
  dat <- agg %>% select(depth_mid, salinita_mean, salinita_se,
                        all_of(cov_mean), all_of(cov_se)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  y_raw <- dat$salinita_mean; sy<-sd(y_raw); my<-mean(y_raw)
  y_obs <- (y_raw-my)/sy
  tau_y <- pmax(dat$salinita_se/sy, 1e-6)
  X_raw <- as.matrix(dat[,cov_mean]); TauX<-as.matrix(dat[,cov_se])
  mx <- colMeans(X_raw); sx <- apply(X_raw,2,sd)
  X_obs <- sweep(sweep(X_raw,2,mx,"-"),2,sx,"/")
  tau_x <- pmax(sweep(TauX,2,sx,"/"),1e-6)
  d_raw <- dat$depth_mid; d <- (d_raw-min(d_raw))/(max(d_raw)-min(d_raw))
  list(stan=list(S=nrow(X_obs), K=ncol(X_obs),
                 y_obs=as.vector(y_obs), X_obs=X_obs,
                 tau_y=as.vector(tau_y), tau_x=tau_x, d=as.vector(d)),
       meta=list(my=my, sy=sy, depth=d_raw, y_raw=y_raw, tau_y=as.vector(tau_y)))
}

sample_u <- function(post, sd_obj) {
  S <- sd_obj$stan$S; tau <- sd_obj$stan$tau_y; d <- sd_obj$stan$d
  ndraw <- length(post$sigma_u)
  U <- matrix(NA_real_, ndraw, S)
  Iy <- diag(tau^2)
  for (i in seq_len(ndraw)) {
    su <- post$sigma_u[i]; ell <- post$ell[i]
    K <- outer(d,d,function(a,b) su^2 * exp(-abs(a-b)/ell))
    Xstar <- if ("X_true" %in% names(post)) post$X_true[i,,]
             else sd_obj$stan$X_obs
    mu <- as.vector(Xstar %*% post$beta[i,])
    A  <- K %*% solve(K + Iy)
    m_u <- A %*% (sd_obj$stan$y_obs - mu)
    V_u <- K - A %*% K
    V_u <- (V_u + t(V_u))/2 + diag(1e-8, S)
    U[i,] <- mvtnorm::rmvnorm(1, mean=as.vector(m_u), sigma=V_u)
  }
  U
}

fit_and_collect <- function(loader, label) {
  sd_obj <- build_data(loader)
  sm <- stan_model(here::here("Stan","eiv_gp.stan"))
  fit <- sampling(sm, data=sd_obj$stan, chains=4, iter=2000, warmup=1000,
                  seed=SEED, control=list(adapt_delta=0.95, max_treedepth=12),
                  refresh=100)
  post <- rstan::extract(fit)
  y_rep <- post$y_rep * sd_obj$meta$sy + sd_obj$meta$my
  idx <- sample(seq_along(post$sigma_u), min(200, length(post$sigma_u)))
  post_sub <- list(beta=post$beta[idx,,drop=FALSE],
                   sigma_u=post$sigma_u[idx], ell=post$ell[idx],
                   X_true=if ("X_true" %in% names(post)) post$X_true[idx,,,drop=FALSE] else NULL)
  u <- sample_u(post_sub, sd_obj)
  list(depth=sd_obj$meta$depth, y_obs=sd_obj$meta$y_raw,
       y_rep=y_rep[idx,], u=u)
}

samp <- list(jan = fit_and_collect(load_jan, "Jan-26"),
             jul = fit_and_collect(load_jul, "Jul-25"))
saveRDS(samp, here::here("Results", "mfull_posterior_samples.rds"))
