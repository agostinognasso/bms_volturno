## ============================================================
## 01_run_analysis.R  —  Full Bayesian model comparison
## Step 1 of the main analysis pipeline.
## Fits all 4 models × 2 campaigns, computes WAIC/LOO/beta
## posteriors and cross-month predictions. Saves results.RData.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).
## ============================================================
rm(list = ls())

library(dplyr)
library(readxl)
library(rstan)
library(loo)
library(scoringRules)
library(mvtnorm)
library(here)

source(here::here("R", "00_data_config.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

STAN_DIR <- here::here("Stan")
DATA_DIR <- here::here("Data")
SEED     <- 1234L
check_data(DATA_DIR)

## ============================================================
## 1. Helper functions
## ============================================================

safe_se <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

build_bins <- function(data, depth_col = "profondita", id_col = "ID",
                       min_depth = 0.01, max_depth = 10,
                       base_step = 0.01, min_sonde = 3) {
  df <- data %>%
    filter(!is.na(.data[[depth_col]]),
           .data[[depth_col]] >= min_depth,
           .data[[depth_col]] <= max_depth) %>%
    mutate(base_bin = floor(.data[[depth_col]] / base_step) * base_step)

  base_counts <- df %>%
    group_by(base_bin) %>%
    summarise(sonde   = list(sort(unique(.data[[id_col]]))),
              n_sonde = n_distinct(.data[[id_col]]), .groups = "drop") %>%
    arrange(base_bin)

  bins <- list(); i <- 1L; k <- 1L
  while (i <= nrow(base_counts)) {
    left_edge <- base_counts$base_bin[i]
    cur       <- base_counts$sonde[[i]]
    j         <- i
    while (length(cur) < min_sonde && j < nrow(base_counts)) {
      j <- j + 1L; cur <- union(cur, base_counts$sonde[[j]])
    }
    if (length(cur) >= min_sonde) {
      right_edge <- base_counts$base_bin[j] + base_step
      bins[[k]]  <- data.frame(bin_id = k, left = left_edge,
                                right = right_edge,
                                mid   = (left_edge + right_edge) / 2,
                                width = right_edge - left_edge,
                                n_sonde = length(cur))
      k <- k + 1L; i <- j + 1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col = "profondita") {
  out <- data
  out$bin_id <- out$depth_mid <- NA_real_
  for (r in seq_len(nrow(bins))) {
    idx <- which(out[[depth_col]] >= bins$left[r] & out[[depth_col]] < bins$right[r])
    out$bin_id[idx]    <- bins$bin_id[r]
    out$depth_mid[idx] <- bins$mid[r]
  }
  out
}

## ============================================================
## 2. Load and clean both campaigns
## ============================================================

load_jan <- function() {
  f  <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx")
  d  <- read_excel(f)
  df <- d %>%
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (C°)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph            = `pH`,
           trx_chla      = `Trx-chl(a)`,
           phycocyanin   = `Phycocyanin`,
           phycoerythrin = `Phycoerythrin`,
           dist          = `Distanza dalla foce (km)`) %>%
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10)
  df$dist <- df$dist - min(df$dist, na.rm = TRUE)
  df
}

load_jul <- function() {
  f  <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx")
  d  <- read_excel(f)
  df <- d %>%
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (°C)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph            = `pH`,
           trx_chla      = `Trx-chl(a)`,
           phycocyanin   = `Phycocyanin`,
           phycoerythrin = `Phycoerythrin`,
           dist          = `Distanza dalla Foce(km)`) %>%
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10)
  df$dist <- df$dist - min(df$dist, na.rm = TRUE)
  df
}

df_jan <- load_jan()
df_jul <- load_jul()


## ============================================================
## 3. Aggregate to depth bins
## ============================================================

vars_phys <- c("temperatura","salinita","ossigeno_mgkg","ph",
               "trx_chla","phycocyanin","phycoerythrin","dist")

aggregate_to_bins <- function(df) {
  bins   <- build_bins(df)
  df_bin <- assign_bins(df, bins)

  by_probe <- df_bin %>%
    filter(!is.na(bin_id)) %>%
    group_by(ID, bin_id, depth_mid) %>%
    summarise(across(all_of(vars_phys), ~ mean(.x, na.rm=TRUE)), .groups="drop")

  by_probe %>%
    group_by(bin_id, depth_mid) %>%
    summarise(
      n_sonde = n(),
      across(all_of(vars_phys),
             list(mean=~mean(.x, na.rm=TRUE), se=~safe_se(.x)),
             .names="{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(depth_mid)
}

agg_jan <- aggregate_to_bins(df_jan)
agg_jul <- aggregate_to_bins(df_jul)


## ============================================================
## 4. Build Stan data
## ============================================================

cov_vars <- setdiff(vars_phys, "salinita")   # 7 covariates

build_stan_data <- function(agg) {
  cov_mean <- paste0(cov_vars, "_mean")
  cov_se   <- paste0(cov_vars, "_se")

  dat <- agg %>%
    select(depth_mid, salinita_mean, salinita_se,
           all_of(cov_mean), all_of(cov_se)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  y_raw <- dat$salinita_mean
  sy    <- sd(y_raw); my <- mean(y_raw)
  y_obs <- (y_raw - my) / sy
  tau_y <- pmax(dat$salinita_se / sy, 1e-6)

  X_raw  <- as.matrix(dat[, cov_mean])
  TauX   <- as.matrix(dat[, cov_se])
  mx     <- colMeans(X_raw); sx <- apply(X_raw, 2, sd)
  X_obs  <- sweep(sweep(X_raw, 2, mx, "-"), 2, sx, "/")
  tau_x  <- pmax(sweep(TauX, 2, sx, "/"), 1e-6)

  d_raw  <- dat$depth_mid
  d      <- (d_raw - min(d_raw)) / (max(d_raw) - min(d_raw))

  list(
    stan = list(S=nrow(X_obs), K=ncol(X_obs),
                y_obs=as.vector(y_obs), X_obs=X_obs,
                tau_y=as.vector(tau_y), tau_x=tau_x,
                d=as.vector(d)),
    meta = list(my=my, sy=sy, mx=mx, sx=sx, d_raw=d_raw,
                y_raw=y_raw, cov_names=cov_vars,
                n_bins=nrow(dat))
  )
}

sd_jan <- build_stan_data(agg_jan)
sd_jul <- build_stan_data(agg_jul)


## ============================================================
## 5. Compile Stan models (once)
## ============================================================

sm_full  <- stan_model(file.path(STAN_DIR, "eiv_gp.stan"),      model_name="EIV_GP")
sm_noEIV <- stan_model(file.path(STAN_DIR, "no_eiv_gp.stan"),   model_name="noEIV_GP")
sm_noGP  <- stan_model(file.path(STAN_DIR, "eiv_no_gp.stan"),   model_name="EIV_noGP")
sm_base  <- stan_model(file.path(STAN_DIR, "no_eiv_no_gp.stan"),model_name="base")

## ============================================================
## 6. Fit all models on both campaigns
## ============================================================

stan_ctrl <- list(adapt_delta=0.95, max_treedepth=12)

fit_model <- function(sm, stan_data, seed_offset=0) {
  sampling(sm, data=stan_data, chains=4L, iter=2000L, warmup=1000L,
           thin=1L, seed=SEED+seed_offset, refresh=200L,
           control=stan_ctrl)
}

no_eiv_data <- function(sd) {
  d <- sd$stan; d$tau_x <- NULL; d
}

fit_full_jan  <- fit_model(sm_full,  sd_jan$stan, 0)
fit_full_jul  <- fit_model(sm_full,  sd_jul$stan, 1)
fit_noEIV_jan <- fit_model(sm_noEIV, no_eiv_data(sd_jan), 2)
fit_noEIV_jul <- fit_model(sm_noEIV, no_eiv_data(sd_jul), 3)
fit_noGP_jan  <- fit_model(sm_noGP,  sd_jan$stan, 4)
fit_noGP_jul  <- fit_model(sm_noGP,  sd_jul$stan, 5)
fit_base_jan  <- fit_model(sm_base,  no_eiv_data(sd_jan), 6)
fit_base_jul  <- fit_model(sm_base,  no_eiv_data(sd_jul), 7)

## ============================================================
## 7. MCMC diagnostics
## ============================================================

diag_summary <- function(fit, name) {
  s        <- summary(fit)$summary
  rhat_max <- max(s[,"Rhat"], na.rm=TRUE)
  neff_min <- min(s[,"n_eff"], na.rm=TRUE)
  divs     <- get_num_divergent(fit)
  bfmi     <- min(get_bfmi(fit))
  data.frame(model=name, rhat_max=rhat_max, neff_min=neff_min,
             divergences=divs, bfmi=bfmi)
}

diag_df <- bind_rows(
  diag_summary(fit_full_jan,  "full_jan"),
  diag_summary(fit_full_jul,  "full_jul"),
  diag_summary(fit_noEIV_jan, "noEIV_jan"),
  diag_summary(fit_noEIV_jul, "noEIV_jul"),
  diag_summary(fit_noGP_jan,  "noGP_jan"),
  diag_summary(fit_noGP_jul,  "noGP_jul"),
  diag_summary(fit_base_jan,  "base_jan"),
  diag_summary(fit_base_jul,  "base_jul")
)

## ============================================================
## 8. LOO-CV / WAIC
## ============================================================

waic_summary <- function(fit, name) {
  ll   <- extract_log_lik(fit, parameter_name="log_lik")
  w    <- waic(ll)
  data.frame(model=name,
             waic     = w$estimates["waic","Estimate"],
             waic_se  = w$estimates["waic","SE"],
             elpd_waic= w$estimates["elpd_waic","Estimate"])
}

waic_df <- bind_rows(
  waic_summary(fit_full_jan,  "full_jan"),
  waic_summary(fit_noEIV_jan, "noEIV_jan"),
  waic_summary(fit_noGP_jan,  "noGP_jan"),
  waic_summary(fit_base_jan,  "base_jan"),
  waic_summary(fit_full_jul,  "full_jul"),
  waic_summary(fit_noEIV_jul, "noEIV_jul"),
  waic_summary(fit_noGP_jul,  "noGP_jul"),
  waic_summary(fit_base_jul,  "base_jul")
)

loo_noGP_jan  <- loo(extract_log_lik(fit_noGP_jan))
loo_base_jan  <- loo(extract_log_lik(fit_base_jan))
loo_noGP_jul  <- loo(extract_log_lik(fit_noGP_jul))
loo_base_jul  <- loo(extract_log_lik(fit_base_jul))

    "se", round(loo_noGP_jan$estimates["elpd_loo","SE"],1), "\n")
    "se", round(loo_base_jan$estimates["elpd_loo","SE"],1), "\n")

## ============================================================
## 9. Beta posteriors
## ============================================================

beta_posterior <- function(fit, cov_names, model_name) {
  post   <- rstan::extract(fit)
  K      <- ncol(post$beta)
  out    <- lapply(seq_len(K), function(k) {
    b <- post$beta[, k]
    data.frame(
      model    = model_name,
      covariate= cov_names[k],
      median   = median(b),
      lo95     = quantile(b, 0.025),
      hi95     = quantile(b, 0.975),
      P_pos    = mean(b > 0)
    )
  })
  bind_rows(out)
}

cov_labels <- c("Temperature","O2","pH","Chl-a","Phycocyanin","Phycoerythrin","Dist")

beta_jan <- bind_rows(
  beta_posterior(fit_full_jan,  cov_labels, "M_full"),
  beta_posterior(fit_noEIV_jan, cov_labels, "M_noEIV"),
  beta_posterior(fit_noGP_jan,  cov_labels, "M_noGP"),
  beta_posterior(fit_base_jan,  cov_labels, "M_base")
)
print(beta_jan, digits=3)

beta_jul <- bind_rows(
  beta_posterior(fit_full_jul,  cov_labels, "M_full"),
  beta_posterior(fit_noEIV_jul, cov_labels, "M_noEIV"),
  beta_posterior(fit_noGP_jul,  cov_labels, "M_noGP"),
  beta_posterior(fit_base_jul,  cov_labels, "M_base")
)
print(beta_jul, digits=3)

## ============================================================
## 10. GP hyperparameters (M_full)
## ============================================================

gp_summary <- function(fit, model_name) {
  post <- rstan::extract(fit)
  data.frame(
    model    = model_name,
    su_med   = median(post$sigma_u),
    su_lo    = quantile(post$sigma_u, 0.025),
    su_hi    = quantile(post$sigma_u, 0.975),
    ell_med  = median(post$ell),
    ell_lo   = quantile(post$ell, 0.025),
    ell_hi   = quantile(post$ell, 0.975)
  )
}

gp_df <- bind_rows(
  gp_summary(fit_full_jan, "full_jan"),
  gp_summary(fit_full_jul, "full_jul")
)
print(gp_df, digits=3)

## ============================================================
## 11. Cross-month predictions (M_full)
## ============================================================

gp_predict <- function(fit_train, sd_train, sd_new, n_draws=400L) {
  post  <- rstan::extract(fit_train)
  n_all <- length(post$sigma_u)
  idx   <- sample(n_all, min(n_draws, n_all))

  S_tr  <- sd_train$stan$S
  S_new <- sd_new$stan$S

  y_pred_std <- matrix(NA_real_, length(idx), S_new)

  for (ii in seq_along(idx)) {
    it      <- idx[ii]
    su      <- post$sigma_u[it]
    el      <- post$ell[it]
    beta_it <- post$beta[it, ]
    Xt_tr   <- post$X_true[it, , ]

    mu_tr   <- as.numeric(Xt_tr %*% beta_it)

    X_star_new <- matrix(NA_real_, S_new, sd_new$stan$K)
    for (k in seq_len(sd_new$stan$K)) {
      mu_k  <- post$mu_x[it, k]
      sg_k  <- post$sigma_x[it, k]
      tau_k  <- sd_new$stan$tau_x[, k]
      v_prior <- sg_k^2
      v_meas  <- tau_k^2
      w_prior <- v_meas / (v_prior + v_meas)
      w_data  <- v_prior / (v_prior + v_meas)
      post_mu  <- w_prior * mu_k + w_data * sd_new$stan$X_obs[, k]
      post_var <- (v_prior * v_meas) / (v_prior + v_meas)
      X_star_new[, k] <- rnorm(S_new, post_mu, sqrt(post_var))
    }

    mu_new  <- as.numeric(X_star_new %*% beta_it)

    d_tr    <- sd_train$stan$d
    d_new   <- sd_new$stan$d
    tau_tr  <- sd_train$stan$tau_y
    tau_new <- sd_new$stan$tau_y
    y_tr    <- sd_train$stan$y_obs

    K_tr    <- su^2 * exp(-abs(outer(d_tr,  d_tr,  "-"))/el) +
               diag(tau_tr^2 + 1e-8)
    K_cross <- su^2 * exp(-abs(outer(d_new, d_tr,  "-"))/el)
    K_nn    <- su^2 * diag(S_new) +
               diag(tau_new^2 + 1e-8)

    L       <- chol(K_tr)
    alpha   <- backsolve(L, forwardsolve(t(L), y_tr - mu_tr))
    cond_mu <- mu_new + as.numeric(K_cross %*% alpha)
    cond_K  <- K_nn - K_cross %*% backsolve(L, forwardsolve(t(L), t(K_cross)))
    cond_K  <- (cond_K + t(cond_K)) / 2
    diag(cond_K) <- pmax(diag(cond_K), 1e-8)

    y_pred_std[ii, ] <- MASS::mvrnorm(1, cond_mu, cond_K)
  }

  y_pred_std * sd_train$meta$sy + sd_train$meta$my
}

pred_jul2jan <- gp_predict(fit_full_jul, sd_jul, sd_jan)
pred_jan2jul <- gp_predict(fit_full_jan, sd_jan, sd_jul)

eval_pred <- function(pred_mat, y_obs_raw, label) {
  y_mean <- colMeans(pred_mat)
  rmspe  <- sqrt(mean((y_mean - y_obs_raw)^2))
  crps_v <- scoringRules::crps_sample(y=y_obs_raw, dat=t(pred_mat))
  data.frame(direction=label, rmspe=rmspe, crps=mean(crps_v))
}

cross_df <- bind_rows(
  eval_pred(pred_jul2jan, sd_jan$meta$y_raw, "Jul->Jan"),
  eval_pred(pred_jan2jul, sd_jul$meta$y_raw, "Jan->Jul")
)

## ============================================================
## 12. Descriptive statistics
## ============================================================

desc_jan <- df_jan %>% group_by(ID) %>%
  summarise(n=n(), sal_mean=mean(salinita,na.rm=T),
            temp_mean=mean(temperatura,na.rm=T),
            o2_mean=mean(ossigeno_mgkg,na.rm=T),
            depth_max=max(profondita,na.rm=T),
            dist=first(dist), .groups="drop") %>% arrange(dist)

desc_jul <- df_jul %>% group_by(ID) %>%
  summarise(n=n(), sal_mean=mean(salinita,na.rm=T),
            temp_mean=mean(temperatura,na.rm=T),
            o2_mean=mean(ossigeno_mgkg,na.rm=T),
            depth_max=max(profondita,na.rm=T),
            dist=first(dist), .groups="drop") %>% arrange(dist)


stats_jan <- df_jan %>% summarise(
  n=n(), sal_mean=mean(salinita,na.rm=T), sal_sd=sd(salinita,na.rm=T),
  sal_min=min(salinita,na.rm=T), sal_max=max(salinita,na.rm=T),
  temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T),
  depth_max=max(profondita,na.rm=T))

stats_jul <- df_jul %>% summarise(
  n=n(), sal_mean=mean(salinita,na.rm=T), sal_sd=sd(salinita,na.rm=T),
  sal_min=min(salinita,na.rm=T), sal_max=max(salinita,na.rm=T),
  temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T),
  depth_max=max(profondita,na.rm=T))


## ============================================================
## 13. Save everything
## ============================================================

save(
  fit_full_jan, fit_full_jul,
  fit_noEIV_jan, fit_noEIV_jul,
  fit_noGP_jan,  fit_noGP_jul,
  fit_base_jan,  fit_base_jul,
  sd_jan, sd_jul, agg_jan, agg_jul,
  df_jan, df_jul,
  diag_df, waic_df, beta_jan, beta_jul, gp_df,
  cross_df,
  desc_jan, desc_jul, stats_jan, stats_jul,
  file = here::here("Results", "results.RData")
)

