## ============================================================
## 08_model_comparison.R
## Fit all four Bayesian model specifications on both campaigns,
## compare via k-fold ELPD and cross-month RMSPE / CRPS.
## Uses adaptive binning (build_adaptive_bins) — a standalone
## alternative to 01_run_analysis.R.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).
## ============================================================

rm(list = ls())

library(dplyr)
library(rstan)
library(loo)
library(scoringRules)
library(here)

source(here::here("R", "00_data_config.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

STAN_DIR  <- here::here("Stan")
DATA_DIR  <- here::here("Data")
SEED_BASE <- 1234L
check_data(DATA_DIR)

## ============================================================
## 1. Helper functions
## ============================================================

build_adaptive_bins <- function(data,
                                depth_col = "profondita",
                                id_col    = "ID",
                                min_depth = 0.1,
                                max_depth = 7.7,
                                base_step = 0.01,
                                min_sonde = 3) {
  df <- data %>%
    filter(!is.na(.data[[depth_col]]),
           .data[[depth_col]] >= min_depth,
           .data[[depth_col]] <= max_depth) %>%
    mutate(base_bin = floor(.data[[depth_col]] / base_step) * base_step)

  base_counts <- df %>%
    group_by(base_bin) %>%
    summarise(
      sonde   = list(sort(unique(.data[[id_col]]))),
      n_sonde = n_distinct(.data[[id_col]]),
      .groups = "drop"
    ) %>%
    arrange(base_bin)

  bins <- list(); i <- 1L; k <- 1L
  while (i <= nrow(base_counts)) {
    left_edge     <- base_counts$base_bin[i]
    current_sonde <- base_counts$sonde[[i]]
    j <- i
    while (length(current_sonde) < min_sonde && j < nrow(base_counts)) {
      j <- j + 1L
      current_sonde <- union(current_sonde, base_counts$sonde[[j]])
    }
    if (length(current_sonde) >= min_sonde) {
      right_edge <- base_counts$base_bin[j] + base_step
      bins[[k]] <- data.frame(
        bin_id = k, left = left_edge, right = right_edge,
        mid    = (left_edge + right_edge) / 2,
        width  = right_edge - left_edge,
        n_sonde = length(current_sonde)
      )
      k <- k + 1L; i <- j + 1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col = "profondita") {
  out <- data
  out$bin_id <- out$depth_left <- out$depth_right <- out$depth_mid <- NA_real_
  for (r in seq_len(nrow(bins))) {
    idx <- which(out[[depth_col]] >= bins$left[r] & out[[depth_col]] < bins$right[r])
    out$bin_id[idx]     <- bins$bin_id[r]
    out$depth_left[idx] <- bins$left[r]
    out$depth_right[idx]<- bins$right[r]
    out$depth_mid[idx]  <- bins$mid[r]
  }
  out
}

safe_se <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

## ============================================================
## 2. Load and pre-process both campaigns
## ============================================================

load_campaign <- function(file) {
  dat <- rio::import(file)
  df  <- dat %>%
    rename(
      profondita    = `Profondità (m)`,
      temperatura   = `temperatura (C°)`,
      salinita      = `Salinità (‰)`,
      ossigeno_mgkg = `Oxygen (mg/Kg)`,
      ph            = `pH`,
      trx_chla      = `Trx-chl(a)`,
      phycocyanin   = `Phycocyanin`,
      phycoerythrin = `Phycoerythrin`,
      dist          = `Distanza dalla foce (km)`
    ) %>%
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10)
  df$dist <- df$dist + abs(min(df$dist, na.rm = TRUE))
  df
}

df_jan <- load_campaign(file.path(DATA_DIR,
  "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx"))

df_jul <- load_campaign(file.path(DATA_DIR,
  "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx"))

vars_depth <- c("temperatura", "salinita", "ossigeno_mgkg",
                "ph", "trx_chla", "phycocyanin", "phycoerythrin", "dist")

aggregate_campaign <- function(df) {
  bins   <- build_adaptive_bins(df, min_depth = 0.01, max_depth = 10,
                                 base_step = 0.01, min_sonde = 3)
  df_bin <- assign_bins(df, bins)

  by_probe <- df_bin %>%
    filter(!is.na(bin_id)) %>%
    group_by(ID, bin_id, depth_left, depth_right, depth_mid) %>%
    summarise(across(all_of(vars_depth), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  by_probe %>%
    group_by(bin_id, depth_left, depth_right, depth_mid) %>%
    summarise(
      n_sonde = n(),
      across(all_of(vars_depth),
             list(mean = ~ mean(.x, na.rm = TRUE), se = ~ safe_se(.x)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(depth_left)
}

agg_jan <- aggregate_campaign(df_jan)
agg_jul <- aggregate_campaign(df_jul)

## ============================================================
## 3. Build Stan data
## ============================================================

build_stan_data <- function(agg) {
  cov_mean <- paste0(setdiff(vars_depth, "salinita"), "_mean")
  cov_sd   <- paste0(setdiff(vars_depth, "salinita"), "_se")

  dat <- agg %>%
    select(depth_mid, salinita_mean, salinita_se, all_of(cov_mean), all_of(cov_sd)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  y_raw <- dat$salinita_mean
  sy    <- sd(y_raw); my <- mean(y_raw)
  y_obs <- (y_raw - my) / sy
  tau_y <- pmax(dat$salinita_se / sy, 1e-6)

  X_raw   <- as.matrix(dat[, cov_mean])
  TauX    <- as.matrix(dat[, cov_sd])
  mx      <- colMeans(X_raw); sx <- apply(X_raw, 2, sd)
  X_obs   <- sweep(sweep(X_raw, 2, mx, "-"), 2, sx, "/")
  tau_x   <- pmax(sweep(TauX, 2, sx, "/"), 1e-6)

  d_raw <- dat$depth_mid
  d     <- (d_raw - min(d_raw)) / (max(d_raw) - min(d_raw))

  list(
    stan   = list(S = nrow(X_obs), K = ncol(X_obs),
                  y_obs = as.vector(y_obs), X_obs = X_obs,
                  tau_y = as.vector(tau_y), tau_x = tau_x,
                  d = as.vector(d)),
    meta   = list(my = my, sy = sy, mx = mx, sx = sx,
                  d_raw = d_raw, cov_names = cov_mean, y_raw = y_raw)
  )
}

sd_jan <- build_stan_data(agg_jan)
sd_jul <- build_stan_data(agg_jul)

## ============================================================
## 4. Stan model files
## ============================================================

model_files <- c(
  full    = file.path(STAN_DIR, "eiv_gp.stan"),
  noEIV   = file.path(STAN_DIR, "no_eiv_gp.stan"),
  noGP    = file.path(STAN_DIR, "eiv_no_gp.stan"),
  base    = file.path(STAN_DIR, "no_eiv_no_gp.stan")
)

strip_tau_x <- function(sd) {
  sd$tau_x <- NULL
  sd
}

get_stan_data <- function(model_name, sd) {
  if (model_name %in% c("noEIV", "base")) strip_tau_x(sd$stan) else sd$stan
}

## ============================================================
## 5. Fit all models × both campaigns
## ============================================================

stan_control <- list(adapt_delta = 0.95, max_treedepth = 12)

fits <- list()
for (camp in c("jan", "jul")) {
  sd_obj <- if (camp == "jan") sd_jan else sd_jul
  for (mod in names(model_files)) {
    key  <- paste0(mod, "_", camp)
    message("Fitting ", key, " ...")
    fits[[key]] <- stan(
      file    = model_files[mod],
      data    = get_stan_data(mod, sd_obj),
      chains  = 4L, iter = 2000L, warmup = 1000L, thin = 1L,
      seed    = SEED_BASE,
      control = stan_control
    )
  }
}

## ============================================================
## 6. MCMC diagnostics
## ============================================================

for (key in names(fits)) {
  s <- summary(fits[[key]])$summary
  rhat_max <- max(s[, "Rhat"], na.rm = TRUE)
  neff_min <- min(s[, "n_eff"], na.rm = TRUE)
  divs     <- get_num_divergent(fits[[key]])
  message(sprintf("%-20s  Rhat_max=%.3f  neff_min=%.0f  divergences=%d",
              key, rhat_max, neff_min, divs))
}

## ============================================================
## 7. Within-campaign k-fold ELPD
## ============================================================

kfold_elpd <- function(fit, stan_data_obj, model_name, k = 10L, seed = 42L) {
  S    <- stan_data_obj$stan$S
  flds <- kfold_split_stratified(k, z = seq_len(S))

  elpd_vec <- numeric(k)
  for (fold in seq_len(k)) {
    test_idx  <- which(flds == fold)
    train_idx <- setdiff(seq_len(S), test_idx)

    make_sub <- function(idx, sd) {
      sd_sub        <- sd
      sd_sub$S      <- length(idx)
      sd_sub$y_obs  <- sd$y_obs[idx]
      sd_sub$X_obs  <- sd$X_obs[idx, , drop = FALSE]
      sd_sub$tau_y  <- sd$tau_y[idx]
      if (!is.null(sd$tau_x)) sd_sub$tau_x <- sd$tau_x[idx, , drop = FALSE]
      sd_sub$d      <- sd$d[idx]
      sd_sub
    }

    sd_train <- make_sub(train_idx, stan_data_obj$stan)
    sd_test  <- make_sub(test_idx,  stan_data_obj$stan)

    if (model_name %in% c("noEIV", "base")) {
      sd_train <- strip_tau_x(list(stan = sd_train))$stan
      sd_test  <- strip_tau_x(list(stan = sd_test))$stan
    }

    fit_fold <- suppressMessages(stan(
      file    = model_files[model_name],
      data    = sd_train,
      chains  = 2L, iter = 1500L, warmup = 750L,
      seed    = seed + fold,
      refresh = 0L,
      control = list(adapt_delta = 0.92)
    ))

    post <- rstan::extract(fit_fold)

    S_test <- length(test_idx)
    log_pd <- numeric(nrow(post$beta))

    for (ii in seq_along(log_pd)) {
      if (model_name %in% c("full", "noGP")) {
        X_t <- post$X_true[ii, seq_len(S_test), , drop = FALSE]
        mu  <- as.numeric(X_t[1, , ] %*% post$beta[ii, ])
      } else {
        mu <- as.numeric(sd_test$X_obs %*% post$beta[ii, ])
      }

      if (model_name %in% c("full", "noEIV")) {
        su  <- post$sigma_u[ii]; el <- post$ell[ii]
        Ktt <- su^2 * exp(-abs(outer(sd_test$d, sd_test$d, "-")) / el)
        Sy  <- Ktt + diag(sd_test$tau_y^2 + 1e-8)
        log_pd[ii] <- mvtnorm::dmvnorm(sd_test$y_obs, mu, Sy, log = TRUE)
      } else {
        log_pd[ii] <- sum(dnorm(sd_test$y_obs, mu, sd_test$tau_y, log = TRUE))
      }
    }
    elpd_vec[fold] <- log(mean(exp(log_pd - max(log_pd)))) + max(log_pd)
  }
  sum(elpd_vec)
}

elpd_results <- list()
for (camp in c("jan", "jul")) {
  sd_obj <- if (camp == "jan") sd_jan else sd_jul
  for (mod in names(model_files)) {
    key <- paste0(mod, "_", camp)
    message("k-fold ELPD for ", key)
    elpd_results[[key]] <- kfold_elpd(fits[[key]], sd_obj, mod)
  }
}

elpd_df <- data.frame(
  model    = names(elpd_results),
  elpd_cv  = unlist(elpd_results)
) %>%
  tidyr::separate(model, into = c("spec", "campaign"), sep = "_")

print(elpd_df)

## ============================================================
## 8. Cross-month posterior predictive RMSPE and CRPS
## ============================================================

cross_predict <- function(fit, sd_train_obj, sd_new_obj, model_name,
                          n_draws = 500L) {

  post   <- rstan::extract(fit)
  n_post <- length(post$beta[, 1])
  idx    <- sample(n_post, min(n_draws, n_post))

  sd_new  <- sd_new_obj$stan
  meta_tr <- sd_train_obj$meta

  S_new <- sd_new$S
  y_pred_std <- matrix(NA_real_, nrow = length(idx), ncol = S_new)

  for (ii in seq_along(idx)) {
    it <- idx[ii]

    if (model_name %in% c("full", "noGP")) {
      X_star_new   <- matrix(NA_real_, S_new, sd_new$K)
      for (k in seq_len(sd_new$K)) {
        X_star_new[, k] <- rnorm(S_new, sd_new$X_obs[, k], sd_new$tau_x[, k])
      }
    } else {
      X_star_new <- sd_new$X_obs
    }

    beta_it <- post$beta[it, ]
    mu_new  <- as.numeric(X_star_new %*% beta_it)

    if (model_name %in% c("full", "noEIV")) {
      su  <- post$sigma_u[it]; el <- post$ell[it]
      d_tr   <- sd_train_obj$stan$d
      y_tr   <- sd_train_obj$stan$y_obs
      tau_tr <- sd_train_obj$stan$tau_y

      if (model_name == "full") {
        Xt_tr <- post$X_true[it, , ]
        mu_tr <- as.numeric(Xt_tr %*% beta_it)
      } else {
        mu_tr <- as.numeric(sd_train_obj$stan$X_obs %*% beta_it)
      }

      K_tr  <- su^2 * exp(-abs(outer(d_tr, d_tr, "-")) / el) +
               diag(tau_tr^2 + 1e-8)
      K_new_tr <- su^2 * exp(-abs(outer(sd_new$d, d_tr, "-")) / el)

      v           <- solve(K_tr, y_tr - mu_tr)
      mu_pred     <- mu_new + as.numeric(K_new_tr %*% v)
      K_pred      <- su^2 * diag(S_new) -
                     K_new_tr %*% solve(K_tr, t(K_new_tr))
      K_pred      <- K_pred + diag(sd_new$tau_y^2 + 1e-8)
      diag(K_pred)<- pmax(diag(K_pred), 1e-8)

      y_pred_std[ii, ] <- MASS::mvrnorm(1, mu_pred, K_pred)
    } else {
      y_pred_std[ii, ] <- rnorm(S_new, mu_new, sd_new$tau_y)
    }
  }

  y_pred_std * meta_tr$sy + meta_tr$my
}

eval_preds <- function(y_pred_mat, y_obs_raw) {
  y_mean <- colMeans(y_pred_mat)
  rmspe  <- sqrt(mean((y_mean - y_obs_raw)^2))
  crps_v <- scoringRules::crps_sample(y = y_obs_raw,
                                      dat = t(y_pred_mat))
  list(rmspe = rmspe, crps = mean(crps_v))
}

cross_results <- list()
for (mod in names(model_files)) {
  preds_jul2jan <- cross_predict(fits[[paste0(mod, "_jul")]],
                                  sd_jul, sd_jan, mod)
  r1 <- eval_preds(preds_jul2jan, sd_jan$meta$y_raw)

  preds_jan2jul <- cross_predict(fits[[paste0(mod, "_jan")]],
                                  sd_jan, sd_jul, mod)
  r2 <- eval_preds(preds_jan2jul, sd_jul$meta$y_raw)

  message(sprintf("%s  Jul->Jan  RMSPE=%.3f  CRPS=%.3f", mod, r1$rmspe, r1$crps))
  message(sprintf("%s  Jan->Jul  RMSPE=%.3f  CRPS=%.3f", mod, r2$rmspe, r2$crps))
  cross_results[[mod]] <- list(jul2jan = r1, jan2jul = r2)
}

## ============================================================
## 9. Summary tables
## ============================================================

post_full_jan <- rstan::extract(fits[["full_jan"]])
cov_names     <- sd_jan$meta$cov_names

beta_summary <- t(apply(post_full_jan$beta, 2, function(b) {
  c(median = median(b),
    lo95   = quantile(b, 0.025),
    hi95   = quantile(b, 0.975),
    P_pos  = mean(b > 0))
}))
rownames(beta_summary) <- cov_names
print(round(beta_summary, 3))

## ============================================================
## 10. Save all fitted objects
## ============================================================

save(fits, elpd_df, cross_results, beta_summary,
     sd_jan, sd_jul, agg_jan, agg_jul,
     file = here::here("Results", "results_model_comparison.RData"))

message("Done. Results saved to Results/results_model_comparison.RData")
