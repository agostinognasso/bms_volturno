## 02_fix_predictions.R — Cross-month predictions (robust GP version).
## Loads results.RData from Step 1 and recomputes cross-month predictions
## using a numerically robust GP solver (nearPD). Saves predictions.RData.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

rm(list = ls())

library(dplyr)
library(rstan)
library(scoringRules)
library(mvtnorm)
library(Matrix)
library(here)

load(here::here("Results", "results.RData"))

## ----------------------------------------------------------------
## Helper: robust GP conditional prediction with nearPD jitter
## ----------------------------------------------------------------
gp_predict_robust <- function(fit_train, sd_train, sd_new,
                               model_name, n_draws = 400L) {
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
    Xt_tr   <- if (!is.null(post$X_true)) post$X_true[it, , ] else NULL

    mu_tr <- if (!is.null(Xt_tr)) {
      as.numeric(Xt_tr %*% beta_it)
    } else {
      as.numeric(sd_train$stan$X_obs %*% beta_it)
    }

    if (!is.null(post$mu_x)) {
      mu_k_all  <- post$mu_x[it, ]
      sg_k_all  <- post$sigma_x[it, ]
      X_star_new <- matrix(NA_real_, S_new, sd_new$stan$K)
      for (k in seq_len(sd_new$stan$K)) {
        tau_k   <- sd_new$stan$tau_x[, k]
        v_pr    <- sg_k_all[k]^2
        v_me    <- tau_k^2
        post_mu  <- (v_me * mu_k_all[k] + v_pr * sd_new$stan$X_obs[, k]) / (v_pr + v_me)
        post_var <- v_pr * v_me / (v_pr + v_me)
        X_star_new[, k] <- rnorm(S_new, post_mu, sqrt(post_var))
      }
      mu_new <- as.numeric(X_star_new %*% beta_it)
    } else {
      mu_new <- as.numeric(sd_new$stan$X_obs %*% beta_it)
    }

    d_tr    <- sd_train$stan$d
    d_new   <- sd_new$stan$d
    tau_tr  <- sd_train$stan$tau_y
    tau_new <- sd_new$stan$tau_y
    y_tr    <- sd_train$stan$y_obs

    K_tr    <- su^2 * exp(-abs(outer(d_tr,  d_tr,  "-")) / el) +
               diag(tau_tr^2 + 1e-8)
    K_cross <- su^2 * exp(-abs(outer(d_new, d_tr,  "-")) / el)
    K_nn    <- su^2 * exp(-abs(outer(d_new, d_new, "-")) / el) +
               diag(tau_new^2 + 1e-8)

    K_tr_pd <- as.matrix(nearPD(K_tr, ensureSymmetry = TRUE)$mat)
    L        <- chol(K_tr_pd)
    alpha    <- backsolve(L, forwardsolve(t(L), y_tr - mu_tr))
    cond_mu  <- mu_new + as.numeric(K_cross %*% alpha)

    V        <- forwardsolve(t(L), t(K_cross))
    cond_K   <- K_nn - t(V) %*% V
    cond_K   <- (cond_K + t(cond_K)) / 2
    cond_K_pd <- as.matrix(nearPD(cond_K, ensureSymmetry = TRUE,
                                   eig.tol = 1e-10)$mat)

    tryCatch({
      y_pred_std[ii, ] <- MASS::mvrnorm(1, cond_mu, cond_K_pd)
    }, error = function(e) {
      y_pred_std[ii, ] <<- rnorm(S_new, cond_mu, sqrt(pmax(diag(cond_K_pd), 1e-8)))
    })
  }

  y_pred_std * sd_train$meta$sy + sd_train$meta$my
}

predict_no_gp <- function(fit_train, sd_train, sd_new,
                           use_eiv = TRUE, n_draws = 400L) {
  post  <- rstan::extract(fit_train)
  n_all <- nrow(post$beta)
  idx   <- sample(n_all, min(n_draws, n_all))
  S_new <- sd_new$stan$S

  y_pred_std <- matrix(NA_real_, length(idx), S_new)
  for (ii in seq_along(idx)) {
    it      <- idx[ii]
    beta_it <- post$beta[it, ]

    if (use_eiv) {
      mu_k_all  <- post$mu_x[it, ]
      sg_k_all  <- post$sigma_x[it, ]
      X_star_new <- matrix(NA_real_, S_new, sd_new$stan$K)
      for (k in seq_len(sd_new$stan$K)) {
        tau_k    <- sd_new$stan$tau_x[, k]
        v_pr     <- sg_k_all[k]^2; v_me <- tau_k^2
        post_mu  <- (v_me * mu_k_all[k] + v_pr * sd_new$stan$X_obs[, k]) / (v_pr + v_me)
        post_var <- v_pr * v_me / (v_pr + v_me)
        X_star_new[, k] <- rnorm(S_new, post_mu, sqrt(post_var))
      }
      mu_new <- as.numeric(X_star_new %*% beta_it)
    } else {
      mu_new <- as.numeric(sd_new$stan$X_obs %*% beta_it)
    }
    y_pred_std[ii, ] <- rnorm(S_new, mu_new, sd_new$stan$tau_y)
  }
  y_pred_std * sd_train$meta$sy + sd_train$meta$my
}

eval_pred <- function(pred_mat, y_obs_raw, label) {
  y_mean <- colMeans(pred_mat)
  rmspe  <- sqrt(mean((y_mean - y_obs_raw)^2))
  crps_v <- scoringRules::crps_sample(y = y_obs_raw, dat = t(pred_mat))
  data.frame(direction = label, rmspe = rmspe, crps = mean(crps_v))
}

## ----------------------------------------------------------------
## Run cross-month predictions
## ----------------------------------------------------------------

perf_rows <- list()

p1 <- gp_predict_robust(fit_full_jul, sd_jul, sd_jan, "M_full")
perf_rows[["full_j2j"]] <- eval_pred(p1, sd_jan$meta$y_raw, "M_full   Jul->Jan")

p2 <- gp_predict_robust(fit_full_jan, sd_jan, sd_jul, "M_full")
perf_rows[["full_j2u"]] <- eval_pred(p2, sd_jul$meta$y_raw, "M_full   Jan->Jul")

noEIV_predict <- function(fit_train, sd_train, sd_new, n_draws = 400L) {
  post  <- rstan::extract(fit_train)
  idx   <- sample(nrow(post$beta), min(n_draws, nrow(post$beta)))
  S_new <- sd_new$stan$S

  y_pred_std <- matrix(NA_real_, length(idx), S_new)
  for (ii in seq_along(idx)) {
    it      <- idx[ii]; beta_it <- post$beta[it, ]
    su      <- post$sigma_u[it]; el <- post$ell[it]
    mu_tr   <- as.numeric(sd_train$stan$X_obs %*% beta_it)
    mu_new  <- as.numeric(sd_new$stan$X_obs %*% beta_it)
    d_tr    <- sd_train$stan$d; d_new <- sd_new$stan$d
    tau_tr  <- sd_train$stan$tau_y; tau_new <- sd_new$stan$tau_y
    y_tr    <- sd_train$stan$y_obs

    K_tr    <- su^2*exp(-abs(outer(d_tr,d_tr,"-"))/el)+diag(tau_tr^2+1e-8)
    K_cross <- su^2*exp(-abs(outer(d_new,d_tr,"-"))/el)
    K_nn    <- su^2*exp(-abs(outer(d_new,d_new,"-"))/el)+diag(tau_new^2+1e-8)

    K_tr_pd  <- as.matrix(nearPD(K_tr, ensureSymmetry=TRUE)$mat)
    L        <- chol(K_tr_pd)
    alpha    <- backsolve(L, forwardsolve(t(L), y_tr-mu_tr))
    cond_mu  <- mu_new + as.numeric(K_cross %*% alpha)
    V        <- forwardsolve(t(L), t(K_cross))
    cond_K   <- K_nn - t(V) %*% V
    cond_K   <- (cond_K+t(cond_K))/2
    cond_K_pd<- as.matrix(nearPD(cond_K, ensureSymmetry=TRUE, eig.tol=1e-10)$mat)

    tryCatch({
      y_pred_std[ii,] <- MASS::mvrnorm(1, cond_mu, cond_K_pd)
    }, error = function(e) {
      y_pred_std[ii,] <<- rnorm(S_new, cond_mu, sqrt(pmax(diag(cond_K_pd), 1e-8)))
    })
  }
  y_pred_std * sd_train$meta$sy + sd_train$meta$my
}

sd_jan_noEIV <- sd_jan; sd_jan_noEIV$stan$tau_x <- NULL
sd_jul_noEIV <- sd_jul; sd_jul_noEIV$stan$tau_x <- NULL

p3 <- noEIV_predict(fit_noEIV_jul, sd_jul_noEIV, sd_jan_noEIV)
perf_rows[["noEIV_j2j"]] <- eval_pred(p3, sd_jan$meta$y_raw, "M_noEIV Jul->Jan")
p4 <- noEIV_predict(fit_noEIV_jan, sd_jan_noEIV, sd_jul_noEIV)
perf_rows[["noEIV_j2u"]] <- eval_pred(p4, sd_jul$meta$y_raw, "M_noEIV Jan->Jul")

p5 <- predict_no_gp(fit_noGP_jul, sd_jul, sd_jan, use_eiv=TRUE)
perf_rows[["noGP_j2j"]] <- eval_pred(p5, sd_jan$meta$y_raw, "M_noGP  Jul->Jan")
p6 <- predict_no_gp(fit_noGP_jan, sd_jan, sd_jul, use_eiv=TRUE)
perf_rows[["noGP_j2u"]] <- eval_pred(p6, sd_jul$meta$y_raw, "M_noGP  Jan->Jul")

p7 <- predict_no_gp(fit_base_jul, sd_jul, sd_jan, use_eiv=FALSE)
perf_rows[["base_j2j"]] <- eval_pred(p7, sd_jan$meta$y_raw, "M_base  Jul->Jan")
p8 <- predict_no_gp(fit_base_jan, sd_jan, sd_jul, use_eiv=FALSE)
perf_rows[["base_j2u"]] <- eval_pred(p8, sd_jul$meta$y_raw, "M_base  Jan->Jul")

perf_df_full <- bind_rows(perf_rows)
print(perf_df_full)

## ----------------------------------------------------------------
## Save predictions
## ----------------------------------------------------------------
predictions <- list(
  full_j2j  = p1, full_j2u  = p2,
  noEIV_j2j = p3, noEIV_j2u = p4,
  noGP_j2j  = p5, noGP_j2u  = p6,
  base_j2j  = p7, base_j2u  = p8
)

save(predictions, perf_df_full, file = here::here("Results", "predictions.RData"))
