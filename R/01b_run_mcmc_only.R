## 01b_run_mcmc_only.R — Alternative Step 1: Fit all models and save before predictions.
## Use this instead of 01_run_analysis.R if you want to separate
## the MCMC fitting step from the prediction step (02_fix_predictions.R).
##
## Run from the Analysis/ directory (open Analysis.Rproj first).
rm(list = ls())

library(dplyr)
library(readxl)
library(rstan)
library(loo)
library(here)

source(here::here("R", "00_data_config.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

STAN_DIR <- here::here("Stan")
DATA_DIR <- here::here("Data")
SEED     <- 1234L
check_data(DATA_DIR)

safe_se <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

build_bins <- function(data, depth_col="profondita", id_col="ID",
                       min_depth=0.01, max_depth=10, base_step=0.01, min_sonde=3) {
  df <- data %>%
    filter(!is.na(.data[[depth_col]]), .data[[depth_col]]>=min_depth,
           .data[[depth_col]]<=max_depth) %>%
    mutate(base_bin=floor(.data[[depth_col]]/base_step)*base_step)
  base_counts <- df %>% group_by(base_bin) %>%
    summarise(sonde=list(sort(unique(.data[[id_col]]))),
              n_sonde=n_distinct(.data[[id_col]]), .groups="drop") %>% arrange(base_bin)
  bins <- list(); i <- 1L; k <- 1L
  while (i <= nrow(base_counts)) {
    left_edge <- base_counts$base_bin[i]; cur <- base_counts$sonde[[i]]; j <- i
    while (length(cur) < min_sonde && j < nrow(base_counts)) {
      j <- j+1L; cur <- union(cur, base_counts$sonde[[j]])
    }
    if (length(cur) >= min_sonde) {
      right_edge <- base_counts$base_bin[j]+base_step
      bins[[k]] <- data.frame(bin_id=k, left=left_edge, right=right_edge,
                               mid=(left_edge+right_edge)/2,
                               width=right_edge-left_edge, n_sonde=length(cur))
      k <- k+1L; i <- j+1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col="profondita") {
  out <- data; out$bin_id <- out$depth_mid <- NA_real_
  for (r in seq_len(nrow(bins))) {
    idx <- which(out[[depth_col]]>=bins$left[r] & out[[depth_col]]<bins$right[r])
    out$bin_id[idx] <- bins$bin_id[r]; out$depth_mid[idx] <- bins$mid[r]
  }
  out
}

vars_phys <- c("temperatura","salinita","ossigeno_mgkg","ph",
               "trx_chla","phycocyanin","phycoerythrin","dist")

load_jan <- function() {
  d <- read_excel(file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx"))
  df <- d %>%
    rename(profondita=`Profondità (m)`, temperatura=`temperatura (C°)`,
           salinita=`Salinità (‰)`, ossigeno_mgkg=`Oxygen (mg/Kg)`, ph=`pH`,
           trx_chla=`Trx-chl(a)`, phycocyanin=Phycocyanin,
           phycoerythrin=Phycoerythrin, dist=`Distanza dalla foce (km)`) %>%
    filter(!is.na(profondita), profondita>=0.1, profondita<=10)
  df$dist <- df$dist - min(df$dist, na.rm=TRUE); df
}

load_jul <- function() {
  d <- read_excel(file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx"))
  df <- d %>%
    rename(profondita=`Profondità (m)`, temperatura=`temperatura (°C)`,
           salinita=`Salinità (‰)`, ossigeno_mgkg=`Oxygen (mg/Kg)`, ph=`pH`,
           trx_chla=`Trx-chl(a)`, phycocyanin=Phycocyanin,
           phycoerythrin=Phycoerythrin, dist=`Distanza dalla Foce(km)`) %>%
    filter(!is.na(profondita), profondita>=0.1, profondita<=10)
  df$dist <- df$dist - min(df$dist, na.rm=TRUE); df
}

df_jan <- load_jan(); df_jul <- load_jul()

aggregate_to_bins <- function(df) {
  bins <- build_bins(df); df_bin <- assign_bins(df, bins)
  by_probe <- df_bin %>% filter(!is.na(bin_id)) %>%
    group_by(ID, bin_id, depth_mid) %>%
    summarise(across(all_of(vars_phys), ~mean(.x,na.rm=TRUE)), .groups="drop")
  by_probe %>% group_by(bin_id, depth_mid) %>%
    summarise(n_sonde=n(),
              across(all_of(vars_phys),
                     list(mean=~mean(.x,na.rm=TRUE), se=~safe_se(.x)),
                     .names="{.col}_{.fn}"), .groups="drop") %>% arrange(depth_mid)
}

agg_jan <- aggregate_to_bins(df_jan)
agg_jul <- aggregate_to_bins(df_jul)

cov_vars <- setdiff(vars_phys, "salinita")

build_stan_data <- function(agg) {
  cov_mean <- paste0(cov_vars,"_mean"); cov_se <- paste0(cov_vars,"_se")
  dat <- agg %>%
    select(depth_mid, salinita_mean, salinita_se, all_of(cov_mean), all_of(cov_se)) %>%
    filter(if_all(everything(), ~!is.na(.x)))
  y_raw <- dat$salinita_mean; sy <- sd(y_raw); my <- mean(y_raw)
  y_obs <- (y_raw-my)/sy; tau_y <- pmax(dat$salinita_se/sy, 1e-6)
  X_raw <- as.matrix(dat[,cov_mean]); TauX <- as.matrix(dat[,cov_se])
  mx <- colMeans(X_raw); sx <- apply(X_raw,2,sd)
  X_obs <- sweep(sweep(X_raw,2,mx,"-"),2,sx,"/")
  tau_x <- pmax(sweep(TauX,2,sx,"/"), 1e-6)
  d_raw <- dat$depth_mid; d <- (d_raw-min(d_raw))/(max(d_raw)-min(d_raw))
  list(stan=list(S=nrow(X_obs), K=ncol(X_obs),
                 y_obs=as.vector(y_obs), X_obs=X_obs,
                 tau_y=as.vector(tau_y), tau_x=tau_x, d=as.vector(d)),
       meta=list(my=my, sy=sy, mx=mx, sx=sx, d_raw=d_raw,
                 y_raw=y_raw, cov_names=cov_vars, n_bins=nrow(dat)))
}

sd_jan <- build_stan_data(agg_jan)
sd_jul <- build_stan_data(agg_jul)

sm_full  <- stan_model(file.path(STAN_DIR,"eiv_gp.stan"),      model_name="EIV_GP")
sm_noEIV <- stan_model(file.path(STAN_DIR,"no_eiv_gp.stan"),   model_name="noEIV_GP")
sm_noGP  <- stan_model(file.path(STAN_DIR,"eiv_no_gp.stan"),   model_name="EIV_noGP")
sm_base  <- stan_model(file.path(STAN_DIR,"no_eiv_no_gp.stan"),model_name="base")

stan_ctrl <- list(adapt_delta=0.95, max_treedepth=12)
no_eiv_d  <- function(sd) { d <- sd$stan; d$tau_x <- NULL; d }

fit_it <- function(sm, data, seed_off) {
  sampling(sm, data=data, chains=4L, iter=2000L, warmup=1000L,
           seed=SEED+seed_off, refresh=200L, control=stan_ctrl)
}


diag_row <- function(fit, name) {
  s <- summary(fit)$summary
  data.frame(model=name, rhat_max=max(s[,"Rhat"],na.rm=TRUE),
             neff_min=min(s[,"n_eff"],na.rm=TRUE),
             divergences=get_num_divergent(fit),
             bfmi=min(get_bfmi(fit)))
}
diag_df <- bind_rows(
  diag_row(fit_full_jan,"full_jan"), diag_row(fit_full_jul,"full_jul"),
  diag_row(fit_noEIV_jan,"noEIV_jan"), diag_row(fit_noEIV_jul,"noEIV_jul"),
  diag_row(fit_noGP_jan,"noGP_jan"),  diag_row(fit_noGP_jul,"noGP_jul"),
  diag_row(fit_base_jan,"base_jan"),  diag_row(fit_base_jul,"base_jul"))

waic_row <- function(fit, name) {
  ll <- extract_log_lik(fit, "log_lik")
  w  <- suppressWarnings(waic(ll))
  data.frame(model=name,
             elpd_waic=w$estimates["elpd_waic","Estimate"],
             waic_se=w$estimates["waic","SE"])
}
waic_df <- bind_rows(
  waic_row(fit_full_jan,"full_jan"), waic_row(fit_noEIV_jan,"noEIV_jan"),
  waic_row(fit_noGP_jan,"noGP_jan"), waic_row(fit_base_jan,"base_jan"),
  waic_row(fit_full_jul,"full_jul"), waic_row(fit_noEIV_jul,"noEIV_jul"),
  waic_row(fit_noGP_jul,"noGP_jul"), waic_row(fit_base_jul,"base_jul"))

loo_noGP_jan <- suppressWarnings(loo(extract_log_lik(fit_noGP_jan)))
loo_base_jan <- suppressWarnings(loo(extract_log_lik(fit_base_jan)))
loo_noGP_jul <- suppressWarnings(loo(extract_log_lik(fit_noGP_jul)))
loo_base_jul <- suppressWarnings(loo(extract_log_lik(fit_base_jul)))
loo_df <- data.frame(
  model = c("noGP_jan","base_jan","noGP_jul","base_jul"),
  elpd_loo = c(loo_noGP_jan$estimates["elpd_loo","Estimate"],
               loo_base_jan$estimates["elpd_loo","Estimate"],
               loo_noGP_jul$estimates["elpd_loo","Estimate"],
               loo_base_jul$estimates["elpd_loo","Estimate"]),
  loo_se = c(loo_noGP_jan$estimates["elpd_loo","SE"],
             loo_base_jan$estimates["elpd_loo","SE"],
             loo_noGP_jul$estimates["elpd_loo","SE"],
             loo_base_jul$estimates["elpd_loo","SE"]))

cov_labels <- c("Temperatura","Ossigeno","pH","Chl-a","Ficocianina","Ficoeritrina","Distanza")
beta_post <- function(fit, cov_names, model_name) {
  post <- rstan::extract(fit)
  bind_rows(lapply(seq_len(ncol(post$beta)), function(k) {
    b <- post$beta[,k]
    data.frame(model=model_name, covariate=cov_names[k],
               median=median(b), lo95=quantile(b,0.025), hi95=quantile(b,0.975),
               P_pos=mean(b>0))
  }))
}
beta_jan <- bind_rows(
  beta_post(fit_full_jan,cov_labels,"M_full"),
  beta_post(fit_noEIV_jan,cov_labels,"M_noEIV"),
  beta_post(fit_noGP_jan,cov_labels,"M_noGP"),
  beta_post(fit_base_jan,cov_labels,"M_base"))
beta_jul <- bind_rows(
  beta_post(fit_full_jul,cov_labels,"M_full"),
  beta_post(fit_noEIV_jul,cov_labels,"M_noEIV"),
  beta_post(fit_noGP_jul,cov_labels,"M_noGP"),
  beta_post(fit_base_jul,cov_labels,"M_base"))

desc_jan <- df_jan %>% group_by(ID) %>%
  summarise(n=n(), sal_mean=mean(salinita,na.rm=T), sal_sd=sd(salinita,na.rm=T),
            temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T),
            depth_max=max(profondita,na.rm=T), dist=first(dist), .groups="drop") %>%
  arrange(dist)
desc_jul <- df_jul %>% group_by(ID) %>%
  summarise(n=n(), sal_mean=mean(salinita,na.rm=T), sal_sd=sd(salinita,na.rm=T),
            temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T),
            depth_max=max(profondita,na.rm=T), dist=first(dist), .groups="drop") %>%
  arrange(dist)
stats_jan <- df_jan %>% summarise(n=n(), sal_mean=mean(salinita,na.rm=T),
  sal_sd=sd(salinita,na.rm=T), sal_med=median(salinita,na.rm=T),
  sal_min=min(salinita,na.rm=T), sal_max=max(salinita,na.rm=T),
  temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T))
stats_jul <- df_jul %>% summarise(n=n(), sal_mean=mean(salinita,na.rm=T),
  sal_sd=sd(salinita,na.rm=T), sal_med=median(salinita,na.rm=T),
  sal_min=min(salinita,na.rm=T), sal_max=max(salinita,na.rm=T),
  temp_mean=mean(temperatura,na.rm=T), o2_mean=mean(ossigeno_mgkg,na.rm=T))


save(fit_full_jan, fit_full_jul,
     fit_noEIV_jan, fit_noEIV_jul,
     fit_noGP_jan, fit_noGP_jul,
     fit_base_jan, fit_base_jul,
     sd_jan, sd_jul, agg_jan, agg_jul, df_jan, df_jul,
     diag_df, waic_df, loo_df, beta_jan, beta_jul,
     desc_jan, desc_jul, stats_jan, stats_jul,
     file=here::here("Results", "results.RData"))
