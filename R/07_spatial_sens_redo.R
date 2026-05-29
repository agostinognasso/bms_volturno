## 07_spatial_sens_redo.R --- M4 fix (correct station IDs)
##
## The original 06_spatial_and_tau_sens.R used the wrong filter strings
## ("FV-CeA" etc.) but the raw data uses "FV - CeA" (with spaces and
## "Grazz" not "Grazzi"). This script re-runs M4 with the correct IDs
## and overwrites Results/spatial_sens.rds and .tex.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(here)
})
rstan_options(auto_write = TRUE); options(mc.cores = 4L)

SEED <- 1234L
source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jul(), build_bins(), aggregate_to_bins(), vars_phys

build_stan_data <- function(df) {
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
  list(S=nrow(X_obs), K=ncol(X_obs),
       y_obs=as.vector(y_obs), X_obs=X_obs,
       tau_y=as.vector(tau_y), tau_x=tau_x, d=as.vector(d))
}

sm <- stan_model(here::here("Stan","eiv_gp.stan"))

fit_summary <- function(sd_obj, tag) {
  fit <- sampling(sm, data=sd_obj, chains=4, iter=2000, warmup=1000,
                  seed=SEED, control=list(adapt_delta=0.95, max_treedepth=12),
                  refresh=0)
  post <- rstan::extract(fit)
  q <- function(x) c(median(x), quantile(x,.025), quantile(x,.975))
  data.frame(tag=tag, S=sd_obj$S,
             bT_med=q(post$beta[,1])[1], bT_lo=q(post$beta[,1])[2], bT_hi=q(post$beta[,1])[3],
             bO_med=q(post$beta[,2])[1], bO_lo=q(post$beta[,2])[2], bO_hi=q(post$beta[,2])[3],
             su_med=q(post$sigma_u)[1], su_lo=q(post$sigma_u)[2], su_hi=q(post$sigma_u)[3])
}

jul_full <- load_jul()
print(unique(jul_full$ID))

# Real IDs in the data (with spaces!)
far_up <- c("FV - CeA","FV - Grazz","FV - Capua")
jul_near <- jul_full %>% filter(!(ID %in% far_up))
print(unique(jul_near$ID))

sd_full <- build_stan_data(jul_full)
sd_near <- build_stan_data(jul_near)

res_full <- fit_summary(sd_full, "Jul-25 (all stations)")
res_near <- fit_summary(sd_near, "Jul-25 (FV-CeA/Grazz/Capua excluded)")

out <- rbind(res_full, res_near)
saveRDS(out, here::here("Results","spatial_sens.rds"))

fmt_ci <- function(m,lo,hi) sprintf("%.2f [%.2f,\\,%.2f]", m, lo, hi)
rows <- character(nrow(out))
for (i in seq_len(nrow(out))) {
  rows[i] <- sprintf("    %s & %s & %s & %s \\\\",
    out$tag[i],
    fmt_ci(out$bT_med[i], out$bT_lo[i], out$bT_hi[i]),
    fmt_ci(out$bO_med[i], out$bO_lo[i], out$bO_hi[i]),
    fmt_ci(out$su_med[i], out$su_lo[i], out$su_hi[i]))
}
writeLines(rows, here::here("Results","spatial_sens.tex"))
print(out, digits=3)
