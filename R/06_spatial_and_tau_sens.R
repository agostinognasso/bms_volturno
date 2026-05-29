## 06_spatial_and_tau_sens.R --- M4 + M5 of the SERRA revision
##
## M4. Refit M_full on Jul-25 excluding the three far-upstream stations
##     (FV-CeA, FV-Grazzi, FV-Capua) that lie well beyond the marine
##     influence. Check stability of the identified tracer coefficients.
##
## M5. Refit M_full on Jan-26 and Jul-25 with tau_y and tau_x uniformly
##     scaled by 0.75 and 1.25 to verify robustness of the EIV layer to
##     misspecification of the plug-in measurement-error scales.
##
## Output:
##   Results/spatial_sens.rds, Results/spatial_sens.tex
##   Results/tau_sens.rds,     Results/tau_sens.tex
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(here)
})
rstan_options(auto_write = TRUE); options(mc.cores = 4L)

SEED <- 1234L
source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jan(), load_jul(), build_bins(), aggregate_to_bins(), vars_phys

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
  data.frame(tag=tag,
             bT_med=q(post$beta[,1])[1], bT_lo=q(post$beta[,1])[2], bT_hi=q(post$beta[,1])[3],
             bO_med=q(post$beta[,2])[1], bO_lo=q(post$beta[,2])[2], bO_hi=q(post$beta[,2])[3],
             su_med=q(post$sigma_u)[1], su_lo=q(post$sigma_u)[2], su_hi=q(post$sigma_u)[3])
}

# =============== M4: spatial sensitivity (Jul-25) =======================
jul_full <- load_jul()
sd_full  <- build_stan_data(jul_full)
res_full <- fit_summary(sd_full, "Jul-25 (all stations)")

far_up <- c("FV-CeA","FV-Grazzi","FV-Capua")
jul_near <- jul_full %>% filter(!(ID %in% far_up))
sd_near  <- build_stan_data(jul_near)
res_near <- fit_summary(sd_near, "Jul-25 (FV-CeA/Grazzi/Capua excluded)")

spatial <- rbind(res_full, res_near)
saveRDS(spatial, here::here("Results","spatial_sens.rds"))

fmt_ci <- function(m,lo,hi) sprintf("%.2f [%.2f,\\,%.2f]", m, lo, hi)
rows_sp <- character(nrow(spatial))
for (i in seq_len(nrow(spatial))) {
  rows_sp[i] <- sprintf("    %s & %s & %s & %s \\\\",
    spatial$tag[i],
    fmt_ci(spatial$bT_med[i], spatial$bT_lo[i], spatial$bT_hi[i]),
    fmt_ci(spatial$bO_med[i], spatial$bO_lo[i], spatial$bO_hi[i]),
    fmt_ci(spatial$su_med[i], spatial$su_lo[i], spatial$su_hi[i]))
}
writeLines(rows_sp, here::here("Results","spatial_sens.tex"))
print(spatial, digits=3)

# =============== M5: tau sensitivity (Jan-26 + Jul-25) ==================
sd_jan <- build_stan_data(load_jan())
sd_jul <- sd_full
scales <- c(0.75, 1.00, 1.25)
res_tau <- list()
for (cmp in c("Jan-26","Jul-25")) {
  sdb <- if (cmp=="Jan-26") sd_jan else sd_jul
  for (s in scales) {
    sd_i <- sdb
    sd_i$tau_y <- pmax(sdb$tau_y * s, 1e-6)
    sd_i$tau_x <- pmax(sdb$tau_x * s, 1e-6)
    res_tau[[paste(cmp,s)]] <- fit_summary(sd_i, sprintf("%s, tau x %.2f", cmp, s))
  }
}
tau_df <- do.call(rbind, res_tau)
saveRDS(tau_df, here::here("Results","tau_sens.rds"))
rows_tau <- character(nrow(tau_df))
for (i in seq_len(nrow(tau_df))) {
  rows_tau[i] <- sprintf("    %s & %s & %s & %s \\\\",
    tau_df$tag[i],
    fmt_ci(tau_df$bT_med[i], tau_df$bT_lo[i], tau_df$bT_hi[i]),
    fmt_ci(tau_df$bO_med[i], tau_df$bO_lo[i], tau_df$bO_hi[i]),
    fmt_ci(tau_df$su_med[i], tau_df$su_lo[i], tau_df$su_hi[i]))
}
writeLines(rows_tau, here::here("Results","tau_sens.tex"))
print(tau_df, digits=3)
