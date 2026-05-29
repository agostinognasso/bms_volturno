## 05_prior_sensitivity.R  --  M11 of the SERRA revision
##
## Refit M_full on Jan-26 under a grid of (s_beta, s_sigma_u) prior
## scales; report posterior median and 95% CI for the two identified
## tracer coefficients and for sigma_u. Output:
##   Results/prior_sensitivity_jan.rds   (raw posterior summaries)
##   Results/prior_sensitivity_jan.tex   (LaTeX rows for Table prior_sens)
##
## Wall-time budget: ~25-35 min on M4 Pro (4 chains, 1000 warmup,
## 1000 sampling per fit, 5 prior settings, only Jan-26).
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(here)
})
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

STAN_FILE <- here::here("Stan", "eiv_gp_priors.stan")
SEED     <- 1234L

source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jan(), build_bins(), aggregate_to_bins(), vars_phys

df  <- load_jan()
bins<- build_bins(df)
agg <- aggregate_to_bins(df, bins)
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

base_data <- list(S=nrow(X_obs), K=ncol(X_obs),
                  y_obs=as.vector(y_obs), X_obs=X_obs,
                  tau_y=as.vector(tau_y), tau_x=tau_x,
                  d=as.vector(d))

grid <- data.frame(
  s_beta    = c(0.25, 0.5, 1.0, 0.5, 0.5),
  s_sigma_u = c(1.0,  1.0, 1.0, 0.5, 2.0)
)

sm <- stan_model(STAN_FILE)

fmt <- function(x) sprintf("%.3f", x)
fmt_ci <- function(med, lo, hi) sprintf("%.2f [%.2f,\\,%.2f]", med, lo, hi)
results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  sd_i <- base_data
  sd_i$prior_s_beta    <- grid$s_beta[i]
  sd_i$prior_s_sigma_u <- grid$s_sigma_u[i]
  fit <- sampling(sm, data=sd_i, chains=4, iter=2000, warmup=1000,
                  seed=SEED+i, control=list(adapt_delta=0.95, max_treedepth=12),
                  refresh=0)
  post <- rstan::extract(fit)
  betaT  <- post$beta[,1]; betaO2 <- post$beta[,2]; sigmau <- post$sigma_u
  results[[i]] <- data.frame(
    s_beta=grid$s_beta[i], s_sigma_u=grid$s_sigma_u[i],
    bT_med=median(betaT),  bT_lo=quantile(betaT,.025),  bT_hi=quantile(betaT,.975),
    bO_med=median(betaO2), bO_lo=quantile(betaO2,.025), bO_hi=quantile(betaO2,.975),
    su_med=median(sigmau), su_lo=quantile(sigmau,.025), su_hi=quantile(sigmau,.975)
  )
}
out <- do.call(rbind, results)
saveRDS(out, here::here("Results", "prior_sensitivity_jan.rds"))

prior_lbl <- function(s, kind) {
  if (kind == "beta") sprintf("$\\mathcal{N}(0,%.2f^2)$", s)
  else sprintf("$\\mathrm{HN}(0,%.1f)$", s)
}
ref_idx <- which(out$s_beta == 0.5 & out$s_sigma_u == 1.0)
rows <- character(nrow(out))
for (i in seq_len(nrow(out))) {
  pre_open <- if (i == ref_idx) "    \\textbf{" else "    "
  pre_close <- if (i == ref_idx) "}" else ""
  row <- sprintf("%s%s & %s & %s & %s & %s%s \\\\",
    pre_open,
    prior_lbl(out$s_beta[i], "beta"),
    prior_lbl(out$s_sigma_u[i], "su"),
    fmt_ci(out$bT_med[i], out$bT_lo[i], out$bT_hi[i]),
    fmt_ci(out$bO_med[i], out$bO_lo[i], out$bO_hi[i]),
    fmt_ci(out$su_med[i], out$su_lo[i], out$su_hi[i]),
    pre_close)
  rows[i] <- row
}
writeLines(rows, here::here("Results", "prior_sensitivity_jan.tex"))
print(out, digits=3)
