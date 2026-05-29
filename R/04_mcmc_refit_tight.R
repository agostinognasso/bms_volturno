## 04_mcmc_refit_tight.R --- M3 of the SERRA revision
##
## Refit M_full on Jan-26 with tighter HMC tuning:
##   adapt_delta = 0.99, max_treedepth = 15, 4000 post-warmup iterations.
## Report updated diagnostics (R-hat, n_eff, divergences, BFMI) and
## produce a pair plot of the parameters involved in the divergences.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(bayesplot); library(here)
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


diag_one <- function(df, label) {
  sd_obj <- build_stan_data(df)
  fit <- sampling(sm, data=sd_obj, chains=4, iter=5000, warmup=1000,
                  seed=SEED, control=list(adapt_delta=0.99, max_treedepth=15),
                  refresh=200)
  sumfit <- summary(fit)$summary
  rhat_max <- max(sumfit[,"Rhat"], na.rm=TRUE)
  neff_min <- min(sumfit[,"n_eff"], na.rm=TRUE)
  div <- sum(sapply(get_sampler_params(fit), function(x) sum(x[,"divergent__"])))
  bfmi <- mean(rstan::get_bfmi(fit), na.rm=TRUE)
  if (label == "Jan-26") {
    pdf(here::here("Results", "fig_pairs_mfull_jan.pdf"), width=7.5, height=7.5)
    print(bayesplot::mcmc_pairs(as.array(fit),
              pars=c("beta[1]","beta[2]","sigma_u","ell"),
              np=nuts_params(fit),
              off_diag_args=list(size=0.4, alpha=0.5)))
    dev.off()
  }
  list(rhat_max=rhat_max, neff_min=neff_min, div=div, bfmi=bfmi, label=label)
}

dj <- diag_one(load_jan(), "Jan-26")
du <- diag_one(load_jul(), "Jul-25")
saveRDS(list(jan=dj, jul=du),
        here::here("Results","mcmc_tight_diag.rds"))
