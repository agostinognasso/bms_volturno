## ============================================================
## sim/run_sim.R — replication harness
## Loops scenarios x replicates x 4 models, fits Stan, scores
## beta / sigma_u / ell against the known truth, records MCMC
## diagnostics. Incremental rds save (crash-safe).
##
## Config via env vars (pilot defaults shown):
##   SIM_R         number of replicates per scenario   [3]
##   SIM_SCEN      comma list of scenarios   [S1_jan,S2_jul,S3_nogp,S4_miss]
##   SIM_CHAINS    MCMC chains                          [2]
##   SIM_ITER      MCMC iterations                      [600]
##   SIM_WARMUP    MCMC warmup                          [300]
##   SIM_TAG       output tag (pilot|full|...)          [pilot]
## Full run example:
##   SIM_R=200 SIM_CHAINS=4 SIM_ITER=1500 SIM_WARMUP=750 SIM_TAG=full \
##     Rscript sim/run_sim.R
## ============================================================
rm(list = ls())
library(here)
SIM_DIR  <- here::here("Simulation")
STAN_DIR <- here::here("Stan")
source(file.path(SIM_DIR, "dgp.R"))

suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)

geti  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
gets  <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }

R_REP   <- geti("SIM_R", 3L)
CHAINS  <- geti("SIM_CHAINS", 2L)
ITER    <- geti("SIM_ITER", 600L)
WARMUP  <- geti("SIM_WARMUP", 300L)
TAG     <- gets("SIM_TAG", "pilot")
SCEN    <- strsplit(gets("SIM_SCEN", "S1_jan,S2_jul,S3_nogp,S4_miss"), ",")[[1]]
## Sharding: NSHARD independent processes split the replicate indices
## (r mod NSHARD == SHARD). Each writes its own crash-safe rds; metrics.R
## rbinds all shards. SHARD is 0-based.
NSHARD  <- geti("SIM_NSHARD", 1L)
SHARD   <- geti("SIM_SHARD", 0L)
options(mc.cores = min(CHAINS, 4L))
BASE_SEED <- 20260518L

truth <- readRDS(file.path(SIM_DIR, "truth_params.rds"))

cat(sprintf("RUN tag=%s  R=%d  scen=%s  chains=%d iter=%d warmup=%d\n",
            TAG, R_REP, paste(SCEN, collapse="+"), CHAINS, ITER, WARMUP))

models <- list(
  M_full  = list(stan = "eiv_gp.stan",       eiv = TRUE,  gp = TRUE),
  M_noEIV = list(stan = "no_eiv_gp.stan",    eiv = FALSE, gp = TRUE),
  M_noGP  = list(stan = "eiv_no_gp.stan",    eiv = TRUE,  gp = FALSE),
  M_base  = list(stan = "no_eiv_no_gp.stan", eiv = FALSE, gp = FALSE)
)
sm <- lapply(models, function(m)
  stan_model(file.path(STAN_DIR, m$stan), model_name = m$stan))

## Strip tau_x for non-EIV models (matches run_analysis.R interface).
make_data <- function(stan_data, eiv) {
  d <- stan_data
  if (!eiv) d$tau_x <- NULL
  d
}

fit_one <- function(model_id, stan_data, seed) {
  m <- models[[model_id]]
  dat <- make_data(stan_data, m$eiv)
  fit <- tryCatch(
    sampling(sm[[model_id]], data = dat, chains = CHAINS, iter = ITER,
             warmup = WARMUP, thin = 1L, seed = seed, refresh = 0L,
             control = list(adapt_delta = 0.95, max_treedepth = 12)),
    error = function(e) { message("  fit error ", model_id, ": ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  post <- rstan::extract(fit)
  smry <- summary(fit)$summary
  list(
    beta_med = apply(post$beta, 2, median),
    beta_lo  = apply(post$beta, 2, quantile, 0.025),
    beta_hi  = apply(post$beta, 2, quantile, 0.975),
    sigma_u_med = if (m$gp) median(post$sigma_u) else NA_real_,
    ell_med     = if (m$gp) median(post$ell)     else NA_real_,
    rhat_max = max(smry[, "Rhat"], na.rm = TRUE),
    neff_min = min(smry[, "n_eff"], na.rm = TRUE),
    divergences = get_num_divergent(fit)
  )
}

rows <- list(); ri <- 1L
out_file <- if (NSHARD > 1L) {
  file.path(SIM_DIR, sprintf("sim_results_%s_shard%d.rds", TAG, SHARD))
} else {
  file.path(SIM_DIR, sprintf("sim_results_%s.rds", TAG))
}
my_reps <- seq_len(R_REP)[(seq_len(R_REP) - 1L) %% NSHARD == SHARD]
cat(sprintf("shard %d/%d handles %d reps\n", SHARD, NSHARD, length(my_reps)))
t0 <- Sys.time()

for (scn in SCEN) {
  for (r in my_reps) {
    seed_r <- BASE_SEED + 1000L * which(SCEN == scn) + r
    ds <- simulate_dataset(truth, scn, seed = seed_r)
    tb <- ds$truth$beta
    K  <- length(tb)
    cn <- ds$truth$cov_names
    for (mid in names(models)) {
      res <- fit_one(mid, ds$stan, seed = seed_r + 7L)
      if (is.null(res)) next
      for (k in seq_len(K)) {
        rows[[ri]] <- data.frame(
          tag = TAG, scenario = scn, base = ds$truth$base,
          rep = r, model = mid, cov = cn[k], k = k,
          beta_true = tb[k],
          beta_med  = res$beta_med[k],
          beta_lo   = res$beta_lo[k],
          beta_hi   = res$beta_hi[k],
          covered   = as.integer(tb[k] >= res$beta_lo[k] &
                                 tb[k] <= res$beta_hi[k]),
          sign_match = as.integer(sign(res$beta_med[k]) == sign(tb[k])),
          sigma_u_true = ds$truth$sigma_u,
          sigma_u_med  = res$sigma_u_med,
          ell_true     = ds$truth$ell,
          ell_med      = res$ell_med,
          rhat_max = res$rhat_max,
          neff_min = res$neff_min,
          divergences = res$divergences,
          stringsAsFactors = FALSE)
        ri <- ri + 1L
      }
    }
    saveRDS(do.call(rbind, rows), out_file)   # incremental, crash-safe
    cat(sprintf("[%s] scen=%s rep=%d/%d done (%.1f min elapsed)\n",
                TAG, scn, r, R_REP,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

cat("\nDONE. Wrote", out_file, "with", ri - 1L, "rows\n")
