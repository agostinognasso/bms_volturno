## ============================================================
## sim/metrics.R — aggregate replication results into the
## summary tables for the "Simulation study" section.
## Usage: Rscript sim/metrics.R [tag]   (default tag = full)
## ============================================================
suppressMessages(library(dplyr))
suppressMessages(library(here))
SIM_DIR <- here::here("Simulation")
args <- commandArgs(trailingOnly = TRUE)
TAG  <- if (length(args) >= 1) args[1] else "full"

# Merge ALL shard files of ANY tag (e.g. 'full' for S1/S2/S4 plus a
# separate 's3' top-up run), then de-duplicate on the replicate key.
# Deterministic per-(scenario,rep) seeds guarantee identical rows for
# any overlap, so keeping the first occurrence is exact.
all_shards <- list.files(SIM_DIR,
  pattern = "^sim_results_.*_shard[0-9]+\\.rds$", full.names = TRUE)
single <- file.path(SIM_DIR, sprintf("sim_results_%s.rds", TAG))
if (length(all_shards) > 0) {
  res <- do.call(rbind, lapply(all_shards, readRDS))
  key <- paste(res$scenario, res$rep, res$model, res$k, sep = "|")
  res <- res[!duplicated(key), ]
  cat("Merged", length(all_shards), "shard files;",
      nrow(res), "unique rows\n")
} else {
  res <- readRDS(single)
}
res$err <- res$beta_med - res$beta_true

## --- Table A: beta recovery per scenario x model (all covariates) ---
tabA <- res %>%
  group_by(scenario, model) %>%
  summarise(
    n_reps   = n_distinct(rep),
    bias     = mean(err),
    rmse     = sqrt(mean(err^2)),
    coverage = mean(covered),
    sign_ok  = mean(sign_match),
    .groups  = "drop")

## --- Table B: the two identified tracers (Temperature, O2) ---
tabB <- res %>%
  filter(cov %in% c("temperatura", "ossigeno_mgkg")) %>%
  group_by(scenario, model, cov) %>%
  summarise(
    beta_true = mean(beta_true),
    bias      = mean(err),
    rmse      = sqrt(mean(err^2)),
    coverage  = mean(covered),
    sign_rev  = mean(1 - sign_match),
    .groups   = "drop")

## --- Table C: GP hyperparameter recovery (M_full only) ---
tabC <- res %>%
  filter(model == "M_full", k == 1) %>%
  group_by(scenario) %>%
  summarise(
    sigma_u_true = mean(sigma_u_true),
    sigma_u_med  = mean(sigma_u_med, na.rm = TRUE),
    ell_true     = mean(ell_true),
    ell_med      = mean(ell_med, na.rm = TRUE),
    .groups      = "drop")

## --- Table D: MCMC health per scenario x model ---
tabD <- res %>%
  filter(k == 1) %>%
  group_by(scenario, model) %>%
  summarise(
    rhat_max_med = median(rhat_max),
    neff_min_med = median(neff_min),
    div_mean     = mean(divergences),
    .groups      = "drop")

cat("\n===== TABLE A: beta recovery (all covariates) =====\n")
print(as.data.frame(tabA), digits = 3)
cat("\n===== TABLE B: tracers (Temperature, O2) =====\n")
print(as.data.frame(tabB), digits = 3)
cat("\n===== TABLE C: GP hyperparameters (M_full) =====\n")
print(as.data.frame(tabC), digits = 3)
cat("\n===== TABLE D: MCMC diagnostics =====\n")
print(as.data.frame(tabD), digits = 3)

saveRDS(list(A = tabA, B = tabB, C = tabC, D = tabD),
        file.path(SIM_DIR, sprintf("sim_summary_%s.rds", TAG)))
cat("\nSaved sim_summary_", TAG, ".rds\n", sep = "")
