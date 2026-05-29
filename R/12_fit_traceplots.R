## 12_fit_traceplots.R
##
## Fits M_full on both campaigns with the long-run MCMC settings used
## in the paper (4 chains x 24000 iterations, 4000 warmup) and produces
## traceplots for the 9 scalar parameters (beta[1:7], sigma_u, ell).
##
## Output:
##   Results/fig_traceplots_full_jan.pdf
##   Results/fig_traceplots_full_jul.pdf
##   ../Journal_Submission/Figures/fig_traceplots_full_jan.pdf
##   ../Journal_Submission/Figures/fig_traceplots_full_jul.pdf
##
## NOTE: This script runs 2 x 4 chains x 24000 iterations and takes
## several hours on a standard laptop. Run it once and the output is
## saved. The traceplots do not need to be regenerated unless the
## model or data change.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(rstan); library(bayesplot); library(here)
})
rstan_options(auto_write = TRUE)
options(mc.cores = 4L)

source(here::here("R", "00_data_config.R"))

SEED     <- 1234L
DATA_DIR <- here::here("Data")
RES_DIR  <- here::here("Results")
check_data(DATA_DIR)
FIG_DIR  <- file.path(dirname(here::here()), "Journal_Submission", "Figures")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jan(), load_jul(), build_bins(), aggregate_to_bins(), vars_phys

## ── Build Stan data ───────────────────────────────────────────────────────────
build_stan_data <- function(loader) {
  df   <- loader()
  bins <- build_bins(df)
  agg  <- aggregate_to_bins(df, bins)
  cov_mean <- paste0(vars_phys, "_mean")
  cov_se   <- paste0(vars_phys, "_se")
  dat <- agg |>
    select(depth_mid, salinita_mean, salinita_se,
           all_of(cov_mean), all_of(cov_se)) |>
    filter(if_all(everything(), ~ !is.na(.x)))
  y_raw <- dat$salinita_mean; sy <- sd(y_raw); my <- mean(y_raw)
  y_obs <- (y_raw - my) / sy
  tau_y <- pmax(dat$salinita_se / sy, 1e-6)
  X_raw <- as.matrix(dat[, cov_mean]); TauX <- as.matrix(dat[, cov_se])
  mx <- colMeans(X_raw); sx <- apply(X_raw, 2, sd)
  X_obs <- sweep(sweep(X_raw, 2, mx, "-"), 2, sx, "/")
  tau_x <- pmax(sweep(TauX, 2, sx, "/"), 1e-6)
  d_raw <- dat$depth_mid
  d     <- (d_raw - min(d_raw)) / (max(d_raw) - min(d_raw))
  list(S = nrow(X_obs), K = ncol(X_obs),
       y_obs = as.vector(y_obs), X_obs = X_obs,
       tau_y = as.vector(tau_y), tau_x = tau_x,
       d = as.vector(d))
}

## ── Compile Stan model ────────────────────────────────────────────────────────
sm <- stan_model(here::here("Stan", "eiv_gp.stan"))

## ── Fit function ──────────────────────────────────────────────────────────────
fit_campaign <- function(loader, label) {
  message(sprintf("\nFitting M_full — %s (4 chains x 24000 iter, 4000 warmup) ...", label))
  sd_obj <- build_stan_data(loader)
  sampling(
    sm, data = sd_obj,
    chains  = 4L, iter = 24000L, warmup = 4000L,
    seed    = SEED,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 500L
  )
}

## ── Traceplot function ────────────────────────────────────────────────────────
pars_to_plot <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                  "beta[5]", "beta[6]", "beta[7]", "sigma_u", "ell")

strip_labels <- c(
  "beta[1]" = "beta[Temp]",
  "beta[2]" = "beta[O[2]]",
  "beta[3]" = "beta[pH]",
  "beta[4]" = "beta[Chl-a]",
  "beta[5]" = "beta[PC]",
  "beta[6]" = "beta[PE]",
  "beta[7]" = "beta[Dist]",
  "sigma_u" = "sigma[u]",
  "ell"     = "ell"
)

make_traceplots <- function(fit, label, outname) {
  arr <- as.array(fit, pars = pars_to_plot)
  np  <- nuts_params(fit)
  dimnames(arr)[[3]] <- strip_labels[pars_to_plot]

  p <- mcmc_trace(
    arr, np = np,
    facet_args = list(ncol = 3, labeller = label_parsed)
  ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title    = paste0("MCMC traceplots — M_full (", label, ")"),
      subtitle = "Divergent transitions highlighted in red"
    ) +
    theme_default(base_size = 8) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position  = "bottom",
      plot.title       = element_text(size = 9, face = "bold"),
      plot.subtitle    = element_text(size = 7.5)
    )

  for (path in c(file.path(RES_DIR, outname), file.path(FIG_DIR, outname))) {
    pdf(path, width = 7.2, height = 8.0); print(p); dev.off()
  }
  message("Wrote ", outname)

  # Quick convergence summary
  s    <- summary(fit, pars = pars_to_plot)$summary
  rhat <- max(s[, "Rhat"], na.rm = TRUE)
  neff <- min(s[, "n_eff"], na.rm = TRUE)
  ndiv <- sum(sapply(get_sampler_params(fit),
                     function(x) sum(x[, "divergent__"])))
  message(sprintf("  %s — Rhat_max=%.3f | neff_min=%.0f | divergences=%d",
                  label, rhat, neff, ndiv))
}

## ── Run ───────────────────────────────────────────────────────────────────────
fit_jan <- fit_campaign(load_jan, "Jan-26")
make_traceplots(fit_jan, "Jan-26", "fig_traceplots_full_jan.pdf")

fit_jul <- fit_campaign(load_jul, "Jul-25")
make_traceplots(fit_jul, "Jul-25", "fig_traceplots_full_jul.pdf")

message("\nTraceplots done. Saved to Results/ and Journal_Submission/Figures/")
