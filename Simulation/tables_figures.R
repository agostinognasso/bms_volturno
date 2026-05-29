## ============================================================
## sim/tables_figures.R
## Turns sim_summary_<tag>.rds into LaTeX fragments + a figure
## for the "Simulation study" section. Run AFTER metrics.R.
## Usage: Rscript sim/tables_figures.R [tag]   (default: full)
## Outputs (in sim/):
##   tab_sim_tracers.tex   beta recovery for Temperature & O2
##   tab_sim_gp.tex        sigma_u / ell recovery (M_full)
##   tab_sim_diag.tex      MCMC diagnostics per scenario x model
##   fig_sim_recovery.pdf  bias +/- and coverage, the two tracers
## ============================================================
suppressMessages(library(dplyr))
suppressMessages(library(here))
SIM_DIR <- here::here("Simulation")
args <- commandArgs(trailingOnly = TRUE)
TAG  <- if (length(args) >= 1) args[1] else "full"

S <- readRDS(file.path(SIM_DIR, sprintf("sim_summary_%s.rds", TAG)))

scen_lab <- c(S1_jan = "S1 (Jan-like, $\\sigma_u\\!\\approx\\!0.02$)",
              S2_jul = "S2 (Jul-like, $\\sigma_u\\!\\approx\\!0.52$)",
              S3_nogp = "S3 ($\\sigma_u\\!=\\!0$)",
              S4_miss = "S4 (misspecified)")
mod_lab  <- c(M_full = "$\\mathrm{M_{full}}$", M_noEIV = "$\\mathrm{M_{noEIV}}$",
              M_noGP = "$\\mathrm{M_{noGP}}$", M_base = "$\\mathrm{M_{base}}$")
cov_lab  <- c(temperatura = "Temp.", ossigeno_mgkg = "O$_2$")
f2 <- function(x) sprintf("%.2f", x)
f3 <- function(x) sprintf("%.3f", x)

## ---- Table: tracer recovery (the central claim) ----
B <- S$B %>%
  mutate(scen = scen_lab[scenario], mod = mod_lab[model], cv = cov_lab[cov]) %>%
  arrange(match(scenario, names(scen_lab)),
          match(cov, names(cov_lab)),
          match(model, names(mod_lab)))

con <- file(file.path(SIM_DIR, "tab_sim_tracers.tex"), "w")
writeLines(c(
"\\begin{table}[h!]",
"  \\caption{Simulation study: recovery of the two identified tracer",
"           coefficients under each scenario. Truth calibrated on the real",
"           $\\mathrm{M_{full}}$ posteriors. Bias and RMSE on the standardised",
"           scale; Cov.\\ is empirical coverage of the nominal 95\\,\\% credible",
"           interval; SR is the sign-reversal frequency.}",
"  \\label{tab:sim_tracers}",
"  \\centering",
"  \\small",
"  \\begin{tabular}{lllrrrr}",
"    \\toprule",
"    Scenario & Coef. & Model & Bias & RMSE & Cov. & SR \\\\",
"    \\midrule"), con)
prev <- ""
for (i in seq_len(nrow(B))) {
  r <- B[i, ]
  key <- paste(r$scenario, r$cov)
  sc  <- if (key != prev) r$scen else ""
  cvc <- if (key != prev) r$cv else ""
  prev <- key
  writeLines(sprintf("    %s & %s & %s & %s & %s & %s & %s \\\\",
                      sc, cvc, r$mod, f3(r$bias), f3(r$rmse),
                      f2(r$coverage), f2(r$sign_rev)), con)
  if (r$model == "M_base" && i < nrow(B)) writeLines("    \\addlinespace", con)
}
writeLines(c("    \\bottomrule", "  \\end{tabular}", "\\end{table}"), con)
close(con)

## ---- Table: GP hyperparameter recovery ----
C <- S$C %>% mutate(scen = scen_lab[scenario]) %>%
  arrange(match(scenario, names(scen_lab)))
con <- file(file.path(SIM_DIR, "tab_sim_gp.tex"), "w")
writeLines(c(
"\\begin{table}[h!]",
"  \\caption{Simulation study: recovery of the GP hyperparameters by",
"           $\\mathrm{M_{full}}$ (posterior-median averaged over replicates).}",
"  \\label{tab:sim_gp}",
"  \\centering",
"  \\begin{tabular}{lcccc}",
"    \\toprule",
"    Scenario & $\\sigma_u$ true & $\\hat{\\sigma}_u$ & $\\ell$ true & $\\hat{\\ell}$ \\\\",
"    \\midrule"), con)
for (i in seq_len(nrow(C))) {
  r <- C[i, ]
  writeLines(sprintf("    %s & %s & %s & %s & %s \\\\", r$scen,
                      f3(r$sigma_u_true), f3(r$sigma_u_med),
                      f2(r$ell_true), f2(r$ell_med)), con)
}
writeLines(c("    \\bottomrule", "  \\end{tabular}", "\\end{table}"), con)
close(con)

## ---- Table: MCMC diagnostics ----
D <- S$D %>% mutate(scen = scen_lab[scenario], mod = mod_lab[model]) %>%
  arrange(match(scenario, names(scen_lab)), match(model, names(mod_lab)))
con <- file(file.path(SIM_DIR, "tab_sim_diag.tex"), "w")
writeLines(c(
"\\begin{table}[h!]",
"  \\caption{Simulation study: MCMC diagnostics (median over replicates).}",
"  \\label{tab:sim_diag}",
"  \\centering",
"  \\begin{tabular}{llccc}",
"    \\toprule",
"    Scenario & Model & $\\hat{R}_{\\max}$ & $n_{\\mathrm{eff,min}}$ & Div \\\\",
"    \\midrule"), con)
prev <- ""
for (i in seq_len(nrow(D))) {
  r <- D[i, ]
  sc <- if (r$scenario != prev) r$scen else ""; prev <- r$scenario
  writeLines(sprintf("    %s & %s & %s & %.0f & %.1f \\\\", sc, r$mod,
                      f3(r$rhat_max_med), r$neff_min_med, r$div_mean), con)
  if (r$model == "M_base" && i < nrow(D)) writeLines("    \\addlinespace", con)
}
writeLines(c("    \\bottomrule", "  \\end{tabular}", "\\end{table}"), con)
close(con)

## ---- Figure: bias (with RMSE band) and coverage, two tracers ----
pdf(file.path(SIM_DIR, "fig_sim_recovery.pdf"), width = 7.2, height = 4.6)
op <- par(mfrow = c(1, 2), mar = c(7, 4, 2, 1), las = 2)
Bf <- B %>% mutate(grp = paste(scenario, cov, model))
cols <- c(M_full = "black", M_noEIV = "#1f77b4",
          M_noGP = "#2ca02c", M_base = "#d62728")
## panel 1: bias +/- RMSE
xs <- seq_len(nrow(Bf))
plot(xs, Bf$bias, pch = 19, col = cols[Bf$model],
     ylim = range(c(Bf$bias - Bf$rmse, Bf$bias + Bf$rmse)),
     xaxt = "n", xlab = "", ylab = "Bias (+/- RMSE)",
     main = "Coefficient bias")
abline(h = 0, lty = 2, col = "grey50")
arrows(xs, Bf$bias - Bf$rmse, xs, Bf$bias + Bf$rmse,
       angle = 90, code = 3, length = 0.02, col = cols[Bf$model])
axis(1, at = xs, labels = paste(Bf$scenario, Bf$cov, Bf$model), cex.axis = 0.5)
## panel 2: coverage
plot(xs, Bf$coverage, pch = 19, col = cols[Bf$model], ylim = c(0, 1),
     xaxt = "n", xlab = "", ylab = "95% CI coverage",
     main = "Credible-interval coverage")
abline(h = 0.95, lty = 2, col = "grey50")
axis(1, at = xs, labels = paste(Bf$scenario, Bf$cov, Bf$model), cex.axis = 0.5)
par(op); dev.off()

cat("Wrote tab_sim_tracers.tex, tab_sim_gp.tex, tab_sim_diag.tex,",
    "fig_sim_recovery.pdf into sim/\n")
