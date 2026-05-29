## 11_generate_paper_figures.R
##
## Generates all R-produced figures appearing in the manuscript:
##   fig_raw_profiles.pdf   — Fig 2: raw vertical profiles
##   fig_cross_month.pdf    — Fig 5: cross-month posterior predictive
##   fig_gp_component.pdf   — Fig 4: GP residual depth process u_s
##   fig_sim_recovery.pdf   — Fig 3: simulation bias/RMSE + coverage
##
## Reads from:
##   Data/          — Excel campaign files
##   Results/       — mfull_posterior_samples.rds, predictions.RData
##   Simulation/    — sim_summary_full.rds
##
## Writes figures to Results/ AND copies to
##   ../Journal_Submission/Figures/
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(ggplot2)
  library(tidyr); library(patchwork); library(here)
})

source(here::here("R", "00_data_config.R"))
DATA_DIR <- here::here("Data")
RES_DIR  <- here::here("Results")
check_data(DATA_DIR)
SIM_DIR  <- here::here("Simulation")
FIG_DIR  <- file.path(dirname(here::here()), "Journal_Submission", "Figures")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

save_fig <- function(p, name, width, height) {
  path_res <- file.path(RES_DIR, name)
  path_fig <- file.path(FIG_DIR, name)
  if (inherits(p, "gg") || inherits(p, "patchwork")) {
    ggsave(path_res, p, width = width, height = height, device = cairo_pdf)
    ggsave(path_fig, p, width = width, height = height, device = cairo_pdf)
  } else {
    pdf(path_res, width = width, height = height); print(p); dev.off()
    pdf(path_fig, width = width, height = height); print(p); dev.off()
  }
  message("Wrote ", name)
}

## ─────────────────────────────────────────────────────────────────────────────
## Data helpers
## ─────────────────────────────────────────────────────────────────────────────

vars_phys <- c("temperatura","ossigeno_mgkg","ph","trx_chla",
               "phycocyanin","phycoerythrin","dist")

safe_se <- function(x){
  n <- sum(!is.na(x)); if (n <= 1) NA_real_ else sd(x, na.rm = TRUE) / sqrt(n)
}

load_jan <- function(){
  read_excel(file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx")) |>
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (C°)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph = `pH`, trx_chla = `Trx-chl(a)`,
           phycocyanin = `Phycocyanin`, phycoerythrin = `Phycoerythrin`,
           dist = `Distanza dalla foce (km)`) |>
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10) |>
    mutate(dist = dist - min(dist, na.rm = TRUE))
}

load_jul <- function(){
  read_excel(file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx")) |>
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (°C)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph = `pH`, trx_chla = `Trx-chl(a)`,
           phycocyanin = `Phycocyanin`, phycoerythrin = `Phycoerythrin`,
           dist = `Distanza dalla Foce(km)`) |>
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10) |>
    mutate(dist = dist - min(dist, na.rm = TRUE))
}

build_bins <- function(data, depth_col = "profondita", id_col = "ID",
                       min_depth = 0.01, max_depth = 10,
                       base_step = 0.01, min_sonde = 3){
  df <- data |>
    filter(!is.na(.data[[depth_col]]),
           .data[[depth_col]] >= min_depth, .data[[depth_col]] <= max_depth) |>
    mutate(base_bin = floor(.data[[depth_col]] / base_step) * base_step)
  bc <- df |> group_by(base_bin) |>
    summarise(sonde = list(sort(unique(.data[[id_col]]))),
              n_sonde = n_distinct(.data[[id_col]]), .groups = "drop") |>
    arrange(base_bin)
  bins <- list(); i <- 1L; k <- 1L
  while (i <= nrow(bc)){
    le <- bc$base_bin[i]; cur <- bc$sonde[[i]]; j <- i
    while (length(cur) < min_sonde && j < nrow(bc)){
      j <- j + 1L; cur <- union(cur, bc$sonde[[j]])
    }
    if (length(cur) >= min_sonde){
      re <- bc$base_bin[j] + base_step
      bins[[k]] <- data.frame(bin_id = k, left = le, right = re,
                              mid = (le + re)/2, width = re - le,
                              n_sonde = length(cur))
      k <- k + 1L; i <- j + 1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col = "profondita"){
  out <- data; out$bin_id <- out$depth_mid <- NA_real_
  for (r in seq_len(nrow(bins))){
    idx <- which(out[[depth_col]] >= bins$left[r] & out[[depth_col]] < bins$right[r])
    out$bin_id[idx] <- bins$bin_id[r]; out$depth_mid[idx] <- bins$mid[r]
  }
  out
}

aggregate_to_bins <- function(df, bins){
  d <- assign_bins(df, bins)
  d |> filter(!is.na(bin_id)) |>
    group_by(bin_id, depth_mid) |>
    summarise(across(c(salinita, all_of(vars_phys)),
                     list(mean = ~mean(.x, na.rm = TRUE), se = ~safe_se(.x))),
              .groups = "drop") |>
    arrange(depth_mid)
}

## ─────────────────────────────────────────────────────────────────────────────
## Fig 2 — Raw vertical profiles (fig_raw_profiles.pdf)
## ─────────────────────────────────────────────────────────────────────────────

message("Generating fig_raw_profiles.pdf ...")

jan_raw <- load_jan()
jul_raw <- load_jul()

make_long <- function(df, label){
  df |>
    select(ID, profondita, salinita, temperatura, ossigeno_mgkg) |>
    rename(Depth = profondita, Station = ID,
           `Salinity (PSU)` = salinita,
           `Temperature (°C)` = temperatura,
           `O2 (mg/kg)` = ossigeno_mgkg) |>
    pivot_longer(c(`Salinity (PSU)`, `Temperature (°C)`, `O2 (mg/kg)`),
                 names_to = "Variable", values_to = "Value") |>
    mutate(Campaign = label,
           Variable = factor(Variable,
             levels = c("Salinity (PSU)", "Temperature (°C)", "O2 (mg/kg)")))
}

keep_jan <- c("FV-C01", "FV-Foce Volturno", "FV-03", "FV-06", "FV-10")
keep_jul <- c("FV-C01", "FV-Foce Volturno", "FV-03", "FV-06", "FV-10", "FV-Capua")

sub <- bind_rows(
  make_long(jan_raw, "Jan-26") |> filter(Station %in% keep_jan),
  make_long(jul_raw, "Jul-25") |> filter(Station %in% keep_jul)
)

p_raw <- ggplot(sub, aes(x = Value, y = Depth, colour = Station, group = Station)) +
  geom_path(alpha = 0.85, linewidth = 0.45) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_y_reverse() +
  facet_grid(Campaign ~ Variable, scales = "free_x") +
  labs(x = NULL, y = "Depth (m)") +
  theme_bw(base_size = 9) +
  theme(legend.position  = "bottom",
        legend.title     = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

save_fig(p_raw, "fig_raw_profiles.pdf", 7.2, 5.0)

## ─────────────────────────────────────────────────────────────────────────────
## Fig 5 — Cross-month posterior predictive (fig_cross_month.pdf)
## ─────────────────────────────────────────────────────────────────────────────

message("Generating fig_cross_month.pdf ...")

pred_file <- file.path(RES_DIR, "predictions.RData")
if (!file.exists(pred_file)) pred_file <- file.path(DATA_DIR, "predictions.RData")
if (!file.exists(pred_file))
  stop("predictions.RData not found. Run 02_fix_predictions.R first.")
load(pred_file)   # loads: predictions, perf_df_full

build_obs <- function(df){
  bins <- build_bins(df); agg <- aggregate_to_bins(df, bins)
  agg |> mutate(depth = depth_mid) |> select(depth, salinita_mean, salinita_se)
}
obs_jan <- build_obs(jan_raw)
obs_jul <- build_obs(jul_raw)

post_summary <- function(P, obs, dir){
  q <- apply(P, 2, quantile, probs = c(.025, .5, .975))
  data.frame(depth = obs$depth, obs = obs$salinita_mean,
             med = q[2,], lo = q[1,], hi = q[3,], dir = dir)
}

pred_df <- bind_rows(
  post_summary(predictions$full_j2j, obs_jan, "Jul-25 → Jan-26 (Mᵂᵈᵁᴸ)"),
  post_summary(predictions$full_j2u, obs_jul, "Jan-26 → Jul-25 (Mᵂᵈᵁᴸ)")
)
# Use plain ASCII for reliability
pred_df$dir <- factor(pred_df$dir,
  levels = c("Jul-25 → Jan-26 (Mᵂᵈᵁᴸ)",
             "Jan-26 → Jul-25 (Mᵂᵈᵁᴸ)"))

pred_df2 <- bind_rows(
  post_summary(predictions$full_j2j, obs_jan,
               "Jul-25 to Jan-26 (M_full)"),
  post_summary(predictions$full_j2u, obs_jul,
               "Jan-26 to Jul-25 (M_full)")
)

p_pred <- ggplot(pred_df2, aes(x = depth)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = med), colour = "steelblue3", linewidth = 0.7) +
  geom_point(aes(y = obs), colour = "black", size = 1.2, alpha = 0.8) +
  facet_wrap(~dir, scales = "free", ncol = 2) +
  labs(x = "Depth (m)", y = "Salinity (PSU)",
       caption = paste("Black: observed bin means.",
                       "Blue line: posterior median.",
                       "Band: 95% posterior predictive interval.")) +
  theme_bw(base_size = 9) +
  theme(plot.caption    = element_text(hjust = 0, size = 7),
        strip.background = element_rect(fill = "grey92", colour = NA))

save_fig(p_pred, "fig_cross_month.pdf", 7.2, 3.5)

## ─────────────────────────────────────────────────────────────────────────────
## Fig 4 — GP residual depth process (fig_gp_component.pdf)
## ─────────────────────────────────────────────────────────────────────────────

message("Generating fig_gp_component.pdf ...")

samp_file <- file.path(RES_DIR, "mfull_posterior_samples.rds")
if (!file.exists(samp_file))
  stop("mfull_posterior_samples.rds not found. Run 03_refit_save_samples.R first.")
samp <- readRDS(samp_file)

gp_panel <- function(s, lbl){
  q  <- apply(s$u, 2, quantile, c(.025, .5, .975))
  df <- data.frame(depth = s$depth, med = q[2,], lo = q[1,], hi = q[3,])
  ggplot(df, aes(x = depth)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkorange", alpha = 0.25) +
    geom_line(aes(y = med), colour = "darkorange3", linewidth = 0.7) +
    labs(x = "Depth (m)",
         y = expression(u[s]~"(posterior, standardised)"),
         subtitle = lbl) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA))
}

p_gp <- gp_panel(samp$jan, "Jan-26: GP residual (near-zero amplitude)") |
        gp_panel(samp$jul, "Jul-25: GP residual (substantial amplitude)")

save_fig(p_gp, "fig_gp_component.pdf", 7.2, 3.2)

## ─────────────────────────────────────────────────────────────────────────────
## Fig 3 — Simulation recovery (fig_sim_recovery.pdf)
## ─────────────────────────────────────────────────────────────────────────────

message("Generating fig_sim_recovery.pdf ...")

sim_sum_file <- file.path(SIM_DIR, "sim_summary_full.rds")
if (!file.exists(sim_sum_file))
  stop("sim_summary_full.rds not found. Run Simulation/metrics.R first.")

S <- readRDS(sim_sum_file)

scen_lab <- c(
  S1_jan  = "S1\n(Jan-like)",
  S2_jul  = "S2\n(Jul-like)",
  S3_nogp = "S3\n(σu=0)",
  S4_miss = "S4\n(misspec.)"
)
mod_colours <- c(
  M_full   = "#1a1a1a",
  M_noEIV  = "#1f77b4",
  M_noGP   = "#2ca02c",
  M_base   = "#d62728"
)
mod_labels <- c(
  M_full   = expression(M[full]),
  M_noEIV  = expression(M[noEIV]),
  M_noGP   = expression(M[noGP]),
  M_base   = expression(M[base])
)
cov_labels  <- c(temperatura = "Temperature", ossigeno_mgkg = expression(O[2]))

Bf <- S$B |>
  filter(cov %in% c("temperatura", "ossigeno_mgkg")) |>
  mutate(
    scen_f = factor(scen_lab[scenario], levels = scen_lab),
    mod_f  = factor(model, levels = names(mod_colours)),
    cov_f  = factor(cov,
               levels = c("temperatura", "ossigeno_mgkg"),
               labels = c("Temperature", "O[2]")),
    bias_lo = bias - rmse,
    bias_hi = bias + rmse
  )

p_bias <- ggplot(Bf, aes(x = scen_f, y = bias, colour = mod_f,
                          shape = mod_f, group = mod_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = bias_lo, ymax = bias_hi),
                position = position_dodge(width = 0.6), width = 0.2, linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  facet_wrap(~cov_f, ncol = 1, labeller = label_parsed) +
  scale_colour_manual(values = mod_colours, labels = mod_labels, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), labels = mod_labels, name = NULL) +
  labs(x = NULL, y = "Bias (± RMSE)",
       title = "Coefficient bias") +
  theme_bw(base_size = 8) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 8, face = "bold"))

p_cov <- ggplot(Bf, aes(x = scen_f, y = coverage, colour = mod_f,
                          shape = mod_f, group = mod_f)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.3, alpha = 0.5) +
  facet_wrap(~cov_f, ncol = 1, labeller = label_parsed) +
  scale_colour_manual(values = mod_colours, labels = mod_labels, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), labels = mod_labels, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "95% CI coverage",
       title = "Credible-interval coverage") +
  theme_bw(base_size = 8) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 8, face = "bold"))

p_sim <- p_bias + p_cov +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig(p_sim, "fig_sim_recovery.pdf", 7.2, 4.6)

message("\nAll paper figures generated.")
message("Saved to: ", RES_DIR)
message("Copied to: ", FIG_DIR)
