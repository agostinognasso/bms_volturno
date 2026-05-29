## 09_make_figures.R  --  Generate all manuscript figures.
##
## Requires:
##   Results/predictions.RData      (from 02_fix_predictions.R)
##   Results/mfull_posterior_samples.rds  (from 03_refit_save_samples.R)
##
## Outputs (all in Results/):
##   fig_raw_profiles.pdf    Fig 2: raw vertical profiles
##   fig_cross_month.pdf     Fig 5: cross-month observed vs predicted
##   fig_ppc.pdf             Fig 3: posterior predictive checks (if samples available)
##   fig_gp_component.pdf    Fig 4: GP residual u_s (if samples available)
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(ggplot2); library(tidyr); library(patchwork); library(here)
})

source(here::here("R", "00_extract_scale_factors.R"), local = TRUE)
## ^ provides load_jan(), load_jul(), build_bins(), aggregate_to_bins(), vars_phys

# ---------------------------------------------------------------- Fig 2 ----
make_raw <- function(df, label) {
  df %>% filter(profondita >= 0.1, profondita <= 10) %>%
    select(ID, profondita, salinita, temperatura, ossigeno_mgkg) %>%
    rename(Salinity=salinita, Temperature=temperatura, `O2`=ossigeno_mgkg,
           Depth=profondita, Station=ID) %>%
    pivot_longer(c(Salinity, Temperature, `O2`), names_to="Variable", values_to="Value") %>%
    mutate(Campaign=label)
}

jan_long <- make_raw(load_jan(), "Jan-26")
jul_long <- make_raw(load_jul(), "Jul-25")
both <- bind_rows(jan_long, jul_long) %>%
  mutate(Variable=factor(Variable, levels=c("Salinity","Temperature","O2"),
                         labels=c("Salinity (PSU)","Temperature (deg C)","O2 (mg/kg)")))

keep_jan <- c("FV-C01","FV-Foce Volturno","FV-03","FV-06","FV-10")
keep_jul <- c("FV-C01","FV-Foce Volturno","FV-03","FV-06","FV-10","FV-Capua")
sub <- both %>% filter((Campaign=="Jan-26" & Station %in% keep_jan) |
                       (Campaign=="Jul-25" & Station %in% keep_jul))

p_raw <- ggplot(sub, aes(x=Value, y=Depth, colour=Station, group=Station)) +
  geom_path(alpha=0.85, linewidth=0.45) + geom_point(size=0.6, alpha=0.7) +
  scale_y_reverse() +
  facet_grid(Campaign ~ Variable, scales="free_x") +
  labs(x=NULL, y="Depth (m)") +
  theme_bw(base_size=9) +
  theme(legend.position="bottom", legend.title=element_blank(),
        strip.background=element_rect(fill="grey92", colour=NA),
        panel.grid.minor=element_blank())
ggsave(here::here("Results","fig_raw_profiles.pdf"), p_raw,
       width=7.2, height=5.0, device=cairo_pdf)

# ---------------------------------------------------------------- Fig 5 ----
load(here::here("Results", "predictions.RData"))

build_obs <- function(df) {
  bins <- build_bins(df); agg <- aggregate_to_bins(df, bins)
  agg %>% mutate(depth=depth_mid) %>% select(depth, salinita_mean, salinita_se)
}
obs_jan <- build_obs(load_jan())
obs_jul <- build_obs(load_jul())

post_summary <- function(P, obs, dir) {
  q <- apply(P, 2, quantile, probs=c(.025,.5,.975))
  data.frame(depth=obs$depth, obs=obs$salinita_mean,
             med=q[2,], lo=q[1,], hi=q[3,], dir=dir)
}
pf_j2j <- post_summary(predictions$full_j2j, obs_jan, "Jul-25 -> Jan-26 (M_full)")
pf_j2u <- post_summary(predictions$full_j2u, obs_jul, "Jan-26 -> Jul-25 (M_full)")
pred_df <- bind_rows(pf_j2j, pf_j2u)

p_pred <- ggplot(pred_df, aes(x=depth)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="steelblue", alpha=0.25) +
  geom_line(aes(y=med), colour="steelblue3", linewidth=0.7) +
  geom_point(aes(y=obs), colour="black", size=1.2, alpha=0.8) +
  facet_wrap(~dir, scales="free", ncol=2) +
  labs(x="Depth (m)", y="Salinity (PSU)",
       caption="Black: observed bin means. Blue line: posterior median. Band: 95% posterior predictive interval.") +
  theme_bw(base_size=9) +
  theme(plot.caption=element_text(hjust=0, size=7),
        strip.background=element_rect(fill="grey92", colour=NA))
ggsave(here::here("Results","fig_cross_month.pdf"), p_pred,
       width=7.2, height=3.5, device=cairo_pdf)

# ---------------------------------------- Figs 3-4 (if samples available) --
samp_file <- here::here("Results", "mfull_posterior_samples.rds")
if (file.exists(samp_file)) {
  samp <- readRDS(samp_file)

  ppc_panel <- function(s, lbl){
    yrep <- as.data.frame(s$y_rep[sample(nrow(s$y_rep), 50),]) |> setNames(s$depth) |>
      mutate(.draw=row_number()) |> pivot_longer(-.draw, names_to="depth", values_to="y_rep") |>
      mutate(depth=as.numeric(depth))
    obs <- data.frame(depth=s$depth, y=s$y_obs)
    ggplot() +
      geom_line(data=yrep, aes(x=depth, y=y_rep, group=.draw), alpha=0.08, colour="steelblue") +
      geom_point(data=obs, aes(x=depth, y=y), size=1) +
      labs(x="Depth (m)", y="Salinity (PSU)", subtitle=lbl) +
      theme_bw(base_size=9)
  }
  p_ppc <- ppc_panel(samp$jan, "Jan-26 PPC") | ppc_panel(samp$jul, "Jul-25 PPC")
  ggsave(here::here("Results","fig_ppc.pdf"), p_ppc, width=7.2, height=3.2, device=cairo_pdf)

  gp_panel <- function(s, lbl){
    q <- apply(s$u, 2, quantile, c(.025,.5,.975))
    df <- data.frame(depth=s$depth, med=q[2,], lo=q[1,], hi=q[3,])
    ggplot(df, aes(x=depth)) +
      geom_hline(yintercept=0, colour="grey60", linewidth=0.3) +
      geom_ribbon(aes(ymin=lo, ymax=hi), fill="darkorange", alpha=0.25) +
      geom_line(aes(y=med), colour="darkorange3", linewidth=0.7) +
      labs(x="Depth (m)", y=expression(u[s]~"(posterior, standardised)"), subtitle=lbl) +
      theme_bw(base_size=9)
  }
  p_gp <- gp_panel(samp$jan, "Jan-26: GP residual") | gp_panel(samp$jul, "Jul-25: GP residual")
  ggsave(here::here("Results","fig_gp_component.pdf"), p_gp, width=7.2, height=3.2, device=cairo_pdf)
}
