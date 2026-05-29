## 10_regen_figures_psu.R
## Regenerates fig_raw_profiles.pdf, fig_cross_month.pdf, fig_ppc.pdf
## with "PSU" labels instead of "‰". Run after 02 and 03.
##
## Run from the Analysis/ directory (open Analysis.Rproj first).

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(ggplot2)
  library(tidyr);  library(patchwork); library(here)
})

source(here::here("R", "00_data_config.R"))
DATA_DIR <- here::here("Data")
check_data(DATA_DIR)

vars_phys <- c("temperatura","ossigeno_mgkg","ph","trx_chla",
               "phycocyanin","phycoerythrin","dist")

safe_se <- function(x){ n <- sum(!is.na(x)); if(n<=1) NA_real_ else sd(x,na.rm=TRUE)/sqrt(n) }

load_jan <- function(){
  f <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx")
  read_excel(f) |>
    rename(profondita=`Profondità (m)`,
           temperatura=`temperatura (C°)`,
           salinita=`Salinità (‰)`,
           ossigeno_mgkg=`Oxygen (mg/Kg)`,
           ph=`pH`,
           trx_chla=`Trx-chl(a)`,
           phycocyanin=`Phycocyanin`,
           phycoerythrin=`Phycoerythrin`,
           dist=`Distanza dalla foce (km)`) |>
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10) |>
    mutate(dist = dist - min(dist, na.rm=TRUE))
}

load_jul <- function(){
  f <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx")
  read_excel(f) |>
    rename(profondita=`Profondità (m)`,
           temperatura=`temperatura (°C)`,
           salinita=`Salinità (‰)`,
           ossigeno_mgkg=`Oxygen (mg/Kg)`,
           ph=`pH`,
           trx_chla=`Trx-chl(a)`,
           phycocyanin=`Phycocyanin`,
           phycoerythrin=`Phycoerythrin`,
           dist=`Distanza dalla Foce(km)`) |>
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10) |>
    mutate(dist = dist - min(dist, na.rm=TRUE))
}

build_bins <- function(data, depth_col="profondita", id_col="ID",
                       min_depth=0.01, max_depth=10,
                       base_step=0.01, min_sonde=3){
  df <- data |> filter(!is.na(.data[[depth_col]]),
                       .data[[depth_col]] >= min_depth,
                       .data[[depth_col]] <= max_depth) |>
    mutate(base_bin = floor(.data[[depth_col]] / base_step) * base_step)
  bc <- df |> group_by(base_bin) |>
    summarise(sonde   = list(sort(unique(.data[[id_col]]))),
              n_sonde = n_distinct(.data[[id_col]]), .groups="drop") |>
    arrange(base_bin)
  bins <- list(); i <- 1L; k <- 1L
  while(i <= nrow(bc)){
    le <- bc$base_bin[i]; cur <- bc$sonde[[i]]; j <- i
    while(length(cur) < min_sonde && j < nrow(bc)){ j <- j+1L; cur <- union(cur, bc$sonde[[j]]) }
    if(length(cur) >= min_sonde){
      re <- bc$base_bin[j] + base_step
      bins[[k]] <- data.frame(bin_id=k, left=le, right=re,
                              mid=(le+re)/2, width=re-le, n_sonde=length(cur))
      k <- k+1L; i <- j+1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col="profondita"){
  out <- data; out$bin_id <- out$depth_mid <- NA_real_
  for(r in seq_len(nrow(bins))){
    idx <- which(out[[depth_col]] >= bins$left[r] & out[[depth_col]] < bins$right[r])
    out$bin_id[idx]  <- bins$bin_id[r]
    out$depth_mid[idx] <- bins$mid[r]
  }
  out
}

aggregate_to_bins <- function(df, bins){
  d <- assign_bins(df, bins)
  d |> filter(!is.na(bin_id)) |>
    group_by(bin_id, depth_mid) |>
    summarise(across(c(salinita, all_of(vars_phys)),
                     list(mean=~mean(.x,na.rm=TRUE), se=~safe_se(.x))),
              .groups="drop") |>
    arrange(depth_mid)
}

# ── Fig 2 : raw profiles ──────────────────────────────────────────────────────

jan_raw <- load_jan()
jul_raw <- load_jul()

make_raw <- function(df, label){
  df |>
    select(ID, profondita, salinita, temperatura, ossigeno_mgkg) |>
    rename(Salinity=salinita, Temperature=temperatura, O2=ossigeno_mgkg,
           Depth=profondita, Station=ID) |>
    pivot_longer(c(Salinity, Temperature, O2),
                 names_to="Variable", values_to="Value") |>
    mutate(Campaign=label)
}

both <- bind_rows(make_raw(jan_raw,"Jan-26"), make_raw(jul_raw,"Jul-25")) |>
  mutate(Variable = factor(Variable,
                           levels = c("Salinity","Temperature","O2"),
                           labels = c("Salinity (PSU)","Temperature (deg C)","O2 (mg/kg)")))

keep_jan <- c("FV-C01","FV-Foce Volturno","FV-03","FV-06","FV-10")
keep_jul <- c("FV-C01","FV-Foce Volturno","FV-03","FV-06","FV-10","FV-Capua")
sub <- both |> filter((Campaign=="Jan-26" & Station %in% keep_jan) |
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

# ── Fig 5 : cross-month predictions ──────────────────────────────────────────

obs_jan <- aggregate_to_bins(jan_raw, build_bins(jan_raw)) |>
  mutate(depth=depth_mid) |> select(depth, salinita_mean, salinita_se)
obs_jul <- aggregate_to_bins(jul_raw, build_bins(jul_raw)) |>
  mutate(depth=depth_mid) |> select(depth, salinita_mean, salinita_se)

load(here::here("Results","predictions.RData"))

post_summary <- function(P, obs, dir){
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

# ── Fig 3 : posterior predictive checks ──────────────────────────────────────

samp <- readRDS(here::here("Results","mfull_posterior_samples.rds"))

ppc_panel <- function(s, lbl){
  set.seed(42)
  yrep <- as.data.frame(s$y_rep[sample(nrow(s$y_rep), 50), ]) |>
    setNames(s$depth) |>
    mutate(.draw = row_number()) |>
    pivot_longer(-.draw, names_to="depth", values_to="y_rep") |>
    mutate(depth = as.numeric(depth))
  obs <- data.frame(depth=s$depth, y=s$y_obs)
  ggplot() +
    geom_line(data=yrep, aes(x=depth, y=y_rep, group=.draw),
              alpha=0.08, colour="steelblue") +
    geom_point(data=obs, aes(x=depth, y=y), size=1) +
    labs(x="Depth (m)", y="Salinity (PSU)", subtitle=lbl) +
    theme_bw(base_size=9)
}

p_ppc <- ppc_panel(samp$jan, "Jan-26 PPC") | ppc_panel(samp$jul, "Jul-25 PPC")
ggsave(here::here("Results","fig_ppc.pdf"), p_ppc,
       width=7.2, height=3.2, device=cairo_pdf)
