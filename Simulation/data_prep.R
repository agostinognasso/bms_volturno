## ============================================================
## sim/data_prep.R — self-contained data loading & binning
## Copied verbatim (logic) from R/run_analysis.R so the simulation
## study is reproducible without re-running the full model pipeline.
## Original run_analysis.R is left untouched (Minimal Impact).
## ============================================================

suppressMessages({
  library(dplyr)
  library(readxl)
  library(here)
})

source(here::here("R", "00_data_config.R"))
DATA_DIR <- here::here("Data")
check_data(DATA_DIR)

safe_se <- function(x) {
  n <- sum(!is.na(x))
  if (n <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n)
}

build_bins <- function(data, depth_col = "profondita", id_col = "ID",
                       min_depth = 0.01, max_depth = 10,
                       base_step = 0.01, min_sonde = 3) {
  df <- data %>%
    filter(!is.na(.data[[depth_col]]),
           .data[[depth_col]] >= min_depth,
           .data[[depth_col]] <= max_depth) %>%
    mutate(base_bin = floor(.data[[depth_col]] / base_step) * base_step)

  base_counts <- df %>%
    group_by(base_bin) %>%
    summarise(sonde   = list(sort(unique(.data[[id_col]]))),
              n_sonde = n_distinct(.data[[id_col]]), .groups = "drop") %>%
    arrange(base_bin)

  bins <- list(); i <- 1L; k <- 1L
  while (i <= nrow(base_counts)) {
    left_edge <- base_counts$base_bin[i]
    cur       <- base_counts$sonde[[i]]
    j         <- i
    while (length(cur) < min_sonde && j < nrow(base_counts)) {
      j <- j + 1L; cur <- union(cur, base_counts$sonde[[j]])
    }
    if (length(cur) >= min_sonde) {
      right_edge <- base_counts$base_bin[j] + base_step
      bins[[k]]  <- data.frame(bin_id = k, left = left_edge,
                                right = right_edge,
                                mid   = (left_edge + right_edge) / 2,
                                width = right_edge - left_edge,
                                n_sonde = length(cur))
      k <- k + 1L; i <- j + 1L
    } else break
  }
  bind_rows(bins)
}

assign_bins <- function(data, bins, depth_col = "profondita") {
  out <- data
  out$bin_id <- out$depth_mid <- NA_real_
  for (r in seq_len(nrow(bins))) {
    idx <- which(out[[depth_col]] >= bins$left[r] & out[[depth_col]] < bins$right[r])
    out$bin_id[idx]    <- bins$bin_id[r]
    out$depth_mid[idx] <- bins$mid[r]
  }
  out
}

load_jan <- function() {
  f  <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx")
  d  <- read_excel(f)
  df <- d %>%
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (C°)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph            = `pH`,
           trx_chla      = `Trx-chl(a)`,
           phycocyanin   = `Phycocyanin`,
           phycoerythrin = `Phycoerythrin`,
           dist          = `Distanza dalla foce (km)`) %>%
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10)
  df$dist <- df$dist - min(df$dist, na.rm = TRUE)
  df
}

load_jul <- function() {
  f  <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx")
  d  <- read_excel(f)
  df <- d %>%
    rename(profondita    = `Profondità (m)`,
           temperatura   = `temperatura (°C)`,
           salinita      = `Salinità (‰)`,
           ossigeno_mgkg = `Oxygen (mg/Kg)`,
           ph            = `pH`,
           trx_chla      = `Trx-chl(a)`,
           phycocyanin   = `Phycocyanin`,
           phycoerythrin = `Phycoerythrin`,
           dist          = `Distanza dalla Foce(km)`) %>%
    filter(!is.na(profondita), profondita >= 0.1, profondita <= 10)
  df$dist <- df$dist - min(df$dist, na.rm = TRUE)
  df
}

vars_phys <- c("temperatura","salinita","ossigeno_mgkg","ph",
               "trx_chla","phycocyanin","phycoerythrin","dist")
cov_vars  <- setdiff(vars_phys, "salinita")   # 7 covariates

aggregate_to_bins <- function(df) {
  bins   <- build_bins(df)
  df_bin <- assign_bins(df, bins)

  by_probe <- df_bin %>%
    filter(!is.na(bin_id)) %>%
    group_by(ID, bin_id, depth_mid) %>%
    summarise(across(all_of(vars_phys), ~ mean(.x, na.rm=TRUE)), .groups="drop")

  by_probe %>%
    group_by(bin_id, depth_mid) %>%
    summarise(
      n_sonde = n(),
      across(all_of(vars_phys),
             list(mean=~mean(.x, na.rm=TRUE), se=~safe_se(.x)),
             .names="{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(depth_mid)
}

build_stan_data <- function(agg) {
  cov_mean <- paste0(cov_vars, "_mean")
  cov_se   <- paste0(cov_vars, "_se")

  dat <- agg %>%
    select(depth_mid, salinita_mean, salinita_se,
           all_of(cov_mean), all_of(cov_se)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  y_raw <- dat$salinita_mean
  sy    <- sd(y_raw); my <- mean(y_raw)
  y_obs <- (y_raw - my) / sy
  tau_y <- pmax(dat$salinita_se / sy, 1e-6)

  X_raw  <- as.matrix(dat[, cov_mean])
  TauX   <- as.matrix(dat[, cov_se])
  mx     <- colMeans(X_raw); sx <- apply(X_raw, 2, sd)
  X_obs  <- sweep(sweep(X_raw, 2, mx, "-"), 2, sx, "/")
  tau_x  <- pmax(sweep(TauX, 2, sx, "/"), 1e-6)

  d_raw  <- dat$depth_mid
  d      <- (d_raw - min(d_raw)) / (max(d_raw) - min(d_raw))

  list(
    stan = list(S=nrow(X_obs), K=ncol(X_obs),
                y_obs=as.vector(y_obs), X_obs=X_obs,
                tau_y=as.vector(tau_y), tau_x=tau_x,
                d=as.vector(d)),
    meta = list(my=my, sy=sy, mx=mx, sx=sx, d_raw=d_raw,
                y_raw=y_raw, cov_names=cov_vars,
                n_bins=nrow(dat))
  )
}

## Convenience: build both campaigns' standardised Stan data.
load_campaign_stan_data <- function() {
  list(
    jan = build_stan_data(aggregate_to_bins(load_jan())),
    jul = build_stan_data(aggregate_to_bins(load_jul()))
  )
}
