## extract_scale_factors.R
## Compute mean and sd (on the bin-averaged scale) for the response
## (salinity) and the seven covariates in each of the two campaigns.
## These factors are needed to back-convert the standardised regression
## coefficients to physical units (M14 of the review).
##
## Output: scale_factors.csv (long format) and a LaTeX snippet
## scale_factors.tex with the table body for the manuscript appendix.
##
## This script is also sourced (local = TRUE) by several other scripts
## to provide load_jan(), load_jul(), build_bins(), aggregate_to_bins(),
## and vars_phys without re-running the CSV/TeX export.

suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(here)
})

source(here::here("R", "00_data_config.R"))

DATA_DIR <- here::here("Data")
OUT_DIR  <- here::here("Results")
check_data(DATA_DIR)

vars_phys <- c("temperatura", "ossigeno_mgkg", "ph", "trx_chla",
               "phycocyanin", "phycoerythrin", "dist")

safe_se <- function(x){ n <- sum(!is.na(x)); if (n<=1) NA_real_ else sd(x,na.rm=TRUE)/sqrt(n) }

build_bins <- function(data, depth_col="profondita", id_col="ID",
                       min_depth=0.01, max_depth=10,
                       base_step=0.01, min_sonde=3){
  df <- data %>% filter(!is.na(.data[[depth_col]]),
                        .data[[depth_col]]>=min_depth,
                        .data[[depth_col]]<=max_depth) %>%
    mutate(base_bin=floor(.data[[depth_col]]/base_step)*base_step)
  bc <- df %>% group_by(base_bin) %>%
    summarise(sonde=list(sort(unique(.data[[id_col]]))),
              n_sonde=n_distinct(.data[[id_col]]),.groups="drop") %>%
    arrange(base_bin)
  bins <- list(); i<-1L; k<-1L
  while(i<=nrow(bc)){
    le <- bc$base_bin[i]; cur <- bc$sonde[[i]]; j <- i
    while(length(cur)<min_sonde && j<nrow(bc)){ j<-j+1L; cur<-union(cur,bc$sonde[[j]]) }
    if(length(cur)>=min_sonde){
      re <- bc$base_bin[j]+base_step
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
    idx <- which(out[[depth_col]]>=bins$left[r] & out[[depth_col]]<bins$right[r])
    out$bin_id[idx]<-bins$bin_id[r]; out$depth_mid[idx]<-bins$mid[r]
  }
  out
}

load_jan <- function(){
  f <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Gennaio 2026_01042026.xlsx")
  d <- read_excel(f)
  d %>% rename(profondita=`Profondità (m)`,
               temperatura=`temperatura (C°)`,
               salinita=`Salinità (‰)`,
               ossigeno_mgkg=`Oxygen (mg/Kg)`,
               ph=`pH`,
               trx_chla=`Trx-chl(a)`,
               phycocyanin=`Phycocyanin`,
               phycoerythrin=`Phycoerythrin`,
               dist=`Distanza dalla foce (km)`) %>%
    filter(!is.na(profondita), profondita>=0.1, profondita<=10) %>%
    mutate(dist=dist-min(dist,na.rm=TRUE))
}

load_jul <- function(){
  f <- file.path(DATA_DIR,
    "Progetto SALWE-CNR_profili di salinità fiume Volturno_Campagna Luglio 2025_01042026.xlsx")
  d <- read_excel(f)
  d %>% rename(profondita=`Profondità (m)`,
               temperatura=`temperatura (°C)`,
               salinita=`Salinità (‰)`,
               ossigeno_mgkg=`Oxygen (mg/Kg)`,
               ph=`pH`,
               trx_chla=`Trx-chl(a)`,
               phycocyanin=`Phycocyanin`,
               phycoerythrin=`Phycoerythrin`,
               dist=`Distanza dalla Foce(km)`) %>%
    filter(!is.na(profondita), profondita>=0.1, profondita<=10) %>%
    mutate(dist=dist-min(dist,na.rm=TRUE))
}

aggregate_to_bins <- function(df, bins){
  d <- assign_bins(df, bins)
  d %>% filter(!is.na(bin_id)) %>%
    group_by(bin_id, depth_mid) %>%
    summarise(across(c(salinita, all_of(vars_phys)),
                     list(mean=~mean(.x,na.rm=TRUE), se=~safe_se(.x))),
              .groups="drop") %>% arrange(depth_mid)
}

extract_factors <- function(df, label){
  bins <- build_bins(df)
  agg  <- aggregate_to_bins(df, bins)
  y_mean <- mean(agg$salinita_mean, na.rm=TRUE)
  y_sd   <- sd(agg$salinita_mean,   na.rm=TRUE)
  res <- data.frame(campaign=label, variable="salinita",
                    mu=y_mean, sigma=y_sd)
  for(v in vars_phys){
    col <- paste0(v,"_mean")
    res <- rbind(res, data.frame(campaign=label, variable=v,
                                 mu=mean(agg[[col]],na.rm=TRUE),
                                 sigma=sd(agg[[col]],na.rm=TRUE)))
  }
  res
}


jan <- load_jan()
jul <- load_jul()

sf <- rbind(extract_factors(jan,"Jan-26"),
            extract_factors(jul,"Jul-25"))
write.csv(sf, file.path(OUT_DIR,"scale_factors.csv"), row.names=FALSE)

## Format LaTeX table body
pretty <- c(salinita="Salinity (\\textperthousand)",
            temperatura="Temperature ($^{\\circ}$C)",
            ossigeno_mgkg="O$_2$ (mg~kg$^{-1}$)",
            ph="pH",
            trx_chla="Trx-Chl-a",
            phycocyanin="Phycocyanin",
            phycoerythrin="Phycoerythrin",
            dist="Distance (km)")

j <- sf[sf$campaign=="Jan-26",]; u <- sf[sf$campaign=="Jul-25",]
rows <- character(nrow(j))
for(i in seq_len(nrow(j))){
  rows[i] <- sprintf("    %s & %.3f & %.3f & %.3f & %.3f \\\\",
                     pretty[j$variable[i]], j$mu[i], j$sigma[i],
                     u$mu[u$variable==j$variable[i]],
                     u$sigma[u$variable==j$variable[i]])
}
writeLines(rows, file.path(OUT_DIR,"scale_factors.tex"))
