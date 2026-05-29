## 00_data_config.R
##
## Central configuration for raw data files.
## Sourced automatically by every script that loads campaign data.
##
## ── IMPORTANT ─────────────────────────────────────────────────────────────────
## The two raw .xlsx files are NOT included in the GitHub repository.
## They are archived on Zenodo. Update ZENODO_DOI below once the dataset
## is published, then run download_data.R to fetch them automatically.
## ─────────────────────────────────────────────────────────────────────────────

## Zenodo DOI — update this when the dataset record is published
ZENODO_DOI <- "10.5281/zenodo.XXXXXXX"
ZENODO_URL <- paste0("https://doi.org/", ZENODO_DOI)

## Canonical file names (do not change; must match the Zenodo deposit)
DATA_FILE_JAN <- paste0(
  "Progetto SALWE-CNR_profili di salinita fiume Volturno",
  "_Campagna Gennaio 2026_01042026.xlsx"
)
DATA_FILE_JUL <- paste0(
  "Progetto SALWE-CNR_profili di salinita fiume Volturno",
  "_Campagna Luglio 2025_01042026.xlsx"
)

## check_data(data_dir) ─────────────────────────────────────────────────────────
## Call this at the top of any script that reads campaign data.
## Stops with a clear Zenodo message if files are missing.
check_data <- function(data_dir) {
  jan_path <- file.path(data_dir, DATA_FILE_JAN)
  jul_path <- file.path(data_dir, DATA_FILE_JUL)

  ## Also accept the original Italian filenames (with accent)
  jan_alt <- file.path(data_dir, sub("salinita", "salinità", DATA_FILE_JAN))
  jul_alt <- file.path(data_dir, sub("salinita", "salinità", DATA_FILE_JUL))

  jan_ok <- file.exists(jan_path) || file.exists(jan_alt)
  jul_ok <- file.exists(jul_path) || file.exists(jul_alt)

  if (!jan_ok || !jul_ok) {
    missing <- c(
      if (!jan_ok) DATA_FILE_JAN,
      if (!jul_ok) DATA_FILE_JUL
    )
    stop(
      "\nData file(s) not found in: ", data_dir, "\n\n",
      "  Missing:\n",
      paste0("    - ", missing, collapse = "\n"), "\n\n",
      "  The raw data are archived on Zenodo:\n",
      "    ", ZENODO_URL, "\n\n",
      "  Run Analysis/download_data.R to fetch them automatically,\n",
      "  or download manually and place the .xlsx files in Analysis/Data/.\n",
      call. = FALSE
    )
  }

  ## Return actual paths (with or without accent)
  invisible(list(
    jan = if (file.exists(jan_path)) jan_path else jan_alt,
    jul = if (file.exists(jul_path)) jul_path else jul_alt
  ))
}
