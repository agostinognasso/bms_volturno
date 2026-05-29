## download_data.R
##
## Downloads the raw campaign data from Zenodo into Analysis/Data/.
## Run this script ONCE before running any analysis script.
##
## Usage (from the Analysis/ directory):
##   Rscript download_data.R
##   # or in RStudio: source("download_data.R")
##
## Requirements: an internet connection.
## The data will be saved to Analysis/Data/.

suppressPackageStartupMessages(library(here))

source(here::here("R", "00_data_config.R"))

DATA_DIR <- here::here("Data")
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR)

## ── Update ZENODO_DOI in R/00_data_config.R before running ───────────────────
if (grepl("XXXXXXX", ZENODO_DOI)) {
  stop(
    "\nZenodo DOI not yet set.\n\n",
    "  Open R/00_data_config.R and replace the placeholder:\n",
    "    ZENODO_DOI <- \"10.5281/zenodo.XXXXXXX\"\n",
    "  with the actual DOI once the dataset is published on Zenodo.\n\n",
    "  Alternatively, download the files manually from Zenodo and\n",
    "  place them in Analysis/Data/.\n",
    call. = FALSE
  )
}

## Zenodo direct-download URL pattern
zenodo_record <- sub("10.5281/zenodo.", "", ZENODO_DOI)
base_url <- paste0("https://zenodo.org/records/", zenodo_record, "/files/")

files <- c(DATA_FILE_JAN, DATA_FILE_JUL)

for (fname in files) {
  dest <- file.path(DATA_DIR, fname)
  if (file.exists(dest)) {
    message("Already present, skipping: ", fname)
    next
  }
  url <- paste0(base_url, utils::URLencode(fname))
  message("Downloading: ", fname, " ...")
  tryCatch(
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
    error = function(e) {
      message("\nDownload failed for: ", fname)
      message("URL tried: ", url)
      message("Please download manually from:\n  ", ZENODO_URL)
      message("and save to:\n  ", DATA_DIR)
    }
  )
}

## Verify
result <- tryCatch(check_data(DATA_DIR), error = function(e) e)
if (inherits(result, "error")) {
  message("\nSome files could not be downloaded automatically.")
  message("Please fetch them manually from:\n  ", ZENODO_URL)
  message("and place them in:\n  ", DATA_DIR)
} else {
  message("\nAll data files are present in Data/. You can now run the analysis.")
}
