# Bayesian modelling of salinity profiles in the Volturno River mouth under depth dependence

**Agostino Gnasso, Lucia Rita Pacifico, Raffaele Mattera, Raffaele D'Adamo, Fabio Matano, Germana Scepi**

---

## Overview

This repository contains the Stan models, R code and simulation scripts used in the paper. The raw data are archived separately on Zenodo (see [Data](#data) below). The study proposes a Bayesian errors-in-variables (EIV) model augmented with a Gaussian Process (GP) component to estimate associations between salinity and physico-chemical covariates from depth-resolved CTD profiles in the Volturno River mouth (Campania, Italy).

Two synoptic field campaigns are analysed:

- **Jan-26** — winter high-discharge campaign (16 January 2026, 13 stations)
- **Jul-25** — summer low-discharge campaign (July 2025, 20 stations)

The full model (`M_full`) is compared against three reduced specifications obtained by removing the EIV layer, the GP layer, or both.

---

## Repository structure

```
Submission/
├── Analysis/                    # All reproducible code and data
│   ├── Analysis.Rproj           # Open this in RStudio before running scripts
│   ├── Data/                    # Raw input data (fetched from Zenodo — see below)
│   │   └── .gitkeep             # placeholder; xlsx files downloaded via download_data.R
│   ├── Stan/                    # Stan model files
│   │   ├── eiv_gp.stan          # M_full: EIV + GP
│   │   ├── eiv_gp_priors.stan   # M_full with tunable priors (prior sensitivity)
│   │   ├── no_eiv_gp.stan       # M_noEIV: GP only
│   │   ├── eiv_no_gp.stan       # M_noGP: EIV only
│   │   └── no_eiv_no_gp.stan    # M_base: neither EIV nor GP
│   ├── download_data.R          # Download raw data from Zenodo (run once)
│   ├── R/                       # Analysis scripts (run in numbered order)
│   │   ├── 00_data_config.R             # Zenodo DOI, file names, check_data()
│   │   ├── 00_extract_scale_factors.R   # Helper sourced by other scripts
│   │   ├── 01_run_analysis.R            # Fit all 4 models × 2 campaigns
│   │   ├── 01b_run_mcmc_only.R          # Alternative step 1 (MCMC only)
│   │   ├── 02_fix_predictions.R         # Cross-month predictions
│   │   ├── 03_refit_save_samples.R      # Posterior samples for figures
│   │   ├── 04_mcmc_refit_tight.R        # Sensitivity: tighter HMC
│   │   ├── 05_prior_sensitivity.R       # Sensitivity: prior grid
│   │   ├── 06_spatial_and_tau_sens.R    # Sensitivity: spatial + tau scaling
│   │   ├── 07_spatial_sens_redo.R       # Corrected spatial sensitivity
│   │   ├── 08_model_comparison.R        # Alternative comparison (adaptive bins)
│   │   ├── 09_make_figures.R            # Intermediate figure generation
│   │   ├── 10_regen_figures_psu.R       # Figures with PSU labels
│   │   ├── 11_generate_paper_figures.R  # All paper figures (main entry point)
│   │   └── 12_fit_traceplots.R          # MCMC traceplots (long run, ~hours)
│   ├── Simulation/              # Simulation study
│   │   ├── data_prep.R          # Data loading helpers
│   │   ├── dgp.R                # Data-generating process
│   │   ├── 00_extract_truth.R   # Extract ground truth from field fits
│   │   ├── run_sim.R            # Main simulation harness
│   │   ├── metrics.R            # Aggregate replicate results
│   │   └── tables_figures.R     # LaTeX tables + simulation figure
│   └── Results/                 # Generated outputs (figures, RDS files)
└── Journal_Submission/
    ├── Figures/                 # Paper-ready figures
    └── Manuscript/              # LaTeX source of the paper
```

---

## Requirements

**R** (≥ 4.2) with the following packages:

| Package | Purpose |
|---------|---------|
| `rstan` (≥ 2.21) | Bayesian inference via Stan |
| `bayesplot` | MCMC diagnostics and traceplots |
| `loo` | WAIC and LOO-CV |
| `scoringRules` | CRPS for predictive evaluation |
| `dplyr`, `tidyr` | Data wrangling |
| `ggplot2`, `patchwork` | Figures |
| `readxl` | Reading Excel data files |
| `mvtnorm`, `MASS`, `Matrix` | Matrix computations |
| `here` | Portable file paths |

Install all at once:

```r
install.packages(c(
  "rstan", "bayesplot", "loo", "scoringRules",
  "dplyr", "tidyr", "ggplot2", "patchwork",
  "readxl", "mvtnorm", "MASS", "Matrix", "here"
))
```

**Stan** is installed automatically with `rstan`. For best performance, follow the [RStan getting started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

---

## How to reproduce

All scripts use `here::here()` for path resolution. Open `Analysis/Analysis.Rproj` in RStudio before running any script — this sets the working directory automatically. Alternatively, set your working directory to `Analysis/` before running `Rscript`.

### Step 0 — Download the data

The raw `.xlsx` files are **not** in this repository. Download them from Zenodo before running anything else:

```r
# From the Analysis/ directory
source("download_data.R")
```

This fetches both campaign files into `Analysis/Data/`. If the DOI placeholder in `R/00_data_config.R` has not yet been updated (i.e. `ZENODO_DOI` still reads `"10.5281/zenodo.20447565"`), the script will print instructions for a manual download.

### Main analysis pipeline

Run scripts in numbered order. Each step produces outputs consumed by later steps.

```r
# Step 1 — Fit all 4 model specifications on both campaigns (~2–4 h)
source("R/01_run_analysis.R")

# Step 2 — Recompute cross-month predictions with robust GP solver
source("R/02_fix_predictions.R")

# Step 3 — Posterior samples for figure generation (~1 h)
source("R/03_refit_save_samples.R")

# Step 4 — Generate all paper figures (instant, reads existing outputs)
source("R/11_generate_paper_figures.R")
```

Figures are written to `Results/` and copied to `Journal_Submission/Figures/`.

### Sensitivity analyses (independent, any order after Step 1)

```r
source("R/04_mcmc_refit_tight.R")     # Tighter HMC: adapt_delta=0.99
source("R/05_prior_sensitivity.R")     # Prior grid on (s_beta, s_sigma_u)
source("R/06_spatial_and_tau_sens.R")  # Spatial extent + tau scaling
source("R/07_spatial_sens_redo.R")     # Corrected station IDs for spatial sens.
```

### MCMC traceplots (Appendix B)

```r
# Runs 2 x 4 chains x 24000 iterations — expect several hours
source("R/12_fit_traceplots.R")
```

### Simulation study

The simulation is designed to run in parallel shards (3 × 200 replicates). From the `Analysis/` directory:

```bash
# Launch 3 parallel shards (uses all results from 00_extract_truth.R)
bash Simulation/launch_full.sh

# After all shards complete, aggregate and generate the simulation figure
Rscript Simulation/metrics.R full
Rscript Simulation/tables_figures.R full
```

`00_extract_truth.R` must be run first to generate `truth_params.rds`:

```r
source("Simulation/00_extract_truth.R")
```

---

## Model specifications

| Model | EIV layer | GP layer | Description |
|-------|-----------|----------|-------------|
| `M_full` | ✓ | ✓ | Full model (proposed) |
| `M_noEIV` | ✗ | ✓ | GP without measurement-error correction |
| `M_noGP` | ✓ | ✗ | EIV regression without depth dependence |
| `M_base` | ✗ | ✗ | Weighted OLS on bin means |

All models are implemented in Stan (`Stan/`) using the exponential (Ornstein–Uhlenbeck) covariance kernel for the GP component.

---

## Data

The raw CTD profiles are **not** hosted in this repository. They are archived on Zenodo under the SAL.WE. (SAline WEdge intrusion) project (CNR):

> **Zenodo DOI:** [10.5281/zenodo.20447565](https://doi.org/10.5281/zenodo.20447565) *(DOI will be activated upon dataset publication)*

Two Excel files are provided in the Zenodo record:

| File | Campaign |
|------|---------|
| `...Campagna Gennaio 2026_01042026.xlsx` | January 2026 (Jan-26) |
| `...Campagna Luglio 2025_01042026.xlsx` | July 2025 (Jul-25) |

Each file contains one row per CTD reading with columns for depth (m), temperature (°C), salinity (PSU), dissolved oxygen (mg/kg), pH, chlorophyll-a (Trx-ChlA), phycocyanin, phycoerythrin and distance from the most seaward station (FV-C01).

**To update the DOI** once the Zenodo record is live, edit the single line in `R/00_data_config.R`:

```r
ZENODO_DOI <- "10.5281/zenodo.20447565"   # replace with actual DOI
```

Then run `download_data.R` to fetch the files automatically.

---

## Outputs

After running the full pipeline, `Results/` contains:

| File | Script | Description |
|------|--------|-------------|
| `results.RData` | `01_run_analysis.R` | All Stan fit objects + diagnostics |
| `predictions.RData` | `02_fix_predictions.R` | Cross-month predictive draws |
| `mfull_posterior_samples.rds` | `03_refit_save_samples.R` | Posterior samples for figures |
| `mcmc_tight_diag.rds` | `04_mcmc_refit_tight.R` | Tight-HMC diagnostics |
| `prior_sensitivity_jan.rds` | `05_prior_sensitivity.R` | Prior sensitivity posteriors |
| `spatial_sens.rds` | `07_spatial_sens_redo.R` | Spatial sensitivity posteriors |
| `tau_sens.rds` | `06_spatial_and_tau_sens.R` | Tau scaling sensitivity |
| `scale_factors.csv` | `00_extract_scale_factors.R` | Standardisation factors (Appendix A) |
| `fig_raw_profiles.pdf` | `11_generate_paper_figures.R` | Figure 2 |
| `fig_sim_recovery.pdf` | `11_generate_paper_figures.R` | Figure 3 |
| `fig_gp_component.pdf` | `11_generate_paper_figures.R` | Figure 4 |
| `fig_cross_month.pdf` | `11_generate_paper_figures.R` | Figure 5 |
| `fig_traceplots_full_jan.pdf` | `12_fit_traceplots.R` | Appendix B (Jan-26) |
| `fig_traceplots_full_jul.pdf` | `12_fit_traceplots.R` | Appendix B (Jul-25) |

---

## Citation

If you use this code or data, please cite:

> Gnasso A., Pacifico L.R., Mattera R., D'Adamo R., Matano F., Scepi G. (2026).
> *Bayesian modelling of salinity profiles in the Volturno River mouth under depth dependence. (WP)*

---

## License

The code is released under the [MIT License](https://opensource.org/licenses/MIT).  
The data are made available for reproducibility purposes under the terms of the SAL.WE. project (CNR).
