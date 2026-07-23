## =========================================================================
## config.R — Shared configuration for the thermal-extremes-CCE pipeline
## =========================================================================
## Sourced at the top of every analysis/figure script. Edit the paths in the
## USER SETTINGS block to point at your local copies of the input data
## (see DATA.md for where to obtain each dataset). Nothing else in this file
## needs to change.
##
## This repository contains CODE ONLY. The input data (SDM habitat-suitability
## output and ROMS/ESM physical projections) are archived externally and are
## described in DATA.md.
## =========================================================================

## ------------------------------------------------------------------------
## USER SETTINGS — edit these to match your machine
## ------------------------------------------------------------------------

## Root folder that holds the downloaded input data and receives outputs.
## Everything below is derived from this. Set it either by editing the path
## below (use an absolute path) or by exporting the TE_CCE_DATA_ROOT
## environment variable before running. See DATA.md for the expected layout.
data_root <- Sys.getenv(
  "TE_CCE_DATA_ROOT",
  unset = "/path/to/project-data-root"
)

## ------------------------------------------------------------------------
## DERIVED PATHS — normally no need to edit
## ------------------------------------------------------------------------
output_data_dir  <- file.path(data_root, "output data")   # NetCDF, CSVs, composite .grd
output_plots_dir <- file.path(data_root, "output plots")  # saved figures
shared_dir       <- file.path(data_root, "shared code")   # CDO mask + grid weights
models_dir       <- file.path(data_root, "EcoROMS models") # SDM BRT .rds objects

# Composite habitat-suitability rasters (built by 03_composite_sdm_output.R)
composites_dir   <- file.path(output_data_dir, "extremes_ccs_ssta", "composites")

# CDO-produced supporting files (see 01_preprocess_sst.sh)
mask_file         <- file.path(shared_dir, "lme-mask-out-cropped.nc")
grid_weights_file <- file.path(shared_dir, "gridweights.nc")

# California Current = region ID 3 in the LME mask
ccs_region_id <- 3

## ------------------------------------------------------------------------
## STUDY CONSTANTS
## ------------------------------------------------------------------------

# Earth System Models (ESMs) under RCP 8.5 from CMIP5
models <- c("GFDL", "Had", "IPSL")

# Target species (codes used throughout)
species_list <- c("casl", "hbwh", "lbst", "swor")

species_labels <- c(
  casl = "California sea lion",
  hbwh = "Humpback whale",
  lbst = "Leatherback sea turtle",
  swor = "Swordfish"
)

# Paths to the original Boosted Regression Tree (BRT) ensemble objects (.rds),
# used for partial dependence plots (Supp Fig 1) and variable contributions
# (Supp Table 1). See DATA.md.
sdm_model_paths <- list(
  casl = file.path(models_dir, "casl.res1.tc4.lr.1.10models.noLat.rds"),
  hbwh = file.path(models_dir, "HBWH.res1.tc3.lr01.10models.rds"),
  lbst = file.path(models_dir, "lbst_noSSH.res1.tc3.lr01.10models (1).rds"),
  swor = file.path(models_dir, "SWOR.res1.tc3.lr03.10models (1).rds")
)

# Analysis periods
historical_period <- as.Date(c("1991-01-01", "2020-12-31"))
future_period     <- as.Date(c("2071-01-01", "2100-12-31"))

## ------------------------------------------------------------------------
## Scenario codes and comparisons
## ------------------------------------------------------------------------
## Habitat-suitability composites are labelled by the months they average:
##   warm     == warm extreme  (90-100th percentile SST anomaly months)
##   cool     == cool extreme  (0-10th  percentile SST anomaly months)
##   average  == period mean (all months in the 30-year period)
## The three manuscript comparisons:
##   F[ave]-H[ave]    "Anthropogenic (long-term) warming"   -> Fig 2
##   H[warm]-H[cool]  "Historical thermal extremes"         -> Fig 3
##   F[warm]-F[cool]  "Future thermal extremes"             -> Fig 4
comp_levels <- c("F[ave]-H[ave]", "H[warm]-H[cool]", "F[warm]-F[cool]")

comp_labels_short <- c(
  "F[ave]-H[ave]" = "Anthropogenic warming",
  "H[warm]-H[cool]" = "Historical thermal extremes",
  "F[warm]-F[cool]" = "Future thermal extremes"
)

## ------------------------------------------------------------------------
## Composite filename tokens
## ------------------------------------------------------------------------
## The pipeline uses the canonical signal codes above (average/warm/cool)
## everywhere. But the composite .grd files on disk may embed a different
## token in their names. This map translates each canonical code to the token
## used in the FILENAMES, for file discovery only (find_composite_grd()).
##
## Set the RIGHT-HAND side to whatever token appears in your .grd filenames.
## Currently set for composites written by the original run, whose names use
## "amplify" (warm) and "minimize" (cool):
signal_file_tokens <- c(average = "average", warm = "amplify", cool = "minimize")
##
## If you (re)generate composites with 03_composite_sdm_output.R, that script
## writes "warm"/"cool", so use the identity mapping instead:
# signal_file_tokens <- c(average = "average", warm = "warm", cool = "cool")
