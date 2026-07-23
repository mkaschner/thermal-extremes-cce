#!/usr/bin/env Rscript
## =========================================================================
## 03_composite_sdm_output.R — Composite daily SDM output by scenario
## =========================================================================
## RUN ON THE SERVER where the daily SDM habitat-suitability .grd files live.
## This is intentionally SELF-CONTAINED (it does not source config.R), because
## it runs in a different environment from the rest of the pipeline. Edit the
## CONFIGURATION block below to match the server's directory layout.
##
## Purpose: for each species x ESM x period x signal, average the daily SDM
## habitat-suitability rasters for the extreme-warming months identified in
## 02_identify_extremes.R, producing the composite .grd files consumed by
## 04_centroid_area_metrics.R and 05_habitat_maps.R.
##
## Inputs:
##   - Daily SDM rasters:  {MODEL_DIR}/{species}*mean.grd  (one per day)
##   - Extreme-dates CSV:  extreme_signal_dates_ccs_ssta_{MODEL}_{date}.csv
##                         (copied from 02's output)
## Outputs (to OUTPUT_DATA_DIR):
##   - {MODEL}_{species}_{period}_{signal}_mean_mlk_{date}.grd
##     (also {period}_average composites for the period-mean scenarios)
##
## Usage:
##   Rscript 03_composite_sdm_output.R                     # all models
##   RUN_SDM_MODEL=GFDL Rscript 03_composite_sdm_output.R  # single model
## =========================================================================

library(dplyr)

# Prefer raster for .grd files; fall back to terra; else dry-run.
use_terra  <- requireNamespace("terra",  quietly = TRUE)
use_raster <- requireNamespace("raster", quietly = TRUE)
if (use_raster) {
  library(raster); message("Using raster package (better for .grd files)")
} else if (use_terra) {
  library(terra);  message("Using terra package")
} else {
  message("WARNING: Neither raster nor terra available — DRY RUN (no output written).")
  use_terra <- FALSE; use_raster <- FALSE
}

## --- CONFIGURATION (edit to match the server) ----------------------------
sdm_base_dir    <- "/data/Meghan_SDMs"                    # holds GFDL/ IPSL/ HAD/ daily .grd
output_data_dir <- "/data/Meghan_SDMs/extremes_ccs_ssta"  # extreme-dates CSVs in, composites out

models       <- c("GFDL", "IPSL", "HAD")
species_list <- c("casl", "hbwh", "lbst", "swor")
sig_types    <- c("warm", "cool")   # warm / cool extremes

period_ranges <- list(
  historical = c(as.Date("1991-01-01"), as.Date("2020-12-31")),
  future     = c(as.Date("2071-01-01"), as.Date("2100-12-31"))
)

# Optional single-model override
if (nzchar(Sys.getenv("RUN_SDM_MODEL"))) {
  models <- toupper(trimws(unlist(strsplit(Sys.getenv("RUN_SDM_MODEL"), ","))))
}

## --- Helpers -------------------------------------------------------------
get_model_dir <- function(model) {
  d <- file.path(sdm_base_dir, model)
  if (!dir.exists(d)) {
    alt <- file.path(sdm_base_dir, tools::toTitleCase(tolower(model)))
    if (dir.exists(alt)) return(alt)
    stop("Model directory not found: ", d)
  }
  d
}

find_extremes_csv <- function(model, data_dir) {
  candidates <- list.files(
    data_dir,
    pattern = paste0("extreme_signal_dates_ccs_ssta_", model, "_\\d{4}-\\d{2}-\\d{2}\\.csv$"),
    full.names = TRUE, ignore.case = TRUE
  )
  model_subdir <- file.path(data_dir, model)
  if (dir.exists(model_subdir)) {
    candidates <- c(candidates, list.files(
      model_subdir,
      pattern = "extreme_signal_dates_.*\\d{4}-\\d{2}-\\d{2}\\.csv$",
      full.names = TRUE))
  }
  if (length(candidates) == 0) return(NA_character_)
  info <- file.info(candidates)
  rownames(info)[which.max(info$mtime)]
}

## =========================================================================
## MAIN
## =========================================================================
cat("\n", strrep("=", 70), "\n")
cat("COMPOSITING DAILY SDM OUTPUT BY WARMING SCENARIO\n")
cat(strrep("=", 70), "\n\n")

dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)

for (model in models) {
  cat(strrep("-", 60), "\n"); cat("Model:", model, "\n"); cat(strrep("-", 60), "\n")
  model_dir <- get_model_dir(model)

  csv_path <- find_extremes_csv(model, output_data_dir)
  if (is.na(csv_path)) {
    cat("  x No extreme signal dates CSV found for", model, "— skipping\n\n"); next
  }
  cat("  Extremes file:", basename(csv_path), "\n")

  signals_df <- read.csv(csv_path) %>%
    mutate(date = as.Date(date), yr_mo = format(date, "%Y-%m"))
  if ("extreme_type" %in% names(signals_df) & !"signal" %in% names(signals_df)) {
    signals_df <- signals_df %>% rename(signal = extreme_type)
  }
  cat("  Signals loaded:", nrow(signals_df), "rows\n\n")

  for (sp in species_list) {
    all_sp_files <- list.files(model_dir, pattern = paste0(sp, ".*mean\\.grd$"),
                               full.names = TRUE)
    if (length(all_sp_files) == 0) { cat("  ", sp, ": no daily .grd files — skipping\n"); next }
    cat("  ", sp, ": found", length(all_sp_files), "daily files\n")

    for (period in names(period_ranges)) {
      date_range <- period_ranges[[period]]

      # --- Warm / cool extreme composites ---
      for (sig in sig_types) {
        target_ym <- signals_df %>%
          filter(signal == sig, date >= date_range[1] & date <= date_range[2]) %>%
          pull(yr_mo) %>% unique()
        if (length(target_ym) == 0) { cat("    ", period, sig, ": no matching dates\n"); next }

        ym_pattern     <- paste(target_ym, collapse = "|")
        matching_files <- all_sp_files[grepl(ym_pattern, all_sp_files)]
        if (length(matching_files) == 0) {
          cat("    ", period, sig, ": no matching files for",
              length(target_ym), "year-months\n"); next
        }

        if (use_raster || use_terra) {
          r_stack <- if (use_raster) stack(matching_files) else rast(matching_files)
          r_mean  <- mean(r_stack)
        }
        out_name <- paste0(model, "_", sp, "_", period, "_", sig, "_mean_mlk_", Sys.Date())
        out_path <- file.path(output_data_dir, out_name)
        if (use_raster || use_terra) {
          writeRaster(r_mean, out_path, overwrite = TRUE)
          cat("    v", period, sig, ":", length(matching_files), "files ->",
              basename(out_path), "\n")
        } else {
          cat("    DRY RUN:", period, sig, ":", length(matching_files),
              "files -> would save", basename(out_path), "\n")
        }
      }

      # --- Period-average composite (all months in the 30-year window) ---
      all_ym         <- format(seq(date_range[1], date_range[2], by = "month"), "%Y-%m")
      ym_pattern     <- paste(all_ym, collapse = "|")
      matching_files <- all_sp_files[grepl(ym_pattern, all_sp_files)]
      if (length(matching_files) > 0) {
        if (use_raster || use_terra) {
          r_stack <- if (use_raster) stack(matching_files) else rast(matching_files)
          r_mean  <- mean(r_stack)
        }
        out_name <- paste0(model, "_", sp, "_", period, "_average_mean_mlk_", Sys.Date())
        out_path <- file.path(output_data_dir, out_name)
        if (use_raster || use_terra) {
          writeRaster(r_mean, out_path, overwrite = TRUE)
          cat("    v", period, "average:", length(matching_files), "files ->",
              basename(out_path), "\n")
        } else {
          cat("    DRY RUN:", period, "average:", length(matching_files),
              "files -> would save", basename(out_path), "\n")
        }
      }
    }
    cat("\n")
  }
}

cat(strrep("=", 70), "\n")
cat("COMPOSITING COMPLETE — output:", output_data_dir, "\n")
cat(strrep("=", 70), "\n")
