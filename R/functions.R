## =========================================================================
## R/functions.R — Shared helper functions
## =========================================================================
## Portable helpers used across multiple pipeline scripts. Source AFTER
## config.R (several functions reference paths/constants defined there:
## composites_dir, models). Requires the `raster` and `dplyr` packages.
## =========================================================================

## --- Calendar-aware NetCDF date conversion -------------------------------
## Converts NetCDF time values to R Date objects, handling the non-standard
## calendars used by CMIP5 ESMs (GFDL/IPSL: 365_day noleap; HadGEM2: 360_day).
nc_time_to_dates <- function(time_raw, time_units, calendar) {
  # Parse base date from units string (e.g., "days since 1861-01-01 00:00:00")
  base_str   <- sub(".*since\\s+", "", time_units)
  base_str   <- sub("\\s+\\d{2}:\\d{2}.*", "", base_str)
  base_parts <- as.numeric(strsplit(base_str, "-")[[1]])
  base_year  <- base_parts[1]
  base_month <- base_parts[2]
  base_day   <- base_parts[3]

  calendar <- tolower(trimws(calendar))

  if (grepl("360", calendar)) {
    ## 360_day: 12 months x 30 days = 360 days/year
    base_total_days <- base_year * 360 + (base_month - 1) * 30 + (base_day - 1)
    total_days  <- base_total_days + time_raw
    years       <- floor(total_days / 360)
    day_of_year <- total_days - years * 360
    months      <- pmin(floor(day_of_year / 30) + 1, 12)
    dates       <- as.Date(paste(years, months, 15, sep = "-"))

  } else if (grepl("365|noleap", calendar)) {
    ## 365_day (noleap): standard month lengths, no Feb 29
    days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    cum_days      <- c(0, cumsum(days_in_month))
    base_doy        <- cum_days[base_month] + (base_day - 1)
    base_total_days <- base_year * 365 + base_doy
    total_days  <- base_total_days + time_raw
    years       <- floor(total_days / 365)
    day_of_year <- total_days - years * 365
    months      <- pmin(sapply(day_of_year, function(d) max(1, sum(d >= cum_days))), 12)
    dates       <- as.Date(paste(years, months, 15, sep = "-"))

  } else {
    ## Standard/Gregorian calendar
    dates <- as.Date(base_str) + time_raw
  }

  return(dates)
}

## --- Grid cell area on the WGS84 ellipsoid (manuscript S2.3.3) -------------
## Returns cell area in km^2 for a cell centred at `lat_deg`.
calc_cell_area <- function(lat_deg, res_lat = 0.1, res_lon = 0.1) {
  R <- 6371  # km (approx. Earth radius)
  to_rad <- pi / 180
  R^2 * (to_rad * res_lat) * (to_rad * res_lon) * cos(to_rad * lat_deg)
}

## --- Composite raster discovery ------------------------------------------
## Finds a composite habitat-suitability .grd (built by
## 03_composite_sdm_output.R). Handles two naming conventions:
##   New: {MODEL}_{species}_{period}_{signal}_mean_mlk_{date}.grd
##   Old: {period}_{signal}_{species}.grd  (in model/sdm_dir/)
find_composite_grd <- function(model, species, period, signal,
                               data_dir = composites_dir,
                               sdm_dir = "sdm output_may25") {
  # Translate the canonical signal code (average/warm/cool) to the token used
  # in the .grd FILENAMES (see signal_file_tokens in config.R). Falls back to
  # the signal itself if no map is defined.
  file_signal <- if (exists("signal_file_tokens") &&
                     signal %in% names(signal_file_tokens)) {
    signal_file_tokens[[signal]]
  } else {
    signal
  }

  # Signal name mapping: "average" -> "{species}_mlk" for old convention
  signal_old <- if (file_signal == "average") paste0(species, "_mlk") else file_signal

  # Try new naming first (flat output directory)
  new_pattern <- paste0(model, "_", species, "_", period, "_", file_signal,
                        "_mean_mlk.*\\.grd$")
  new_files <- list.files(data_dir, pattern = new_pattern, full.names = TRUE,
                          ignore.case = TRUE)
  if (length(new_files) > 0) return(new_files[length(new_files)])  # most recent

  # Try old naming (in model/sdm_dir/)
  old_dir <- file.path(data_dir, model, sdm_dir)
  if (dir.exists(old_dir)) {
    old_files <- list.files(old_dir, full.names = TRUE)
    matched <- old_files[
      grepl(period, old_files, ignore.case = TRUE) &
        grepl(signal_old, old_files) &
        grepl(species, old_files) &
        grepl("\\.grd$", old_files)
    ]
    if (length(matched) > 0) return(matched[1])
  }

  return(NA_character_)
}

## --- Load one composite raster -------------------------------------------
load_grd <- function(model, species, period, signal) {
  fn <- find_composite_grd(model, species, period, signal, composites_dir)
  if (is.na(fn)) return(NULL)
  raster(fn)
}

## --- Ensemble mean of a list of rasters ----------------------------------
ensemble_mean <- function(raster_list) {
  calc(stack(raster_list), fun = mean, na.rm = TRUE)
}

## --- 75th-percentile core-habitat threshold (manuscript S2.3.3) ----------
## Threshold defined per species x ESM from the HISTORICAL AVERAGE composite,
## then held constant across all scenario composites for that species x ESM.
get_presence_thresholds <- function(model_esm, species_vec) {
  thresholds <- data.frame(species = species_vec, core_thresh_75 = NA)
  for (i in seq_along(species_vec)) {
    sp <- species_vec[i]
    r  <- load_grd(model_esm, sp, "historical", "average")
    if (is.null(r)) {
      warning("No historical mean file found for: ", sp, " in ESM: ", model_esm)
      next
    }
    thresholds$core_thresh_75[i] <- quantile(r[], probs = 0.75, na.rm = TRUE)
  }
  return(thresholds)
}

## --- Raster -> tidy data frame -------------------------------------------
raster_to_df <- function(r) {
  df <- rasterToPoints(r) %>% as.data.frame()
  colnames(df) <- c("lon", "lat", "value")
  df
}

## --- Spatial extremes (max/min values + coordinates) ---------------------
## Records the highest and lowest cell value (and location) of a map, used to
## report spatial statistics quoted in the Results.
get_spatial_extremes <- function(df, species, model, label) {
  dplyr::bind_rows(
    df %>% dplyr::slice_max(value, n = 1, with_ties = TRUE) %>%
      dplyr::mutate(max_or_min = "max"),
    df %>% dplyr::slice_min(value, n = 1, with_ties = TRUE) %>%
      dplyr::mutate(max_or_min = "min")
  ) %>%
    dplyr::mutate(species = toupper(species), model = model,
                  id = label, value = round(value, 3))
}
