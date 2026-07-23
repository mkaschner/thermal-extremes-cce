## =========================================================================
## 04_centroid_area_metrics.R — Core-habitat centroids & area (S2.3.3)
## =========================================================================
## Inputs  : composite habitat-suitability .grd files (03_composite_sdm_output.R)
## Outputs : centroid_results.csv  (per species x ESM x period x signal)
##             consumed by 06_figure_summary.R
##
## For each species x ESM the core-habitat threshold is the 75th percentile of
## the HISTORICAL AVERAGE composite (get_presence_thresholds), held constant
## across all scenario composites. Within each composite, cells >= threshold
## are "core habitat". Reported per composite:
##   - centroid  : habitat-suitability-weighted (HS^1) mean lon/lat of core cells
##   - core_area_km2 : unweighted sum of core-cell areas on the WGS84 ellipsoid
## These are the definitions used in the manuscript (Fig. 5 metrics). Legacy
## HS^2-weighted centroid and suitability-weighted area variants were dropped.
## =========================================================================

source("config.R")
source(file.path("R", "functions.R"))

library(raster)
library(dplyr)
library(geosphere)

generate_centroid_results <- function(model, species_list) {
  results_all  <- data.frame()
  threshold_df <- get_presence_thresholds(model, species_list)

  signals <- c("average", "warm", "cool")
  periods <- c("historical", "future")

  for (sp in species_list) {
    message(sp, " ...")
    thresh <- threshold_df$core_thresh_75[threshold_df$species == sp]
    if (is.na(thresh)) next

    for (prd in periods) {
      for (sig in signals) {
        r <- load_grd(model, sp, prd, sig)
        if (is.null(r)) { message("  Missing: ", model, " ", sp, " ", prd, " ", sig); next }

        if (is.na(crs(r))) {
          crs(r) <- CRS("+proj=longlat +datum=WGS84")
        } else {
          r <- projectRaster(r, crs = CRS("+proj=longlat +datum=WGS84"))
        }

        # Mask to core habitat (>= threshold)
        r_masked <- r
        r_masked[r_masked < thresh] <- NA

        df <- rasterToPoints(r_masked) %>% as.data.frame() %>%
          setNames(c("x", "y", "value"))
        if (nrow(df) == 0) next

        df$cell_area_km2 <- mapply(calc_cell_area, df$y)

        # Habitat-suitability-weighted (HS^1) centroid
        w              <- df$value
        centroid_x_hs1 <- sum(df$x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
        centroid_y_hs1 <- sum(df$y * w, na.rm = TRUE) / sum(w, na.rm = TRUE)

        # Core area: unweighted sum of core-cell areas
        core_area_km2 <- sum(df$cell_area_km2, na.rm = TRUE)

        results_all <- rbind(results_all, data.frame(
          model          = model,
          species        = sp,
          period         = prd,
          signal         = sig,
          core_thresh    = thresh,
          centroid_x_hs1 = centroid_x_hs1,
          centroid_y_hs1 = centroid_y_hs1,
          core_area_km2  = core_area_km2,
          stringsAsFactors = FALSE
        ))
      }
    }
    message("  ok ", toupper(sp))
  }
  results_all
}

# Run for all ESMs
centroid_results <- bind_rows(lapply(models, generate_centroid_results,
                                     species_list = species_list))

# Canonical output (loaded by 06_figure_summary.R) + a date-stamped copy
write.csv(centroid_results,
          file.path(output_data_dir, "centroid_results.csv"), row.names = FALSE)
write.csv(centroid_results,
          file.path(output_data_dir,
                    paste0("centroid_results_", Sys.Date(), ".csv")), row.names = FALSE)
message("Saved centroid_results.csv")
