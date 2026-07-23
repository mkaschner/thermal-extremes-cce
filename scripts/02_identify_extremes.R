## =========================================================================
## 02_identify_extremes.R — Thermal-extreme month identification (S2.3.2)
## =========================================================================
## Inputs  : full_anoms_cropped_{MODEL}.nc  (from 01_preprocess_sst.sh)
##           lme-mask-out-cropped.nc, gridweights.nc
## Outputs : ccs_ssta_{MODEL}_{date}.csv                (regional SSTa series)
##           extreme_signal_dates_ccs_ssta_{MODEL}_{date}.csv (extreme months;
##             consumed by 03_composite_sdm_output.R)
##           suppFig2_ssta_timeSeries_{date}.png        (Supplementary Fig. 2)
##
## Method: apply the California Current LME mask to the gridded SST anomalies,
## take the area-weighted monthly regional mean, then within each 30-year
## period classify the warmest 10% of months (>= 90th percentile anomaly) as
## "warm" extremes and the coolest 10% (<= 10th percentile) as "cool" extremes.
## =========================================================================

source("config.R")
source(file.path("R", "functions.R"))

library(ncdf4)
library(tidyverse)

## =========================================================================
## Stage A — Area-weighted regional SST anomaly time series
## =========================================================================
compute_ccs_ssta <- function(model, data_dir, mask_file, grid_weights_file,
                             ccs_region_id = 3) {
  message("Computing CCS SSTa for: ", model)

  # 1. Read SST anomaly array
  ssta_path <- file.path(data_dir, model, paste0("full_anoms_cropped_", model, ".nc"))
  ssta_nc   <- nc_open(ssta_path)

  sst_array <- ncvar_get(ssta_nc, "thetao")
  fill_val  <- ncatt_get(ssta_nc, "thetao", "_FillValue")$value
  sst_array[sst_array == fill_val] <- NA
  # Note: zero is a valid anomaly — do NOT set 0 to NA

  time_raw   <- ncvar_get(ssta_nc, "time")
  time_units <- ncatt_get(ssta_nc, "time", "units")$value
  calendar   <- ncatt_get(ssta_nc, "time", "calendar")$value
  nc_close(ssta_nc)

  # 2. Calendar-aware date conversion
  dates <- nc_time_to_dates(time_raw, time_units, calendar)
  message(sprintf("  Calendar: %s | %d timesteps | %s to %s",
                  calendar, length(dates),
                  as.character(min(dates)), as.character(max(dates))))

  # 3. Read LME mask — CCS is region ID 3 only
  mask_nc    <- nc_open(mask_file)
  mask_array <- ncvar_get(mask_nc, "mask")
  mask_fill  <- ncatt_get(mask_nc, "mask", "_FillValue")$value
  mask_array[mask_array == mask_fill | mask_array == 0] <- NA
  nc_close(mask_nc)
  mask_ccs <- ifelse(mask_array == ccs_region_id, 1, 0)

  # 4. Read grid cell area weights
  wts_nc       <- nc_open(grid_weights_file)
  grid_weights <- ncvar_get(wts_nc, "cell_weights")
  wts_fill     <- ncatt_get(wts_nc, "cell_weights", "_FillValue")$value
  grid_weights[grid_weights == wts_fill | grid_weights == 0] <- NA
  nc_close(wts_nc)

  # 5. Area-weighted regional mean
  region_weights <- mask_ccs * grid_weights
  norm_factor    <- mean(region_weights, na.rm = TRUE)

  n_times <- dim(sst_array)[3]
  ssta_regional <- vapply(seq_len(n_times), function(t) {
    mean(sst_array[, , t] * region_weights, na.rm = TRUE) / norm_factor
  }, numeric(1))

  # 6. Verify: historical mean must be ~0 (anomalies baselined to 1991-2020)
  hist_idx  <- which(dates >= historical_period[1] & dates <= historical_period[2])
  hist_mean <- mean(ssta_regional[hist_idx], na.rm = TRUE)
  message(sprintf("  Historical (1991-2020): n=%d, mean=%.6f C",
                  length(hist_idx), hist_mean))
  if (length(hist_idx) != 360)
    warning("  Expected 360 timesteps in 1991-2020, got ", length(hist_idx), "!")
  if (abs(hist_mean) > 0.001)
    warning("  Historical mean is not zero — check CDO baseline!")

  fut_idx <- which(dates >= future_period[1] & dates <= future_period[2])
  message(sprintf("  Future (2071-2100): n=%d, mean=%.2f C",
                  length(fut_idx), mean(ssta_regional[fut_idx], na.rm = TRUE)))

  # 7. Save with date stamp
  ccs_ssta <- data.frame(date = dates, ssta_reg = ssta_regional)
  out_path <- file.path(data_dir, model,
                        paste0("ccs_ssta_", model, "_", Sys.Date(), ".csv"))
  write.csv(ccs_ssta, out_path, row.names = FALSE)
  message("  Saved: ", out_path)

  return(ccs_ssta)
}

ccs_ssta_list <- setNames(
  lapply(models, compute_ccs_ssta,
         data_dir          = output_data_dir,
         mask_file         = mask_file,
         grid_weights_file = grid_weights_file,
         ccs_region_id     = ccs_region_id),
  models
)

## =========================================================================
## Stage B — Classify months as warm or cool thermal extremes
## =========================================================================
## Percentile thresholds computed independently within each 30-year period.
identify_extreme_dates <- function(model, ccs_ssta) {
  message("Classifying extreme dates for: ", model)

  ccs_ssta <- ccs_ssta %>% mutate(date = as.Date(date))

  hist_data <- ccs_ssta %>%
    filter(date >= historical_period[1] & date <= historical_period[2])
  fut_data <- ccs_ssta %>%
    filter(date >= future_period[1] & date <= future_period[2])

  categorize_extremes <- function(ssta_df, period_label) {
    p10 <- quantile(ssta_df$ssta_reg, probs = 0.10, na.rm = TRUE)
    p90 <- quantile(ssta_df$ssta_reg, probs = 0.90, na.rm = TRUE)
    message(sprintf("  %s %s: p10 = %.4f C, p90 = %.4f C",
                    model, period_label, p10, p90))

    ssta_df %>%
      mutate(
        extreme_type = case_when(
          ssta_reg <= p10 ~ "cool",   # cool extreme
          ssta_reg >= p90 ~ "warm",    # warm extreme
          TRUE ~ NA_character_
        ),
        period = period_label
      ) %>%
      filter(!is.na(extreme_type))
  }

  all_extremes <- bind_rows(
    categorize_extremes(hist_data, "historical"),
    categorize_extremes(fut_data,  "future")
  ) %>%
    mutate(
      period = factor(period, levels = c("historical", "future")),
      yr_mo  = format(date, "%Y-%m")
    )

  cat("\n"); print(table(all_extremes$period, all_extremes$extreme_type))

  # Name matches the pattern expected by 03_composite_sdm_output.R
  out_path <- file.path(output_data_dir, model,
                        paste0("extreme_signal_dates_ccs_ssta_", model, "_",
                               Sys.Date(), ".csv"))
  write.csv(all_extremes, out_path, row.names = FALSE)
  message("  Saved: ", out_path, "\n")

  return(all_extremes)
}

extremes_list <- setNames(
  mapply(identify_extreme_dates, models, ccs_ssta_list, SIMPLIFY = FALSE),
  models
)

# Quick verification: expect ~36 months per group (10% of 360)
for (m in names(extremes_list)) {
  counts <- extremes_list[[m]] %>%
    group_by(period, extreme_type) %>%
    summarise(n_months = n_distinct(yr_mo), .groups = "drop")
  cat(m, ":\n"); print(as.data.frame(counts)); cat("\n")
}

## =========================================================================
## SUPPLEMENTARY FIGURE 2 — SSTa time series with extreme months highlighted
## =========================================================================
ssta_all <- bind_rows(lapply(names(ccs_ssta_list), function(m) {
  ccs_ssta_list[[m]] %>% mutate(date = as.Date(date), model = m)
}))
extremes_all <- bind_rows(lapply(names(extremes_list), function(m) {
  extremes_list[[m]] %>% mutate(date = as.Date(date), model = m)
}))

# Display labels: Had -> HAD
model_display <- c("Had" = "HAD", "IPSL" = "IPSL", "GFDL" = "GFDL")
ssta_all$model     <- factor(model_display[ssta_all$model],     levels = c("HAD", "IPSL", "GFDL"))
extremes_all$model <- factor(model_display[extremes_all$model], levels = c("HAD", "IPSL", "GFDL"))

# Period mean lines (~0 historical, positive future)
average_lines <- bind_rows(
  ssta_all %>% filter(date >= historical_period[1] & date <= historical_period[2]) %>%
    group_by(model) %>%
    summarise(y = mean(ssta_reg, na.rm = TRUE),
              x_start = historical_period[1], x_end = historical_period[2], .groups = "drop"),
  ssta_all %>% filter(date >= future_period[1] & date <= future_period[2]) %>%
    group_by(model) %>%
    summarise(y = mean(ssta_reg, na.rm = TRUE),
              x_start = future_period[1], x_end = future_period[2], .groups = "drop")
)

highlight_periods <- data.frame(
  start = c(historical_period[1], future_period[1]),
  end   = c(historical_period[2], future_period[2]),
  fill  = c("#a2d6fb", "#f96f66")
)
panel_labels <- c("HAD" = "A. HAD", "IPSL" = "B. IPSL", "GFDL" = "C. GFDL")

supp_fig2 <- ggplot() +
  geom_rect(data = highlight_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_line(data = ssta_all, aes(x = date, y = ssta_reg),
            linewidth = 0.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_segment(data = average_lines,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               color = "#ffaf03", linewidth = 1.25) +
  geom_point(data = extremes_all,
             aes(x = date, y = ssta_reg, color = extreme_type, shape = extreme_type),
             size = 2, alpha = 0.9) +
  scale_color_manual(
    name   = "Thermal extreme",
    values = c("warm" = "#d21919", "cool" = "#1f609e"),
    labels = c("warm" = "Warm extreme", "cool" = "Cool extreme")
  ) +
  scale_shape_manual(
    name   = "Thermal extreme",
    values = c("warm" = 2, "cool" = 6),
    labels = c("warm" = "Warm extreme", "cool" = "Cool extreme")
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years", expand = c(0.01, 0),
               limits = c(historical_period[1], future_period[2])) +
  labs(x = "Time", y = "SST Anomaly (°C)") +
  guides(
    color = guide_legend(override.aes = list(size = 4, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 4, alpha = 1))
  ) +
  facet_wrap(~ model, strip.position = "top", ncol = 1, scales = "free_y",
             labeller = as_labeller(panel_labels)) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", color = "white"),
    axis.line        = element_line(color = "black"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing    = unit(0.8, "lines"),
    strip.background = element_rect(fill = NA, color = NA),
    strip.placement  = "outside",
    strip.text       = element_text(hjust = 0, face = "bold", size = 16),
    legend.position  = "bottom",
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(size = 12)
  )

ggsave(
  filename = file.path(output_plots_dir,
                       paste0("suppFig2_ssta_timeSeries_", Sys.Date(), ".png")),
  plot = supp_fig2, width = 9, height = 9.5, units = "in", dpi = 600, bg = "white"
)
message("Supplementary Figure 2 saved.")
