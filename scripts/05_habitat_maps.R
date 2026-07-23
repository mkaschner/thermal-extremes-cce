## =========================================================================
## 05_habitat_maps.R — Habitat-suitability difference maps (Figures 2-4)
## =========================================================================
## Inputs  : composite habitat-suitability .grd files (03_composite_sdm_output.R)
## Outputs : individual map panels (PNG, transparent background) + a table of
##           per-map spatial extremes. Panels are assembled into the manuscript
##           multi-panel figures externally.
##
## Figure mapping (manuscript numbering):
##   Fig 2  Long-term warming        F[ave] - H[ave]   -> generate_longterm_change_maps
##   Fig 3  Historical thermal extr. H[warm] - H[cool]   -> generate_interannual_maps (historical)
##   Fig 4  Future thermal extremes  F[warm] - F[cool]   -> generate_interannual_maps (future)
##   Mean habitat-suitability maps (H[ave], F[ave], ...) -> generate_mean_maps
## =========================================================================

source("config.R")
source(file.path("R", "functions.R"))

library(ncdf4)
library(raster)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)
library(ragg)

## --- Base map + shared theme ---------------------------------------------
wrld2 <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -134, xmax = -114, ymin = 30, ymax = 48) %>%
  st_make_valid()

coast_linewidth <- 0.6

theme_map_common <- list(
  coord_sf(xlim = c(-130, -116), ylim = c(31, 46.5), expand = FALSE),
  scale_x_continuous(breaks = seq(-130, -116, by = 2),
                     labels = function(x) sprintf("%d", x)),
  scale_y_continuous(breaks = seq(32, 46, by = 2),
                     labels = function(y) sprintf("%d", y)),
  xlab("Longitude (°E)"), ylab("Latitude (°N)"),
  theme_minimal(),
  theme(
    legend.position  = c(0.73, 0.66),
    legend.spacing.x = unit(0.5, "cm"),
    legend.key.size  = unit(1, "cm"),
    legend.text      = element_text(size = 20),
    legend.title     = element_text(size = 18),
    axis.text        = element_text(size = 18, color = "black"),
    axis.title       = element_text(size = 20),
    panel.grid.major = element_line(linewidth = 0.5, color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = NA, color = NA),
    axis.line        = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.5),
    panel.background = element_rect(fill = "grey95")
  )
)

## --- Colour palettes ------------------------------------------------------
hs_colors   <- colorRampPalette(c("#f4f8fa", "#4689d4", "#0a192a"))(4)  # 0-1
diff_colors <- c("#164542", "#359b9b", "#f6f3f4", "#d21919", "#720e0e") # diverging
abs_colors  <- c("#1f609e", "#a0b3c7", "#f6f3f4")                       # sequential

guide_vertical <- guides(fill = guide_colorbar(
  title.position = "top", title.hjust = 0.5, title.vjust = 3,
  direction = "vertical", barwidth = 3, barheight = 14, ticks = TRUE))

## --- Plot builders --------------------------------------------------------
plot_hs_map <- function(df) {
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value), interpolate = TRUE) +
    geom_sf(data = wrld2, fill = "grey", color = "black", linewidth = coast_linewidth) +
    scale_fill_gradientn(
      "Habitat suitability",
      colours = hs_colors, limits = c(0, 1), na.value = "black",
      values = scales::rescale(seq(0, 1, 0.1)), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_map_common + guide_vertical
}

plot_diff_map <- function(df, title_expr, range_limits, legend_title_size = 18) {
  breaks5 <- seq(range_limits[1], range_limits[2], length.out = 5)
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value), interpolate = TRUE) +
    geom_sf(data = wrld2, fill = "grey", color = "black", linewidth = coast_linewidth) +
    scale_fill_gradientn(
      title_expr, colours = colorRampPalette(rev(diff_colors))(300),
      limits = range_limits, values = scales::rescale(c(range_limits[1], 0, range_limits[2])),
      breaks = breaks5,
      labels = sapply(signif(breaks5, 2), function(x) if (x == 0) "0" else sprintf("%+g", x)),
      na.value = "black") +
    theme_map_common + guide_vertical +
    theme(legend.title = element_text(size = legend_title_size))
}

plot_abs_map <- function(df, title_expr, range_limits) {
  breaks5 <- seq(range_limits[1], range_limits[2], length.out = 5)
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = value), interpolate = TRUE) +
    geom_sf(data = wrld2, fill = "grey", color = "black", linewidth = coast_linewidth) +
    scale_fill_gradientn(
      title_expr, colours = colorRampPalette(rev(abs_colors))(300),
      limits = range_limits, values = scales::rescale(range_limits),
      breaks = breaks5, labels = signif(breaks5, 2), na.value = "black") +
    theme_map_common + guide_vertical
}

## =========================================================================
## A) Mean habitat-suitability maps (H[ave], F[ave], warm/cool composites)
## =========================================================================
generate_mean_maps <- function(models, species_list) {
  message("\n--- Mean habitat suitability maps ---")
  plots <- list(); stats <- list()
  signals <- c("average", "warm", "cool")
  periods <- c("historical", "future")

  for (sp in species_list) {
    for (prd in periods) for (sig in signals) {
      model_rasters <- list()
      for (m in models) {
        r <- load_grd(m, sp, prd, sig)
        if (!is.null(r)) model_rasters[[m]] <- r
      }
      if (length(model_rasters) == 0) next

      for (m in names(model_rasters)) {
        df  <- raster_to_df(model_rasters[[m]])
        key <- paste0(sp, "_", m, "_", prd, "_", sig)
        plots[[key]] <- plot_hs_map(df)
        stats[[key]] <- get_spatial_extremes(df, sp, m, paste(sp, m, prd, sig, sep = "_"))
      }
      if (length(model_rasters) == length(models)) {
        ens <- ensemble_mean(model_rasters); df <- raster_to_df(ens)
        key <- paste0(sp, "_ensemble_", prd, "_", sig)
        plots[[key]] <- plot_hs_map(df)
        stats[[key]] <- get_spatial_extremes(df, sp, "ensemble", paste(sp, "ensemble", prd, sig, sep = "_"))
      }
    }
    message("  ok ", toupper(sp))
  }
  list(plots = plots, stats = bind_rows(stats))
}

## =========================================================================
## B) Figure 2 — Long-term change maps: F[ave] - H[ave]
## =========================================================================
generate_longterm_change_maps <- function(models, species_list) {
  message("\n--- Long-term change maps (F[ave] - H[ave]) [Fig 2] ---")
  plots <- list(); stats <- list(); rasters <- list()

  for (sp in species_list) {
    model_diffs <- list()
    for (m in models) {
      r_fut  <- load_grd(m, sp, "future", "average")
      r_hist <- load_grd(m, sp, "historical", "average")
      if (is.null(r_fut) | is.null(r_hist)) next
      model_diffs[[m]] <- r_fut - r_hist
    }
    if (length(model_diffs) == 0) next

    ens <- if (length(model_diffs) == length(models)) ensemble_mean(model_diffs) else NULL

    model_vals <- unlist(lapply(model_diffs, values))
    model_max  <- quantile(abs(model_vals), 1, na.rm = TRUE)
    model_rng  <- c(-model_max, model_max)
    ens_rng <- if (!is.null(ens)) { em <- quantile(abs(values(ens)), 1, na.rm = TRUE); c(-em, em) } else NULL

    for (m in names(model_diffs)) {
      df <- raster_to_df(model_diffs[[m]]); key <- paste0(sp, "_", m, "_longterm")
      plots[[key]]   <- plot_diff_map(df, expression(F[ave] - H[ave]), model_rng)
      stats[[key]]   <- get_spatial_extremes(df, sp, m, key)
      rasters[[key]] <- model_diffs[[m]]
    }
    if (!is.null(ens)) {
      df <- raster_to_df(ens); key <- paste0(sp, "_ensemble_longterm")
      plots[[key]]   <- plot_diff_map(df, expression(F[ave] - H[ave]), ens_rng)
      stats[[key]]   <- get_spatial_extremes(df, sp, "ensemble", key)
      rasters[[key]] <- ens
    }
    message("  ok ", toupper(sp))
  }
  list(plots = plots, stats = bind_rows(stats), rasters = rasters)
}

## =========================================================================
## C) Figures 3 & 4 — Interannual variability maps
##    Fig 3: H[warm] - H[cool] (historical)   Fig 4: F[warm] - F[cool] (future)
## =========================================================================
generate_interannual_maps <- function(models, species_list) {
  message("\n--- Interannual variability maps (amp - min) [Figs 3 & 4] ---")
  plots <- list(); stats <- list(); rasters <- list()

  for (sp in species_list) {
    model_diffs <- list()
    for (m in models) {
      diffs_m <- list()
      for (prd in c("historical", "future")) {
        r_amp <- load_grd(m, sp, prd, "warm")
        r_min <- load_grd(m, sp, prd, "cool")
        if (is.null(r_amp) | is.null(r_min)) next
        diffs_m[[prd]] <- r_amp - r_min
      }
      if (length(diffs_m) > 0) model_diffs[[m]] <- diffs_m
    }
    if (length(model_diffs) == 0) next

    model_vals <- c()
    for (m in names(model_diffs)) for (prd in names(model_diffs[[m]]))
      model_vals <- c(model_vals, values(model_diffs[[m]][[prd]]))
    model_max <- quantile(abs(model_vals), 1, na.rm = TRUE)
    model_rng <- c(-model_max, model_max)

    title_exprs <- list(historical = expression(H[warm] - H[cool]),
                        future     = expression(F[warm] - F[cool]))

    for (prd in c("historical", "future")) {
      prd_rasters <- lapply(model_diffs, function(x) x[[prd]])
      prd_rasters <- prd_rasters[!sapply(prd_rasters, is.null)]

      ens <- NULL; ens_rng <- NULL
      if (length(prd_rasters) == length(models)) {
        ens <- ensemble_mean(prd_rasters)
        em  <- quantile(abs(values(ens)), 1, na.rm = TRUE); ens_rng <- c(-em, em)
      }
      for (m in names(model_diffs)) {
        if (!prd %in% names(model_diffs[[m]])) next
        df <- raster_to_df(model_diffs[[m]][[prd]])
        key <- paste0(sp, "_", m, "_", prd, "_interannual")
        plots[[key]] <- plot_diff_map(df, title_exprs[[prd]], model_rng)
        stats[[key]] <- get_spatial_extremes(df, sp, m, key)
      }
      if (!is.null(ens)) {
        df <- raster_to_df(ens); key <- paste0(sp, "_ensemble_", prd, "_interannual")
        plots[[key]]   <- plot_diff_map(df, title_exprs[[prd]], ens_rng)
        stats[[key]]   <- get_spatial_extremes(df, sp, "ensemble", key)
        rasters[[key]] <- ens
      }
    }
    message("  ok ", toupper(sp))
  }
  list(plots = plots, stats = bind_rows(stats), rasters = rasters)
}

## =========================================================================
## RUN ALL MAP GENERATION
## =========================================================================
mean_maps        <- generate_mean_maps(models, species_list)
longterm_maps    <- generate_longterm_change_maps(models, species_list)   # Fig 2
interannual_maps <- generate_interannual_maps(models, species_list)       # Figs 3 & 4

all_map_stats <- bind_rows(mean_maps$stats, longterm_maps$stats,
                           interannual_maps$stats)
write.csv(all_map_stats,
          file.path(output_data_dir, paste0("all_map_spatial_stats_", Sys.Date(), ".csv")),
          row.names = FALSE)
message("\nAll map stats saved.")

## --- Save individual map panels (transparent background) ------------------
save_map_plots <- function(plot_list, output_dir, width_in = 5, height_in = 7, dpi = 750) {
  fig_dir <- file.path(output_dir, "figures", "maps")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(plot_list)) {
    p <- plot_list[[nm]]
    if (!inherits(p, "ggplot")) next
    out_path <- file.path(fig_dir, paste0(nm, "_", Sys.Date(), ".png"))
    ragg::agg_png(filename = out_path, width = width_in, height = height_in,
                  units = "in", res = dpi, background = NA)
    print(p); dev.off()
  }
  message("Maps saved to: ", fig_dir)
}

save_map_plots(mean_maps$plots,        output_plots_dir)
save_map_plots(longterm_maps$plots,    output_plots_dir)
save_map_plots(interannual_maps$plots, output_plots_dir)
