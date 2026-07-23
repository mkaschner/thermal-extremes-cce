## =========================================================================
## 06_figure_summary.R — Figure 5: displacement distance | area | direction
## =========================================================================
## Input : centroid_results.csv (04_centroid_area_metrics.R)
## Output: figure5_summary_{date}.png  (arrows | circles | unit vectors)
##         Final panel placement was refined externally (PowerPoint).
##
## Panel A = displacement distance (arrows, mean ± SD)
## Panel B = change in core habitat area (circles: size = |absolute change|,
##           fill = % change)
## Panel C = displacement direction (unit vectors with circular-SD cones)
##
## Displacement distance = geosphere::distHaversine between core-habitat
## centroids; bearing = geosphere::bearing; ensemble mean bearing = circular
## mean via atan2; ensemble spread = circular SD (Fisher 1993). Metrics are
## computed per ESM, then averaged across the three-ESM ensemble.
## =========================================================================

source("config.R")

library(geosphere)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)

centroid_results <- read.csv(file.path(output_data_dir, "centroid_results.csv"),
                             stringsAsFactors = FALSE)

## --- Labels & colours -----------------------------------------------------
comp_labels_plain <- c(
  "F[ave]-H[ave]"   = "Anthropogenic\nwarming",
  "H[warm]-H[cool]" = "Historical\nthermal extremes",
  "F[warm]-F[cool]" = "Future\nthermal extremes"
)
comp_colors <- c(
  "F[ave]-H[ave]"   = "#4e9c81",
  "H[warm]-H[cool]" = "#7d87b2",
  "F[warm]-F[cool]" = "#b5622a"
)
species_order_list <- c("casl", "hbwh", "lbst", "swor")

## =========================================================================
## DATA PREPARATION
## =========================================================================
# Ensemble-mean centroids (HS^1 centroid + unweighted core area)
ensemble_centroids <- centroid_results %>%
  group_by(species, period, signal) %>%
  summarise(
    centroid_x_hs1 = mean(centroid_x_hs1, na.rm = TRUE),
    centroid_y_hs1 = mean(centroid_y_hs1, na.rm = TRUE),
    core_area_km2  = mean(core_area_km2,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(model = "ensemble")

all_centroids <- bind_rows(centroid_results, ensemble_centroids)

# Per-model displacement between two (period, signal) scenarios
per_model_displacements <- function(from_period, from_signal,
                                    to_period,   to_signal,   comp_label) {
  from <- centroid_results %>%
    filter(period == from_period, signal == from_signal) %>%
    select(model, species, x_from = centroid_x_hs1, y_from = centroid_y_hs1,
           area_from = core_area_km2)
  to <- centroid_results %>%
    filter(period == to_period, signal == to_signal) %>%
    select(model, species, x_to = centroid_x_hs1, y_to = centroid_y_hs1,
           area_to = core_area_km2)

  left_join(from, to, by = c("model", "species")) %>%
    mutate(
      comparison      = comp_label,
      dist_km         = mapply(function(x1, y1, x2, y2)
        distHaversine(c(x1, y1), c(x2, y2)) / 1000, x_from, y_from, x_to, y_to),
      bearing_deg     = (mapply(function(x1, y1, x2, y2)
        bearing(c(x1, y1), c(x2, y2)), x_from, y_from, x_to, y_to) + 360) %% 360,
      area_abs_change = area_to - area_from,
      area_pct_change = (area_to - area_from) / area_from * 100
    )
}

per_model_data <- bind_rows(
  per_model_displacements("historical", "average", "future",     "average", "F[ave]-H[ave]"),
  per_model_displacements("historical", "cool",    "historical", "warm",    "H[warm]-H[cool]"),
  per_model_displacements("future",     "cool",    "future",     "warm",    "F[warm]-F[cool]")
)

# Ensemble means + circular statistics for bearing (Fisher 1993)
summary_data_uncertainty <- per_model_data %>%
  group_by(species, comparison) %>%
  summarise(
    dist_km_mean  = mean(dist_km, na.rm = TRUE),
    dist_km_sd    = sd(dist_km,   na.rm = TRUE),
    bearing_mean  = (atan2(mean(sin(bearing_deg * pi / 180), na.rm = TRUE),
                           mean(cos(bearing_deg * pi / 180), na.rm = TRUE)) *
                       180 / pi + 360) %% 360,
    bearing_sd    = sqrt(-2 * log(sqrt(
      mean(sin(bearing_deg * pi / 180), na.rm = TRUE)^2 +
        mean(cos(bearing_deg * pi / 180), na.rm = TRUE)^2))) * 180 / pi,
    area_pct_mean = mean(area_pct_change, na.rm = TRUE),
    area_pct_sd   = sd(area_pct_change,   na.rm = TRUE),
    area_abs_mean = mean(area_abs_change, na.rm = TRUE),
    area_abs_sd   = sd(area_abs_change,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    comparison   = factor(comparison, levels = comp_levels),
    species_name = factor(species, levels = c("casl", "hbwh", "lbst", "swor"),
                          labels = c("California\nsea lion", "Humpback\nwhale",
                                     "Leatherback\nsea turtle", "Swordfish")),
    bearing_rad  = bearing_mean * pi / 180
  )

## --- Shared elements ------------------------------------------------------
size_limit <- max(abs(summary_data_uncertainty$area_abs_mean) +
                    summary_data_uncertainty$area_abs_sd, na.rm = TRUE)

# Species image column (falls back to text if illustrations are absent)
species_drawings <- lapply(species_order_list, function(sp) {
  img_path <- file.path(output_plots_dir, "images", paste0(toupper(sp), "_img.jpg"))
  if (file.exists(img_path)) {
    ggdraw() + draw_image(img_path, scale = 0.95) + theme(plot.margin = margin(0, 0, 0, 0))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = toupper(sp),
                        fontface = "bold", size = 4) + theme_void()
  }
})
species_column <- plot_grid(plotlist = species_drawings, ncol = 1, align = "v", axis = "tb")

# Unit-vector cone data (Panel C)
unit_data <- summary_data_uncertainty %>%
  mutate(x_unit = sin(bearing_rad), y_unit = cos(bearing_rad))

unit_cone_data <- do.call(rbind, lapply(seq_len(nrow(summary_data_uncertainty)), function(i) {
  row     <- summary_data_uncertainty[i, ]
  bear_lo <- (row$bearing_mean - row$bearing_sd) * pi / 180
  bear_hi <- (row$bearing_mean + row$bearing_sd) * pi / 180
  arc_seq <- seq(bear_lo, bear_hi, length.out = 40)
  data.frame(x = c(0, sin(arc_seq), 0), y = c(0, cos(arc_seq), 0),
             species_name = row$species_name, comparison = as.character(row$comparison))
})) %>% mutate(comparison = factor(comparison, levels = comp_levels))

## =========================================================================
## FIGURE — Arrows | circles | unit vectors
## =========================================================================
arrow_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position  = "none", strip.text = element_blank(),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5)
  )

# --- Panel A: displacement distance (arrows) -----------------------------
arrow_x_pos <- c("F[ave]-H[ave]" = 1, "H[warm]-H[cool]" = 2, "F[warm]-F[cool]" = 3)
p_A_dat <- summary_data_uncertainty %>%
  mutate(x_pos_arrow = arrow_x_pos[as.character(comparison)])

p_A <- ggplot(p_A_dat, aes(x = x_pos_arrow, color = comparison)) +
  geom_linerange(aes(ymin = dist_km_mean, ymax = dist_km_mean + dist_km_sd),
                 linewidth = 7, alpha = 0.4) +
  geom_linerange(aes(ymin = pmax(0, dist_km_mean - dist_km_sd), ymax = dist_km_mean),
                 linewidth = 7, alpha = 0.4) +
  geom_segment(aes(y = 0, yend = dist_km_mean, xend = x_pos_arrow), linewidth = 1.8,
               arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  scale_color_manual(values = comp_colors, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("Anthropogenic\nwarming", "Historical\nthermal extremes",
                                "Future\nthermal extremes"), limits = c(0.5, 3.5)) +
  facet_grid(species_name ~ ., switch = "y", scales = "free_y") +
  labs(x = NULL, y = "Displacement (km)", title = "A. Displacement\ndistance") +
  arrow_theme +
  theme(axis.text.x = element_text(size = 9, hjust = 0.5),
        axis.text.y = element_text(size = 9), strip.text.y.left = element_blank())

# --- Panel B: change in core habitat area (circles) ----------------------
p_B <- ggplot() +
  geom_point(data = summary_data_uncertainty,
             aes(y = species_name, x = comparison, size = abs(area_abs_mean),
                 fill = area_pct_mean), shape = 21, color = "black", alpha = 1, stroke = 1) +
  scale_fill_gradient2(high = "#024d72", mid = "#f6f3f4", low = "#750c00", midpoint = 0,
                       limits = c(-100, 100), oob = scales::squish, name = "% area\nchange",
                       guide = guide_colorbar(barwidth = 6, barheight = 0.8,
                                              title.vjust = 0.75, direction = "horizontal")) +
  scale_size_area(max_size = 50, limits = c(0, size_limit), guide = "none") +
  scale_x_discrete(labels = comp_labels_plain) +
  facet_grid(species_name ~ ., scales = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "B. Change in core\nhabitat area (km²)") +
  arrow_theme +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 9),
        panel.grid.major = element_blank(), legend.position = "bottom",
        legend.title = element_text(size = 9), legend.text = element_text(size = 8))

# --- Panel C: displacement direction (unit vectors + cones) --------------
make_unit_panel <- function(sp) {
  sp_unit  <- unit_data      %>% filter(species_name == sp)
  sp_cones <- unit_cone_data %>% filter(species_name == sp)
  ggplot() +
    geom_polygon(data = sp_cones, aes(x = x, y = y, fill = comparison, group = comparison),
                 alpha = 0.2, color = NA) +
    geom_segment(data = sp_unit, aes(x = 0, y = 0, xend = x_unit, yend = y_unit, color = comparison),
                 linewidth = 1.3, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    geom_point(aes(x = 0, y = 0), size = 2, color = "black", fill = "white", shape = 21) +
    annotate("text", x = 0, y = 1.25, label = "N", size = 3, color = "grey50", fontface = "bold") +
    annotate("text", x = 0, y = -1.25, label = "S", size = 3, color = "grey50", fontface = "bold") +
    annotate("text", x = 1.25, y = 0, label = "E", size = 3, color = "grey50", fontface = "bold") +
    annotate("text", x = -1.25, y = 0, label = "W", size = 3, color = "grey50", fontface = "bold") +
    scale_fill_manual(values = comp_colors, guide = "none") +
    scale_color_manual(values = comp_colors, guide = "none") +
    coord_equal(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    arrow_theme +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          plot.margin = margin(1, 1, 1, 1))
}

sp_levels   <- levels(summary_data_uncertainty$species_name)
unit_panels <- lapply(sp_levels, make_unit_panel)
p_C <- plot_grid(
  ggdraw() + draw_label("C. Displacement direction", fontface = "bold", size = 13),
  plot_grid(plotlist = unit_panels, ncol = 1, align = "v", axis = "lr"),
  ncol = 1, rel_heights = c(0.05, 1))

# --- Assemble ------------------------------------------------------------
legend_B <- get_legend(p_B + theme(legend.position = "right",
                                   legend.title = element_text(face = "bold", size = 10),
                                   legend.text = element_text(size = 9)))
p_A <- p_A + theme(strip.text = element_blank(), strip.text.y.left = element_blank())
p_B <- p_B + theme(strip.text = element_blank(), strip.text.y.left = element_blank(),
                   legend.position = "none",
                   panel.background = element_rect(fill = "transparent", color = NA),
                   plot.background = element_rect(fill = "transparent", color = NA),
                   panel.border = element_blank())

panels_hyp <- plot_grid(p_A, p_B, p_C, nrow = 1, align = "h", axis = "tb",
                        rel_widths = c(0.5, 0.8, 0.7))
full_hyp <- plot_grid(species_column, panels_hyp, nrow = 1, align = "h", axis = "tb",
                      rel_widths = c(0.15, 3))

# Scale bar for the area circles
scale_breaks <- c(100000, 300000, 500000, 700000)
scale_bar <- ggplot(data.frame(x = c(1, 1.5, 2.3, 3.5), size = scale_breaks),
                    aes(x = x, y = 1, size = size)) +
  geom_point(shape = 21, fill = NA, color = "black", alpha = 1) +
  scale_size_area(max_size = 50, limits = c(0, size_limit)) +
  annotate("text", x = c(1, 1.5, 2.3, 3.5), y = 1,
           label = c("100k", "300k", "500k", "700k"), size = 3, hjust = 0.5) +
  coord_cartesian(xlim = c(0.3, 5), ylim = c(0.1, 2)) +
  theme_void() + theme(legend.position = "none")

legends_hyp <- plot_grid(scale_bar, legend_B, nrow = 1, rel_widths = c(1.5, 0.5))
figure5 <- plot_grid(full_hyp, legends_hyp, ncol = 1, rel_heights = c(0.8, 0.2))

ggsave(file.path(output_plots_dir, paste0("figure5_summary_", Sys.Date(), ".png")),
       figure5, width = 14, height = 12, dpi = 600, bg = "transparent")
message("Figure 5 saved.")
