## =========================================================================
## 07_partial_dependence.R — Supp Fig 1 + Supp Table 1 (variable importance)
## =========================================================================
## Inputs : BRT ensemble .rds objects (config.R: sdm_model_paths)
## Outputs: suppFig1_partial_dependence_{date}.png     (Supplementary Fig. 1)
##          suppTable1_variable_contributions_{date}.csv (Supplementary Table 1)
##
## For each species, loads the 10-model BRT ensemble, extracts relative
## variable contributions (summary.gbm) and partial dependence (plot.gbm),
## selects the three highest-contribution predictors, and plots their partial
## dependence curves (grey = individual models; colored = ensemble smooth).
## =========================================================================

source("config.R")

library(gbm)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# Species display config: name + accent colour (paths come from config.R)
species_colors <- c(casl = "#6f5400", hbwh = "#5b264a",
                    lbst = "#314731", swor = "#274f58")
species_config <- lapply(species_list, function(sp) {
  list(name = species_labels[[sp]], color = species_colors[[sp]],
       rds = sdm_model_paths[[sp]])
})
names(species_config) <- species_list

## --- Readable variable labels + units ------------------------------------
var_labels <- c(
  sst = "Sea Surface Temperature", sst_sd = "SST Variability",
  ild = "Isothermal Layer Depth", EKE = "Eddy Kinetic Energy",
  ssh = "Sea Surface Height", ssh_sd = "SSH Variability",
  bath = "Depth", bath_sd = "Depth, SD", rugosity = "Rugosity",
  z = "Depth", z_sd = "Depth Variability", l.chl = "log(Chlorophyll)",
  l.chl_sd = "log(Chl) Variability", mld = "Mixed Layer\nDepth",
  bv = "Brunt-Väisälä\nFrequency", curl = "Wind Stress Curl",
  su = "Surface Eastward Current", sv = "Surface Northward Current",
  sustr = "Surface Eastward Wind Stress", svstr = "Surface Northward Wind Stress",
  lunar = "Lunar Illumination"
)
get_var_label <- function(v) ifelse(v %in% names(var_labels), var_labels[v], v)

var_units <- c(
  sst = "°C", ild = "m", EKE = "log m² s⁻²", ssh = "m",
  bath = "m", bath_sd = "m", z = "m", z_sd = "m",
  l.chl = "log mg m⁻³", mld = "m", curl = "N/m³",
  su = "m/s", sv = "m/s", sustr = "N/m²", svstr = "N/m²", lunar = "%"
)

## --- Core functions -------------------------------------------------------
extract_partial_dependence <- function(brt_models) {
  all_pd <- list()
  for (model_id in seq_along(brt_models)) {
    mod <- brt_models[[model_id]]
    for (var in mod$var.names) {
      pd <- as.data.frame(plot.gbm(mod, i.var = var, n.trees = mod$n.trees,
                                   return.grid = TRUE))
      colnames(pd) <- c("var_value", "predicted_HS")
      pd$variable <- var; pd$model <- model_id
      all_pd[[paste0(var, "_", model_id)]] <- pd
    }
  }
  do.call(rbind, all_pd)
}

extract_contributions <- function(brt_models) {
  contribs <- list()
  for (model_id in seq_along(brt_models)) {
    scores <- as.data.frame(summary.gbm(brt_models[[model_id]], plotit = FALSE))
    colnames(scores) <- c("variable", "contribution")
    scores$model_id <- model_id
    contribs[[model_id]] <- scores
  }
  do.call(rbind, contribs) %>%
    group_by(variable) %>%
    summarise(mean_contribution = mean(contribution, na.rm = TRUE),
              sd_contribution   = sd(contribution, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_contribution))
}

plot_partial_dependence <- function(pd_data, top_vars, sp_color, sp_code) {
  pd_subset <- pd_data %>% filter(variable %in% top_vars)
  pd_subset$variable <- factor(pd_subset$variable, levels = top_vars)
  facet_labels <- sapply(top_vars, function(v) {
    unit <- if (v %in% names(var_units)) paste0("\n(", var_units[v], ")") else ""
    paste0(get_var_label(v), unit)
  })
  names(facet_labels) <- top_vars
  ggplot() +
    geom_line(data = pd_subset, aes(x = var_value, y = predicted_HS, group = model),
              color = "grey70", linewidth = 0.6, alpha = 0.7) +
    geom_smooth(data = pd_subset, aes(x = var_value, y = predicted_HS),
                color = sp_color, linewidth = 1.2, method = "loess", se = FALSE) +
    labs(y = "Predicted Habitat Suitability", x = NULL) +
    facet_wrap(~ variable, scales = "free_x", strip.position = "top", nrow = 1,
               labeller = as_labeller(facet_labels)) +
    theme_minimal(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(size = 10, face = "bold"),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          axis.title.x = element_blank(), axis.text = element_text(size = 12),
          panel.grid.minor = element_blank(), axis.title = element_text(size = 12),
          plot.margin = margin(5, 10, 5, 10))
}

## --- Generate per species -------------------------------------------------
sp_plots <- list(); sp_contribs <- list()
for (sp in names(species_config)) {
  cfg <- species_config[[sp]]
  message("Processing: ", cfg$name, " (", sp, ")")
  brt_models <- readRDS(cfg$rds)
  message("  Loaded ", length(brt_models), " BRT models")

  contribs <- extract_contributions(brt_models)
  pd_data  <- extract_partial_dependence(brt_models)
  top3     <- head(contribs$variable, 3)
  message("  Top 3 predictors: ", paste(top3, collapse = ", "))

  sp_contribs[[sp]] <- contribs %>% mutate(species = sp)
  sp_plots[[sp]]    <- plot_partial_dependence(pd_data, top3, cfg$color, sp)
}

## --- Supplementary Table 1: variable contributions -----------------------
all_contribs <- do.call(rbind, sp_contribs) %>%
  mutate(species_name   = sapply(species, function(sp) species_config[[sp]]$name),
         variable_label = sapply(variable, get_var_label),
         mean_contribution = round(mean_contribution, 2),
         sd_contribution   = round(sd_contribution, 2)) %>%
  select(species_name, variable_label, mean_contribution, sd_contribution) %>%
  arrange(species_name, desc(mean_contribution))

write.csv(all_contribs,
          file.path(output_plots_dir,
                    paste0("suppTable1_variable_contributions_", Sys.Date(), ".csv")),
          row.names = FALSE)
message("Supplementary Table 1 saved.")

## --- Assemble Supplementary Figure 1 -------------------------------------
add_species_label <- function(sp) {
  cfg <- species_config[[sp]]
  img_path <- file.path(output_plots_dir, "images", paste0(toupper(sp), "_img.jpg"))
  if (file.exists(img_path)) {
    ggdraw() + draw_image(img_path, scale = 0.85) +
      draw_label(cfg$name, size = 9, hjust = 0.5, vjust = 6.5, fontface = "bold") +
      theme(plot.margin = margin(0, 0, 0, 0))
  } else {
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = cfg$name,
                        size = 4, fontface = "bold") + theme_void()
  }
}

panel_labels <- c("A", "B", "C", "D")
rows <- lapply(seq_along(species_list), function(i) {
  sp <- species_list[i]
  add_species_label(sp) +
    (sp_plots[[sp]] + ggtitle(paste0(panel_labels[i], ".")) +
       theme(plot.title = element_text(size = 13, face = "bold", hjust = 0))) +
    plot_layout(widths = c(0.18, 1))
})

supp_fig1 <- wrap_plots(rows, ncol = 1) +
  plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))

ggsave(filename = file.path(output_plots_dir,
                            paste0("suppFig1_partial_dependence_", Sys.Date(), ".png")),
       plot = supp_fig1, width = 10, height = 12, units = "in", dpi = 600, bg = "white")
message("Supplementary Figure 1 saved.")
