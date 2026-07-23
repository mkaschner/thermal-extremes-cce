# Thermal extremes shape future highly mobile species redistribution in the California Current

<!-- DOI badge: add after publishing the reserved Zenodo record -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

Code accompanying the manuscript:

> Kaschner, M. L., Lezama-Ochoa, N., Welch, H., Cluett, A., Pozo Buil, M.,
> Bograd, S. J., Jacox, M. G., Hazen, E. L., & Barton, A. D. *Thermal extremes
> shape future highly mobile species redistribution in the California Current.*
> (planning to submit to *Diversity and Distributions*).

This repository reproduces the analysis and figures. It assesses how habitat
suitability for four highly mobile marine species — California sea lion
(`casl`), humpback whale (`hbwh`), leatherback sea turtle (`lbst`), and
swordfish (`swor`) — responds to long-term warming versus transient thermal
extremes in the historical (1991–2020) and future (2071–2100) periods, using
published species distribution models forced by three downscaled Earth System
Models (GFDL, IPSL, HadGEM2) under RCP 8.5.

## Repository layout

```
thermal-extremes-cce/
├── config.R                       # paths + study constants (EDIT THIS FIRST)
├── R/functions.R                  # shared helper functions
├── scripts/
│   ├── 01_preprocess_sst.sh       # CDO: regrid, anomalies, crop  (S2.3.1)
│   ├── 02_identify_extremes.R     # regional SSTa + extreme months (S2.3.2) + Supp Fig 2
│   ├── 03_composite_sdm_output.R  # server-side: composite daily SDM output
│   ├── 04_centroid_area_metrics.R # core-habitat centroids & area  (S2.3.3)
│   ├── 05_habitat_maps.R          # Figures 2–4 (difference maps)
│   ├── 06_figure_summary.R        # Figure 5 (distance | direction | area)
│   └── 07_partial_dependence.R    # Supp Fig 1 + Supp Table 1 (gbm)
├── DATA.md                        # data availability / where to get inputs
└── README.md
```

Figure 1 is a conceptual schematic of the study design and hypotheses,
assembled externally (its SST-time-series element derives from the outputs of
`02_identify_extremes.R`); it is not produced by a standalone script here.

## Requirements

- **R 4.5.1** with: `ncdf4`, `raster`, `terra`, `sf`, `tidyverse` (`dplyr`,
  `tidyr`, `ggplot2`, `purrr`), `rnaturalearth`, `rnaturalearthdata`,
  `patchwork`, `cowplot`, `geosphere`, `gbm`, `viridis`, `scales`, `scico`,
  `ragg`, `zoo`.
- **CDO 2.6.0** (Climate Data Operators) for `01_preprocess_sst.sh`.

## How to run

1. Obtain the input data (see [DATA.md](DATA.md)) and set `data_root` in
   [`config.R`](config.R).
2. Run the pipeline in order. R scripts are run from the repository root
   (they `source("config.R")`):

| Step | Script | Produces |
|------|--------|----------|
| 1 | `bash scripts/01_preprocess_sst.sh` | SST anomalies, grid weights, LME mask |
| 2 | `Rscript scripts/02_identify_extremes.R` | regional SSTa, extreme-month CSVs, **Supp Fig 2** |
| 3 | `Rscript scripts/03_composite_sdm_output.R` | composite habitat-suitability `.grd` *(run on SDM server)* |
| 4 | `Rscript scripts/04_centroid_area_metrics.R` | `centroid_results.csv` |
| 5 | `Rscript scripts/05_habitat_maps.R` | **Figures 2–4** |
| 6 | `Rscript scripts/06_figure_summary.R` | **Figure 5** |
| 7 | `Rscript scripts/07_partial_dependence.R` | **Supp Fig 1**, **Supp Table 1** |

Steps 5–7 are independent of one another once steps 1–4 have produced the
composites and `centroid_results.csv`. Figure 5 panel placement was refined
externally (PowerPoint) after export.

## Notes on terminology

Internally the code labels the two thermal extremes `amplify` (warm; 90–100th
percentile SST-anomaly months) and `minimize` (cool; 0–10th percentile). The
three manuscript comparisons are:

| Code | Manuscript term | Figure |
|------|-----------------|--------|
| `F[ave] − H[ave]` | Anthropogenic (long-term) warming | Fig 2 |
| `H[amp] − H[min]` | Historical thermal extremes | Fig 3 |
| `F[amp] − F[min]` | Future thermal extremes | Fig 4 |

Multi-panel figures (e.g. the 4-species × 3-column layouts of Figs 2–4) are
assembled from the individual panels that the scripts export.

## License

Released under the [MIT License](LICENSE).

## Citation

If you use this code, please cite the article and the archived software release
(see [CITATION.cff](CITATION.cff)):

> Kaschner, M. L. et al. (2026). Analysis code for "Thermal extremes shape
> future highly mobile species redistribution in the California Current".
> Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

_Replace `XXXXXXX` with the reserved Zenodo DOI once the record is published._
