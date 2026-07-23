# Data availability

This repository contains **code only**. The input datasets are archived
externally and are described below. Set `data_root` in [`config.R`](config.R)
(or the `TE_CCE_DATA_ROOT` environment variable) to the folder where you have
downloaded them.

## Input data

### 1. Species habitat-suitability projections (SDM output)
Projected habitat suitability for the four focal species (California sea lion,
humpback whale, leatherback sea turtle, swordfish) from published boosted
regression tree species distribution models, forced by three downscaled ESMs
under RCP 8.5.

- Highly mobile species projections:
  <https://oceanview.pfeg.noaa.gov/erddap/files/HMS_habitat_suitability/>
- Source SDMs: Becker et al. (2020); Brodie et al. (2018);
  Lezama-Ochoa et al. (2025); Welch et al. (2023).

The daily habitat-suitability rasters (`{species}*mean.grd`) are the input to
`03_composite_sdm_output.R`.

### 2. Boosted regression tree model objects (`.rds`)
The fitted 10-model BRT ensembles are required only for the partial-dependence
plots (Supp. Fig. 1) and variable-contribution table (Supp. Table 1),
produced by `08_partial_dependence.R`. Paths are set in `config.R`
(`sdm_model_paths`).

### 3. Physical ocean projections (ROMS / ESM SST)
Regionally downscaled ROMS output and the CMIP5 ESM SST fields
(GFDL-ESM2M, IPSL-CM5A-MR, HadGEM2-ES) used to identify thermal extremes.

- ROMS physical variables:
  <https://oceanview.pfeg.noaa.gov/erddap/search/index.html?searchFor=roms>
- ESM SST are the raw input to `01_preprocess_sst.sh`
  (Pozo Buil et al., 2021).

### 4. Supporting spatial files
- `lme-mask-out.nc` — Large Marine Ecosystem mask (California Current = region 3).
- Grid weights and the cropped mask are produced by `01_preprocess_sst.sh`.

## Expected local layout

```
<data_root>/
├── output data/
│   ├── GFDL/  IPSL/  Had/          # raw + preprocessed SST NetCDFs per ESM
│   └── extremes_ccs_ssta/
│       └── composites/             # composite .grd (from script 03)
├── output plots/                   # figures are written here
├── shared code/                    # lme-mask-out-cropped.nc, gridweights.nc
└── EcoROMS models/                 # BRT .rds ensembles
```

> Note: `01_preprocess_sst.sh` writes `gridweights.nc` and
> `lme-mask-out-cropped.nc` into `output data/`; copy them into `shared code/`
> (or adjust `mask_file` / `grid_weights_file` in `config.R`) before running
> `02_identify_extremes.R`.
