#!/bin/bash
## =========================================================================
## 01_preprocess_sst.sh — SST preprocessing with CDO (manuscript S2.3.1)
## =========================================================================
## For each ESM (GFDL, IPSL, HadGEM2), this pipeline:
##   (a) Regrids historical (1976-2005) and future/RCP8.5 (2006-2100) monthly
##       SST to a common 1 deg x 1 deg grid via bilinear interpolation.
##   (b) Merges the two into one continuous 1976-2100 time series.
##   (c) Computes monthly anomalies relative to the 1991-2020 climatology.
##   (d) Crops to the study region bounding box (150-260E, 20-60N).
## It also produces the shared grid-cell area weights and the cropped
## California Current LME mask used in 02_identify_extremes.R.
##
## Requirements: CDO 2.6.0 (Climate Data Operators).
## Input NetCDFs are described in DATA.md and are NOT included in this repo.
##
## Usage:
##   OUTPUT_DATA_DIR=/path/to/output_data bash 01_preprocess_sst.sh
## The output data directory must contain one subfolder per model (GFDL/,
## IPSL/, Had/) holding that model's raw SST NetCDFs, and an `lme-mask-out.nc`.
## =========================================================================

set -euo pipefail

# Directory holding per-model SST NetCDFs and the LME mask (edit or override).
OUTPUT_DATA_DIR="${OUTPUT_DATA_DIR:-../output data}"
cd "$OUTPUT_DATA_DIR" || { echo "Cannot cd to $OUTPUT_DATA_DIR"; exit 1; }

# ==============================================================================
# 1. REUSABLE PIPELINE (anomalies & cropping)
# ==============================================================================
run_analysis() {
    local PREFIX=$1
    local INPUT=$2
    local MODEL=$3
    local ANOM_CMD=$4

    echo "---- Running ${PREFIX} pipeline for ${MODEL} ----"

    # Monthly anomalies relative to the 1991-2020 monthly climatology
    eval "cdo -L -ymonsub ${INPUT} ${ANOM_CMD} ${PREFIX}_anoms_${MODEL}.nc"

    # Crop to study region bounding box
    cdo sellonlatbox,150,260,20,60 ${PREFIX}_anoms_${MODEL}.nc ${PREFIX}_anoms_cropped_${MODEL}.nc
}

# ==============================================================================
# 2. LOOP THROUGH MODELS
# ==============================================================================
MODELS=("GFDL" "IPSL" "Had")

for MODEL in "${MODELS[@]}"; do
    cd "$MODEL" || continue

    if [ "$MODEL" == "GFDL" ]; then
        M_STR="GFDL-ESM2M"; ENS="r1i1p1"
    elif [ "$MODEL" == "IPSL" ]; then
        M_STR="IPSL-CM5A-MR"; ENS="r1i1p1"
    elif [ "$MODEL" == "Had" ]; then
        M_STR="HadGEM2-ES"; ENS="r2i1p1"
    fi

    # File names
    HIST_RAW="sst_Omon_${M_STR}_historical_${ENS}_1976-2005.nc"
    FUT_RAW="sst_Omon_${M_STR}_rcp85_${ENS}_2006-2100.nc"
    HIST_REMAP="remap_sst_Omon_${M_STR}_historical_${ENS}_1976-2005.nc"
    FUT_REMAP="remap_sst_Omon_${M_STR}_rcp85_${ENS}_2006-2100.nc"
    FULL_REMAP="remap_sst_Omon_${M_STR}_full_1976-2100.nc"

    # Step A: Remap to 1-degree grid (bilinear)
    cdo remapbil,global_1 "$HIST_RAW" "$HIST_REMAP"
    cdo remapbil,global_1 "$FUT_RAW" "$FUT_REMAP"

    # Step B: Merge historical + future into the continuous "full" series
    cdo mergetime "$HIST_REMAP" "$FUT_REMAP" "$FULL_REMAP"

    # Step C: Anomalies (vs 1991-2020 monthly climatology) + crop
    run_analysis "full" "$FULL_REMAP" "$MODEL" "-ymonmean -selyear,1991/2020 ${FULL_REMAP}"

    cd ..
done

# ==============================================================================
# 3. SHARED GRID WEIGHTS & LME MASK
# ==============================================================================
cdo gridweights GFDL/remap_sst_Omon_GFDL-ESM2M_full_1976-2100.nc gridweights.nc
cdo sellonlatbox,150,260,20,60 lme-mask-out.nc lme-mask-out-cropped.nc
cdo gridweights lme-mask-out-cropped.nc gridweights-cropped.nc

echo "Preprocessing complete."
