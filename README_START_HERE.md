# Supplementary R code: 1981–2024 Farming Systems Trial yield stability

Manuscript: **Crop-specific grain yield stability across 44 years of organic and conventional cropping systems**
Journal: *Farming System*

## Contents

- `01_Primary_Model/`: Nine R source files for the annual-yield Gaussian REML, weather effects, model checks, 1,999-replicate primary parametric bootstrap, omnibus tests, matched-year FT/RT comparisons, and model figures. The file `07_FIGURES_AND_MANUSCRIPT_TABLES.R` contains the previously prepared correction that converts weather-regression effects to **one observed-reference SD** before labelling them as such. The unpatched file is preserved in `04_Archived_Unpatched_Original/` for provenance only.
- `02_Manuscript_Figures/`: Four modular figure scripts, a standalone figure script, and 12 small derived plotting CSVs. The plotting datasets are included to reproduce descriptive figures without re-distributing the raw plot-level dataset. The archived REML point estimates in `Derived_data/12_primary_model_point_estimates.csv` should be refreshed from the primary model if that model is re-estimated.
- `03_SPI3_Drought_Figures/`: R code and three input CSVs for the descriptive SPI-3 drought-year/system-yield figure and supplementary drought-year wheat tillage comparison. These figures are descriptive, not inferential tests of climate resilience.
- `05_Notes/FIGURE_CAPTIONS.md`: Reference captions from the figure-authoring stage. Reconcile figure numbering with the final manuscript.

## Required original data (not included)

The primary analysis requires a source Excel workbook with tabs `corn_81_24`, `soy_81_24`, `wheat_81_24` and the site-specific PRISM daily precipitation/temperature CSV.

Place the two original raw input files at **the root of this extracted package**, i.e. immediately above `01_Primary_Model/` and `02_Manuscript_Figures/`.

For the original statistical pipeline, name the yield workbook exactly `FST_1981_2024_3_crops_BG.xlsx`; if your upload is named `FST_1981_2024_3_crops_BG(1).xlsx`, rename it (or edit `01_Primary_Model/01_SETTINGS.R`). The PRISM name expected is the current `WEATHER_FILE` in `01_SETTINGS.R`.

## Running the primary model

Requires R >= 4.2 and packages `readxl` and `ggplot2`. Execute from the extracted package, for example:

```
cd 01_Primary_Model
Rscript 00_RUN_EVERYTHING_V9.R
```

The initial settings request `RUN_MODE <- 'FULL'` with 1,999 primary bootstrap runs and 999 tillage bootstrap runs; use `SMOKE` for a preliminary functional check before obtaining final results. Outputs are saved to `FST_V9_FINAL_ANALYSIS_outputs` in the input data directory. This code is the original reproducibility model from the authors' files, with the weather unit correction noted above.

## Running manuscript figures

```
cd 02_Manuscript_Figures
Rscript 00_RUN_ALL_FIGURES.R --plots-only
```

The `--plots-only` option draws graphics using the supplied Derived_data CSV files and requires only R/ggplot2. To recalculate the derived datasets, run without `--plots-only` after placing the original input files in the package root; requires `readxl`. Alternatively, run `Rscript FST_All_Figures_Standalone.R` from within `02_Manuscript_Figures`; this script automatically searches its own and parent directory for the raw files.

## Running SPI-3 figure

```
cd 03_SPI3_Drought_Figures
Rscript plot_SPI3_drought_FST.R
```

Requires R packages `ggplot2` and `patchwork`. Outputs PNG, TIFF, and SVG files.

## Important methodological boundaries

- A crop absent from a rotation is **unobserved**, not a zero-yield event.
- Annual means average entry points within replicate before equal replicate weighting.
- Long-term series combines FT before 2008 and equal FT/RT arm means after 2008 only when both treatments were present.
- Yield gaps compare systems over **matched harvest years**, not total rotation output.
- The original model's 2008 breakpoint describes a management transition; it does not isolate a causal effect of reduced tillage.
- The corn stability contrast is sensitive to the covariance model (B versus UN). Supplementary sensitivity results should be considered with the primary estimates.
- The drought-year visualizations classify meteorological precipitation deficits with SPI-3; they do not directly quantify crop physiological drought stress.

## Submission and verification

The archive provides the author-supplied R programs and derived plotting inputs. **The full R pipeline has not been executed in this authoring environment.** The recipient should run the scripts in R, inspect generated graphics and output tables, and verify that all manuscript statistics agree with the final model outputs before public deposition.

Input data are **not** bundled; raw-yield access and PRISM redistribution should be handled separately under the authors' applicable data-sharing policy. Some R filenames still contain an internal historical version suffix; these suffixes are code-provenance identifiers, not manuscript version labels.
