# Configuration for the FST grain yield and stability analysis.
# Input files are located in the analysis data directory.
ROOT_OVERRIDE <- NULL  # Parent data directory is detected when this is NULL
YIELD_FILE <- 'FST_1981_2024_3_crops_BG.xlsx'
WEATHER_FILE <- 'PRISM_ppt_tmin_tmean_tmax_tdmean_soltotal_Unlikely_to_change_800m_19810101_20260101_40.5509_-75.7301.csv'
OUTPUT_FOLDER <- 'FST_V9_FINAL_ANALYSIS_outputs'
YEAR_RANGE <- 1981:2024
ERA_BREAK <- 2008L
CROPS <- c('corn','soybean','wheat')
SYSTEMS <- c('CNV','LEG','MNR')
# Crop-specific trend and covariance specifications.
PRIMARY <- list(
  corn=list(trend='M1_STEP',covariance='B',serial='AR1'),
  soybean=list(trend='M4_PIECEWISE',covariance='B',serial='NONE'),
  wheat=list(trend='M2_LINEAR',covariance='B',serial='NONE'))
# Bootstrap replication counts for diagnostic and full runs.
RUN_MODE <- 'FULL'  # Diagnostic run: SMOKE; full analysis: FULL
SEED <- 20261112L
MAIN_BOOT_SMOKE <- 49L
MAIN_BOOT_FULL <- 1999L
OMNIBUS_BOOT_SMOKE <- 19L
OMNIBUS_BOOT_FULL <- 499L
TILLAGE_BOOT_SMOKE <- 39L
TILLAGE_BOOT_FULL <- 999L
CHECKPOINT_EVERY <- 25L
# Optional analyses and output generation.
RUN_SENSITIVITY <- TRUE
RUN_LEAVE_ONE_YEAR_OUT <- TRUE
RUN_POST2008_TILLAGE <- TRUE
RUN_BOOTSTRAP <- TRUE
# Reference weather uses calendar years with a harvest in any system.
# Mean yields are standardized across the full 1981–2024 period.
REFERENCE_WEATHER <- 'OBSERVED_CROP_YEARS'
# Weather-reference sensitivity analysis.
ALTERNATIVE_REFERENCE <- 'ALL_ELIGIBLE_WEATHER_YEARS'
# Weather indices are derived from PRISM precipitation and temperature records.
HOT_DAY_C <- 32
SENSITIVITY_BOOT_SMOKE <- 49L
SENSITIVITY_BOOT_FULL <- 999L
# Additional bootstrap analysis for the corn covariance specification.
RUN_CORN_UN_BOOTSTRAP <- TRUE
# Confidence level for simultaneous within-crop system comparisons.
FAMILYWISE_CONFIDENCE <- 0.95
