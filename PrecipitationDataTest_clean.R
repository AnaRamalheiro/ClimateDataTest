###############################################################################
# CLIMATE DATA TEST
# Bias adjustment and example threshold analysis of precipitation projections
###############################################################################

# Purpose:
# This script is a methodological test developed to understand the workflow
# required to analyse climate-model precipitation projections.
#
# The analysis:
#   1. Reads historical and projected precipitation from IPSL-CM6A-LR
#   2. Uses ERA5-Land as the reference dataset
#   3. Evaluates historical precipitation bias
#   4. Tests Quantile Delta Mapping (QDM) using calibration and validation periods
#   5. Applies QDM to SSP5-8.5 precipitation projections
#   6. Tests a hypothetical precipitation indicator and threshold
#
# IMPORTANT:
# The precipitation indicator (>15 mm/day for 5 consecutive days) and its
# threshold are hypothetical and are used only to test the methodology.
#
# This is a single-GCM methodological exercise and does not include spatial
# downscaling.


# ---------------------------------------------------------------------------
# 1. SETUP
# ---------------------------------------------------------------------------

rm(list = ls())

library(ncdf4)
library(ggplot2)
library(dplyr)
library(MBC)
library(ggrepel)


# ---------------------------------------------------------------------------
# 2. LOAD CLIMATE DATA
# ---------------------------------------------------------------------------

# Historical IPSL-CM6A-LR precipitation (1981-2014)
hist_file <- "CMIP6_Historical_1981_2014/pr_day_IPSL-CM6A-LR_historical_r1i1p1f1_gr_19810101-20141231.nc"

# ERA5-Land precipitation
era_file <- "ERA5 Total Precipitation/reanalysis-era5-land-timeseries-sfc-pressure-precipitation3rt9rcwn.nc"

# IPSL-CM6A-LR SSP5-8.5 precipitation projection (2015-2030)
ssp_file <- "CMIP6_SSP5-8.5/pr_day_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr_20150101-20301231.nc"


# Open NetCDF files
nc_hist <- nc_open(hist_file)
nc_era  <- nc_open(era_file)
nc_ssp  <- nc_open(ssp_file)


# ---------------------------------------------------------------------------
# 3. PROCESS HISTORICAL IPSL DATA
# ---------------------------------------------------------------------------

# Extract daily precipitation.
# CMIP6 precipitation (pr) is expressed in kg m-2 s-1.
pr_hist <- ncvar_get(nc_hist, "pr")

# Convert kg m-2 s-1 to mm/day.
# 1 kg m-2 water = 1 mm water, therefore multiply by 86,400 seconds/day.
pr_hist_mm <- pr_hist * 86400


# Convert longitude from 0-360 to -180-180.
lon_hist <- nc_hist$dim$lon$vals
lon_hist[lon_hist > 180] <- lon_hist[lon_hist > 180] - 360

lat_hist <- nc_hist$dim$lat$vals


# Select the IPSL grid cell closest to the Ria de Aveiro reference point
# (approximately 40.6 N, -8.7 E).
#
# The closest available IPSL grid cell is:
# latitude  = 40.56338
# longitude = -7.5
#
# Note: this coarse GCM grid cell does not coincide exactly with the
# ERA5-Land reference point.
hist_point <- pr_hist_mm[2, 2, ]


# ---------------------------------------------------------------------------
# 4. PROCESS SSP5-8.5 IPSL DATA
# ---------------------------------------------------------------------------

pr_ssp <- ncvar_get(nc_ssp, "pr")

# Convert kg m-2 s-1 to mm/day.
pr_ssp_mm <- pr_ssp * 86400

lat_ssp <- nc_ssp$dim$lat$vals
lon_ssp <- nc_ssp$dim$lon$vals
lon_ssp[lon_ssp > 180] <- lon_ssp[lon_ssp > 180] - 360


# Extract the same IPSL grid cell used for the historical period.
ssp_point <- pr_ssp_mm[2, 2, ]


# ---------------------------------------------------------------------------
# 5. PROCESS ERA5-LAND DATA
# ---------------------------------------------------------------------------

# Extract ERA5-Land time information.
era_time <- nc_era$dim$valid_time$vals

era_dates <- as.POSIXct(
  era_time * 3600,
  origin = "1970-01-01",
  tz = "UTC"
)


# Extract total precipitation (tp), provided in metres.
tp_era <- ncvar_get(nc_era, "tp")

# Convert metres to millimetres.
tp_era_mm <- tp_era * 1000


era_df <- data.frame(
  datetime = era_dates,
  precip_mm = tp_era_mm
)


# Convert hourly ERA5-Land precipitation to daily precipitation totals.
era_df$date <- as.Date(era_df$datetime)

era_daily <- era_df %>%
  group_by(date) %>%
  summarise(
    precip_mm = sum(precip_mm, na.rm = TRUE),
    .groups = "drop"
  )


# ---------------------------------------------------------------------------
# 6. MATCH HISTORICAL PERIODS
# ---------------------------------------------------------------------------

# Restrict ERA5-Land to the IPSL historical period (1981-2014).
era_hist <- era_daily %>%
  filter(
    date >= as.Date("1981-01-01"),
    date <= as.Date("2014-12-31")
  )


# Convert IPSL historical time to calendar dates.
hist_time <- nc_hist$dim$time$vals

hist_datetime <- as.POSIXct(
  hist_time * 86400,
  origin = "1850-01-01",
  tz = "UTC"
)

hist_dates <- as.Date(hist_datetime)


# Check that ERA5-Land and IPSL dates correspond.
all(hist_dates == era_hist$date)


# Combine historical ERA5-Land and IPSL precipitation.
comparison_hist <- data.frame(
  date = hist_dates,
  era_mm = era_hist$precip_mm,
  ipsl_mm = hist_point
)


# Remove very small negative ERA5-Land precipitation values caused by
# numerical/encoding effects.
comparison_hist$era_mm[comparison_hist$era_mm < 0] <- 0


# ---------------------------------------------------------------------------
# 7. HISTORICAL MODEL EVALUATION
# ---------------------------------------------------------------------------

# 7.1 Annual precipitation --------------------------------------------------

annual_hist <- comparison_hist %>%
  mutate(year = format(date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    era_total = sum(era_mm, na.rm = TRUE),
    ipsl_total = sum(ipsl_mm, na.rm = TRUE),
    .groups = "drop"
  )


# Mean absolute annual bias (mm/year).
mean_bias <- mean(annual_hist$ipsl_total) -
  mean(annual_hist$era_total)

mean_bias


# Mean relative annual bias (%).
relative_bias <- (
  mean(annual_hist$ipsl_total) /
    mean(annual_hist$era_total) - 1
) * 100

relative_bias


# Plot annual precipitation.
annual_long <- annual_hist %>%
  tidyr::pivot_longer(
    cols = c(era_total, ipsl_total),
    names_to = "dataset",
    values_to = "precip_mm"
  )


ggplot(
  annual_long,
  aes(x = as.numeric(year), y = precip_mm, color = dataset)
) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Annual precipitation (1981-2014)",
    x = "Year",
    y = "Annual precipitation (mm)",
    color = "Dataset"
  ) +
  theme_minimal()


# 7.2 Monthly precipitation climatology ------------------------------------

comparison_hist <- comparison_hist %>%
  mutate(
    year = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m"))
  )


monthly_hist <- comparison_hist %>%
  group_by(year, month) %>%
  summarise(
    era_total = sum(era_mm, na.rm = TRUE),
    ipsl_total = sum(ipsl_mm, na.rm = TRUE),
    .groups = "drop"
  )


monthly_climatology <- monthly_hist %>%
  group_by(month) %>%
  summarise(
    era_mean = mean(era_total),
    ipsl_mean = mean(ipsl_total),
    .groups = "drop"
  )


monthly_plot <- monthly_climatology %>%
  tidyr::pivot_longer(
    cols = c(era_mean, ipsl_mean),
    names_to = "dataset",
    values_to = "precip_mm"
  )


ggplot(
  monthly_plot,
  aes(x = month, y = precip_mm, color = dataset)
) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  labs(
    title = "Mean monthly precipitation (1981-2014)",
    x = "Month",
    y = "Mean monthly precipitation (mm)",
    color = "Dataset"
  ) +
  theme_minimal()


# 7.3 Wet-day frequency and intensity --------------------------------------

# For this diagnostic, a wet day is defined as precipitation >= 1 mm/day.
comparison_hist$era_wet <- comparison_hist$era_mm >= 1
comparison_hist$ipsl_wet <- comparison_hist$ipsl_mm >= 1


# Mean number of wet days per year.
sum(comparison_hist$era_wet) / length(unique(comparison_hist$year))
sum(comparison_hist$ipsl_wet) / length(unique(comparison_hist$year))


# Mean precipitation intensity on wet days.
mean(comparison_hist$era_mm[comparison_hist$era_wet])
mean(comparison_hist$ipsl_mm[comparison_hist$ipsl_wet])


# Median precipitation intensity on wet days.
median(comparison_hist$era_mm[comparison_hist$era_wet])
median(comparison_hist$ipsl_mm[comparison_hist$ipsl_wet])


# 95th percentile of wet-day precipitation.
quantile(
  comparison_hist$era_mm[comparison_hist$era_wet],
  probs = 0.95
)

quantile(
  comparison_hist$ipsl_mm[comparison_hist$ipsl_wet],
  probs = 0.95
)


# Historical evaluation indicates differences between IPSL-CM6A-LR and
# ERA5-Land in precipitation amount, wet-day frequency and wet-day intensity.


# ---------------------------------------------------------------------------
# 8. TEST QDM USING CALIBRATION AND VALIDATION PERIODS
# ---------------------------------------------------------------------------

# Split historical data:
# calibration = 1981-1999
# validation  = 2000-2014
#
# The validation period is kept separate to evaluate QDM on data that were
# not used to estimate the bias adjustment.

calibration <- comparison_hist[comparison_hist$year < 2000, ]
validation  <- comparison_hist[comparison_hist$year >= 2000, ]


era_cal  <- calibration$era_mm
ipsl_cal <- calibration$ipsl_mm
ipsl_val <- validation$ipsl_mm


# Apply Quantile Delta Mapping.
# ratio = TRUE is used because precipitation is treated multiplicatively.
qdm_val <- QDM(
  o.c = era_cal,
  m.c = ipsl_cal,
  m.p = ipsl_val,
  ratio = TRUE
)


# Extract bias-adjusted IPSL validation precipitation.
ipsl_val_corrected <- qdm_val$mhat.p


validation_qdm <- data.frame(
  date = validation$date,
  ERA = validation$era_mm,
  IPSL_raw = ipsl_val,
  IPSL_corrected = ipsl_val_corrected
)


# ---------------------------------------------------------------------------
# 9. EVALUATE QDM PERFORMANCE
# ---------------------------------------------------------------------------

# Compare mean daily precipitation.
mean(validation_qdm$ERA)
mean(validation_qdm$IPSL_raw)
mean(validation_qdm$IPSL_corrected)


# Compare wet-day frequency (>=1 mm/day).
validation_qdm$ERA_wet <- validation_qdm$ERA >= 1
validation_qdm$IPSL_raw_wet <- validation_qdm$IPSL_raw >= 1
validation_qdm$IPSL_corrected_wet <- validation_qdm$IPSL_corrected >= 1

sum(validation_qdm$ERA_wet)
sum(validation_qdm$IPSL_raw_wet)
sum(validation_qdm$IPSL_corrected_wet)


# Compare the 95th percentile of wet-day precipitation.
quantile(
  validation_qdm$ERA[validation_qdm$ERA_wet],
  probs = 0.95
)

quantile(
  validation_qdm$IPSL_raw[validation_qdm$IPSL_raw_wet],
  probs = 0.95
)

quantile(
  validation_qdm$IPSL_corrected[validation_qdm$IPSL_corrected_wet],
  probs = 0.95
)


# Check QDM behaviour during the calibration period.
quantile(
  calibration$era_mm[calibration$era_wet],
  probs = 0.95
)

quantile(
  calibration$ipsl_mm[calibration$ipsl_wet],
  probs = 0.95
)

quantile(
  qdm_val$mhat.c[qdm_val$mhat.c >= 1],
  probs = 0.95
)


# QDM closely reproduces the reference P95 during calibration.
# Performance differs during validation because the changes in precipitation
# distributions between the calibration and validation periods differ between
# ERA5-Land and IPSL. QDM corrects systematic distributional biases while
# preserving modelled changes; it cannot force the model to reproduce the
# temporal evolution observed in the reference dataset.


# ---------------------------------------------------------------------------
# 10. APPLY QDM TO SSP5-8.5 PROJECTION
# ---------------------------------------------------------------------------

# After testing QDM using independent calibration and validation periods,
# use the complete historical period (1981-2014) to estimate the adjustment
# applied to the future SSP5-8.5 projection.

ipsl_hist <- comparison_hist$ipsl_mm


qdm_hist <- QDM(
  o.c = comparison_hist$era_mm,
  m.c = ipsl_hist,
  m.p = ssp_point,
  ratio = TRUE
)


ssp5_corrected <- qdm_hist$mhat.p


# Create dates for SSP5-8.5 projection.
ssp_time <- nc_ssp$dim$time$vals
ssp_dates <- as.Date("2015-01-01") + ssp_time


ssp5_qdm <- data.frame(
  Date = ssp_dates,
  IPSL_SSP5_raw = ssp_point,
  IPSL_SSP5_corrected = ssp5_corrected
)


# Compare raw and bias-adjusted projected precipitation.
summary(ssp5_qdm$IPSL_SSP5_raw)
summary(ssp5_qdm$IPSL_SSP5_corrected)


# ---------------------------------------------------------------------------
# 11. HYPOTHETICAL THRESHOLD ANALYSIS
# ---------------------------------------------------------------------------

# Example indicator:
# occurrence of >=5 consecutive days with precipitation >15 mm/day.
#
# IMPORTANT: both 15 mm/day and 5 consecutive days are hypothetical values
# used solely to test the analysis workflow. They should not be interpreted
# as ecological thresholds for the Ria de Aveiro salt marshes.


ssp5_qdm$Wet_days <- ssp5_qdm$IPSL_SSP5_corrected > 15


# Identify consecutive runs of TRUE/FALSE values.
wet_runs <- rle(ssp5_qdm$Wet_days)


# Identify runs containing at least 5 consecutive high-precipitation days.
precip_threshold <- wet_runs$values == TRUE &
  wet_runs$lengths >= 5


# Number of threshold events.
sum(precip_threshold)


# Duration of each threshold event.
wet_runs$lengths[precip_threshold]


# ---------------------------------------------------------------------------
# 12. CALCULATE CONSECUTIVE HIGH-PRECIPITATION DAYS
# ---------------------------------------------------------------------------

# Create a counter that resets to zero whenever daily precipitation falls
# to <=15 mm/day.

ssp5_qdm$Consecutive_wet_days <- 0


# Treat the first observation separately because it has no preceding day.
if (ssp5_qdm$Wet_days[1] == TRUE) {
  ssp5_qdm$Consecutive_wet_days[1] <- 1
}


# For each subsequent day:
# - if precipitation >15 mm/day, add 1 to the previous day's counter;
# - otherwise, reset the counter to zero.
for (i in 2:nrow(ssp5_qdm)) {
  
  if (ssp5_qdm$Wet_days[i] == TRUE) {
    
    ssp5_qdm$Consecutive_wet_days[i] <-
      1 + ssp5_qdm$Consecutive_wet_days[i - 1]
    
  } else {
    
    ssp5_qdm$Consecutive_wet_days[i] <- 0
  }
}


# Check the longest projected sequence.
max(ssp5_qdm$Consecutive_wet_days)


# Identify the dates on which the hypothetical 5-day threshold is reached.
threshold_events <- ssp5_qdm[
  ssp5_qdm$Consecutive_wet_days == 5,
]


threshold_events


# ---------------------------------------------------------------------------
# 13. PLOT HYPOTHETICAL THRESHOLD EVENTS
# ---------------------------------------------------------------------------

ggplot(
  data = ssp5_qdm,
  aes(x = Date, y = Consecutive_wet_days)
) +
  geom_point() +
  geom_hline(
    aes(
      yintercept = 5,
      linetype = "Hypothetical threshold (5 consecutive days)"
    ),
    colour = "red",
    linewidth = 0.8
  ) +
  geom_point(
    data = threshold_events
  ) +
  geom_text_repel(
    data = threshold_events,
    aes(label = Date),
    vjust = 1
  ) +
  scale_y_continuous(
    breaks = 0:5
  ) +
  labs(
    title = "Projected consecutive high precipitation (>15 mm) events",
    subtitle = "IPSL-CM6A-LR, SSP5-8.5 (2015-2030)",
    x = "Year",
    y = "Number of consecutive days",
    linetype = NULL
  )


# ---------------------------------------------------------------------------
# END OF TEST
# ---------------------------------------------------------------------------
