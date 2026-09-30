
# Precipitation data analysis ---------------------------------------------
## This is an example of a threshold analysis done on precipitation data under
## SSP5-8.5 using IPSL-CM6A-LR and ERA5 data. I am practicing understanding this
## type of data, doing a bias adjustment, downscaling the data and then analysing
## a random precipitation indicator and its threshold


# -------------------------------------------------------------------------

rm(list = ls())
library(terra)
library(ncdf4)
library(ggplot2)
library(dplyr)   
library(stars)
library(RNetCDF)


########################### Preparing and checking data ####################

#Upload data

hist_file <-"CMIP6_Historical_1981_2014/pr_day_IPSL-CM6A-LR_historical_r1i1p1f1_gr_19810101-20141231.nc"
era_file <- "ERA5 Total Precipitation/reanalysis-era5-land-timeseries-sfc-pressure-precipitation3rt9rcwn.nc"
ssp_file <- "CMIP6_SSP5-8.5/pr_day_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr_20150101-20301231.nc"

#Comparing data

nc_hist <- nc_open(hist_file)
nc_era <- nc_open(era_file)
nc_ssp <- nc_open(ssp_file)

#names(nc_era$var)
#names(nc_era$dim)
#names(nc_ssp$var)
#names(nc_ssp$dim)
#nc_close(nc_hist)
#nc_close(nc_era)

#nc_ssp$var$pr
#nc_hist$var$pr

## Processing historical data
###Put historical precipitation data  into R format
pr_hist <- ncvar_get(nc_hist, "pr")
#dim(pr_hist)
pr_hist[, , 1]

### Put data from km/m2/sec into mm/day 
pr_hist_mm <- pr_hist * 86400
pr_hist_mm[, , 1]

### Put spatial data in same format as ERA
lon_hist <- nc_hist$dim$lon$vals
lon_hist[lon_hist > 180] <- lon_hist[lon_hist > 180] - 360
#lon_hist


## Processing SSP data
###Put SSP data into R format and mm/day
pr_ssp <- ncvar_get(nc_ssp, "pr")
pr_ssp_mm <- pr_ssp * 86400
pr_ssp_mm[, , 1]
#dim(pr_ssp_mm)


## Processing ERA data
nc_era$var$tp
era_time <- nc_era$dim$valid_time$vals
era_dates <- as.POSIXct(
  era_time * 3600,
  origin = "1970-01-01",
  tz = "UTC"
)
#range(era_dates)
#head(diff(era_dates))
#head(ncvar_get(nc_era, "latitude"))
#head(ncvar_get(nc_era, "longitude"))

## Matching CMIP6 with ERA geographic point 
lat_hist <- nc_hist$dim$lat$vals
#which.min(abs(lat_hist - 40.6))
#lat_hist
#which.min(abs(lon_hist + 8.7))

hist_point <- pr_hist_mm[2, 2, ]
#length(hist_point)

# Checking SSP5 data
lat_ssp <- nc_ssp$dim$lat$vals
lon_ssp <- nc_ssp$dim$lon$vals
lon_ssp[lon_ssp > 180] <- lon_ssp[lon_ssp > 180] - 360

#lat_ssp
#lon_ssp


## Matching point with ERA geographical point
#which.min(abs(lat_ssp - 40.6))
#which.min(abs(lon_ssp + 8.7))

ssp_point <- pr_ssp_mm[2, 2, ]
#length(ssp_point)
#head(ssp_point)

##Convert ERA5 hourly to daily total precipitation
tp_era <- ncvar_get(nc_era, "tp")
tp_era_mm <- tp_era * 1000
era_df <- data.frame(
  datetime = era_dates,
  precip_mm = tp_era_mm
)
#head(era_df)
era_df$date <- as.Date(era_df$datetime)
era_daily <- era_df %>%
  group_by(date) %>%
  summarise(
    precip_mm = sum(precip_mm, na.rm = TRUE)
  )
#head(era_daily)
#tail(era_daily)

### Sense checking precipitation values on ERA5
era_2000 <- era_daily %>%
  filter(date >= as.Date("2000-01-01"),
         date <= as.Date("2000-12-31"))

sum(era_2000$precip_mm)
summary(era_2000$precip_mm)
max(era_2000$precip_mm)


##################### BIAS CORRECTION ##############################

## Subset ERA5 to 1981-2014 to match CMIP6 historical
era_hist <- era_daily %>%
  filter(
    date >= as.Date("1981-01-01"),
    date <= as.Date("2014-12-31")
  )
#head(era_hist)
#tail(era_hist)
#nrow(era_hist)
#nc_hist$dim$time

## Saving ERA5 time in yyyy-mm-dd format
hist_time <- nc_hist$dim$time$vals
hist_datetime <- as.POSIXct(
  hist_time * 86400,
  origin = "1850-01-01",
  tz = "UTC"
)
head(hist_datetime)
tail(hist_datetime)

hist_dates <- as.Date(hist_datetime)
#length(hist_dates)
all(hist_dates == era_hist$date)

## Putting ERA5 and CMIP6 Historical into same dataframe
comparison_hist <- data.frame(
  date = hist_dates,
  era_mm = era_hist$precip_mm,
  ipsl_mm = hist_point
)

#head(comparison_hist)
#summary(comparison_hist)

comparison_hist$era_mm[comparison_hist$era_mm < 0] <- 0

## Comparing annual precipitation data

annual_hist <- comparison_hist %>%
  mutate(year = format(date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    era_total = sum(era_mm, na.rm = TRUE),
    ipsl_total = sum(ipsl_mm, na.rm = TRUE)
  )

#head(annual_hist)

#mean(annual_hist$era_total)
#mean(annual_hist$ipsl_total)

## Calculate the bias
mean_bias <- mean(annual_hist$ipsl_total) - mean(annual_hist$era_total)
mean_bias

relative_bias <- (
  mean(annual_hist$ipsl_total) /
    mean(annual_hist$era_total) - 1
) * 100
relative_bias

## Plot the bias
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
    x = "Year",
    y = "Annual precipitation (mm)",
    color = "Dataset"
  ) +
  theme_minimal()

## Comparing monthly precipitation data
comparison_hist <- comparison_hist %>%
  mutate(
    year = as.numeric(format(date, "%Y")),
    month = as.numeric (format(date, "%m"))
  )
comparison_hist

monthly_hist <- comparison_hist %>%
  group_by(year, month) %>%
  summarise(
    era_total = sum(era_mm, na.rm = TRUE),
    ipsl_total = sum(ipsl_mm, na.rm = TRUE),
    .groups = "drop"
  )
#head(monthly_hist)

monthly_climatology <- monthly_hist %>%
  group_by(month) %>%
  summarise(
    era_mean = mean(era_total),
    ipsl_mean = mean(ipsl_total)
  )
monthly_climatology

monthly_plot <- monthly_climatology %>%
  tidyr::pivot_longer(
    cols = c(era_mean, ipsl_mean),
    names_to = "dataset",
    values_to = "precip_mm"
  )
head(monthly_plot)

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
    x = "Month",
    y = "Mean monthly precipitation (mm)",
    color = "Dataset",
    title = "Mean monthly precipitation (1981–2014)"
  ) +
  theme_minimal()

head(comparison_hist)
comparison_hist$era_wet <- comparison_hist$era_mm >= 1 #create a column where I check if era5 precipitation is >=1mm
comparison_hist$ipsl_wet <- comparison_hist$ipsl_mm >= 1 #create a column where I check if ipsl precipitation is >=1mm
comparison_hist
sum(comparison_hist$era_wet) #sums all values (TRUE=1 FALSE=0) so I can compare era5 and ipsl data
sum(comparison_hist$ipsl_wet)

tail(comparison_hist)



sum(comparison_hist$era_wet)/length(unique(comparison_hist$year)) # calculate average wet days in era5
sum(comparison_hist$ipsl_wet)/length(unique(comparison_hist$year)) #calculate average wet days in ipsl
  ### This has confirmed that the ipsl models a lot more wet days than era5, confirming a frequency bias on the precipitation data


mean(comparison_hist$era_mm[comparison_hist$era_wet == TRUE]) ## calculate average precipitation in wet days for era5
mean(comparison_hist$ipsl_mm[comparison_hist$ipsl_wet == TRUE]) ## calculate average precipitation in wet days for ipsl
 ### This test confirmed that there is an intensity bias as well

median(comparison_hist$era_mm[comparison_hist$era_wet == TRUE])
median(comparison_hist$ipsl_mm[comparison_hist$ipsl_wet == TRUE])
 ### The median is also bigger in ipsl, but it;s quite lower than the median, showing that there are a lot of low values and a few heavy precipitation events

quantile(comparison_hist$era_mm, probs = 0.95)
quantile(comparison_hist$ipsl_mm, probs = 0.95)

era_wet_days <- comparison_hist[comparison_hist$era_wet == TRUE,]
ipsl_wet_days <- comparison_hist[comparison_hist$ipsl_wet == TRUE, ]
max(era_wet_days$era_mm)
max(ipsl_wet_days$ipsl_mm)

ggplot(era_wet_days, aes( x = era_mm)) + geom_histogram(binwidth = 5) + scale_x_continuous(limits = c(0, 160), breaks = seq(0,160, by =10))
ggplot(ipsl_wet_days, aes(x = ipsl_mm)) + geom_histogram(binwidth = 5) + scale_x_continuous(limits = c(0, 160), breaks = seq(0,160, by =10))

calibration <- comparison_hist[comparison_hist$year < 2000,]
validation <- comparison_hist[comparison_hist$year >=2000,]
nrow(calibration)
nrow(validation)


# Bias-Correction using Quantile delta mapping ----------------------------

## Creating distinct objects so I can run the quantile delta mapping and confirm 
## it's working properly using validation data
era_cal <- calibration$era_mm
ipsl_cal <- calibration$ipsl_mm
ipsl_val <- validation$ipsl_mm

library(MBC)

qdm_val <- QDM(o.c = era_cal, m.c = ipsl_cal, m.p = ipsl_val, ratio = TRUE)

## Use the mhat.p data that has been created with the QDM, representing
## the adjusted validation data (IPSL data for 200-2014)

ipsl_val_corrected <- qdm_val$mhat.p

validation_qdm <- data.frame(date = validation$date, ERA = validation$era_mm, IPSL_raw = ipsl_val, IPSL_corrected = ipsl_val_corrected)
head(validation_qdm)

## Testing that the mean intensity of precipitation has been corrected 
## with the QDM, which we could see it was.
mean(validation_qdm$ERA)
mean(validation_qdm$IPSL_raw)
mean(validation_qdm$IPSL_corrected)

## Testing that the frequency of "wet days" (here classified as >=1mm)
## is closer to the observed data
validation_qdm$ERA_wet <- validation_qdm$ERA >=1
validation_qdm$IPSL_raw_wt <-  validation_qdm$IPSL_raw >=1
validation_qdm$IPSL_corrected_wet <- validation_qdm$IPSL_corrected >=1
head(validation_qdm)

sum(validation_qdm$ERA_wet)
sum(validation_qdm$IPSL_raw_wt)
sum(validation_qdm$IPSL_corrected_wet)

## Testing to see the quantiles of corrected data. With this analysis
## we noted that the p95 is actually quite different. This might be due
## due to a trend in climatic changes. 
quantile(validation_qdm$ERA[validation_qdm$ERA_wet], probs = 0.95)
quantile(validation_qdm$IPSL_raw[validation_qdm$IPSL_raw_wt], probs = 0.95)
quantile(validation_qdm$IPSL_corrected[validation_qdm$IPSL_corrected_wet], probs = 0.95)

## we can see that there was a big change in the quantiles for the 1981-1999
## period between ERA5 and IPSL data, meaning that the QDM worked, 
## but the difference in the IPSL 2000-2014 quantiles comes from a climatic trend.
## This shows that QDM can correct systematic statistical biases, but it cannot
## reproduce climatic trends that have not been simulated. 

calibration
quantile(calibration$era_mm[calibration$era_wet], probs = 0.95)
quantile(calibration$ipsl_mm[calibration$ipsl_wet], probs = 0.95)

quantile(qdm_val$mhat.c[qdm_val$mhat.c >= 1], probs = 0.95)

###################################################################
# We have just confirmed that QDM corrected the bias in the historical
# period, so now we do the QDM to the entire historical period
###################################################################

ipsl_hist <- comparison_hist$ipsl_mm
qdm_hist <- QDM(o.c = era_hist$precip_mm, m.c = ipsl_hist, m.p =  ssp_point, ratio = TRUE)
ssp5_corrected <- qdm_hist$mhat.p

ssp_time <- nc_ssp$dim$time$vals
ssp_dates <- as.Date("2015-01-01") + ssp_time

ssp5_qdm <- data.frame(Date = ssp_dates, IPSL_SSP5_raw = ssp_point, IPSL_SSP5_corrected = ssp5_corrected)

summary(ssp5_qdm$IPSL_SSP5_raw)
summary(ssp5_qdm$IPSL_SSP5_corrected)


# Threshold analysis ------------------------------------------
## For this analysis, a threshold of 5 consecutive days with 
## precipitation above 15mm/day was used

ssp5_qdm$Wet_days <- ssp5_qdm$IPSL_SSP5_corrected>15
wet_runs <- rle(ssp5_qdm$Wet_days)
precip_threshold <- wet_runs$values == TRUE & wet_runs$lengths >=5
sum(precip_threshold)
wet_runs$lengths[precip_threshold]

run_end <- cumsum(wet_runs$lengths)
event_end <- run_end[precip_threshold]
event_end
event_start <- event_end - wet_runs$lengths[precip_threshold] + 1
event_start

ssp5_qdm$Date[event_start]
ssp5_qdm$Date[event_end]

ssp5_qdm$Threshold_event <- FALSE
ssp5_qdm

for (i in 1:length(event_start)) {
  
  ssp5_qdm$Threshold_event[event_start[i]:event_end[i]] <- TRUE
  
}

### Creating a graph with consecutive days of heavy precipitation

ssp5_qdm$Wet_days <- ssp5_qdm$IPSL_SSP5_corrected >15
ssp5_qdm$Consecutive_wet_days <- 0

if (ssp5_qdm$Wet_days[1] == TRUE){
  ssp5_qdm$Consecutive_wet_days[1] = 1
}

for (i in 2:nrow(ssp5_qdm)) {
    if (ssp5_qdm$Wet_days[i] == TRUE) {
      ssp5_qdm$Consecutive_wet_days[i] = 1 + ssp5_qdm$Consecutive_wet_days[i-1]
    }
      else {ssp5_qdm$Consecutive_wet_days[i] = 0} 
}



threshold_events <- ssp5_qdm[ssp5_qdm$Consecutive_wet_days==5, ]
threshold_events

library(ggrepel)

ggplot(data= ssp5_qdm, aes(x= Date, y = Consecutive_wet_days)) +
  scale_y_continuous(breaks = 0:5)+
  geom_point() +
  geom_hline( aes(yintercept = 5, linetype  = "Hypothetical threshold"),colour = "red", linewidth = 0.8) + 
  geom_point(data= threshold_events) + 
  geom_text_repel(data = threshold_events, aes(label = Date), vjust = +1) +
  labs (
    title = "Projected consecutive high precipitation (>15mm) events",
    subtitle = "IPSL-CM6A-LR, SSP5-8.5 (2015–2030)",
    x = "Year",
    y = "Number of days"
  )
