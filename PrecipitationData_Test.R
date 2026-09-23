rm(list = ls())
library(terra)
library(ncdf4)
library(ggplot2)
library(dplyr)   
library(stars)
library(RNetCDF)


########################### Preparing and checking data ####################

#Upload data

hist_file <- "C:/Users/anafi/OneDrive - Wageningen University & Research/Desktop/R/Climate data_test/CMIP6_Historical_1981_2014/pr_day_IPSL-CM6A-LR_historical_r1i1p1f1_gr_19810101-20141231.nc"
era_file <- "C:/Users/anafi/OneDrive - Wageningen University & Research/Desktop/R/Climate data_test/ERA5 Total Precipitation/reanalysis-era5-land-timeseries-sfc-pressure-precipitation3rt9rcwn.nc"
ssp_file <- "C:/Users/anafi/OneDrive - Wageningen University & Research/Desktop/R/Climate data_test/CMIP6_SSP5-8.5/pr_day_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr_20150101-20301231.nc"

#Comparing data

nc_hist <- nc_open(hist_file)
nc_era <- nc_open(era_file)
nc_ssp <- nc_open(ssp_file)

names(nc_era$var)
names(nc_era$dim)
names(nc_ssp$var)
names(nc_ssp$dim)
#nc_close(nc_hist)
#nc_close(nc_era)

nc_ssp$var$pr
nc_hist$var$pr

## Processing historical data
###Put historical precipitation data  into R format
pr_hist <- ncvar_get(nc_hist, "pr")
dim(pr_hist)
pr_hist[, , 1]

### Put data from km/m2/sec into mm/day 
pr_hist_mm <- pr_hist * 86400
pr_hist_mm[, , 1]

### Put spatial data in same format as ERA
lon_hist <- nc_hist$dim$lon$vals
lon_hist[lon_hist > 180] <- lon_hist[lon_hist > 180] - 360
lon_hist


## Processing SSP data
###Put SSP data into R formar and mm/day
pr_ssp <- ncvar_get(nc_ssp, "pr")
pr_ssp_mm <- pr_ssp * 86400
pr_ssp_mm[, , 1]
dim(pr_ssp_mm)


## Processing ERA data
nc_era$var$tp
era_time <- nc_era$dim$valid_time$vals
era_dates <- as.POSIXct(
  era_time * 3600,
  origin = "1970-01-01",
  tz = "UTC"
)
range(era_dates)
head(diff(era_dates))
head(ncvar_get(nc_era, "latitude"))
head(ncvar_get(nc_era, "longitude"))

## Matching CMIP6 with ERA geographic point 
lat_hist <- nc_hist$dim$lat$vals
which.min(abs(lat_hist - 40.6))
lat_hist
which.min(abs(lon_hist + 8.7))

hist_point <- pr_hist_mm[2, 2, ]
length(hist_point)

# Checking SSP5 data
lat_ssp <- nc_ssp$dim$lat$vals
lon_ssp <- nc_ssp$dim$lon$vals
lon_ssp[lon_ssp > 180] <- lon_ssp[lon_ssp > 180] - 360

lat_ssp
lon_ssp


## Matching point with ERA geographical point
which.min(abs(lat_ssp - 40.6))
which.min(abs(lon_ssp + 8.7))

ssp_point <- pr_ssp_mm[2, 2, ]
length(ssp_point)
head(ssp_point)

##Convert ERA5 hourly to daily total precipitation
tp_era <- ncvar_get(nc_era, "tp")
tp_era_mm <- tp_era * 1000
era_df <- data.frame(
  datetime = era_dates,
  precip_mm = tp_era_mm
)
head(era_df)
era_df$date <- as.Date(era_df$datetime)
era_daily <- era_df %>%
  group_by(date) %>%
  summarise(
    precip_mm = sum(precip_mm, na.rm = TRUE)
  )
head(era_daily)
tail(era_daily)

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
head(era_hist)
tail(era_hist)
nrow(era_hist)
nc_hist$dim$time

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
length(hist_dates)
all(hist_dates == era_hist$date)

## Putting ERA5 and CMIP6 Historical into same dataframe
comparison_hist <- data.frame(
  date = hist_dates,
  era_mm = era_hist$precip_mm,
  ipsl_mm = hist_point
)

head(comparison_hist)
summary(comparison_hist)

comparison_hist$era_mm[comparison_hist$era_mm < 0] <- 0

## Comparing annual precipitation data

annual_hist <- comparison_hist %>%
  mutate(year = format(date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    era_total = sum(era_mm, na.rm = TRUE),
    ipsl_total = sum(ipsl_mm, na.rm = TRUE)
  )

head(annual_hist)

mean(annual_hist$era_total)
mean(annual_hist$ipsl_total)

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
head(monthly_hist)

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
 ###