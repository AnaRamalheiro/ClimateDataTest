rm(list = ls())
setwd("C:\\Users\\anafi\\OneDrive - Wageningen University & Research\\Desktop\\R\\Climate data_test")
library(terra)
library(ncdf4)

cmip_hist <- nc_open ("CMIP6_Historical_1981_2014\\pr_day_IPSL-CM6A-LR_historical_r1i1p1f1_gr_19810101-20141231.nc")
cmip_SSP5 <- nc_open("CMIP6_SSP5-8.5\\pr_day_IPSL-CM6A-LR_ssp585_r1i1p1f1_gr_20150101-20301231.nc")
era5_hist <- nc_open("ERA5 Total Precipitation\\reanalysis-era5-land-timeseries-sfc-pressure-precipitation3rt9rcwn.nc")


length(cmip_hist)
names(cmip_hist)
mode(cmip_hist)
cmip_hist
head(cmip_hist)
tail(cmip_hist)
cmip_hist$var$pr
mode(cmip_hist$var$pr)
cmip_hist$var
