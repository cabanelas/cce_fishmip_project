################################################################################
#############                FEISTY                #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################
# FEISTY data downloaded from ISIMIP repository on 24-APR-25

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(ncdf4) #to open netcdf
library(tidyverse)
library(lubridate)
library(terra) # for bilinear interpolation
library(sf)
library(stars)
library(gganimate)
library(gifski)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
feisty <- nc_open("raw/feisty_gfdl-esm4_nobasd_historical_nat_default_tcb_global_monthly_1950_2014.nc")

## ------------------------------------------ ##
#     Explore file -----
## ------------------------------------------ ##
print(feisty)
names(feisty$var)  
feisty$var$tcb

## ------------------------------------------ ##
#    1) Temporal aggregation -----
## ------------------------------------------ ##
lon <- ncvar_get(feisty, "lon")     # 360
lat <- ncvar_get(feisty, "lat")     # 180
time <- ncvar_get(feisty, "time")   # 780

ncatt_get(feisty, "time", "units")$value # "months since 1601-01-01"
years <- floor(1601 + time / 12) 
table(floor(years))

annual_list_native <- list()

for (yr in unique(years)) { 
  message("Processing year: ", yr)
  
  idx <- which(years == yr)
  
  monthly_stack <- lapply(idx, function(i) {
    slice <- ncvar_get(feisty, "tcb", start = c(1, 1, i), count = c(-1, -1, 1))
    slice[slice > 1e19] <- NA  
    slice_t <- t(slice) #slice[, rev(seq_along(lat))]
      ## t(slice) needed bc NetCDF == [longitude, latitude]
      ##                    terra:rast == [latiude, longitude]
    rast(slice_t) # convert 2D matrix to SpatRaster 
  })
  
  # combine the list of monthly rasters into a single SpatRaster stack
  r_stack <- rast(monthly_stack) 
  
  ext(r_stack) <- ext(min(lon), max(lon), min(lat), max(lat))
  crs(r_stack) <- "EPSG:4326"
  
  # mean across 12 months
  r_mean <- app(r_stack, mean, na.rm = TRUE)

  annual_list_native[[as.character(yr)]] <- r_mean
}
#sanity check
plot(annual_list_native[["2000"]], main = "FEISTY Biomass — Year 2000")


## ------------------------------------------ ##
#    2) Spatial aggregation -----
## ------------------------------------------ ##
feisty_stack <- rast(annual_list_native) # annual feisty data

zooms <- rast("raw/zoomss_gfdl-esm4_nobasd_historical_nat_default_tcb_global_annual_1950_2014.nc")
names(zooms)
crs(zooms)

feisty_stars <- st_as_stars(feisty_stack) 
zooms_stars <- st_as_stars(zooms)

feisty_resampled <- st_warp(feisty_stars, #annual time
                            dest = zooms_stars, 
                            method = "bilinear", 
                            use_gdal = TRUE)


#spatial aggregate
src_resampled <- st_warp(src_nc_annual, dest = target_nc, method = "bilinear")
ecocean_5deg <- st_warp(src = ecocean, dest = zooms, method = "average")
nnual_weighted_means <- lapply(1:dim(src_resampled)[3], function(i) {
  exact_extract(src_resampled[,,,i], shape, 'weighted_mean')
})


#clip