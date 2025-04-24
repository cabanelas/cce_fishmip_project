################################################################################
#############                FEISTY                #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(ncdf4) #to open netcdf
#library(raster) #for brick()
#library(ocedata) #for plot()
library(tidyverse)
#library(easyNCDF)
library(lubridate)
library(terra) # for bilinear interpolation
library(gganimate)
library(ggOceanMaps)
library(rnaturalearth)
library(gifski)
library(sf)
library(stars)

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
years <- floor(1601 + time / 12) 

annual_list_native <- list()

for (yr in unique(years)) { #table(floor(years))
  message("Processing year: ", yr)
  
  idx <- which(years == yr)
  
  monthly_stack <- lapply(idx, function(i) {
    slice <- ncvar_get(feisty, "tcb", start = c(1, 1, i), count = c(-1, -1, 1))
    slice[slice > 1e19] <- NA  
    rast(slice)
  })
  
  r_stack <- rast(monthly_stack)
  ext(r_stack) <- c(min(lon), max(lon), min(lat), max(lat))
  crs(r_stack) <- "EPSG:4326"
  
  # flip lat
  r_stack <- flip(r_stack, direction = "vertical")
  
  # mean across 12 months
  r_mean <- app(r_stack, mean, na.rm = TRUE)

  annual_list_native[[as.character(yr)]] <- r_mean
}


# Get variable names
var_name <- names(ecocean$var)[1]
# Read in the variables
data_array <- ncvar_get(ecocean, var_name)
# Read time units and convert to dates
time_units <- ncatt_get(ecocean, "time", "units")$value
time_origin <- sub(".*since ", "", time_units)
time_dates <- as.Date(time, origin = time_origin)

## ------------------------------------------ ##
#    2) Spatial aggregation -----
## ------------------------------------------ ##
feisty_stack <- rast(annual_list_native) 

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