################################################################################
#############                FEISTY                #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################
# FEISTY data downloaded from ISIMIP repository on 24-APR-25
# Looking at total consumer biomass = tcb

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(ncdf4) #to open netcdf; v1.21
library(tidyverse) #v2.0.0
library(lubridate) #v1.9.3
library(terra) # for bilinear interpolation; #v1.7.29
library(sf) #v1.0.15
library(stars) #v0.6.2
library(exactextractr) #v0.10.0

library(gganimate) #v1.0.8
library(gifski) #v1.32.0.2

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
dim(ncvar_get(feisty, "tcb"))

## ------------------------------------------ ##
#     Extract variables and coordinates -----
## ------------------------------------------ ##
lon <- ncvar_get(feisty, "lon")     # 360
lat <- ncvar_get(feisty, "lat")     # 180
time <- ncvar_get(feisty, "time")   # 780

ncatt_get(feisty, "time", "units")$value # "months since 1601-01-01"
years <- floor(1601 + time / 12) 
table(floor(years))

## ------------------------------------------ ##
#    1) Temporal aggregation -----
## ------------------------------------------ ##
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
#    2) ReGrid : bilinear interpolation ----- not needed, skip to 3
## ------------------------------------------ ##
feisty_stack <- rast(annual_list_native) # annual feisty data

zooms <- rast("raw/zoomss_gfdl-esm4_nobasd_historical_nat_default_tcb_global_annual_1950_2014.nc")
names(zooms)
crs(zooms)
res(zooms)

zoomss <- nc_open("raw/zoomss_gfdl-esm4_nobasd_historical_nat_default_tcb_global_annual_1950_2014.nc")
lonzoo <- ncvar_get(zoomss, "lon")
latzoo <- ncvar_get(zoomss, "lat")
timezoo <- ncvar_get(zoomss, "time") 

# check grid
range(lonzoo) # 1 degree spacing ;  -179.5  179.5
length(lonzoo)
range(latzoo) # 1 degree spacing ;-89.5  89.5
length(latzoo)

zoo_slice <- ncvar_get(zoomss, "tcb", start = c(1, 1, 1), count = c(-1, -1, 1))
image.plot(lonzoo, latzoo, zoo_slice)
## ZOOMSS IS ONE DEGREEEEEEEE

feisty_stars <- st_as_stars(feisty_stack) 
zooms_stars <- st_as_stars(zooms)

# warp FEISTY to match ZooMSS using bilinear interpolation
feisty_5deg <- st_warp(feisty_stars, #annual time; another opt resample() fnction
                       dest = zooms_stars, 
                       method = "bilinear", 
                       use_gdal = TRUE)

plot(feisty_5deg[,,,1], main = "FEISTY Year 1950 (Warped to ZooMSS Grid)")
#still one degree

## ------------------------------------------ ##
#    3) LMEs -----
## ------------------------------------------ ##
lme <- st_read("raw/lme66/LME66.shp")

names(lme)
head(lme)

# we are using LME #s
## 10 = Hawaii
## 14 = Patagonian Shelf
## 21 = Norwegian Shelf
## 28 = Guinea Current 

lme_patagonia <- lme %>% 
  filter(LME_NUMBER %in% c(14))

patagonia_mean <- lapply(1:nlyr(feisty_stack), function(i) {
  exact_extract(feisty_stack[[i]], lme_patagonia, 
                fun = 'weighted_mean', weights = "area")
})
annual_weighted_means <- lapply(1:dim(src_resampled)[3], function(i) {
  exact_extract(src_resampled[,,,i], shape, 'weighted_mean')
})

####################
# 3a) All 66 LMEs - loop -----
####################
# empty list
lme_results <- vector("list", length = nrow(lme))

for (i in seq_len(nrow(lme))) {
  message("Processing LME: ", lme$LME_NAME[i])
  
  lme_poly <- lme[i, ]  # single LME as sf
  
  means <- lapply(1:nlyr(feisty_stack), function(j) {
    exact_extract(feisty_stack[[j]], lme_poly, 
                  fun = "weighted_mean", weights = "area")
  })
  
  lme_results[[i]] <- tibble(
    LME_NUMBER = lme$LME_NUMBER[i],
    LME_NAME = lme$LME_NAME[i],
    year = 1950:2014,
    tcb = unlist(means)
  )
}

lme_all_df <- bind_rows(lme_results)


