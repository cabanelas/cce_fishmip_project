# Climate Change Ecology Final Project
  # GH repo: https://github.com/cabanelas/cce_fishmip_project

library(ncdf4)
library(lubridate)
library(sf)
library(stars)
library(dplyr)
library(exactextractr)
library(PCICt)
library(raster)
library(exactextractr) #v0.10.0
library(gganimate) #v1.0.8
library(gifski) #v1.32.0.
library(SpatRaster)

# open NetDC file
ecocean <- nc_open("/Users/kjehickman/Documents/Classes/Spring_25/7-430_CCData/ccecology_final_project/ecoocean_gfdl-esm4_nobasd_historical_nat_default_tcb_global_monthly_1950_2014.nc")
zooms <- nc_open("/Users/kjehickman/Documents/Classes/Spring_25/7-430_CCData/ccecology_final_project/zoomss_gfdl-esm4_nobasd_picontrol_nat_default_tcb_global_annual_2015_2100.nc")

# open shape file
lme <- st_read("/Users/kjehickman/Documents/Classes/Spring_25/7-430_CCData/ccecology_final_project/lme66/LME66.shp")

print(ecocean)

# Get variable names
var_name <- names(ecocean$var)[1]
ecocean$var[[var_name]]$dimids 
sapply(ecocean$var[[var_name]]$dim, function(x) x$name)  # lon, lat, time

# Read in the variables
data_array <- ncvar_get(ecocean, var_name)
dim(data_array)  # lon:360, lat:180, time:780

# Get dimensions
lon <- ncvar_get(ecocean, "lon")
lat <- ncvar_get(ecocean, "lat")
time <- ncvar_get(ecocean, "time")

# Aggregate by time ####

# Read time units and convert to dates--ecocean
time_units <- ncatt_get(ecocean, "time", "units")$value
time_origin <- sub(".*since ", "", time_units)
time_dates <- as.Date(time, origin = time_origin)

# Extract years
years <- year(time_dates)
unique_years <- sort(unique(years))
n_years <- length(unique_years)

# Create an array to hold the annual averages
  # If your array is [lon, lat, time]
annual_array <- array(NA, dim = c(length(lon), length(lat), n_years))
# annual_array <- array(NA, dim = c(dim(data_array)[1], dim(data_array)[2], n_years))

for (i in seq_along(unique_years)) {
  year_idx <- which(years == unique_years[i])
  annual_array[,,i] <- apply(data_array[,,year_idx], c(1, 2), mean, na.rm = TRUE)
}

# Define dimensions for output
lon_dim <- ncdim_def("lon", "degrees_east", lon)
lat_dim <- ncdim_def("lat", "degrees_north", lat)
time_dim <- ncdim_def("time", paste("years since", min(unique_years)), unique_years)

# Define variable
var_def <- ncvar_def(var_name,
                     units = ecocean$var[[var_name]]$units,
                     dim = list(lon_dim, lat_dim, time_dim),
                     missval = NA, prec = "float")

# Create and write to a new file
ecocean_out <- nc_create("ecocean_annual.nc", list(var_def))
ncvar_put(ecocean_out, var_def, annual_array)


# Aggregate Spatially ####

# aggregate to 5-degree resolution
  # make sure aggregation matches ZooMS layer

# 1. Load the ecocean file as a RasterBrick (multi-band = time)
ecocean_r <- brick("ecocean_annual.nc")

# 2. Load the zooms NetCDF and use it to define the new resolution grid
zooms_r <- brick("/Users/kjehickman/Documents/Classes/Spring_25/7-430_CCData/ccecology_final_project/zoomss_gfdl-esm4_nobasd_historical_nat_default_tcb_global_annual_1950_2014.nc")  # assuming it has lon/lat structure

  # Check projection and match if needed
crs(ecocean_r) <- crs(zooms_r)

# 3. Resample using bilinear interpolation
ecocean_resampled <- resample(ecocean_r, zooms_r, method = "bilinear")

# 4. Save the resampled object as a new NetCDF
writeRaster(ecocean_resampled, filename = "ecocean_annual_5deg.nc", format = "CDF", overwrite = TRUE)


# Clip netCDF to shp file ####

# selet LMEs to target:
  ## 10 = HOTS
  ## 14 = Patagonian Shelf
  ## 21 = Norwegian Shelf
  ## 28 = Guinea Current

# read in spatially & temporally aggregated CDF as Raster
ecocean_full <- brick("ecocean_annual_5deg.nc")

names(lme)
head(lme)

# subset to 4 LME of interest
lme_subset <- lme %>% 
  filter(LME_NUMBER==10 | LME_NUMBER==14| LME_NUMBER==21 | LME_NUMBER==28)

extract_lme_timeseries <- function(rast, shape, start_year = 1950) {
  n_years <- nlayers(rast)
  years <- start_year:(start_year + n_years - 1)
  # 
  out <- lapply(1:n_years, function(i) {
    layer <- rast[[i]]
    # Extract values and coverage fractions as a dataframe
    extract_df <- exact_extract(layer, shape)[[1]]
    # Remove NA values
    extract_df <- extract_df[!is.na(extract_df$value), ]
    if (nrow(extract_df) == 0) {
      return(data.frame(Year = years[i], TCB_mean = NA, TCB_median = NA, TCB_sd = NA))
    }
    
    values <- extract_df$value
    weights <- extract_df$coverage_fraction
    
    # Weighted mean
    w_mean <- sum(values * weights) / sum(weights)
    
    # Weighted median (approximate method)
    ord <- order(values)
    values_ord <- values[ord]
    weights_ord <- weights[ord]
    cum_weights <- cumsum(weights_ord) / sum(weights_ord)
    w_median <- values_ord[which.min(abs(cum_weights - 0.5))]
    
    # Weighted standard deviation
    w_var <- sum(weights * (values - w_mean)^2) / sum(weights)
    w_sd <- sqrt(w_var)
    
    data.frame(Year = years[i], TCB_mean = w_mean, TCB_median = w_median, TCB_sd = w_sd)
  })
  
  do.call(rbind, out)
}

# call for each LME
ts_10 <- extract_lme_timeseries(ecocean_full, lme_subset %>% filter(LME_NUMBER == 10)) %>% mutate(LME = "10")
ts_14 <- extract_lme_timeseries(ecocean_full, lme_subset %>% filter(LME_NUMBER == 14)) %>% mutate(LME = "14")
ts_21 <- extract_lme_timeseries(ecocean_full, lme_subset %>% filter(LME_NUMBER == 21)) %>% mutate(LME = "21")
ts_28 <- extract_lme_timeseries(ecocean_full, lme_subset %>% filter(LME_NUMBER == 28)) %>% mutate(LME = "28")

# binda all LME output into df
tcb_by_lme <- bind_rows(ts_10, ts_14, ts_21, ts_28)

# # Plot masked result
summary(tcb_by_lme$TCB_mean)
tcb_by_lme <- tcb_by_lme %>% 
  mutate(LME_name = case_when(LME == 10 ~ "HOTS",
                   LME == 14 ~ "Patagonian Shelf",
                   LME == 21 ~ "Norwegian Shelf",
                   TRUE ~ "Guinean Current"))

# export to csv

ggplot(tcb_by_lme, aes(x = Year, y = TCB_mean)) +
  geom_ribbon(aes(ymin = TCB_mean - TCB_sd, ymax = TCB_mean + TCB_sd, fill = LME), alpha = 0.2) +
  geom_line(aes(color = LME), size = 1) +
  facet_wrap(~LME_name, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "EcoOcean Annual Total Community Biomass (Mean ± 1 SD) by LME",
    x = "Year", y = "TCB"
  ) +
  ylim(-100, 700)+
  theme(legend.position = "none")

# graph map of 2014 - 1950 to show difference over time
  # use the CDF for this


# Clean up
nc_close(ecocean)
nc_close(ecocean_out)

# make a .csv containing: year, tcb mean, tcb median, tcb stdev, tcb max, tcb min, LME number
  # stdev = (data_point - mean)/weighted_grid_cell (as n)
  # coverage fraction comes from exact_extract, this is spatial weight
# figures:
  # Olivia: overlay LME shape file onto full map, explain why these LMEs (basemap)
  # All: each person does time-series for each LME (4-panels)
  # Jo: then show 4 panels (1 for each LME) of line graph tcb over time
    # have comparison of 3 models in one line graph for each LME
  # Olivia: 3-panel showing difference in model predictions
  # Alex: gif of most interesting LME, 3 panels 1 for each 
    # LME of interest is either 14 or 21 (shelves)
  # ?: for case study, show static time difference for case study LME
  # subtract 2014 from 1950 and show on a map for each model



