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

#library(gganimate) #v1.0.8
#library(gifski) #v1.32.0.2

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
feisty <- nc_open("raw/feisty_gfdl-esm4_nobasd_historical_nat_default_tcb_global_monthly_1950_2014.nc")

## ------------------------------------------ ##
#     Explore file -----
## ------------------------------------------ ##
print(feisty)
names(feisty$var)  
#feisty$var$tcb
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

feisty_stack <- rast(annual_list_native) # annual feisty data
#writeRaster(feisty_stack, "feisty_stack_output.tif", overwrite = TRUE)

## ------------------------------------------ ##
#    2) LMEs -----
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

####################
# 2a) All 66 LMEs - loop -----
####################

# empty list
lme_results <- vector("list", length = nrow(lme))

library(Hmisc)  # for wtd.quantile

for (i in seq_len(nrow(lme))) {
  message("Processing LME: ", lme$LME_NAME[i])
  
  lme_poly <- lme[i, ]  # single LME
  
  extracted_data <- lapply(1:nlyr(feisty_stack), function(j) {
    exact_extract(feisty_stack[[j]], lme_poly)
  })
  
  # weighted stats
  processed_stats <- lapply(extracted_data, function(df_list) {
    df <- df_list[[1]] # extract actual data.frame
    
    df <- df[!is.na(df$value) & !is.nan(df$value), ]  # clean
    
    if (nrow(df) > 0 && all(c("value", "coverage_fraction") %in% names(df))) {
      mean_val <- weighted.mean(df$value, df$coverage_fraction, na.rm = TRUE)
      median_val <- Hmisc::wtd.quantile(df$value, df$coverage_fraction, probs = 0.5, na.rm = TRUE)
      var_val <- sum(df$coverage_fraction * (df$value - mean_val)^2) / sum(df$coverage_fraction)
      sd_val <- sqrt(var_val)
      
      list(mean = mean_val, median = median_val, sd = sd_val)
    } else {
      list(mean = NA, median = NA, sd = NA)
    }
  })
  
  lme_results[[i]] <- tibble(
    LME_NUMBER = lme$LME_NUMBER[i],
    LME_NAME = lme$LME_NAME[i],
    year = 1950:2014,
    tcb_mean = sapply(processed_stats, `[[`, "mean"),
    tcb_median = sapply(processed_stats, `[[`, "median"),
    tcb_sd = sapply(processed_stats, `[[`, "sd")
  )
}

lme_all_df <- bind_rows(lme_results)
head(lme_all_df)

#write.csv(lme_all_df, "output/FEISTY_tcb_66lme.csv")

## ------------------------------------------ ##
#    3) Plot -----
## ------------------------------------------ ##
tcb_by_lme <- lme_all_df %>%
  filter(LME_NUMBER %in% c(10, 14, 21, 28))

colnames(tcb_by_lme) <- c("LME", "LME_name", "Year", "TCB_mean", 
                          "TCB_median", "TCB_sd")

(p <- ggplot(tcb_by_lme, aes(x = Year, y = TCB_mean)) +
  geom_ribbon(aes(ymin = TCB_mean - TCB_sd, ymax = TCB_mean + TCB_sd, fill = LME), 
              alpha = 0.2) +
  geom_line(aes(color = LME), size = 1) +
  facet_wrap(~ LME_name, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Annual Total Community Biomass (Mean ± 1 SD)\nFEISTY",
    x = "Year", y = expression("TCB (g/m"^2*")")
  ) +
  #ylim(-100, 50)+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    strip.text = element_text(size = 13, face = "bold"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.text.x = element_text(size = 12, color = "black"),  
    axis.text.y = element_text(size = 12, color = "black")   
  )
)

#ggsave("plot_FEISTY_4LMEtcb_line.png", plot = p, width = 10, height = 8, dpi = 300,
#       bg = "white")

#write.csv(tcb_by_lme, "output/FEISTY_tcb_by_lme.csv")

## ------------------------------------------ ##
#   THE END 
## ------------------------------------------ ##
tcb_nes <- lme_all_df %>%
  filter(LME_NAME %in% c("Northeast U.S. Continental Shelf"))

colnames(tcb_nes) <- c("LME", "LME_name", "Year", "TCB_mean", 
                          "TCB_median", "TCB_sd")
#NES
(nesp <- ggplot(tcb_nes, aes(x = Year, y = TCB_mean)) +
  geom_ribbon(aes(ymin = TCB_mean - TCB_sd, ymax = TCB_mean + TCB_sd, fill = "skyblue"), 
              alpha = 0.3) +
  #geom_ribbon(aes(ymin = q05, ymax = q95), fill = "skyblue", alpha = 0.3) +
  geom_line(color = "black", linewidth = 1) +
  labs(title = "Mean TCB (NES LME, FEISTY ISIMIP)",
       y = expression("TCB (g/m"^2*")"),
       x = NULL) +
  theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
      strip.text = element_text(size = 13, face = "bold"),  
      axis.title.x = element_text(size = 14, face = "bold"),  
      axis.title.y = element_text(size = 15, face = "bold"), 
      axis.text.x = element_text(size = 12, color = "black"),  
      axis.text.y = element_text(size = 12, color = "black")   
    ) +
  ylim(0, 30)
)
#ggsave("FEISTY_NES_isimipdata.png", plot = nesp, width = 10, height = 6, dpi = 300, bg = "white")

# caterpillar plot 
lme1 <- st_make_valid(lme)

# Extract the latitude of the centroid of each LME (central latitude)
lme1 <- lme1 %>%
  mutate(lat = st_coordinates(st_centroid(geometry))[, 2])

#lme_sorted <- lme1 %>% arrange(lat)
lme_sorted <- lme1 %>% 
  arrange(desc(lat)) %>%
  mutate(position = row_number()) %>%
  select(LME_NUMBER, LME_NAME, position, lat)

lme_all_df1 <- left_join(lme_all_df, lme_sorted, by = c("LME_NUMBER", "LME_NAME"))

(catplot <- ggplot(lme_all_df1, aes(x = tcb_mean , y = reorder(LME_NAME, -position))) +
  geom_point() +  # Plot the mean
  geom_errorbarh(aes(xmin = tcb_mean - tcb_sd,
                     xmax = tcb_mean + tcb_sd), height = 0) +  
  labs(
    title = "FEISTY TCB (Mean ± 1 SD) 1950-2014",
    x = expression("TCB (g/m"^2*")"),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) 
)
#ggsave("plot_FEISTY_caterpillar.png", plot = catplot, width = 8, height = 10, dpi = 300, bg = "white")

## ------------------------------------------ ##
#   DID NOT DO -- CODE BELOW IF REGRID NEEDED 
## ------------------------------------------ ##
## ------------------------------------------ ##
#    2) ReGrid : bilinear interpolation ----- not needed, skip to 3
## ------------------------------------------ ##
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


###ORIGINAL loop


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