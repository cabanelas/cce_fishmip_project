################################################################################
#############             ZooMSS  GIF MAPS TS      #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################
#ZooMSS

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(ncdf4) #to open netcdf
library(tidyverse)
library(gifski)
library(sf)
library(gganimate)
library(terra) #rast()
library(ggspatial)
library(marmap) # for bathymetry
library(scico)
library(rnaturalearth)
library(metR) # for geom text contour

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
zoo <- rast("netCDF/ZooMSS_annual_1deg.nc")

nlyr(zoo)

names(zoo)
names(zoo) <- 1950:2014

# LME
lme <- st_read("raw/lme66/LME66.shp") 

subset_lmes <- lme %>% 
  filter(LME_NUMBER %in% c(10, 14, 21, 28))

subset_lmes %>% 
  group_by(LME_NUMBER) %>% 
  group_split() %>% 
  purrr::map(~ list(
    lme_number = unique(.x$LME_NUMBER),
    lme_name = unique(.x$LME_NAME), 
    bbox = st_bbox(.x)
  ))

land <- ne_countries(scale = "medium", returnclass = "sf")

## ----------- ##
#   LME10  ---
## ----------- ##
## 10 = Hawaii
bbox_lme10 <- c(xmin = -179.9, ymin = 18, #-179.77818, ymin = 17.39859
                xmax = -152, ymax = 30) #xmax = -153.47676, ymax = 29.14307
bbox_lme10_ext <- ext(bbox_lme10[c("xmin", "xmax", "ymin", "ymax")])

# download and process the bathymetry data
bathyHawaii <- getNOAA.bathy(
  lon1 = -180, lon2 = -152, 
  lat1 = 16, lat2 = 30, 
  resolution = 4
)
bathy_Hawaii <- fortify.bathy(bathyHawaii)

## ----------- ##
#   LME14  ---
## ----------- ##
## 14 = Patagonian Shelf
bbox_lme14 <- c(xmin = -69.61389, xmax = -51.74283, 
                ymin = -55.04733, ymax = -32.44833)
bbox_lme14_ext <- ext(bbox_lme14[c("xmin", "xmax", "ymin", "ymax")])

# bathymetry data
bathyPatagon <- getNOAA.bathy(
  lon1 = -71, lon2 = -50, 
  lat1 = -57, lat2 = -30, 
  resolution = 4
)
bathy_Patagon <- fortify.bathy(bathyPatagon)

## ----------- ##
#   LME21  ---
## ----------- ##
## 21 = Norwegian Shelf
bbox_lme21 <- c(xmin = -7.952473, ymin = 62.005292,
                xmax = 17.046112, ymax = 78) #77.004961
bbox_lme21_ext <- ext(bbox_lme21[c("xmin", "xmax", "ymin", "ymax")])

# bathymetry data
bathyNorw <- getNOAA.bathy(
  lon1 = -9, lon2 = 19, 
  lat1 = 60, lat2 = 79, 
  resolution = 4
)
bathy_Norw <- fortify.bathy(bathyNorw)

## ----------- ##
#   LME28  ---
## ----------- ##
## 28 = Guinea Current 
bbox_lme28 <- c(xmin = -20, ymin = -8,
                xmax = 13, ymax = 12) 
bbox_lme28_ext <- ext(bbox_lme28[c("xmin", "xmax", "ymin", "ymax")])

# bathymetry data
bathyGuinea <- getNOAA.bathy(
  lon1 = -21, lon2 = 14, 
  lat1 = -9, lat2 = 13, 
  resolution = 4
)
bathy_Guinea <- fortify.bathy(bathyGuinea)

## ------------------------------------------ ##
#            ZooMSS map -----
## ------------------------------------------ ##
plot(zoo$"1950")

## ----------- ##
#   LME10  --- Hawaii
## ----------- ##
zoo_lme10 <- crop(zoo, bbox_lme10_ext)
zoo_df_lme10 <- as.data.frame(zoo_lme10, xy = TRUE)

# pivot long
zoo_l_lme10 <- zoo_df_lme10 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

zoo_10_p <- ggplot() +
  geom_raster(data = zoo_l_lme10, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"), 
                   na.value = "transparent") +
  # VIRIDIS color Version
  #scale_fill_viridis_c(
  #  option = "viridis",  
  #  direction = 1,       
  #  begin = 0,
  #  end = 1,
  #  name = expression("TCB Biomass (g/m"^2*")"),
  #  na.value = "transparent"
  #)
  # DIFFERENT LANDMASS AES
  #geom_sf(data = land, fill = "gray70", color = "black", size = 0.3) +
  borders("world", colour = "gray85", fill = "#575151") +
  geom_contour(
    data = bathy_Hawaii,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -4000),
    color = "gray50",
    size = 0.3,
    alpha = 0.3
  ) +
  geom_text_contour(
    data = bathy_Hawaii,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -5000),
    size = 2.5,           
    color = "gray60",    
    skip = 1,
    min.size = 1e2,
    rotate = TRUE
  ) +
  # limiting the map to focus on the region
  coord_sf(
    xlim = range(zoo_l_lme10$x, na.rm = TRUE),
    ylim = range(zoo_l_lme10$y, na.rm = TRUE),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "ZooMSS Pacific-Hawaiian LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
zoo_10_anim <- gganimate::animate(zoo_10_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("zoomss_hawaii_biomass.gif", animation = zoo_10_anim)

## ----------- ##
#   LME14  --- Patagonia
## ----------- ##
zoo_lme14 <- crop(zoo, bbox_lme14_ext)
zoo_df_lme14 <- as.data.frame(zoo_lme14, xy = TRUE)

# pivot long
zoo_l_lme14 <- zoo_df_lme14 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))


zoo_14_p <- ggplot() +
  geom_raster(data = zoo_l_lme14, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"),
                   na.value = "transparent") + 
  borders("world", colour = "gray85", fill = "#575151") +
  geom_contour(
    data = bathy_Patagon,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -4000),
    color = "gray50",
    size = 0.3,
    alpha = 0.3
  ) +
  geom_text_contour(
    data = bathy_Patagon,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -5000),
    size = 2.5,           
    color = "gray30",    
    skip = 1,
    min.size = 1e2,
    rotate = TRUE
  ) +
  coord_sf(
    xlim = range(zoo_l_lme14$x, na.rm = TRUE),
    ylim = range(zoo_l_lme14$y, na.rm = TRUE),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "ZooMSS Patagonian Shelf LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
zoo_14_anim <- gganimate::animate(zoo_14_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("zoomss_patagonia_biomass.gif", animation = zoo_14_anim)

## ----------- ##
#   LME21  --- Norwegian
## ----------- ##
zoo_lme21 <- crop(zoo, bbox_lme21_ext)
zoo_df_lme21 <- as.data.frame(zoo_lme21, xy = TRUE)


# pivot long
zoo_l_lme21 <- zoo_df_lme21 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

zoo_21_p <- ggplot() +
  geom_raster(data = zoo_l_lme21, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"), #of WW ??
                   na.value = "transparent") + 
  borders("world", colour = "gray85", fill = "#575151") +
  geom_contour(
    data = bathy_Norw,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -4000),
    color = "gray50",
    size = 0.3,
    alpha = 0.3
  ) +
  geom_text_contour(
    data = bathy_Norw,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -5000),
    size = 2.5,           
    color = "gray30",    
    skip = 1,
    min.size = 1e2,
    rotate = TRUE
  ) +
  coord_sf(
    xlim = range(zoo_l_lme21$x, na.rm = TRUE),
    ylim = range(zoo_l_lme21$y, na.rm = TRUE),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "ZooMSS Norwegian Sea LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
zoo_21_anim <- gganimate::animate(zoo_21_p, fps = 10, 
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("zoomss_norwegian_biomass.gif", animation = zoo_21_anim)

## ----------- ##
#   LME28  --- Guinea
## ----------- ##
zoo_lme28 <- crop(zoo, bbox_lme28_ext)
zoo_df_lme28 <- as.data.frame(zoo_lme28, xy = TRUE)

# pivot long
zoo_l_lme28 <- zoo_df_lme28 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

zoo_28_p <- ggplot() +
  geom_raster(data = zoo_l_lme28, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"), #of WW ??
                   na.value = "transparent") +
  borders("world", colour = "gray85", fill = "#575151") +
  geom_contour(
    data = bathy_Guinea,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -4000),
    color = "gray50",
    size = 0.3,
    alpha = 0.3
  ) +
  geom_text_contour(
    data = bathy_Guinea,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -5000),
    size = 2.5,           
    color = "gray50",    
    skip = 1,
    min.size = 1e3,
    rotate = TRUE
  ) +
  coord_sf(
    xlim = range(zoo_l_lme28$x, na.rm = TRUE),
    ylim = range(zoo_l_lme28$y, na.rm = TRUE),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "ZooMSS Guinea Current LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
zoo_28_anim <- gganimate::animate(zoo_28_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("zoomss_guinea_biomass.gif", animation = zoo_28_anim)
