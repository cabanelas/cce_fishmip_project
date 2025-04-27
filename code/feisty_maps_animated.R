################################################################################
#############             FEISTY  GIF MAPS TS      #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################
# FEISTY ISIMIP data downloaded APR 2025

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
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
#library(rnaturalearthdata)
#library(ggnewscale)
#library(cowplot)
#library(usmap)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
feisty <- rast("feisty_stack_output.tif")

# LME
lme <- st_read("raw/lme66/LME66.shp") %>%
  mutate(LME_focus = case_when(
    LME_NUMBER == 10 ~ "Insular Pacific-Hawaiian",
    LME_NUMBER == 14 ~ "Patagonian Shelf",
    LME_NUMBER == 21 ~ "Norwegian Sea",
    LME_NUMBER == 28 ~ "Guinea Current",
    TRUE ~ NA_character_
  ))

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
#            FEISTY map -----
## ------------------------------------------ ##
plot(feisty$"1950")

## ----------- ##
#   LME10  --- Hawaii
## ----------- ##
feisty_lme10 <- crop(feisty, bbox_lme10_ext)
feisty_df_lme10 <- as.data.frame(feisty_lme10, xy = TRUE)

# pivot long
feisty_long_lme10 <- feisty_df_lme10 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

feisty_10_p <- ggplot() +
  geom_raster(data = feisty_long_lme10, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"), #of WW ??
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
    #stroke = 0.1,          # thin white outline for text
    size = 2.5,           
    color = "gray60",    
    skip = 1,
    min.size = 1e2,
    rotate = TRUE
  ) +
  # limiting the map to focus on the region
  coord_sf(xlim = c(-178, -152), ylim = c(16, 30), expand = FALSE) +
  #coord_sf(
  #  xlim = c(bbox_lme10["xmin"], bbox_lme10["xmax"]),
  #  ylim = c(bbox_lme10["ymin"], bbox_lme10["ymax"])
  #)
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(title = "FEISTY Pacific-Hawaiian LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + #breaks = seq(-72, -69, 1))
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
feisty_10_anim <- gganimate::animate(feisty_10_p, fps = 10, #frames per second
                #width = 800, height = 600, 
                width = 1600, 
                height = 1200, 
                res = 150, 
                renderer = gifski_renderer())

anim_save("FEISTY_hawaii_biomass.gif", animation = feisty_10_anim)

# video file instead of GIF -- have to install ffmpeg to computer - so didnt do
#gganimate::ffmpeg_renderer()) 

## ----------- ##
#   LME14  --- Patagonia
## ----------- ##
feisty_lme14 <- crop(feisty, bbox_lme14_ext)
feisty_df_lme14 <- as.data.frame(feisty_lme14, xy = TRUE)

# pivot long
feisty_long_lme14 <- feisty_df_lme14 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

feisty_14_p <- ggplot() +
  geom_raster(data = feisty_long_lme14, aes(x = x, y = y, fill = biomass)) +
  scale_fill_scico(palette = "cork", 
                   direction = -1, 
                   begin = 0,
                   end = 0.4,
                   name = expression("TCB Biomass (g/m"^2*")"), #of WW ??
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
    xlim = range(feisty_long_lme14$x, na.rm = TRUE),
    ylim = range(feisty_long_lme14$y, na.rm = TRUE),
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
  labs(title = "FEISTY Patagonian Shelf LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
feisty_14_anim <- gganimate::animate(feisty_14_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("FEISTY_patagonia_biomass.gif", animation = feisty_14_anim)

## ----------- ##
#   LME21  --- Norwegian
## ----------- ##
feisty_lme21 <- crop(feisty, bbox_lme21_ext)
feisty_df_lme21 <- as.data.frame(feisty_lme21, xy = TRUE)

# pivot long
feisty_long_lme21 <- feisty_df_lme21 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

feisty_21_p <- ggplot() +
  geom_raster(data = feisty_long_lme21, aes(x = x, y = y, fill = biomass)) +
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
    xlim = range(feisty_long_lme21$x, na.rm = TRUE),
    ylim = range(feisty_long_lme21$y, na.rm = TRUE),
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
  labs(title = "FEISTY Norwegian Sea LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
feisty_21_anim <- gganimate::animate(feisty_21_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("FEISTY_norwegian_biomass.gif", animation = feisty_21_anim)

## ----------- ##
#   LME28  --- Guinea
## ----------- ##
feisty_lme28 <- crop(feisty, bbox_lme28_ext)
feisty_df_lme28 <- as.data.frame(feisty_lme28, xy = TRUE)

# pivot long
feisty_long_lme28 <- feisty_df_lme28 %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "biomass"
  ) %>% 
  mutate(year = as.integer(year))

feisty_28_p <- ggplot() +
  geom_raster(data = feisty_long_lme28, aes(x = x, y = y, fill = biomass)) +
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
    xlim = range(feisty_long_lme28$x, na.rm = TRUE),
    ylim = range(feisty_long_lme28$y, na.rm = TRUE),
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
  labs(title = "FEISTY Guinea Current LME Biomass {frame_time}",
       x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + 
  scale_y_continuous(labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N"))) +
  transition_time(year) +
  ease_aes('linear') 

# GIF
feisty_28_anim <- gganimate::animate(feisty_28_p, fps = 10, #frames per second
                                     width = 1600, 
                                     height = 1200, 
                                     res = 150, 
                                     renderer = gifski_renderer())

anim_save("FEISTY_guinea_biomass.gif", animation = feisty_28_anim)


## ------------------------------------------ ##
#            All LMEs world map -----
## ------------------------------------------ ##
ggplot() +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#cce6ff", color = NA),  # light blue ocean
    plot.background = element_rect(fill = "#cce6ff", color = NA)
  ) +
  borders("world", fill = "gray80", color = "gray60") +   # land
  geom_sf(data = lme, fill = NA, color = "black", size = 0.3)

bathy_world <- getNOAA.bathy(
  lon1 = -180, lon2 = 180,
  lat1 = -80, lat2 = 90,
  resolution = 10  
)
bathy_world_df <- fortify.bathy(bathy_world)

(all_lme_p <- ggplot() +
  geom_raster(data = bathy_world_df, aes(x = x, y = y, fill = z),
              show.legend = FALSE) +
  ggnewscale::new_scale_fill() +  # reset fills
  borders("world", fill = "gray80", colour = "gray60") +
  geom_sf(data = filter(lme, !is.na(LME_focus)), aes(fill = LME_focus), color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(
      "Insular Pacific-Hawaiian" = "#FFD700",
      "Patagonian Shelf" = "#FF4500",
      "Norwegian Sea" = "#00FF7F",
      "Guinea Current" = "#984EA3"
    ),
    breaks = c(
      "Insular Pacific-Hawaiian",
      "Patagonian Shelf",
      "Guinea Current",
      "Norwegian Sea"
    )) +
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12))
)
#ggsave("figures/all_lme_world.png", plot = all_lme_p, width = 10, height = 8, dpi = 300,bg = "white")

# same but showing all other LME polygons
(all_lme_p_v2 <-ggplot() +
  geom_raster(data = bathy_world_df, aes(x = x, y = y, fill = z),
              show.legend = FALSE) +
  ggnewscale::new_scale_fill() +  # reset fills
  borders("world", fill = "gray80", colour = "gray60") +
  geom_sf(data = lme, fill = NA, color = "black", size = 0.2) +
  geom_sf(data = filter(lme, !is.na(LME_focus)), aes(fill = LME_focus), color = "black", size = 0.3) +
  scale_fill_manual(
    values = c(
      "Insular Pacific-Hawaiian" = "#FFD700",
      "Patagonian Shelf" = "#FF4500",
      "Norwegian Sea" = "#00FF7F",
      "Guinea Current" = "#984EA3"
    ),
    breaks = c(
      "Insular Pacific-Hawaiian",
      "Patagonian Shelf",
      "Guinea Current",
      "Norwegian Sea"
    )) +
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 90), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12))
)
#ggsave("figures/all_lme_worldv2.png", plot = all_lme_p_v2, width = 10, height = 8, dpi = 300,bg = "white")


## ------------------------------------------ ##
#            THE END -----
## ------------------------------------------ ##



################################################################
########  TRASHHH     ##########################################
################################################################
ggplot(feisty_long_lme10, aes(x = x, y = y, fill = biomass)) +
  geom_tile() +
  borders("world", colour = "gray30", fill = "gray90", size = 0.2) +
  scale_fill_viridis_c(name = "Biomass (g m²)", option = "C", na.value = "white") +
  scale_x_continuous(breaks = seq(-180, 180, by = 30), name = "Longitude") +
  scale_y_continuous(breaks = seq(-90, 90, by = 30), name = "Latitude") +
  coord_sf(xlim = c(-178, -152), ylim = c(16, 30), expand = FALSE) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  labs(title = "FEISTY Biomass (g WW/m²) — Year: {frame_time}") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold")
  ) 

library(ggOceanMaps)
basemap_lme10 <- ggOceanMaps::basemap(
  limits = c(bbox_lme10["xmin"], bbox_lme10["xmax"], bbox_lme10["ymin"], bbox_lme10["ymax"]),
  land.col = "gray",
  land.border.col = "black",  # no border
  bathymetry = FALSE
)

# works but meh
# bounding box for the main map area
bbox_lme10 <- st_as_sfc(st_bbox(c(xmin = -179.9, ymin = 18,
                                  xmax = -152, ymax = 30), 
                                crs = 4326))

basemap_lme10 +
  geom_tile(data = feisty_long_lme10, aes(x = x, y = y, fill = biomass), 
            inherit.aes = FALSE) +
  geom_sf(data = land, fill = "gray70", color = NA) +
  geom_sf(data = bbox_lme10, fill = NA, color = "black", linewidth = 1.5) +
  coord_sf(ylim = c(18, 30),
           xlim = c(-179.9, -152),
           expand = TRUE)

# raster covers land ; no land
basemap_lme10 +
  geom_raster(data = feisty_long, aes(x = x, y = y, fill = biomass), 
              inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "C") +
  labs(title = "LME 10 Biomass {frame_time}", x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  transition_time(year) +
  ease_aes('linear')

# tile covers land ; no land
basemap_lme10 +
  geom_tile(data = feisty_long, aes(x = x, y = y, fill = biomass), inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "C", na.value = "transparent")

# raster covers land ; no land
basemap_lme10 +
  geom_raster(data = feisty_long, aes(x = x, y = y, fill = biomass), inherit.aes = FALSE) +
  coord_sf(
    xlim = c(bbox_lme10["xmin"], bbox_lme10["xmax"]),
    ylim = c(bbox_lme10["ymin"], bbox_lme10["ymax"]),
    expand = FALSE
  ) +
  scale_fill_viridis_c(option = "C") +
  theme_minimal()