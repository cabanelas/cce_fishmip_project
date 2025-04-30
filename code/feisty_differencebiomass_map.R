################################################################################
#############          FEISTY Difference MAP      #############################
#############             APR-2025                 #############################
#############     Climate Change Ecology Class     #############################
## by: Alexandra Cabanelas 
################################################################################
# FEISTY ISIMIP data downloaded APR 2025

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse)
library(terra) #rast()
library(rnaturalearth)
library(sf)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
feisty <- rast("feisty_stack_output.tif")
land <- ne_countries(scale = "medium", returnclass = "sf")
print(feisty)

## ------------------------------------------ ##
#            Calculation -----
## ------------------------------------------ ##
biomass_1950 <- feisty[[1]]
biomass_2014 <- feisty[[65]]

# compute the difference: 2014 - 1950
biomass_diff <- biomass_2014 - biomass_1950

# convert to df for ggplot2
diff_df <- as.data.frame(biomass_diff, xy = TRUE, na.rm = TRUE)
colnames(diff_df)[3] <- "diff"

## ------------------------------------------ ##
#            Plot -----
## ------------------------------------------ ##
(worldfeisty <- ggplot(diff_df, aes(x = x, y = y, fill = diff)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "gray70"
  ) +
  borders("world", colour = "gray40", fill = "gray70") +
  coord_fixed() +
  labs(title = "Change in TCB (2014 − 1950)",
       fill = "Biomass Δ (g/m²)", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black", 
                                 margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_x_continuous(expand = c(0, 0),
                     labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + #breaks = seq(-72, -69, 1))
  scale_y_continuous(expand = c(0, 0),
                     labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N")))
)

#ggsave("worldmapfeisty_2014minus1950.png", plot = worldfeisty, width = 10, height = 8, dpi = 300, bg = "white")

## ------------------------------------------ ##
#            Plot with 4 LME polygons -----
## ------------------------------------------ ##
lme <- st_read("raw/lme66/LME66.shp") 

subset_lmes <- lme %>% 
  filter(LME_NUMBER %in% c(10, 14, 21, 28))

(worldmapfeisty <- ggplot() +
  geom_raster(data = diff_df, aes(x = x, y = y, fill = diff)) +
  geom_sf(data = subset_lmes, fill = NA, color = "black", size = 0.5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, na.value = "gray70"
  ) +
  borders("world", colour = "gray40", fill = "gray70") +
  #coord_fixed() +
  labs(title = "Change in TCB (2014 − 1950)",
       fill = "Biomass Δ (g/m²)", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15, color = "black", 
                                 margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  scale_x_continuous(expand = c(0, 0),
                     labels = function(x) paste0(abs(x), "°", 
                                                 ifelse(x < 0, "W", "E"))) + #breaks = seq(-72, -69, 1))
  scale_y_continuous(expand = c(0, 0),
                     labels = function(y) paste0(abs(y), "°", 
                                                 ifelse(y < 0, "S", "N")))
)

#ggsave("worldmapfeisty_2014minus1950_lmepolyg.png", plot = worldmapfeisty, width = 10, height = 8, dpi = 300, bg = "white")