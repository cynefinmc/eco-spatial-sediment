####################################################################################################
###                                                                                              ###
###             Figure 1: Map of extensive sediment sampling locations                           ###
###                                                                                              ###
###                                                                                              ###
####################################################################################################


### 1. Load packages and data -----------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tmaptools)
library(tmap)
library(sf)
library(maptiles) #' Supplies basemaps for `tmap` package
library(viridis)


#' Dataset with coordinates (both as UTM and lat/long for plotting with basemap)
my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))


### 2. Wetland locations for plotting -------------------------------------------------------------------------------

### i) Summarize to a single lat/long per wetland ----
#' Doesn't matter which sampling site is used as map is so zoomed out
ext.locs <- ext.sed %>% 
  group_by(wtld) %>% 
  summarize(long = first(longitude),
            lat = first(latitude)) %>% 
  ungroup()

ext.locs <- st_as_sf(x = ext.locs,
                     coords = c("long", "lat"),
                     crs = "+proj=longlat +datum=WGS84")


### ii) Make a wider bounding box for WETLANDS to zoom plot out ----
# Summarize the extent of all sediment UTM coordinates and zoom out by 10KM
utm_box <- ext.sed %>% 
  summarize(min_x_east = min(east) - 10000,
            min_y_north = min(north) - 10000,
            max_x_east = max(east) + 10000,
            max_y_north = max(north) + 10000)
utm_box
# min_x_east  min_y_north  max_x_east  max_y_north
# 476376      5418716      601475      5476163   (with +10km)

# Now need to turn this into a dataframe to be able to convert it
utm_df <- data.frame(x = c(476376, 601475),
                     y = c(5418716, 5476163))

### iii) Turn into spatial points then bounding box
#' I measured WGS84 (in metric units) and the UTM zone is 10N whose's EPSG code is 32610
sp_utm <- st_as_sf(x = utm_df,
                   coords = c("x", "y"),
                   crs = "+proj=utm +zone=10 +datum=WGS84")

# Now transform the coordinate system into lat/long used by google maps
sp_LL <- st_transform(x = sp_utm,
                         crs = "+proj=longlat +datum=WGS84")

#' Turn lat long spatial object into a character string that matches tmap requirements
#'   format: c(xmin, ymin, xmax, ymax)
bbox_str <- sf::st_bbox(sp_LL)


### 3. Make map ------------------------------------------------------------------------------------------------------
base.map <- tm_shape(bbox = bbox_str) +
              tm_basemap("OpenStreetMap", alpha = 0.5)
base.map

ext.map <- tm_shape(bbox = bbox_str) +
            tm_basemap("OpenStreetMap", alpha = 0.6) +
            tm_shape(ext.locs) +
            tm_dots(fill = "#414487FF",
                    col = "black",
                    shape = 21,
                    lwd = 0.8,
                    size = 1) +
            tm_scalebar(position = c("right", "bottom")) +
            tm_compass(position = c("left", "top")) #+
            #tm_layout(legend.position = c("left", "bottom"),
            #          legend.bg.color = "white",
            #          legend.bg.alpha = 0.5) 
          
ext.map

png(filename = paste0(my.wd, "/Figures/Fig1_Locations_of_Extensive_Sediment_Wetlands.png"),
    width = 10.5,
    height = 5,
    units = "in",
    res = 300)
ext.map
dev.off()










