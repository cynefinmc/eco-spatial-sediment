###########################################################################################################
###                                                                                                     ###
###               FIGURE 7: PLOT KULLDORFF'S SPATIAL SCAN CLUSTERS OVERLAYING STUDY AREA                ### 
###                                                                                                     ###
###                                                                                                     ###
###########################################################################################################


### 1. Load packages and data ------------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggmap)
library(ggplot2)
library(ggrepel)
library(sf)
library(viridis)
library(scales)

my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

# Read in shapefile results file from SatScan for plotting
col.shp <- st_read(paste0(my.wd,
                          "/Results/SaTScan-results/Ext-sed-all-satscan-results.col.shp"))
table(st_geometry_type(col.shp)) #'Should have 12 polygons


### 2. Set up for plotting -------------------------------------------------------------------------------------------------

# Get extra wide basemap
remake.basemap <- F
if(remake.basemap){
  bbox.mat <- data.frame(min = c(-123.5500, 48.8500),
                         max = c(-121.6005, 49.42967))
  bbox.mat <- as.matrix(bbox.mat)
  rownames(bbox.mat) <- c("x", "y")
  
  basemap.satscan <- get_map(bbox.mat, maptype = "terrain")
  save(basemap.satscan, file = paste0(my.wd, "/Figures/Basemap_terrain_ggmap_obj_for_satscan_map.Rdata"))
}else{
  load(paste0(my.wd, "/Figures/Basemap_terrain_ggmap_obj_for_satscan_map.Rdata")) #basemap.satscan
}

# Set up location data for wetlands - plot with lat/long to match SatScan results
wtld.locs <- ext.sed %>% group_by(wtld) %>% summarize(X = first(longitude),
                                                     Y = first(latitude))

# Set up shapefile with extra columns for visualizing on plot
shp.sign <- col.shp %>% 
              mutate("significance" = if_else(P_VALUE >= 0.05,
                                            "not_sign",
                                            "sign"),
                     "direction" = if_else(REL_RISK > 1,
                                         "high",
                                         "low")) %>% 
              filter(significance == "sign")

# Transform shp.sign to ensure it has the same names for lat and lon as the basemap
attributes(shp.sign)$names <-str_replace(attributes(shp.sign)$names, "LATITUDE", "lat")
attributes(shp.sign)$names <-str_replace(attributes(shp.sign)$names, "LONGITUDE", "lon")

attributes(shp.sign)$names

# Set up colors
rr.cols <- c("high" = "red",
             "low" = "blue")



### 3. Plot ------------------------------------------------------------------------------------------------------------------------

# Set random seed so the layout of the text on the plot is repeatable
set.seed(42)

# Plot
cluster.plot <- ggmap(basemap.satscan, darken = c(0.5, "white")) +
  geom_jitter(data = wtld.locs,      #jitter slightly so it's easier to distinguish wetland centroids from cluster centroids
              aes(x = X, y = Y),
              shape = 17, color = "grey27") +
  geom_text_repel(data = wtld.locs,
                  aes(x = X, y = Y,
                      label = wtld),
                  color = "grey27",
                  size = 3.5,
                  nudge_y = -0.02) +
  geom_sf(data = shp.sign,  #polygons of the significant clusters
          aes(color = direction, fill = direction),
          alpha = 0.2) +   
  geom_label_repel(data = shp.sign,
                   aes(x = lon, y = lat, 
                       label = CLUSTER),
                   box.padding = 0.5,
                   min.segment.length = 0,
                   nudge_y = 0.03,
                   nudge_x = 0.01,
                   size = 4.5) +
  theme_bw() +
  scale_fill_manual(values = rr.cols) +
  scale_color_manual(values = rr.cols) +
  guides(size = guide_legend(override.aes = aes(label = ""))) + 
  labs(x = "Longitude",
       y = "Latitude")


cluster.plot

# Save
png(filename = paste0(my.wd, 
                      "/Figures/Fig7_Extensive sediment significant spatial scan clusters with wtld names and basemap.png"),
    width = 10.5,
    height = 5,
    units = "in",
    res = 300)
cluster.plot
dev.off()





