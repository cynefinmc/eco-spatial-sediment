############################################################################################
###                                                                                      ###
###   FIGURE 6: VISUALIZE FINAL INLA MODEL - SPATIAL RANDOM FIELD ZOOMED TO WETLANDS     ###
###                                                                                      ###
###                                                                                      ###
############################################################################################


### 1. Set up environment -------------------------------------------------------------------------------------

### i. Load Packages ----
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(viridis)
library(scales)    
library(INLA)
library(fields)    

my.wd <- getwd()
source(paste0(my.wd, "/00_aiv_ext_functions.R"))

#' Load data
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

# Load final model and associated INLA objects for the SRF 
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_gam_spatial_meshsp_n21958_finalmodel.rdata")) #M.final

N <- nrow(ext.sed) #Sample size

load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_meshsp_obj_fromdesktop_n21958.rdata"))   #mesh.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_Asp_obj_fromdesktop_n21958.rdata"))      #A.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_SPDEsp_obj_fromdesktop_n21958.rdata"))   #spde.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_windexsp_obj_fromdesktop_n21958.rdata")) #w.index.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_Stack1_obj_fromdesktop_n21958.rdata"))   #Stack1


### 2. Make plot for posterior mean of SRF, zoomed in to area around sampling sites from each wetland ---------------------------------------

### (i) Get bbox for each wetland, with some wiggle room around collection locations ----
bbox.list <- list()
wtld.sites <- sort(unique(ext.sed$wtld))

for(i in seq_along(wtld.sites)){
  
  mybbox <- m.MakeBbox(Dataset = ext.sed,
                       CoordNames = c("east", "north"),
                       Buffer = 0.2,
                       LocCol = "wtld",
                       Location = wtld.sites[i],
                       MakeSpatial = F,
                       MyCRS = NULL)
  
  coord.df <- data.frame("min_X" = mybbox$X[1],
                         "max_X" = mybbox$X[2],
                         "min_Y" = mybbox$Y[1],
                         "max_Y" = mybbox$Y[2],
                         "wtld"  = wtld.sites[i])

  bbox.list[[i]] <- coord.df
}

bbox.df <- bind_rows(bbox.list)

rm(bbox.list, coord.df, i)


### (ii) Get w's for the entire spatial field ----
w.pm <- M.final$summary.random$w$mean


### (iii) Make data frame of sampling location points for each wetland ----
sampling.locs <- select(ext.sed,
                        wtld, stID, east, north)


### (iv) Make the spatial projector field, by wetland ----
field.list <- list()

for(i in seq_along(wtld.sites)){
  # Subset coords by wetland
  mycoords <- bbox.df[bbox.df$wtld == wtld.sites[i], ]
  
  myproj <- inla.mesh.projector(mesh.sp, 
                                xlim = c(mycoords$min_X, mycoords$max_X), 
                                ylim = c(mycoords$min_Y, mycoords$max_Y), 
                                dims = c(300, 300)) 
  
  field.proj <- inla.mesh.project(myproj, w.pm)
  
  grid.df <- expand.grid(x = myproj$x, y = myproj$y)
  grid.df$w.pm <- as.vector(field.proj)
  grid.df$wtld <- wtld.sites[i]
  
  field.list[[i]] <- grid.df
}

field.df <- bind_rows(field.list)

rm(field.list, mycoords, myproj, field.proj, grid.df, i)


### (v) Plot a faceted version of whole study area, with free y and x scales ----

# Determine limits for color gradient bar
scalemin <- min(field.df$w.pm) 
scalemax <- max(field.df$w.pm)

# Make plot
#' Note that it throws a warning with geom_raster but still works
p.srf <- ggplot() +
            geom_raster(data = field.df, aes(x = x, y = y, fill = w.pm)) +
            geom_point(data = sampling.locs, aes(x = east, y = north),
                       pch = 21,
                       col = "grey93",
                       size = 2) +
            facet_wrap(~ wtld, scales = "free") +
            #coord_fixed(ratio = 1) +  #' doesn't work with scales = "free"
            theme_bw() + 
            scale_fill_viridis(name = "Spatial \nrandom \neffect",
                               na.value = "transparent",
                               option = "D",
                               limits = c(scalemin, scalemax)) +
            theme(axis.text.x = element_text(angle = 90), 
                  text = element_text(size = 14),
                  strip.text = element_text(size = 15)) +  
            labs(x = "Easting", y = "Northing")

p.srf


### (vi) Save image ----
png(filename = paste0(mywd, "/Figures/",
                      "Fig6_INLA posterior mean of spatial random field by wetland.png"),
    width = 10,
    height = 8,
    units = "in",
    bg = "white",
    res = 300)
p.srf
dev.off()



