#################################################################################
###                                                                           ###
###   FIGURE 4: VISUALIZE FINAL INLA MODEL - ENTIRE SPATIAL RANDOM FIELD      ###
###                                                                           ###
###                                                                           ###
#################################################################################


### 1. Set up environment -------------------------------------------------------------------------------------

### i. Load Packages ----
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(viridis)
library(scales)    
library(INLA)
library(fields)    #' Need this for the PlotField2 function to work

#' Load data
my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

# Load final model and associated INLA objects for the SRF 
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_gam_spatial_meshsp_n21958_finalmodel.rdata")) #M.final

N <- nrow(ext.sed) #Sample size

# Load in the INLA objects need to plot the spatial field and to transform to log scale
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_meshsp_obj_fromdesktop_n21958.rdata"))   #mesh.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_Asp_obj_fromdesktop_n21958.rdata"))      #A.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_SPDEsp_obj_fromdesktop_n21958.rdata"))   #spde.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_windexsp_obj_fromdesktop_n21958.rdata")) #w.index.sp
load(paste0(my.wd, "/Results/Final-inla-objects/", "Extsed_Stack1_obj_fromdesktop_n21958.rdata"))   #Stack1



### 2. Make plot -------------------------------------------------------------------------------------

w.pm <- M.final$summary.random$w$mean

p.srf <- m.gg2PlotField(field = w.pm,
               mesh = mesh.sp,
               ColourScheme = 3,                #You can use the values 1 to 6; note the 'u' in color!!!
               xlim = range(mesh.sp$loc[ ,1]), #Either need to add xlim and ylim, OR add a ContourMap
               ylim = range(mesh.sp$loc[ ,2]),
               scalemin = min(w.pm),    #Provide an input for min value for colorbar, to better compare between different plots
               scalemax = max(w.pm),    #Provide an input for max value for colorbar, to better compare between different plots
               Xlabel = "Easting",
               Ylabel = "Northing",
               LegendName = "Spatial \nrandom \neffect",
               #AddPoints = T,                 #Default is 'F'
               TextSize = 14,
               MyTitle = "") 

p.srf


### 3. Save with a specified size and resolution ------------------------------------------------------
png(filename = paste0(mywd, "/Figures/",
                      "Fig4_INLA posterior mean of spatial random field of entire study area.png"),
    width = 10.5,
    height = 6,
    units = "in",
    bg = "white",
    res = 300)
p.srf
dev.off()



