###########################################################################################################
###                                                                                                     ###
###                         FIGURE 5: INLA MATERN CORRELATION FUNCTION                                  ###
###           Visualize the extent of spatial autocorrelation in the INLA Bernoulli GAM                 ###
###                                                                                                     ###
###########################################################################################################

#' The majority of this is heavily based on code written by Highland Statistics Ltd.
#' Code comes from Volume I and II of Beginner's Guide to Spatial, Temporal and Spatial-Temporal Ecological
#'   Data Analysis with R-INLA by Zuur, Ieno and Saveliev, and associated courses.
#' Books and courses can be found at: https://www.highstat.com


### 1. Load packages and data ------------------------------------------------------------------------------
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(INLA)
library(brinla)
library(fields)   

my.wd <- getwd()

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



### 2. Set up for plot  ---------------------------------------------------------------------------------

### i) Report the hyper-parameters on scale of the response variable  ----
bri.hyperpar.summary(M.final) 

### ii) Get all the Spatial parameters for the model  ---- 
str(spde.sp, max.level = 1) #Already have spde loaded

# By hand if want to ensure are doing transformation correctly
SpatField <- inla.spde2.result(inla = M.final,
                               name = "w",
                               spde = spde.sp,
                               do.transfer = TRUE)

Kappa <- inla.emarginal(function(x) x, SpatField$marginals.kappa[[1]] )

Sigma.u <- inla.emarginal(function(x) sqrt(x), 
                          SpatField$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x, 
                        SpatField$marginals.range.nominal[[1]] )

Range #295.7996...remember this is in METERS


### iii) Do the calculation to determine the matern correlation ----
LocMeshExt <- mesh.sp$loc[ ,1:2]

D.ext <- as.matrix(dist(LocMeshExt))

# Calculate the Matern correlation values based on parameters from model
dvec.ext <- seq(0, max(D.ext), length = 10000)      
#The length is quite large as the difference between distances in our study are so large

CorM.ext <- (Kappa * dvec.ext) * besselK(Kappa * dvec.ext, 1) 
CorM.ext[1] <- 1  #Because the first element is NaN instead of 100% correlation


### iv) Add interpretation to thea graph  -----
# STRONG spatial correlation = sites separated by distance at which is corr b/w 0.8 to 1.0
stg_corr <- which(CorM.ext >= 0.8) 
( stg_dist <- dvec.ext[max(stg_corr)] ) #~44m; 44.1848

# MODERATE spatial correlation = sites separated by distance at which is corr b/w 0.5 to 0.8
mod_corr <- which(CorM.ext >= 0.5) 
( mod_dist <- dvec.ext[max(mod_corr)] ) #~110m; 110.462

# WEAK spatial correlation = sites separated by distance at which is corr b/w 0.1 to 0.5 
weak_corr <- which(CorM.ext >= 0.1) 
( weak_dist <- dvec.ext[max(weak_corr)] ) #~298m; 298.2474
Range  #295.7996 <- this of course is the ACTUAL value the model uses

# Determine the max distance between sampling sites within an wetland, averaged across all wetlands (for reference)
wtld_dists <- ext.sed %>% 
  select(wtld, east, north) %>% 
  group_by(wtld) %>% 
  nest()  %>% 
  transmute(max_dist = map_dbl(data, m.getmaxdist)) %>%   #m.getmaxdist is in the aiv_ext_func.R script
  ungroup()

range(wtld_dists$max_dist) #Max distance between sampling sites within a wetland ranges from 144.8m to 1132.7m

( wtld_avgD <- mean(wtld_dists$max_dist) ) #~430m (429.84m) is the average maximum distance between sampling locations within a wetland



### 3. Plot Figure 3  ---------------------------------------------------------------------------------

# Make df to use with ggplot
extcorr_df <- tibble("Distance"     = dvec.ext,
                     "Correlation"  = CorM.ext)

# Distances on X-axis are defined above
Range # Range = extent of weak correlation
mod_dist
stg_dist

# Make color palates to make it easier to quickly see the strength of spatial correlation
cor.colors <- c("Strong" = "red",
                "Moderate" = "darkorange",
                "Weak" = "lightgoldenrod")

cor.grayscale <- c("Strong" = "gray50",
                   "Moderate" = "gray70",
                   "Weak" = "gray90")

# Plot
mcor_plot <- ggplot(data = extcorr_df, 
                    aes(x = Distance, y = Correlation)) +
  geom_rect(aes(xmin = 0, xmax = stg_dist, ymin = 0, ymax = 1, 
                fill = "Strong")) +  
  geom_rect(aes(xmin = stg_dist, xmax = mod_dist, ymin = 0, ymax = 1, 
                fill = "Moderate")) +
  geom_rect(aes(xmin = mod_dist, xmax = Range, ymin = 0, ymax = 1, 
                fill = "Weak")) +
  geom_line(size = 2) +
  xlim(0, 1000) +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
  scale_fill_manual(name = "Strength of \nspatial correlation",
                    values = cor.grayscale,
                    limits = c("Strong", "Moderate", "Weak")) + #To adjust the order of colors in the legend
  labs(x = "Distance between sediment sampling sites (m)") 

mcor_plot

# Save
png(filename = paste0(mywd, "/Figures/",
                      "Fig5_Matern correlation function in final model.png"),
    width = 8.5,
    height = 8,
    units = "in",
    bg = "white",
    res = 300)
p.srf
dev.off()




