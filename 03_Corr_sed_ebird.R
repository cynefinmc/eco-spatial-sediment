###########################################################################################################
###                                                                                                     ###
###            CORRELATION OF EXTENSIVE SEDIMENT MATRIX RESULTS AND EBIRD OBSERVATIONS                  ###
###                                                                                                     ###
###     This script examines the Pearson's correlation between                                          ###
###     rate of matrix gene PCR-positive samples by wetland and month-year of collection to             ###
###     the species richness and relative abundance of 32 species of local wild waterfowl               ###
###                                                                                                     ###
###########################################################################################################

### 1. Load packages and data ------------------------------------------------------------------------------
library(dplyr)
library(tibble)

my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))
ext.sed$study_yr <- factor(ext.sed$study_yr,
                           levels = c("yr1", "yr2"))

### 2. Summarize by wetland-study year combination  ------------------------------------------------------------
#' As eBird data was summarized at level of wetland within the same month/year combination as sediment collection

ebird.summary <- ext.sed %>% 
                    group_by(wtld, mon_yr, study_yr) %>% 
                    summarize("Matrix_positivity" = sum(matrix01 == 1) / n(),
                              "latitude_wtld" = mean(latitude),   #To plot with ggmaps
                              "longitude_wtld" = mean(latitude),  #To plot with ggmaps
                              "easting_wtld" = mean(east),
                              "northing_wtld" = mean(north),
                              "N_checklist" = unique(N_checklist),
                              "Total_abundance" = unique(Total_abundance),
                              "Relative_abundance" = unique(Relative_abundance),
                              "Species_richness" = unique(Species_richness)) 


#' Remove those without eBird observations
ebird.summary <- filter(ebird.summary, !is.na(N_checklist))



### 3. Calculate Pearson's correlation  -------------------------------------------------------------------------------------------------
plot(ebird.summary$Relative_abundance, ebird.summary$Matrix_positivity) #I don't see a non-linear relationship, so ok to use Pearsons
plot(ebird.summary$Species_richness, ebird.summary$Matrix_positivity)  #I don't see a non-linear relationship, so ok to use Pearsons

cor.test(ebird.summary$Relative_abundance, ebird.summary$Matrix_positivity,
         method = "pearson")
#Correlation = 0.197 (-0.245, 0.571); not correlated

cor.test(ebird.summary$Species_richness, ebird.summary$Matrix_positivity,
         method = "pearson")
#Correlation = 0.374 (-0.056, 0.687); not correlated


