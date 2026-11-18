###########################################################################################################
###                                                                                                     ###
###               SPATIAL SCAN MODELING OF EXTENSIVE SEDIMENT MATRIX GENE PCR RESULTS                   ###
###                                                                                                     ###
###     This script prepares data of AIV detection in sediment for spatial clustering scan              ###
###     using SatScan software and then summarizing output                                              ###
###                                                                                                     ###
###########################################################################################################

### 1. Load the files and packages --------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(stringr)

my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

### 2. PRE-ANALYSIS: Prep data to run in SaTScan ------------------------------------------------------------

#' Add month and study year covariates as numeric values
ext.sed <- ext.sed %>% 
              mutate("Nmon" = lubridate::month(date_col),
                     "Nstudy_yr" = as.numeric(factor(study_yr)))

#' Separate data by study year and then into cases (matrix01 == 1) vs "controls" (matrix01 == 0) 
#' Keep lat/long for coordinates so can plot output using ggmaps
sed.cases <- ext.sed %>% 
                  filter(matrix01 == 1) %>% 
                  select(wtld, stID, matrix01, date_col, Nmon, Nstudy_yr, east, north)

sed.controls <- ext.sed %>% 
                  filter(matrix01 == 0) %>% 
                  select(wtld, stID, matrix01, date_col, Nmon, Nstudy_yr, east, north)

sed.spatlocs <- ext.sed %>% 
                  select(wtld, stID, east, north, longitude, latitude)

#' Save as tab delimited file
write.table(sed.cases, paste0(my.wd, "/Extensive_sediment_satscan_data_all_cases.txt"),
            row.names = F, quote = F)
write.table(sed.controls, paste0(my.wd, "/Extensive_sediment_satscan_data_all_controls.txt"),
            row.names = F, quote = F)
write.table(sed.spatlocs, paste0(my.wd, "/Extensive_sediment_satscan_data_all_coords.txt"),
            row.names = F, quote = F)

### 3. Run SatScan software -----------------------------------------------------------------------------------------

#' Cluster analysis performed in SatScan software version 10.0
#' Summary of cluster results from `Ext-sed-all-satscan-results.txt` output file
#' Plot clusters from shapefile `Ext-sed-all-satscan-results.gis.shp`





