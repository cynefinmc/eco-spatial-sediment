###########################################################################################################
###                                                                                                     ###
###        FIGURE 3: EBIRD DIVERSITY AND AIV MATRIX POSITIVITY BY WETLAND PER STUDY YEAR                ###
###                                                                                                     ###
###########################################################################################################


### 1. Load packages and data ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel) 
library(scales)
library(viridis)

#' Data
my.wd <- getwd()
ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))
ext.sed$study_yr <- factor(ext.sed$study_yr,
                           levels = c("yr1", "yr2"))



### 2. Prepare data for plotting ------------------------------------------------------------------------------

# Subset to ebird data and change unit of analysis to wetland-level
#   as no ebird data at level of sampling sites
ebird.summary <- ext.sed %>% 
  group_by(wtld, mon_yr, study_yr) %>% 
  summarize("Matrix_positivity" = sum(matrix01 == 1) / n(),
            #"latitude_wtld" = mean(latitude),   #For plotting with ggmaps
            #"longitude_wtld" = mean(latitude),  #For plotting with ggmaps
            "easting_wtld" = mean(east),
            "northing_wtld" = mean(north),
            "N_checklist" = unique(N_checklist),
            "Total_abundance" = unique(Total_abundance),
            "Relative_abundance" = unique(Relative_abundance),
            "Species_richness" = unique(Species_richness)) %>% 
  #' Format units for visualization on plot
  mutate("Matrix_positivity" = format(round(Matrix_positivity, digits = 2), nsmall = 2)) %>% 
  #' Add column to indicate if wetland has ebird data or not
  mutate("Bird_data" = case_when(!is.na(N_checklist) ~ "yes",
                                 is.na(N_checklist) ~ "no"))


# Turn wetlands into factor
ebird.summary$wtld <- factor(ebird.summary$wtld,
                             levels = c("BHER", "BRY", "BURN", "LMT", "MIL",
                                        "MNP", "NIC", "PIT", "REI", "SAR",
                                        "SDW", "SFEN", "SIS", "STA", "TRL",
                                        "TUY", "TWU", "WCK"))

# Replace missing values for wetlands without bird data
ebird.summary$Species_richness[is.na(ebird.summary$Species_richness) == T] <- 0
ebird.summary$Relative_abundance[is.na(ebird.summary$Relative_abundance) == T] <- 0

# Make a wetland-matrix positivity label
ebird.summary <- mutate(ebird.summary,
                   "wtld_label" = paste0(wtld, " (", Matrix_positivity,")"))


### 3. Plot  ---------------------------------------------------------------------------------------------------

# Set up objects for graph
bird_pch <- c("yes" = 21, 
              "no" = 8)

# Set random seed to make sure ggrepel makes the same graph consistently
set.seed(63) 

# Plot
bird.plot <- ggplot(ebird.summary, aes(x = easting_wtld, y = northing_wtld)) + 
  geom_point(alpha = 0.9,  
             aes(size = Species_richness,
                 fill = Relative_abundance,
                 shape = Bird_data)) +
  geom_text_repel(aes(label = wtld_label,
                      size = 5),
                  box.padding = 0.5) +  
  facet_grid(study_yr ~.) +
  theme_bw() +
  scale_fill_viridis(discrete = F, option = "A",
                     name = "Relative abundance") +
  scale_shape_manual(values = bird_pch,
                     name = "Species richness") +
  theme(text = element_text(size = 14)) +
  guides(size = guide_legend(override.aes = aes(label = "")),  #To remove the 'a' in the 'size' legend
         shape = guide_legend("Has eBird data")) + 
  labs(x = "Easting", y = "Northing")

  
bird.plot

# Save
png(filename = paste0(my.wd, "/Figures/",
                      "Fig3_Ebird_diversity_richness_rel_abund_and_pos_rate.png"),
    width = 10.5,
    height = 7,
    units = "in",
    bg = "white",
    res = 300)
bird.plot
dev.off()






