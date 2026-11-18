####################################################################################
###                                                                             ###
###      ECOLOGICAL REGRESSION MODEL OF EXTENSIVE SEDIMENT MATRIX RESULTS       ###
###                                                                             ###
###      Table 1: Summarize wild waterfowl observations from eBird              ###
###                                                                             ###
###################################################################################


### 1. Set up working environment -------------------------------------------------
library(dplyr)
library(tibble)
library(tidyr)

my.wd <- getwd()
source(paste0(my.wd, "/00_aiv_ext_functions.R"))
# loads species.key that links 32 waterfowl species common names, scientific names and abbreviations
# Exclude hybrids they aren't listed in the filtered ebird object


# Load subset ebird data; only has aggregated observations for ALL month-years of extensive wetlands
# within a 2.5km radius of the centroid of sediment collection for that wetland
# Object contains subset of proprietary data
# To request access to raw data, contact eBird directly 
# https://science.ebird.org/en/use-ebird-data/download-ebird-data-products
load(paste0(ebird.dir, "/Ebird_data_subset_to_Extsed_2500m_buffer.rdata")) # ebird_ext
ebird_ext <- select(ebird_ext, -geometry)


### 2. Summarize by species ------------------------------------------------------------------------------

# Summarize ebird observations by study year and bird species
ebird.tab <- ebird_ext %>% 
                filter(!is.na(checklist_id)) %>%   #To remove mon_yr-wtld combos without ebird obs
                group_by(scientific_name, year) %>% 
                summarize(n = sum(observation_count)) %>% 
                pivot_wider(names_from = year,
                            values_from = n)

# Add common names for species
species.tab <- select(species.key, -spp)
ebird.tab <- left_join(ebird.tab, species.tab, by = "scientific_name")

ebird.tab <- select(ebird.tab, common_name, everything())
ebird.tab <- arrange(ebird.tab, common_name)

# Add a row with totals at the bottom
totals.obs <- c("Total bird observations", "", sum(ebird.tab$`2016`), sum(ebird.tab$`2017`))

# Add a row with total number of checklists by year (i.e., for relative abundance)
chk.tab <- ebird_ext %>% 
              filter(!is.na(checklist_id)) %>%   #To remove mon_yr-wtld combos without ebird obs
              group_by(year) %>% 
              summarize(n_checklist = n_distinct(checklist_id)) %>% 
              pivot_wider(names_from = year,
                          values_from = n_checklist)

totals.chk <- c("Total number of checklists", "", chk.tab$`2016`, chk.tab$`2017`)


### 3. Formating for final table ------------------------------------------------------------------------------

# Make column names legible for final table
names(ebird.tab) <- c("Common name", "Scientific name", "Fall 2016", "Fall 2017")
names(totals.obs) <- c("Common name", "Scientific name", "Fall 2016", "Fall 2017")
names(totals.chk) <- c("Common name", "Scientific name", "Fall 2016", "Fall 2017")

# Change to same class for binding
ebird.tab$`Fall 2016` <- as.character(ebird.tab$`Fall 2016`)
ebird.tab$`Fall 2017` <- as.character(ebird.tab$`Fall 2017`)

# Ungroup ebird.tab
ebird.tab <- ungroup(ebird.tab)

# Add rows for totals
ebird.tab <- rbind(ebird.tab, totals.chk)
ebird.tab <- rbind(ebird.tab, totals.obs)

# Change numbers back to numeric class
ebird.tab$`Fall 2016` <- as.numeric(ebird.tab$`Fall 2016`)
ebird.tab$`Fall 2017` <- as.numeric(ebird.tab$`Fall 2017`)

# Save for exporting as a table to Word
write.csv(ebird.tab, 
          file = paste0(my.wd, "/Tables/Table1-summary-of-ebird-observations-for-extensive-sediment"), #Exclude .csv at end
          row.names = F, quote = F) 

#' 1. Save table object with write_csv() but DO NOT put the '.csv' format at the end (prevents formatting for Excel)
#' 2. Open with TextEdit
#' 3. Copy all the text
#' 4. Open word and paste text
#' 5. Highlight all the text in Word and from top toolbar, select 'table' > 'convert' > 'convert text to table'
#' 6. Ensure the number of columns is correct and that 'comma' is selected as the appropriate delim


