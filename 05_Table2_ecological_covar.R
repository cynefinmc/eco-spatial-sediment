####################################################################################
###                                                                             ###
###      ECOLOGICAL REGRESSION MODEL OF EXTENSIVE SEDIMENT MATRIX RESULTS       ###
###                                                                             ###
###      Table 2: Summarize ecological covariates used in regression model      ###
###                                                                             ###
###################################################################################


### 1. Set up working environment -------------------------------------------------
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

my.wd <- getwd()
source(paste0(my.wd, "/00_aiv_ext_functions.R"))

ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

#' Shorten some covariate names for easier coding
ext.sed <- rename(ext.sed,
                  "sal" = sal_ppt,       #' salinity
                  "temp" = temp_C,       #' temperature
                  "fdepth" = fdepth_cm)  #' sediment depth as a factor



### 2) Apply summary table functions to determine summary values for each covariate -----------------------------------------------------
### i) pH ----
pH.summary <- mSumTabCont(Dataset = ext.sed, 
                          Variable = "pH", 
                          PrettyVarName = "pH",
                          Ndigits = 2)

### ii) salinity ----
sal.summary <- mSumTabCont(Dataset = ext.sed,
                           Variable = "sal",
                           PrettyVarName = "Salinity (ppt)",
                           Ndigits = 2)


### iii) temperature ----
temp.summary <- mSumTabCont(Dataset = ext.sed,
                            Variable = "temp",
                            PrettyVarName = "Temperature (C)",
                            Ndigits = 2)


### iv) sediment characteristics ----
sedclass.summary <- mSumTabCat(Dataset = ext.sed,
                               Variable = "sedclass",
                               PrettyVarName = "Sediment quality",
                               Ndigits = 2)

### v) water depth ----
fdepth.summary <- mSumTabCat(Dataset = ext.sed,
                             Variable = "fdepth",
                             PrettyVarName = "Water depth (cm)",
                             Ndigits = 2)



### vi) feces -----
feces.summary <- mSumTabCat(Dataset = ext.sed,
                            Variable = "feces",
                            PrettyVarName = "Feces present",
                            Ndigits = 2)


### vii) wetland ----
wtld.summary <- mSumTabCat(Dataset = ext.sed,
                           Variable = "wtld",
                           PrettyVarName = "Wetland",
                           Ndigits = 2)

wtld.summary$category <- str_to_upper(wtld.summary$category)


### viii) spatial coordinates ----
east.summary <- mSumTabCont(Dataset = ext.sed,
                            Variable = "east",
                            PrettyVarName = "Easting (m)",
                            Ndigits = 0)

north.summary <- mSumTabCont(Dataset = ext.sed,
                             Variable = "north",
                             PrettyVarName = "Northing (m)",
                             Ndigits = 0)

### 4. Combine into final table -----------------------------------------------------------------------------
tab2 <- bind_rows(pH.summary,
                  sal.summary,
                  temp.summary,
                  fdepth.summary,
                  sedclass.summary,
                  feces.summary,
                  wtld.summary,
                  east.summary,
                  north.summary)

names(tab2) <- c("Covariate", "Category", "Total (%)", "AIV positive (%)", "AIV negative (%)")


# Clean up workspace
rm(pH.summary, sal.summary, temp.summary, sedclass.summary, fdepth.summary, feces.summary,
   wtld.summary, east.summary, north.summary)


### 5. Save as .csv for export to Word ----------------------------------------------------------------------
write.csv(tab2, 
          file = paste0(my.wd, "/Tables/Table2-summary-of-ecological-covariates"), #Exclude .csv at end
          row.names = F, quote = F) 

#' 1. Save table object with write_csv() but DO NOT put the '.csv' format at the end (prevents formatting for Excel)
#' 2. Open with TextEdit
#' 3. Copy all the text
#' 4. Open word and paste text
#' 5. Highlight all the text in Word and from top toolbar, select 'table' > 'convert' > 'convert text to table'
#' 6. Ensure the number of columns is correct and that 'comma' is selected as the appropriate delim

