#################################################################################################################################
###                                                                                                                           ###
###                             ECOLOGICAL REGRESSION MODEL OF EXTENSIVE SEDIMENT MATRIX RESULTS                              ###
###                                                                                                                           ###
###                             This script holds the FUNCTIONS and OBJECTS used for the analysis                             ###
###                                                                                                                           ###
#################################################################################################################################

#' The majority of these function were written by or heavily based on code written by Highland Statistics Ltd.
#' Code comes from Volume I and II of Beginner's Guide to Spatial, Temporal and Spatial-Temporal Ecological
#'   Data Analysis with R-INLA by Zuur, Ieno and Saveliev, and associated courses.
#' Books and courses can be found at: https://www.highstat.com


### 1.Identification keys   --------------------------------------------------------------------------------------------------------------------------

# Extensive sites with abbreviations
ext.sites <- c("BHER", "BRY", "BURN", "LMT",
               "MIL", "MNP", "NIC", "PIT", 
               "REI", "SAR", "SDW", "SFEN", 
               "SIS", "STA", "TRL", "TUY",
               "TWU", "WCK")

### Species key that lines up the scientific name, common name, and abbreviation
#Exclude hybrids they aren't listed in the filtered ebird object
species.key <- data.frame(common_name = c("American Wigeon", 
                                          "Barrow's Goldeneye", 
                                          "Blue-winged Teal", 
                                          "Brant",
                                          "Bufflehead", 
                                          "Cackling Goose", 
                                          "Canada Goose", 
                                          "Canvasback", 
                                          "Cinnamon Teal", 
                                          "Common Goldeneye", 
                                          "Common Merganser", 
                                          "Eurasian Wigeon", 
                                          "Gadwall", 
                                          "Greater Scaup", 
                                          "Greater White-fronted Goose", 
                                          "Green-winged Teal",  
                                          "Hooded Merganser",
                                          "Lesser Scaup", 
                                          "Mallard", 
                                          "Mute Swan", 
                                          "Northern Pintail", 
                                          "Northern Shoveler", 
                                          "Red-breasted Merganser", 
                                          "Redhead",  
                                          "Ring-necked Duck", 
                                          "Ruddy Duck", 
                                          "Snow Goose",  
                                          "Surf Scoter", 
                                          "Trumpeter Swan", 
                                          "Tundra Swan", 
                                          "White-winged Scoter", 
                                          "Wood Duck"),
                          spp = c("amewig",
                                  "bargol",
                                  "buwtea",
                                  "brant",
                                  "buffle",
                                  "cackgoo",
                                  "cangoo",
                                  "canvas",
                                  "cintea",
                                  "comgol",
                                  "commer",
                                  "eurwig",
                                  "gadwall",
                                  "gresca",
                                  "gwfgoo",
                                  "gnwtea",
                                  "hoomer",
                                  "lessca",
                                  "mallard",
                                  "mutswa",
                                  "norpin",
                                  "norsho",
                                  "rebmer",
                                  "redhead",
                                  "rinduc",
                                  "rudduc",
                                  "snogoo",  
                                  "sursco",
                                  "truswa",
                                  "tunswa", 
                                  "whwsco",
                                  "wooduc"),
                          scientific_name = c("Mareca americana",
                                              "Bucephala islandica",
                                              "Spatula discors",
                                              "Branta bernicla",
                                              "Bucephala albeola",
                                              "Branta hutchinsii",
                                              "Branta canadensis",
                                              "Aythya valisineria",
                                              "Spatula cyanoptera",
                                              "Bucephala clangula",
                                              "Mergus merganser", 
                                              "Mareca penelope",
                                              "Mareca strepera",
                                              "Aythya marila",
                                              "Anser albifrons",
                                              "Anas crecca",
                                              "Lophodytes cucullatus",
                                              "Aythya affinis",
                                              "Anas platyrhynchos",
                                              "Cygnus olor",
                                              "Anas acuta",
                                              "Spatula clypeata",
                                              "Mergus serrator",
                                              "Aythya americana",
                                              "Aythya collaris",
                                              "Oxyura jamaicensis",
                                              "Anser caerulescens",
                                              "Melanitta perspicillata",
                                              "Cygnus buccinator",
                                              "Cygnus columbianus",
                                              "Melanitta deglandi",
                                              "Aix sponsa"),
                          stringsAsFactors = F)


### 3. Functions  ----------------------------------------------------------------------------------------------------------------------------

#######################################################################
#                     %nin% function                                  #
#       this function acts as the opposite to the function '%in%'     #
#######################################################################

'%nin%' <- Negate('%in%')


#############################################################################
#                     m.getresid.bern FUNCTION                              #
#                                                                           #
#    Extracts fitted values, calculates Pearson residuals and returns       #
#    a data frame with the Observed, fitted and Pearson residuals from an   #
#    INLA object. Currently only set up for BERNOULLI models                #
# Works with the following call, but not when multiple models in the INLA list
#  MyINLAlist <- list("M.glmm2" = M.glmm2)
#  Does NOT work directly with MyINLAlist = M.glmm2
#############################################################################

m.getresid.bern <- function(MyData, MyYvar, MyINLAlist){
  
  ResidList = list()
  
  for(i in seq_along(MyINLAlist)){
    N <- nrow(MyData)
    Pi <- purrr::map(MyINLAlist[i], list("summary.fitted.values", "mean")) %>% flatten_dbl() 
    Pi <- Pi[1:N]
    
    ExpY <- Pi  
    VarY <- (1 - Pi) * Pi 
    E1 <- (MyYvar - ExpY) / sqrt(VarY) #Pearson residuals
    
    ObsY <- MyYvar
    
    # Put Pearson residuals, Expected values, Observed values, model name AND covariates into a tibble
    df <- tibble("E1" = E1, 
                 "ExpY" = ExpY, 
                 "ObsY" = ObsY,
                 "model" = rep(names(MyINLAlist)[i], N)) 
    
    df <- cbind(df, MyData)
    df <- as_tibble(df)
    
    ResidList[[i]] <- df
    
    # Just for trouble-shooting a loop
    #print(names(MyINLAlist)[i])
  }
  resid.df <- bind_rows(ResidList)
}

#############################################################################
#                         m.getws FUNCTION                                  #
#                                                                           #
#    Extracts posterior means and standard deviations of w's from           #
#    multiple INLA objects with a SPATIAL field, and returns a data frame   #
#############################################################################

m.getws <- function(myspatialmodels){
  
  wList = list()
  
  for(i in seq_along(myspatialmodels)){
    w.pm <- pluck(myspatialmodels[i], 1, "summary.random", "w", "mean")   #Posterior means of w's
    w.sd <- pluck(myspatialmodels[i], 1, "summary.random", "w", "sd")     #Standard deviations of w's
    
    # Put posterior means model name into a tibble
    df <- tibble("w.pm" = w.pm, 
                 "w.sd" = w.sd, 
                 "model" = rep(names(myspatialmodels)[i], length(w.pm))) 
    
    wList[[i]] <- df
    names(wList)[i] <- names(myspatialmodels)[i]
    
    # Just for trouble-shooting a loop
    #print(names(wList)[i])
  }
  w.df <- bind_rows(wList)
}

#######################################################################################
#                     m.gg2PlotField FUNCTION                                         #
#   Allows you to plot the spatial mean and sd generated in INLA                      #
#                                                                                     #
#   This is a ggplot version of the PlotFields function written by Hakkon Baaka       #
#   Written by Highland Statistics Ltd and I've slightly modified.                    #
#                                                                                     #
#   field = mean or sd of random spatial field                                        #
#         eg) mymodel$summary.random$w$mean or mymodel$summary.random$w$sd            #
#   mesh = inla mesh object used within your model                                    #
#         eg) mesh1                                                                   #
#   scalemin and scalemax = numeric values for min and max of colorgradient legend    #
#         to ensure that multiple graphs with different values of w.pm are comparable #
#   xlim and ylim coords OR ContourMap (something that has a bbox) to                 #
#         provide the x and y axes for the graph                                      #
#   You can also change the color scheme or supply x, y, and title labels             #
#######################################################################################

m.gg2PlotField <- function(field, 
                           mesh, 
                           ContourMap, 
                           xlim, 
                           ylim, 
                           MyTitle = "",
                           ColourScheme = 1, 
                           Xlabel = "", 
                           Ylabel = "", 
                           LegendName = "",
                           TextSize = NULL,
                           AddPoints = FALSE, #if want to add points to the plot
                           PointSize = 0.1,
                           MyData = "",  #need specify a dataset for points if AddPoints = T
                           scalemin,        #scalemin and scalemax are numeric values so can compare between graphs
                           scalemax, ...){
  stopifnot(length(field) == mesh$n)
  # Plotting region to be the same as the study area polygon
  if (missing(xlim)) xlim <- ContourMap@bbox[1, ] 
  if (missing(ylim)) ylim <- ContourMap@bbox[2, ]
  
  # inla.mesh.projector creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, xlim = xlim, ylim = ylim,
                              dims = c(300, 300))
  
  # The function inla.mesh.project can then be used to project the w's on this grid.
  field.proj <- inla.mesh.project(proj, field)
  
  # Make data frame from grid
  df <- expand.grid(x = proj$x, y = proj$y)
  df$w.pm <- as.vector(field.proj)
  
  # Determine limits for color gradient bar
  if (missing(scalemin)) scalemin <- min(df$w.pm) 
  if (missing(scalemax)) scalemax <- max(df$w.pm)
  
  # Make plot
  w <- ggplot()
  w <- w + geom_raster(data = df, aes(x = x, y = y, fill = w.pm))
  w <- w + coord_fixed(ratio = 1) + theme_bw() + ggtitle(MyTitle)
  w <- w + xlab(Xlabel) + ylab(Ylabel)
  
  if(ColourScheme == 1){
    w <- w + scale_fill_gradient2(name = LegendName,
                                  midpoint = 0, 
                                  low = "blue", 
                                  mid = "white", 
                                  high = "red",
                                  na.value = "transparent",
                                  limits = c(scalemin, scalemax))}
  
  if(ColourScheme == 2){
    w <- w + scale_fill_gradient(name = LegendName,
                                 low = "blue", 
                                 high = "orange",
                                 na.value = "transparent",
                                 limits = c(scalemin, scalemax))}
  if(ColourScheme == 3){
    w <- w + scale_fill_viridis(name = LegendName,
                                na.value = "transparent",
                                option = "D",
                                limits = c(scalemin, scalemax))}
  
  if(ColourScheme == 4){
    w <- w + scale_fill_viridis(name = LegendName,
                                na.value = "transparent",
                                option = "G",
                                limits = c(scalemin, scalemax))}
  
  if(ColourScheme == 5){
    w <- w + scale_fill_viridis(name = LegendName,
                                na.value = "transparent",
                                option = "E",
                                limits = c(scalemin, scalemax))}
  if(ColourScheme == 6){
    w <- w + scale_fill_gradientn(name = LegendName,
                                  colours = rainbow(5),
                                  na.value = "transparent",
                                  limits = c(scalemin, scalemax))}
  
  if (AddPoints){
    if(!is.data.frame(MyData)){
      stop("If AddPoints = TRUE, you need to supply a data frame for MyData", call. = F)
    }
    w <- w + geom_point(data = MyData,
                        aes(x = east, y = north),
                        size = PointSize) +
      #Need to reset coords limits for extent of plot again
      xlim(xlim) +  
      ylim(ylim)}
  
  if(!is.null(TextSize)){
    w <- w + theme(text = element_text(size = TextSize))
  }
  
  w
}


##################################################################################################
#                             m.GetPredictedpH function                                          #
#                                                                                                #
#   This function predicts pH values for different wetland or site locations, with and without   #
#   the presences of feces covariate, using the final INLA model                                 #
#   There are two ways to run this function:                                                     #
#       1. With an example site: use LocType = "stID" and providing an stID                      #
#       2. Using the mean location for a specific wetland: use LocType = "mean" and provide wtld #
#                                                                                                #
#                                                                                                #
#   LocType = character string indicating how function should obtain the location info           #
#             can be either 'stID' or 'mean'                                                     #       
#   stID = default is NULL; otherwise if LocType is "stID", need character string in the         #
#          format of stID column (e.g., "BURN-S14_DEC-2016")                                     #
#   MyPredWtld = default is NULL; otherwise if LocType is "mean", need to provide a wtld that    #
#                exactly matches one of the 18 extensive wetlands (e.g., BURN, MIL)              #
#                                                                                                #
##################################################################################################

#' To run, need to have loaded INLA objects used in making model
#'  eg) mesh.sp, A.sp, spde.sp, w.index.sp, Stack1, Lego.pH

m.GetPredictedpH <- function(LocType, MystID = NULL, MyPredWtld = NULL){
  
  ### i) Obtain columns with location information ----
  if(LocType == "stID"){
    LocDF <- ext.sed %>%
      filter(stID == MystID) %>%
      select(wtld, site, stID, east, north, longitude, latitude)
    
    print(paste("Starting to generate predicted data for:", unique(LocDF$stID)))
  } 
  
  if(LocType == "mean"){
    LocDF <- ext.sed %>% 
      filter(wtld %in% MyPredWtld) %>% 
      group_by(wtld) %>% 
      summarize(wtld = unique(wtld),
                east = mean(east),
                north = mean(north),
                longitude = mean(longitude),
                latitude = mean(latitude)) %>% 
      mutate(site = "avg_site",
             stID = paste0(MyPredWtld, "_avg-site")) %>% 
      select(wtld, site, stID, east, north, longitude, latitude) #So columns are in the same order
    
    print(paste("Starting to generate predicted data for:", "the mean location within", unique(LocDF$wtld)))
  }
  
  if(nrow(LocDF) != 1){
    print("Warning: Location coordinates are not equal to 1. Check issue manually.")
  }
  
  
  ### ii) Create dataframe of artifical covariate values ----
  ### Feces = Yes
  #' should only predict within pH range that was observed for presence of feces
  MyDataYes <- data.frame(pH = seq(from = min(ext.sed$pH[ext.sed$feces == "yes"]), 
                                   to = max(ext.sed$pH[ext.sed$feces == "yes"]),
                                   length = 100),
                          sal = mean(ext.sed$sal),
                          temp = mean(ext.sed$temp),
                          feces = factor("yes", levels = c("no", "yes")),
                          fdepth = factor("1 -25", levels = c("1 -25", "over 25")),
                          sedclass = factor("floats", levels = c("floats", "mud", "grit", "rock")))
  
  MyDataYes$east <- unique(LocDF$east)
  MyDataYes$north <- unique(LocDF$north)
  MyDataYes$wtld <- unique(LocDF$wtld)
  MyDataYes$stID <- unique(LocDF$stID)
  
  # Set Y variable (matrix01) to NA to predict for this variable
  MyDataYes$matrix01 <- NA_character_
  
  # Determine Nrow for data
  Nyes <- nrow(MyDataYes) #100
  
  ### Feces = No
  #' should only predict within pH range that was observed for absence of feces 
  MyDataNo <- data.frame(pH = seq(from = min(ext.sed$pH[ext.sed$feces == "no"]), 
                                  to = max(ext.sed$pH[ext.sed$feces == "yes"]),
                                  length = 100),
                         sal = mean(ext.sed$sal),
                         temp = mean(ext.sed$temp),
                         feces = factor("no", levels = c("no", "yes")),
                         fdepth = factor("1 -25", levels = c("1 -25", "over 25")),
                         sedclass = factor("floats", levels = c("floats", "mud", "grit", "rock")))
  
  MyDataNo$east <- unique(LocDF$east)
  MyDataNo$north <- unique(LocDF$north)
  MyDataNo$wtld <- unique(LocDF$wtld)
  MyDataNo$stID <- unique(LocDF$stID)
  
  # Set Y variable (matrix01) to NA to predict for this variable
  MyDataNo$matrix01 <- NA_character_
  
  # Determine Nrow for data
  Nno <- nrow(MyDataNo) #100
  
  
  ### iii) Create model matrix of covariates using predicted data ----
  #' Use the model.matrix() function for INLA models using spatial correlation
  
  ### Feces = Yes
  Xyes <- model.matrix(~ sal + temp + feces + fdepth + sedclass,  #do NOT include pH as it is being modeled as a smoother
                       data = MyDataYes)
  Xyes <- as.data.frame(Xyes)
  
  colnames(Xyes) 
  OldColnames  <- colnames(Xyes)
  NewColNames  <- gsub(pattern = "fdepthover 25", replacement = "fdepth_over25", x = OldColnames)
  NewColNames  <- gsub(pattern = "fecesyes", replacement = "feces_yes", x = NewColNames)
  NewColNames  <- gsub(pattern = "sedclass", replacement = "sedclass_", x = NewColNames)
  NewColNames
  colnames(Xyes) <- NewColNames
  
  
  ### Feces = No
  Xno <- model.matrix(~ sal + temp + feces + fdepth + sedclass,  #do NOT include pH as it is being modeled as a smoother
                      data = MyDataNo)
  Xno <- as.data.frame(Xno)
  
  if(identical(colnames(Xno), OldColnames) == F){
    print("WARNING: need to fix Xno column names manually")
  } 
  
  if(identical(colnames(Xno), OldColnames) == T){
    colnames(Xno) <- NewColNames
  } 
  
  ### iv) Create basis function for the predicted pH values ----
  #' `Lego.pH` is the basis function used with the observed data
  
  ### Feces = Yes
  XpBasisYes <- PredictMat(Lego.pH, MyDataYes)
  colnames(XpBasisYes) <- paste0("XpH", 1:ncol(XpBasisYes)) #Need colnames() not names() as is a matrix not df
  
  ncy <- ncol(XpBasisYes)
  nry <- ncol(XpBasisYes)
  
  ### Feces = No
  XpBasisNo <- PredictMat(Lego.pH, MyDataNo)
  colnames(XpBasisNo) <- paste0("XpH", 1:ncol(XpBasisNo)) #Need colnames() not names() as is a matrix not df
  
  ncn <- ncol(XpBasisNo)
  nrn <- ncol(XpBasisNo)
  
  
  ### v) Determine the spatial location and make A matrix for it ----
  #' Needs to be in *c(easting, northing)* format to match observed A matrix
  LocPred <- cbind(MyDataYes$east, MyDataYes$north) #Will be same locations for feces Yes or No so will just make one A
  
  A.Pred <- inla.spde.make.A(mesh = mesh.sp, loc = LocPred)
  
  
  ### vi) Create a stack for the data for which predictions are needed ----
  #' Note that NA is needed for Y variable to make INLA predict for that variable 
  
  ### Feces = Yes
  StackPredYes <- inla.stack(
    tag = "Predict_yes",
    data = list(matrix01 = NA),  # <- need to use NA for response variable!
    A = list(1, 1, 1, A.Pred),     # <- need to use Predicted A matrix!
    effects = list(
      Intercept = rep(1, Nyes),  # <- need to use N row in predicted data (1 in this case)
      XpHYes = as.data.frame(XpBasisYes),  #Need to make sure this the data frame object with column names as in the formula call
      Xyes = Xyes[ ,-1],
      w = w.index.sp))
  
  
  ### Feces = No
  StackPredNo <- inla.stack(
    tag = "Predict_no",
    data = list(matrix01 = NA),  # <- need to use NA for response variable!
    A = list(1, 1, 1, A.Pred),     # <- need to use Predicted A matrix!
    effects = list(
      Intercept = rep(1, Nno),  # <- need to use N row in predicted data (1 in this case)
      XpHNo = as.data.frame(XpBasisNo),  #Need to make sure this the data frame object with column names as in the formula call
      Xno = Xno[ ,-1],
      w = w.index.sp))
  
  
  ### vii) Combine Stacks ----
  All.stacks <- inla.stack(Stack1, StackPredYes, StackPredNo)
  
  
  ### viii) Run INLA ----
  print(paste("Running predictions in INLA for:", unique(LocDF$stID)))
  
  M.Pred <- inla(matrix01 ~ -1 + Intercept + XpH1 + XpH2 + feces_yes + 
                   f(w, model = spde.sp),
                 family = "binomial",
                 Ntrials = 1,  
                 data = inla.stack.data(All.stacks),                  # <- Use Combined Stack here
                 control.predictor = list(link = 1,                   # <- results will be on scale of response variable
                                          A = inla.stack.A(All.stacks))) # <- Use Combined Stack here
  
  
  ### ix) Extract the relevant objects ----
  
  # Determine the indices for the predicted values 
  index.Yes <- inla.stack.index(All.stacks,
                                tag = "Predict_yes")$data
  
  index.No <- inla.stack.index(All.stacks,
                               tag = "Predict_no")$data
  
  # Extract the predicted matrix01 values
  PredYesValues <- M.Pred$summary.fitted.values[index.Yes, c("mean", "0.025quant", "0.975quant")]
  PredNoValues <- M.Pred$summary.fitted.values[index.No, c("mean", "0.025quant", "0.975quant")]
  
  # Add them to the PredData object
  MyYesData <- cbind(MyDataYes, PredYesValues)
  MyNoData <- cbind(MyDataNo, PredNoValues)
  
  MyPredData <- bind_rows(MyYesData, MyNoData)
  
  MyPredData <- rename(MyPredData,
                       mu = "mean",
                       selow = "0.025quant",
                       seup = "0.975quant")
  
  if(LocType == "stID"){
    print(paste("Finished adding predicted matrix01 value to data for:", 
                unique(LocDF$stID)))
  }
  
  if(LocType == "mean"){
    print(paste("Finished adding predicted matrix01 value to data for:", 
                "the mean location within", unique(LocDF$wtld)))
  }
  
  # Return this object to the Global environment and exit function
  return(MyPredData)
  
} # End of m.GetPredictedpH function



##################################################################################################
#                                 m.MakeBbox function                                            #
#                                                                                                #
#   This function makes a bounding box from a dataframe with coordinate columns                  #
#   It allows you to add some wiggle room around the provided coordinates                        #
#   The output can then be used to make a basemap, or otherwise define plotting parameters       #
#                                                                                                #
##################################################################################################

m.MakeBbox <- function(Dataset, 
                       CoordNames = c("X", "Y"),
                       Buffer = NULL,
                       LocCol = NULL,
                       Location = NULL,
                       MakeSpatial = T,
                       MyCRS = NULL){
  
  require(dplyr)
  require(sf) # To turn bbox into a spatial object
  
  # Filter dataset to location of interest and ensure there are standard coordinate column names
  if(!is.null(Location)){
    loc.i    <- which(names(Dataset) == LocCol)
    xcoord.i <- which(names(Dataset) == CoordNames[1])
    ycoord.i <- which(names(Dataset) == CoordNames[2])
    
    names(Dataset)[loc.i]    <- "location"
    names(Dataset)[xcoord.i] <- "X"
    names(Dataset)[ycoord.i] <- "Y"
    
    mydf <- Dataset[ , c(loc.i, xcoord.i, ycoord.i)]
    mydf <- mydf[mydf$location == Location, ]
    mydf$location <- droplevels(mydf$location)
    
  } else {
    xcoord.i <- which(names(Dataset) == CoordNames[1])
    ycoord.i <- which(names(Dataset) == CoordNames[2])
    
    names(Dataset)[xcoord.i] <- "X"
    names(Dataset)[ycoord.i] <- "Y"
    
    mydf <- Dataset[, c(xcoord.i, ycoord.i)]
  }
  
  # Calculate the min and max values of coordinates, with a buffer if desired
  if(!is.null(Buffer)){
    xadjust <- (max(mydf$X) -min(mydf$X)) * Buffer
    yadjust <- (max(mydf$Y) -min(mydf$Y)) * Buffer
    
    minmax <- mydf %>% 
      summarize(min_x = min(X) - xadjust,
                min_y = min(Y) - yadjust,
                max_x = max(X) + xadjust,
                max_y = max(Y) + yadjust)
  } else {
    minmax <- mydf %>% 
      summarize(min_x = min(X),
                min_y = min(Y),
                max_x = max(X),
                max_y = max(Y))
  }
  
  # Turn into a data frame
  bbox <- data.frame(X = c(minmax$min_x, minmax$max_x),
                     Y = c(minmax$min_y, minmax$max_y))
  rownames(bbox) <- c("min", "max")
  
  # Make into a spatial object recognized by `sf` package
  if(MakeSpatial){
    if(is.null(MyCRS)) stop("MyCRS argument needs to be specified to return an sf boundary box")
    
    bbox <- st_as_sf(x = bbox,
                     coords = c("X", "Y"),
                     crs = "+proj=longlat +datum=WGS84")
  }
  
  return(bbox)
}


##################################################################################################
#                                 mSumTabCont function                                           #
#                                                                                                #
#   This function generates summary values of CONTINOUS variables in ext.sed dataframe           #
#   so they can be combined into a table suitable for publication                                #
#   Summary outputs are median, min and max values                                               #
#                                                                                                #
##################################################################################################

mSumTabCont <- function(Dataset, Variable, PrettyVarName, Ndigits){
  DF <- as_tibble(Dataset[ , Variable])
  names(DF) <- "variable"
  
  DF.pos <- filter(Dataset, matrix01 == 1)
  DF.pos <- as_tibble(DF.pos[ , Variable])
  names(DF.pos) <- "variable"
  
  DF.neg <- filter(Dataset, matrix01 == 0)
  DF.neg <- as_tibble(DF.neg[ , Variable])
  names(DF.neg) <- "variable"
  
  sum.tab <- DF %>% 
    summarize(median = median(variable),
              min = min(variable),
              max = max(variable)) %>% 
    mutate("var" = PrettyVarName,
           "category" = "Median (Range)",
           "sum_stat" = paste0(format(round(median, digits = Ndigits), nsmall = Ndigits),
                               " (",
                               format(round(min, digits = Ndigits), nsmall = Ndigits),
                               " - ",
                               format(round(max, digits = Ndigits), nsmall = Ndigits),
                               ")")) %>% 
    select(var, category, sum_stat)
  
  sum.pos <- DF.pos %>% 
    summarize(median = median(variable),
              min = min(variable),
              max = max(variable)) %>% 
    mutate("var" = PrettyVarName,
           "category" = "Median (Range)",
           "AIV_pos" = paste0(format(round(median, digits = Ndigits), nsmall = Ndigits),
                              " (",
                              format(round(min, digits = Ndigits), nsmall = Ndigits),
                              " - ",
                              format(round(max, digits = Ndigits), nsmall = Ndigits),
                              ")")) %>% 
    select(var, category, AIV_pos)
  
  sum.neg <- DF.neg %>% 
    summarize(median = median(variable),
              min = min(variable),
              max = max(variable)) %>% 
    mutate("var" = PrettyVarName,
           "category" = "Median (Range)",
           "AIV_neg" = paste0(format(round(median, digits = Ndigits), nsmall = Ndigits),
                              " (",
                              format(round(min, digits = Ndigits), nsmall = Ndigits),
                              " - ",
                              format(round(max, digits = Ndigits), nsmall = Ndigits),
                              ")")) %>% 
    select(var, category, AIV_neg)
  
  sum.tab <- full_join(sum.tab, sum.pos, by = c("var", "category"))
  sum.tab <- full_join(sum.tab, sum.neg, by = c("var", "category"))
}



##################################################################################################
#                                 mSumTabCat function                                            #
#                                                                                                #
#   This function generates summary values of CATEGORICAL variables in ext.sed dataframe         #
#   Where the Y variable is positive (1) or negative (0)                                         #
#   so they can be combined into a table suitable for publication                                #
#   Summary outputs are Number and Percent of all, positive and negative observations            #
#                                                                                                #
##################################################################################################

mSumTabCat <- function(Dataset, Variable, PrettyVarName, Ndigits){
  DF <- as_tibble(Dataset[ , c(Variable, "matrix01")])
  names(DF) <- c("variable", "matrix01")
  
  sum.tab <- DF %>% 
    group_by(variable) %>% 
    summarize(n_samples = n(),
              perc_samples = n()/nrow(DF) *100,
              n_pos = sum(matrix01 == 1),
              perc_pos = n_pos / n() *100,
              n_neg = sum(matrix01 == 0),
              perc_neg = n_neg / n() *100) %>% 
    mutate("var" = PrettyVarName,
           "category" = variable,
           "sum_stat" = paste0(n_samples,
                               " (",
                               format(round(perc_samples, digits = Ndigits), nsmall = Ndigits),
                               ")"),
           "AIV_pos" = paste0(n_pos,
                              " (",
                              format(round(perc_pos, digits = Ndigits), nsmall = Ndigits),
                              ")"),
           "AIV_neg" = paste0(n_neg,
                              " (",
                              format(round(perc_neg, digits = Ndigits), nsmall = Ndigits),
                              ")")) %>% 
    select("var", "category", "sum_stat", "AIV_pos", "AIV_neg")
  
  sum.tab$category <- as.character(sum.tab$category)
  sum.tab$category <- str_to_title(sum.tab$category)
  
  return(sum.tab)
  
}


















