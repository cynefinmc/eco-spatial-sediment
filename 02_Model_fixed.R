###########################################################################################################
###                                                                                                     ###
###                ECOLOGICAL REGRESSION MODEL OF EXTENSIVE SEDIMENT MATRIX RESULTS                     ###
###                                                                                                     ###                                                                                                     ###
###          This script determines the most appropriate FIXED structure for the model                  ###
###                     using backwards and forwards model selection                                    ###                  
###                                                                                                     ###
###########################################################################################################

#' Much of this code is heavily based on code written by Highland Statistics Ltd.
#' Code comes from Volume I and II of Beginner's Guide to Spatial, Temporal and Spatial-Temporal Ecological
#'   Data Analysis with R-INLA by Zuur, Ieno and Saveliev, and associated courses.
#' Books and courses can be found at: https://www.highstat.com


### (1) Load packages and data   --------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(scales) 
library(dplyr)
library(tibble)
library(lubridate)
library(purrr)
library(lme4)      
library(gamm4)     
library(INLA)
library(brinla)
library(fields) 

#Set seed so have reproducible results
set.seed(858)

my.wd <- getwd()
source(paste0(my.wd, "/00_aiv_ext_functions.R"))

ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

#' Shorten some covariate names for easier coding
ext.sed <- rename(ext.sed,
                  "sal" = sal_ppt,       #' salinity
                  "temp" = temp_C,       #' temperature
                  "fdepth" = fdepth_cm)  #' sediment depth as a factor

#' Remove the ebird summary covariates, as at the level of wetland and have missing data, 
#' and MODIS non-water dominant landcover due to collinearity
ext.sed <- select(ext.sed,
                  -c("Matrix_positivity", "N_checklist", "Total_abundance", "Relative_abundance", "Species_richness",
                     "dominant_MODIS_landcover_class"))


### (2) Set up objects for the 'full' model  --------------------------------------------------------------------------------------------------------------------------------------

# i) Make smoother for the pH covariate  ----
#The smoother is defined as f(pH) = XpH * beta and does NOT include the intercept

# Get basis function using smoothCon
  #This is a unpenalized (fx = T) cubic regression spline (bs = "cr") smoother
  #k is dim of basis to represent the smooth term, see ?s
Lego.pH <- smoothCon(s(pH, bs = "cr", k = 3, fx = TRUE), #choosing 3 knots as expect 2 edf based on gamm4 model
                     data = ext.sed,  
                     absorb.cons = TRUE)[[1]]
XpH <- Lego.pH$X 

# Turn into data frame
XpH.df <- data.frame(XpH)
names(XpH.df) <- paste0("XpH", 1:ncol(XpH))
head(XpH.df)

# Need the lincomb step to be able to plot the smoother
#This code calculates the sum of the basis functions: XpH * beta = f(pH)
lcs.pH <- inla.make.lincombs(XpH.df)


# ii) Make an 'X' matrix for linear covariates  ----

# Xm 
Xm <- model.matrix(~ sal + temp + feces + fdepth + sedclass,  #do NOT include pH as it is being modeled as a smoother
                   data = ext.sed)
Xm <- as.data.frame(Xm)

colnames(Xm)
OldColnames  <- colnames(Xm)
NewColNames  <- gsub(pattern = "fdepthover 25", replacement = "fdepth_over25", x = OldColnames)
NewColNames  <- gsub(pattern = "fecesyes", replacement = "feces_yes", x = NewColNames)
NewColNames  <- gsub(pattern = "sedclass", replacement = "sedclass_", x = NewColNames)
NewColNames
colnames(Xm) <- NewColNames



# iii) Make meshes for spatial models  ----
RangeGuess <- 100 #in meters, estimate based on observations of waterfowl feeding in wetlands while in the field
MaxEdge <- RangeGuess / 5 
ELoc <- cbind(ext.sed$east, ext.sed$north) 
OuterBound <- 1500
bound <- inla.nonconvex.hull(ELoc, convex = -0.003, 
                             resolution = c(340, 115)) #Increased resolution as per INLA warning


mesh.sp <- inla.mesh.2d(ELoc,
                       boundary = bound,
                       max.edge = c(2* 1, 100 * 5) * MaxEdge,
                       cutoff = 2* MaxEdge / 5,
                       offset = c(MaxEdge, OuterBound))
plot(mesh.sp, asp = 1)
points(ELoc, col = 2, pch = 16, cex = 1)
mesh.sp$n  #n = 21,958
#save(mesh.sp, file = "Results/Final-inla-objects/Extsed_meshsp_obj_fromdesktop_n21958.rdata")

# iv) Make projector matrices for spatial model  ----
A.sp <- inla.spde.make.A(mesh.sp, loc = ELoc)
#save(A.sp, file = "Results/Final-inla-objects/Extsed_Asp_obj_fromdesktop_n21958.rdata")


# v) Define SPDE and spatial field (w.index) for spatial models  ----
    #see notes in 'Mesh_building.R' for how chose priors
spde.sp <- inla.spde2.pcmatern(mesh.sp,
                              prior.range = c(100, 0.05),
                              prior.sigma = c(2, 0.05))

w.index.sp <- inla.spde.make.index(name = "w",
                                  n.spde = spde.sp$n.spde)
#save(spde.sp, file = "Results/Final-inla-objects/Extsed_SPDEsp_obj_fromdesktop_n21958.rdata")
#save(w.index.sp, file = "Results/Final-inla-objects/Extsed_windexsp_obj_fromdesktop_n21958.rdata")


# vii) Define the number of observations ----
N <- nrow(ext.sed)



### (3) Run the 'full' model: Bernoulli GAM  + spatial field (Msp.full = "spatial_full model")  --------------------------------------------------------------------------------------------------------------------------------------
  #Using the spatial random field to account for small scale variation
mesh.sp$n                       #Done; 21958
dim(A.sp)                       #Done; 701 21958
str(spde.sp, max.level = 1)     #Already done
str(w.index.sp, max.level = 1)  #Already done
head(Xm)                       #Already done
head(XpH)                      #Already done
head(XpH.df)                   #Already done
N                              #Already done

Stack1 <- inla.stack(tag = "Fit_full",
                     data = list(matrix01 = ext.sed$matrix01), 
                     A = list(1, 1, 1, A.sp),
                     effects = list(
                       Intercept = rep(1, N), 
                       XpH = XpH.df,  #Need to make sure this the data frame object with column names as in the formula call
                       Xm = Xm[ ,-1], 
                       w = w.index.sp)) 
#save(Stack1, file = "Results/Final-inla-objects/Extsed_Stack1_obj_fromdesktop_n21958.rdata")

f.full <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                 XpH1 + XpH2 + #Need to manually enter the column names from XpH.df 
                 sal + temp +     
                 feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                 f(w, model = spde.sp))     #Spatial Random Field

Msp.full <- inla(f.full,
              family = "binomial",
              Ntrials = 1,  
              lincomb = lcs.pH,
              data = inla.stack.data(Stack1),
              control.inla = list(strategy = 'simplified.laplace'), #default option
              control.compute = list(dic = TRUE, waic = TRUE),
              control.predictor = list(A = inla.stack.A(Stack1)))

summary(Msp.full)
print(Msp.full$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(Msp.full)
Msp.full$dic$dic
Msp.full$waic$waic
#save(Msp.full, file = "Results/Final-inla-objects/Extsed_gam_spatial_meshsp_n21958_FULL.rdata")


### (4) Backwards Model Selection  --------------------------------------------------------------------------------------------------------------------------------------

### i) Round 1 ----
# Remove sal
f.01 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                     XpH1 + XpH2 + 
                     #sal + 
                     temp + feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                     f(w, model = spde.sp))     

M.01 <- inla(f.01,
                 family = "binomial",
                 Ntrials = 1,  
                 lincomb = lcs.pH,
                 data = inla.stack.data(Stack1),
                 control.inla = list(strategy = 'simplified.laplace'), #default option
                 control.compute = list(dic = TRUE, waic = TRUE),
                 control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.01)
print(M.01$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.01)
M.01$dic$dic
M.01$waic$waic

# Remove temp
f.02 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + 
                   #temp + 
                   feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                   f(w, model = spde.sp))     

M.02 <- inla(f.02,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.02)
print(M.02$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.02)
M.02$dic$dic
M.02$waic$waic

# Remove feces
f.03 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + 
                   #feces_yes + 
                   fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                   f(w, model = spde.sp))     

M.03 <- inla(f.03,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.03)
print(M.03$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.03)
M.03$dic$dic
M.03$waic$waic


# Remove fdepth
f.04 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + feces_yes + 
                   #fdepth_over25 + 
                   sedclass_mud + sedclass_grit + sedclass_rock +    
                   f(w, model = spde.sp))     

M.04 <- inla(f.04,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.04)
print(M.04$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.04)
M.04$dic$dic
M.04$waic$waic


# Remove sedclass
f.05 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + feces_yes + fdepth_over25 + 
                   #sedclass_mud + sedclass_grit + sedclass_rock +    
                   f(w, model = spde.sp))     

M.05 <- inla(f.05,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.05)
print(M.05$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.05)
M.05$dic$dic
M.05$waic$waic


# Remove pH
Stack.nopH <- inla.stack(tag = "Fit",
                     data = list(matrix01 = ext.sed$matrix01), 
                     A = list(1, 1, A.sp),
                     effects = list(
                       Intercept = rep(1, N), 
                       #XpH = XpH.df,  
                       Xm = Xm[ ,-1], 
                       w = w.index.sp))

f.06 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   #XpH1 + XpH2 + 
                   sal + temp + feces_yes + fdepth_over25 + 
                   sedclass_mud + sedclass_grit + sedclass_rock +    
                   f(w, model = spde.sp))     

M.06 <- inla(f.06,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Change the Stack
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #change the stack

summary(M.06)
print(M.06$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.06)
M.06$dic$dic
M.06$waic$waic


# Round 1 model selection ----
Rnd1.Models <- list("All covar" = Msp.full,
                    "no_sal" = M.01, 
                    "no_temp" = M.02, 
                    "no_feces" = M.03, 
                    "no_fdepth" = M.04, 
                    "no_sedclass" = M.05, 
                    "no_pH" = M.06)

Rnd1.Names <- c("All covar",
                "no_sal", "no_temp", "no_feces", "no_fdepth", "no_sedclass","no_pH")

m.plotbetas(mymodlist = Rnd1.Models, 
            mymodnames = Rnd1.Names)

m.plotHyperpar(mymodlist = Rnd1.Models,
               mymodnames = Rnd1.Names)

# Select model
dic.r1  <- c(Msp.full$dic$dic, 
             M.01$dic$dic, M.02$dic$dic, M.03$dic$dic, M.04$dic$dic, M.05$dic$dic, M.06$dic$dic)
waic.r1 <- c(Msp.full$waic$waic, 
             M.01$waic$waic, M.02$waic$waic, M.03$waic$waic, M.04$waic$waic, M.05$waic$waic, M.06$waic$waic)
Out.r1 <- cbind(dic.r1, waic.r1)
rownames(Out.r1) <- c("All covar", 
                      "no_sal", "no_temp", "no_feces", "no_fdepth", "no_sedclass","no_pH")
Out.r1
#              dic.r1   waic.r1
# All covar   461.3181 460.0314
# no_sal      460.1909 458.7334
# no_temp     460.0762 458.5177
# no_feces    464.4012 462.3584
# no_fdepth   458.7541 457.0665
# no_sedclass 456.8232 453.3314  <- lowest DIC and WAIC; use this for full model in next round (i.e., drop 'sed_class')
# no_pH       476.7387 474.3460


### ii) Round 2 ----
# Assign new full model for this next round of selection
Mrnd2.full <- M.05
Mrnd2.full$waic$waic 

# Remove sal
f.07 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   #sal + 
                   temp + feces_yes + fdepth_over25 + 
                   f(w, model = spde.sp))     

M.07 <- inla(f.07,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.07)
print(M.07$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.07)
M.07$dic$dic
M.07$waic$waic


# Remove temp
f.08 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + 
                   #temp + 
                   feces_yes + fdepth_over25 + 
                   f(w, model = spde.sp))     

M.08 <- inla(f.08,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.08)
print(M.08$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.08)
M.08$dic$dic
M.08$waic$waic


# Remove feces
f.09 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + 
                   #feces_yes + 
                   fdepth_over25 + 
                   f(w, model = spde.sp))     

M.09 <- inla(f.09,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.09)
print(M.09$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.09)
M.09$dic$dic
M.09$waic$waic


# Remove fdepth
f.10 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + feces_yes + 
                   #fdepth_over25 + 
                   f(w, model = spde.sp))     

M.10 <- inla(f.10,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.10)
print(M.10$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.10)
M.10$dic$dic
M.10$waic$waic


# Remove pH
str(Stack.nopH, max.level = 2) #Use a different stack

f.11 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   #XpH1 + XpH2 + 
                   sal + temp + feces_yes + fdepth_over25 + 
                   f(w, model = spde.sp))     

M.11 <- inla(f.11,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                                    # <- Don't forget to remove the lincomb for smoother
             data = inla.stack.data(Stack.nopH),                   # <- Don't forget to change the stack!!!
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH))) #<- Don't forget to change the stack

summary(M.11)
print(M.11$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.11)
M.11$dic$dic
M.11$waic$waic


# Round 2 model selection ----
Rnd2.Models <- list("Full" = Mrnd2.full,
                    "no_sal" = M.07, 
                    "no_temp" = M.08, 
                    "no_feces" = M.09, 
                    "no_fdepth" = M.10, 
                    "no_pH" = M.11)

Rnd2.Names <- c("Full",
                "no_sal", "no_temp", "no_feces", "no_fdepth", "no_pH")

m.plotbetas(mymodlist = Rnd2.Models, 
            mymodnames = Rnd2.Names)

m.plotHyperpar(mymodlist = Rnd2.Models,
               mymodnames = Rnd2.Names)

# Select model
dic.r2  <- c(Mrnd2.full$dic$dic, 
             M.07$dic$dic, M.08$dic$dic, M.09$dic$dic, M.10$dic$dic, M.11$dic$dic)
waic.r2 <- c(Mrnd2.full$waic$waic, 
             M.07$waic$waic, M.08$waic$waic, M.09$waic$waic, M.10$waic$waic, M.11$waic$waic)
Out.r2 <- cbind(dic.r2, waic.r2)
rownames(Out.r2) <- c("Full", 
                      "no_sal", "no_temp", "no_feces", "no_fdepth", "no_pH")
Out.r2
#           dic.r2   waic.r2
#Full      456.8232 453.3314
#no_sal    456.0475 452.4480
#no_temp   455.8593 452.1228
#no_feces  460.1897 455.9761
#no_fdepth 454.1015 450.2617  <- lowest DIC and WAIC; use this for full model in next round (i.e., drop 'fdepth') 
#no_pH     472.9864 467.9727


### iii) Round 3 ----
# Assign new full model for this next round of selection
Mrnd3.full <- M.10 #Remove fdepth
Mrnd3.full$waic$waic 

# Remove sal
f.12 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   #sal + 
                   temp + feces_yes + 
                   f(w, model = spde.sp))     

M.12 <- inla(f.12,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.12)
print(M.12$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.12)
M.12$dic$dic
M.12$waic$waic


# Remove temp
f.13 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + 
                   #temp + 
                   feces_yes + 
                   f(w, model = spde.sp))     

M.13 <- inla(f.13,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.13)
print(M.13$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.13)
M.13$dic$dic
M.13$waic$waic


# Remove feces
f.14 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + temp + 
                   #feces_yes + 
                   f(w, model = spde.sp))     

M.14 <- inla(f.14,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.14)
print(M.14$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.14)
M.14$dic$dic
M.14$waic$waic


# Remove pH
str(Stack.nopH, max.level = 2) #Use a different stack

f.15 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   #XpH1 + XpH2 + 
                   sal + temp + feces_yes + 
                   f(w, model = spde.sp))     

M.15 <- inla(f.15,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                                    # <- Don't forget to remove the lincomb for smoother
             data = inla.stack.data(Stack.nopH),                   # <- Don't forget to change the stack!!!
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH))) #<- Don't forget to change the stack

summary(M.15)
print(M.15$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.15)
M.15$dic$dic
M.15$waic$waic


# Round 3 model selection ----
# Double check nothing funky going on when adding/removing different covariates
Rnd3.Models <- list("Full" = Mrnd3.full,
                    "no_sal" = M.12, 
                    "no_temp" = M.13, 
                    "no_feces" = M.14, 
                    "no_pH" = M.15)

Rnd3.Names <- c("Full",
                "no_sal", "no_temp", "no_feces", "no_pH")

m.plotbetas(mymodlist = Rnd3.Models, 
            mymodnames = Rnd3.Names)

m.plotHyperpar(mymodlist = Rnd3.Models,
               mymodnames = Rnd3.Names)

# Select model
dic.r3  <- c(Mrnd3.full$dic$dic, 
             M.12$dic$dic, M.13$dic$dic, M.14$dic$dic, M.15$dic$dic)
waic.r3 <- c(Mrnd3.full$waic$waic, 
             M.12$waic$waic, M.13$waic$waic, M.14$waic$waic, M.15$waic$waic)
Out.r3 <- cbind(dic.r3, waic.r3)
rownames(Out.r3) <- c("Full", 
                      "no_sal", "no_temp", "no_feces", "no_pH")
Out.r3
#          dic.r3   waic.r3
#Full     454.1015 450.2617
#no_sal   453.5217 449.4441
#no_temp  453.2793 449.1375  <- Lowest DIC/WAIC of the equivalent models
#no_feces 457.4783 452.7579  
#no_pH    470.3556 465.0249  


### iv) Round 4 ----
# Assign new full model for this next round of selection
Mrnd4.full <- M.13 #Remove temp
Mrnd4.full$waic$waic 

# Remove sal (M.16) ----
f.16 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   #sal + 
                   feces_yes + 
                   f(w, model = spde.sp))     

M.16 <- inla(f.16,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.16)
print(M.16$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.16)
M.16$dic$dic
M.16$waic$waic


# Remove feces
f.17 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal + 
                   #feces_yes + 
                   f(w, model = spde.sp))     

M.17 <- inla(f.17,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.17)
print(M.17$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.17)
M.17$dic$dic
M.17$waic$waic


# Remove pH
str(Stack.nopH, max.level = 2) #Use a different stack

f.18 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   #XpH1 + XpH2 + 
                   sal + feces_yes + 
                   f(w, model = spde.sp))     

M.18 <- inla(f.18,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                                    # <- Don't forget to remove the lincomb for smoother
             data = inla.stack.data(Stack.nopH),                   # <- Don't forget to change the stack!!!
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH))) #<- Don't forget to change the stack

summary(M.18)
print(M.18$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.18)
M.18$dic$dic
M.18$waic$waic


# Round 4 model selection ----
Rnd4.Models <- list("Full" = Mrnd4.full,
                    "no_sal" = M.16, 
                    "no_feces" = M.17, 
                    "no_pH" = M.18)

Rnd4.Names <- c("Full",
                "no_sal", "no_feces", "no_pH")

m.plotbetas(mymodlist = Rnd4.Models, 
            mymodnames = Rnd4.Names)

m.plotHyperpar(mymodlist = Rnd4.Models,
               mymodnames = Rnd4.Names)

# Select model
dic.r4  <- c(Mrnd4.full$dic$dic, 
             M.16$dic$dic, M.17$dic$dic, M.18$dic$dic)
waic.r4 <- c(Mrnd4.full$waic$waic, 
             M.16$waic$waic, M.17$waic$waic, M.18$waic$waic)
Out.r4 <- cbind(dic.r4, waic.r4)
rownames(Out.r4) <- c("Full", 
                      "no_sal", "no_feces", "no_pH")
Out.r4
#          dic.r4   waic.r4
#Full     453.2793 449.1375
#no_sal   452.6718 448.3825 <- DIC/WAIC equivalent (within +/-2) of the full model but has more parameters
#no_feces 456.4960 451.5363
#no_pH    469.3660 463.7569


### v) Round 5 ----
# Assign new full model for this next round of selection
Mrnd5.full <- M.16 #Remove salinity
Mrnd5.full$waic$waic 

# Remove feces
f.19 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   #feces_yes + 
                   f(w, model = spde.sp))     

M.19 <- inla(f.19,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.19)
print(M.19$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.19)
M.19$dic$dic
M.19$waic$waic

# Remove pH
str(Stack.nopH, max.level = 2) #Use a different stack

f.20 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   #XpH1 + XpH2 + 
                   feces_yes + 
                   f(w, model = spde.sp))     

M.20 <- inla(f.20,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                                    # <- Don't forget to remove the lincomb for smoother
             data = inla.stack.data(Stack.nopH),                   # <- Don't forget to change the stack!!!
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH))) #<- Don't forget to change the stack

summary(M.20)
print(M.20$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.20)
M.20$dic$dic
M.20$waic$waic


# Round 5 model selection ----
Rnd5.Models <- list("Full" = Mrnd5.full,
                    "no_feces" = M.19, 
                    "no_pH" = M.20)

Rnd5.Names <- c("Full",
                "no_feces", "no_pH")

m.plotbetas(mymodlist = Rnd5.Models, 
            mymodnames = Rnd5.Names)

m.plotHyperpar(mymodlist = Rnd5.Models,
               mymodnames = Rnd5.Names)

# Select model
dic.r5  <- c(Mrnd5.full$dic$dic, 
             M.19$dic$dic, M.20$dic$dic)
waic.r5 <- c(Mrnd5.full$waic$waic, 
             M.19$waic$waic, M.20$waic$waic)
Out.r5 <- cbind(dic.r5, waic.r5)
rownames(Out.r5) <- c("Full", 
                      "no_feces", "no_pH")
Out.r5
#          dic.r5   waic.r5
#Full     452.6718 448.3825 <- lowest DIC/WAIC 
#no_feces 455.8298 450.6320
#no_pH    468.7447 463.2006

# Stop backwards model selection 
# The final model is Mrnd5.full (= M.16) with only feces and pH smoother as covariates
summary(Mrnd5.full)
print(Mrnd5.full$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(Mrnd5.full)
Mrnd5.full$dic$dic
Mrnd5.full$waic$waic



### (5) Forwards Model Selection  --------------------------------------------------------------------------------------------------------------------------------------

# i) Round 1 forward model selection ----

# Intercept only
str(Stack.nopH, max.level = 2) #Use a different stack

f.ionly <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   f(w, model = spde.sp))     

M.ionly <- inla(f.ionly,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.ionly)
print(M.ionly$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.ionly)
M.ionly$dic$dic
M.ionly$waic$waic

# Add sal
f.21 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                  sal +
                  f(w, model = spde.sp))     

M.21 <- inla(f.21,
                family = "binomial",
                Ntrials = 1,  
                #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
                data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
                control.inla = list(strategy = 'simplified.laplace'), 
                control.compute = list(dic = TRUE, waic = TRUE),
                control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.21)
print(M.21$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.21)
M.21$dic$dic
M.21$waic$waic


# Add temp
f.22 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   temp +
                   f(w, model = spde.sp))     

M.22 <- inla(f.22,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.22)
print(M.22$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.22)
M.22$dic$dic
M.22$waic$waic


# Add feces
f.23 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   feces_yes +
                   f(w, model = spde.sp))     

M.23 <- inla(f.23,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.23)
print(M.23$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.23)
M.23$dic$dic
M.23$waic$waic

# Add fdepth
f.24 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   fdepth_over25 +
                   f(w, model = spde.sp))     

M.24 <- inla(f.24,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.24)
print(M.24$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.24)
M.24$dic$dic
M.24$waic$waic


# Add sedclass
f.25 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   sedclass_mud + sedclass_grit + sedclass_rock +
                   f(w, model = spde.sp))     

M.25 <- inla(f.25,
             family = "binomial",
             Ntrials = 1,  
             #lincomb = lcs.pH,                        #No smoother for pH, so remove the lincomb argument
             data = inla.stack.data(Stack.nopH),       #Use the Stack without effects for pH
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack.nopH)))  #Stack without effect for pH

summary(M.25)
print(M.25$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.25)
M.25$dic$dic
M.25$waic$waic


# Add pH as a smoother (need different stack!)
f.26 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   f(w, model = spde.sp))     

M.26 <- inla(f.26,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.26)
print(M.26$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.26)
M.26$dic$dic
M.26$waic$waic


# Model selection: Round 1 (forwards) ----
Rnd1f.Models <- list("intercept_only" = M.ionly,
                    "add_sal" = M.21,
                    "add_temp" = M.22,
                    "add_feces" = M.23, 
                    "add_fdepth" = M.24,
                    "add_sedclass" = M.25,
                    "add_pH" = M.26)

Rnd1f.Names <- c("intercept_only",
                "add_sal", "add_temp","add_feces", "add_fdepth", "add_sedclass", "add_pH")

m.plotbetas(mymodlist = Rnd1f.Models, 
            mymodnames = Rnd1f.Names)

m.plotHyperpar(mymodlist = Rnd1f.Models,
               mymodnames = Rnd1f.Names)

# Select model
dic.r1f  <- c(M.ionly$dic$dic, 
             M.21$dic$dic, M.22$dic$dic, M.23$dic$dic, M.24$dic$dic, M.25$dic$dic, M.26$dic$dic)
waic.r1f <- c(M.ionly$waic$waic, 
             M.21$waic$waic, M.22$waic$waic, M.23$waic$waic, M.24$waic$waic, M.25$waic$waic, M.26$waic$waic)
Out.r1f <- cbind(dic.r1f, waic.r1f)
rownames(Out.r1f) <- c("intercept_only",
                       "add_sal", "add_temp","add_feces", "add_fdepth", "add_sedclass", "add_pH")
Out.r1f
#                dic.r1f  waic.r1f
#intercept_only 474.4392 467.8125
#add_sal        475.4424 468.8320
#add_temp       475.5255 469.1224
#add_feces      468.7455 463.2006
#add_fdepth     477.2887 471.0821
#add_sedclass   478.3853 474.4800
#add_pH         455.8306 450.6324  <- Lowest DIC/WAIC, select this for model for the next round


# ii) Round 2 forward model selection ----

# Optimal model selected from last round (add_pH smoother)
Mfrnd2.opt <- M.26
M.26$waic$waic  

# Add sal
f.27 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sal +
                   f(w, model = spde.sp))     

M.27 <- inla(f.27,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.27)
print(M.27$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.27)
M.27$dic$dic
M.27$waic$waic


# Add temp
f.28 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   temp +
                   f(w, model = spde.sp))     

M.28 <- inla(f.28,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.28)
print(M.28$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.28)
M.28$dic$dic
M.28$waic$waic

# Add feces
f.29 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   feces_yes +
                   f(w, model = spde.sp))     

M.29 <- inla(f.29,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.29)
print(M.29$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.29)
M.29$dic$dic
M.29$waic$waic

# Add fdepth
f.30 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   fdepth_over25 +
                   f(w, model = spde.sp))     

M.30 <- inla(f.30,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.30)
print(M.30$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.30)
M.30$dic$dic
M.30$waic$waic

# Add sedclass
f.31 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + 
                   sedclass_mud + sedclass_grit + sedclass_rock +
                   f(w, model = spde.sp))     

M.31 <- inla(f.31,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.31)
print(M.31$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.31)
M.31$dic$dic
M.31$waic$waic


# Model selection: Round 2 (forwards) ----
Rnd2f.Models <- list("null_model" = Mfrnd2.opt,
                     "add_sal" = M.27,
                     "add_temp" = M.28,
                     "add_feces" = M.29, 
                     "add_fdepth" = M.30,
                     "add_sedclass" = M.31)

Rnd2f.Names <- c("null_model",
                 "add_sal", "add_temp","add_feces", "add_fdepth", "add_sedclass")

m.plotbetas(mymodlist = Rnd2f.Models, 
            mymodnames = Rnd2f.Names)

m.plotHyperpar(mymodlist = Rnd2f.Models,
               mymodnames = Rnd2f.Names)

# Select model
dic.r2f  <- c(Mfrnd2.opt$dic$dic, 
              M.27$dic$dic, M.28$dic$dic, M.29$dic$dic, M.30$dic$dic, M.31$dic$dic)
waic.r2f <- c(Mfrnd2.opt$waic$waic, 
              M.27$waic$waic, M.28$waic$waic, M.29$waic$waic, M.30$waic$waic, M.31$waic$waic)
Out.r2f <- cbind(dic.r2f, waic.r2f)
rownames(Out.r2f) <- c("null_model",
                       "add_sal", "add_temp","add_feces", "add_fdepth", "add_sedclass")
Out.r2f
#                dic.r2f  waic.r2f
#null_model     455.8306 450.6324
#add_sal        456.4925 451.5343
#add_temp       456.7675 451.7841
#add_feces      452.6694 448.3769 <- lower DIC and just barely WAIC; use for next round
#add_fdepth     458.4714 453.7724
#add_sedclass   459.8106 456.9845


# iii) Round 3 forward model selection ----

# Optimal model selected from last round (add_feces)
Mfrnd3.opt <- M.29
M.29$waic$waic  

# Add sal
f.32 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + feces_yes +
                   sal +
                   f(w, model = spde.sp))     

M.32 <- inla(f.32,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.32)
print(M.32$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.32)
M.32$dic$dic
M.32$waic$waic

# Add temp
f.33 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + feces_yes +
                   temp +
                   f(w, model = spde.sp))     

M.33 <- inla(f.33,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.33)
print(M.33$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.33)
M.33$dic$dic
M.33$waic$waic

# Add fdepth
f.34 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + feces_yes +
                   fdepth_over25 +
                   f(w, model = spde.sp))     

M.34 <- inla(f.34,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.34)
print(M.34$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.34)
M.34$dic$dic
M.34$waic$waic

# Add sedclass
f.35 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                   XpH1 + XpH2 + feces_yes +
                   sedclass_mud + sedclass_grit + sedclass_rock +
                   f(w, model = spde.sp))     

M.35 <- inla(f.35,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack1),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.35)
print(M.35$summary.fixed[,c("mean", "0.025quant", "0.975quant")] , digits = 3)
bri.hyperpar.summary(M.35)
M.35$dic$dic
M.35$waic$waic


# Model selection: Round 3 (forwards) ----
Rnd3f.Models <- list("null_model" = Mfrnd3.opt,
                     "add_sal" = M.32,
                     "add_temp" = M.33,
                     "add_fdepth" = M.34,
                     "add_sedclass" = M.35)

Rnd3f.Names <- c("null_model",
                 "add_sal", "add_temp", "add_fdepth", "add_sedclass")

m.plotbetas(mymodlist = Rnd3f.Models, 
            mymodnames = Rnd3f.Names)

m.plotHyperpar(mymodlist = Rnd3f.Models,
               mymodnames = Rnd3f.Names)

# Select model
dic.r3f  <- c(Mfrnd3.opt$dic$dic, 
              M.32$dic$dic, M.33$dic$dic, M.34$dic$dic, M.35$dic$dic)
waic.r3f <- c(Mfrnd3.opt$waic$waic, 
              M.32$waic$waic, M.33$waic$waic, M.34$waic$waic, M.35$waic$waic)
Out.r3f <- cbind(dic.r3f, waic.r3f)
rownames(Out.r3f) <- c("null_model",
                       "add_sal", "add_temp", "add_fdepth", "add_sedclass")
Out.r3f
#              dic.r3f  waic.r3f
#null_model   452.6694 448.3769  <- lowest DIC/WAIC of the equivalent models and has less parameters; this is final model
#add_sal      453.2810 449.1395
#add_temp     453.5210 449.4431
#add_fdepth   455.0957 451.3139
#add_sedclass 456.9489 454.9612

identical(print(Mfrnd3.opt$summary.fixed, digits = 3), print(Mrnd5.full$summary.fixed, digits = 3))
# The results of the two final models from backwards vs forward selection
# are very slightly different at the third decimal point
# So I will chose to use the model generated from backwards selection in model validation, for consistency sake

# FINAL MODEL (M.16 = Mrnd5.full ~ Mfrnd3.opt) -----
f.16
M.final <- Mrnd5.full
#save(M.final, file = "Results/Final-inla-objects/Extsed_gam_spatial_meshsp_n21958_finalmodel.rdata")


### (6) Model Validation --------------------------------------------------------------------------------------------------------------------------------------


### i) Obtain the Pearson residuals ----
N  #Sample size

# Obtain fitted values from INLA model
Pi1 <- M.final$summary.fitted.values[1:N, "mean"]  #includes fixed effects + random effects


# Calculate Bernoulli residuals by hand
ExpY <- Pi1
varY <- Pi1 *(1 - Pi1)
E1 <- (ext.sed$matrix01 - ExpY) / (sqrt(varY))


# Or for easier plotting
fmodlist <- list("Final Model" = M.final)
myresids <- m.getresid.bern(MyData = ext.sed, MyYvar = ext.sed$matrix01, MyINLAlist = fmodlist)


### ii) Plot residuals vs fitted values ----
ggplot(myresids, aes(x = ExpY, y = E1)) +
  geom_point() +
  labs(x = "Fitted values", y = "Residuals")


### iii) Plot residuals vs. observed values ----
ggplot(myresids, aes(x = ObsY, y = E1)) +
  geom_point() +
  labs(x = "Observed values", y = "Residuals")


### iv) Plot residuals vs every covariate in the model ----

# pH
ggplot(myresids, aes(x = pH, y = E1)) +
  geom_point() +
  geom_smooth() +
  labs(x = "pH", y = "Residuals")
  #No non-linear pattern remaining

# feces (categorical)
myresids <- mutate(myresids,
                      nfeces = as.numeric(feces))

ggplot(myresids, aes(x = feces, y = E1)) +
  geom_boxplot() +
  labs(x = "feces", y = "Residuals")

ggplot(myresids, aes(x = nfeces, y = E1)) +
  geom_jitter(alpha = 0.4, width = 0.2, height = 0.2) +
  labs(x = "feces", y = "Residuals")


### v) Plot residuals vs covariates not in the model ----

# Salinity
ggplot(myresids, aes(x = sal, y = E1)) +
  geom_point() +
  geom_smooth() +
  labs(x = "salinity", y = "Residuals")

# Temperature
ggplot(myresids, aes(x = temp, y = E1)) +
  geom_point() +
  geom_smooth() +
  labs(x = "temperature", y = "Residuals")

# Water Depth (above or below 25cm)
ggplot(myresids, aes(x = fdepth, y = E1)) +
  geom_boxplot() +
  labs(x = "water depth", y = "Residuals")

# Sediment Type 
ggplot(myresids, aes(x = sedclass, y = E1)) +
  geom_boxplot() +
  labs(x = "sediment type", y = "Residuals")

# Wetland
ggplot(myresids, aes(x = wtld, y = E1)) +
  geom_boxplot() +
  labs(x = "wetland", y = "Residuals")



### vii) Check residuals for temporal, spatial and spatial-temporal dependency ----

### TEMPORAL
# Study Month
ggplot(myresids, aes(x = mon_yr, y = E1)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Month of Study", y = "Residuals")

# Study Year
ggplot(myresids, aes(x = study_yr, y = E1)) +
  geom_jitter() +
  geom_smooth() +
  labs(x = "Year of Study", y = "Residuals")


### SPATIAL (i.e., averaged over both study years)

#Hard to interpret a variogram due to different spatial scales
  #within and between wetlands

#Large scale: Make dataframe that summarizes residuals by wetland location
res.cols <- c("2" = "darkred",
              "1" = "seagreen")

t <- myresids %>% 
        group_by(wtld) %>% 
        summarize(res = mean(E1),
                  east = mean(east),
                  north = mean(north),
                  Sign = as.factor(ifelse(res > 0, 1, 2))) %>% 
        mutate(MyCex = abs(res) / max(abs(res)))  #Note this is scaled, not just the absolute value...

ggplot(t, aes(x = east, y = north, size = MyCex, fill = Sign)) +
        geom_point(alpha = 0.7, pch = 21) +
        theme_light() +
        scale_fill_manual(values = res.cols, labels = c("Neg", "Pos")) +
        labs(x = "Easting", y = "Northing",
             title = paste("Spatial distribution of residuals"))

#Small-scale: plot the INLA spatial field
wt.pm <- M.final$summary.random$w$mean

m.gg2PlotField(field = wt.pm,
               mesh = mesh.sp,
               ColorScheme = 1,                #You can use the values 1 to 4
               xlim = range(mesh.sp$loc[ ,1]), #Either need to add xlim and ylim, OR add a ContourMap
               ylim = range(mesh.sp$loc[ ,2]),
               Xlabel = "Easting",
               Ylabel = "Northing",
               MyTitle = "Final model - Posterior mean values w")


### SPATIAL-TEMPORAL

#Large scale: Make dataframe that summarizes residuals by wetland location, by study year
res.cols <- c("2" = "darkred",
              "1" = "seagreen")

st <- myresids %>% 
          group_by(study_yr, wtld) %>%  
          summarize(res = mean(E1),
                    east = mean(east),
                    north = mean(north),
                    Sign = as.factor(ifelse(res > 0, 1, 2))) %>% 
          mutate(MyCex = abs(res) / max(abs(res)))  #Note this is scaled, not just the absolute value...

ggplot(st, aes(x = east, y = north, size = MyCex, fill = Sign)) +
      geom_point(alpha = 0.7, pch = 21) +
      theme_light() +
      facet_grid(study_yr ~.) +    
      scale_fill_manual(values = res.cols, labels = c("Neg", "Pos")) +
      labs(x = "Easting", y = "Northing",
           title = "Spatial-temporal distribution of residuals",
           subtitle = "plotted by year of study")


