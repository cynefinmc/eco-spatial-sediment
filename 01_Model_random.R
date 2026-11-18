###########################################################################################################
###                                                                                                     ###
###                ECOLOGICAL REGRESSION MODEL OF EXTENSIVE SEDIMENT MATRIX RESULTS                     ###
###                                                                                                     ###
###          This script determines the most appropriate RANDOM structure for the model                 ###
###                                                                                                     ###
###########################################################################################################

#' Much of this code is heavily based on code written by Highland Statistics Ltd.
#' Code comes from Volume I and II of Beginner's Guide to Spatial, Temporal and Spatial-Temporal Ecological
#'   Data Analysis with R-INLA by Zuur, Ieno and Saveliev, and associated courses.
#' Books and courses can be found at: https://www.highstat.com


### (1) Load packages and data   ------------------------------------------------------------------------------
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

my.wd <- getwd()
source(paste0(my.wd, "/00_aiv_ext_functions"))

ext.sed <- read.csv(file = paste0(my.wd,
                                  "/Extensive_sediment_with_ebird_MODIS_summaries_final.csv"))

### (2) Select final data used for spatial model ------------------------------------------------------------

#' Shorten some covariate names
ext.sed <- rename(ext.sed,
                  "sal" = sal_ppt,       #' salinity
                  "temp" = temp_C,       #' temperature
                  "fdepth" = fdepth_cm)  #' sediment depth as a factor

#' Remove the ebird summary covariates, as at the level of wetland and have missing data, 
#' and MODIS non-water dominant landcover due to collinearity
ext.sed <- select(ext.sed,
                  -c("Matrix_positivity", "N_checklist", "Total_abundance", "Relative_abundance", "Species_richness",
                     "dominant_MODIS_landcover_class"))


### (3) Frequentist models for comparison  -------------------------------------------------------------------
#' Helps to estimate what to use for sigma priors

#' GLM only without covariates
Mfrq1 <- glm(matrix01 ~ 1,          
              data = ext.sed,
              family = "binomial")
summary(Mfrq1)
#' the se around the intercept is ~0.1, so sigma of ~1 to 2 with a logistic model is reasonable to capture variation

#' GLMM with only wtld and without covariates 
Mfrq2 <- glmer(matrix01 ~ 1 + (1 | wtld), #wtld only
               data = ext.sed,
               family = "binomial")
summary(Mfrq2)
#' the se around the intercept is ~0.4, so sigma of ~1 to 2 with a logistic model is still reasonable to capture variation


#' Nested GLMM without covariates 
Mfrq3 <- glmer(matrix01 ~ 1 + (1 | wtld/site), #site nested within wtld
               data = ext.sed,
               family = "binomial")
summary(Mfrq3)

#' The full GLMM model with all covariates 
Mfrq4 <- glmer(matrix01 ~ pH + sal + temp + feces + fdepth + sedclass +
                  (1 | wtld/site), 
                  data = ext.sed,
                  family = "binomial")
summary(Mfrq4)

#' The full GLMM model with wetland as only random effect
Mfrq5 <- glmer(matrix01 ~ pH + sal + temp + feces + fdepth + sedclass +
                 (1 | wtld), 
               data = ext.sed,
               family = "binomial")
summary(Mfrq5)


#' GAMM with wetland as random effect and covariates
Mfrq6 <- gamm4(matrix01 ~ s(pH) + sal + temp + feces + fdepth + sedclass,
               random = ~ (1 | wtld),
               data = ext.sed,
               family = binomial)
summary(Mfrq6$gam)
summary(Mfrq6$mer)
plot(Mfrq6$gam, pages = 1)


### (4) INLA Model formulation: choosing the dependency structure  ---------------------------------------------------------------------------------------------------------------

### Set up coding for models:
# i) Make smoother for the pH covariate  ----

#' Get basis function using smoothCon
  #' This is a unpenalized (fx = T) cubic regression spline (bs = "cr") smoother
  #' k is dim of basis to represent the smooth term, see ?s
Lego.pH <- smoothCon(s(pH, bs = "cr", k = 3, fx = TRUE), #choosing 3 knots as expect 2 edf based on gamm4 model
                     data = ext.sed,  
                     absorb.cons = TRUE)[[1]]
XpH <- Lego.pH$X 

#' Turn into data frame for lincomb function
XpH.df <- data.frame(XpH)
names(XpH.df) <- paste0("XpH", 1:ncol(XpH))

lcs.pH <- inla.make.lincombs(XpH.df)


# ii) Make an 'X' matrix for linear covariates  ----

#' Xm for most models
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
RangeGuess <- 100 #in meters, best guess based on where I observed birds congregating within the wetlands
MaxEdge <- RangeGuess / 5 
ELoc <- cbind(ext.sed$east, ext.sed$north) 
OuterBound <- 1500
bound <- inla.nonconvex.hull(ELoc, convex = -0.003, 
                               resolution = c(340, 115)) #Increased resolution as per INLA warning

# Small mesh for troubleshooting
mesh.sm <- inla.mesh.2d(ELoc,
                       boundary = bound, 
                       max.edge = c(5* 1, 1000 * 5) * MaxEdge, #Try to adjust for fast processing time
                       cutoff = 5* MaxEdge / 5,                #Try to adjust for fast processing time
                       offset = c(MaxEdge, OuterBound))  
plot(mesh.sm, asp = 1) 
points(ELoc, col = 2, pch = 16, cex = 1) 
mesh.sm$n  #n = 4,804 


# Final mesh 
mesh.f <- inla.mesh.2d(ELoc,
                       boundary = bound, 
                       max.edge = c(2* 1, 100 * 5) * MaxEdge,
                       cutoff = 2* MaxEdge / 5,         
                       offset = c(MaxEdge, OuterBound))  
plot(mesh.f, asp = 1) 
points(ELoc, col = 2, pch = 16, cex = 1) 
mesh.f$n  #n = 21,958



# iv) Make projector matrices for spatial models  ----
A.f <- inla.spde.make.A(mesh.f, loc = ELoc)


# v) Define SPDE and spatial field (w.index) for spatial models  ----
spde.f <- inla.spde2.pcmatern(mesh.f,
                              prior.range = c(100, 0.05),   
                              prior.sigma = c(2, 0.05))

w.index.f <- inla.spde.make.index(name = "w",
                                  n.spde = spde.f$n.spde)


# vii) Define penalized complexity priors for random effects (intercepts)  ----
  #' pc priors for the Spatial Random Field will be defined while making the SPDE
summary(Mfrq2)
summary(Mfrq3)
  #As sd for 'site' is ~0.8 and sd for 'wtld' is ~1.5 in both models, sigma less than 2 is reasonable
  # P(sigma_site > 2) = 0.05  and P(sigma_wtld > 2) = 0.05 
priorpc <- list(prec = list(prior = "pc.prec", param = c(2, 0.05)))

# viii) Define the number of observations
N <- nrow(ext.sed)



### 1. Bernoulli GAMM with wtld only ---- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# Use this formula so can use variables as listed in the Xm Stack
f1 <-  formula(matrix01 ~ -1 + Intercept +  #Add -1 and + Intercept, as defined Intercept separately in Stack above
                 XpH1 + XpH2 + 
                 sal + temp +     
                 feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +   #Needs to be the same as in the 'Xm' object 
                 f(wtld, model = "iid", hyper = priorpc))    #Random intercept for 'wtld'

head(Xm)
N <- nrow(ext.sed)
Stack1 <- inla.stack(tag = "Fit1",
                     data = list(matrix01 = ext.sed$matrix01), 
                     A = list(1, 1, 1, 1),
                     effects = list(
                       Intercept = rep(1, N), 
                       XpH = XpH.df,
                       Xm = Xm[ ,-1],
                       wtld = ext.sed$wtld)) 

M.gamm1 <- inla(f1,
                family = "binomial",
                Ntrials = 1,     
                lincomb = lcs.pH,
                data = inla.stack.data(Stack1), 
                control.inla = list(strategy = 'simplified.laplace'),  
                control.compute = list(dic = TRUE, waic = TRUE),       
                #control.predictor = list(compute = T),                 #To get fitted values
                control.predictor = list(A = inla.stack.A(Stack1)))

summary(M.gamm1) 

mbetas <- M.gamm1$summary.fixed[,c("mean", "0.025quant", "0.975quant")] 
print(mbetas, digits = 3) 
bri.hyperpar.summary(M.gamm1)  
M.gamm1$dic$dic
M.gamm1$waic$waic


# Plot the smoother (on scale of linear covariates)
Edata <- ext.sed

SmoothpH   <- M.gamm1$summary.lincomb.derived
Edata$fpH  <- SmoothpH[, "mean"]
Edata$SeLo <- SmoothpH[, "0.025quant"]
Edata$SeUp <- SmoothpH[, "0.975quant"]

p <- ggplot()
p <- p + xlab("pH") + ylab("Smoother")
p <- p + theme(text = element_text(size=15)) 
p <- p + geom_line(data = Edata, 
                   aes(x = pH, 
                       y = fpH),
                   size = 1)
p <- p + geom_ribbon(data = Edata, 
                     aes(x = pH, 
                         ymax = SeUp, 
                         ymin = SeLo),
                     alpha = 0.6)
p



### 2. Bernoulli GAM + spatial field (.so = "spatial only") ---- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

mesh.f$n                       #Done; 21958
dim(A.f)                       #Done; 701 21958
str(spde.f, max.level = 1)     #Already done
str(w.index.f, max.level = 1)  #Already done
head(Xm)                       #Already done
head(XpH)                      #Already done
head(XpH.df)                   #Already done
N                              #Already done

Stack4 <- inla.stack(tag = "Fit4",
                     data = list(matrix01 = ext.sed$matrix01), 
                     A = list(1, 1, 1, A.f),
                     effects = list(
                       Intercept = rep(1, N), 
                       XpH = XpH.df, 
                       Xm = Xm[ ,-1], 
                       w = w.index.f)) 

f4 <-  formula(matrix01 ~ -1 + Intercept + #Add -1 and + Intercept
                 XpH1 + XpH2 + 
                 sal + temp +     
                 feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                 f(w, model = spde.f))     #Spatial Random Field

M.soF <- inla(f4,
             family = "binomial",
             Ntrials = 1,  
             lincomb = lcs.pH,
             data = inla.stack.data(Stack4),
             control.inla = list(strategy = 'simplified.laplace'), 
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack4)))

summary(M.soF)
bri.hyperpar.summary(M.soF)


### 3. Bernoulli GAMM with wetland + spatial field (.ws = "wetland and space") ---- ~~~~~~~~~~~~~~~~~~~~~~~~~~~###

mesh.f$n                       #Done; 21958
dim(A.f)                       #Done; 701 21958
str(spde.f, max.level = 1)     #Already done
str(w.index.f, max.level = 1)  #Already done
head(Xm)                       #Already done
head(XpH)                      #Already done
head(XpH.df)                   #Already done
N                              #Already done

Stack5 <- inla.stack(tag = "Fit5",
                     data = list(matrix01 = ext.sed$matrix01), 
                     A = list(1, 1, 1, 1, A.f),
                     effects = list(
                       Intercept = rep(1, N), 
                       XpH = XpH.df,
                       Xm = Xm[ ,-1], 
                       wtld = ext.sed$wtld,
                       w = w.index.f)) 

f5 <-  formula(matrix01 ~ -1 + Intercept +      #Add -1 and + Intercept
                 XpH1 + XpH2 +
                 temp + sal +
                 feces_yes + fdepth_over25 + sedclass_mud + sedclass_grit + sedclass_rock +    
                 f(wtld, model = "iid", hyper = priorpc) +   #Random intercept for 'wtld'
                 f(w, model = spde.f))                        #Spatial Random Field

M.wsF <- inla(f5,
             family = "binomial",
             Ntrials = 1,    
             lincomb = lcs.pH,
             data = inla.stack.data(Stack5),
             control.inla = list(strategy = 'simplified.laplace'), #default option
             control.compute = list(dic = TRUE, waic = TRUE),
             control.predictor = list(A = inla.stack.A(Stack5)))

summary(M.wsF)
bri.hyperpar.summary(M.wsF)




### (5) Model comparision via DIC/WAIC (and parameter estimates) between potential dependency structures -----------------------------------------------------

### i) DIC/WAIC ----
dic.f  <- c(M.gamm1$dic$dic, 
            M.soF$dic$dic, M.wsF$dic$dic)  

waic.f <- c(M.gamm1$waic$waic,
            M.soF$waic$waic, M.wsF$waic$waic) 

Out.f <- cbind(dic.f, waic.f)
rownames(Out.f) <- c("GAMM wtld",
                     "GAM + spatial", "GAMM wtld + spatial") 
Out.f


### ii) Compare posterior means of hyper-parameters ----
bri.hyperpar.summary(M.gamm1)
bri.hyperpar.summary(M.soF)
bri.hyperpar.summary(M.wsF)


### iii) Compare posterior means of betas ----
AllModels <- list(
  M.gamm1$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
  M.soF$summary.fixed[, c("mean", "0.025quant", "0.975quant")],
  M.wsF$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
  
ModelNames <- c("GAMM wtld",
                "GAM + spatial", "GAMM wtld + spatial")

mCompareBetasofModels(AllModels = AllModels, 
                      ModelNames = ModelNames)


### (6) Model validation for comparison between potential dependency structures --------------------------------------------------------

### i) Get Pearson residuals ----
MyfinalINLAs <- list("M.gamm1" = M.gamm1,
                     "M.soF" = M.soF,
                     "M.wsF" = M.wsF) 

# Extract the residuals and bind to variables in the original data for plotting
myf.resids <- m.getresid.bern(MyData = ext.sed, MyYvar = ext.sed$matrix01, MyINLAlist = MyfinalINLAs)

# Make a vector of model names
myfinalnames <- unique(myf.resids$model)

### ii) Plot Pearson residuals vs. Fitted values -----
resid_vs_fitted <- purrr::map(.x = myfinalnames,
                              .f = ~{
                               myf.resids %>% 
                                 filter(model == .x) %>% 
                                 ggplot() +
                                 geom_point(aes(x = ExpY, y = E1)) +
                                 labs(x = "Fitted values",
                                      y = "Pearson residuals",
                                      title = .x)
                              })

resid_vs_fitted[1]
resid_vs_fitted[2]
resid_vs_fitted[3]
resid_vs_fitted[4]
resid_vs_fitted[5]


### iii) Plot fitted values vs. observed values (to check for overfitting) ----
obs_vs_fitted <- purrr::map(.x = myfinalnames,
                            .f = ~{
                                myf.resids %>% 
                                  filter(model == .x) %>% 
                                  ggplot() +
                                  geom_point(aes(x = ExpY, y = ObsY)) +
                                  labs(x = "Fitted values",
                                       y = "Observed values",
                                       title = .x)
                              })

obs_vs_fitted[1]
obs_vs_fitted[2]
obs_vs_fitted[3]
obs_vs_fitted[4]
obs_vs_fitted[5]


### iv) Plot Pearson residual vs. covariates in or not in model ----

# Boxplot with wtld as that is the one I am the most concerned about
wtld_vs_resid <- purrr::map(.x = myfinalnames,
                            .f = ~{
                                myf.resids %>% 
                                  filter(model == .x) %>% 
                                  ggplot() +
                                  geom_boxplot(aes(x = wtld, y = E1)) +
                                  labs(x = "Wetland",
                                       y = "Pearson Residuals",
                                       title = .x)
                              })

wtld_vs_resid[1]
wtld_vs_resid[2]
wtld_vs_resid[3]
wtld_vs_resid[4]
wtld_vs_resid[5]



# Study year
year_vs_resid <- purrr::map(.x = myfinalnames,
                            .f = ~{
                              myf.resids %>% 
                                filter(model == .x) %>% 
                                ggplot() +
                                geom_boxplot(aes(x = study_yr, y = E1)) +
                                labs(x = "Year of study",
                                     y = "Pearson Residuals",
                                     title = .x)
                            })

year_vs_resid[1]
year_vs_resid[2]
year_vs_resid[3]
year_vs_resid[4]
year_vs_resid[5]



### v) Check for Pearson residuals for remaining spatial dependency ----

# Plot SPATIAL FIELDS
  #First need to get the posterior means and sds of the w's into one data frame
myspatialmodels = list("M.soF" = M.soF,
                       "M.wsF" = M.wsF)

myspatialnames <- unique(names(myspatialmodels))

w.df <- m.getws(myspatialmodels = myspatialmodels)
range(w.df$w.pm) #-2.044170  2.831429 #to get on a scale that can compare between all models

#If want to zoom in on a few wetlands
xlim.multi <- c(556878, 580992)    #LMT, BHER, and SAR
ylim.multi <- c(5436383, 5444603)  #LMT, BHER, and SAR


#Plot "M.soF" with 'm.gg2PlotField' function
w.pm <- M.soF$summary.random$w$mean

m.gg2PlotField(field = w.pm,
               mesh = mesh.f,
               ColorScheme = 1,                #You can use the values 1 to 4
               #xlim = range(mesh.f$loc[ ,1]), #Either need to add xlim and ylim, OR add a ContourMap
               #ylim = range(mesh.f$loc[ ,2]),
               xlim = xlim.multi,
               ylim = ylim.multi,
               scalemin = min(w.df$w.pm),    #Provide an input for min value for colorbar, to better compare between different plots
               scalemax = max(w.df$w.pm),    #Provide an input for max value for colorbar, to better compare between different plots
               Xlabel = "Easting",
               Ylabel = "Northing",
               #AddPoints = T,                 #Default is 'F'
               MyTitle = "M.soF - Posterior mean values w")


# Plot "M.wsF" with 'm.gg2PlotField' function
w.pm <- M.wsF$summary.random$w$mean

m.gg2PlotField(field = w.pm,
               mesh = mesh.f,
               ColorScheme = 1,                #You can use the values 1 to 4
               xlim = range(mesh.f$loc[ ,1]), #Either need to add xlim and ylim, OR add a ContourMap
               ylim = range(mesh.f$loc[ ,2]),
               #xlim = xlim.multi,
               #ylim = ylim.multi,
               scalemin = min(w.df$w.pm),    #Provide an input for min value for colorbar, to better compare between different plots
               scalemax = max(w.df$w.pm),    #Provide an input for max value for colorbar, to better compare between different plots
               Xlabel = "Easting",
               Ylabel = "Northing",
               #AddPoints = T,                 #Default is 'F'
               MyTitle = "M.wsF - Posterior mean values w")


# Plot residuals vs. wetland location for all the models
spatial_resids <- purrr::map(.x = myfinalnames,
                             .f = ~{
                                  #Need to make an average residual per wetland
                                  t <- myf.resids %>% 
                                        filter(model == .x) %>% 
                                        group_by(wtld) %>%  #' add 'study_yr' covariate to check for potential of temporal differences
                                        summarize(res = mean(E1),
                                                  east = mean(east),
                                                  north = mean(north),
                                                  Sign = as.factor(ifelse(res > 0, 1, 2))) %>% 
                                        mutate(MyCex = abs(res) / max(abs(res)))  #Note this is scaled, not just the absolute value...
                                 #Define colors for plotting      
                                 res.cols <- c("2" = "darkred",
                                               "1" = "seagreen")
                                 #Make the plot
                                 ggplot(t, aes(x = east, y = north, size = MyCex, fill = Sign)) +
                                   geom_point(alpha = 0.7, pch = 21) +
                                   theme_light() +
                                   facet_grid(study_yr ~.) +    #Use if want ST distribution
                                   scale_fill_manual(values = res.cols, labels = c("Neg", "Pos")) +
                                   labs(x = "Easting", y = "Northing",
                                        title = paste(.x, "spatial distribution of residuals", sep = ": "))
                            })

spatial_resids[1]
spatial_resids[2]
spatial_resids[3]
spatial_resids[4]
spatial_resids[5]



