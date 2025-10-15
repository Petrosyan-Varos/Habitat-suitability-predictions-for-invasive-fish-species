##====================================================================================##
##  Script name: CBI_EcoInf.R
##  Authors: Varos Petrosyan, Fedor Osipov

##  Purpose: Bootstrap assessment of the Continuous Boyce Index (CBI)
##  
##  Purpose: Assessing model performance using Continuous Boyce's index
##
##  Description:
##    This script performs a bootstrap-based evaluation of model performance
##    using the Continuous Boyce Index (CBI). It estimates the mean and confidence
##    intervals of CBI values from habitat suitability raster data and species
##    occurrence records (SORs). The script randomly generates background points,
##    computes CBI values for multiple bootstrap replicates, and saves both
##    summary statistics and individual replicate results.
##
##  Input data:
##    pred_rast — raster layer of predicted habitat suitability (*.tif format)
##    pres_pts  — table of species occurrence records (SORs, *.csv format)
##
##  Output files:
##    CBI_confIntervals.csv — summary table including:
##         - number of replicates and background points
##         - estimation method
##         - mean CBI
##         - lower and upper confidence interval limits
##    CBI_bootstrap_results.csv — individual bootstrap replicates
##
##  File structure:
##    Input files location:  Working_Directory/Source_Data/
##    Output files location: Working_Directory/Results/
##
##  Requirements:
##    - R version = 4.2
##    - Packages: raster, ecospat, sp
##
##  Example usage:
##    source("CBI_EcoInf.R")
##
##  Authors: [Petrosyan Varos, Fedor Osipov ]
##  Affiliation: [A.N. Severtsov Institute of Ecology and Evolution of the Russian Academy of Sciences, Moscow, Russia]
##  Date: [October 14, 2025]
##====================================================================================##


library(raster)
library(ecospat)
library(sp)

# Setting up working directories 
PathToWorkingDirectory <-"Working_Directory" 
setwd(PathToWorkingDirectory)
PathToSources<-paste0(getwd(),"/","Source_Data")
PartToResults <-paste0(getwd(),"/","Results")

##=== Input data ===##
# Specify the raster file containing predicted habitat suitability
SDM_raster_FileName <- file.path(PathToSources, "Maxent_P_parva_SDM.tif") 

##=== Read the raster layer  
pred_rast <- raster(SDM_raster_FileName)

##=== Specify the file containing species occurrence records
presence_points_File<-file.path(PathToSources, "P_Parva_Armenia_and_Adjreg.csv") #*

##=== Read the file that contains two columns: lon and lat
pres_pts <- read.csv(presence_points_File)

##=== Set the SpatialPoints projection to match the raster projection
coordinates(pres_pts) <- ~lon + lat
crs(pres_pts) <- crs(pred_rast)


##=== Function for estimating the CBI index from a single run

boyce_once <- function(pred_rast, pres_pts, n_bg = 1000) {
  
  ##=== The suitability  values is extracted at the points of presence
  obs_val <- raster::extract(pred_rast, pres_pts)
  obs_val <- obs_val[!is.na(obs_val)]  # 
  
  ##=== All background pixels in the raster
  fit_all <- raster::getValues(pred_rast)
  fit_all <- fit_all[!is.na(fit_all)]
  
  ##=== It is suggested to additionally select 1000 random background points
  ##===  (not necessary, but increases stability), this is important for this case 
  
  if (length(fit_all) > n_bg) {
    fit_all <- sample(fit_all, n_bg)
  }
  ##=== CBI assessment 
  boy <- ecospat.boyce(fit = fit_all, obs = obs_val, nclass = 0, PEplot = FALSE)
  return(boy$Spearman.cor)
}


set.seed(123)  
n_reps <- 10  # Number of repetitions: 100

##=== Random generation of background points and CBI estimation 
boyce_values <- replicate(n_reps, boyce_once(pred_rast, pres_pts, n_bg = 1000))

##=== Bootstrap estimation of the mean CBI
mean_boyce <- mean(boyce_values, na.rm = TRUE)
cat("CBI mean =", round(mean_boyce, 3), "\n")

##=== Confidence interval estimation - 0.025, 0.975
ci <- quantile(boyce_values, probs = c(0.025, 0.975), na.rm = TRUE)
cat("95% Confidence interval: [", round(ci[1], 3), ";", round(ci[2], 3), "]\n")

##=== Creating a data frame containing the mean and CBI confidence interval
CBI_Val <-data.frame(
  Repitation = n_reps,
  BackgrPoints = n_bg,
  Methods= c("bootstrap"),
  CBI_mean = mean_boyce ,
  CBILow_Lim = round(ci[[1]], 3),
  CBIUpper_Lim= round(ci[[2]], 3)
  )


##  Saving results
# Contains:
#   - number of replicates and background points
#   - name of the estimation method
#   - mean CBI
#   - lower and upper confidence interval limits

CBI_ConfInt<-file.path(PartToResults, "CBI_confIntervals.csv")
write.csv(CBI_Val, CBI_ConfInt, row.names = FALSE)

boyce_df <- data.frame(run = 1:n_reps, Boyce = boyce_values)
CBI_bootstrap_FileName <- file.path(PartToResults, "CBI_bootstrap_results.csv")

# Individual bootstrap replicates
write.csv(boyce_df, CBI_bootstrap_FileName, row.names = FALSE)

