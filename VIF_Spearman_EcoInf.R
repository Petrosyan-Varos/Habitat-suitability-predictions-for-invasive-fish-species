## =========================================================================================##
##============  Assessment of multicollinearity of predictor variables  ====================##
# 
## Script name: VIF_Spearman_EcoInf.R
##
## Authors: Varos Petrosyan, Fedor Osipov
##
## Purpose: Selection of predictor variables for constructing SDMs
## 
## Description: The VIF_Spearman_EcoInf.R script estimates Spearman’s pairwise rank correlation coefficients 
## and assesses multicollinearity of predictor variables for constructing SDMs. 
## Variables with an absolute Spearman’s pairwise rank correlation coefficient greater than 0.70 
## are excluded from the model. 
## Multicollinearity of the remaining variables is further assessed using the 
## Variance Inflation Factor (VIF) implemented in the R package 'usdm'. 
## A variable is considered multicollinear and excluded from the model if VIF > 10 
## (Chatterjee and Hadi, 2006).
## Example: This script demonstrates that the selected variables meet the chosen 
## VIF threshold for modeling.
## Requirements:
## R packages: raster, usdm, ENMTools, corrplot
## Reference:
## Chatterjee, S., Hadi, A. S. 2006. Regression analysis by example, fourth ed., New Jersey
### ==============================================================================================###

library(raster)
library(usdm)
library(ENMTools) 
library(corrplot)

# Setting up working directories 
PathToWorkingDirectory <- "Working_Directory"  
setwd(PathToWorkingDirectory)
PathToSources<-paste0(getwd(),"/","Source_Data/Predictors/")
PartToResults <-paste0(getwd(),"/","Results/Predictors")

# Loading predictor variables 
# All predictors must have the same resolution and spatial extent

# Input predictors:    
#  Bio_03.asc  - Isothermality (BIO2/BIO7 ?? 100)
#  Bio_08.asc  - Mean temperature of the wettest quarter (°C)
#  Bio_12.asc  - Annual precipitation (mm)
#  Bio_18.asc  - Precipitation of the warmest quarter (mm)
#  Windmean.asc - Mean wind speed (m/s)
#  Elevation.asc - Altitude (m)
#  Slope.asc - Slope inclination (°)
#  WaterPoly.asc - Distance to water bodies (lakes, reservoirs, ponds) (m)
#  WaterLine.asc - Distance to rivers or streams (m)

# Input data 
files <- list.files(PathToSources, pattern = "\\.asc$", full.names = TRUE)
stack_clim <- stack(files)   # create a raster stack
names(stack_clim) # Let's check the layer names

# To reduce computational load, select 10,000 random points (can be adjusted)
# Extracting random points ===
set.seed(123)
sample_points <- sampleRandom(stack_clim, size = 10000, na.rm = TRUE, sp = FALSE)

# === Spearman correlation ===
cor_spearman <- cor(sample_points, method = "spearman", use = "complete.obs")
round (cor_spearman, 2)

# Creating a Spearman correlation matrix 

matr<-as.matrix(cor_spearman)
corrplot(matr, method="color",type = c("upper"),addCoef.col=" darkgreen", 
         addgrid.col = "gray33", tl.col = "black")

# === VIF assessment ===
# Specifying the VIF threshold value (10 is recommended) – Chatterjee and Hadi, 2006 

vif_results <- vifstep(sample_points, th = 10)  
vif_results@results     # all values
vif_results@excluded    # excluded variables
vif_results@kept        # variables that remain

# === Saving results ===

## Spearman’s correlation matrix heatmap of predictor variables
CorrPlotFile <- file.path(PartToResults, "CorrSpearmanPlot.png")
png(CorrPlotFile , width = 800, height = 600)
corrplot(matr, method="color",type = c("upper"),addCoef.col=" darkgreen", 
         addgrid.col = "gray33", tl.col = "black")
dev.off()

## Spearman’s correlation matrix of predictor variables
cor_spearman_file <- file.path(PartToResults, "CorrSpearmanMatrix.csv")
write.csv(cor_spearman, "correlation_spearman.csv", row.names = TRUE)

## List of excluded variables
Vif_Value_File <- file.path(PartToResults, "vif_results.csv")
write.csv(vif_results@results, Vif_Value_File, row.names = FALSE)



