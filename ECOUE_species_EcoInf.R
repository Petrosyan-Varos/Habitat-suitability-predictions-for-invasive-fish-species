##====================================================================================##
##  Script name: ECOUE_species_EcoInf.R
##  Author: Varos Petrosyan, Fedor Osipov
##  Purpose: Estimation of seven climatic niche metrics under the ECOUE framework
##
##  Description:
##    This script calculates seven climatic niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim)
##    for invasive species under the extended ECOUE (Ecological Niche Overlap and
##    Unfilled–Expansion) framework. The ECOUE approach integrates PCA-env analysis
##    to quantify and compare the realized climatic niches between native and invasive
##    ranges of a species.
##
##    The extended ECOUE framework (Liu et al., 2020; Petrosyan et al., 2023b) expands
##    traditional niche overlap metrics by including additional niche breadth indices:
##      - Sm (stability index): proportion of the invasive niche overlapping with
##        the native niche.
##      - Em (expansion index): proportion of the invasive niche not occupied
##        in the native range.
##      - Um (unfilled index): proportion of the native niche unoccupied in
##        the invasive range.
##      - Bn (native niche breadth) = Sm + Um
##      - Bi (invasive niche breadth) = Sm + Em
##      - BR (breadth ratio) = ln(Bn / Bi)
##      - Sim (niche similarity) = 2Sm / (Bn + Bi)
##
##  Input:
##    Working_Directory/Source_Data/
##      - Pseudo_parva_BackgrNative.csv        : P. parva SORs in the native range and pseudo-absence points
##      - Pseudo_parva_BackgrInvasive.csv      : P. parva SORs in the invasive range and pseudo-absence points
##      - Carassius_gibelio_BackgrNative.csv   : C. gibelio SORs in the native range and pseudo-absence points
##      - Carassius_gibelio_BackgrInvasive.csv : C. gibelio SORs in the invasive range and pseudo-absence points
##
##    Each input file must include the following columns:
##      OBJECTID   – Record ID
##      Species    – Species name
##      Dd_long    – Longitude in decimal degrees
##      Dd_lat     – Latitude in decimal degrees
##      Bio_03     – Isothermality (BIO2/BIO7 × 100)
##      Bio_08     – Mean temperature of the wettest quarter (°C)
##      Bio_12     – Annual precipitation (mm)
##      Bio_18     – Precipitation of the warmest quarter (mm)
##      Elevation  – Altitude (m)
##      Slope      – Slope angle (°)
##      WindMean   – Mean wind speed (m/s)
##      Spec_Occ   – 1 = presence point; 0 = background point
##      Region     – Invasive region name or pseudo-absence area name
##
##  Output:
##    Working_Directory/Results/
##      - P_parva_NicheComparison_AllRegions.csv : Estimates of the seven niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim)
##                                                 for P. parva under the ECOUE framework.
##
##  Example:
##    The provided dataset demonstrates the ECOUE analysis using P. parva data
##    (Pseudo_parva_BackgrNative.csv and Pseudo_parva_BackgrInvasive.csv).
##    Additional files for C. gibelio are included for running the script on another species.
##
##  Requirements:
##    R packages: ecospat, raster, sp
##
##  Reference:
##    Liu, C. et al. (2020). https://doi.org/10.1073/pnas.2004289117.
##    Petrosyan, V. et al. (2023b). https://doi.org/10.1134/S106235902360126X.
##
##====================================================================================##

library("ecospat") 
library(raster)

# Setting up working directories 
PathToWorkingDirectory <-"E:/Global_Change_2024" 
setwd(PathToWorkingDirectory)
PathToSources<-paste0(getwd(),"/","Source_Data")
PartToResults <-paste0(getwd(),"/","Results")

## List of different invasive regions 

  regions <- list(
  "N_Eurasia" = c("Russia"),
  "Europe"    = c("Europe"),
  "Armenia"   = c("Armenia")
  )
  
  ## Files for PCA-env analysis of P. parva using the ECOUE concept  
  Species<-"P_parva" #  
  Native_Range_FileName <- file.path(PathToSources, "Pseudo_parva_BackgrNative.csv") 
  Inv_Range_FileName <- file.path(PathToSources, "Pseudo_parva_BackgrInvasive.csv")
  
  ## Read data – species occurrence records and background points in the native range 
natData <- read.csv(Native_Range_FileName,fill = TRUE, header = TRUE  )

##    Each input file must include the following columns:
##      OBJECTID   – Record ID
##      Species    – Species name
##      Dd_long    – Longitude in decimal degrees
##      Dd_lat     – Latitude in decimal degrees
##      Bio_03     – Isothermality (BIO2/BIO7 × 100)
##      Bio_08     – Mean temperature of the wettest quarter (°C)
##      Bio_12     – Annual precipitation (mm)
##      Bio_18     – Precipitation of the warmest quarter (mm)
##      Elevation  – Altitude (m)
##      Slope      – Slope angle (°)
##      WindMean   – Mean wind speed (m/s)
##      Spec_Occ   – 1 = presence point; 0 = background point
##      Region     – Invasive region name or pseudo-absence area 


### Select predictor variables from the tables for the native range
selCol<-c(5:13)

## Selection from the table of species occurrence records using 
## Spec_Occ – column    
presence<-1
# Read Presence Data
nat1<-natData[which(natData[,12]==presence),selCol]

# Read pseudo-absence records
presence<-0
nat2<-natData[which(natData[,12]==presence),selCol]
## 
nat0<- rbind(nat1,nat2)

# Delete NA
nat <-nat0[complete.cases(nat0), ] # Remove missing values (NA)

### Load data – species occurrence records and background points in the invasive range 
invDat<-read.csv(Inv_Range_FileName , fill = TRUE, header = TRUE)

results_all <- data.frame()

for (region_name in names(regions)) {
  region_data <- regions[[region_name]]
   print(region_data)
   
   
# Select records  From the Species Table
PointType<-1
inv1<-invDat[which(invDat[,13]==region_data & invDat[,12]==PointType),selCol] 
presence<-0
inv2<-invDat[which(invDat[,12]==presence),selCol]

inv0<- rbind(inv1,inv2)
# Delete NA
inv<-inv0[complete.cases(inv0), ] # Remove missing values (NA)# Remove missing values (NA)

## Niche quantification and comparison
### PCA-ENVIRONMENT
# The PCA is calibrated on all sites within the study area


ColumnselCol<-c(1,2,3,4,7)  ##

## PCA-env 
par(mfrow=c(1,1))
pca.env<-dudi.pca(rbind(nat,inv)[,ColumnselCol],scannf=FALSE,nf=2)

# Plot	Variables	Contribution	with	ecospat.plot.contrib()
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

# Predict the scores on the axes
scores.globclim <-pca.env$li

#	PCA	scores	for	the	species	native range
scores.sp.nat<-suprow(pca.env,nat[which(nat[,8]==1),ColumnselCol])$li

#	PCA	scores	for	the	species	invasive range 
scores.sp.inv<-suprow(pca.env,inv[which(inv[,8]==1),ColumnselCol])$li


#	PCA	scores	for	the	whole	native	area
scores.clim.nat<-suprow(pca.env,nat[,ColumnselCol])$li 


#	PCA	scores	for	the	whole	invasive area
scores.clim.inv<-suprow(pca.env,inv[,ColumnselCol])$li 


#	Calculate	the	Occurrence	Densities	Grid	with	ecospat.grid.clim.dyn()
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                     glob1=scores.clim.nat, sp=scores.sp.nat, R=100, th.sp=0)

grid.clim.inv<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv, sp=scores.sp.inv, R=100, th.sp=0)

#	Calculate	Niche	Overlap	with	
overlap<-ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=TRUE)

# Evaluation of E,S, U metrics
Metrics<-ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection=0.1)

Expantion<-Metrics$dynamic.index.w [["expansion"]]
Stability <-Metrics$dynamic.index.w [["stability"]]
Unfilling <-Metrics$dynamic.index.w [["unfilling"]]

# Evaluation of an extended set of metrics (Em, Sm, Um, Bn, Bi, BR, Sim) within the ECOUE concept
Ut<-(Stability/(1-Unfilling))/Unfilling
scal<-1+Ut
ModExpantion<-Expantion/(Expantion+Stability+Ut)
ModStability<-Stability/(Expantion+Stability+Ut)
ModUnfilling<-Ut/(Expantion+Stability+Ut)
Bnval<-(ModStability+ModUnfilling)
Bival<-(ModExpantion+ModStability)
LnBRval<-log(Bnval/Bival)
SimIndval<-2*ModStability/(Bnval+Bival)

  new_row <-data.frame(
  Region = region_data,
  Em = ModExpantion,
  Sm = ModStability,
  Um = ModUnfilling,
  Bn=Bnval,
  Bi=Bival,
  LnBR=LnBRval,
  SimInd=SimIndval
   )
  results_all<- rbind(results_all, new_row)
} # for 

Result_Table_FileName <- paste0(PartToResults,"/","NicheComparison_AllRegions.csv")

write.csv(results_all, Result_Table_FileName, row.names = FALSE)

cat("\n=== All results are saved in NicheComparison_AllRegions NicheComparison_AllRegions.csv ===\n")
print(results_all)

