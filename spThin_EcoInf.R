## =============================================================================================###
##  Spatial Thinning of Species Occurrence Records (SORs)
##  Script name: spThin_EcoInf.R
##
##  Authors: Varos Petrosyan, Fedor Osipov
##
##  Purpose: Reduce spatial autocorrelation in species occurrence records
### ============================================================================================###

# Description:
# This R script uses the spThin R package (Aiello-Lammens et al., 2015) 
# to generate seven thinned SOR subsets, each with a minimum inter-record distance 
# ranging from 10 to 70 km (in 10 km increments). 
# The script is executed separately for each geographic region—
# North America, Europe, Northern Eurasia, Armenia and adjacent regions, and Southeast Asia—
# and for each species individually. 
# This approach accounts for uneven sampling intensity across regions and reduces sampling bias.


# Input data for Carassius gibelio:
# Carassius_gibelio_N_America.csv            - North America
# Carassius_gibelio_Europe.csv               - Europe
# Carassius_gibelio_Northern_Eurasia.csv     - Northern Eurasia
# Carassius_gibelio_Armenia_and_adjacent_regions.csv - Armenia and adjacent regions
# Carassius_gibelio_Southeast_Asia.csv      - Southeast Asia

# Input data for Pseudorasbora parva:
# Pseudorasbora_parva_Europe.csv             - Europe
# Pseudorasbora_parva_Northern_Eurasia.csv   - Northern Eurasia
# Pseudorasbora_parva_Armenia_and_adjacent_regions.csv - Armenia and adjacent regions
# Pseudorasbora_parva_Southeast_Asia.csv    - Southeast Asia


## Input file location:
## Working_Directory/Source_Data/ALL_SORS/CSV_Full/

## All .csv files contain three columns: Species, Dd_long, and Dd_lat.
## After selecting a species, the R script generates seven files, each representing a dataset 
## with a specified minimum distance between occurrence records (10–70 km).
### ===============================================================================================###

# Output:
# Working_Directory/Source_Data/ALL_SORS/CSV_Full/
# All .csv files contain three columns: Species, Dd_long, and Dd_lat.
# After selecting a species, the script generates seven files, 
# each representing a dataset with a specified minimum distance between occurrence records (10–70 km).

# Output file location:
# Working_Directory/Results/Reduced_SORS/

# Next steps:
# It is recommended to evaluate the resulting subsets using the Average Nearest Neighbor Index (ANNI),
# for example with ArcGIS Desktop 10.8.1 (ESRI, 2020).

# Example:
# The results of record thinning (SORs_Full_csv.zip) using this script for C. gibelio and P. parva,
# including both invasive and native ranges, are provided in the SORs_Reduced_csv.zip archive.


library("spThin")

PathToWorkingDirectory <-"Working_Directory" 
setwd(PathToWorkingDirectory)
PathToSources<-paste0(getwd(),"/","Source_Data/ALL_SORS/CSV_Full")
PartToResults <-paste0(getwd(),"/","Results/Reduced_SORS")


## Select the species SORs for Thinning (optional) - Carassius gibelio

Thin_FileName<-file.path(PathToSources, "Carassius_gibelio_N_America.csv"); Region<-"N_America"; speciesName<-"Carassius_gibelio"
Thin_FileName<-file.path(PathToSources, "Carassius_gibelio_Europe.csv"); Region<-"Europe"; speciesName<-"Carassius_gibelio"
Thin_FileName<-file.path(PathToSources, "Carassius_gibelio_Northern_Eurasia.csv"); Region<-"N_Eurasia"; speciesName<-"Carassius_gibelio"
Thin_FileName<-file.path(PathToSources, "Carassius_gibelio_Armenia_and_adjacent_regions.csv"); Region<-"Armenia_AdjReg"; speciesName<-"Carassius_gibelio"
Thin_FileName<-file.path(PathToSources, "Carassius_gibelio_Southeast_Asia.csv"); Region<-"SE_Asia"; speciesName<-"Carassius_gibelio"

## Select the species SORs for Thinning (optional) - Pseudorasbora parva
Thin_FileName<-file.path(PathToSources, "Pseudorasbora_parva_Europe.csv"); Region<-"Europe";speciesName<-"Pseudorasbora_parva"
Thin_FileName<-file.path(PathToSources, "Pseudorasbora_parva_Northern_Eurasia.csv"); Region<-"N_Eurasia";speciesName<-"Pseudorasbora_parva"
Thin_FileName<-file.path(PathToSources, "Pseudorasbora_parva_Armenia_and_adjacent_regions.csv"); Region<-"Armenia_AdjReg"; speciesName<-"Pseudorasbora_parva"
Thin_FileName<-file.path(PathToSources, "Pseudorasbora_parva_Southeast_Asia.csv"); Region<-"SE_Asia"; speciesName<-"Pseudorasbora_parva"

ThinDataSet <- read.csv(Thin_FileName)
ThinDataSet$Species <- as.character(ThinDataSet$Species)

for (Distance in seq(from = 10, to = 70, by = 10)) {
  print(Distance )

  outputFile <- paste0(speciesName,"_",Region,"_",Distance,'km', '_')

Thin_dataset_full_20km <- thin( loc.data = ThinDataSet, 
                                     lat.col = "Dd_lat", long.col = "Dd_long", 
                                     spec.col = "Species", 
                                     thin.par = Distance, reps = 1, 
                                     locs.thinned.list.return = TRUE, 
                                     write.files = TRUE, 
                                     max.files = 1, 
                                     out.dir =PartToResults, out.base = outputFile, 
                                     write.log.file = FALSE,
                                     log.file = "Thinned_full_log_file.txt")
}


