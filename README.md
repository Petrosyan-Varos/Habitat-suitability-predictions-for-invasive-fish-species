This document describes the R scripts used for processing species occurrence records, reducing spatial autocorrelation, evaluating and selecting predictor variables for species distribution modeling (SDMs), analyzing climatic niches using PCA-env, and conducting bootstrap assessments of the Continuous Boyce Index (CBI).
1. Spatial Thinning of Species Occurrence Records (SORs)
Script name: spThin_EcoInf.R
Authors: Varos Petrosyan, Fedor Osipov
Purpose: Reduce spatial autocorrelation in species occurrence records.
Description: This R script uses the spThin R package (Aiello-Lammens et al., 2015) to generate seven thinned SOR subsets, each with a minimum inter-record distance ranging from 10 to 70 km (in 10 km increments). The script is executed separately for each geographic region—North America, Europe, Northern Eurasia, Armenia and adjacent regions, and Southeast Asia—and for each species individually. This approach accounts for uneven sampling intensity across regions and reduces sampling bias.
Input data for Carassius gibelio:
•	Carassius_gibelio_N_America.csv — North America
•	Carassius_gibelio_Europe.csv — Europe
•	Carassius_gibelio_Northern_Eurasia.csv — Northern Eurasia
•	Carassius_gibelio_Armenia_and_adjacent_regions.csv — Armenia and adjacent regions
•	Carassius_gibelio_Southeast_Asia.csv — Southeast Asia
Input data for Pseudorasbora parva:
•	Pseudorasbora_parva_Europe.csv — Europe
•	Pseudorasbora_parva_Northern_Eurasia.csv — Northern Eurasia
•	Pseudorasbora_parva_Armenia_and_adjacent_regions.csv — Armenia and adjacent regions
•	Pseudorasbora_parva_Southeast_Asia.csv — Southeast Asia
Output:
•	Location: Working_Directory/Source_Data/ALL_SORS/CSV_Full/
•	All *.csv files contain three columns: Species, Dd_long, and Dd_lat.
•	After selecting a species, the script generates seven files, each representing a dataset with a specified minimum distance between occurrence records (10–70 km).
Output file location:
Working_Directory/Results/Reduced_SORS/
Next steps:
It is recommended to evaluate the resulting subsets using the Average Nearest Neighbor Index (ANNI), for example with ArcGIS Desktop 10.8.1 (ESRI, 2020).
Example:
The results of record thinning (SORs_Full_csv.zip) using this script for C. gibelio and P. parva, including both invasive and native ranges, are provided in the SORs_Reduced_csv.zip archive.
 
2. Assesment of multicollinearity of predictor variables
Script name: VIF_Sperman_EcoInf.R
Authors: Varos Petrosyan, Fedor Osipov
Purpose: Selection of predictor variables for constructing SDMs
Description: The VIF_Sperman_EcoInf.R script estimates Spearman's pairwise rank correlation coefficients and multicollinearity of variables for constructing SDMs. Variables with an absolute Spearman's pairwise rank correlation coefficient greater than 0.70 were excluded from the model. Multicollinearity of the remaining variables was further assessed using the variance inflation factor (VIF) in the R usdm package. A variable was considered multicollinear and excluded from the model if VIF > 10 (Chatterjee and Hadi, 2006).

Input:
Working_Directory/Source_Data/Predictors/
Bio_03.asc - Isothermality (BIO2/BIO7 (x100))
Bio_08.asc - Mean temperature in wettest quarter (°C)
Bio12.asc - Annual precipitation (mm)
Bio_18.asc - Precipitation in warmest quarter (mm)
Windmean.asc - Wind speed (m/s)
Elevation.asc - Altitude, m
Slope.asc - Inclination angles, °
WaterPoly.asc - Distance to water polygons (lakes, reserves, ponds), m
WaterLine.asc - Distance to river, stream, m)

Output:
Working_Directory/Results/Predictors/
CorrSpearmanPlot.png- Spearman’s correlation matrix heatmap of predictor variables
CorrSpearmanMatrix.csv - Spearman’s correlation matrix of predictor variables
Vif_results.csv - excluded variables

Requirements:
R packages: raster, usdm, ENMTools, corrplot

Example: This script demonstrates that the selected variables meet the chosen VIF threshold for modeling.

Reference:
Chatterjee, S., Hadi, A. S. 2006. Regression analysis by example, fourth ed., New Jersey
 
3. Estimates of climatic niche metrics for the species under the ECOUE concept 
Script name: ECOUE_species_EcoInf.R
Authors: Varos Petrosyan, Fedor Osipov
Purpose: Estimation of seven climatic niche metrics under the ECOUE framework
Description: This script calculates seven climatic niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim) for invasive species under the extended ECOUE (Ecological Niche Overlap and Unfilled–Expansion) framework. The ECOUE approach integrates PCA-env analysis to quantify and compare the realized climatic niches between native and invasive ranges of a species.
The extended ECOUE framework (Liu et al., 2020; Petrosyan et al., 2023) expands traditional niche overlap metrics by including additional niche breadth indices:
 - Sm (stability index): proportion of the invasive niche overlapping with the native niche.
 - Em (expansion index): proportion of the invasive niche not occupied in the native range.
 - Um (unfilled index): proportion of the native niche unoccupied in the invasive range.
 - Bn (native niche breadth) = Sm + Um
 - Bi (invasive niche breadth) = Sm + Em
 - BR (breadth ratio) = ln(Bn / Bi)
 - Sim (niche similarity) = 2Sm / (Bn + Bi)

Input:
Working_Directory/Source_Data/
Pseudo_parva_BackgrNative.csv : P. parva SORs in the native range and pseudo-absence points
Pseudo_parva_BackgrInvasive.csv : P. parva SORs in the invasive range and pseudo-absence points
Carassius_gibelio_BackgrNative.csv : C. gibelio SORs in the native range and pseudo-absence points
Carassius_gibelio_BackgrInvasive.csv : C. gibelio SORs in the invasive range and pseudo-absence points
Each input file must include the following columns:##    OBJECTID	- Record ID
OBJECTID   – Record ID
Species    – Species name
Dd_long    – Longitude in decimal degrees
Dd_lat     – Latitude in decimal degrees
Bio_03     – Isothermality (BIO2/BIO7 × 100)
Bio_08     – Mean temperature of the wettest quarter (°C)
Bio_12     – Annual precipitation (mm)
Bio_18     – Precipitation of the warmest quarter (mm)
Elevation  – Altitude (m)
Slope      – Slope angle (°)
WindMean   – Mean wind speed (m/s)
Spec_Occ   – 1 = presence point; 0 = background point
Region     – Invasive region name or pseudo-absence area name

 Output:
 Working_Directory/Results/
 - P_parva_NicheComparison_AllRegions.csv : Estimates of the seven niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim) for P. parva under the ECOUE framework.

Example:
The provided dataset demonstrates the ECOUE analysis using P. parva data   (Pseudo_parva_BackgrNative.csv and Pseudo_parva_BackgrInvasive.csv).
Additional files for C. gibelio are included for running the script on another species.
Requirements:
R packages: ecospat, raster, sp

Reference:
Liu, C. et al. (2020). https://doi.org/10.1073/pnas.2004289117.
Petrosyan, V. et al. (2023). https://doi.org/10.1134/S106235902360126X.

 
4. Bootstrap assessment of the Continuous Boyce Index (CBI) 
Script name: CBI_EcoInf.R
Authors: Varos Petrosyan, Fedor Osipov
Purpose: Assessing model performance using Continuous Boyce's index
Description: This script performs a bootstrap-based evaluation of model performance using the Continuous Boyce Index (CBI). It estimates the mean and confidence intervals of CBI values from habitat suitability raster data and species occurrence records (SORs). The script randomly generates background points, computes CBI values for multiple bootstrap replicates, and saves both  summary statistics and individual replicate results.
Input files location:
Working_Directory/Source_Data/
Input files (see SDMs_SORs.zip):
This archive contains four files:
1.	Maxent_C_gibelio_SDM.tif – predicted suitable habitats for C. gibelio in Armenia and adjacent regions
2.	C_gibelio_SORs_for_ArmeniaAdReg.csv – species occurrence records of C. gibelio in Armenia and adjacent regions
3.	Maxent_P_parva_SDM.tif – predicted suitable habitats for P. parva in Armenia and adjacent regions
4.	P_parva_SORs_for_ArmeniaAdReg.csv – species occurrence records of P. parva in Armenia and adjacent regions
Output files location:
Working_Directory/Results/
Output files:
•	CBI_confIntervals.csv — contains the number of replicates and background points, the estimation method, the average CBI, and the lower and upper limits of the confidence interval.
•	CBI_bootstrap_results.csv — contains individual bootstrap replicates.
This example uses C. gibelio data (Maxent_C_gibelio_SDM.tif and C_gibelio_SORs_for_ArmeniaAdReg.csv) to demonstrate how the script works. 
However, the archive also includes data (Maxent_P_parva_SDM.tif and P_parva_SORs_for_ArmeniaAdReg.csv) for running the script on a different species.
 
List of files
1. Spatial Thinning of Species Occurrence Records (SORs)
1.1. spThin_EcoInf.R
1.2. SORs_Full_csv.zip
1.3. SORs_Reduced_csv.zip
2. Assesment of multicollinearity of predictor variables
2.1.VIF_Sperman_EcoInf.R
2.2. Predictors.zip
3. Estimates of climatic niche metrics for the species under the ECOUE concept
3.1. ECOUE_species_EcoInf.R
3.2. Species_NatInv-PsAb_points.zip
4. Bootstrap assessment of the Continuous Boyce Index (CBI)
4.1. CBI_EcoInf.R
4.2. SDMs_SORs.zip


## License
This project is licensed under the MIT License — see the [LICENSE](./LICENSE) file for details.
Copyright (c) [2025] [Varos Petrosyan and Fedor Osipov]
For questions or collaboration inquiries, please contact:
Varos Petrosyan
Email: vgpetrosyan@gmail.com
