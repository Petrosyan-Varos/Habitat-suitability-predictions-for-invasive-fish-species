# This document describes the R scripts used for:

## 1. Spatial Thinning of Species Occurrence Records (SORs). Script name: spThin_EcoInf.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Reduce spatial autocorrelation in species occurrence records.  
Description: This R script uses the spThin R package (Aiello-Lammens et al., 2015) to generate seven thinned SOR subsets, each with a minimum inter-record distance ranging from 10 to 70 km (in 10 km increments). The script is executed separately for each geographic region—North America, Europe, Northern Eurasia, Armenia and adjacent regions, and Southeast Asia—and for each species individually. This approach accounts for uneven sampling intensity across regions and reduces sampling bias.  

Input data for Carassius gibelio:  
• Carassius_gibelio_N_America.csv — North America  
• Carassius_gibelio_Europe.csv — Europe  
• Carassius_gibelio_Northern_Eurasia.csv — Northern Eurasia  
• Carassius_gibelio_Armenia_and_adjacent_regions.csv — Armenia and adjacent regions  
• Carassius_gibelio_Southeast_Asia.csv — Southeast Asia  

Input data for Pseudorasbora parva:  
Pseudorasbora_parva_Europe.csv — Europe; Pseudorasbora_parva_Northern_Eurasia.csv — Northern Eurasia; Pseudorasbora_parva_Armenia_and_adjacent_regions.csv — Armenia and adjacent regions; Pseudorasbora_parva_Southeast_Asia.csv — Southeast Asia.  

Output Data:  
Location: Working_Directory/Source_Data/ALL_SORS/CSV_Full/  
All *.csv files contain three columns: Species, Dd_long, and Dd_lat. After selecting a species, the script generates seven files, each representing a dataset with a specified minimum distance between occurrence records (10–70 km).  

Output file location: Working_Directory/Results/Reduced_SORS/  

Next steps: It is recommended to evaluate the resulting subsets using the Average Nearest Neighbor Index (ANNI), for example with ArcGIS Desktop 10.8.1 (ESRI, 2020).  

Example: The results of record thinning (SORs_Full_csv.zip) using this script for C. gibelio and P. parva, including both invasive and native ranges, are provided in the SORs_Reduced_csv.zip archive.  

________________________________________

## 2. Assessment of multicollinearity of predictor variables. Script name: VIF_Spearman_EcoInf.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Selection of predictor variables for constructing SDMs  
Description: The VIF_Spearman_EcoInf.R script estimates Spearman's pairwise rank correlation coefficients and multicollinearity of variables for constructing SDMs. Variables with an absolute Spearman's pairwise rank correlation coefficient greater than 0.70 were excluded from the model. Multicollinearity of the remaining variables was further assessed using the variance inflation factor (VIF) in the R usdm package. A variable was considered multicollinear and excluded from the model if VIF > 10 (Chatterjee and Hadi, 2006).  

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

Requirements: R packages: raster, usdm, ENMTools, corrplot  

Example: This script demonstrates that the selected variables meet the chosen VIF threshold for modeling.  

Reference: Chatterjee, S., Hadi, A. S. 2006. Regression Analysis by Example, fourth ed., New Jersey.  

________________________________________

## 3. Estimates of climatic niche metrics for the species under the ECOUE concept. Script name: ECOUE_species_EcoInf.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Estimation of seven climatic niche metrics under the ECOUE framework  
Description: This script calculates seven climatic niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim) for invasive species under the extended ECOUE (Ecological Niche Overlap and Unfilled–Expansion) framework. The ECOUE approach integrates PCA-env analysis to quantify and compare the realized climatic niches between native and invasive ranges of a species.  

The extended ECOUE framework (Liu et al., 2020; Petrosyan et al., 2023) expands traditional niche overlap metrics by including additional niche breadth indices:  
• Sm (stability index): proportion of the invasive niche overlapping with the native niche.  
• Em (expansion index): proportion of the invasive niche not occupied in the native range.  
• Um (unfilled index): proportion of the native niche unoccupied in the invasive range.  
• Bn (native niche breadth) = Sm + Um  
• Bi (invasive niche breadth) = Sm + Em  
• BR (breadth ratio) = ln(Bn / Bi)  
• Sim (niche similarity) = 2Sm / (Bn + Bi)  

Input:  
Working_Directory/Source_Data/  
Pseudo_parva_BackgrNative.csv: P. parva SORs in the native range and pseudo-absence points  
Pseudo_parva_BackgrInvasive.csv : P. parva SORs in the invasive range and pseudo-absence points  
Carassius_gibelio_BackgrNative.csv : C. gibelio SORs in the native range and pseudo-absence points  
Carassius_gibelio_BackgrInvasive.csv : C. gibelio SORs in the invasive range and pseudo-absence points  

Each input file must include the following columns:## OBJECTID - Record ID OBJECTID – Record ID Species – Species name Dd_long – Longitude in decimal degrees Dd_lat – Latitude in decimal degrees Bio_03 – Isothermality (BIO2/BIO7 × 100) Bio_08 – Mean temperature of the wettest quarter (°C) Bio_12 – Annual precipitation (mm) Bio_18 – Precipitation of the warmest quarter (mm) Elevation – Altitude (m) Slope – Slope angle (°) WindMean – Mean wind speed (m/s) Spec_Occ – 1 = presence point; 0 = background point Region – Invasive region name or pseudo-absence area name  

Output:  
Working_Directory/Results/  
P_parva_NicheComparison_AllRegions.csv : Estimates of the seven niche metrics (Em, Sm, Um, Bn, Bi, BR, Sim) for P. parva under the ECOUE framework.  

Example: The provided dataset demonstrates the ECOUE analysis using P. parva data (Pseudo_parva_BackgrNative.csv and Pseudo_parva_BackgrInvasive.csv). Additional files for C. gibelio are included for running the script on another species.  

Requirements: R packages: ecospat, raster, sp  

References:  
Liu, C. et al. (2020). https://doi.org/10.1073/pnas.2004289117.  
Petrosyan, V. et al. (2023). https://doi.org/10.1134/S106235902360126X.  

________________________________________

## 4. Bootstrap assessment of the Continuous Boyce Index (CBI). Script name: CBI_EcoInf.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Assessing model performance using Continuous Boyce's index.  
Description: This script performs a bootstrap-based evaluation of model performance using the Continuous Boyce Index (CBI). It estimates the mean and confidence intervals of CBI values from habitat suitability raster data and species occurrence records (SORs). The script randomly generates background points, computes CBI values for multiple bootstrap replicates, and saves both summary statistics and individual replicate results.  

Input files location: Working_Directory/Source_Data/  
Input files (see SDMs_SORs.zip):  
This archive contains two files:  
1. Maxent_P_parva_SDM.tif – predicted suitable habitats for P. parva in Armenia and adjacent regions  
2. P_parva_SORs_for_ArmeniaAdReg.csv – species occurrence records of P. parva in Armenia and adjacent regions  

Output files location: Working_Directory/Results/  
Output files:  
• CBI_confIntervals.csv — contains the number of replicates and background points, the estimation method, the average CBI, and the lower and upper limits of the confidence interval.  
• CBI_bootstrap_results.csv — contains individual bootstrap replicates.  

Example: This example uses P. parva data (Maxent_P_parva_SDM.tif and P_parva_SORs_for_ArmeniaAdReg.csv) to demonstrate how the script works.  

________________________________________

## 5. Estimation of accessible area (M) radius for SDM. Script name: Calculate_M_Area_Radius.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Estimation of accessible area (M) radius for SDM.  
Description: This script estimates the radius of the accessible area (M) used for background point generation in SDM workflows. The radius is derived from nearest-neighbor distances between species occurrence points.  

Method:  
1. Convert geographic coordinates to kilometers  
2. Compute all pairwise distances  
3. Estimate nearest-neighbor distance per point  
4. Define M radius as: mean(NN) + SD(NN)  

Input data for Carassius gibelio: Carasus_Global_All.csv  
Input data for Pseudorasbora parva: Parva_Global_All.csv  

Output data:  
Mean nearest-neighbor distance (km),  
Standard deviation (km),  
Recommended M radius (km),  

Example: Estimation of accessible area (M) radius for SDM of  Carassius gibelio and Pseudorasbora parva using SORS for Carasus_Global_All.csv and Parva_Global_All.csv (Carasus_ Parva_global_Sors.zip) respectively. This example uses to demonstrate how the script works.  

______________________________________

## 6. MaxEnt Hyperparameter Optimization using biomod2. Script name: Maxent_parameters_using_Biomod2_tuning.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Determining the optimal MaxEnt parameters using BIOMOD_tuning functions.  
Description: Ranking of MaxEnt setting of global SDM (G-SDM) for Carassius gibelio based on a composite overall score. The overall score was calculated as a weighted sum of normalized tuning metrics: Mean.AUC.DIFF (overfitting indicator), OR10 (10% omission rate), ORmin (minimum omission rate), and Mean.AUC (model discrimination). All metrics were normalized using min–max scaling; for Mean.AUC, the scale was inverted so that higher AUC values corresponded to a lower penalty. The best-performing MaxEnt model with optimal parameterization, including the regularization multiplier (RM) and the selected feature classes (FC: Linear, Quadratic, Product, Threshold, and Hinge), is located at the bottom of the Y-axis  

Input data:  
Elevation.tif, Parva_SDMglobClim.tif, Slope.tif, WaterLine.tif, WaterPoly.tif  from (Layers_of_ Armenia.Zip) and species occurrence record - Parva_Local_Full.csv  

Output data:  
MaxEnt_tuning_scored_results.csv  
MaxEnt_settings_score_plot.png  

Example: This example uses Armenia_TopWat_Layers.zip and P. parva data (Parva_Local_Full.csv) to demonstrate how the script works.  

________________________________________

## 7. Statistical comparison of SDM model performance. Script name: Compare_SDM_Metrics.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Statistical comparison of SDM model performance.  
Description:  This script performs: - Standardization of SDM evaluation metrics (AUC, TSS, CBI), Construction of combined performance indices - COMB_equal, COMB_expert, COMB_invvar; Statistics for combined indices, Model ranking, Boxplots and density plots, ANOVA, Kruskals-Wallis, pairwise Wilcoxon tests (BH correction) of evaluation metrics, Correlation analysis and heatmap visualization,  

Input: CSV file with columns:  model, fold, AUC, TSS, CBI - Parva_Compare_SDM_V2.csv, Caras_Compare_SDM_V2.csv  

Output:  Summary CSV tables, Metrics_Boxplot.png, Metrics_Density.png,  statistical_reports.txt, Figure_Composite indices across SDM frameworks, ANOVA_results.txt, Correlation_Heatmap.png, KW_results.txt, Wilcoxon_results.txt, Correlation_Table.csv, ANOVA_COMB_results.txt, KW_COMB_results.txt,  

Example: This example uses Parva_Compare_SDM_V2.csv and Caras_Compare_SDM_V2.csv (Compare_SDM.zip) to demonstrate how the script works.  

________________________________________

## 8. Global SDM modelling (G-SDM) with MaxEnt. Script name: GSDM_Global_Model_and_Projection.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Global SDM modelling (G-SDM) with MaxEnt and projection to South Caucasus / Armenia  
This script performs: 1. Global-scale SDM calibration using MaxEnt, 2. Spatial block cross-validation (BlockCV), 3. Background generation (three alternative strategies), 4. Calculation of AUC, TSS and Boyce index (global only), 5. Final model training on all global occurrences, 6. Projection to South Caucasus / Armenia, 7. Uncertainty estimation (SD across folds),  8. Variable importance and response curves  

Input:  Global environmental rasters (30 arc-sec), Regional environmental rasters (30 arc-sec), Global occurrence CSV with columns: Lon, Lat  

Output: SDM_Global_final.tif,  SDM_Regional_projection.tif,  CV metrics summary (AUC, TSS, Boyce),  Ensemble SD and mean rasters,  Variable_Importance.csv,  Response curves, Final Boyce curve  

Example: This example uses Global environmental rasters (from https://www.worldclim.org/data/worldclim21.html - Bio_03.tif, Bio_08.tif, Bio_12.tif, Bio_18.tif, Windmean.tif),  Regional Layers (Bio_03.tif, Bio_08.tif, Bio_12.tif, Bio_18.tif, Windmean.tif)  and Parva_Global_All.csv from the Armenia_Climate.zip  to demonstrate how the script works  

________________________________________

## 10. Regional SDM modelling (R-SDM) with MaxEnt. Script name: R-SDM_Creation_and_Projection.R

Authors: Varos Petrosyan, Fedor Osipov  
Purpose: Regional SDM modelling (R-SDM) with MaxEnt  
This script performs: 1. Regional-scale SDM calibration using MaxEnt,  Spatial block cross-validation (BlockCV),  Background generation (three alternative strategies),  Calculation of AUC, TSS and Boyce index,  Final model training on all regional occurrences,  Projection to South Caucasus,  Uncertainty estimation (SD across folds),  MESS raster estimation (across folds) , Variable importance and response curves  

IMPORTANT: All evaluation metrics are computed ONLY at the regional scale  

Input: Regional environmental rasters (30 arc-sec), Regional occurrence CSV with columns: Lon, Lat  

Output: results_summary.csv, SDM_Regional_final.tif, SDM_Regional_mean.tif, SDM_Regional_SD.tif, SDM_Regional_MESS.tif,  Variable_Contribution_PC.csv, Variable_Importance_PI.csv, Response_Curves.png  

Example: This example uses Regional Layers (Bio_03.tif, Bio_08.tif, Bio_12.tif, Bio_18.tif, Windmean.tif - Armenia_Climate_Layers.zip) and (Elevation.tif, Slope.tif, WaterLine.tif, WaterPoly.tif ), and Parva_local_Full.csv from Armenia_TopWat_Layers.zip  to demonstrate how the script works  

________________________________________

## List of files

1. Spatial Thinning of Species Occurrence Records (SORs)  
1.1. spThin_EcoInf.R  
1.2. SORs_Full_csv.zip  
1.3. SORs_Reduced_csv.zip  

2. Assessment of multicollinearity of predictor variables  
2.1. VIF_Spearman_EcoInf.R  
2.2. Predictors.zip (These layers are presented at a resolution of 30 arc seconds due to a file size limitation.)  

3. Estimates of climatic niche metrics for the species under the ECOUE concept  
3.1. ECOUE_species_EcoInf.R  
3.2. Species_NatInv-PsAb_points.zip  

4. Bootstrap assessment of the Continuous Boyce Index (CBI)  
4.1. CBI_EcoInf.R  
4.2. SDMs_SORs.zip (The SDM raster layers of suitable habitats for P. parva in the study region are presented at a resolution of 30 arc seconds due to a file size limitation.)  

5. Estimation of accessible area (M) radius for SDM, Script name: Calculate_M_Area_Radius.R  
5.1. Species occurrence records - Carasus_ Parva_global_Sors.zip(Carassius gibelio: Carasus_Global_All.csv, Pseudorasbora parva: Parva_Global_All.csv)  

6. MaxEnt Hyperparameter Optimization using biomod2, Script name: Maxent_parameters_using_Biomod2_tuning.R  
6.1. Layers_of_ Armenia.Zip including Elevation.tif, Parva_SDMglobClim.tif, Slope.tif, WaterLine.tif, WaterPoly.tif , and species occurrence record - Parva_Local_Full.csv  

7. Statistical comparison of SDM model performance, Script name: Compare_SDM_Metrics.R  
7.1. CSV files with columns:  model, fold, AUC, TSS, CBI included in Compare_SDM.zip  

8. Global SDM modelling (G-SDM) with MaxEnt, Script name: GSDM_Global_Model_and_Projection.R  
8.1. Armenia_Climate_Layers.zip  

9. Regional SDM modelling (R-SDM) with MaxEnt, Script name: R-SDM_Creation_and_Projection.R  
9.1. Armenia_Climate_Layers.zip and Armenia_TopWat_Layers.zip  

## License

This project is licensed under the MIT License — see the [LICENSE](./LICENSE) file for details.  
Copyright (c) 2025 Varos Petrosyan and Fedor Osipov  

For questions or collaboration inquiries, please contact:  
Varos Petrosyan  
Email: vgpetrosyan@gmail.com
