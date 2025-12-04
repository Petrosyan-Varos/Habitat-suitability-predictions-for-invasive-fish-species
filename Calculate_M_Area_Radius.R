###############################################################
# Script name: Calculate_M_Area_Radius.R
# Purpose: Estimation of accessible area (M) radius for SDM
# Author: Varos Petrosyan, Fedor Osipov
# Project: Multi-scale SDM in Armenia
# Species: Pseudorasbora parva / Carassius gibelio
# Date: 2025
#
# Description:
# This script estimates the radius of the accessible area (M)
# used for background point generation in SDM workflows.
# The radius is derived from nearest-neighbor distances
# between species occurrence points.
#
# Method:
# 1. Convert geographic coordinates to kilometers
# 2. Compute all pairwise distances
# 3. Estimate nearest-neighbor distance per point
# 4. Define M radius as: mean(NN) + SD(NN)
####  License
# This code is released under the MIT License.
# You are free to use, modify, and redistribute with attribution.
###############################################################

# ----------------------------
# 1. Input data
# ----------------------------

PathToSources <- "Working_directory" #  ## 

## Loading data from the archive Carasus_Parva_global_Sors.zip
# Choose one species file
fileSpecName <- paste0(PathToSources, "/Parva_Global_All.csv")
fileSpecName <- paste0(PathToSources, "/Carasus_Global_All.csv")

# CSV must contain columns: Lon, Lat
coords <- read.csv(fileSpecName, header = TRUE)

# Inspect structure
head(coords)

# Safety check
stopifnot(all(c("Lon", "Lat") %in% names(coords)))

# ----------------------------
# 2. Coordinate conversion (degrees → kilometers)
# ----------------------------
# Local planar approximation for Armenia

mean_lat <- mean(coords$Lat, na.rm = TRUE)

coords$x_km <- coords$Lon * 111 * cos(mean_lat * pi / 180)  # longitude → km
coords$y_km <- coords$Lat * 111                              # latitude  → km

# ----------------------------
# 3. Pairwise distance matrix
# ----------------------------

dist_matrix <- as.matrix(
  dist(coords[, c("x_km", "y_km")], method = "euclidean")
)

# ----------------------------
# 4. Nearest-neighbor distance per point
# ----------------------------

nearest_dist <- apply(dist_matrix, 1, function(d){
  min(d[d > 0], na.rm = TRUE)
})

# ----------------------------
# 5. Summary statistics
# ----------------------------

mean_nn <- mean(nearest_dist, na.rm = TRUE)
sd_nn   <- sd(nearest_dist, na.rm = TRUE)

# ----------------------------
# 6. Accessible area radius (M)
# ----------------------------
# Primary definition:
radius_M <- mean_nn + sd_nn

# Alternative option (commented):
# radius_M <- 1.5 * mean_nn

# ----------------------------
# 7. Output
# ----------------------------

cat("========================================\n")
cat("Accessible area (M) radius estimation\n")
cat("========================================\n")
cat("Mean nearest-neighbor distance (km):", round(mean_nn, 2), "\n")
cat("Standard deviation (km):", round(sd_nn, 2), "\n")
cat("Recommended M radius (km):", round(radius_M, 2), "\n")
cat("========================================\n")
