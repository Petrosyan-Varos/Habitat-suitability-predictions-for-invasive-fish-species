############################################################
#  MaxEnt Hyperparameter Optimization using biomod2
#  Multi-metric weighted tuning approach
#  Author: Varos Petrosyan, Fedor Osipov
#  Project: Armenia 2025 SDM
#  Script: Maxent_parameters_using_Biomod2_tuning.R
####  License
# This code is released under the MIT License.
# You are free to use, modify, and redistribute with attribution.
###############################################################

### =========================================================
### 1. Load required libraries
### =========================================================

library(biomod2)
library(raster)
library(ENMeval)
library(dismo)
library(rJava)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(doParallel)

### =========================================================
### 2. Define working directories and input data
### =========================================================

driveName <- "working_drive_name"
# Loading data from the archive
### Select environmental dataset - Elevation.tif, Parva_SDMglobClim.tif, Slope.tif, WaterLine.tif, WaterPoly.tif
## from Armenia_TopWat_Layers.Zip

MainDirectory <- paste0(driveName, "Armenia_2025")
## Parva_Local_Full.csv
# Select environmental dataset from 
CurrentClimDir <- "Layers_directory"
#                
# Input species file (select one)
PathToSources <- "Sors_directory"   
## Select Sors  dataset 
fileSpecName  <- paste0(PathToSources, "/Parva_Local_Full.csv")

# Load species presence data
DataSpecies <- read.csv(fileSpecName, fill = TRUE, header = TRUE)
NN <- nrow(DataSpecies)

myResp    <- rep(1, NN)
myRespXY  <- DataSpecies
myRespName <- "Parva"   

### =========================================================
### 3. Load predictor variables (select ONE SDM type)
### =========================================================

# Example: ON-SDM / MN-SDM (regional multi-source predictors)
myExpl <- raster::stack(
  raster(paste0(CurrentClimDir, "Elevation.tif")),
  raster(paste0(CurrentClimDir, "Parva_SDMglobClim.tif")),
  raster(paste0(CurrentClimDir, "Slope.tif")),
  raster(paste0(CurrentClimDir, "WaterLine.tif")),
  raster(paste0(CurrentClimDir, "WaterPoly.tif"))
)

### =========================================================
### 4. Format data for biomod2
### =========================================================

myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,
  expl.var = myExpl,
  resp.xy  = myRespXY,
  resp.name = myRespName,
  PA.nb.rep = 1,
  PA.nb.absences = 5000,
  PA.strategy = "random"
)

### =========================================================
### 5. Parallel computation setup
### =========================================================

cl <- makeCluster(6)
registerDoParallel(cl)

time.seq <- system.time({
  Biomod.tuning <- BIOMOD_tuning(
    myBiomodData,
    models   = c("MAXENT.Phillips"),
    Yweights = NULL,
    metric.ME = "TSS"
  )
})

stopCluster(cl)

### =========================================================
### 6. Visual inspection of tuning results
### =========================================================

par(mfrow = c(1, 2))
eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results,
          "Mean.AUC",
          legend.position = "bottomright")

eval.plot(Biomod.tuning$tune.MAXENT.Phillips@results,
          "Mean.AUC.DIFF",
          variance = "Var.AUC.DIFF")

### =========================================================
### 7. Extract and prepare tuning metrics
### =========================================================

res <- Biomod.tuning$tune.MAXENT.Phillips@results

df <- res %>%
  dplyr::select(settings, features, rm,
                Mean.AUC.DIFF, Mean.OR10, Mean.ORmin, Mean.AUC) %>%
  filter(
    !is.na(Mean.AUC.DIFF),
    !is.na(Mean.OR10),
    !is.na(Mean.ORmin),
    !is.na(Mean.AUC)
  )

### =========================================================
### 8. Min–max normalization of metrics
### =========================================================

minmax <- function(x){
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

df_norm <- df %>%
  mutate(
    s_AUC.DIFF = minmax(Mean.AUC.DIFF),   # lower is better
    s_OR10     = minmax(Mean.OR10),
    s_ORmin    = minmax(Mean.ORmin),
    s_AUC_raw  = minmax(Mean.AUC),
    s_AUC      = 1 - s_AUC_raw            # invert AUC
  )

### =========================================================
### 9. Weighted Overall Score
### =========================================================

# Final weights used in the study
w <- c(AUC.DIFF = 0.60, OR10 = 0.15, ORmin = 0.10, AUC = 0.15)

df_scored <- df_norm %>%
  mutate(
    Score = w["AUC.DIFF"] * s_AUC.DIFF +
      w["OR10"]     * s_OR10 +
      w["ORmin"]    * s_ORmin +
      w["AUC"]      * s_AUC
  ) %>%
  arrange(Score) %>%
  mutate(rank = row_number())

### =========================================================
### 10. Top-10 optimal MaxEnt configurations
### =========================================================

top10 <- df_scored %>% slice(1:10)
print(top10 %>%
        dplyr::select(rank, settings, features, rm,
                      Mean.AUC, Mean.AUC.DIFF,
                      Mean.OR10, Mean.ORmin, Score))

### =========================================================
### 11. Visualization of tuning metrics
### =========================================================

df_long <- df_scored %>%
  select(settings, features, rm,
         Mean.AUC.DIFF, Mean.OR10, Mean.ORmin, Mean.AUC,
         Score, rank) %>%
  pivot_longer(
    cols = c(Mean.AUC.DIFF, Mean.OR10, Mean.ORmin, Mean.AUC),
    names_to = "metric", values_to = "value"
  )

df_long$metric <- factor(
  df_long$metric,
  levels = c("Mean.AUC.DIFF", "Mean.OR10", "Mean.ORmin", "Mean.AUC")
)

# Metric distributions
ggplot(df_long, aes(x = forcats::fct_inorder(settings), y = value)) +
  geom_point(aes(color = metric), alpha = 0.6) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(
    x = "Model configuration",
    y = "Metric value",
    title = "Tuning metrics for MaxEnt configurations"
  )

# Overall Score barplot
ggplot(df_scored, aes(x = reorder(settings, Score), y = Score)) +
  geom_col(fill = "#1f77b4", width = 0.7) +
  coord_flip() +
  labs(
    x = "Settings (ordered by Score)",
    y = "Overall Score (lower is better)",
    title = "Ranked MaxEnt configurations by combined score"
  ) +
  theme_minimal(base_size = 14)

# Boxplots by RM and feature class
ggplot(df_scored, aes(x = as.factor(rm), y = Mean.AUC.DIFF)) +
  geom_boxplot() +
  facet_wrap(~features, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Regularization multiplier (RM)",
    y = "Mean.AUC.DIFF",
    title = "Overfitting patterns by feature classes and RM"
  )

### =========================================================
### 12. Export final results
### =========================================================

write.csv(df_scored,
          "MaxEnt_tuning_scored_results.csv",
          row.names = FALSE)

ggsave("MaxEnt_settings_score_plot.png",
       width = 8, height = 6, dpi = 300)

############################################################
# End of script
############################################################
