###############################################################
# Script name: GSDM_Global_Model_and_Projection.R
# Purpose: Global SDM modelling (G-SDM) with MaxEnt
# Author: Varos Petrosyan, Fedor Osipov
# Project: Multi-scale SDM analysis (South Caucasus)
# R version: 3.6.2+
#
# Description:
# This script performs:
# 1. Global-scale SDM calibration using MaxEnt
# 2. Spatial block cross-validation (BlockCV)
# 3. Background generation (three alternative strategies)
# 4. Calculation of AUC, TSS and Boyce index (global only)
# 5. Final model training on all global occurrences
# 6. Projection to South Caucasus / Armenia
# 7. Uncertainty estimation (SD across folds)
# 8. Variable importance and response curves
#
# IMPORTANT:
# - All evaluation metrics are computed ONLY at the global scale
# - Regional layers are used exclusively for projection
#
# INPUT:
# - Global environmental rasters (30 arc-sec)
# - Regional environmental rasters (30 arc-sec)
# - Global occurrence CSV with columns: Lon, Lat
#
# OUTPUT:
# - SDM_Global_final.tif
# - SDM_Regional_projection.tif
# - CV metrics summary (AUC, TSS, Boyce)
# - Ensemble SD and mean rasters
# - Variable_Importance.csv
# - Response curves
# - Final Boyce curve
####  License
# This code is released under the MIT License.
# You are free to use, modify, and redistribute with attribution.
###############################################################

###############################################################

# =============================================================
# 0. Java memory
# =============================================================

options(java.parameters = "-Xmx4g")

# =============================================================
# 1. Libraries
# =============================================================

library(rJava)
library(dismo)
library(raster)
library(blockCV)
library(ecospat)
library(ggplot2)
library(dplyr)
library(sf)
sf_use_s2(FALSE)

# =============================================================
# 2. Paths and folders
# =============================================================

maxent_path <- "maxent.jar_directory" ## ... E:/maxent/Maxent_old/maxent.jar

output_folder <- "Results_output_directory" ##
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

logfile <- file.path(output_folder, "processing_log.txt")

# =============================================================
# 3. Environmental predictors
# =============================================================

# ---- Global predictors (used for ALL evaluations)
### Need load Global environmental rasters -
### from https://www.worldclim.org/data/worldclim21.html 
env_glob <- stack(c(
  "Gl_data_Directory_name/Bio_03.tif",
  "Gl_data_Directory_name/Bio_08.tif",
  "Gl_data_Directory_name/Bio_12.tif",
  "Gl_data_Directory_name/Bio_18.tif",
  "Gl_data_Directory_name/Windmean.tif"
))

# ---- Regional predictors 
### Loading data from the archive -  Armenia_Climate.zip
env_arm <- stack(c(
  "Reg_data_Directory_name/Bio_03.tif",
  "Reg_data_Directory_name/Bio_08.tif",
  "Reg_data_Directory_name/Bio_12.tif",
  "Reg_data_Directory_name/Bio_18.tif",
  "Reg_data_Directory_name/Windmean.tif"
))

# =============================================================
# 4. Occurrence data (GLOBAL)
# =============================================================
### Loading Parva_Global_All.csv data from the archive Armenia_Climate.zip
occ_glob <- read.csv("../Parva_Global_All.csv")
if(!"presence" %in% names(occ_glob)) occ_glob$presence <- 1

coordinates(occ_glob) <- ~Lon+Lat
proj4string(occ_glob) <- CRS("+proj=longlat +datum=WGS84")
occ_sf <- st_as_sf(occ_glob, coords = c("Lon", "Lat"), crs = 4326)

# =============================================================
# 5. Spatial CV (BlockCV) – GLOBAL
# =============================================================

set.seed(123)

sb <- spatialBlock(
  speciesData = occ_sf,
  rasterLayer = env_glob,
  theRange    = 25000,
  k           = 10,
  selection   = "random",
  iteration   = 100,
  biomod2Format = FALSE
)

folds_vec <- sb$foldID

# =============================================================
# 6. Background generation (GLOBAL)
# =============================================================

bg_method <- 3  # 1 = full raster, 2 = random, 3 = buffer-based M
bg_size   <- 10000

if(bg_method == 1){
  all_coords <- coordinates(env_glob)
}

if(bg_method == 2){
  set.seed(123)
  bg_all <- sampleRandom(env_glob, size = bg_size, sp = TRUE, na.rm = TRUE)
  all_coords <- coordinates(bg_all)
}

if(bg_method == 3){
  occ_sf_m <- st_transform(occ_sf, 3857)
  buffer_radius_m <- 155000
  buf_m <- st_union(occ_sf_m) %>% st_buffer(dist = buffer_radius_m)
  buf_wgs84 <- st_transform(buf_m, 4326)
  buf_sp <- as(buf_wgs84, "Spatial")
  
  env_buf <- mask(env_glob, buf_sp)
  set.seed(123)
  bg_all <- sampleRandom(env_buf, size = bg_size, sp = TRUE, na.rm = TRUE)
  all_coords <- coordinates(bg_all)
}

# =============================================================
# 7. Cross-validation loop (GLOBAL METRICS ONLY)
# =============================================================

produce_rasters_in_cv <- FALSE

n_folds <- 10
auc_vec   <- rep(NA, n_folds)
tss_vec   <- rep(NA, n_folds)
boyce_vec <- rep(NA, n_folds)

preds_list  <- list()
models_list <- list()

cat("Run started: ", format(Sys.time()), "\n", file = logfile)

# Добавим лог-файл
logfile <- file.path(output_folder, "processing_log.txt")
cat("Run started: ", format(Sys.time()), "
", file = logfile, append = FALSE)

for(i in 1:n_folds){
  cat(sprintf("[%s] Fold %d - start
", format(Sys.time(), "%H:%M:%S"), i), file = logfile, append = TRUE)
  cat(sprintf("[%s] Fold %d - start
", format(Sys.time(), "%H:%M:%S"), i))
  
  train_idx <- which(folds_vec != i)
  test_idx  <- which(folds_vec == i)
  
  if(length(test_idx)==0){
    cat(sprintf("Fold %d skipped (no test)
", i), file = logfile, append = TRUE)
    next
  }
  
  occ_train <- occ_glob[train_idx,]
  occ_test  <- occ_glob[test_idx,]
  
  tr_coords <- coordinates(occ_train)
  te_coords <- coordinates(occ_test)
  
  # ---- MODEL TRAINING (global env_glob) ----
  t0_train <- Sys.time()
  me <- maxent(
    x = env_glob,
    p = tr_coords,
    path = output_folder,
    args = c("betamultiplier=1", "responsecurves=true", "jackknife=true")
  )
  # args = c("betamultiplier=1", "linear=true", "quadratic=true", "hinge=true",
  #         "product=false", "threshold=false")
  
  
  models_list[[i]] <- me
  dt_train <- as.numeric(difftime(Sys.time(), t0_train, units = "secs"))
  cat(sprintf("Model trained in %.1f s
", dt_train), file = logfile, append = TRUE)
  cat(sprintf("Model trained in %.1f s
", dt_train))
  
  # ---- optionally produce full raster for this fold (slow) ----
  if(produce_rasters_in_cv){
    t0_pred <- Sys.time()
    pred_rast <- predict(env_glob, me)
    preds_list[[i]] <- pred_rast
    dt_pred <- as.numeric(difftime(Sys.time(), t0_pred, units = "secs"))
    cat(sprintf("Full raster predicted in %.1f s
", dt_pred), file = logfile, append = TRUE)
    cat(sprintf("Full raster predicted in %.1f s
", dt_pred))
  } else {
    pred_rast <- NULL
    preds_list[[i]] <- NULL
  }
  
  # ---- BACKGROUND FOR METRICS (sample from global all_coords) ----
  # use moderate n_bg for CV speed; 
  n_bg_cv <- 2000
  bg_pool <- all_coords
  if(is.null(bg_pool) || nrow(bg_pool) == 0){ stop("bg pool is empty") }
  if(nrow(bg_pool) <= n_bg_cv){ bg_coords <- bg_pool } else { bg_coords <- bg_pool[sample(1:nrow(bg_pool), n_bg_cv), , drop=FALSE] }
  
  # ---- Extract env variables only at points (fast) ----
  env_te <- extract(env_glob, te_coords)
  env_bg <- extract(env_glob, bg_coords)
  
  # remove NA rows
  ok_te <- which(!apply(env_te, 1, function(x) any(is.na(x))))
  ok_bg <- which(!apply(env_bg, 1, function(x) any(is.na(x))))
  if(length(ok_te) == 0){ warning(sprintf("Fold %d: no valid env at test points - skipping" , i)); next }
  env_te_df <- as.data.frame(env_te[ok_te, , drop=FALSE])
  te_coords_ok <- te_coords[ok_te, , drop=FALSE]
  
  if(length(ok_bg) == 0){
    warning(sprintf("Fold %d: no valid env at sampled bg points - AUC/TSS not computed", i))
    env_bg_df <- NULL
  } else {
    env_bg_df <- as.data.frame(env_bg[ok_bg, , drop=FALSE])
    bg_coords_ok <- bg_coords[ok_bg, , drop=FALSE]
  }
  
  # ---- Predict at points using model (fast) ----
  pred_te_vals <- tryCatch({ predict(me, env_te_df) }, error=function(e){ message("Predict test error: ", e$message); rep(NA, nrow(env_te_df)) })
  pred_bg_vals <- NULL
  if(!is.null(env_bg_df) && nrow(env_bg_df)>0){ pred_bg_vals <- tryCatch({ predict(me, env_bg_df) }, error=function(e){ message("Predict bg error: ", e$message); rep(NA, nrow(env_bg_df)) }) }
  
  # ---- AUC and TSS using numeric preds ----
  ev <- NULL
  if(!is.null(pred_bg_vals) && length(pred_te_vals)>0){
    ev <- tryCatch({ evaluate(p = pred_te_vals, a = pred_bg_vals) }, error=function(e){ message("evaluate error: ", e$message); NULL })
    eval_list[[i]] <- ev
    auc_vec[i] <- if(!is.null(ev)) ev@auc else NA
    
    # TSS via spec_sens threshold
    tss_vec[i] <- NA
    if(!is.null(ev)){
      th <- tryCatch(threshold(ev, "spec_sens"), error=function(e) NA)
      if(!is.na(th)){
        test_pres <- ifelse(pred_te_vals >= th, 1, 0)
        bg_abs    <- ifelse(pred_bg_vals < th, 1, 0)
        TPR <- sum(test_pres==1, na.rm=TRUE) / sum(!is.na(test_pres))
        TNR <- sum(bg_abs==1, na.rm=TRUE) / sum(!is.na(bg_abs))
        tss_vec[i] <- TPR + TNR - 1
      }
    }
  } else {
    eval_list[[i]] <- NULL
    auc_vec[i] <- NA
    tss_vec[i] <- NA
  }
  
  # ---- Boyce for fold: use larger sampled fit background (fast) ----
  boyce_vec[i] <- NA
  n_fit_bg <- 10000
  if(nrow(all_coords) <= n_fit_bg){ fit_bg_coords <- all_coords } else { fit_bg_coords <- all_coords[sample(1:nrow(all_coords), n_fit_bg), , drop=FALSE] }
  env_fit_bg <- extract(env_glob, fit_bg_coords)
  ok_fit <- which(!apply(env_fit_bg,1,function(x) any(is.na(x))))
  if(length(ok_fit) > 100){
    env_fit_bg_df <- as.data.frame(env_fit_bg[ok_fit, , drop=FALSE])
    pred_fit_sample <- tryCatch({ predict(me, env_fit_bg_df) }, error=function(e) NULL)
    if(!is.null(pred_fit_sample) && length(pred_te_vals) > 5){
      boyce_obj <- tryCatch(ecospat.boyce(fit = pred_fit_sample, obs = pred_te_vals, nclass = 0), error=function(e) NULL)
      boyce_vec[i] <- if(!is.null(boyce_obj) && !is.null(boyce_obj$cor)) boyce_obj$cor else NA
    }
  }
  
  dt_fold <- as.numeric(difftime(Sys.time(), t0_train, units = "secs"))
  cat(sprintf("[%s] Fold %d done — elapsed=%.1f s  AUC=%.3f  TSS=%s  Boyce=%s
", format(Sys.time(), "%H:%M:%S"),
              i, dt_fold,
              ifelse(is.na(auc_vec[i]), NA, auc_vec[i]),
              ifelse(is.na(tss_vec[i]), "NA", format(round(tss_vec[i],3), nsmall=3)),
              ifelse(is.na(boyce_vec[i]), "NA", format(round(boyce_vec[i],3), nsmall=3)) ), file = logfile, append = TRUE)
  cat(sprintf("[%s] Fold %d done — elapsed=%.1f s  AUC=%.3f  TSS=%s  Boyce=%s
", format(Sys.time(), "%H:%M:%S"),
              i, dt_fold,
              ifelse(is.na(auc_vec[i]), NA, auc_vec[i]),
              ifelse(is.na(tss_vec[i]), "NA", format(round(tss_vec[i],3), nsmall=3)),
              ifelse(is.na(boyce_vec[i]), "NA", format(round(boyce_vec[i],3), nsmall=3)) ))
}

# =============================================================
# SUMMARY
# =============================================================
# =============================================================

results_df <- data.frame(fold=1:n_folds, AUC=auc_vec, TSS=tss_vec, Boyce=boyce_vec)
write.csv(results_df, file.path(output_folder,"cv_results_summary.csv"), row.names=FALSE)

# =============================================================
# FINAL MODEL -  Global prediction
# =============================================================
final_mod <- maxent(x = env_glob, p = coordinates(occ_glob),
                    path = output_folder,
                    args = c("betamultiplier=1", "responsecurves=true", "jackknife=true"))

final_pred <- predict(env_glob, final_mod)
#writeRaster(final_pred, file.path(output_folder,"SDM_Global_final.tif"), overwrite=TRUE)

writeRaster(
  final_pred,
  file.path(output_folder,"SDM_Global_final.tif"),
  overwrite=TRUE,
  options=c("COMPRESS=LZW", "BIGTIFF=YES")
)



# ---- Regional projection ----
final_pred_regional <- predict(env_arm, final_mod)
writeRaster(final_pred_regional, file.path(output_folder,"SDM_Regional_projection.tif"), overwrite=TRUE)
plot(final_pred_regional )

##### Modification: save the final SD of the fold forecast
sd_pred   <- calc(stack(preds_list), fun = sd)
writeRaster(sd_pred,    filename = file.path(output_folder, "SDM_Regional_SD.tif"),    format = "GTiff", overwrite = TRUE)

# ---- Ensemble mean/sd ----
#replace the direct stack(preds_list) with a safe version that takes into account NULL/non-uniform rasters
is_raster_obj <- function(x){ inherits(x, c("RasterLayer","RasterStack","RasterBrick")) }
valid_preds <- preds_list[ unlist(lapply(preds_list, is_raster_obj)) ]


if(length(valid_preds)==0){
  message("No per-fold rasters produced — skipping mean/sd.")
  mean_pred <- NULL
  sd_pred <- NULL
} else {
  if(!exists("final_pred") || is.null(final_pred)) stop("final_pred missing!")
  
  message("Aligning per-fold rasters to final_pred grid...")
  aligned_preds <- vector("list", length(valid_preds))
  
  for(j in seq_along(valid_preds)){
    rj <- valid_preds[[j]]
    
    if(!compareRaster(rj, final_pred, stopiffalse=FALSE)){
      rj <- try(resample(rj, final_pred, method="bilinear"), silent=TRUE)
      if(inherits(rj, "try-error")){
        message(sprintf("Resample failed for raster %d — skipped", j))
        next
      }
    }
    aligned_preds[[j]] <- rj
  }
  
  aligned_preds <- aligned_preds[ sapply(aligned_preds, is_raster_obj) ]
  
  if(length(aligned_preds)==0){
    message("All rasters failed alignment — skipping mean/sd.")
    mean_pred <- NULL
    sd_pred <- NULL
  } else {
    pred_stack <- try(stack(aligned_preds), silent=TRUE)
    if(inherits(pred_stack, "try-error")){
      message("Stack creation error — skipping mean/sd.")
      mean_pred <- NULL
      sd_pred <- NULL
    } else {
      mean_pred <- calc(pred_stack, fun=mean)
      sd_pred   <- calc(pred_stack, fun=sd)
      writeRaster(mean_pred, file.path(output_folder,"SDM_Global_mean.tif"), overwrite=TRUE)
      writeRaster(sd_pred, file.path(output_folder,"SDM_Global_SD.tif"), overwrite=TRUE)
    }
  }
}

# =============================================================
# VARIABLE IMPORTANCE
# =============================================================
pc_rows <- grep("contribution", rownames(final_mod@results), ignore.case=TRUE)
pi_rows <- grep("permutation.importance", rownames(final_mod@results), ignore.case=TRUE)
pc <- if(length(pc_rows)) final_mod@results[pc_rows,1] else rep(NA,nlayers(env_glob))
pi <- if(length(pi_rows)) final_mod@results[pi_rows,1] else rep(NA,nlayers(env_glob))
vars_df <- data.frame(Variable=names(env_glob), PercentContribution=pc, PermutationImportance=pi)
write.csv(vars_df, file.path(output_folder, "Variable_Importance.csv"), row.names=FALSE)

# =============================================================
# RESPONSE CURVES
# =============================================================
png(file.path(output_folder, "Response_Curves.png"), width=1200, height=900)
try(response(final_mod), silent=TRUE)
dev.off()

# =============================================================
# FINAL BOYCE INDEX (on global final predictior!)
# =============================================================
all_vals <- values(final_pred); all_vals <- all_vals[!is.na(all_vals)]
occ_vals <- extract(final_pred, coordinates(occ_glob)); occ_vals <- occ_vals[!is.na(occ_vals)]

final_boyce <- NA
boyce_obj <- NULL
if(length(occ_vals)>5 && length(all_vals)>100){
  boyce_obj <- tryCatch(ecospat.boyce(fit=all_vals, obs=occ_vals, nclass=0), error=function(e) NULL)
  if(!is.null(boyce_obj)) final_boyce <- boyce_obj$cor
}
###


png(file.path(output_folder,"Final_Boyce.png"), width=700, height=600)
if(!is.null(boyce_obj) &&
   !is.null(boyce_obj$classes) &&
   !is.null(boyce_obj$F.ratio) &&
   length(boyce_obj$classes) > 1 &&
   length(boyce_obj$F.ratio) > 1){
  
  plot(boyce_obj$F.ratio ~ boyce_obj$classes, type="b",
       xlab="Suitability class", ylab="P/E (Boyce)",
       main=paste0("Final Boyce = ", round(final_boyce,3)))
  abline(h=1, lty=2)
  
} else {
  plot.new()
  text(0.5,0.5,"Final Boyce unavailable (insufficient data)")
}

dev.off()


cat("DONE. All metrics are correctly calculated at the GLOBAL level.")

######################## End ###########################################################

