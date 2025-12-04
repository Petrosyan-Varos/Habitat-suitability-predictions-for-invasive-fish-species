###############################################################
# Script name: R-SDM_Creation_and_Projection.R
# Purpose: Regional SDM modelling (R-SDM) with MaxEnt
# Author: Varos Petrosyan, Fedor Osipov
# Project: Multi-scale SDM analysis (South Caucasus)
# R version: 3.6.2+
#
# Description:
# This script performs:
# 1. Reginal -scale SDM calibration using MaxEnt
# 2. Spatial block cross-validation (BlockCV)
# 3. Background generation (three alternative strategies)
# 4. Calculation of AUC, TSS and Boyce index 
# 5. Final model training on all reginal occurrences
# 6. Projection to South Caucasus 
# 7. Uncertainty estimation (SD across folds)
# 8. MESS raster estimation (across folds) 
# 9. Variable importance and response curves
#
# IMPORTANT:
# - All evaluation metrics are computed ONLY at the reginal scale
#
# INPUT:
# - Regional environmental rasters (30 arc-sec)
# - Global occurrence CSV with columns: Lon, Lat
#
# OUTPUT:
# - CV metrics summary (AUC, TSS and Boyce) 
# - SDM_Regional_final_30s.tif
# - SDM_Regional_mean_30s.tif
# - SDM_Regional_SD_30s.tif
# - SDM_Regional_MESS_30s.tif
# - Ensemble SD and mean rasters
# - Variable_Contribution_PC_30s.csv
# - Variable_Importance_PI_30s.csv
# - Response curves
####  License
# This code is released under the MIT License.
# You are free to use, modify, and redistribute with attribution.
###############################################################
###############################################################

options(java.parameters = "-Xmx8g")
library(rJava)
library(dismo)
library(raster)
library(terra)
library(blockCV)
library(ecospat)
library(ggplot2)
library(dplyr)
library(sf)
sf_use_s2(FALSE)

# ---------------- User parameters --------------------------------------------
maxent_path <- "E:/maxent/Maxent_old/maxent.jar"
maxent_path <- "maxent.jar_directory" ## ... maxent.jar
output_folder <- "Results_output_directory" ##
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# ---- Regional predictors 
### Loading data from the archives -  Armenia_Climate.zip +Armenia_TopWat_Layers.zip
env_files_30s <- c(
  "Reg_data_Directory_name/Bio_03.tif",
  "Reg_data_Directory_name/Bio_08.tif",
  "Reg_data_Directory_name/Bio_12.tif",
  "Reg_data_Directory_name/Bio_18.tif",
  "Reg_data_Directory_name/Windmean.tif",
  "Reg_data_Directory_nam/Elevation.tif",
  "Reg_data_Directory_nam/Slope.tif",
  "Reg_data_Directory_nam/WaterLine.tif",
  "Reg_data_Directory_nam/WaterPoly.tif"
)


### Loading Parva_Global_All.csv data from the archive Armenia_Climate.zip
 occ_file <- "../Parva_Global_All.csvs/Parva_Local_Full.csv"

# Background options
bg_method <- 3  # 2 = sampleRandom, 3 = buffer
bg_pool_sample_size <- 1000
n_bg_eval <- 1000
buffer_radius_m <- 21000 # 

# blockCV param (theRange chosen by you)
theRange_val <- 5000   # as requested
k_test_min <- 2
k_test_max <- 10
min_test_per_fold <- 10  # minimal acceptable number of presences in each test-fold

# If you want to force a k (skip auto-search), set force_k to a number, else NA
force_k <- NA  # e.g. force_k <- 5
k_folds_default <- 5  # fallback

set.seed(123)
use_crop <- FALSE

# MaxEnt args
 maxent_args_cv <- c("betamultiplier=1", "responsecurves=false", "jackknife=false")
 maxent_args_final <- c("betamultiplier=1", "responsecurves=true", "jackknife=true")

 ## Parva -V2  
 maxent_args_cv = c("betamultiplier=1.5", "linear=true", "quadratic=true", "hinge=False",
                    "product=false", "threshold=false", "responsecurves=true", "jackknife=true")
 maxent_args_final<-maxent_args_cv
 

# processing - Log file
logfile <- file.path(output_folder, "processing_log.txt")
cat("Run started: ", format(Sys.time()), "\n", file=logfile, append=FALSE)

# ---------------- Load data --------------------------------------------
cat("Reading rasters via terra...\n"); cat("Reading rasters via terra...\n", file=logfile, append=TRUE)
env_terra <- rast(env_files_30s)
names(env_terra) <- make.names(names(env_terra), allow_ = TRUE)

if(use_crop){
  occ_tmp <- read.csv(occ_file)
  bb <- ext(min(occ_tmp$Lon), max(occ_tmp$Lon), min(occ_tmp$Lat), max(occ_tmp$Lat))
  bb <- extend(bb, 0.1)
  env_terra <- crop(env_terra, bb)
}

# Occurrences ===
occ_local <- read.csv(occ_file, stringsAsFactors = FALSE)
if(!all(c("Lon","Lat") %in% names(occ_local))) stop("occ_local must have columns Lon and Lat")
if(!"presence" %in% names(occ_local)) occ_local$presence <- 1

occ_sf <- st_as_sf(occ_local, coords = c("Lon","Lat"), crs = 4326, remove = FALSE)
occ_sp <- as(occ_sf, "Spatial")

# For dismo/blockCV we use raster::stack built from terra object
env_raster <- tryCatch({
  stack(env_terra)        # convert SpatRaster -> RasterStack
}, error = function(e){
  stop("Cannot convert env_terra to raster::stack: ", e$message)
})
# Ensure identical names/order
names(env_raster) <- names(env_terra)

cat("CRS checks:\n", file=logfile, append=TRUE)
# cat(" - env_terra CRS: ", crs(env_terra), "\n", file=logfile, append=TRUE)
# cat(" - env_raster CRS: ", crs(env_raster), "\n", file=logfile, append=TRUE)
# cat(" - occ_sf CRS: ", st_crs(occ_sf)$wkt[1], "\n", file=logfile, append=TRUE)
cat("Number of presence records:", nrow(occ_local), "\n"); 
cat("Number of presence records:", nrow(occ_local), "\n", file=logfile, append=TRUE)

# ---------------- Auto-select k (2..10) using spatialBlock --------------------
evaluate_k <- function(k_value, occ_sf, env_raster, theRange, min_test){
  set.seed(123)
  sb_tmp <- spatialBlock(
    speciesData = occ_sf,
    rasterLayer = env_raster,
    theRange = theRange,
    k = k_value,
    selection = "random",
    iteration = 100,
    biomod2Format = FALSE
  )
  freq <- table(sb_tmp$foldID)
  ok <- all(as.integer(freq) >= min_test)
  list(k = k_value, freq = freq, ok = ok, sb = sb_tmp)
}

chosen_k <- NA
if(!is.na(force_k)){
  chosen_k <- force_k
  cat("force_k provided -> using k =", chosen_k, "\n"); cat("force_k provided -> using k =", chosen_k, "\n", file=logfile, append=TRUE)
} else {
  cat("Searching best k in range", k_test_min, ":", k_test_max, " (min test per fold =", min_test_per_fold, ")\n")
  cat("Searching best k in range", k_test_min, ":", k_test_max, " (min test per fold =", min_test_per_fold, ")\n", file=logfile, append=TRUE)
  res_list <- list()
  for(kv in seq(k_test_min, k_test_max)){
    r <- evaluate_k(kv, occ_sf, env_raster, theRange_val, min_test_per_fold)
    res_list[[as.character(kv)]] <- r
    cat(sprintf("k=%d : min per fold = %d, ok=%s\n", kv, min(as.integer(r$freq)), ifelse(r$ok,"YES","NO")), file = logfile, append = TRUE)
  }
  ok_k <- as.integer(names(Filter(function(x) x$ok, res_list)))
  if(length(ok_k) == 0){
    cat("No k in tested range satisfied min_test_per_fold. Choosing conservative k.\n", file=logfile, append=TRUE)
    N <- nrow(occ_local)
    k_cand <- seq(k_test_min, k_test_max)
    test_per_fold <- floor(N / k_cand)
    k_ok2 <- k_cand[test_per_fold >= floor(min_test_per_fold/2)]
    chosen_k <- if(length(k_ok2)>0) min(k_ok2) else k_folds_default
  } else {
    chosen_k <- min(ok_k)  # choose smallest accepted k for more stable folds
  }
  cat("Selected k =", chosen_k, "\n", file=logfile, append=TRUE)
}

# Final spatialBlock with chosen_k
set.seed(123)
sb <- spatialBlock(
  speciesData = occ_sf,
  rasterLayer = env_raster,
  theRange = theRange_val,
  k = chosen_k,
  selection = "random",
  iteration = 100,
  biomod2Format = FALSE
)


set.seed(123)
sb <- spatialBlock(
  speciesData = occ_sf,
  rasterLayer = env_raster,
  theRange = theRange_val,
  k = 5, # 5
  selection = "random",
  iteration = 100,
  biomod2Format = FALSE
)

folds_vec <- sb$foldID
cat("Final fold distribution:\n"); print(table(folds_vec)); cat("Final fold distribution:\n", file=logfile, append=TRUE); capture.output(print(table(folds_vec)), file=logfile, append=TRUE)

# ---------------- Background pool ------------------------------------------------
bg_method == 2

if(bg_method == 2){
  set.seed(123)
  samp <- terra::spatSample(env_terra, size = bg_pool_sample_size, method = "random", xy = TRUE, as.points = TRUE)
  all_coords <- coordinates(as(samp, "Spatial"))
}
bg_method <- 3
buffer_radius_m <- 21000 # 21000- Parva
if(bg_method == 3){
  occ_sf_m <- st_transform(occ_sf, 3857)
  buf_m <- st_union(occ_sf_m) %>% st_buffer(dist = buffer_radius_m)
  buf_wgs84 <- st_transform(buf_m, 4326)
  buf_sp <- as(buf_wgs84, "Spatial")
  env_buf <- mask(env_raster, buf_sp)
  if(!all(names(env_buf) %in% names(env_raster))){
    env_buf <- env_buf[[names(env_raster)]]
  }
  set.seed(123)
  samp <- sampleRandom(env_buf, size = bg_pool_sample_size, sp = TRUE, na.rm = TRUE)
  all_coords <- coordinates(samp)
}

# ---------------- Prepare CV storage ----------------------------------------
n_folds <- length(unique(folds_vec))
models_list <- vector("list", n_folds)
preds_list <- vector("list", n_folds)   # store per-fold prediction SpatRasters
eval_list <- vector("list", n_folds)
auc_vec <- rep(NA_real_, n_folds)
tss_vec <- rep(NA_real_, n_folds)
boyce_vec <- rep(NA_real_, n_folds)

# ---------------- CV loop ----------------------------------------------------
for(i in seq_len(n_folds)){
  fold_start <- Sys.time()
  cat(sprintf("\n=== Fold %d started at %s ===\n", i, format(fold_start)))
  cat(sprintf("\n=== Fold %d started at %s ===\n", i, format(fold_start)), file=logfile, append=TRUE)
  
  train_idx <- which(folds_vec != i)
  test_idx  <- which(folds_vec == i)
  if(length(test_idx) == 0){
    message("Fold ", i, " skipped (no test points)"); next
  }
  
  occ_train <- occ_sp[train_idx, ]
  occ_test  <- occ_sp[test_idx, ]
  tr_coords <- coordinates(occ_train)
  te_coords <- coordinates(occ_test)
  
  # Fit MaxEnt on training coords
  me <- tryCatch({
    maxent(x = env_raster, p = tr_coords,
           path = file.path(output_folder, paste0("fold_", i)),
           args = maxent_args_cv)
  }, error = function(e){
    message("MaxEnt error fold ", i, ": ", e$message)
    cat("MaxEnt error fold ", i, ": ", e$message, "\n", file=logfile, append=TRUE)
    return(NULL)
  })
  if(is.null(me)) next
  models_list[[i]] <- me
  
  # ---- Safe full-area prediction for this fold: raster::predict -> convert to SpatRaster ----
  pred_rast_r <- tryCatch({
    predict(env_raster, me)   # raster::predict (returns RasterLayer/Stack)
  }, error = function(e){
    message("predict(env_raster, me) error fold ", i, ": ", e$message)
    cat("predict(env_raster, me) error fold ", i, ": ", e$message, "\n", file=logfile, append=TRUE)
    return(NULL)
  })
  if(!is.null(pred_rast_r)){
    pred_rast_spat <- tryCatch({
      rast(pred_rast_r)  # convert to SpatRaster
    }, error = function(e){
      message("Conversion to SpatRaster failed fold ", i, ": ", e$message)
      cat("Conversion to SpatRaster failed fold ", i, ": ", e$message, "\n", file=logfile, append=TRUE)
      return(NULL)
    })
    preds_list[[i]] <- pred_rast_spat
  } else {
    preds_list[[i]] <- NULL
  }
  
  # Prepare background for evaluate: remove pres coords
  key_all <- paste(round(all_coords[,1],6), round(all_coords[,2],6), sep="_")
  pres_xy <- coordinates(occ_sp)
  key_pres <- paste(round(pres_xy[,1],6), round(pres_xy[,2],6), sep="_")
  keep_mask <- !(key_all %in% key_pres)
  bg_pool <- all_coords[keep_mask, , drop = FALSE]
  if(nrow(bg_pool) == 0){ message("Fold ", i, ": no bg after filtering -> skip"); next }
  if(nrow(bg_pool) <= n_bg_eval) bg_coords <- bg_pool else bg_coords <- bg_pool[sample(1:nrow(bg_pool), n_bg_eval), , drop = FALSE]
  
  # Extract predictor values (raster::extract) and filter NA
  vals_test <- raster::extract(env_raster, te_coords)
  if(is.null(vals_test) || nrow(vals_test)==0) { message("No test values extracted -> skip"); next }
  vals_test <- vals_test[complete.cases(vals_test), , drop = FALSE]
  vals_bg <- raster::extract(env_raster, bg_coords)
  vals_bg <- vals_bg[complete.cases(vals_bg), , drop = FALSE]
  if(nrow(vals_test) == 0 || nrow(vals_bg) == 0){
    message("Fold ", i, ": insufficient non-NA test/bg values -> skip"); next
  }
  
  # Ensure column names match training raster names (if dims match)
  if(ncol(vals_test) == nlayers(env_raster)) colnames(vals_test) <- names(env_raster)
  if(ncol(vals_bg) == nlayers(env_raster)) colnames(vals_bg) <- names(env_raster)
  
  # Predict on tabular data via maxent object (safe)
  pred_test_vals <- tryCatch({ predict(me, vals_test) }, error = function(e){
    message("Predict-on-test error fold ", i, ": ", e$message); cat("Predict-on-test error fold ", i, e$message, "\n", file=logfile, append=TRUE); return(NULL)
  })
  pred_bg_vals <- tryCatch({ predict(me, vals_bg) }, error = function(e){
    message("Predict-on-bg error fold ", i, ": ", e$message); cat("Predict-on-bg error fold ", i, e$message, "\n", file=logfile, append=TRUE); return(NULL)
  })
  if(is.null(pred_test_vals) || is.null(pred_bg_vals)){ next }
  
  # Evaluate
  ev <- tryCatch({ evaluate(p = pred_test_vals, a = pred_bg_vals) }, error = function(e){
    message("Evaluate error fold ", i, ": ", e$message); cat("Evaluate error fold ", i, ": ", e$message, "\n", file=logfile, append=TRUE); return(NULL)
  })
  eval_list[[i]] <- ev
  auc_vec[i] <- if(!is.null(ev)) ev@auc else NA
  
  # TSS
  tss_vec[i] <- NA
  if(!is.null(ev)){
    th <- tryCatch({ threshold(ev, 'spec_sens') }, error = function(e) NA)
    if(!is.na(th)){
      test_pres <- ifelse(pred_test_vals >= th, 1, 0)
      bg_abs <- ifelse(pred_bg_vals < th, 1, 0)
      TPR <- ifelse(length(test_pres)>0, sum(test_pres==1)/length(test_pres), NA)
      TNR <- ifelse(length(bg_abs)>0, sum(bg_abs==1)/length(bg_abs), NA)
      tss_vec[i] <- if(!is.na(TPR) && !is.na(TNR)) TPR + TNR - 1 else NA
    }
  }
  
  # Boyce (sample predictions over landscape)
  boyce_vec[i] <- NA
  if(length(pred_test_vals) > 5){
    samp_vals <- tryCatch({
      samp_pts <- terra::spatSample(env_terra, size = 5000, method = "random", xy = TRUE, as.points = TRUE)
      samp_tab <- raster::extract(env_raster, coordinates(as(samp_pts, "Spatial")))
      samp_tab <- samp_tab[complete.cases(samp_tab), , drop = FALSE]
      if(nrow(samp_tab) == 0) return(NULL)
      predict(me, samp_tab)
    }, error = function(e) NULL)
    if(!is.null(samp_vals)){
      boyce_obj <- tryCatch({ ecospat.boyce(fit = samp_vals, obs = pred_test_vals, nclass = 0, window.w = "default") }, error = function(e) NULL)
      if(!is.null(boyce_obj) && !is.null(boyce_obj$cor)) boyce_vec[i] <- boyce_obj$cor
    }
  }
  
  cat(sprintf("Fold %d: AUC=%.3f, TSS=%.3f, Boyce=%.3f\n", i,
              ifelse(is.na(auc_vec[i]),NA,auc_vec[i]),
              ifelse(is.na(tss_vec[i]),NA,tss_vec[i]),
              ifelse(is.na(boyce_vec[i]),NA,boyce_vec[i])))
  fold_end <- Sys.time()
  cat(sprintf("Fold %d finished. Duration: %s\n", i, format(difftime(fold_end, fold_start, units = "mins"))))
  cat(sprintf("Fold %d finished. Duration: %s\n", i, format(difftime(fold_end, fold_start, units = "mins"))), file=logfile, append=TRUE)
}

# --------------- Summarize CV results ---------------------------------------
results_df <- data.frame(
  fold = 1:n_folds,
  n_test = as.integer(sapply(1:n_folds, function(k) sum(folds_vec == k))),
  AUC = auc_vec,
  TSS = tss_vec,
  Boyce = boyce_vec
)
print(results_df)
write.csv(results_df, file.path(output_folder, "results_summary.csv"), row.names = FALSE)
cat("CV summary written to", file.path(output_folder, "cv_results_summary_30s_v6.csv"), "\n", file=logfile, append=TRUE)

# ---------------- Final model on all occurrences -----------------------------
cat("Fitting final MaxEnt on all occurrences...\n"); cat("Fitting final MaxEnt on all occurrences...\n", file=logfile, append=TRUE)
final_mod <- maxent(x = env_raster, p = coordinates(occ_sp), path = output_folder, args = maxent_args_final)

# Final prediction (raster::predict -> convert to SpatRaster)
final_pred_r <- tryCatch({
  predict(env_raster, final_mod)
}, error = function(e){
  message("Final predict (raster::predict) error: ", e$message); cat("Final predict (raster::predict) error: ", e$message, "\n", file=logfile, append=TRUE); return(NULL)
})

if(!is.null(final_pred_r)){
  final_pred_spat <- tryCatch({ rast(final_pred_r) }, error = function(e){
    message("Conversion final prediction to SpatRaster failed: ", e$message); cat("Conversion final prediction to SpatRaster failed: ", e$message, "\n", file=logfile, append=TRUE); return(NULL)
  })
  if(!is.null(final_pred_spat)){
    terra::writeRaster(final_pred_spat, file.path(output_folder, "SDM_Regional_final_30s.tif"), overwrite = TRUE)
    # mean and sd across folds (use preds_list; ensure only non-NULL)
    valid_preds <- preds_list[!sapply(preds_list, is.null)]
    # convert to SpatRaster if any raster-type slipped in
    valid_preds2 <- lapply(valid_preds, function(x){
      if(inherits(x, "Raster")) {
        tryCatch(rast(x), error=function(e) NULL)
      } else if(inherits(x, "SpatRaster")) x else NULL
    })
    valid_preds2 <- valid_preds2[!sapply(valid_preds2, is.null)]
    if(length(valid_preds2) > 0){
      pred_stack <- rast(valid_preds2)
      mean_pred <- terra::app(pred_stack, fun = mean, na.rm = TRUE)
      sd_pred   <- terra::app(pred_stack, fun = sd, na.rm = TRUE)
      terra::writeRaster(mean_pred, file.path(output_folder, "SDM_Regional_mean_30s.tif"), overwrite = TRUE)
      terra::writeRaster(sd_pred,   file.path(output_folder, "SDM_Regional_SD_30s.tif"), overwrite = TRUE)
      cat("Saved mean and SD rasters.\n", file=logfile, append=TRUE)
    } else {
      cat("No per-fold prediction rasters available to compute mean/SD.\n", file=logfile, append=TRUE)
    }
    
    # ---------------- MESS (dismo::mess) --------------------------------
    ref_vals <- raster::extract(env_raster, coordinates(occ_sp))
    ref_vals <- ref_vals[complete.cases(ref_vals), , drop = FALSE]
    if(nrow(ref_vals) > 0){
      mess_rast <- tryCatch({ dismo::mess(env_raster, ref_vals) }, error = function(e){
        message("MESS error: ", e$message); cat("MESS error: ", e$message, "\n", file=logfile, append=TRUE); return(NULL)
      })
      if(!is.null(mess_rast)){
        writeRaster(mess_rast, file.path(output_folder, "SDM_Armenia_local_MESS_30s_v6.tif"), overwrite = TRUE)
        cat("Saved MESS raster.\n", file=logfile, append=TRUE)
      }
    } else {
      cat("No valid training reference values for MESS.\n", file=logfile, append=TRUE)
    }
  }
}

# ---------------- Variable importance & Response curves -----------------------
if(exists("final_mod")){
  pc_rows <- grep("contribution", rownames(final_mod@results), ignore.case = TRUE)
  pi_rows <- grep("permutation.importance", rownames(final_mod@results), ignore.case = TRUE)
  pc <- if(length(pc_rows)) final_mod@results[pc_rows,1] else rep(NA, nlayers(env_raster))
  pi <- if(length(pi_rows)) final_mod@results[pi_rows,1] else rep(NA, nlayers(env_raster))
  vars_df <- data.frame(Variable = names(env_raster), PercentContribution = pc, PermutationImportance = pi)
  write.csv(vars_df, file.path(output_folder, "Variable_Contribution_PC_30s_v6.csv"), row.names = FALSE)
  write.csv(vars_df, file.path(output_folder, "Variable_Importance_PI_30s_v6.csv"), row.names = FALSE)
  
  png(file.path(output_folder, "Response_Curves_30s_v6.png"), width = 1200, height = 900)
  tryCatch({ response(final_mod) }, error = function(e){ plot.new(); text(0.5,0.5,"Response curves error") })
  dev.off()
  cat("Saved PC/PI and Response curves.\n", file=logfile, append=TRUE)
}

cat("Run finished: ", format(Sys.time()), "\n", file=logfile, append=TRUE)
cat("Finished. Results in:", output_folder, "\n")
