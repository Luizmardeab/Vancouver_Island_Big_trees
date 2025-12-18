# Load required libraries
rm(list = ls()); gc()
library(pacman) 

p_load(terra,exactextractr,sf,ggplot2,dplyr,patchwork,tidyverse,nngeo)


setwd('G:/Big_tree/')
# Prediction from Model with all Sentinel 1 and 2 bands
#Our prediction
prod_path = 'Pred_Smooth_S12_10GLOBMAE_ensemble.tif'
Prediction <- terra::rast(prod_path)
plot(Prediction)

#Water mask
water_mask_path = '/water_mask.tif'
water_mask <- terra::rast(water_mask_path)[[4]]
plot(water_mask)

#Study area
shapepath<-"Van_Is.shp"
vc<-terra::vect(shapepath)
vc$AOI=1

#Canopy cover model
ccover_path = "Z:/CDFCP/CC&CH/CCvoer.tif"

r_clamped <- clamp(ccover, lower = 0, upper = 100, values=TRUE)
plot(r_clamped)

r_clamped = crop(r_clamped, ext(Prediction))
plot(r_clamped)

gc()
# Create an empty raster template (extent + resolution)
template <- rast(ext(Prediction), res = 10, crs = "EPSG:3005")
r1 <- rasterize(vc, template, field = "AOI") 
plot(r1)

# Make sure all rasters match r1 resolution/projection before masking
r <- mask(Prediction[[1]], r1, maskvalues = NA, updatevalue = NA) # keep only r1==1
plot(r)

# Make sure all rasters match r1 resolution/projection before masking
water_mask = crop(water_mask,ext(r1))
wm <- mask(water_mask, r1, maskvalues = NA, updatevalue = NA) # keep only r1==1
plot(wm)

#Keep only validy canopy cover
ccover <- mask(r_clamped, r1, maskvalues = NA, updatevalue = NA) # keep only r1==1
plot(ccover)

# Save Vancouver Island Mask, water Mask, Canopy Height, and Cover.
getwd()
pathout <- file.path(getwd(), "msk_rast")

# Create output folder if it does not exist
dir.create(pathout, recursive = TRUE, showWarnings = FALSE)

# File paths
vi_msk <- file.path(pathout, "VI_Mask.tif")
wt_msk <- file.path(pathout, "water_Mask.tif")
ch_msk <- file.path(pathout, "Cheight.tif")    
cc_msk <- file.path(pathout, "CCover.tif")

# Write rasters
writeRaster(r1,      vi_msk, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(wm,      wt_msk, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW"))
writeRaster(r,       ch_msk, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW"))   
writeRaster(ccover,  cc_msk, filetype = "GTiff", overwrite = TRUE, gdal = c("COMPRESS=LZW"))

### Clear memory ###
rm(list = ls()); gc()
###
path <- file.path(getwd(), "msk_rast")
vi_msk <- rast(file.path(path, "VI_Mask.tif"))
wt_msk <- rast(file.path(path, "water_Mask.tif"))
ch_msk <- rast(file.path(path, "Cheight.tif"))    
fc_msk <- rast(file.path(path, "FCover.tif"))

## BEC #
BEC = terra::vect("BEC_VI.shp")
crs(BEC)==crs(fc_msk)

fc_msk = resample(fc_msk,ch_msk)
# Create an empty raster template (extent + resolution)
template <- rast(ext(fc_msk), res = 10, crs = crs(fc_msk))
r1 <- rasterize(BEC, template, field = "ZONE") 
plot(r1)

cats(r1)
# Keep mask: remove CMAunp
keep_mask <- terra::ifel(r1 %in% c(0, 2, 3), 1, NA)

vi_msk <- ifel(vi_msk == 1, 1, NA)
fc_msk <- ifel(fc_msk == 1, 1, NA)
keep_mask <- ifel(keep_mask == 1, 1, NA)

# Start from canopy height raster
prediction_band1 <- ch_msk

# Apply masks to our prediction
prediction_band1 <- mask(prediction_band1, wt_msk,   maskvalue = 1)  # remove water
prediction_band1 <- mask(prediction_band1, vi_msk) # keep only VI
prediction_band1 <- mask(prediction_band1, fc_msk)  # remove FC=NA
prediction_band1 <- mask(prediction_band1, keep_mask) #remove KeepM=NA

# Remove values < 5
prediction_band1 <- terra::ifel(prediction_band1 < 5, NA, prediction_band1)

plot(prediction_band1)

hist(values(prediction_band1)[!is.na(values(prediction_band1))], main = "Histogram of Predicted Canopy Height", xlab = "Values", col = "lightblue", breaks = 50)

gc()
# Count the number of valid (non-NA) pixels
valid_pixel_count <- sum(!is.na(values(prediction_band1)))

# Print the result
print(paste("Number of valid pixels:", valid_pixel_count))

# Get the top 1% highest values
top10 <- quantile(values(prediction_band1), 0.90, na.rm = TRUE)
top5 <- quantile(values(prediction_band1), 0.95, na.rm = TRUE)
top1 <- quantile(values(prediction_band1), 0.99, na.rm = TRUE)

################################################################################
## Create Big-treed forest maps
outdir <- "Big_tree_raster_ex_CMA"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "top_trees"), recursive = TRUE, showWarnings = FALSE)

thresholds <- c(top10, top5, top1)

for (i in thresholds) {
  # --- 1. Create binary raster: canopy > threshold 
  top <- prediction_band1 > i
  top <- classify(top, cbind(NA, 0), include.lowest = TRUE) # ensure 0/1 output
  names(top) <- paste0("Top_", i)
  plot(top, main = paste("Top", i, "percent"))
  
  # --- 2. Write result
  out_path <- file.path(outdir, "top_trees", paste0("Top", i, ".tif"))
  writeRaster(top, out_path, overwrite = TRUE)
}

## BEC #
BEC = terra::vect("BEC_VI.shp")
crs(BEC)==crs(prediction_band1)

# Create an empty raster template (extent + resolution)
template <- rast(ext(prediction_band1), res = 10, crs = crs(prediction_band1))
r1 <- rasterize(BEC, template, field = "ZONE") 
plot(r1)
r2 = rasterize(BEC, template, field = "MAP_LABEL") 
plot(r2)

gc()

### Create big-treed forest maps with defintion set by BEC zones
# Table S4
## BEC Threshold
z <- zonal(prediction_band1, r1, fun = quantile, probs = c(0.90,0.95,0.99), na.rm = TRUE)
lut_r1 <- as.data.frame(cats(r1)[[1]])
names(lut_r1)[1] <- "ID"
z <- z%>%left_join(lut_r1, by=c("ZONE"))
z$top10 = z$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,1]
z$top5 = z$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,2]
z$top1 = z$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,3]
z <- dplyr::select(z, dplyr::all_of(c( "ID","ZONE", "top10", "top5", "top1")))

hvals <- values(prediction_band1)[values(r1) == 2]
hvals <- hvals[!is.na(hvals)]-10
hist(hvals, breaks = 50)

gc()

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "top_trees_BEC"), recursive = TRUE, showWarnings = FALSE)

for(i in colnames(z)[3:5]){
  
  # --- 2. Create raster with matching percentile value at each pixel
  prast <- classify(r1, cbind(z$ID, z[[i]]))  # <- use [[i]] to extract vector
  
  # --- 3. Create binary raster: canopy > regional percentile
  top <- prediction_band1 > prast
  top <- classify(top, cbind(NA, 0), include.lowest = TRUE) # ensure 0/1 output
  names(top) <- paste(i,"pct")
  
  # --- 4. Write result
  out_path <- file.path(outdir, "top_trees_BEC", paste0(i, ".tif"))
  writeRaster(top, out_path, overwrite = TRUE)
  gc()
}


# Create big-treed forest maps with defitinion set for each BEC variant
# Table S4
## RDs Threshold
zrd <- zonal(prediction_band1, r2, fun = quantile, probs = c(0.90,0.95,0.99), na.rm = TRUE)
lut_r1 <- as.data.frame(cats(r2)[[1]])
zrd <- zrd%>%left_join(lut_r1, by=c("MAP_LABEL"))
zrd$top10 = zrd$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,1]
zrd$top5 = zrd$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,2]
zrd$top1 = zrd$Pred_Smooth_S12_10GLOBMAE_ensemble_1[,3]

names(zrd)[3] = "ID"
zrd <- dplyr::select(zrd, dplyr::all_of(c("ID", "MAP_LABEL", "top10", "top5", "top1")))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "top_trees_variant"), recursive = TRUE, showWarnings = FALSE)

for(i in colnames(zrd)[3:5]){
  prast <- classify(r2, cbind(zrd$ID, zrd[[i]]))  # <- use [[i]] to extract vector
  top <- prediction_band1 > prast
  top <- classify(top, cbind(NA, 0), include.lowest = TRUE) # ensure 0/1 output
  names(top) <- paste(i,"pct")
  
  #Write result
  out_path <- file.path(outdir, "top_trees_variant", paste0(i, ".tif"))
  writeRaster(top, out_path, overwrite = TRUE)
  gc()
}

############### Create Raster Sum ###################
# used in Figure 3
# 1. List all folders directly inside outdir
folders <- list.dirs(outdir, full.names = TRUE, recursive = FALSE)

# 2. For each folder, list all .tif files
tif_files <- unlist(lapply(folders, function(f) {
  list.files(f, pattern = "\\.tif$", full.names = TRUE)
}))

# Check
print(tif_files)

# 3. Stack all rasters together
r_stack <- rast(tif_files)

# 4. Sum all layers
r_sum <- sum(r_stack, na.rm = TRUE)
plot(r_sum)

r_sum <- mask(r_sum, r, maskvalues = NA, updatevalue = NA) # keep only r1==1
plot(r_sum)
# 5. Save output
out_file <- file.path(outdir, "sum_all_rasters.tif")
writeRaster(r_sum, out_file, overwrite = TRUE)

gc()

# r_sum is your summed raster
vals <- freq(r_sum, digits = 0)   # digits=0 ensures integer grouping

print(vals)

vals$percent <- 100 * vals$count / sum(vals$count)

vals = vals%>%filter(value!=0)%>%mutate(Perc_Ov = 100 * count / sum(count))

print(vals)

# List all folders inside outdir
folders <- list.dirs(outdir, full.names = TRUE, recursive = FALSE)

for (f in folders) {
  
  message("Processing folder: ", f)
  
  # List all .tif files in the folder
  tif_files <- list.files(f, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(tif_files) == 0) {
    message("  No tif files found, skipping.")
    next
  }
  
  # Load rasters
  r_stack <- rast(tif_files)
  
  # Pixel-wise sum
  r_sum <- sum(r_stack, na.rm = TRUE)
  
  # Output file name
  out_file <- file.path(f, "sum_raster.tif")
  
  writeRaster(r_sum, out_file, overwrite = TRUE)
  gc()
  message("  Saved: ", out_file)
}

## Create big-treed forest maps for each alternative canopy height product ####
##############  Canopy height Products ##########################################
prod_path = "Height_Products.tif"
Hghgt_prod<-terra::rast(prod_path)

Hghgt_prod <- crop(Hghgt_prod, ext(fc_msk))
Hghgt_prod <- mask(Hghgt_prod, fc_msk, maskvalues = NA, updatevalue = NA)
Hghgt_prod[[5]] = Hghgt_prod[[5]]*100
plot(Hghgt_prod[[5]])
get_top_percentiles <- function(band, water_mask, min_height = 5) {
  
  # Remove water (set to NA)
  v <- values(band)
  v[values(water_mask) == 1] <- NA
  
  # Remove very small heights
  v[v < min_height] <- NA  
  
  # Remove NA before computing percentiles
  v_clean <- v[!is.na(v)]
  
  if (length(v_clean) == 0) {
    return(c(top10 = NA, top5 = NA, top1 = NA))
  }
  
  # Compute thresholds
  t10 <- quantile(v_clean, 0.90, na.rm = TRUE)
  t5  <- quantile(v_clean, 0.95, na.rm = TRUE)
  t1  <- quantile(v_clean, 0.99, na.rm = TRUE)
  
  return(c(top10 = t10, top5 = t5, top1 = t1))
}

bands_to_extract <- c(3, 4, 5)

results <- lapply(bands_to_extract, function(i) {
  get_top_percentiles(Hghgt_prod[[i]], wt_msk)
})


names(results) <- paste0("Band_", bands_to_extract)


# Thresholds per band
cutoffs <- list(
  Band_3 = c(top10 = 34, top5 = 36, top1 = 40),
  Band_4 = c(top10 = 42, top5 = 45, top1 = 49),
  Band_5 = c(top10 = 32, top5 = 34, top1 = 39)
)

# Output directory
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "top_trees_GProd"), recursive = TRUE, showWarnings = FALSE)

# Loop over bands and cutoffs
for(band_idx in 3:5){
  band_name <- paste0("Band_", band_idx)
  thresholds <- cutoffs[[band_name]]
  
  band_rast <- Hghgt_prod[[band_idx]]  # Select the band
  
  for(th_name in names(thresholds)){
    thresh <- thresholds[[th_name]]
    
    # Create binary raster
    top_rast <- band_rast > thresh
    top_rast <- classify(top_rast, cbind(NA, 0), include.lowest = TRUE)  # ensure 0/1
    
    # Set name
    names(top_rast) <- paste0(band_name, "_", th_name)
    
    # Save raster
    out_path <- file.path(outdir, "top_trees_GProd", paste0(names(top_rast), ".tif"))
    writeRaster(top_rast, out_path, overwrite = TRUE)
    gc()
    print(paste("Saved:", out_path))
  }
}

#################################################################
######### VRI top 1, 5, and 10% ##################
VRI = vect("VRI_VI.shp")
names(VRI)

template <- rast(ext(fc_msk), res = 10, crs = crs(fc_msk))
VRI <- rasterize(VRI, template, field = "PROJ_HEIGH") 
plot(VRI)
VRI <-  mask(VRI, fc_msk, maskvalues = NA, updatevalue = NA)
# Get the top 1% highest values
top10_threshold <- quantile(values(VRI), 0.90, na.rm = TRUE)
top5_threshold <- quantile(values(VRI), 0.95, na.rm = TRUE)
top1_threshold <- quantile(values(VRI), 0.99, na.rm = TRUE)

outdir <- "Big_tree_raster_ex_CMA"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "top_trees_VRI"), recursive = TRUE, showWarnings = FALSE)

thresholds <- c(top10_threshold, top5_threshold, top1_threshold)

for (i in thresholds) {
  # --- 1. Create binary raster: canopy > threshold
  top <- VRI > i
  top <- classify(top, cbind(NA, 0), include.lowest = TRUE) # ensure 0/1 output
  names(top) <- paste0("Top_", i)
  plot(top, main = paste("Top", i, "percent"))
  
  # --- 2. Write result
  out_path <- file.path(outdir, "top_trees_VRI", paste0("Top", i, "_canopy.tif"))
  writeRaster(top, out_path, overwrite = TRUE)
  gc()
}
