# Load required libraries
rm(list = ls()); gc()
library(pacman) 

p_load(terra,exactextractr,sf,ggplot2,dplyr,patchwork,tidyverse,nngeo)

setwd('G:/Big_tree/')

# 1. List all folders directly inside outdir
foldir = "Big_tree_raster_ex_CMA"
folders <- list.dirs(foldir, full.names = TRUE, recursive = FALSE)

# 2. For each folder, list all .tif files
tif_files <- unlist(lapply(folders, function(f) {
  list.files(f, pattern = "\\.tif$", full.names = TRUE)
}))

# Check
print(tif_files)

##############################################################################
## BEC #
BEC = terra::vect("BEC_VI.shp")
crs(BEC)==crs(prediction_band1)

# Create an empty raster template (extent + resolution)
template <- rast(ext(prediction_band1), res = 10, crs = crs(prediction_band1))
r1 <- rasterize(BEC, template, field = "ZONE") 
plot(r1)

r2 = rasterize(BEC, template, field = "MAP_LABEL") 
plot(r2)

###################### Calculate spatial overlap ###############################
# List all folders inside outdir
cats(r1)
keep_mask <- terra::ifel(r1 %in% c(0, 2, 3), 1, NA) #Remove CMAunp

plot(keep_mask)

compute_overlap <- function(ref_raster, target_raster,keep_mask=keep_mask){
  r_ref <- rast(ref_raster)
  r_target <- rast(target_raster)

  
  r_ref    <- mask(r_ref, keep_mask)
  r_target <- mask(r_target, keep_mask)
  
  
  # Create masks for cells == 1
  r_ref_mask <- r_ref[] == 1
  r_target_mask <- r_target[] == 1
  
  # Count overlapping cells (both == 1)
  overlap_cells <- sum(r_ref_mask & r_target_mask, na.rm = TRUE)
  
  # Count total valid cells in reference
  total_ref_cells <- sum(r_ref_mask, na.rm = TRUE)
  
  # Fraction of reference cells that overlap
  overlap_fraction <- overlap_cells / total_ref_cells
  
  # Return overlap count and fraction
  return(c(overlap_cells, overlap_fraction))
}

################################################################################
# Table 1
folders <- list.dirs(foldir, full.names = TRUE, recursive = FALSE)
target_folders = folders[c(1,3,5)]

# Initialize a list to store results
overlap_data <- list()

# Percentages you want to loop over
perc_values <- c(1, 5, 10)

for (j in perc_values){
  
  ref_files <- list.files(target_folders[1], pattern = paste0("top", j, ".tif$"), full.names = TRUE)
  target_files <- list.files(target_folders[2:3], pattern = paste0("top", j, ".tif$"), full.names = TRUE)
  
  # Initialize vectors to store overlap cells and fraction
  overlap_cells <- numeric(4)
  overlap_frac <- numeric(4)
  
  for (i in 1:4){
    ov <- compute_overlap(ref_files, target_files[i], keep_mask=keep_mask)
    overlap_cells[i] <- ov[1]
    overlap_frac[i] <- ov[2]
    gc()
  }
  
  # Store in the list
  overlap_data[[paste0(j, "%")]] <- overlap_cells
  overlap_data[[paste0(j, "%_frac")]] <- overlap_frac
}

# Convert to a table
df <- data.frame(
  perc = rep(names(overlap_data)[c(TRUE, FALSE)], each = 4),
  model = rep(paste0("Model_", c("Potapov","Lang","Pauls","VRI")), length(perc_values)),
  overlap_cells = unlist(overlap_data[seq(1, length(overlap_data), 2)]),
  overlap_fraction = unlist(overlap_data[seq(2, length(overlap_data), 2)])
)

df <- df %>%
  mutate(overlap_info = paste0(round(overlap_cells,2)*100/10000, " (", round(overlap_fraction,2), ")")) %>%
  dplyr::select(perc, model, overlap_info) %>%
  pivot_wider(names_from = model, values_from = overlap_info)%>%data.frame()

df

 
### Overlap between total, BEC, and variant level analysis #####
# Table 2
# List all .tif files 
tif_files <- unlist(lapply(folders[c(1,2,4)], function(x) {
  list.files(x, pattern = "\\.tif$", full.names = TRUE)
}))

tif_files

# Extract group name and the top percentage from filename
file_info <- data.frame(
  path = tif_files,
  group = ifelse(grepl("top_trees_variant", tif_files), "Variant",
                 ifelse(grepl("top_trees_BEC", tif_files), "BEC", "Total")),
  top = sub(".*top([0-9]+).*", "\\1", tif_files),
  stringsAsFactors = FALSE
)

file_info[1:3,3] = c(1, 10, 5)
# Only keep top1, top5, top10
file_info <- file_info[file_info$top %in% c("1","5","10"), ]

# ---- Overlap Function ----
compute_overlap <- function(ref_file, target_file){
  r_ref <- rast(ref_file)
  r_tgt <- rast(target_file)
  
  # Align target to reference
  r_tgt_res <- resample(r_tgt, r_ref, method="near")
  
  # Masks: valid = cell == 1
  m_ref <- r_ref[] == 1
  m_tgt <- r_tgt_res[] == 1
  
  overlap <- sum(m_ref & m_tgt, na.rm=TRUE)
  ref_total <- sum(m_ref, na.rm=TRUE)
  
  return(c(overlap_cells = overlap,
           overlap_fraction = overlap / ref_total))
}

# ---- Run all comparisons ----
results <- data.frame()

for (top_level in c("1","5","10")) {
  
  base_file <- file_info$path[file_info$group=="Total" & file_info$top==top_level]
  
  bec_file      <- file_info$path[file_info$group=="BEC" & file_info$top==top_level]
  variant_file  <- file_info$path[file_info$group=="Variant" & file_info$top==top_level]
  
  # Compare Base vs BEC
  ov1 <- compute_overlap(base_file, bec_file)
  
  # Compare Base vs Variant
  ov2 <- compute_overlap(base_file, variant_file)
  
  # Store results
  results <- rbind(results,
                   data.frame(
                     top = paste0("top", top_level),
                     compare = "Base vs BEC",
                     overlap_cells = ov1["overlap_cells"],
                     overlap_fraction = ov1["overlap_fraction"]
                   ),
                   data.frame(
                     top = paste0("top", top_level),
                     compare = "Base vs Variant",
                     overlap_cells = ov2["overlap_cells"],
                     overlap_fraction = ov2["overlap_fraction"]
                   )
  )
}

results%>%data.frame()


#### BEC variant representation for each BIg tree definition and analysis level
### Overlap between total, BEC, and variant level analysis #####
## Figure S3
# List all .tif files 
# outdir <- "Chapter_3/Big_tree_raster_ex_CMA"
folders <- list.dirs(outdir, full.names = TRUE, recursive = FALSE)

tif_files <- unlist(lapply(folders[c(1,2,4)], function(x) {
  list.files(x, pattern = "\\.tif$", full.names = TRUE)
}))

tif_files

# Extract group name and the top percentage from filename
file_info <- data.frame(
  path = tif_files,
  group = ifelse(grepl("top_trees_variant", tif_files), "Variant",
                 ifelse(grepl("top_trees_BEC", tif_files), "BEC", "Total")),
  top = sub(".*top([0-9]+).*", "\\1", tif_files),
  stringsAsFactors = FALSE
)

# Only keep top1, top5, top10
file_info <- file_info[file_info$top %in% c("1","5","10"), ]


r2 = rasterize(BEC, template, field = "MAP_LABEL")
plot(r2)


# ---- Run all BEC variant representation ----
results <- data.frame()

for (top_level in c("1","5","10")) {
  for (type in c("Total", "BEC", "Variant")) {
    
    # Load raster
    base_file <- file_info$path[file_info$group == type & file_info$top == top_level]
    file <- rast(base_file)
    
    file[file == 0] <- NA
    z <- zonal(file, r2, fun = "sum", na.rm = TRUE)
    names(z)[2] <- "Count"
    
    # Store results
    results <- rbind(
      results,
      data.frame(
        top      = paste0("top_", top_level),
        Level    = type,
        BEC_Var  = z$MAP_LABEL,
        Count    = z$Count,
        Area     = z$Count * 100 / 10000000
      )
    )
    
    # Explicit cleanup
    rm(file, z)
    gc()
  }
}


results

df_perc <- results %>%
  group_by(top, Level) %>%
  mutate(
    TotalArea = sum(Area, na.rm = TRUE),
    Perc = round((Area / TotalArea) * 100, 2)   # <- put 2 *inside* round()
  ) %>%
  ungroup() %>%
  data.frame()


# Prepare data for plotting
plot_data <- df_perc %>%filter(BEC_Var!="CMAunp")%>%
  group_by(Level, BEC_Var) %>%
  reframe(
    Perc_top5 = Area[top == "top_5"],           # bar value
    Perc_min = min(Area[top %in% c("top_1", "top_10")]), # lower whisker
    Perc_max = max(Area[top %in% c("top_1", "top_10")])  # upper whisker
  ) %>%
  ungroup()%>%data.frame()

plot_data%>%filter(Level == "Total")

# Plot
B = plot_data %>%
  mutate(
    Level = factor(Level, levels = c("Total", "BEC", "Variant")),
    Level = recode(Level,
                   "BEC" = "BEC zone",
                   "Variant" = "BEC variant")
  )%>%
ggplot(aes(x = BEC_Var, y = Perc_top5, fill =Level, group = Level)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Perc_min, ymax = Perc_max), 
                width = 0.3, color = "black", 
                position = position_dodge(width = 0.8), alpha=0.75) +
  scale_fill_manual(values = c("#1E4F26","#2F6DAD","#9AA6B4")) +
  # facet_wrap(~Level)+
  labs(x = "BEC Variant", y = "Area (1000 ha)", fill = "BEC Variant") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.position=c(0.75,0.15),plot.margin = unit(c(0,30,0,0), "pt"),
                 strip.text.x = element_text( size = 12, color = "black", face = "bold"),
                 strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
                 legend.direction='vertical', legend.title=element_text( size = 12, color = "black", face = "bold"),
                 legend.background = element_rect(fill="transparent", linetype="solid", linewidth=0.5, colour="black"),
                 legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
                 legend.text=element_text(size=12, colour="black"),
                 panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
                 panel.background=element_rect(fill="white", colour="black", linewidth=1, linetype="solid"),
                 panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
                 panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
                 plot.title=element_text(color="black", size=14, hjust=0,face='bold'),
                 axis.title.x=element_text(size=12, colour="black"),#axis.ticks.x=element_blank(),
                 axis.text.x=element_text(size=12, colour="black"),
                 axis.title.y=element_text(size=12, colour="black"),
                 axis.text.y=element_text(size=12, colour="black"), legend.key=element_rect(fill="transparent", colour="transparent"))

ggsave(paste("Paper_Images/BEC_ForestCover_Big_treed.jpeg",sep=""), height=7, width=9, units="in", dpi=300)
## Get total forest cover per BEC variant
# count pixels per MAP_LABEL
r2 = mask(r2, fc_msk)
plot(r2)
tab <- as.data.frame(freq(r2, digits = 0))
tab$area_ha <- tab$count * 100/10000

total_area <- sum(tab$area_ha, na.rm = TRUE)

tab$proportion <- tab$area_ha / total_area*100

A = tab %>%filter(value!="CMAunp")%>%
  ggplot(aes(x = value, y = area_ha/1000)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8), alpha=0.75, fill="#1E4F26") +
  labs(x = "BEC Variant", y = "Area (1000 ha)", fill = "BEC Variant") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(legend.position=c(0.75,0.15),plot.margin = unit(c(0,30,0,0), "pt"),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text( size = 12, color = "black", face = "bold"),
        legend.background = element_rect(fill="transparent", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, linewidth=0.5),
        legend.text=element_text(size=12, colour="black"),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="black", linewidth=1, linetype="solid"),
        panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'),
        axis.title.x=element_text(size=12, colour="black"),#axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=12, colour="black"),
        axis.title.y=element_text(size=12, colour="black"),
        axis.text.y=element_text(size=12, colour="black"), legend.key=element_rect(fill="transparent", colour="transparent"))


final_plot <- (A|B)+
  plot_annotation(tag_levels = "A") 

#can change 14 to any desired size


ggsave(paste("Chapter_3/Paper_Images/BIg_tree_zones.jpeg",sep=""), height=7, width=12, units="in", dpi=300)


################################################################################
## Canopy height and AGE
# Figure 4
VRI = vect("VRI_VI.shp")
names(VRI)
hist(VRI$PROJ_AGE_1)
VRI$STAge <- ifelse(
  VRI$PROJ_AGE_1 <= 20, "Unaged",
  ifelse(VRI$PROJ_AGE_1 > 40 & VRI$PROJ_AGE_1 < 80, "Young",
         ifelse(VRI$PROJ_AGE_1 >= 80 & VRI$PROJ_AGE_1 < 250, "Mature",
                ifelse(VRI$PROJ_AGE_1 >= 250, "Old", NA)))
)


VRI$STAge <- factor(VRI$STAge, levels = c("Unaged","Young", "Mature", "Old"))
# Add numeric code column
VRI$STAge_code <- as.numeric(VRI$STAge)  # if you want 0,1,2
table(VRI$STAge,VRI$STAge_code)

template <- rast(ext(fc_msk), res = 10, crs = crs(fc_msk))
AGE <- rasterize(VRI, template, field = "STAge_code") 
plot(AGE)
AGE <-  mask(AGE, fc_msk, maskvalues = NA, updatevalue = NA)
plot(AGE)
gc()
freq(fc_msk)

cats <- data.frame(
  value = c(1, 2, 3, 4),
  label = c("Unaged","Young", "Mature", "Old")
)

levels(AGE) <- cats

cats(AGE)
cats(Btres_cat)
gc()
# Create a combined raster code
### Figure 4
Btres = terra::rast("Chapter_3/Big_tree_raster_ex_CMA/sum_raster_top.tif")
Btres_cat <- classify(Btres_cat, rcl = matrix(c(0, 0, NA), ncol = 3, byrow = TRUE))

combined <- (AGE - 1) * 10 + Btres_cat
plot(combined)

# Make sure raster is integer
combined <- as.factor(combined)

cats(combined)
f <- freq(combined)
f = f%>%filter(value!=0)
f$perc <- f$count / sum(f$count) * 100
f

# Define categories matching the actual raster values
cats_df <- data.frame(
  value = c(0,1,2,3,10, 11,12,13,20,21,22,23),
  label = c(
    NA, "Young_Top10","Young_Top5","Young_Top1",
    NA,"Mature_Top10","Mature_Top5","Mature_Top1",
    NA,"Old_Top10","Old_Top5","Old_Top1"
  )
)

levels(combined) <- cats_df
rcl <- matrix(c(
  1,  1,
  2,  2,
  3,  3,
  11, 4,
  12, 5,
  13, 6,
  21, 7,
  22, 8,
  23, 9
), ncol = 2, byrow = TRUE)

combined <- classify(combined, rcl, others = NA)

# Check categories
combined = as.factor(combined)
cats(combined)

plot(combined)
# Save the classified raster
writeRaster(combined, "Chapter_3/Big_tree_raster_ex_CMA/Combined_9class.tif", overwrite=TRUE)


################### PAs ###################################
# Table S5
PA_OECM = st_read("VI_PA.shp")
if(crs(PA_OECM)!=crs(ch_msk)){
  PA_OECM <- st_transform(PA_OECM, crs(ch_msk))
}

template <- rast(ext(fc_msk), res = 10, crs = crs(fc_msk))
PAs <- rasterize(PA_OECM, template, field = "PA_OECM") 
plot(PAs)
PAs <-  mask(PAs, fc_msk, maskvalues = NA, updatevalue = NA)
plot(PAs)
PAs <- terra::ifel(is.na(PAs), 2, PAs)

cats(PAs)

# 1. List all folders directly inside outdir
folders <- list.dirs(outdir, full.names = TRUE, recursive = FALSE)

# 2. For each folder, list all .tif files
tif_files <- unlist(lapply(folders, function(f) {
  list.files(f, pattern = "\\.tif$", full.names = TRUE)
}))

# Check
print(tif_files)

# Extract group name and the top percentage from filename
file_info <- data.frame(
  path = tif_files,
  group = ifelse(grepl("top_trees_variant", tif_files), "Variant",
                 ifelse(grepl("top_trees_BEC", tif_files), "BEC",
                        ifelse(grepl("Band_3", tif_files), "Potapov",
                               ifelse(grepl("Band_4", tif_files), "Lang",
                                      ifelse(grepl("Band_5", tif_files), "Pauls",
                                             ifelse(grepl("top_trees_VRI", tif_files), "VRI", "Total")))))),
  top = sub(".*top([0-9]+).*", "\\1", tif_files),
  stringsAsFactors = FALSE
)


# Only keep top1, top5, top10
file_info <- file_info[file_info$top %in% c("1","5","10"), ]


results <- data.frame()

for (top_level in c("1","5","10")) {
  for (type in c("Total", "BEC", "Variant", "Potapov", "Lang", "Pauls", "VRI")) {
    
    # Get the file
    base_file <- file_info$path[file_info$group == type & file_info$top == top_level]
    
    # Skip missing files
    if (length(base_file) == 0) {
      message("No file found for: type = ", type, " | top = ", top_level)
      next
    }
    
    # Load raster
    file <- rast(base_file)
    file[file == 0] <- NA
    
    # Calculate zonal stats
    z <- zonal(file, PAs, fun = "sum", na.rm = TRUE)
    names(z)[2] <- "Count"
    
    # Convert pixel count ??? hectares
    z$Area <- z$Count * 0.01   # 0.01 ha per 10 m pixel
    totalA = sum(z$Area)
    z$Perc = z$Area/totalA*100
    z$PA_OECM[3] = "Not Protected"
    # Append results
    results <- rbind(
      results,
      data.frame(
        top     = paste0("top_", top_level),
        Level   = type,
        PA_OECM = z$PA_OECM,
        Count   = z$Count,
        Area    = z$Area,
        Perc    = z$Perc
      )
    )
    
    rm(file, z)
    gc()
  }
}

results

df <- results %>%mutate(Perc = Area )%>%
  select(top, Level, PA_OECM, Perc) %>%
  mutate(
    top_n = as.numeric(str_remove(top, "top_"))
  ) %>%
  group_by(Level, PA_OECM) %>%
  summarise(
    main = Perc[top_n == 5],
    min  = min(Perc[top_n %in% c(1,10)]),
    max  = max(Perc[top_n %in% c(1,10)]),
    .groups = "drop"
  ) %>%
  mutate(
    value = sprintf("%.1f (%.1f-%.1f)", main, min, max)
  ) %>%
  select(Level, PA_OECM, value) %>%
  pivot_wider(
    names_from = PA_OECM,
    values_from = value
  )%>%data.frame()

df

# Results presented in the text
########################### HFP ################################
HFP = rast('G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/HFP_3005.tif')
Btres = terra::rast("Chapter_3/Big_tree_raster_ex_CMA/top_trees/top5.tif")
HFP = resample(HFP, Btres)

HFP_class <- mask(HFP,  Btres, maskvalue = 0, updatevalue = NA)
plot(HFP_class)
gc()

# Reclassify using nested ifel
HFP_class <- classify(HFP_class, rcl = matrix(c(
  -Inf, 1, 1,      # 1 = Wilderness
  1, 4, 2,        # 2 = Low
  4, 10, 3,       # 3 = Moderate
  10, Inf, 4       # 4 = Highly Modified
), ncol = 3, byrow = TRUE))

# Assign meaningful category names
levels(HFP_class ) <- data.frame(
  ID = 1:4,
  category = c("Wilderness", "Low", "Moderate", "Highly Modified")
)

HFP_class <- as.factor(HFP_class)

# Check
cats(HFP_class)
plot(HFP_class)
gc()


f <- freq(HFP_class)
f = f%>%filter(value!=0)
f$perc <- f$count / sum(f$count) * 100
f

freq(Btres)
#################### HFP vs AGE ###########################
VRI = vect("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_2/VRI_VI.shp")
names(VRI)
hist(VRI$PROJ_AGE_1)
VRI$STAge <- ifelse(
  VRI$PROJ_AGE_1 <= 20, "Unaged",
  ifelse(VRI$PROJ_AGE_1 > 40 & VRI$PROJ_AGE_1 < 80, "Young",
         ifelse(VRI$PROJ_AGE_1 >= 80 & VRI$PROJ_AGE_1 < 250, "Mature",
                ifelse(VRI$PROJ_AGE_1 >= 250, "Old", NA)))
)


VRI$STAge <- factor(VRI$STAge, levels = c("Unaged","Young", "Mature", "Old"))

# Add numeric code column
VRI$STAge_code <- as.numeric(VRI$STAge)  # if you want 0,1,2
table(VRI$STAge,VRI$STAge_code)

template <- rast(ext(Btres), res = 10, crs = crs(fc_msk))
AGE <- rasterize(VRI, template, field = "STAge_code") 
plot(AGE)

cats <- data.frame(
  value = c(1, 2, 3, 4),
  label = c("Unaged","Young", "Mature", "Old")
)

levels(AGE) <- cats

cats(AGE)
cats(HFP_class)

combined_age = (HFP_class - 1) * 4 + AGE

gc()

combined_age <- as.factor(combined_age)
cats(combined_age)
f <- freq(combined_age)
f$perc <- f$count / sum(f$count) * 100
f


############ BIG trees within PAs ##############################
Btres = terra::rast("Chapter_3/Big_tree_raster_ex_CMA/top_trees/top5.tif")
PA_OECM = vect("G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Chapter_3/VI_PA.shp")
if(crs(PA_OECM)!=crs(Btres)){
  PA_OECM <- st_transform(PA_OECM, crs(Btres))
}

# Keep only value == 1 (others become NA)
Btres1 <- Btres
Btres1[Btres1 != 1] <- NA

# Count pixels (sum of 1s) inside each polygon
pixel_counts <- terra::extract(
  Btres1,
  PA_OECM,
  fun = sum,
  na.rm = TRUE
)

# Attach counts back to polygon attribute table
PA_OECM$Btres_pixels <- pixel_counts[, 2]
PA_OECM$Btres_area_ha <- PA_OECM$Btres_pixels *100/10000
# Inspect result
head(PA_OECM)

PA_OECM$ID = c(1:nrow(PA_OECM))
PA_OECM$PA_A <- expanse(PA_OECM, unit = "ha")


PA_OECM %>%
  data.frame() %>%        # important if PA_OECM is sf
  group_by(NAME_E) %>%
  summarise(Area = sum(Btres_area_ha, na.rm = TRUE)) %>%
  arrange(desc(Area)) %>%          # largest area first
  data.frame()

names(PA_OECM)
unique(PA_OECM$PA_OECM)

PA_OECM %>%
  data.frame() %>% 
  filter(PA_OECM == "Protected Area") %>%
  group_by(NAME_E) %>%
  summarise(Area = sum(Btres_area_ha, na.rm = TRUE), .groups = "drop") %>%
  filter(Area > 1) %>%
  nrow()


PA_OECM %>%
  data.frame() %>%
  filter(PA_OECM == "OECM") %>%
  group_by(ID) %>%
  summarise(PA_a = sum(PA_A, na.rm = TRUE), Area = sum(Btres_area_ha, na.rm = TRUE)) %>%
  filter(PA_a >1, Area>1) %>%
  nrow()


## Harvested areas #############
Btres = terra::rast("top5.tif")
CBlocks = vect("VI_CCutBlocks.shp")
if(crs(CBlocks)!=crs(Btres)){
  CBlocks <- st_transform(CBlocks, crs(Btres))
}

# Keep only value == 1 (others become NA)
Btres1 <- Btres
Btres1[Btres1 != 1] <- NA

# Count pixels (sum of 1s) inside each polygon
pixel_counts <- terra::extract(
  Btres1,
  CBlocks,
  fun = sum,
  na.rm = TRUE
)

# Attach counts back to polygon attribute table
CBlocks$Btres_pixels <- pixel_counts[, 2]
CBlocks$Btres_area_ha <- CBlocks$Btres_pixels *100/10000

CBlocks%>%data.frame()%>%
  filter(HARVEST__1>=2018)%>%summarise(Havested = sum(Btres_area_ha, na.rm = TRUE))

