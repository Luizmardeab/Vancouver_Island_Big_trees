rm(list = ls()); gc()
library(pacman) 

p_load(raster,terra,exactextractr,sf,ggplot2,dplyr,patchwork,tidyverse,nngeo,spatialEco, scales)

setwd('G:/Big_tree/')

# ##################################################################################
mytheme =  theme(legend.position="bottom",plot.margin = unit(c(0,30,0,0), "pt"),
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


############################################################################
# Loading all environemtnal descriptors 
#DEM
DEM<-terra::rast('DEM_3005.tif')

if(crs(Prediction)!=crs(DEM)){
  DEM<-terra::project(DEM, crs(Prediction))
}

#Wetness index
Wet_ind<-terra::rast('Wet_ind_3005.tif')
if(crs(Prediction)!=crs(Wet_ind)){
  Wet_ind<-terra::project(Wet_ind, crs(Prediction))
}

#Climate
clim<-terra::rast("clim_3005.tif")
# c('AHM.asc', 'bFFP.asc', 'CMD.asc', 'CMI.asc',"DD_0.asc",'DD_18.asc','DD18.asc', 'DD5.asc',   'eFFP.asc',
#   'EMT.asc', 'Eref.asc', 'EXT.asc', 'FFP.asc', 'MAP.asc', 'MAT.asc', 'MCMT.asc', 'MSP.asc', 'MWMT.asc',
#   'NFFD.asc', 'PAS.asc', 'RH.asc', 'SHM.asc', 'TD.asc')

if(crs(Prediction)!=crs(clim)){
  clim<-terra::project(clim, crs(Prediction))
}

#Human Footprint
HFP = terra::rast("HFP_3005.tif")
if(crs(Prediction)!=crs(HFP)){
  HFP<-terra::project(HFP, crs(Prediction))
}

#Topographical Position Index
TPI = terra::rast("TPI.tif")

# #Slope 
slope = terra::rast("slope.tif")

# Distance from ocean and Distance from wet areas and streams
DFO = terra::rast("dist_Ocean.tif")
DWS = terra::rast("dist_WStream.tif")

# Check CRS of all rasters
all_rasters <- list(Prediction, DEM,slope, Wet_ind,TPI, clim, DWS, HFP,DFO)
crs_check <- sapply(all_rasters, crs)

# Extract the first band/layer from the Prediction SpatRaster
path <- file.path(getwd(), "msk_rast")
vi_msk <- rast(file.path(path, "VI_Mask.tif"))
wt_msk <- rast(file.path(path, "water_Mask.tif"))
ch_msk <- rast(file.path(path, "Cheight.tif"))     # fixed name
fc_msk <- rast(file.path(path, "FCover.tif"))

## Biogeoclimatic zones #
BEC = terra::vect("BEC_VI.shp")
crs(BEC)==crs(fc_msk)

fc_msk = resample(fc_msk,ch_msk)
# Create an empty raster template (extent + resolution)
template <- rast(ext(fc_msk), res = 10, crs = crs(fc_msk))
r1 <- rasterize(BEC, template, field = "ZONE") 
plot(r1)

cats(r1)
# Keep mask: keep only classes 0, 2, 3 → 1, otherwise NA
keep_mask <- terra::ifel(r1 %in% c(0, 2, 3), 1, NA)

vi_msk <- ifel(vi_msk == 1, 1, NA)
fc_msk <- ifel(fc_msk == 1, 1, NA)
keep_mask <- ifel(keep_mask == 1, 1, NA)

# Start from canopy height raster
prediction_band1 <- ch_msk

# Apply masks
prediction_band1 <- mask(prediction_band1, wt_msk,   maskvalue = 1)  # remove water
prediction_band1 <- mask(prediction_band1, vi_msk) # keep only VI
prediction_band1 <- mask(prediction_band1, fc_msk)  # remove FC=NA
prediction_band1 <- mask(prediction_band1, keep_mask) #remove KeepM=NA

# Remove values < 5
prediction_band1 <- terra::ifel(prediction_band1 < 5, NA, prediction_band1)

ch_msk<- prediction_band1

hist(values(prediction_band1)[!is.na(values(prediction_band1))], main = "Histogram of Predicted Canopy Height", xlab = "Values", col = "lightblue", breaks = 50)

# Count the number of valid (non-NA) pixels
valid_pixel_count <- sum(!is.na(values(prediction_band1)))

# Print the result
print(paste("Number of valid pixels:", valid_pixel_count))

# Get the top 1% highest values
top1 <- quantile(values(prediction_band1), 0.99, na.rm = TRUE)
top5<- quantile(values(prediction_band1), 0.95, na.rm = TRUE)
top10 <- quantile(values(prediction_band1), 0.90, na.rm = TRUE)

# Print the result
print(paste("Top 1%: ", top1, "\nTop 5%:", top5, "\nTop 10%:", top10,sep=""))

# Extract values above the >44m threshold
prediction_band1[values(prediction_band1) <= top10 ] <- NA

# Count the number of valid (non-NA) pixels
top_values_count <- sum(!is.na(values(prediction_band1)))

# Extract values above the top 5% threshold
prediction_band1[values(prediction_band1) <= top5] <- NA

# Count the number of valid (non-NA) pixels
top5_values_count <- sum(!is.na(values(prediction_band1)))

# Extract values above the top 1% threshold
prediction_band1[values(prediction_band1) <= top1] <- NA

# Count the number of valid (non-NA) pixels
top1_values_count <- sum(!is.na(values(prediction_band1)))

# Print the result
print(paste("Big Tree top 10%:", top_values_count,"Big Tree top 5%:", top5_values_count,
            "Big Tree top 1%:", top1_values_count ))

# Verify the result (summary of top1_values)
prediction_band1 <- ch_msk
prediction_band1[values(prediction_band1) <= top10] <- NA

# Extract coordinates for the top 1% values as Positive for Big Tree
top_coords <- which(values(prediction_band1) > top10) 

# Retrieve the coordinates of the top 1% values using xyFromCell
top_points <- xyFromCell(prediction_band1, top_coords)

# Convert the coordinates to an sf object
top_sf <- st_as_sf(data.frame(x = top_points[, 1], y = top_points[, 2]), coords = c("x", "y"))

# Set the CRS of top1_sf to be the same as the Prediction raster's CRS
top_sf <- st_set_crs(top_sf, crs(Prediction))

nrow(top_sf) # Positive Value for Big Trees
#####################################

rm(top_points,top_coords,top_values,prediction_band1)

gc()
# ###########################################################################
# ##################################################################################
mytheme =  theme(legend.position="bottom",plot.margin = unit(c(0,30,0,0), "pt"),
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


# # # Use exactextract to get the raster values at the top 1% locations
names(top_sf)
nrow(top_sf)
top_sf$FID <- c(1:nrow(top_sf))

# # Extract Height values at these points
top_sf$Pred_hgt <- raster::extract(ch_msk, top_sf, fun = "mean")


top_sf$Prd_hgt = top_sf$Pred_hgt$Pred_Smooth_S12_10GLOBMAE_ensemble_1
gc()


#### Get BECzone value to tall canopy points
BEC = st_read("BEC_VI.shp")
BEC <- st_transform(BEC, crs(top_sf))

top_sf = st_join(top_sf,BEC[c("ZONE","SUBZONE","VARIANT","PHASE", "NATURAL_DI", "MAP_LABEL" , "BGC_LABEL",  "ZONE_NAME",  "SUBZONE_NA","VARIANT_NA", "PHASE_NAME")])
names(top_sf)

top_sf$PHASE[is.na(top_sf$PHASE)] = ""
top_sf$VARIANT[is.na(top_sf$VARIANT)] = ""
top_sf$BEC_name = top_sf$MAP_LABEL

# Cut Blocks
CBlocks = st_read("VI_CCutBlocks.shp")
CBlocks<- st_transform(CBlocks, crs(top_sf))
top_sf = st_join(top_sf,CBlocks[c("HARVEST__1", "PERCENT_CL")])
duplicate = duplicated(top_sf$geometry)
table(duplicate)
names(top_sf)
summary(top_sf$Pred_hgt)

#### Get Vegetation resource inventory (VRI) height and age values to tall canopy points
VRI = st_read("VRI_VI.shp")
VRI <- st_transform(VRI, crs(top_sf))
names(VRI)

top_sf = st_join(top_sf,VRI[c("PROJ_HEIGH","PROJ_HEI_1", "PROJ_AGE_1","SITE_INDEX")])
names(top_sf)
nrow(top_sf)


duplicate = duplicated(top_sf$geometry)
table(duplicate)
top_sf = top_sf[duplicate==FALSE,]
nrow(top_sf)


####### Update stand age with harvested years and classify forest into age classes
top_sf = top_sf %>%
  mutate(HARVEST = if_else(is.na(HARVEST__1), 0, HARVEST__1), PROJ_AG =PROJ_AGE_1) %>%
  mutate(StandA = case_when(
    HARVEST > 0 & HARVEST <=2019 ~ 2019 - HARVEST,
    HARVEST == 0 ~ PROJ_AG
  ))

top_sf = top_sf %>%
  mutate(
    Prd_hgt = Prd_hgt,
    Age_C = case_when(
      PROJ_AG >= 250 ~ "Old",
      PROJ_AG < 250 & PROJ_AG >=70 ~ "Mature",
      PROJ_AG < 70 ~ "Young",
      TRUE ~ NA_character_
    ),
    Hght_C = case_when(
      Prd_hgt >= 54 ~ "Top 1%",
      Prd_hgt < 54 & Prd_hgt >= 48 ~ "Top 5%", 
      Prd_hgt < 48 ~ "Top 10%",
      TRUE ~ NA_character_
    ),
    Age_Hgt = paste(Age_C, Hght_C, sep = "-"),
    Age_Hgt = factor(Age_Hgt,
                     levels = c("Young-Top 10%", "Young-Top 5%", "Young-Top 1%",
                                "Mature-Top 10%", "Mature-Top 5%", "Mature-Top 1%",
                                "Old-Top 10%", "Old-Top 5%", "Old-Top 1%"))
  )

top_sf = top_sf %>%
  mutate(
    Prd_hgt = Prd_hgt,
    Age_C = case_when(
      StandA >= 250 ~ "Old",
      StandA < 250 & StandA >=70 ~ "Mature",
      StandA < 70 ~ "Young",
      TRUE ~ NA_character_
    ),
    Hght_C = case_when(
      Prd_hgt >= 54 ~ "Top 1%",
      Prd_hgt < 54 & Prd_hgt >= 48 ~ "Top 5%", 
      Prd_hgt < 48 ~ "Top 10%",
      TRUE ~ NA_character_
    ),
    Age_HgtCB = paste(Age_C, Hght_C, sep = "-"),
    Age_HgtCB = factor(Age_HgtCB,
                     levels = c("Young-Top 10%", "Young-Top 5%", "Young-Top 1%",
                                "Mature-Top 10%", "Mature-Top 5%", "Mature-Top 1%",
                                "Old-Top 10%", "Old-Top 5%", "Old-Top 1%"))
  )
table(top_sf$Age_HgtCB)
table(top_sf$Age_Hgt)
gc()


# Extract DEM values at these points
DEM_sf <- raster::extract(DEM, vc, fun = "mean")
nrow(DEM_sf)
nrow(vc)
# Extract slope values at these points (first band of slope_Asp)
slope_sf <- raster::extract(slope, vc, fun = "mean")

# Extract wetness index values at these points
Wet_ind_sf <- raster::extract(Wet_ind, vc, fun = "mean")

# Extract Topographic Position Index values at these points
TPI_sf <- raster::extract(TPI, vc, fun = "mean")
table(is.na(TPI_sf$TPI_500))

# Extract Climate values at these points
clim_sf <- raster::extract(clim, vc, fun = "mean")

# Extract Human Footprint values
HFP_sf<-raster::extract(HFP, vc, fun = "mean")

# Extract Human Footprint values
DFO_sf<-raster::extract(DFO, vc, fun = "mean")

# Extract Human Footprint values
DWS_sf<-raster::extract(DWS, vc, fun = "mean")

names(clim_sf)[2:23] = c("AHM",	"bFFP","CMD","CMI","DD_0",	"DD_18","DD18", "DD5","eFFP","EMT",	"Eref","EXT","FFP",	"MAT",
                          "MCMT","MSP","MWMT","NFFD","PAS","RH","SHM","TD")


# # nrow(top1_Pred2024)
nrow(DEM_sf)
nrow(slope_sf)
nrow(Wet_ind_sf)
nrow(TPI_sf)
nrow(clim_sf)
nrow(DFO_sf)
nrow(DWS_sf)

#### Get PA_OECM value to tall canopy points
# top_sf = st_as_sf(vc)
PA_OECM = st_read("VI_PA.shp")
PA_OECM <- st_transform(PA_OECM, crs(Prediction))


names(PA_OECM)
top_sf = st_join(top_sf,PA_OECM[c("PA_OECM", "PA_TYPE","NAME_E","OBJECTID" )], left = TRUE)
top_sf$PA_OECM[is.na(top_sf$PA_OECM)] = "Not Protected"
top_sf$PA_TYPE[is.na(top_sf$PA_TYPE)] = "Not Protected"
top_sf$NAME_E[is.na(top_sf$NAME_E)] = "Not Protected"
top_sf$OBJECTID[is.na(top_sf$OBJECTID)] = "Not Protected"
# data = st_join(data,PA_OECM[c("Poly")], left = TRUE)
# Ensure only unique points are retained
data <- top_sf%>%
  distinct(geometry, .keep_all = TRUE)

names(data)
nrow(data)*100/10000
data = data[c(1:31,35:36)]
nrow(DEM_sf)
## Create dataset
data = cbind(data,DEM_sf)
data = cbind(data,slope_sf)
data = cbind(data,Dom_sps_sf)
data = cbind(data,Wet_ind_sf)
data = cbind(data,TPI_sf)
data = cbind(data,clim_sf)
data = cbind(data,DFO_sf)
data = cbind(data,DWS_sf)
data = cbind(data,HFP_sf)

###############################################################################
df = data%>%st_drop_geometry()

nrow(df)
names(df)
df = df%>%mutate(HARVEST = if_else(is.na(HARVEST),0,HARVEST))%>%
  mutate(Harv = if_else(HARVEST>=2019,-70,Prd_hgt))%>%filter(Prd_hgt>=44)
df$Harv

df$BEC = df$BEC_name
unique(df$BEC)

table(df$BEC)
nrow(df)
#Calculate the area per stange age group for the 3 definitions
df = df%>%mutate(Age_C = case_when(
                PROJ_AG >= 250 ~ "Old",
                PROJ_AG < 250 & PROJ_AG >=80 ~ "Mature",
                PROJ_AG < 80 ~ "Young",
                TRUE ~ NA_character_ ))
df = df%>%mutate(StandA = if_else(HARVEST>0 & HARVEST<2019, 2019-HARVEST,PROJ_AG))%>%
    mutate( StandAC = case_when(
      StandA >= 250 ~ "Old",
      StandA < 250 & StandA >=80 ~ "Mature",
      StandA < 80 ~ "Young",
      TRUE ~ NA_character_ ))


# Example: define thresholds
thresholds <- c("44m" = top10, "48m" =top5, "54m" = top1)

nrow(df[df$Prd_hgt>=top10,])*100/10000000
nrow(df[df$Prd_hgt>=top5,])*100/10000000
nrow(df[df$Prd_hgt>=top1,])*100/10000000


# Loop over thresholds and calculate area per Age_C

results <- lapply(names(thresholds), function(name) {
  thr <- thresholds[name]
  
  df %>%
    filter(Prd_hgt >= thr) %>%
    group_by(StandAC,BEC) %>%
    summarise(
      Count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Threshold = name,
      Area = Count*100/10000000  # convert pixels to ha if 1 pixel = 100 m2
    ) %>%
    group_by(BEC) %>%
    mutate(
      Proportion = Area / sum(Area)  # proportion of StandAC within each BEC
    ) %>%
    ungroup()
}) %>% bind_rows()

OMY_area = results%>%mutate(Area = Area)%>%
  group_by(StandAC,BEC) %>%#summarise(Area = sum(Area))%>%
  summarise(
    xmin = min(Area),          # lower bound
    xmid = median(Area),       # intermediate value for bar
    xmax = max(Area),          # upper bound
    .groups = "drop"
  )%>%data.frame()

# Set dodge width
dodge <- position_dodge(width = 0.8)

p1 = OMY_area%>%filter(!is.na(StandAC))%>% mutate(
  StandAC = factor(StandAC, levels = c("Young", "Mature", "Old")),
  BEC = factor(BEC, levels = c("CDFmm", "CWHxm1",  "CWHxm2","CWHmm1",
                               "CWHvh1", "CWHmm2", "MHmm1", "CWHvm2" ,"CWHvm1"))) %>%
  ggplot(aes(y = BEC, x = xmid, group=StandAC, fill=StandAC)) +
  geom_bar(stat = "identity", position = dodge, width = 0.7) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), position = dodge, height = 0.4) +
  geom_point(aes(x = xmid), position = dodge, size = 0.5, color = "black", show.legend = FALSE) +
  labs(x = "Area (1000 ha)", y = "Biogeoclimatic (BEC) Variant") +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
    ),
    labels = c(
      "Young","Mature","Old-growth"
    ),
    name = "Big Tree Forests"  # Remove legend title
  ) +
  mytheme +
  theme(legend.position=c(0.8,0.15), axis.text.x = element_text(hjust = 1))

OMY_area = results%>%filter(!is.na(StandAC))%>%
  group_by(StandAC, BEC) %>%
  summarise(
    xmin = min(Proportion),          # lower bound
    xmid = median(Proportion),       # intermediate value for bar
    xmax = max(Proportion),          # upper bound
    .groups = "drop"
  )%>%data.frame()

p2 = OMY_area%>% mutate(
  StandAC = factor(StandAC, levels = c("Young", "Mature", "Old")),
  BEC = factor(BEC, levels = c("CDFmm", "CWHxm1",  "CWHxm2","CWHmm1",
                               "CWHvh1", "CWHmm2", "MHmm1", "CWHvm2" ,"CWHvm1"))) %>%
  ggplot(aes(y = BEC, x = xmid, group=StandAC, fill=StandAC)) +
  geom_bar(stat = "identity", position = dodge, width = 0.7) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), position = dodge, height = 0.4) +
  geom_point(aes(x = xmid), position = dodge, size = 0.5, color = "black", show.legend = FALSE) +
  labs(x = "Class Proportion (%)", y = "Biogeoclimatic (BEC) Variant") +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
    ),
    labels = c(
      "Young","Mature","Old-growth"
    ),
    name = "Big Tree Forests"  # Remove legend title
  ) +
  mytheme +
  theme(legend.position="", axis.text.x = element_text(hjust = 1),
        legend.title=element_text( size = 12, color = "black", face = "bold"),
        axis.text.y = element_blank(), axis.title.y = element_blank())

combined_pt = wrap_plots(p1|p2) +
  plot_annotation(tag_levels = "A")

ggsave("Paper_Images/BEC_Age_BT.jpeg",combined_pt, height=6, width=12, units="in", dpi=300)

################# Figure S4 ##############
results2 <- lapply(names(thresholds), function(name) {
  thr <- thresholds[name]
  
  df %>%
    filter(Prd_hgt >= thr) %>%
    group_by(PA_OECM,StandAC, BEC) %>%
    summarise(
      Count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Threshold = name,
      Area = Count*100 / 10000  # convert pixels to ha if 1 pixel = 100 m2
    ) %>%
    group_by(PA_OECM, BEC) %>%
    mutate(
      Proportion = Area / sum(Area)  # proportion of StandAC within each BEC
    ) %>%
    ungroup()
}) %>% bind_rows()

results2%>%filter(  BEC=="CDFmm")

OMY_area = results2%>%mutate(Area = Area)%>%filter(!is.na(StandAC))%>%
  group_by(StandAC, BEC, Threshold) %>%summarise(Area=sum(Area)/1000)%>%
  group_by(StandAC, BEC) %>%
  summarise(
    xmin = min(Area),          # lower bound
    xmid = median(Area),       # intermediate value for bar
    xmax = max(Area),          # upper bound
    .groups = "drop"
  )%>%data.frame()

OMY_area%>%filter(BEC=="CWHxm1")
# Set dodge width
dodge <- position_dodge(width = 0.8)
unique(df$PA_OECM)
p3 = OMY_area%>% mutate(
  StandAC = fct_relevel(StandAC, c("Young","Mature","Old" )),
  Class = "All Big-tree Forests",
  # PA_OECM = factor(PA_OECM, levels = c("Not Protected", "Protected Area","OECM")),
  BEC = factor(BEC, levels = c("CDFmm", "CWHxm1",  "CWHxm2","CWHmm1",
                               "CWHvh1", "CWHmm2", "MHmm1", "CWHvm2" ,"CWHvm1"))) %>%
  ggplot(aes(y = BEC, x = xmid, group=StandAC, fill=StandAC)) +
  geom_bar(stat = "identity", position = dodge, width = 0.7) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), position = dodge, height = 0.4) +
  geom_point(aes(x = xmid), position = dodge, size = 0.5, color = "black", show.legend = FALSE) +
  facet_wrap(~Class,scales = "free_x")+
  labs(x = "Total Area (1000 ha)", y = "Biogeoclimatic (BEC) Variant") +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    # labels = c(
    #   "Young","Mature","Old-growth"
    #   # "Not Protected" , 
    #   # "OECM",
    #   # "Protected"
    #   
    # ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(legend.position=c(0.65,0.15), axis.text.x = element_text(hjust = 1))
p3

OMY_area = results2%>%filter(!is.na(StandAC))%>%
  group_by(PA_OECM, StandAC, BEC) %>%
  summarise(
    xmin = min(Proportion),          # lower bound
    xmid = median(Proportion),       # intermediate value for bar
    xmax = max(Proportion),          # upper bound
    .groups = "drop"
  )%>%data.frame()

OMY_area <- OMY_area %>%
  mutate(across(where(is.factor), as.character)) %>%
  bind_rows(data.frame(
    PA_OECM = "OECM",
    StandAC = "Old",
    BEC = "CDFmm",
    xmin = 0,
    xmid = 0,
    xmax = 0
  ))


pts = c()
z=0
for (i in unique(OMY_area$PA_OECM)){
  p4 = OMY_area%>% mutate(
    PA_OECM = factor(PA_OECM, levels = c("Not Protected","Protected Area","OECM" )),
    StandAC = factor(StandAC, levels = c("Young", "Mature", "Old")),
    BEC = factor(BEC, levels = c("CDFmm", "CWHxm1",  "CWHxm2","CWHmm1",
                                 "CWHvh1", "CWHmm2", "MHmm1", "CWHvm2" ,"CWHvm1"))) %>%
    filter(PA_OECM==i)%>%
    ggplot(aes(y = BEC, x = xmid, group=StandAC, fill=StandAC)) +
    geom_bar(stat = "identity", position = dodge, width = 0.7) +
    geom_errorbarh(aes(xmin = xmin, xmax = xmax), position = dodge, height = 0.4) +
    facet_wrap(~PA_OECM,scales = "free_x")+
    geom_point(aes(x = xmid), position = dodge, size = 0.5, color = "black", show.legend = FALSE) +
    labs(x = "Class Proportion (%)", y = "Biogeoclimatic (BEC) Variant") +
    scale_fill_manual(
      values = c(
        "Young"="#33CCFF",
        "Mature"="#99CC66",
        "Old" = "#006600"
        # "Not Protected" = "#A0A0A0",
        # "OECM" = "#88AADD",
        # "Protected" = "#1B7837"
      ),
      labels = c(
        "Young","Mature","Old-growth"
        # "Not Protected" , 
        # "OECM",
        # "Protected"
        
      ),
      name = "Stand Age"  # Remove legend title
    ) +
    mytheme +
    theme(legend.position="", axis.text.x = element_text(hjust = 1),
          legend.title=element_text( size = 12, color = "black", face = "bold"),
          axis.text.y = element_blank(), axis.title.y = element_blank())
  
  p4
  z=z+1
  pts[[z]]=p4
}
combined_pt = (p3 | pts[[1]] | pts[[3]] | pts[[2]]) +
  # plot_layout(guides = "collect", widths = c(1, 1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave("Paper_Images/BEC_PAs_BT.jpeg",combined_pt, height=6, width=14, units="in", dpi=300)

################################################################################
####### BEC values
names(df)
df %>%mutate(TPI=scale(TPI_250),DEM=Elevation, OECM=PA_OECM, HFP =HFP )%>%
  filter(!is.na(DEM),Prd_hgt>=top5) %>%
  mutate(
    Binned = case_when(
      HFP < 1   ~ "Wilderness",
      HFP >= 1  & HFP < 4  ~ "Low",
      HFP >= 4  & HFP < 10 ~ "Moderate",
      HFP >= 10 ~ "Highly Modified"
      # Binned = cut(
      #   HFP,
      #   breaks = c(seq(0, 10, by = 1), Inf),
      #   labels = c(paste(seq(0, 9, by = 1), seq(1,10, by = 1), sep = " to "), ">10"),
      #   include.lowest = TRUE,
      #   right = FALSE
    ),
    PA = case_when(
      OECM == "Not Protected" ~ "Not Protected",
      OECM == "Protected Area" ~ "Protected",
      OECM == "OECM" ~ "OECM"
    )
  )%>%
  group_by(BEC,OGC) %>%
  summarise(
    Count = n()*100/10000,
    ord = mean(HFP, na.rm = TRUE)#,.groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(BEC, levels(BEC)),  # preserve order
    Percentage = round((Count / sum(Count)) * 100,2)
  ) %>%
  filter(!is.na(Binned), BEC%in%c("CDFmm","CWHxm1","CWHxm2")) %>% data.frame()

################################################################################
## Figure S5

# Plot distribution for DEM with values > 1200 grouped
plot1 <- df %>%mutate(DEM=Elevation, OECM=PA_OECM)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      DEM,
      breaks = c(seq(0, 1200, by = 100), Inf),  # Define breaks for bins
      labels = c(paste(seq(0, 1200 - 100, by = 100), seq(100, 1200, by = 100), sep = "-"), paste0(">", 1200)),
      include.lowest = TRUE,
      right = FALSE  # Closed interval on the left
    )
  ) %>%
  # mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
  #                       OECM=="Protected Area" ~'Protected',
  #                       OECM=="OECM" ~'OECM')) %>% 
  group_by(Binned,OGC) %>%
  summarize(
    Count = n(),
    ord = mean(DEM),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),
    Percentage = (Count / sum(Count)) * 100# Ensure bins remain ordered correctly
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "DEM", x = "DEM (m)", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"

    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot1+theme(legend.position=c(0.8,0.9))

# Plot for slope using the original function (as no special binning is required)

plot2 <- df %>%mutate(Slope=slope, OECM=PA_OECM, DEM=Elevation)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      Slope,
      breaks = c(seq(0, 80, by = 8), Inf),  # Define breaks for bins
      labels = c(paste(seq(0, 80 - 8, by = 8), seq(8, 80, by = 8), sep = "-"), paste0(">", 80)),
      include.lowest = TRUE,
      right = FALSE  # Closed interval on the left
    )
  ) %>%
  # mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
  #                       OECM=="Protected Area" ~'Protected',
  #                       OECM=="OECM" ~'OECM')) %>%  
  group_by(Binned, OGC) %>%
  summarize(
    Count = n(),
    ord = mean(Slope),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),
    Percentage = (Count / sum(Count)) * 100  # Calculate percentages
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Slope", x = "Slope (degrees)", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot2+theme(legend.position=c(0.8,0.9))

# Helper function to bin data and calculate percentages, including conversion of aspect to direction
# Example for Aspect conversion to direction
plot3 <- df %>%mutate(Aspect=aspect, OECM=PA_OECM, DEM=Elevation)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Rounded = Aspect, 
    Direction = case_when(
      Rounded >= 0   & Rounded < 22.5   ~ "North",
      Rounded >= 22.5  & Rounded < 67.5  ~ "North-East",
      Rounded >= 67.5  & Rounded < 112.5 ~ "East",
      Rounded >= 112.5 & Rounded < 157.5 ~ "South-East",
      Rounded >= 157.5 & Rounded < 202.5 ~ "South",
      Rounded >= 202.5 & Rounded < 247.5 ~ "South-West",
      Rounded >= 247.5 & Rounded < 292.5 ~ "West",
      Rounded >= 292.5 & Rounded < 337.5 ~ "North-West",
      Rounded >= 337.5 & Rounded <= 360 ~ "North"
    )
  ) %>%
  # mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
  #                       OECM=="Protected Area" ~'Protected',
  #                       OECM=="OECM" ~'OECM')) %>%  
  group_by(Direction, OGC) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%  # Convert counts to percentages
  data.frame()%>%
  ggplot(aes(x = Direction, y = Percentage)) +  # Correct placement of aes()
  geom_bar(stat = "identity",  color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Aspect", x = "Aspect Direction", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot3+theme(legend.position=c(0.8,0.9))

# Plot distribution for TPI with 
summary(df$TPI_100)
table(df$OGC)
plot4 <- df %>%mutate(TPI=TPI_500/255,DEM=Elevation, OECM=PA_OECM)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      TPI,
      breaks = c(seq(-1, 1, by = 0.2), Inf),
      labels = c(paste(seq(-1, 1 - 0.2, by = 0.2), seq(-1 + 0.2, 1, by = 0.2), sep = " to "), ">1"),
      include.lowest = TRUE,
      right = FALSE
    ),
    PA = case_when(
      OECM == "Not Protected" ~ "Not Protected",
      OECM == "Protected Area" ~ "Protected",
      OECM == "OECM" ~ "OECM"
    )
  ) %>%
  group_by(Binned, OGC) %>%
  summarise(
    Count = n(),
    ord = mean(TPI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),  # preserve order
    Percentage = (Count / sum(Count)) * 100
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Topographic Posiition Index", x = "TPI", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot4+theme(legend.position=c(0.8,0.9))


# Human Footprint
summary(df$HFP)
summary(df$PA_OECM)
plot5 <- df %>%mutate(TPI=scale(TPI_250),DEM=Elevation, OECM=PA_OECM)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    # Binned = case_when(
      # HFP < 1   ~ "Wilderness",
      # HFP >= 1  & HFP < 4  ~ "Low",
      # HFP >= 4  & HFP < 10 ~ "Moderate",
      # HFP >= 10 ~ "Highly Modified"
    Binned = cut(
      HFP,
      breaks = c(seq(0, 10, by = 1), Inf),
      labels = c(paste(seq(0, 9, by = 1), seq(1,10, by = 1), sep = " to "), ">10"),
      include.lowest = TRUE,
      right = FALSE
    ),
    PA = case_when(
      OECM == "Not Protected" ~ "Not Protected",
      OECM == "Protected Area" ~ "Protected",
      OECM == "OECM" ~ "OECM"
    )
  ) %>%
  group_by(Binned, OGC) %>%
  summarise(
    Count = n(),
    ord = mean(TPI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),  # preserve order
    Percentage = (Count / sum(Count)) * 100
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Human Footprint Index", x = "HFP", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot5+theme(legend.position=c(0.8,0.9))

# Plot for slope using the original function (as no special binning is required)
plot6 <- df %>%mutate(TPI=scale(TPI_250),DEM=Elevation,Slope=slope, OECM=PA_OECM, DFOcean=DFOcean)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      DFOcean,
      breaks = c(seq(0, 35000, by = 5000), Inf),  # Define breaks for bins
      labels = c(paste(seq(0, 35000 - 5000, by = 5000), seq(5000, 35000, by = 5000), sep = "-"), paste0(">", 35000)),
      include.lowest = TRUE,
      right = FALSE  # Closed interval on the left
    )
  ) %>%
  mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
                        OECM=="Protected Area" ~'Protected',
                        OECM=="OECM" ~'OECM')  
  ) %>%
  group_by(Binned, OGC) %>%
  summarise(
    Count = n(),
    ord = mean(TPI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),  # preserve order
    Percentage = (Count / sum(Count)) * 100
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC),alpha=0.75) +
  labs(title = "Distance from the Ocean", x = "Distance (m)", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot6+theme(legend.position=c(0.8,0.9))


# Plot for slope using the original function (as no special binning is required)
plot7 <- df %>%mutate(TPI=scale(TPI_250),DEM=Elevation,Slope=slope, OECM=PA_OECM)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      MAP,
      breaks = c(seq(0, 2400, by = 200), Inf),  # Define breaks for bins
      labels = c(paste(seq(0, 2400 - 200, by = 200), seq(200, 2400, by = 200), sep = "-"), paste0(">", 2400)),
      include.lowest = TRUE,
      right = FALSE  # Closed interval on the left
    )
  ) %>%
  mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
                        OECM=="Protected Area" ~'Protected',
                        OECM=="OECM" ~'OECM')
  ) %>%
  group_by(Binned, OGC) %>%
  summarise(
    Count = n(),
    ord = mean(TPI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),  # preserve order
    Percentage = (Count / sum(Count)) * 100
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Mean Annual Preciptation (mm)", x = "Preciptation (mm)", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot7+theme(legend.position=c(0.8,0.9))

# Plot for slope using the original function (as no special binning is required)
plot8 <- df %>%mutate(TPI=scale(TPI_250),DEM=Elevation,Slope=slope, OECM=PA_OECM, DFOcean=DWS)%>%
  filter(!is.na(DEM), Prd_hgt>=top5, !is.na(OGC)) %>%#     99% = 54.34176
  mutate(
    Binned = cut(
      DFOcean,
      breaks = c(seq(0, 1000, by = 100), Inf),  # Define breaks for bins
      labels = c(paste(seq(0, 1000 - 100, by = 100), seq(100, 1000, by = 100), sep = "-"), paste0(">", 1000)),
      include.lowest = TRUE,
      right = FALSE  # Closed interval on the left
    )
  ) %>%
  mutate(PA = case_when(OECM=="Not Protected" ~'Not Protected',
                        OECM=="Protected Area" ~'Protected',
                        OECM=="OECM" ~'OECM')  
         ) %>%
  group_by(Binned, OGC) %>%
  summarise(
    Count = n(),
    ord = mean(TPI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Binned = fct_relevel(Binned, levels(Binned)),  # preserve order
    Percentage = (Count / sum(Count)) * 100
  ) %>%
  filter(!is.na(Binned)) %>% data.frame()%>%
  ggplot(aes(x = Binned, y = Percentage)) +
  geom_bar(stat = "identity", color='black',aes(fill = OGC), alpha=0.75) +
  labs(title = "Distance from Wetlands and Streams", x = "Distance (m)", y = "Percentage (%)") +
  coord_flip() +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_discrete(labels = function(x) gsub("\\[|\\]|\\(|\\)", "", x)) +
  scale_fill_manual(
    values = c(
      "Young"="#33CCFF",
      "Mature"="#99CC66",
      "Old" = "#006600"
      # "Not Protected" = "#A0A0A0",
      # "OECM" = "#88AADD",
      # "Protected" = "#1B7837"
    ),
    labels = c(
      "Young","Mature","Old-growth"
      # "Not Protected" , 
      # "OECM",
      # "Protected"
      
    ),
    name = "Stand Age"  # Remove legend title
  ) +
  mytheme +
  theme(axis.text.x = element_text(hjust = 1))

plot8+theme(legend.position=c(0.8,0.9))


combined_plot <- (
  plot1 + plot2 + plot3 + plot8 +
    plot5 +plot4+ plot7 + plot6 #plot_spacer()
) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "none")  # Hide legends from each plot

# Extract the legend from one plot
legend <- cowplot::get_legend(plot1 + theme(legend.position = "right", legend.key = element_rect(fill = NA),
                                            legend.background = element_rect(fill = "white"),
                                            legend.text = element_text(size = 16)))
legend_plot <- ggplot() + theme_void() + annotation_custom(legend)

# Insert the legend inside the first plot
final_plot <- combined_plot+
  plot_annotation(tag_levels = "A") + 
  inset_element(legend_plot, left = 3.9, bottom = 0.34, right = 0, top =0) & 
  theme(plot.tag = element_text(size = 14), plot.title = element_text(size = 12))  # You can change 14 to any desired size

ggsave("Paper_Images/DEM_Slope_BT_SAge.png", plot = final_plot, width = 14, height = 14, dpi = 300)
