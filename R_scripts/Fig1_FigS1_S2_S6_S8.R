rm(list = ls()); gc()
library(pacman) 

p_load(raster,terra,exactextractr,sf,ggplot2,dplyr,patchwork,tidyverse,viridis, scales)

setwd('G:/Big_tree/')

# --- Load your rasters ---
ref <- rast("G:/Rasters/Testing_LiDAR.tif")      # ground truth raster
pred <- rast("G:/Rasters/CCvoer.tif")     # predicted raster

# Make sure they align
ref<-project(ref, pred)
pred <- resample(pred, ref)


# Extract coordinates for the top 1% values as Positive for Big Tree
top_coords <- which(values(ref) >=0) 

# Retrieve the coordinates of the top 1% values using xyFromCell
top_points <- xyFromCell(ref, top_coords)

# Convert the coordinates to an sf object
top_sf <- st_as_sf(data.frame(x = top_points[, 1], y = top_points[, 2]), coords = c("x", "y"))

# Set the CRS of top1_sf to be the same as the Prediction raster's CRS
top_sf <- st_set_crs(top_sf, crs(pred))

nrow(top_sf)

# --- Randomly sample 1% of total pixels ---
set.seed(42)
n_bins <- 100

n_sample <- 10407844   # number of points to randomly select (1% of all pixels)

top_sf_sample <- top_sf %>% slice_sample(n = n_sample)

top_sf_sample<-st_buffer(top_sf_sample, 10)

ref_df = exact_extract(ref, top_sf_sample, fun = "mean")
pred_df = exact_extract(pred, top_sf_sample, fun = "mean")

##

df <- data.frame(
  Reference = ref_df,
  Predicted = pred_df
)

df <- na.omit(df)


# --- Compute R2 ---
lm_model <- lm(Reference ~ Predicted, data = df) #
r2_value <- summary(lm_model)$r.squared  # Extract R²

## Figure S2 ######

p <- df%>%#filter(CC_ab5>=0)%>%
  ggplot( aes(x = Reference, y = Predicted)) +
  geom_bin2d(bins = 200) +  # bins resolution
  scale_fill_viridis_c(option = "D", 
                       limits = c(0, 500),  # max value for color scaling
                       oob = scales::squish,  
                       trans = "sqrt"  
  )+
  geom_abline(slope = 1, intercept =  0, linetype = "dashed", color = "red", size = 1) +
  # geom_smooth(method = "lm",  color = "grey", size = 0.5, formula = 'y ~ x') +#fill = "grey",
  theme_minimal() +
  coord_cartesian(ylim = c(0, 100)) +  
  labs(x = "Reference Canopy Cover (%)", y = "Predicted Canopy Cover (%)", fill = "Density")

# Add R� as an annotation
p + annotate("text", x = 1/15, y = 98, 
             label = paste0("R2 = ", round(r2_value, 3)), 
             size = 6, hjust = 0, fontface = "bold")+
  theme(legend.position=c(0.9,0.22),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

ggsave(paste("Paper_Images/CanopyCover.jpeg",sep=""), height=5, width=6, units="in", dpi=300)



################################################################################
########### Get height value from our predicted canopy height, Patopov, Lang, and Pauls' ##
# Predicted (Our)
pred_path = "Pred_Smooth_S12_10GLOBMAE_ensemble.tif"
prd<-terra::rast(pred_path)

## LiDAR reference
ldr_path = "p99height.tif"
ldr<-terra::rast(ldr_path)

## global height products
prods_path = "Height_Products.tif"
prod<-terra::rast(prods_path)

ref_df = exact_extract(ldr[[9]], top_sf_sample, fun = c("mean","median","max"))
pred_df = exact_extract(prd[[1]], top_sf_sample, fun = c("mean","median","max"))
hprod_df= exact_extract(prod, top_sf_sample, fun = c("mean","median","max"))

names(hprod_df) = c(
                    "mean_1m_RES", "mean_Sothe", "mean_Potapov", "mean_Lang", "mean_Pauls",
                    "median_1m_RES", "median_Sothe", "median_Potapov", "median_Lang", "median_Pauls",  
                    "max_1m_RES", "max_Sothe", "max_Potapov", "max_Lang", "max_Pauls"
                    )

### Get Biogeclimatic Variants (BEC)
BEC = st_read("data/VI_BEC.shp")

# Subset BEC columns
BEC_sub <- BEC[, c("ZONE","SUBZONE","VARIANT","PHASE", "MAP_LABEL", "BGC_LABEL", "ZONE_NAME", "PHASE_NAME")]

sf_sam = st_centroid(top_sf_sample)

df_sam = st_join(sf_sam,BEC_sub )

## Data set from comparision (1% of total reference pixels)
df = data.frame(
  ref_mean = ref_df$mean,
  ref_median = ref_df$median,
  ref_max = ref_df$max,
  pred_mean = pred_df$mean,
  pred_median = pred_df$median,
  pred_max = pred_df$max,
  BEC = df_sam$MAP_LABEL
  
)
df = cbind(df,hprod_df)
# Check new column names
names(df)
nrow(df)


## Figure 2 ####
lm_model <- lm(ref_median ~ median_Potapov, data = df) #
r2_value <- summary(lm_model)$r.squared  # Extract R²
r2_value

p <- df%>%#filter(MEDIAN_prdhgt>45&MEDIAN_ldrhgt>24)%>%
  ggplot( aes(x = ref_median, y = pred_median)) +
  geom_bin2d(bins = 90) +  # Adjust bins for resolution
  scale_fill_viridis_c(option = "A", 
                       limits = c(0, 1000),  # Set max value for color scaling
                       oob = scales::squish,  # Prevents values > 40000 from distorting the scale
                       trans = "sqrt"  # Adjusts contrast for better visibility) +  # Use viridis color scale
  )+
  geom_abline(slope = 1, intercept =  0, linetype = "dashed", color = "red", size = 1) +
  # geom_smooth(method = "lm",  color = "grey", size = 0.5, formula = 'y ~ x') +#fill = "grey",
  theme_minimal() +
  coord_cartesian(ylim = c(0, 70), xlim = c(0, 75)) +  
  labs(x = "Reference Canopy Height (m)", y = "Predicted Canopy Height (m)", fill = "Density")
# Add R² as an annotation
p + annotate("text", x = 1/15, y = 65, 
             label = paste0("R2 = ", round(r2_value, 3)), 
             size = 6, hjust = 0, fontface = "bold")+
  theme(legend.position=c(0.9,0.22),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", linewidth=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", linewidth=0.5, linetype="solid"),
        panel.grid.major=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(linewidth=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=12), axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 

ggsave(paste("Paper_Images/pred_ref.jpeg",sep=""), height=5, width=6, units="in", dpi=300)

## Figure S6
# Calculate all metrics for each BEC_VARIAN group
metrics <- df %>%filter(BEC!='CMAunp')%>%
  mutate(MEDIAN_prdhgt =pred_median ,MEDIAN_ldrhgt=ref_median)%>%
  st_drop_geometry() %>%
  group_by(BEC) %>%
  summarise(
    r2 = ifelse(all(is.na(MEDIAN_prdhgt) | is.na(MEDIAN_ldrhgt)), 
                NA_real_, 
                summary(lm(MEDIAN_prdhgt ~ MEDIAN_ldrhgt))$r.squared),
    mae = mean(abs(MEDIAN_prdhgt - MEDIAN_ldrhgt), na.rm = TRUE),
    rmse = sqrt(mean((MEDIAN_prdhgt - MEDIAN_ldrhgt)^2, na.rm = TRUE)),
    range = max(MEDIAN_ldrhgt, na.rm = TRUE)-min(MEDIAN_ldrhgt, na.rm = TRUE),
    nrmse = sqrt(mean((MEDIAN_prdhgt - MEDIAN_ldrhgt)^2, na.rm = TRUE)) /range *100,
    bias = mean(MEDIAN_prdhgt - MEDIAN_ldrhgt, na.rm = TRUE)  # Calculate bias
  ) %>%
  mutate(across(c(mae, rmse, bias), ~ ifelse(is.nan(.), NA_real_, .))) %>%
  mutate(label = sprintf("R2 = %.2f\nMAE = %.2f\nRMSE = %.2f\nBias = %.2f", 
                         r2, mae, rmse, bias))

# Your original plot code with annotation
p <- df %>%filter(BEC!='CMAunp')%>%
  ggplot(aes(y = pred_median, x = ref_median)) +
  geom_bin2d(bins = 90) +
  scale_fill_viridis_c(option = "A", 
                       limits = c(0, 500),
                       oob = scales::squish,
                       trans = "sqrt") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~BEC, ncol=3) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 65), xlim = c(0, 75)) +  
  labs(x = "Reference Canopy Height (m)", y = "VRI Canopy Height (m)", fill = "Density")

# Add annotations to each facet
p + geom_text(data = metrics, 
              aes(x = Inf, y = -Inf, label = label),
              hjust = 1.1, vjust = -0.5, size = 3.5) +
  theme(
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.text.y = element_text(size = 12, color = "black", face = "bold.italic"),
    legend.direction = 'vertical', legend.title = element_text(size = 12),
    legend.background = element_rect(fill = "white", linetype = "solid", size = 0.5, colour = "black"),
    legend.spacing.y = unit(1.5, "mm"), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.spacing.x = unit(0.1, "lines"), panel.spacing.y = unit(0.1, "lines"),
    panel.background = element_rect(fill = "white", colour = "grey", size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'dashed', colour = "lightgrey"),
    panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', colour = "lightgrey"), 
    plot.title = element_text(color = "black", size = 14, hjust = 0, face = 'bold'),
    axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), legend.key = element_rect(fill = "transparent", colour = "transparent"))

ggsave(paste("Paper_images/heightBEC.jpeg",sep=""), height=10, width=9, units="in", dpi=300)



######################################################################################
###################### plotting the error by groups of 5m height #################
## Figure S8
df$hght_res = df$pred_median-df$ref_median

table(df$BEC)
p_ref <- df %>%
  filter(!is.na(ref_median), BEC!="CMAunp") %>%
  ggplot(aes(x = BEC, y = ref_median)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill='darkgreen', alpha=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
  # coord_cartesian(ylim = c(-25, 35)) +
  coord_flip(ylim = c(0, 70))+
  scale_fill_viridis_d(option = "B", begin = 0.1) +  # Skip the very dark color
  theme_minimal() +
  labs(
    x = "BEC Zone",
    y = "LiDAR Reference Height (m)",
    title = "Boxplot of Reference Height by BEC Zone",
    fill = "BEC Zone"
  ) + 
  theme(legend.position="",
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_blank(),#element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),
        axis.ticks = element_line(color="black"),
        axis.title.y=element_blank(), axis.text.x=element_text(color="black",size=12),
        axis.text.y=element_text(color="black",size=12), legend.key=element_rect(fill="transparent", colour="transparent")) 


p_res2 <- df %>%
  filter(!is.na(ref_median), BEC!="CMAunp") %>%
  ggplot(aes(x = BEC, y = hght_res)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill='yellow', alpha=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
  # coord_cartesian(ylim = c(-25, 35)) +
  coord_flip(ylim = c(-12, 25))+
  
  theme_minimal() +
  labs(
    x = "",
    y = "Height Residual (m)",
    title = "Boxplot of Height Residuals by BEC Zone",
    fill = "BEC Zone"
  ) + 
  theme(legend.position="",
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=12),
        legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_blank(),#element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),
        axis.ticks = element_line(color="black"),
        axis.title.y=element_blank(), axis.text.x=element_text(color="black",size=12),
        axis.text.y=element_blank(), legend.key=element_rect(fill="transparent", colour="transparent")) 

print(p_res2)

# Calculate percent per group
percent_data <- df%>% data.frame()%>%
  filter(!is.na(hght_res), BEC!="CMAunp") %>%
  count(BEC) %>%
  mutate(percent = (n / sum(n)) * 100)

# Plot
BEC_dat<- data.frame(
  BEC = c("MHmm1", "CWHvm1", "CWHvm2", "CWHvh1", "CWHxm1", "CWHxm2", "CDFmm", "CWHmm2",  "CWHmm1"),
  total_area = c(10.05, 29.35, 11.39, 10.00, 7.48, 13.84, 5.28, 6.85,  3.97)
)
p_percent <- ggplot(percent_data, aes(x = BEC, y = percent)) +
  geom_bar(data =BEC_dat,color='black', aes( x = BEC, y = total_area, fill = "Total Area"), 
           stat = "identity", alpha = 1) +  # Overlay reference bars in red
  geom_bar(stat = "identity", aes( fill = "Sampled Area"),  color='black', alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "BEC Zone",
    y = "Sample Percentage (%)",
    title = "Percentage of Samples per BEC Zone"
  )  + 
  scale_fill_manual(values = c("Sampled Area" = "darkgrey", 
                               "Total Area" = "red"), 
                    name = "") +guides(fill = guide_legend(ncol = 1))+ 
  theme(legend.position=c(0.78,0.88),
        strip.text.x = element_text( size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(size = 12, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_blank(),
        legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="grey", size=0.5, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_blank(),#element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),
        axis.ticks = element_line(color="black"),
        axis.title.y=element_blank(), axis.text.x=element_text(color="black",size=12),
        axis.text.y=element_blank(), legend.key=element_rect(fill="transparent", colour="transparent")) 


print(p_percent)

(p_ref | p_res2 | p_percent) +
  # plot_layout(widths = c(2, 2, 1)) +
  plot_annotation(tag_levels = "A", tag_suffix = ")")


ggsave(paste("Paper_images/BEC_residuals.jpeg",sep=""), height=10, width=9, units="in", dpi=300)
#################### Height intervals #########################
### Figure S1
df$hght_res_ours = df$pred_median-df$ref_median
df$hght_res_1m = df$median_1m_RES-df$ref_median
df$hght_res_Sothe = df$median_Sothe-df$ref_median
df$hght_res_Potapov = df$median_Potapov-df$ref_median
df$hght_res_Lang = df$median_Lang-df$ref_median
df$hght_res_Pauls = df$median_Pauls/100-df$ref_median
# Pivot data to long format
p_res_long <- df %>%mutate(MAX_ldrhgt=ref_median)%>%
  filter(!is.na(hght_res_ours)) %>%
  mutate(Hght_clss = case_when(
    MAX_ldrhgt >= 0 & MAX_ldrhgt < 5  ~ "0-5m",
    MAX_ldrhgt >= 5 & MAX_ldrhgt < 10 ~ "5-10m",
    MAX_ldrhgt >= 10 & MAX_ldrhgt < 15 ~ "10-15m",
    MAX_ldrhgt >= 15 & MAX_ldrhgt < 20 ~ "15-20m",
    MAX_ldrhgt >= 20 & MAX_ldrhgt < 25 ~ "20-25m",
    MAX_ldrhgt >= 25 & MAX_ldrhgt < 30 ~ "25-30m",
    MAX_ldrhgt >= 30 & MAX_ldrhgt < 35 ~ "30-35m",
    MAX_ldrhgt >= 35 & MAX_ldrhgt < 40 ~ "35-40m",
    MAX_ldrhgt >= 40 & MAX_ldrhgt < 45 ~ "40-45m",
    MAX_ldrhgt >= 45 & MAX_ldrhgt < 50 ~ "45-50m",
    # MAX_ldrhgt >= 50 & MAX_ldrhgt < 55 ~ "50-55m",
    # MAX_ldrhgt >= 55 & MAX_ldrhgt < 60 ~ "55-60m",
    MAX_ldrhgt >= 50 ~ "50+m"
  )) %>%
  mutate(Hght_clss = factor(Hght_clss, levels = c("0-5m", "5-10m", "10-15m", 
                                                  "15-20m", "20-25m", "25-30m",
                                                  "30-35m", "35-40m", "40-45m", "45-50m",
                                                  "50+m"))) %>%#,"55+m", "55-60m","60+m"
  pivot_longer(cols = starts_with("hght_res_"), 
               names_to = "Method", 
               values_to = "Height_Residual") %>%st_drop_geometry()%>%data.frame()%>%
  mutate(Model = case_when(
    # Method == "hght_res_VRI" ~ "VRI",
    Method == "hght_res_ours" ~ "Our",
    # Method == "hght_res_ourNLL" ~ "Our+Wloss+NLL",
    # Method == "hght_res_ourUM" ~ "Our",
    Method == "hght_res_1m" ~ "1m_Res",
    Method == "hght_res_Sothe" ~ "Sothe",
    Method == "hght_res_Potapov" ~ "Potapov",
    Method == "hght_res_Lang" ~ "Lang",
    Method == "hght_res_Pauls" ~ "Pauls",
    TRUE ~ NA_character_  # Handles any unexpected cases
  ))%>%
  mutate(Model = factor(Model, levels = c("Our","1m_Res", "Sothe","Potapov","Lang","Pauls")))


# Plot with different colors for each method
# Generate Viridis colors and skip the first (black) one
num_colors <- length(unique(p_res_long$Model))  # Number of models

viridis_colors <- viridis(num_colors + 1, option = "B")[-c(1,3,4)]  # Skip the first (black) color

p_res <- p_res_long%>%
  filter(Model%in%c("Our","Potapov","Lang","Pauls")&!is.na(Hght_clss))%>%
  ggplot(aes(x = Hght_clss, y = Height_Residual, fill = Model)) +
  geom_hline(yintercept = 0, color = "grey",  size = 0.75) +  # Add grey line at y=0 linetype = "dashed",
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(0.8)) +  # Add whiskers explicitly
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8), size=0.25) +
  # scale_fill_viridis_d(option = "B") +  # Use Viridis color palette for categorical data
  coord_cartesian(ylim = c(-35, 35)) +  
  scale_fill_manual(values = viridis_colors,
                    labels = c(
                      "Our           (R²:0.75, MAE:4.90, RMSE:6.92, BIAS:0.34)",
                      # "VRI             (R²:0.16, MAE:11.2, RMSE:15.68, BIAS:-5.27)",
                      # "Ours           (R²:0.76, MAE:4.73, RMSE:6.79, BIAS:-0.86)",
                      "Potapov    (R²:0.39, MAE:8.09, RMSE:10.85, BIAS:-4.39)", 
                      "Lang          (R²:0.44, MAE:8.11, RMSE:10.20, BIAS:3.52)", 
                      "Pauls         (R²:0.56, MAE:7.06, RMSE:9.27, BIAS:-3.35)"
                    ),
                    name = ""  # Remove legend title
  ) +  # Manually apply Viridis without black
  theme_minimal() +
  labs(x = "", y = "Residuals (m)", title = "Comparison of Height Residuals") 
# +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
size = 14
p_res  + 
  # annotate("text", x = 1/15, y = 0.7, 
  #                 label = paste0("R? = ", round(r2_value, 3)), 
  #                 size = 6, hjust = 0, fontface = "bold")+
  theme(legend.position=c(0.73,0.9),
        strip.text.x = element_text( size = size, color = "black", face = "bold"),
        strip.text.y = element_text(size = size, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=size),legend.text = element_text(size=12),
        # legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        # legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="black", size=1, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=size), axis.text.x=element_text(size=size,color="black"),
        axis.text.y=element_text(color="black",size=size), legend.key=element_rect(fill="transparent", colour="transparent")) 

ggsave(paste("Paper_Images/prehghtvshght_residual_long.jpeg",sep=""), height=7, width=14, units="in", dpi=300)



#################################################################################
# BEC error (Not included in the paper)
p_res <- p_res_long%>%
  filter(Model%in%c("Our","Potapov","Lang","Pauls")&!is.na(Hght_clss), !is.na(BEC), BEC!="CMAunp")%>%
  ggplot(aes(x = BEC, y = Height_Residual, fill = Model)) +
  geom_hline(yintercept = 0, color = "grey",  size = 0.75) +  # Add grey line at y=0 linetype = "dashed",
  stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(0.8)) +  # Add whiskers explicitly
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.8), size=0.25) +
  # scale_fill_viridis_d(option = "B") +  # Use Viridis color palette for categorical data
  coord_flip()+
  coord_cartesian(xlim = c(-35, 35)) +  
  scale_fill_manual(values = viridis_colors,
                    labels = c(
                      "Our (Wloss+M)",
                      # "VRI",
                      "Potapov", 
                      "Lang", 
                      "Pauls"
                    ),
                    name = ""  # Remove legend title
  ) +  # Manually apply Viridis without black
  theme_minimal() +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "", y = "Residuals (m)", title = "Comparison of Height Residuals") 
# +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
size = 14
p_res  + 
  # annotate("text", x = 1/15, y = 0.7, 
  #                 label = paste0("R² = ", round(r2_value, 3)), 
  #                 size = 6, hjust = 0, fontface = "bold")+
  theme(legend.position=c(0.3,0.95),
        strip.text.x = element_text( size = size, color = "black", face = "bold"),
        strip.text.y = element_text(size = size, color = "black", face = "bold.italic" ),
        legend.direction='vertical', legend.title=element_text(size=size),legend.text = element_text(size=12),
        # legend.background = element_rect(fill="white", linetype="solid", size=0.5, colour="black"),
        # legend.spacing.y=unit(1.5, "mm"), panel.border=element_rect(colour="black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.1, "lines"),panel.spacing.y = unit(0.1, "lines"),
        panel.background=element_rect(fill="white", colour="black", size=1, linetype="solid"),
        panel.grid.major=element_line(size=0.25, linetype='dashed', colour="lightgrey"),
        panel.grid.minor=element_line(size=0.25, linetype='dashed', colour="lightgrey"), 
        plot.title=element_text(color="black", size=14, hjust=0,face='bold'), #axis.title.x=element_blank(),axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),#element_text(size=10, angle=45, hjust = 1), 
        axis.title.y=element_text(size=size), axis.text.x=element_text(size=size,color="black"),
        axis.text.y=element_text(color="black",size=size), legend.key=element_rect(fill="transparent", colour="transparent")) 


ggsave(paste("Paper_Images/prehghtvshght_residual_BEC.jpeg",sep=""), height=7, width=9, units="in", dpi=300)
