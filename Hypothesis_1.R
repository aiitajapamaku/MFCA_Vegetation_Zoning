
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
-
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  install.packages("ggh4x")
}
library(ggh4x)

# --- 1. Set working directory ---
setwd("/Users/aiitajoshuaapamaku/Desktop/MFNP/20250708")

# --- 2. Load points (zones) & rasters ---
pts_externalbuf <- read_sf("Buffer_zone_points.shp")
pts_core        <- read_sf("Core_points.shp")
pts_intbuf      <- read_sf("ptos.shp")

mfca2024 <- rast("DynamicWorld_Consensus_Uganda_MFCA2024.tif")
evi2024  <- rast("EVI_50thPercentile_2024_MFCA.tif")
evi2018  <- rast("EVI_50thPercentile_2018_MFCA.tif")
cv.evi   <- rast("EVI_CV_2324_MFCA_Harmonized_CloudMasked.tif")

# --- 3. Reproject points to match raster CRS ---
pts_externalbuf <- st_transform(pts_externalbuf, crs(mfca2024))
pts_core        <- st_transform(pts_core,        crs(mfca2024))
pts_intbuf      <- st_transform(pts_intbuf,      crs(mfca2024))

# Add Zone column
pts_externalbuf$Zone <- "ExternalBuffer"
pts_core$Zone        <- "Core"
pts_intbuf$Zone      <- "InternalBuffer"

# Combine all zones
pts <- bind_rows(pts_externalbuf, pts_core, pts_intbuf)
pts$ID <- 1:nrow(pts)

# --- 4. Preprocess Land Cover raster ---
mfca2024 <- round(mfca2024)
mfca2024 <- ifel(mfca2024 %in% c(1,2,3,5), mfca2024, NA)

# --- 5. Extract raster values for each point ---
lc_vals <- terra::extract(mfca2024, pts, df=TRUE)
names(lc_vals)[2] <- "LC"

evi24 <- terra::extract(evi2024, pts, df=TRUE)
names(evi24)[2] <- "EVI2024"

evi18 <- terra::extract(evi2018, pts, df=TRUE)
names(evi18)[2] <- "EVI2018"

cv_vals <- terra::extract(cv.evi, pts, df=TRUE)
names(cv_vals)[2] <- "EVI_CV"

# Combine all data
df_all <- pts %>%
  dplyr::select(ID, Zone) %>%
  left_join(lc_vals,  by="ID") %>%
  left_join(evi24,    by="ID") %>%
  left_join(evi18,    by="ID") %>%
  left_join(cv_vals,  by="ID")

df_all <- df_all %>% filter(!is.na(LC))

# Recode LC to factor with Dynamic World labels
df_all$LC <- factor(df_all$LC,
                    levels=c(1,2,3,5),
                    labels=c("Forest","Shrub & Scrub","Grass","Flooded Vegetation"))

# Compute ΔEVI
df_all <- df_all %>% mutate(DeltaEVI = EVI2024 - EVI2018)

# --- 6. Summary statistics ---
summary_stats <- df_all %>%
  group_by(Zone, LC) %>%
  summarise(
    mean_EVI2024   = mean(EVI2024, na.rm=TRUE),
    median_EVI2024 = median(EVI2024, na.rm=TRUE),
    sd_EVI2024     = sd(EVI2024, na.rm=TRUE),
    mean_EVI_CV    = mean(EVI_CV, na.rm=TRUE),
    mean_DeltaEVI  = mean(DeltaEVI, na.rm=TRUE),
    .groups="drop"
  )
print(summary_stats)

# --- 7. Visualization settings ---
dw_colors <- c("Forest"="#1C5F2C",
               "Shrub & Scrub"="#DFC35A",
               "Grass"="#88B053",
               "Flooded Vegetation"="#7A87C6")

# Facet strip colors for zones
zone_strip_colors <- list(
  ExternalBuffer="#FF9999",
  InternalBuffer="#99FF99",
  Core="#FFFF99"
)

# Panel background colors — increased opacity (~80%)
zone_panel_colors <- c(
  "ExternalBuffer"="#FF9999CC",         
  "InternalBuffer"="#99FF99CC",         
  "Core"="#FFFF99CC"                     
)

y_limits <- c(0, 0.75)

# --- 8. Plot function using ggh4x and stronger panel overlays ---
plot_zone_overlay <- function(yvar, ylab, title) {
  ggplot(df_all, aes(x=LC, y=.data[[yvar]], fill=LC)) +
    # Panel background overlay
    geom_rect(data=distinct(df_all, Zone),
              aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=Zone),
              inherit.aes=FALSE) +
    geom_boxplot(position="dodge") +
    ggh4x::facet_wrap2(~Zone,
                       strip = ggh4x::strip_themed(
                         background_x = ggh4x::elem_list_rect(zone_strip_colors)
                       )) +
    scale_fill_manual(values=c(dw_colors, zone_panel_colors)) +
    coord_cartesian(ylim=y_limits) +
    theme_minimal() +
    labs(title=title, y=ylab, x="Vegetation Class", fill="Class")
}

# --- 9. Generate plots with updated titles ---
plot_zone_overlay("EVI2024", "EVI", "Vegetation condition")
plot_zone_overlay("EVI_CV", "Coefficient of Variation (EVI)", "Coefficient of variation")
plot_zone_overlay("DeltaEVI", "ΔEVI", "Change in Vegetation condition")

# Remove geometry column (if summary_stats is an sf object)
summary_df <- sf::st_set_geometry(summary_stats, NULL)

# Save as CSV
write.csv(summary_df, "MFCA_summary_stats.csv", row.names = FALSE)

