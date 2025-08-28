
# 1. Install & Load Packages

required_pkgs <- c("terra", "ggplot2", "dplyr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(terra)
library(ggplot2)
library(dplyr)


# 2. Set Working Directory & Load Rasters

setwd("/Users/aiitajoshuaapamaku/Desktop/MFNP/20250708")

# Load rasters
cl.jt_filled <- rast("clssEVIjt_filled.tif")     
pred_delta   <- rast("predicted_DeltaEVI.tif")   


# 3.CRS Alignment

target_crs <- crs(cl.jt_filled)  # reference CRS
rasters <- list(cl.jt_filled=cl.jt_filled, pred_delta=pred_delta)

for(i in names(rasters)){
  if(crs(rasters[[i]]) != target_crs){
    rasters[[i]] <- project(rasters[[i]], target_crs)
    cat(i, "reprojected.\n")
  }
}

cl.jt_filled <- rasters$cl.jt_filled
pred_delta   <- rasters$pred_delta


# 4. Define Productivity × Stability Color Ramp

class_colors <- c(
  "11" = "#D9D9D9",  # Low Productivity, High Stability
  "12" = "#D1BCE3",  # Low Productivity, Moderate Stability
  "13" = "#C39BCB",  # Low Productivity, Low Stability
  "21" = "#A3D79C",  # Moderate Productivity, High Stability
  "22" = "#89B4B4",  # Moderate Productivity, Moderate Stability
  "23" = "#9B83B3",  # Moderate Productivity, Low Stability
  "31" = "#49B845",  # High Productivity, High Stability
  "32" = "#43967C",  # High Productivity, Moderate Stability
  "33" = "#833F91"   # High Productivity, Low Stability
)


# 5. Save Separate Legend

png("Legend_ProductivityStability.png", width=600, height=800)
par(mar=c(1,1,1,1))
plot.new()
legend("center", legend=names(class_colors),
       fill=class_colors, title="Productivity × Stability Classes",
       cex=1.2, bty="n")
dev.off()
cat("Legend saved: Legend_ProductivityStability.png\n")


# 6. Composite Map (Classes + Predictor Overlay)

# Interactive view
plot(cl.jt_filled, col=class_colors, legend=FALSE,
     main="Composite: Productivity × Stability + Delta EVI",
     maxpixels=1e6)
plot(pred_delta, col=terrain.colors(100),
     alpha=0.6, legend=TRUE, add=TRUE)

# Save as PNG
png("Map_CompositePredictor.png", width=1200, height=1000)
plot(cl.jt_filled, col=class_colors, legend=FALSE,
     main="Composite: Productivity × Stability + Delta EVI",
     maxpixels=1e6)
plot(pred_delta, col=terrain.colors(100),
     alpha=0.6, legend=TRUE, add=TRUE)
dev.off()
cat("Composite map saved: Map_CompositePredictor.png\n")


# 7. Save Raster Layers

writeRaster(cl.jt_filled, "cljt_filled_raster.tif", overwrite=TRUE)
writeRaster(pred_delta, "predicted_delta_raster.tif", overwrite=TRUE)


# 8. Productivity × Stability Matrix (3x3 grid)

# Axes
stability_levels <- c("High Stability", "Moderate Stability", "Low Stability") # X-axis
productivity_levels <- c("Low Productivity", "Moderate Productivity", "High Productivity") # Y-axis

grid_df <- expand.grid(
  Stability = stability_levels,
  Productivity = productivity_levels
)

# Assign class codes
grid_df$class <- c("11","12","13",
                   "21","22","23",
                   "31","32","33")

# Map colors
grid_df$color <- class_colors[grid_df$class]

# Plot with ggplot
p_matrix <- ggplot(grid_df, aes(x=Stability, y=Productivity, fill=color)) +
  geom_tile(color="white", linewidth=1.2) +
  scale_fill_identity() +
  geom_text(aes(label=class), color="black", size=6) +
  scale_y_discrete(limits=rev(productivity_levels)) + # Low -> High top to bottom
  theme_minimal(base_size=14) +
  labs(title="Productivity × Stability Classes",
       x="Stability", y="Productivity") +
  theme(
    panel.grid=element_blank(),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14, face="bold"),
    plot.title=element_text(hjust=0.5, size=16, face="bold")
  )

# Print matrix
print(p_matrix)

# Save matrix as PNG
ggsave("Matrix_Productivity_Stability.png", p_matrix,
       width=8, height=6, dpi=300)
cat("Matrix grid saved: Matrix_Productivity_Stability.png\n")

