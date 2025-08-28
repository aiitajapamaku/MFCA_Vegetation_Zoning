library(terra)
library(spatialRF)
library(vegan)
library(car)
library(sf)

setwd("/Users/aiitajoshuaapamaku/Desktop/MFNP/20250708")

# 1. Load files 

ptos <- read_sf("ptos.shp")

mfca2024 <- rast("DynamicWorld_Consensus_Uganda_MFCA2024.tif") |>
  project(crs(ptos)) |>
  round()
mfca2024 <- ifel(mfca2024 %in% c(1,2,3,5), mfca2024, NA)

evi2024 <- rast("EVI_50thPercentile_2024_MFCA.tif")
evi2018 <- rast("EVI_50thPercentile_2018_MFCA.tif")
cv.evi  <- rast("EVI_CV_2324_MFCA_Harmonized_CloudMasked.tif")


# 2. Mask to natural vegetation

mask_nat <- ifel(mfca2024 %in% c(1,2,3,5), 1, NA) |>
  project(crs(evi2024)) |>
  resample(evi2024, method="near")

evi2024 <- evi2024 * mask_nat
evi2018 <- evi2018 * mask_nat


# 3. Classification layers (terciles)

qt_evi <- quantile(values(evi2024), c(0.33, 0.66), na.rm=TRUE)
cl.evi24 <- ifel(evi2024 <= qt_evi[1], 1,
                 ifel(evi2024 <= qt_evi[2], 2, 3))

qt_cv <- quantile(values(cv.evi), c(0.33, 0.66), na.rm=TRUE)
cl.cv24 <- ifel(cv.evi <= qt_cv[1], 1,
                ifel(cv.evi <= qt_cv[2], 2, 3))

cl.jt <- (cl.evi24 * 10) + cl.cv24
writeRaster(cl.jt, "clssEVIjt_20250708.tif", overwrite=TRUE)


# 4. Focal percentage layers

foc_layers <- lapply(c(11,12,13,21,22,23,31,32,33), function(i) {
  focal(ifel(cl.jt==i,1,0), matrix(1,5,5), "mean", na.rm=TRUE)
})
names(foc_layers) <- paste0("fc", c(11,12,13,21,22,23,31,32,33))


# 5. Extract values at points

dt.dw24  <- extract(mfca2024, ptos)[,2]
dt.evi24 <- extract(evi2024, ptos)[,2]
dt.evi18 <- extract(evi2018, ptos)[,2]

df_pred <- sapply(foc_layers, function(x) extract(x, ptos)[,2])
df_pred <- decostand(df_pred, "standardize")


# 6. Remove correlated predictors

dt_model_df <- spatialRF::auto_cor(x=as.data.frame(df_pred), cor.threshold=0.7)$selected.variables.df


# 7. Prepare modelling dataframe

valid_idx <- as.numeric(rownames(dt_model_df))
dt_model <- data.frame(
  clss = factor(dt.dw24[valid_idx], levels=c(1,2,3,5)),
  y    = dt.evi24[valid_idx] - dt.evi18[valid_idx],
  dt_model_df
)
dt_model <- dt_model[complete.cases(dt_model),]


# 8. Fit linear model (main effects, clean version)

library(car)
library(lmtest)
library(sandwich)


all_predictors <- setdiff(names(dt_model), c("y", "clss"))
formula_full <- as.formula(paste("y ~ clss +", paste(all_predictors, collapse=" + ")))


m_main <- lm(formula_full, data = dt_model)

# Check for aliased coefficients
aliased <- alias(m_main)$Complete
if(length(aliased) > 0){
  cat("Aliased terms detected, removing:\n")
  aliased_terms <- names(which(rowSums(aliased != 0) > 0))
  print(aliased_terms)
  # Remove aliased predictors
  predictors_clean <- setdiff(all_predictors, aliased_terms)
  formula_clean <- as.formula(paste("y ~ clss +", paste(predictors_clean, collapse=" + ")))
  m_main <- lm(formula_clean, data=dt_model)
}

# Model summary
summary(m_main)

# VIF check
vif_vals <- vif(m_main)
cat("VIF values for clean model:\n")
print(vif_vals)



# 9. Diagnostics with robust standard errors and ready-to-use summary table

library(lmtest)   # for coeftest and bptest
library(sandwich) # for robust SE
library(car)      # for residual plots and VIF

# 9.1 Robust covariance matrix (HC3)
robust_cov <- vcovHC(m_main, type = "HC3")

# 9.2 Coefficients with robust SEs
robust_se <- coeftest(m_main, vcov = robust_cov)
cat("Model coefficients with robust standard errors:\n")
print(robust_se)

# 9.3 Robust confidence intervals
robust_ci <- confint(m_main, vcov. = robust_cov)
cat("Robust 95% confidence intervals:\n")
print(robust_ci)

# 9.4 Build summary table with robust SEs and VIFs
coef_mat <- as.matrix(robust_se)

summary_table <- data.frame(
  Predictor = rownames(coef_mat),
  Estimate  = coef_mat[,1],
  Std.Error = coef_mat[,2],
  t.value   = coef_mat[,3],
  Pr.t      = coef_mat[,4],
  stringsAsFactors = FALSE
)

# Add VIF column initialized as NA
summary_table$VIF <- NA

# Assign VIFs to matching predictors only
for (pred in names(vif_vals[, "GVIF^(1/(2*Df))"])) {
  if (pred %in% summary_table$Predictor) {
    summary_table$VIF[summary_table$Predictor == pred] <- vif_vals[pred, "GVIF^(1/(2*Df))"]
  }
}

# Round numeric columns for reporting
summary_table[,2:6] <- round(summary_table[,2:6], 4)

cat("Summary table of coefficients, robust SEs, and VIFs (ready for reporting):\n")
print(summary_table)

# 9.5 Standard diagnostic plots
par(mfrow = c(2, 2))
plot(m_main)  # Residuals vs fitted, Q-Q, Scale-Location, Cook's distance

# 9.6 Histogram of residuals
hist(residuals(m_main), breaks = 50, col = "lightblue",
     main = "Histogram of Residuals", xlab = "Residuals")

# 9.7 Shapiro-Wilk test for normality
res_sample <- sample(residuals(m_main), min(3000, length(residuals(m_main))))
shapiro_test <- shapiro.test(res_sample)
cat("Shapiro-Wilk normality test:\n")
print(shapiro_test)

# 9.8 Breusch-Pagan test for heteroscedasticity
bp_test <- bptest(m_main)
cat("Breusch-Pagan test for heteroscedasticity:\n")
print(bp_test)

# 9.9 Cook's distance to detect influential points
cook_d <- cooks.distance(m_main)
plot(cook_d, type = "h", main = "Cook's Distance", ylab = "Cook's distance")
abline(h = 4 / nrow(dt_model), col = "red", lty = 2)  # threshold line

# Reset plotting layout
par(mfrow = c(1,1))


# 10. Coefficient plot with robust SEs and 95% CI

library(ggplot2)

# 10.1 Prepare data for plotting
plot_df <- summary_table
plot_df <- plot_df[plot_df$Predictor != "(Intercept)", ]  # Exclude intercept

# Calculate 95% CI using robust SE if not already included
plot_df$CI_lower <- plot_df$Estimate - 1.96 * plot_df$Std.Error
plot_df$CI_upper <- plot_df$Estimate + 1.96 * plot_df$Std.Error

# 10.2 Optional: order predictors by effect size for readability
plot_df$Predictor <- factor(plot_df$Predictor, levels = plot_df$Predictor[order(plot_df$Estimate)])

# 10.3 Create coefficient plot
ggplot(plot_df, aes(x = Estimate, y = Predictor)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Model Coefficients with Robust 95% Confidence Intervals",
    x = "Coefficient Estimate",
    y = "Predictor"
  ) +
  theme_minimal(base_size = 14)


# 11. Prediction of ΔEVI for each vegetation class


# 11.1 Extract mean values of numeric predictors
predictor_means <- colMeans(dt_model[, setdiff(names(dt_model), c("y", "clss"))], na.rm = TRUE)

# 11.2 Build new data.frame for prediction (one row per class)
pred_df <- data.frame(
  clss = factor(c(1, 2, 3, 5), levels = c(1,2,3,5))
)

# Repeat mean values for all predictors
for (pred in names(predictor_means)) {
  pred_df[[pred]] <- predictor_means[pred]
}

# 11.3 Predict ΔEVI using the fitted model
# Use robust SEs for confidence intervals
library(sandwich)
library(lmtest)

robust_cov <- vcovHC(m_main, type = "HC3")  # reuse robust covariance
pred_values <- predict(m_main, newdata = pred_df, se.fit = TRUE)

# Compute robust 95% CI for predictions
pred_df$DeltaEVI <- pred_values$fit
pred_df$CI_lower <- pred_values$fit - 1.96 * pred_values$se.fit
pred_df$CI_upper <- pred_values$fit + 1.96 * pred_values$se.fit

# 11.4 View results
print(pred_df)

# 11.5 Optional: plot predictions with error bars
library(ggplot2)
ggplot(pred_df, aes(x = clss, y = DeltaEVI)) +
  geom_point(size = 3, color = "forestgreen") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, color = "darkgreen") +
  labs(title = "Predicted ΔEVI by Vegetation Class",
       x = "Vegetation Class",
       y = "Predicted ΔEVI") +
  theme_minimal(base_size = 14)


# 12. Visualization of predicted ΔEVI with Vegetation Overlay (Composite Map)

library(ggplot2)
library(terra)
library(tmap)
library(viridis)

# --- 12.0 Dynamic World class info ---
dw_classes <- data.frame(
  clss  = 0:8,
  name  = c("water","trees","grass","flooded_vegetation","crops",
            "shrub_and_scrub","built","bare","snow_and_ice"),
  color = c("#419bdf","#397d49","#88b053","#7a87c6","#e49635",
            "#dfc35a","#c4281b","#a59b8f","#b39fe1") # actual Dynamic World colors
)

# --- 12.1 Merge prediction results with class info ---
dw_pred <- merge(pred_df, dw_classes, by = "clss")

# --- 12.2 Horizontal bar plot ---
ggplot(dw_pred, aes(x = DeltaEVI, y = name, fill = name)) +
  geom_col() +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.2, color = "black") +
  scale_fill_manual(values = setNames(dw_pred$color, dw_pred$name)) +
  labs(title = "Predicted ΔEVI by Vegetation Class",
       x = "Predicted ΔEVI", y = "Vegetation Class") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# --- 12.3 Raster of predicted ΔEVI ---
mfca2024_num <- mfca2024  # numeric DW codes 0–8
lut <- data.frame(from = as.numeric(dw_pred$clss), to = dw_pred$DeltaEVI)
pred_values_raster <- classify(mfca2024_num, lut, others = NA)

# --- 12.4 Save rasters ---
writeRaster(pred_values_raster, "predicted_DeltaEVI.tif", overwrite = TRUE)
writeRaster(mfca2024_num, "vegetation_classes.tif", overwrite = TRUE)

# --- 12.5 Prepare vegetation raster with names ---
used_classes <- sort(unique(values(mfca2024_num)))
used_classes <- used_classes[!is.na(used_classes)]

dw_colors_named <- setNames(
  dw_classes$color[dw_classes$clss %in% used_classes],
  dw_classes$name[dw_classes$clss %in% used_classes]
)

mfca2024_named <- mfca2024_num
values(mfca2024_named) <- factor(values(mfca2024_num),
                                 levels = dw_classes$clss[dw_classes$clss %in% used_classes],
                                 labels = dw_classes$name[dw_classes$clss %in% used_classes])

# --- 12.6 Composite tmap (ΔEVI + Vegetation) ---
tmap_mode("plot")
tmap_options(component.autoscale = FALSE)

tm_composite <- tm_shape(pred_values_raster) +
  tm_raster(
    col.scale = tm_scale_continuous(values = viridis(50, option = "A")), # green viridis ramp
    legend.show = FALSE
  ) +
  tm_shape(mfca2024_named) +
  tm_raster(
    col.scale = tm_scale_categorical(values = dw_colors_named), # Dynamic World colors
    col_alpha = 0.35,
    legend.show = FALSE
  ) +
  tm_layout(
    frame = FALSE,
    main.title = "Composite Map: ΔEVI and Vegetation Classes"
  )

tm_composite

# --- 12.7 Save composite map raster ---
composite_raster <- c(pred_values_raster, mfca2024_num)
writeRaster(composite_raster, "composite_map.tif", overwrite = TRUE)

# --- 12.8 Standalone legends ---
# ΔEVI legend
tm_legend_evi <- tm_shape(pred_values_raster) +
  tm_raster(
    col.scale = tm_scale_continuous(values = viridis(50, option = "A")),
    col.legend = tm_legend("Predicted ΔEVI")
  )
tm_legend_evi

# Vegetation class legend
tm_legend_veg <- tm_shape(mfca2024_named) +
  tm_raster(
    col.scale = tm_scale_categorical(values = dw_colors_named),
    col.legend = tm_legend("Vegetation Class")
  )
tm_legend_veg

#Saving my tables for reporting
write.csv(summary_table, file = "summary_table_hypothesis2.csv", row.names = FALSE)
write.csv(summary_table, file = "/Users/aiitajoshuaapamaku/Desktop/MFNP/summary_table_hypothesis2.csv", row.names = FALSE)

# Convert robust_ci to a data frame
robust_ci_df <- as.data.frame(robust_ci)
robust_ci_df$Predictor <- rownames(robust_ci_df)  # add predictor names as a column
rownames(robust_ci_df) <- NULL  # remove row names

# Reorder columns for clarity
robust_ci_df <- robust_ci_df[, c("Predictor", "2.5 %", "97.5 %")]

# Save as CSV
write.csv(robust_ci_df, file = "robust_CI.csv", row.names = FALSE)

# Assuming pred_df contains your predictions and confidence intervals
# Save as CSV for Excel
write.csv(pred_df, file = "predicted_DeltaEVI.csv", row.names = FALSE)
write.csv(pred_df, file = "/Users/aiitajoshuaapamaku/Desktop/MFNP/predicted_DeltaEVI.csv", row.names = FALSE)


