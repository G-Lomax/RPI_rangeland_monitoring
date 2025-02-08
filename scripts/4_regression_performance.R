## This script reports and visualises the performance of the random forest
## regression models fitted in 3_quantile_regression.R

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 24th September 2024

## 1. Setup ----

source("scripts/load.R")

## 2. Load data and models ----

task_gpp <- read_rds("data/processed/rds/regr_task_gpp.rds")

rf_model_sp <- read_rds("data/processed/rds/rf_tuned_sp.rds")

# Country polygons for maps

ke_tz <- st_read("data/raw/vector/kenya_tanzania.geojson")

## 3. Calculate model performance metrics ----

# Measures to evaluate model performance
measures <- msrs(c("regr.mae", "regr.rmse", "regr.rsq"))

rf_predictions_sp <- rf_model_sp$predict(task_gpp)

rf_predictions_sp$score(measures)

## 4. Visualise model fit ----

# Extract quantile predictions

rf_predictions_sp_qu <- predict(
  rf_model_sp$learner$model,
  task_gpp$data(),
  type = "quantiles",
  quantiles = 0.9
)

data_with_quantiles <- task_gpp$data() %>%
  bind_cols(rf_predictions_sp_qu$predictions) %>%
  rename(quantile_pred = "quantile= 0.9") %>%
  mutate(across(.cols = c("GPP", starts_with("quantile")), .fns = ~ .x / 1000))

# Plot distribution of model residuals (using 0.9 quantile) 
density_breaks <- c(1, 10, 100, 1000, 10000)

density_plot_sp <- ggplot(data_with_quantiles, aes(x = quantile_pred, y = GPP)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(direction = 1, trans = "log", breaks = density_breaks) +
  geom_abline(slope = seq(0.2, 1, 0.2), intercept = 0, colour = "grey", lwd = 0.8, linetype = "longdash") +
  geom_abline(slope = 1, intercept = 0, colour = "grey", lwd = 1.6) +
  theme_classic() +
  theme(legend.position = c(0.18, 0.75), legend.key.height = unit(1.5, "cm"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent")) +
  xlim(0, 16) +
  ylim(0, 16) +
  labs(x = expression(atop("Potential GPP", (kg~C~m^-2~yr^-1))),
       y = expression(atop("Actual GPP", (kg~C~m^-2~yr^-1))),
       fill = "Number of points")

ggsave(
  "results/figures/density_plot_sp.png",
  density_plot_sp,
  width = 24, height = 24, units = "cm", dpi = 300
)

## 5. Variable importance ----

# Variable importance

importance_sp <- rf_model_sp$learner$importance()

static_vars <- c(
  "pptMean",
  "treeCover",
  "ECO_ID",
  "DEM",
  "slope",
  "merit_twi_fd8",
  "dist_to_river",
  "sand"
)

var_labels <- c(
  precipitation = "<i>Annual precipitation</i>",
  mean_ppt = "Mean annual precipitation",
  temperature_2m = "Mean annual temperature",
  pptIntensity = "Precipitation intensity",
  pptMeanDayAnomaly = "Anomaly in<br>mean precipitation day",
  GMT_0900_PAR = "Mean PAR",
  pptAnomaly = "Annual precipitation anomaly",
  potential_evaporation_sum = "Potential evapotranspiration",
  ugi = "Unranked Gini index",
  sand = "Soil sand fraction",
  treeCover = "Tree cover fraction",
  ECO_ID = "<i>Ecoregion</i>",
  DEM = "<i>Elevation</i>",
  slope = "Slope",
  distToRiver = "Distance to river",
  landform = "Landform",
  twi = "Topographic Wetness Index"
)

var_labels_df <- data.frame(variable = names(var_labels), label = unname(var_labels))

importance_df <- tibble(variable = names(importance_sp)) %>%
  mutate(importance_sp = importance_sp[variable],
         type = ifelse(variable %in% static_vars, "Static", "Annual")) %>%
  left_join(var_labels_df)

importance_plot <- importance_df %>%
  mutate(label = reorder(label, importance_sp)) %>%
  pivot_longer(starts_with("importance")) %>%
  ggplot(aes(x = value, y = label, colour = type)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("darkblue", "red")) +
  theme_bw() +
  labs(x = "Variable importance", y = "Variable", colour = "Variable type") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        axis.text.y = ggtext::element_markdown())

ggsave(
  "results/figures/var_importance.png",
  importance_plot,
  width = 20, height = 20, units = "cm", dpi = 300
)

## 6. Comparison with previous versions

# Convert rasters to consistent data frames covering the same spatial extent
rpi_rast <- rast("data/processed/raster/rpi_rast_sp.tif")
rpi_rast_legacy <- rast("data/processed/raster/rpi_legacy.tif")

rpi_rast_legacy_crop <- rpi_rast_legacy %>%
  crop(rpi_rast[[1]], extend = TRUE, mask = TRUE)

rpi_rast_crop <- rpi_rast[[1:nlyr(rpi_rast_legacy_crop)]] %>%
  mask(rpi_rast_legacy_crop)

rpi_df <- rpi_rast_crop %>%
  as.data.frame(cells = TRUE) %>%
  tidy_annual_vars() %>%
  drop_na()

rpi_legacy_df <- rpi_rast_legacy_crop %>%
  as.data.frame(cells = TRUE) %>%
  tidy_annual_vars() %>%
  drop_na()


# Calculate performance metrics for current and legacy model
calc_performance <- function(truth, predicted, na.rm = FALSE) {

  error <- predicted - truth

  mae <- mean(abs(error), na.rm = na.rm)
  rmse <- sqrt(mean(error ^ 2, na.rm = na.rm))
  rsq <- 1 - sum((error ^ 2), na.rm = na.rm) / sum((truth - mean(truth, na.rm = na.rm)) ^ 2, na.rm = na.rm)

  c(mae = mae, rmse = rmse, rsq = rsq)

}

# 
# perf_full <- calc_performance(rpi_df$GPP, rpi_df$gpp_predicted)
# perf_legacy <- calc_performance(rpi_legacy_df$GPP, rpi_legacy_df$gpp_predicted)

# Calculate pixel-wise performance in explaining temporal extent

gpp_actual <- subset(rpi_rast, str_detect(names(rpi_rast), "GPP"))
gpp_pred <- subset(rpi_rast, str_detect(names(rpi_rast), regex("^gpp_predicted")))

gpp_actual_subset <- subset(rpi_rast_crop, str_detect(names(rpi_rast_crop), "GPP"))
gpp_pred_subset <- subset(rpi_rast_crop, str_detect(names(rpi_rast_crop), regex("^gpp_predicted")))

gpp_actual_legacy <- subset(rpi_rast_legacy_crop, str_detect(names(rpi_rast_legacy_crop), "GPP"))
gpp_pred_legacy <- subset(rpi_rast_legacy_crop, str_detect(names(rpi_rast_legacy_crop), regex("^gpp_predicted")))

# (Load RESTREND residuals to compare)
restrend_old <- rast("data/processed/raster/restrend_resids_old.tif") %>%
  crop(rpi_rast[[1]], mask = TRUE, extend = TRUE)
restrend_old_resids <- subset(restrend_old, str_detect(names(restrend_old), "resid"))
restrend_old_preds <- gpp_actual_legacy - restrend_old_resids

# restrend_new <- rast("data/processed/raster/restrend.tif")
# restrend_new_resids <- subset(restrend_new, str_detect(names(restrend_new), "resid"))
# restrend_new_preds <- gpp_actual - restrend_new_resids

calc_performance_raster <- function(truth, predicted, na.rm = FALSE) {
  
  annual_anomaly <- truth - mean(truth, na.rm = na.rm)
  pred_anomaly <- predicted - mean(predicted, na.rm = na.rm)
  error <- pred_anomaly - annual_anomaly
  
  mae <- mean(abs(error), na.rm = na.rm)
  rmse <- sqrt(mean(error ^ 2, na.rm = na.rm))
  rsq <- 1 - (sum(error ^ 2, na.rm = na.rm) / sum(annual_anomaly ^ 2, na.rm = na.rm))
  
  output <- c(mae, rmse, rsq)
  names(output) <- c("mae", "rmse", "rsq")
  output
  
}

perf_new_full_temporal <- calc_performance_raster(gpp_actual, gpp_pred)
perf_new_subset_temporal <- calc_performance_raster(gpp_actual_subset, gpp_pred_subset)
perf_legacy_temporal <- calc_performance_raster(gpp_actual_legacy, gpp_pred_legacy)
perf_restrend_old_temporal <- calc_performance_raster(gpp_actual_legacy, restrend_old_preds)
# perf_restrend_new_temporal <- calc_performance_raster(gpp_actual, restrend_new_preds)

writeRaster(perf_new_full_temporal,
            "data/processed/raster/temporal_performance_full_new.tif",
            overwrite = TRUE)
writeRaster(perf_new_subset_temporal,
            "data/processed/raster/temporal_performance_new.tif",
            overwrite = TRUE)
writeRaster(perf_legacy_temporal,
            "data/processed/raster/temporal_performance_old.tif",
            overwrite = TRUE)
writeRaster(perf_restrend_old_temporal,
            "data/processed/raster/temporal_performance_restrend_old.tif",
            overwrite = TRUE)
# writeRaster(perf_restrend_new_temporal,
#             "data/processed/raster/temporal_performance_restrend_new.tif",
#             overwrite = TRUE)

perf_new_full_temporal <- rast("data/processed/raster/temporal_performance_full_new.tif")
perf_new_subset_temporal <- rast("data/processed/raster/temporal_performance_new.tif")
perf_legacy_temporal <- rast("data/processed/raster/temporal_performance_old.tif")
perf_restrend_old_temporal <- rast("data/processed/raster/temporal_performance_restrend_old.tif")
# perf_restrend_new_temporal <- rast("data/processed/raster/temporal_performance_restrend_new.tif")

# Plot performance metrics

performance_df_old <- as.data.frame(perf_new_subset_temporal, cells = TRUE) %>%
  rename_with(~ paste0(.x, ".new_subset"), -cell) %>%
  left_join(as.data.frame(perf_legacy_temporal, cells = TRUE) %>% rename_with(~paste0(.x, ".old"), -cell)) %>%
  pivot_longer(-cell) %>%
  separate_wider_delim(name, ".", names = c("measure", "version")) %>%
  drop_na()

performance_df_restrend <- as.data.frame(perf_new_full_temporal, cells = TRUE) %>%
  rename_with(~ paste0(.x, ".new_full"), -cell) %>%
  left_join(as.data.frame(perf_restrend_new_temporal, cells = TRUE) %>% rename_with(~paste0(.x, ".restrend"), -cell)) %>%
  pivot_longer(-cell) %>%
  separate_wider_delim(name, ".", names = c("measure", "version")) %>%
  drop_na()

measure_labels <- as_labeller(c(mae = "Mean~Absolute~Error~(g~C~m^-2~yr^-1)", rmse = "RMSE", rsq = "R^2"), label_parsed)

temporal_performance_hist <- performance_df_old %>%
  filter(measure != "rmse") %>%
  mutate(version = ordered(version, levels = c("old", "new_subset"))) %>%
  arrange(version) %>%
  ggplot(aes(x = value, fill = version)) +
  geom_density(alpha = 0.6, kernel = "rectangular") +
  facet_wrap(~measure, ncol = 1, nrow = 3, scales = "free", labeller = measure_labels) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", labels = c("RPI v1", "RPI v2")) +
  ggh4x::facetted_pos_scales(x = list(
    measure == "mae" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rmse" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rsq" ~ scale_x_continuous(limits = c(-0.5, 1))
  )) +
  geom_vline(xintercept = 0, colour = "grey10") +
  labs(x = "Value", y = "Density", fill = "Version")

temporal_performance_hist_restrend <- performance_df_restrend %>%
  filter(measure != "rmse") %>%
  mutate(version = ordered(version, levels = c("restrend", "new_full"))) %>%
  arrange(version) %>%
  ggplot(aes(x = value, fill = version)) +
  geom_density(alpha = 0.6, kernel = "rectangular") +
  facet_wrap(~measure, ncol = 1, nrow = 3, scales = "free", labeller = measure_labels) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", labels = c("RESTREND", "RPI v2")) +
  ggh4x::facetted_pos_scales(x = list(
    measure == "mae" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rmse" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rsq" ~ scale_x_continuous(limits = c(-0.5, 1))
  )) +
  geom_vline(xintercept = 0, colour = "grey10") +
  labs(x = "Value", y = "Density", fill = "Method")

ggsave("results/figures/rpi_version_performance_hist.png",
       temporal_performance_hist,
       width = 16, height = 16, units = "cm", dpi = 250)

ggsave("results/figures/rpi_restrend_performance_hist.png",
       temporal_performance_hist_restrend,
       width = 16, height = 16, units = "cm", dpi = 250)

# Calculate overall performance in explaining spatial patterns in mean GPP

gpp_actual_mean <- mean(gpp_actual)
gpp_pred_mean <- mean(gpp_pred)

gpp_actual_subset_mean <- mean(gpp_actual_subset)
gpp_pred_subset_mean <- mean(gpp_pred_subset)

gpp_actual_legacy_mean <- mean(gpp_actual_legacy)
gpp_pred_legacy_mean <- mean(gpp_pred_legacy)

spatial_performance_new <- calc_performance(values(gpp_actual_subset_mean), values(gpp_pred_subset_mean), na.rm = TRUE)
spatial_performance_old <- calc_performance(values(gpp_actual_legacy_mean), values(gpp_pred_legacy_mean), na.rm = TRUE)


# Maps of performance by pixel

fill <- tm_shape(ke_tz, is.main = FALSE) +
  tm_fill(fill = "grey95")
borders <- tm_shape(ke_tz) +
  tm_borders(col = "black")

rpi_old_mae_map <- fill +
  tm_shape(perf_legacy_temporal$mae, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0, 1500),
      values = ("oranges"),
      outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(
      title = "MAE (kg C m-2 yr-1)",
      reverse = TRUE,
      frame = FALSE
    )
  ) +
  borders +
  tm_layout(frame = FALSE)

rpi_new_mae_map <- fill +
  tm_shape(perf_new_subset_temporal$mae, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0, 1500),
      values = ("oranges"),
      outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(
      title = "MAE (kg C m-2 yr-1)",
      reverse = TRUE,
      frame = FALSE
    )
  ) +
  borders +
  tm_layout(frame = FALSE)

rpi_best_mae <- which.min(c(perf_legacy_temporal$mae, perf_new_subset_temporal$mae))
rpi_best_rsq <- which.max(c(perf_legacy_temporal$rsq, perf_new_subset_temporal$rsq))

generate_performance_map <- function(performance_rast, metric, palette) {
  # Prep data
  layer <- performance_rast[[metric]]
  
  if (metric == "rsq") {
    negative_rsq <- layer < 0
    negative_rsq <- terra::mask(negative_rsq, negative_rsq, maskvalues = c(0, NA))
  }
  
  label_lookup <- c(
    mae = expression(atop("MAE", (kg~C~m^-2~yr^-1))),
    rmse = expression("RMSE ("~kg^2~m^-4~yr^-2~")"),
    rsq = expression(R^2)
  )
  
  limits_lookup <- list(
    mae = c(0, 1000),
    rmse = c(0, 1000),
    rsq = c(0, 1)
  )
  
  label <- label_lookup[metric]
  limits <- limits_lookup[[metric]]

  # Prep map
  
  fill <- tm_shape(ke_tz, is.main = FALSE) +
    tm_fill(fill = "grey95")
  borders <- tm_shape(ke_tz) +
    tm_borders(col = "black")
  
  map <- fill +
    tm_shape(layer, is.main = TRUE) +
    tm_raster(
      col.scale = tm_scale_continuous(
        limits = limits,
        values = palette,
        outliers.trunc = c(TRUE, TRUE),
        midpoint = NA
      ),
      col.legend = tm_legend(
        title = label,
        reverse = TRUE,
        frame = FALSE
      )
    ) +
    borders +
    tm_layout(frame = FALSE, asp = 3/4)
  
  if (metric == "rsq" & min(values(layer), na.rm = TRUE) < 0) {
    map <- map +
      tm_shape(negative_rsq) +
      tm_raster(col.scale = tm_scale_categorical(values = c("lightsalmon"), labels = expression(Negative~R^2)),
                col.legend = tm_legend(title = ""))
  }
  
  map
}

rpi_old_mae_map <- generate_performance_map(perf_legacy_temporal, "mae", "or_rd")
rpi_new_mae_map <- generate_performance_map(perf_new_subset_temporal, "mae", "or_rd")
rpi_old_rsq_map <- generate_performance_map(perf_legacy_temporal, "rsq", "viridis")
rpi_new_rsq_map <- generate_performance_map(perf_new_subset_temporal, "rsq", "viridis")

# Plot differences in performance between two methods

rpi_mae_diff_map <- fill +
  tm_shape(perf_legacy_temporal$mae - perf_new_subset_temporal$mae, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(-300, 300), values = "pu_or", midpoint = 0, outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(
      reverse = TRUE, title = expression(atop("MAE difference", (kg~C~m^-2~yr^-1))), frame = FALSE
    )
  ) +
  borders +
  tm_layout(frame = FALSE, asp = 3/4)
  

rpi_rsq_diff_map <- fill +
  tm_shape(perf_new_subset_temporal$rsq - perf_legacy_temporal$rsq, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(-0.75, 0.75), values = "pu_or", midpoint = 0, outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(
      reverse = TRUE, title = expression(R^2~difference), frame = FALSE
    )
  ) +
  borders +
  tm_layout(frame = FALSE, asp = 3/4)

# Save
tmap_save(rpi_old_mae_map, "results/figures/rpi_old_mae_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)
tmap_save(rpi_new_mae_map, "results/figures/rpi_new_mae_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)
tmap_save(rpi_old_rsq_map, "results/figures/rpi_old_rsq_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)
tmap_save(rpi_new_rsq_map, "results/figures/rpi_new_rsq_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)
tmap_save(rpi_mae_diff_map, "results/figures/rpi_mae_diff_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)
tmap_save(rpi_rsq_diff_map, "results/figures/rpi_rsq_diff_map.png",
          width = 12, height = 12, units = "cm", dpi = 300)


## Plot overlapping study area

legacy_sa <- !is.na(rpi_rast_legacy$potential_gpp_predicted.2000)
new_sa <- !is.na(rpi_rast$potential_gpp_predicted.2000)

new_sa_ext <- extend(new_sa, legacy_sa)
legacy_sa_ext <- extend(legacy_sa, new_sa_ext)
overlap <- legacy_sa_ext + 2 * new_sa_ext
overlap_masked <- mask(overlap, overlap, maskvalues = c(0, NA))

levels(overlap_masked) <- data.frame(ID = 1:3, category = c("RPI v1", "RPI v2", "Both versions"))

shared_study_area_map <- fill +
  tm_shape(overlap_masked, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "set1",
      # labels = c("RPI v1", "RPI v2", "Both versions")
    ),
  col.legend = tm_legend(title = "Study area", frame = FALSE)
  ) +
  borders

