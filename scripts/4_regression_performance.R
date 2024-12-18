## This script reports and visualises the performance of the random forest
## regression models fitted in 3_quantile_regression.R

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 24th September 2024

## 1. Setup ----

source("scripts/load.R")

## 2. Load data and models ----

task_gpp <- read_rds("data/processed/rds/regr_task_gpp.rds")

rf_model_spt <- read_rds("data/processed/rds/rf_tuned_spt.rds")
rf_model_sp <- read_rds("data/processed/rds/rf_tuned_sp.rds")

## 3. Calculate model performance metrics ----

# Measures to evaluate model performance
measures <- msrs(c("regr.mae", "regr.rmse", "regr.rsq"))

rf_predictions_spt <- rf_model_spt$predict(task_gpp)
rf_predictions_sp <- rf_model_sp$predict(task_gpp)

rf_predictions_spt$score(measures)
rf_predictions_sp$score(measures)

## 4. Visualise model fit ----

# Extract quantile predictions
rf_predictions_spt_qu <- predict(
  rf_model_spt$learner$model,
  task_gpp$data(),
  type = "quantiles",
  quantiles = 0.9
)

rf_predictions_sp_qu <- predict(
  rf_model_sp$learner$model,
  task_gpp$data(),
  type = "quantiles",
  quantiles = 0.9
)

data_with_quantiles <- task_gpp$data() %>%
  bind_cols(rf_predictions_spt_qu$predictions, rf_predictions_sp_qu$predictions) %>%
  rename(quantile_pred_spt = "quantile= 0.9...18", quantile_pred_sp = "quantile= 0.9...19") %>%
  mutate(across(.cols = c("GPP", starts_with("quantile")), .fns = ~ .x / 1000))

# Plot distribution of model residuals (using 0.9 quantile) 
density_breaks <- c(1, 10, 100, 1000, 10000)

density_plot_spt <- ggplot(data_with_quantiles, aes(x = quantile_pred_spt, y = GPP)) +
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

density_plot_sp <- ggplot(data_with_quantiles, aes(x = quantile_pred_sp, y = GPP)) +
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
  "results/figures/density_plot_spt.png",
  density_plot_spt,
  width = 24, height = 24, units = "cm", dpi = 300
)

ggsave(
  "results/figures/density_plot_sp.png",
  density_plot_sp,
  width = 24, height = 24, units = "cm", dpi = 300
)

## 5. Variable importance ----

# Variable importance

importance_spt <- rf_model_spt$learner$importance()
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

importance_df <- tibble(variable = unique(c(names(importance_spt), names(importance_sp)))) %>%
  mutate(importance_spt = importance_spt[variable],
         importance_sp = importance_sp[variable],
         type = ifelse(variable %in% static_vars, "Static", "Annual"))

importance_plot <- importance_df %>%
  replace_na(list("importance_sp" = 0, "importance_spt" = 0)) %>%
  mutate(variable = reorder(variable, importance_sp)) %>%
  pivot_longer(starts_with("importance")) %>%
  ggplot(aes(x = value, y = variable, colour = name, shape = type)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_colour_manual(labels = c("Spatial block CV", "Spatiotemporal LLTO-CV"), values = c("darkblue", "red")) +
  labs(x = "Variable importance", y = "Variable", colour = "CV method", shape = "Variable type") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

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
  terra::resample(rpi_rast[[1]]) %>%
  crop(rpi_rast[[1]], mask = TRUE)

rpi_rast_crop <- rpi_rast[[1:nlyr(rpi_rast_legacy_crop)]] %>%
  crop(rpi_rast_legacy_crop[[1]], mask = TRUE)

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

perf_full <- calc_performance(rpi_df$GPP, rpi_df$gpp_predicted)
perf_legacy <- calc_performance(rpi_legacy_df$GPP, rpi_legacy_df$gpp_predicted)

# Calculate pixel-wise performance in explaining temporal extent

gpp_actual <- subset(rpi_rast_crop, str_detect(names(rpi_rast_crop), "GPP"))
gpp_pred <- subset(rpi_rast_crop, str_detect(names(rpi_rast_crop), regex("^gpp_predicted")))

gpp_actual_legacy <- subset(rpi_rast_legacy_crop, str_detect(names(rpi_rast_legacy_crop), "GPP"))
gpp_pred_legacy <- subset(rpi_rast_legacy_crop, str_detect(names(rpi_rast_legacy_crop), regex("^gpp_predicted")))

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

perf_full_temporal <- calc_performance_raster(gpp_actual, gpp_pred)
perf_legacy_temporal <- calc_performance_raster(gpp_actual_legacy, gpp_pred_legacy)

writeRaster(perf_full_temporal, "data/processed/raster/temporal_performance_new.tif")
writeRaster(perf_legacy_temporal, "data/processed/raster/temporal_performance_old.tif")

# Plot performance metrics

performance_df <- as.data.frame(perf_full_temporal, cells = TRUE) %>%
  rename_with(~ paste0(.x, "_new"), -cell) %>%
  left_join(as.data.frame(perf_legacy_temporal, cells = TRUE) %>% rename_with(~paste0(.x, "_old"), -cell)) %>%
  pivot_longer(-cell) %>%
  separate_wider_delim(name, "_", names = c("measure", "version"))

measure_labels <- as_labeller(c(mae = "MAE", rmse = "RMSE", rsq = "R-squared"))
temporal_performance_hist <- ggplot(performance_df, aes(x = value, fill = version)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~measure, ncol = 1, nrow = 3, scales = "free", labeller = measure_labels) +
  theme_classic() +
  ggh4x::facetted_pos_scales(x = list(
    measure == "mae" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rmse" ~ scale_x_continuous(limits = c(0, 1500)),
    measure == "rsq" ~ scale_x_continuous(limits = c(-1, 1))
  )) +
  geom_vline(xintercept = 0, colour = "grey10") +
  labs(x = "Value", y = "Density", fill = "Version")

# Calculate overall performance in explaining spatial patterns in mean GPP

gpp_actual_mean <- mean(gpp_actual)
gpp_pred_mean <- mean(gpp_pred)

gpp_actual_legacy_mean <- mean(gpp_actual_legacy)
gpp_pred_legacy_mean <- mean(gpp_pred_legacy)

spatial_performance_new <- calc_performance(values(gpp_actual_mean), values(gpp_pred_mean), na.rm = TRUE)
spatial_performance_old <- calc_performance(values(gpp_actual_legacy_mean), values(gpp_pred_legacy_mean), na.rm = TRUE)


