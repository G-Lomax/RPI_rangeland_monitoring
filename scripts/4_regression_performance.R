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
  rename(quantile_pred_spt = "quantile= 0.9...17", quantile_pred_sp = "quantile= 0.9...18") %>%
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

importance_spt <- rf_tuned_spt$learner$importance()
importance_sp <- rf_tuned_sp$learner$importance()

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
