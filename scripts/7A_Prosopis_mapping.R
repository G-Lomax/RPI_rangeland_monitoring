## Evaluation against Prosopis-invaded land areas in Southern Kenya
## 25th October 2024
## Guy Lomax
## G.Lomax@exeter.ac.uk

## 1. Setup ----

source("scripts/load.R")

START <- 2009
END <- 2018

## 2. Load data ----

# RPI annual rasters

rpi_map <- rast("data/processed/raster/rpi_rast_sp.tif")

# Regional land cover map

lc_map <- rast("data/raw/raster/hunter2020/svmRadial_Multiple_season_time_series_Raster.tif")

lc_points <- read_csv("data/raw/csv/hunter2020/training_data_acess.csv") %>%
  st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326") %>%
  st_transform(crs = crs(lc_map))

lc_labels <- read_csv("data/raw/raster/hunter2020/hunter_raster_key.csv")

lc_points_labeled <- left_join(lc_points, lc_labels)

# Kenya/Tanzania polygon
ke_tz <- st_read("data/raw/vector/kenya_tanzania.geojson")
ke_tz$ADM0_NAME[ke_tz$ADM0_CODE == 257] <- "Tanzania"

ke_buffer <- filter(ke_tz, ADM0_NAME == "Kenya") %>%
  st_buffer(dist = set_units(50, "km"))

# Plot land cover map

lc_tmap <- tm_shape(lc_map) +
  tm_raster(
    col.scale = tm_scale_categorical(values = lc_labels$hex2, labels = lc_labels$label),
    col.legend = tm_legend(title = "Land cover", frame = FALSE)) +
  tm_shape(ke_tz) +
  tm_borders(lwd = 1.5) +
  tm_shape(lc_points_labeled) +
  tm_symbols(col = "white", fill = "label", fill.scale = tm_scale_discrete(values = "navyblue"),
             fill.legend = tm_legend(title = "", labels = "Training point"), size = 0.3, lwd = 0.6) +
  tm_graticules(lines = FALSE)

inset_map <- tm_shape(ke_buffer) +
  tm_shape(ke_tz) +
  tm_borders(lwd = 2) +
  tm_text("ADM0_NAME") +
  tm_shape(st_bbox(lc_map) %>% st_as_sfc()) +
  tm_borders(col = "red", lwd = 2) +
  tm_layout(bg.color = "white")

tmap_save(lc_tmap, "results/figures/hunter_lc_map.png",
          insets_tm = inset_map, insets_vp = viewport(x=0.65, y=0.8, width=0.15, height=0.25),
          width = 20, height = 16, units = "cm", dpi = 300)

## 3. Prepare data ----

# Create rasters of fractional land cover at RPI resolution
na_frac_map <- lc_map %>%
  is.na() %>%
  project(rpi_map, method = "average") %>%
  trim()

cover_frac_rasters <- map(1:9, function(x) {
  binary_map <- lc_map == x
  binary_map_reproj <- binary_map %>%
    project(rpi_map, method = "average") %>%
    trim() %>%
    mask(na_frac_map <= 0.2, maskvalues = 0)
}) %>% rast()

names(cover_frac_rasters) <- lc_labels$label

# RPI rasters

rpi_only <- subset(rpi_map, str_detect(names(rpi_map), "rpi"))

rpi_ts <- subset(rpi_only, str_detect(names(rpi_only), paste(START:END ,collapse = '|')))

rpi_2018 <- subset(rpi_only, "rpi.2018") %>%
  crop(cover_frac_rasters)

rpi_mean <- rpi_ts %>%
  crop(cover_frac_rasters) %>%
  mean()

# rpi_trend <- rpi_ts %>%
#   crop(cover_frac_rasters) %>%
#   app(function(x) {
#   if (any(is.na(x))) {
#     NA
#   } else {
#     mod <- deming::theilsen(x ~ seq_along(x))
#     
#     coef(mod)[2]
#   }
# })

names(rpi_2018) <- "rpi_2018"
names(rpi_mean) <- "rpi_mean"
# names(rpi_trend) <- "rpi_trend"

## 4. Assess RPI means and trends by Prosopis presence - linear mixture models ----

cover_frac_combined <- c(cover_frac_rasters, rpi_2018, rpi_mean)
cover_frac_df <- as.data.frame(cover_frac_combined, na.rm = TRUE, cells = TRUE, xy = TRUE)

# cover_frac_names <- rename_with(cover_frac_df, \(x) str_replace_all(x, " ", "_"))

# OLS multiple linear regression of RPI values against fractions of different
# land covers (exclude the largest land cover (Grass Control) to avoid
# collinearity - results will be relative to 100% Grass Control mean value)

# cover_frac_rpi_2018 <- lm(
#   rpi_2018 ~ `Vachellia spp.` + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
#     `Mixed Vegetation Control` + `Neltuma spp.` + `Sorghum bicolor` + `Sporobolus cordofanus`,
#   data = cover_frac_df
# )

cover_frac_rpi_mean <- lm(
  rpi_mean ~ 0 + `Vachellia spp.` + `Cynodon plectostachyus` + `Ficus sur` + `General Control` + `Grass Control` +
    `Mixed Vegetation Control` + `Neltuma spp.` + `Sorghum bicolor` + `Sporobolus cordofanus`,
  data = cover_frac_df
)

# cover_frac_rpi_gam <- gam(
#   rpi_mean ~ 0 + Vachellia_spp. + Cynodon_plectostachyus + Ficus_sur + General_Control + Grass_Control +
#     Mixed_Vegetation_Control + Neltuma_spp. + Sorghum_bicolor + Sporobolus_cordofanus,
#   method = "REML",
#   family = Gamma(link = "identity"),
#   data = cover_frac_names
# )

# draw(cover_frac_rpi_gam, parametric = TRUE, scales = "fixed", ci_col = "steelblue", ci_alpha = 0.6, residuals = TRUE) & theme_bw() & ylim(0, 1.5)

# cover_frac_rpi_partial <- lm(
#   rpi_mean ~ `Vachellia spp.` + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
#     `Mixed Vegetation Control` + `Sorghum bicolor` + `Sporobolus cordofanus`,
#   data = cover_frac_df
# )
# 
# cover_frac_rpi_added <- lm(
#   residuals(cover_frac_rpi_partial) ~ `Neltuma spp.`,
#   data = cover_frac_df
# )

# cover_frac_rpi_trend <- lm(
#   rpi_trend ~ `Vachellia spp.` + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
#     `Mixed Vegetation Control` + `Neltuma spp.` + `Sorghum bicolor` + `Sporobolus cordofanus`,
#   data = cover_frac_df
# )

# Partial residual plot showing slope/effect of Neltuma vs. other covers

neltuma_resids <- cover_frac_df %>%
  mutate(resids_mean= residuals(cover_frac_rpi_mean),
         partial_resids_mean =
           resids_mean +
           `Neltuma spp.` * coef(cover_frac_rpi_mean)[["`Neltuma spp.`"]])

partial_resid_plot <- ggplot(neltuma_resids, aes(x = `Neltuma spp.`, y = partial_resids_mean)) +
  geom_point(size = 0.5, alpha = 0.3) +
  # geom_abline(
  #   slope = coef(cover_frac_rpi_mean)[["`Neltuma spp.`"]],
  #   colour = "steelblue") +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlim(0, 1) +
  labs(x = "Land cover fraction - Neltuma spp.",
       y = "Partial residuals")

ggsave("results/figures/neltuma_partial_resids.png",
       width = 12, height = 9, units = "cm", dpi = 250)

# Coefficient plot

model_df <- as.data.frame(summary(cover_frac_rpi_mean)$coefficients) %>%
  select(-"Pr(>|t|)") %>%
  mutate(ci_upper = Estimate + 1.96 * `Std. Error`,
         ci_lower = Estimate - 1.96 * `Std. Error`)

model_df$covariate <- rownames(model_df)
rownames(model_df) <- NULL
model_df$neltuma <- ifelse(model_df$covariate == "`Neltuma spp.`", "Neltuma spp.", "Other")
model_df$label <- str_replace_all(model_df$covariate, "`", "")

parameter_plot <- ggplot(model_df, aes(y = reorder(label, Estimate))) +
  geom_point(aes(x = Estimate, colour = neltuma), size = 3, shape = 18, show.legend = FALSE) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper, colour = neltuma), show.legend = FALSE) +
  geom_vline(xintercept = c(0, 1)) +
  theme_bw() +
  xlim(0, 2) +
  # scale_y_discrete(limits = rev, labels = sort(lc_labels$label, TRUE)) +
  scale_colour_manual(values = c("red", "steelblue")) +
  labs(title = "a)",
       x = "Estimated pure\nRPI value", y = "Land cover")

ggsave("results/figures/lc_mean_parameter_plot.png",
       parameter_plot,
       width = 14, height = 10, units = "cm", dpi = 200)

## 5. Assess mean RPI at "pure" presence points ----

lc_points_sf <- st_as_sf(lc_points_labeled, coords = c("lon", "lat"), crs = "EPSG:4326")

lc_points_rpi <- extract(rpi_mean, lc_points_sf, cells = TRUE, bind = TRUE) %>%
  st_as_sf() %>%
  drop_na()

# Drop cells with more than one type of presence point within
lc_points_unique <- lc_points_rpi %>%
  group_by(cell) %>%
  filter(length(unique(Class)) == 1) %>%
  summarise(label = first(label),
            rpi_mean = mean(rpi_mean)) %>%
  mutate(geometry = st_centroid(geometry))

# Convert label to ordered factor based on mean pure RPI value above

lc_points_fct <- lc_points_unique %>%
  mutate(label = ordered(label, levels = levels(reorder(ordered(model_df$label), model_df$Estimate))))

field_obs_rpi_boxplot <- lc_points_fct %>%
  group_by(label) %>%
  mutate(mean = mean(rpi_mean),
         neltuma = ifelse(label == "Neltuma spp.", "Neltuma spp.", "Other")) %>%
  ungroup() %>%
  ggplot(aes(y = label, x = rpi_mean, colour = neltuma)) +
  geom_boxplot(coef = 6, width = 0.4, show.legend = FALSE) +
  geom_point(aes(x = mean), size = 3, shape = 18, show.legend = FALSE) +
  geom_vline(xintercept = c(0,1)) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 2)) +
  labs(title = "b)",
       x = "RPI value", y = "Land cover", colour = "") +
  # scale_y_discrete(limits = rev, labels = sort(lc_labels$label, TRUE)) +
  scale_colour_manual(values = c("red", "steelblue"))

ggsave("results/figures/field_obs_rpi_boxplot.png", field_obs_rpi_boxplot,
       width = 14, height = 10, units = "cm", dpi = 250)

# Combined figure

combined_rpi_estimates <- parameter_plot + field_obs_rpi_boxplot +
  plot_layout(ncol = 2, axes = "collect", guides = "collect", axis_titles = "collect")

ggsave("results/figures/neltuma_combined_rpi.png", combined_rpi_estimates,
       width = 20, height = 12, units = "cm", dpi = 250)

# ggplot(lc_points_unique, aes(x = Class, y = rpi_trend)) +
#   geom_hline(yintercept = 0, colour = "black") +
#   geom_boxplot() +
#   theme_bw() +
#   coord_cartesian(ylim = c(-0.15, 0.15)) +
#   labs(y = "RPI Trend") +
#   theme(axis.text.x = element_text(size = 10))

# Categorical model

cat_mod_mean <- aov(rpi_mean ~ label, data = lc_points_unique)
# cat_mod_trend <- aov(rpi_trend ~ Class, data = lc_points_unique)

TukeyHSD(cat_mod_mean)
