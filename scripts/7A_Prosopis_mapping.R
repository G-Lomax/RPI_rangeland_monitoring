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

lc_points <- read_csv("data/raw/csv/hunter2020/training_data_acess.csv")

lc_labels <- read_csv("data/raw/raster/hunter2020/hunter_raster_key.csv")

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

rpi_trend <- rpi_ts %>%
  crop(cover_frac_rasters) %>%
  app(function(x) {
  if (any(is.na(x))) {
    NA
  } else {
    mod <- deming::theilsen(x ~ seq_along(x))
    
    coef(mod)[2]
  }
})

names(rpi_2018) <- "rpi_2018"
names(rpi_mean) <- "rpi_mean"
names(rpi_trend) <- "rpi_trend"


## 4. Assess RPI means and trends by Prosopis presence ----

cover_frac_combined <- c(cover_frac_rasters, rpi_2018, rpi_mean, rpi_trend)
cover_frac_df <- as.data.frame(cover_frac_combined, na.rm = TRUE, cells = TRUE, xy = TRUE)

# OLS multiple linear regression of RPI values against fractions of different
# land covers (exclude the largest land cover (Grass Control) to avoid
# collinearity - results will be relative to Grass Control mean value)

cover_frac_rpi_2018 <- lm(
  rpi_2018 ~ Acacia + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
    `Mixed Vegetation Control` + Prosopis + `Sorghum bicolor` + `Sporobolus cordofanus`,
  data = cover_frac_df
)

cover_frac_rpi_mean <- lm(
  rpi_mean ~ Acacia + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
    `Mixed Vegetation Control` + Prosopis + `Sorghum bicolor` + `Sporobolus cordofanus`,
  data = cover_frac_df
)

cover_frac_rpi_trend <- lm(
  rpi_trend ~ Acacia + `Cynodon plectostachyus` + `Ficus sur` + `General Control` +
    `Mixed Vegetation Control` + Prosopis + `Sorghum bicolor` + `Sporobolus cordofanus`,
  data = cover_frac_df
)

# Assess RPI means and trends of >80% pure pixels of different categories



cover_pure <- cover_frac_rasters %>%
  mask(any(cover_frac_rasters >= 0.8), maskvalues = 0) %>%
  which.max()


pure_combined <- c(cover_pure, rpi_mean, rpi_trend)
names(pure_combined) <- c("cover_pure", "rpi_mean", "rpi_trend")

pure_df <- as.data.frame(pure_combined, na.rm = TRUE)
ggplot(pure_df, aes(x = as.factor(cover_pure), y = rpi_mean)) +
  geom_boxplot()


# Assess by presence points

lc_points_sf <- st_as_sf(lc_points, coords = c("lon", "lat"), crs = "EPSG:4326")

lc_points_rpi <- extract(combined_rast[[2:4]], lc_points_sf, cells = TRUE, bind = TRUE) %>%
  st_as_sf() %>%
  drop_na()

# Drop cells with more than one type of presence point within
lc_points_unique <- lc_points_rpi %>%
  group_by(cell) %>%
  filter(length(unique(Class)) == 1) %>%
  summarise(Class = first(Class),
            rpi_2018 = mean(rpi_2018),
            rpi_mean = mean(rpi_mean),
            rpi_trend = mean(rpi_trend)) %>%
  mutate(geometry = st_centroid(geometry))

ggplot(lc_points_unique, aes(x = Class, y = rpi_mean)) +
  geom_boxplot() +
  geom_hline(yintercept = c(0,1), colour = "black") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(y = "Mean RPI") +
  theme(axis.text.x = element_text(size = 10))

ggplot(lc_points_unique, aes(x = Class, y = rpi_trend)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_boxplot() +
  theme_bw() +
  coord_cartesian(ylim = c(-0.15, 0.15)) +
  labs(y = "RPI Trend") +
  theme(axis.text.x = element_text(size = 10))

# Categorical model

cat_mod_mean <- aov(rpi_mean ~ Class, data = lc_points_unique)
cat_mod_trend <- aov(rpi_trend ~ Class, data = lc_points_unique)


TukeyHSD(cat_mod_trend)
