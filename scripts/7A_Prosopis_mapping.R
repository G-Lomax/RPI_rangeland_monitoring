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

## 3. Prepare data ----

lc_map_prosopis <- lc_map == 7

lc_map_prosopis_reproj <- lc_map_prosopis %>%
  project(rpi_map, method = "average") %>%
  trim()

rpi_only <- subset(rpi_map, str_detect(names(rpi_map), "rpi"))

rpi_ts <- subset(rpi_only, str_detect(names(rpi_only), paste(START:END ,collapse = '|')))

rpi_2018 <- subset(rpi_only, "rpi.2018") %>%
  crop(lc_map_prosopis_reproj)

rpi_mean <- rpi_ts %>%
  crop(lc_map_prosopis_reproj) %>%
  mean()

rpi_trend <- rpi_ts %>%
  crop(lc_map_prosopis_reproj) %>%
  app(function(x) {
  if (any(is.na(x))) {
    NA
  } else {
    mod <- deming::theilsen(x ~ seq_along(x))
    
    coef(mod)[2]
  }
})


## 4. Assess RPI means and trends by Prosopis presence

combined_rast <- c(lc_map_prosopis_reproj, rpi_2018, rpi_mean, rpi_trend)
names(combined_rast) <- c("prosopis_frac", "rpi_2018", "rpi_mean", "rpi_trend")

combined_df <- as.data.frame(combined_rast, na.rm = TRUE, cells = TRUE, xy = TRUE)

# Assess RPI means and trends of >80% pure pixels of different categories

na_frac_map <- lc_map %>%
  is.na() %>%
  project(rpi_map, method = "average") %>%
  crop(lc_map_prosopis_reproj)

cover_frac_rasters <- map(1:9, function(x) {
  binary_map <- lc_map == x
  binary_map_reproj <- binary_map %>%
    project(rpi_map, method = "average") %>%
    trim()
  
  # Remove cells with > 20% NA
  mask(binary_map_reproj, na_frac_map <= 0.2, maskvalues = 0)
}) %>% rast()

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
