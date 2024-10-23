## Evaluation against LandPKS field monitoring dataset
## 1st October 2024
## Guy Lomax
## G.Lomax@exeter.ac.uk

## 1. Setup ----

source("scripts/load.R")

START <- date("2000-09-01")
END <- date("2023-08-31")

## 2. Load data ----

# Field land cover data
landpks_raw <- read_csv("data/raw/csv/LandPKS/landCover.csv")

# Lookup table for summarising data
landpks_lookup <- read_csv("data/raw/csv/LandPKS/landpks_lookup.csv")

# RPI and covariate rasters

rpi_map <- rast("data/processed/raster/rpi_rast_sp.tif")

static_vars <- rast("data/raw/raster/covariateMaps/staticVars.tif") %>%
  crop(rpi_map[[1]])

precipitation <- Sys.glob("data/raw/raster/covariateMaps/dynamicVars*.tif") %>%
  map(rast) %>%
  map(function(r) subset(r, str_detect(names(r), "precipitation"))) %>%
  rast() %>%
  crop(rpi_map)

names(precipitation) <- paste0(names(precipitation), ".", 2000:2022)

map <- mean(precipitation)
names(map) <- "map"

## 3. Data cleaning and preparation -----

# Select desired columns and limit to study area
landpks_cols <- select(
  landpks_raw,
  c("Name", "Latitude", "Longitude", "ObservationDate_GMT",
    "Transect", "Direction", starts_with("Intercept"), ends_with("Count"),
    -"SegmentSpecies1Density", -"SegmentSpecies2Density")
)

landpks_study_area <- landpks_cols %>%
  filter(Longitude >= xmin(rpi_map) & Longitude <= xmax(rpi_map)) %>%
  filter(Latitude >= ymin(rpi_map) & Latitude <= ymax(rpi_map)) %>%
  filter(ObservationDate_GMT >= START & ObservationDate_GMT <= END)

# Remove duplicate rows
# Remove sites with inconsistent coordinates
# Remove sites with missing transect data or more than one row per
# transect entry

landpks_clean <- landpks_study_area %>%
  distinct() %>%  # Remove duplicate rows
  group_by(Name, ObservationDate_GMT) %>%
  filter(length(unique(Latitude)) == 1 & length(unique(Longitude)) == 1) %>%  # Remove sites with ambiguous location data
  group_by(Name, ObservationDate_GMT, Transect, Direction) %>%
  mutate(n = n(),
         n_na = rowSums(across(starts_with("Intercept"), is.na)),
         cover_points = rowSums(across(ends_with("Count"), sum))) %>%
  filter(max(n) == 1) %>%  # Remove sites with more than one row per transect entry
  filter(max(n_na) == 0) %>%  # Remove sites with missing data in intercepts
  filter(sum(cover_points) <=  20) %>%  # Remove sites with anomalously large numbers of strikes recorded per intercept (95th percentile)
  group_by(Name, ObservationDate_GMT) %>%
  filter(n() == 20) %>%  # Remove sites/measurements with fewer than 20 rows
  filter(!(str_detect(Name, coll("test", TRUE)) | str_detect(Name, regex("practi.e", TRUE))))

# A. Record total % cover recorded of all cover types (can include multiple
# strikes per segment)
landpks_total <- landpks_clean %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT) %>%
  summarise(across(
    .cols = ends_with("Count"),
    .fns = sum,
    .names = stringr::str_replace("{col}", "Segment", "")
  ))


landpks_all_intercepts <- landpks_clean %>%
  select(-ends_with("Count")) %>%
  pivot_longer(cols = starts_with("Intercept"), names_to = "Intercept", values_to = "cover") %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(cover = str_split(cover, ",")) %>%
  unnest(cover) %>%
  mutate(cover = str_trim(cover))

# Fix data error
landpks_all_intercepts$cover[landpks_all_intercepts$cover == "HPerennial grasses"] <- "Perennial grasses"

landpks_full <- landpks_all_intercepts %>%
  left_join(landpks_lookup) %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  reframe(group = unique(group)) %>%
  mutate(group = factor(group))

# Add new row for "bare" soil when no annual or perennial herbaceous plant listed
# (LandPKS method otherwise considers bare to be only fully exposed)
bare_soil_check <- landpks_full %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(group = check_bare(group, covers = c("annual", "perennial", "shrub"))) %>%
  filter(!is.na(group))

landpks_full_with_bare <- landpks_full %>%
  bind_rows(bare_soil_check) %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT) %>%
  count(group, .drop = FALSE, name = "total_n")

# B. Record the first cover type intersected by a raindrop falling on the point
# (e.g., if a tree canopy is overhanging grass, record as tree)
# Calculate fractional cover as intersected by raindrop
landpks_raindrop <- landpks_all_intercepts %>%
  left_join(landpks_lookup) %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(raindrop_intercept = group[which.min(raindrop_order)]) %>%
  group_by(Name, Latitude, Longitude, ObservationDate_GMT) %>%
  mutate(raindrop_intercept = factor(raindrop_intercept)) %>%
  count(raindrop_intercept, .drop = FALSE, name = "raindrop_n") %>%
  rename(group = raindrop_intercept)

landpks_all <- inner_join(landpks_raindrop, landpks_full_with_bare)

# Convert to sf object

landpks_sf <- st_as_sf(landpks_all, coords = c("Longitude", "Latitude"),  crs = "EPSG:4326")

## 4. Extract raster values for points ----

# Convert Observation Date to hydrological year and day in hYear

landpks_hyear <- landpks_sf %>%
  mutate(time_since_start = interval(START, ObservationDate_GMT),
         hYear = year(START) + time_since_start %/% years(),
         days_in_hyear = time_since_start %% years() %/% days())

# Extract relevant layers

data_years <- unique(landpks_hyear$hYear) %>% sort()

annual_vars <- c(rpi_map, precipitation)
static_vars_with_map <- c(static_vars, map)

landpks_rpi_gpp <- map(data_years, raster_extract, sf = landpks_hyear, raster = annual_vars) %>%
  bind_rows() %>%
  drop_na()

landpks_all_vars <- bind_cols(landpks_rpi_gpp, terra::extract(static_vars_with_map, vect(landpks_rpi_gpp)))

# Take mean of plots in the same grid cell that occur in the same hydrological year

landpks_average <- landpks_all_vars %>%
  group_by(group, cell, hYear) %>%
  summarise(Name = first(Name),
            across(
              .cols = c(ObservationDate_GMT, raindrop_n, total_n,
                        time_since_start, days_in_hyear,
                        GPP, gpp_predicted, potential_gpp_predicted, rpi,
                        precipitation, treeCover),
              .fns = mean
            )) %>%
  group_by(Name, ObservationDate_GMT) %>%
  mutate(woody = sum(raindrop_n * group %in% c("tree", "shrub")),
         woody_tot = sum(total_n * group %in% c("tree", "shrub")),
         herb = sum(raindrop_n * group %in% c("annual", "perennial")),
         herb_tot = sum(total_n * group %in% c("annual", "perennial")),
         bare = sum(raindrop_n * group %in% c("bare")),
         bare_ground = sum(total_n * group %in% c("bare"))) %>%
  ungroup() %>%
  st_centroid()

## 5. EDA ----

# Relationships between RPI/GPP and land cover as viewed from above

rpi_gpp_raindrop <- landpks_average %>%
  filter(group != "base") %>%
  pivot_longer(cols = c("GPP", "rpi", "precipitation")) %>%
  ggplot(aes(x = raindrop_n, y = value)) +
  geom_point(size = 0.2) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(cols = vars(group), rows = vars(name), scales = "free_y") +
  theme_bw()

# Relationships between RPI/GPP and land cover as viewed from the field

rpi_gpp_total <- landpks_average %>%
  pivot_longer(cols = c("GPP", "rpi", "precipitation")) %>%
  ggplot(aes(x = total_n, y = value)) +
  geom_point(size = 0.2) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(cols = vars(group), rows = vars(name), scales = "free_y") +
  theme_bw()

# Correlations

rpi_gpp_cor <- landpks_average %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarise(gpp_raindrop_cor = cor(raindrop_n, GPP),
            gpp_total_cor = cor(total_n, GPP),
            rpi_raindrop_cor = cor(raindrop_n, rpi),
            rpi_total_cor = cor(total_n, rpi),
            ppt_raindrop_cor = cor(raindrop_n, precipitation),
            ppt_total_cor = cor(total_n, precipitation))

# GPP-herbaceous cover relationship for different levels of woody cover
landpks_average %>%
  mutate(woody_bin = cut_interval(woody, n = 5)) %>%
  ggplot(aes(x = herb, y = GPP)) +
  geom_point(size = 0.2) +
  facet_wrap(~woody_bin, nrow = 1, ncol = 5) +
  theme_bw() +
  geom_smooth(se = FALSE, method = "lm")

# GPP and RPI compared to annual vs. perennial balance

rpi_gpp_annual_perennial <- landpks_average %>%
  st_drop_geometry() %>%
  select(-bare) %>%
  pivot_wider(id_cols = c("Name", "ObservationDate_GMT", "herb", "GPP", "rpi"), names_from = "group", values_from = "raindrop_n") %>%
  mutate(herb_bin = cut_interval(herb, n = 5),
         perennial_ratio = log((perennial + 0.5) / (annual + 0.5))) %>%
  filter(herb > 0) %>%
  pivot_longer(cols = c("GPP", "rpi")) %>%
  ggplot(aes(x = perennial_ratio, y = value)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(cols = vars(herb_bin), rows = vars(name), scales = "free") +
  theme_bw()

# Variation through the hydrological year

hYear_gpp_cor <- landpks_average %>%
  mutate(month = month(ObservationDate_GMT)) %>%
  group_by(month, group) %>%
  filter(n() >= 10) %>%
  summarise(gpp_raindrop_cor = cor(raindrop_n, GPP),
            gpp_total_cor = cor(total_n, GPP),
            rpi_raindrop_cor = cor(raindrop_n, rpi),
            rpi_total_cor = cor(total_n, rpi)) %>%
  pivot_longer(ends_with("_cor")) %>%
  ggplot(aes(x = month, y = value, colour = name)) +
  geom_line() +
  facet_wrap(~group) +
  theme_bw()

cover_by_doy <- landpks_average %>%
  mutate(yday = yday(ObservationDate_GMT)) %>%
  ggplot(aes(x = yday, y = total_n)) +
  geom_point(size = 0.5) +
  geom_smooth(se = FALSE, method = "gam") +
  facet_wrap(~group) +
  theme_bw()

## 6. Indices as predictive tools ----

landpks_model <- landpks_average %>%
  group_by(cell, hYear, ObservationDate_GMT, days_in_hyear) %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  summarise(x = mean(x), y = mean(y),
            rpi = mean(rpi), GPP = mean(GPP),
            bare = mean(bare),
            herb = mean(herb),
            herb_tot = sum(total_n * group %in% c("annual", "perennial")),
            per_ann_balance = sum(total_n * group %in% c("perennial")) - sum(total_n * group %in% c("annual"))
  ) %>%
  st_drop_geometry()

# Total bare ground
bare_gam_raw <- gam(bare / 100 ~ s(days_in_hyear, bs = "cc") + 
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = betar(eps = 0.001))

bare_gam_rpi <- gam(bare / 100 ~ s(rpi, bs = "ts") + s(days_in_hyear, bs = "cc") +
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

bare_gam_gpp <- gam(bare/ 100 ~ s(GPP, bs = "ts") + s(days_in_hyear, bs = "cc") +
                     te(x, y, bs = "ts"),
                   data = landpks_model,
                   family = betar(link = "logit", eps = 0.001)
)

# Total herbaceous cover

herb_gam_raw <- gam(herb / 100 ~ s(days_in_hyear, bs = "cc") + 
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = betar(eps = 0.001))

herb_gam_rpi <- gam(herb / 100 ~ s(rpi, bs = "ts") + s(days_in_hyear, bs = "cc") +
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

herb_gam_gpp <- gam(herb/ 100 ~ s(GPP, bs = "ts") + s(days_in_hyear, bs = "cc") +
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

# Perennial vs. annual balance

per_ann_gam_raw <- gam(per_ann_balance / 100 ~ s(days_in_hyear, bs = "cc") + 
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = "gaussian"
)

per_ann_gam_rpi <- gam(per_ann_balance / 100 ~ s(rpi, bs = "ts") + s(days_in_hyear, bs = "cc") +
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = "gaussian"
)

per_ann_gam_gpp <- gam(per_ann_balance/ 100 ~ s(GPP, bs = "ts") + s(days_in_hyear, bs = "cc") +
                      te(x, y, bs = "ts"),
                    data = landpks_model,
                    family = "gaussian"
)
