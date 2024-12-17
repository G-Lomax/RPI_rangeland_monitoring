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
  filter(diff(range(Latitude)) <= 0.005 & diff(range(Longitude)) <= 0.005) %>%  # Remove sites with ambiguous location data
  mutate(Latitude = median(Latitude), Longitude = median(Longitude)) %>%  # Ensure consistent lat-lon
  group_by(Name, ObservationDate_GMT, Transect, Direction) %>%
  mutate(n = n(),
         n_na = rowSums(across(starts_with("Intercept"), is.na)),
         cover_points = rowSums(across(ends_with("Count"), sum))) %>%
  filter(max(n) == 1) %>%  # Remove sites with more than one row per transect entry
  filter(max(n_na) == 0) %>%  # Remove sites with missing data in intercepts
  filter(sum(cover_points) <=  20) %>%  # Remove sites with anomalously large numbers of strikes recorded per intercept (95th percentile)
  group_by(Name, ObservationDate_GMT) %>%
  filter(n() == 20) %>%  # Remove sites/measurements with fewer than 20 rows
  filter(!(str_detect(Name, coll("test", TRUE)) | str_detect(Name, regex("practi.e", TRUE)))) %>%
  mutate(sample_id = cur_group_id())

# A. Record total % cover recorded of all cover types (can include multiple
# strikes per segment)
# landpks_total <- landpks_clean %>%
#   group_by(Name, Latitude, Longitude, ObservationDate_GMT) %>%
#   summarise(across(
#     .cols = ends_with("Count"),
#     .fns = sum,
#     .names = stringr::str_replace("{col}", "Segment", "")
#   ))


landpks_all_intercepts <- landpks_clean %>%
  select(-ends_with("Count")) %>%
  pivot_longer(cols = starts_with("Intercept"), names_to = "Intercept", values_to = "cover") %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(cover = str_split(cover, ",")) %>%
  unnest(cover) %>%
  mutate(cover = str_trim(cover))

# Fix data error
landpks_all_intercepts$cover[landpks_all_intercepts$cover == "HPerennial grasses"] <- "Perennial grasses"

landpks_full <- landpks_all_intercepts %>%
  left_join(landpks_lookup) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  reframe(group = unique(group)) %>%
  mutate(group = factor(group))

# Add new row for "bare" soil when no annual or perennial herbaceous plant listed
# (LandPKS method otherwise considers bare to be only fully exposed)
bare_soil_check <- landpks_full %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(group = check_bare(group, covers = c("annual", "perennial", "shrub"))) %>%
  filter(!is.na(group))

landpks_full_with_bare <- landpks_full %>%
  bind_rows(bare_soil_check) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT) %>%
  count(group, .drop = FALSE, name = "total_n")

# B. Record the first cover type intersected by a raindrop falling on the point
# (e.g., if a tree canopy is overhanging grass, record as tree)
# Calculate fractional cover as intersected by raindrop
landpks_raindrop <- landpks_all_intercepts %>%
  left_join(landpks_lookup) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(raindrop_intercept = group[which.min(raindrop_order)]) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT) %>%
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

# Where more than one plot occurs in an RPI/GPP grid cell for a given
# hydrological year, randomly sample one plot to retain

set.seed(212)

landpks_id_keep <- landpks_all_vars %>%
  group_by(cell, hYear, sample_id) %>%
  slice_head(n = 1) %>%
  group_by(cell, hYear) %>%
  reframe(sample_id = unique(sample_id)) %>%
  group_by(cell, hYear) %>%
  slice_sample(n = 1)

landpks_grid <- landpks_id_keep %>%
  left_join(landpks_all_vars) %>%
  st_as_sf() %>%
  group_by(Name, sample_id, cell, hYear, ObservationDate_GMT) %>%
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

rpi_gpp_raindrop <- landpks_grid %>%
  filter(group != "base") %>%
  rename("Potential GPP" = "potential_gpp_predicted", "RPI" = "rpi") %>%
  pivot_longer(cols = c("Potential GPP", "GPP", "RPI")) %>%
  ggplot(aes(x = raindrop_n, y = value)) +
  geom_point(size = 0.6) +
  geom_smooth(method = "lm") +
  facet_grid(cols = vars(group), rows = vars(name), scales = "free_y") +
  theme_bw()

# Relationships between RPI/GPP and land cover as viewed from the field

rpi_gpp_total <- landpks_grid %>%
  filter(group != "base") %>%
  rename("Potential GPP" = "potential_gpp_predicted", "RPI" = "rpi") %>%
  pivot_longer(cols = c("Potential GPP", "GPP", "RPI")) %>%
  ggplot(aes(x = total_n, y = value)) +
  geom_point(size = 0.6) +
  geom_smooth(method = "lm") +
  facet_grid(cols = vars(group), rows = vars(name), scales = "free_y") +
  theme_bw() +
  labs(x = "Fractional cover (%)",
       y = "Index value",
       title = "Relationship between indices and LandPKS land cover %")

# Correlations

rpi_gpp_cor <- landpks_grid %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarise(gpp_raindrop_cor = cor(raindrop_n, GPP),
            gpp_total_cor = cor(total_n, GPP),
            rpi_raindrop_cor = cor(raindrop_n, rpi),
            rpi_total_cor = cor(total_n, rpi),
            pot_raindrop_cor = cor(raindrop_n, potential_gpp_predicted),
            pot_total_cor = cor(total_n, potential_gpp_predicted))

rpi_gpp_cor_key_params <- landpks_grid %>%
  st_drop_geometry() %>%
  group_by(cell, ObservationDate_GMT) %>%
  summarise(rpi = mean(rpi),
            GPP = mean(GPP),
            pot_gpp = mean(potential_gpp_predicted),
            bare = mean(bare),
            herb_tot = mean(herb_tot),
            perennial = sum(total_n * group %in% c("perennial")),
            perennial_ratio = perennial / herb_tot) %>%
  ungroup() %>%
  summarise(gpp_bare_cor = cor(GPP, bare),
            rpi_bare_cor = cor(rpi, bare),
            pot_bare_cor = cor(pot_gpp, bare),
            gpp_per_cor = cor(GPP, perennial),
            rpi_per_cor = cor(rpi, perennial),
            pot_per_cor = cor(pot_gpp, perennial),
            gpp_per_ratio_cor = cor(GPP[herb_tot >= 20], perennial_ratio[herb_tot >= 20]),
            rpi_per_ratio_cor = cor(rpi[herb_tot >= 20], perennial_ratio[herb_tot >= 20]),
            pot_per_ratio_cor = cor(pot_gpp[herb_tot >= 20], perennial_ratio[herb_tot >= 20]))

# GPP-herbaceous cover relationship for different levels of woody cover
landpks_grid %>%
  mutate(woody_bin = cut_interval(woody, n = 5)) %>%
  ggplot(aes(x = herb, y = GPP)) +
  geom_point(size = 0.2) +
  facet_wrap(~woody_bin, nrow = 1, ncol = 5) +
  theme_bw() +
  geom_smooth(se = FALSE, method = "lm")

# GPP and RPI compared to annual vs. perennial balance

rpi_gpp_annual_perennial <- landpks_grid %>%
  st_drop_geometry() %>%
  select(-bare) %>%
  pivot_wider(id_cols = c("Name", "ObservationDate_GMT", "herb", "GPP", "rpi"), names_from = "group", values_from = "raindrop_n") %>%
  mutate(herb_bin = cut_interval(herb, n = 5),
         perennial_ratio = perennial / (annual + perennial)) %>%
  filter(herb > 0) %>%
  pivot_longer(cols = c("GPP", "rpi")) %>%
  ggplot(aes(x = perennial_ratio, y = value)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_grid(cols = vars(herb_bin), rows = vars(name), scales = "free") +
  theme_bw()

# Variation through the hydrological year

wet_start <- c(3, 10)
wet_end <- c(5, 12)
season_df <- data.frame(start = wet_start, end = wet_end, season = seq_along(wet_start),
                        day_start = wet_start * 30, day_end = wet_end * 30)


hYear_gpp_cor <- landpks_grid %>%
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
  geom_hline(yintercept = 0, colour = "grey20") +
  facet_wrap(~group) +
  theme_bw() +
  scale_x_continuous(breaks = 1:12) +
  geom_rect(inherit.aes = FALSE, data = season_df, aes(xmin = start, xmax = end, group = season),
            ymin = -1, ymax = 1, colour = "transparent", fill = "orange", alpha = 0.3)

cover_by_doy <- landpks_grid %>%
  mutate(yday = yday(ObservationDate_GMT)) %>%
  ggplot(aes(x = yday, y = total_n)) +
  geom_point(size = 0.5) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cc")) +
  facet_wrap(~group) +
  theme_bw() +
  geom_rect(inherit.aes = FALSE, data = season_df, aes(xmin = day_start, xmax = day_end, group = season),
            ymin = 0, ymax = 100, colour = "transparent", fill = "orange", alpha = 0.3)

## 6. Indices as predictive tools ----

landpks_model <- landpks_grid %>%
  group_by(cell, hYear, ObservationDate_GMT, days_in_hyear) %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  summarise(x = mean(x), y = mean(y), Name = first(Name),
            rpi = mean(rpi), GPP = mean(GPP), pot_gpp = mean(potential_gpp_predicted),
            bare = mean(bare),
            herb = mean(herb),
            woody = mean(woody),
            per_cover = sum(total_n * group %in% c("perennial")),
            per_frac = sum(total_n * group %in% c("perennial")) / sum(total_n * group %in% c("perennial", "annual"))
  )

# Total bare ground
bare_gam_raw <- gam(bare / 100 ~ s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(eps = 0.005))

bare_gam_rpi <- gam(bare / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.005)
)

bare_gam_gpp <- gam(bare/ 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                   data = landpks_model,
                   family = betar(link = "logit", eps = 0.005)
)

bare_gam_pot_gpp <- gam(bare/ 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                          te(x, y, k = 6, bs = "tp"),
                        data = landpks_model,
                        family = betar(link = "logit", eps = 0.005)
)

# Total perennial cover

per_gam_raw <- gam(per_cover / 100 ~ s(days_in_hyear, k = 12, bs = "cc") + 
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(eps = 0.001))

per_gam_rpi <- gam(per_cover / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

per_gam_gpp <- gam(per_cover / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

per_gam_pot_gpp <- gam(per_cover / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                     te(x, y, k = 6, bs = "tp"),
                   data = landpks_model,
                   family = betar(link = "logit", eps = 0.001)
)

# Perennial vs. annual balance

landpks_per_frac <- landpks_model %>%
  filter(woody <= 50 & herb >= 20)

per_ann_gam_raw <- gam(per_frac / 100 ~ s(days_in_hyear, bs = "cc") + 
                      te(x, y, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_rpi <- gam(per_frac / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_gpp <- gam(per_frac / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_pot_gpp <- gam(per_frac / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                         te(x, y, k = 6, bs = "tp"),
                       data = landpks_model,
                       family = betar(link = "logit", eps = 0.001)
)



# Perennial fraction total

per_gam_raw <- gam(per_cover / 100 ~ s(days_in_hyear, bs = "cc") + 
                         te(x, y, bs = "tp"),
                       data = landpks_model,
                       family = betar(link = "logit", eps = 0.001)
)

per_gam_rpi <- gam(per_cover / 100 ~ rpi + s(days_in_hyear, bs = "cc") +
                         te(x, y, bs = "tp"),
                       data = landpks_model,
                       family = betar(link = "logit", eps = 0.001)
)

per_gam_gpp <- gam(per_cover / 100 ~ GPP + s(days_in_hyear, bs = "cc") +
                         te(x, y, bs = "tp"),
                       data = landpks_model,
                       family = betar(link = "logit", eps = 0.001)
)

## 7. Trend consistency

# Filter LandPKS data to those with multiple years sampled in the same location

landpks_multitemporal <- landpks_model%>%
  group_by(cell) %>%
  filter(length(unique(hYear)) > 1)

landpks_trend <- landpks_multitemporal %>%
  group_by(cell) %>%
  mutate(centroid = find_centroid(geometry),
         dist_to_centroid = st_distance(geometry, centroid, by_element = TRUE)) %>%
  filter(dist_to_centroid <= set_units(100, "metres")) %>%
  filter(length(unique(hYear)) > 1) %>%
  summarise(names = list(unique(Name)), n = length(unique(hYear)),
            start = first(hYear), end = last(hYear), length = end - start + 1,
            centroid = first(centroid),
            max_distance = max(dist_to_centroid),
            bare_trend = (lm(bare ~ hYear) %>% coef())[2],
            per_trend = (lm(per_cover ~ hYear) %>% coef())[2],
            gpp_trend = (lm(GPP ~ hYear) %>% coef())[2],
            rpi_trend = (lm(rpi ~ hYear) %>% coef())[2]) %>%
  filter(length >= 3)

# GPP vs perennial

ggplot(landpks_trend, aes(x = gpp_trend, y = per_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits)

# GPP vs bare ground
ggplot(landpks_trend, aes(x = gpp_trend, y = bare_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits)

# RPI vs perennial
ggplot(landpks_trend, aes(x = rpi_trend, y = per_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits)

# RPI vs bare ground
ggplot(landpks_trend, aes(x = rpi_trend, y = bare_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits)
