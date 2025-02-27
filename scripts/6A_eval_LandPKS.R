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

# precipitation <- Sys.glob("data/raw/raster/covariateMaps/dynamicVars*.tif") %>%
#   map(rast) %>%
#   map(function(r) subset(r, str_detect(names(r), "precipitation"))) %>%
#   rast() %>%
#   crop(rpi_map)
# 
# names(precipitation) <- paste0(names(precipitation), ".", 2000:2022)
# 
# map <- mean(precipitation)
# names(map) <- "map"

# Country polygons
ke_tz <- st_read("data/raw/vector/kenya_tanzania.geojson")

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
  mutate(cover = str_trim(cover)) %>%
  ungroup()

# Fix data error
landpks_all_intercepts$cover[landpks_all_intercepts$cover == "HPerennial grasses"] <- "Perennial grasses"

landpks_full <- landpks_all_intercepts %>%
  left_join(landpks_lookup) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  reframe(group = unique(group)) %>%
  mutate(group = factor(group)) %>%
  ungroup()

# Add new row for "bare" soil when no annual or perennial herbaceous plant listed
# (LandPKS method otherwise considers bare to be only fully exposed)
bare_soil_check <- landpks_full %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT, Transect, Direction, Intercept) %>%
  summarise(group = check_bare(group, covers = c("annual", "perennial"))) %>%
  filter(!is.na(group)) %>%
  ungroup()

landpks_full_with_bare <- landpks_full %>%
  bind_rows(bare_soil_check) %>%
  group_by(Name, sample_id, Latitude, Longitude, ObservationDate_GMT) %>%
  count(group, .drop = FALSE, name = "total_n") %>%
  ungroup()

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
  rename(group = raindrop_intercept) %>%
  ungroup()

landpks_all <- inner_join(landpks_raindrop, landpks_full_with_bare)

# Convert to sf object

landpks_sf <- st_as_sf(landpks_all, coords = c("Longitude", "Latitude"),  crs = "EPSG:4326")

# Add zero rows for each sample_id for land covers not present
landpks_sf_complete <- landpks_sf %>%
  complete(sample_id, group, fill = list(raindrop_n = 0, total_n = 0)) %>%
  group_by(sample_id) %>%
  mutate(Name = first(Name[!is.na(Name)]),
         ObservationDate_GMT = first(ObservationDate_GMT[!is.na(ObservationDate_GMT)]),
         geometry = st_centroid(st_union(geometry))) %>%
  ungroup() %>%
  st_as_sf(crs = "EPSG:4326")

## 4. Extract raster values for points ----

# Convert Observation Date to hydrological year and day in hYear

landpks_hyear <- landpks_sf_complete %>%
  mutate(time_since_start = interval(START, ObservationDate_GMT),
         hYear = year(START) + time_since_start %/% years(),
         days_in_hyear = time_since_start %% years() %/% days())

# Extract relevant layers

data_years <- unique(landpks_hyear$hYear) %>% sort()

landpks_rpi_gpp <- map(data_years, raster_extract, sf = landpks_hyear, raster = rpi_map) %>%
  bind_rows() %>%
  drop_na() %>%
  select(Name, sample_id, cell,
         ObservationDate_GMT, time_since_start, hYear, days_in_hyear,
         group, raindrop_n, total_n,
         GPP, potential_gpp_predicted, rpi,
         geometry)

# Where more than one plot occurs in an RPI/GPP grid cell for a given
# hydrological year, take the mean

# Ignore warning of non-numeric values - doesn't affect outcome
landpks_id_mean <- landpks_rpi_gpp %>%
  select(-Name, -sample_id) %>%
  group_by(cell, hYear, group) %>%
  summarise(across(everything(), mean)) %>%
  ungroup()

landpks_grid <- landpks_id_mean %>%
  group_by(cell, hYear, ObservationDate_GMT) %>%
  mutate(woody = sum(raindrop_n * group %in% c("tree", "shrub")),
         woody_tot = sum(total_n * group %in% c("tree", "shrub")),
         herb = sum(raindrop_n * group %in% c("annual", "perennial")),
         herb_tot = sum(total_n * group %in% c("annual", "perennial")),
         bare = sum(raindrop_n * group %in% c("bare")),
         bare_ground = sum(total_n * group %in% c("bare"))) %>%
  ungroup() %>%
  st_centroid()

# Plot retained points on map

landpks_locs <- landpks_grid %>%
  filter(group == unique(group)[1])

# landpks_map <- tm_shape(ke_tz, is.main = FALSE) +
#   tm_polygons(fill = "grey95") +
#   tm_shape(rpi_map$GPP.2000, is.main = TRUE) +
#   tm_raster(col.scale = tm_scale_continuous(values = "wheat"),
#             col.legend = tm_legend(show = FALSE)) +
#   tm_shape(landpks_locs) +
#   tm_symbols(
#     size = 0.3,
#     fill = "hYear",
#     fill.scale = tm_scale_discrete(
#       values = "paired", 
#       label.format = list(fun=function(x) formatC(x, digits=0, format="d"))),
#     fill.legend = tm_legend(title = "Year")) +
#   tm_shape(ke_tz) +
#   tm_borders() +
#   tm_layout(frame = FALSE, legend.frame = FALSE)
# 
# `tmap`_save(landpks_map, "results/figures/landpks_site_map.png",
#           height = 16, width = 12, units = "cm", dpi = 250)

# Facet map

rpi_lpks <- rpi_map[[paste0("rpi.", as.character(2014:2022))]]

# landpks_to_map <- landpks_locs %>%
#   mutate(hYear = paste0("rpi.", hYear))
# 
# lpks_facet_map <- tm_shape(rpi_lpks) +
#   tm_raster(
#     col.scale = tm_scale_continuous(
#       values = "viridis", limits = c(0.2, 1.2), outliers.trunc = c(TRUE, TRUE)),
#     col.legend = tm_legend(show = FALSE)) +
#   tm_shape(landpks_to_map) +
#   tm_symbols(col = "white", fill = "red", size = 0.8) +
#   tm_facets_wrap(by = "hYear", ncol = 3, nrow = 3, labels = 2014:2022) +
#   tm_shape(ke_tz) +
#   tm_borders(lwd = 1.5) +
#   tm_layout(panel.label.bg.color = "white", panel.label.frame = FALSE, frame = FALSE,
#             panel.labels = 2014:2022, panel.label.size = 12)
# 
# tmap_save(lpks_facet_map, "results/figures/landpks_facet_map.png",
#           width = 16, height = 24, units = "cm", dpi = 250)

# With grid showing density of points

density_grid <- rpi_map[[1]] %>% aggregate(100, fun = "mean", na.rm = TRUE)

lpks_raster <- map(2014:2022, function(n) {
  
  landpks_year <- filter(landpks_rpi_gpp, hYear == n) %>%
    group_by(sample_id) %>%
    summarise(hYear = mean(hYear))
  
  landpks_density <- rasterize(landpks_year, density_grid, field = "sample_id", fun = "count")
  
  landpks_density_mask <- mask(landpks_density, landpks_density, maskvalues = c(0, NA))
  
  names(landpks_density_mask) <- paste0("Y", n)
  
  landpks_density_mask
}) %>% rast()

year_n <- landpks_rpi_gpp %>%
  st_drop_geometry() %>%
  group_by(sample_id) %>%
  summarise(hYear = mean(hYear)) %>%
  count(hYear)

panel_labels <- paste0(year_n$hYear, " (n=", year_n$n, ")")

lpks_grid_maps <- tm_shape(ke_tz, is.main = FALSE) +
  tm_fill("grey90") +
  tm_shape(rpi_lpks, is.main = TRUE) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0.2, 1.2),
      values = "-yl_gn_bu",
      outliers.trunc = c(TRUE, TRUE)),
    col.legend = tm_legend(show = FALSE),
    col_alpha = 0.6) +
  tm_shape(lpks_raster) +
  tm_raster(
    col.scale = tm_scale_continuous_log10(
      limits = c(1, 200),
      values = "magma",
      outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(show = FALSE)) +
  tm_shape(ke_tz) +
  tm_borders(lwd = 2) +
  tm_facets_wrap(ncol = 3, nrow = 3) +
  tm_layout(frame = FALSE, panel.label.frame = FALSE,
            panel.label.bg.color = "transparent",
            panel.labels = panel_labels, panel.label.size = 10,
            inner.margins = rep(0.05, 4))

tmap_save(lpks_grid_maps, "results/figures/landpks_grid_density_maps.png",
          width = 16, height = 24, units = "cm", dpi = 250)

rpi_legend <- tm_shape(rpi_lpks[[1]]) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0.2, 1.2),
      values = "-yl_gn_bu",
      outliers.trunc = c(TRUE, TRUE)),
    col.legend = tm_legend(reverse = TRUE, title = "RPI", frame = FALSE),
    col_alpha = 0.6) +
  tm_layout(legend.only = TRUE)

lpks_legend <- tm_shape(lpks_raster[[1]]) +
  tm_raster(
    col.scale = tm_scale_continuous_log10(
      limits = c(1, 200),
      ticks = c(1, 10, 100),
      labels = c(1, 10, 100),
      values = "magma",
      outliers.trunc = c(TRUE, TRUE)
    ),
    col.legend = tm_legend(reverse = TRUE, title = "N", frame = FALSE)) +
  tm_layout(legend.only = TRUE)

tmap_save(rpi_legend, "results/figures/landpks_grid_rpi_legend.png",
          width = 4, height = 12, units = "cm", dpi = 250)

tmap_save(lpks_legend, "results/figures/landpks_grid_density_legend.png",
          width = 4, height = 12, units = "cm", dpi = 250)

## 5. EDA ----

# Relationships between RPI/GPP and land cover as viewed from above

cap_labels <- capitalise_string(unique(landpks_grid$group))
names(cap_labels) <- unique(landpks_grid$group)

rpi_gpp_raindrop <- landpks_grid %>%
  filter(group != "base" & group != "litter") %>%
  rename("Potential GPP" = "potential_gpp_predicted", "RPI" = "rpi") %>%
  pivot_longer(cols = c("Potential GPP", "GPP", "RPI")) %>%
  ggplot(aes(x = raindrop_n, y = value)) +
  geom_point(size = 0.6) +
  geom_smooth(method = "lm") +
  facet_grid(name ~ group, scales = "free_y", labeller = labeller(group = as_labeller(cap_labels))) +
  theme_bw() +
  labs(x = "Fractional cover", y = "Metric value")

# Relationships between RPI/GPP and land cover as viewed from the field

rpi_gpp_total <- landpks_grid %>%
  filter(group != "base" & group != "litter") %>%
  rename("Potential GPP" = "potential_gpp_predicted", "RPI" = "rpi") %>%
  pivot_longer(cols = c("Potential GPP", "GPP", "RPI")) %>%
  ggplot(aes(x = total_n, y = value)) +
  geom_point(size = 0.6) +
  geom_smooth(method = "gam", alpha = 0.4, colour = "darkgreen", fill = "green3") +
  facet_grid(name ~ group, scales = "free_y", labeller = labeller(group = as_labeller(cap_labels))) +
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

rpi_gpp_key_params <- landpks_grid %>%
  st_drop_geometry() %>%
  group_by(cell, ObservationDate_GMT) %>%
  summarise(rpi = mean(rpi),
            GPP = mean(GPP),
            pot_gpp = mean(potential_gpp_predicted),
            bare = mean(bare),
            herb_tot = mean(herb_tot),
            perennial = sum(total_n * group %in% c("perennial")),
            perennial_ratio = perennial / herb_tot) %>%
  ungroup()

param_labels = c(bare = "% bare ground", herb_tot = "% herbaceous cover", perennial_ratio = "Perennial fraction")
index_labels = c(rpi = "RPI", GPP = "Observed GPP", pot_gpp = "Modelled potential GPP")
all_labels = c(param_labels, index_labels)

## NEED TO FIX SO WE DON'T GET > 100% TOTAL HERB COVER

key_params_plot <- rpi_gpp_key_params %>%
  pivot_longer(cols = c("bare", "herb_tot", "perennial_ratio"), names_to = "param", values_to = "param_value") %>%
  pivot_longer(cols = c("rpi", "GPP", "pot_gpp"), names_to = "index", values_to = "index_value") %>%
  ggplot(aes(x = param_value, y = index_value)) +
  geom_point(size = 0.6) +
  geom_smooth(method = "gam", alpha = 0.2, col = "darkgreen", fill = "forestgreen") +
  facet_grid(index ~ param, scales = "free", labeller = as_labeller(all_labels), switch = "both") +
  theme_bw() +
  labs(x = "", y = "") +
  theme(strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.placement = "outside",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggh4x::facetted_pos_scales(
    x = list(
      param == "bare" ~ scale_x_continuous(limits = c(0, 100)),
      param == "herb_tot" ~ scale_x_continuous(limits = c(0, 100)),
      param == "perennial_ratio" ~ scale_x_continuous(limits = c(0, 1))),
    y = list(
      index == "rpi" ~ scale_y_continuous(limits = c(0, 1.5)),
      index == "GPP" ~ scale_y_continuous(limits = c(0, 6000)),
      index == "rsq" ~ scale_y_continuous(limits = c(0, 6000))
  ))

ggsave("results/figures/index_field_scatter_plots.png", key_params_plot,
       width = 24, height = 30, units = "cm", dpi = 250)

rpi_gpp_key_params_cor <- rpi_gpp_key_params %>%
  summarise(gpp_bare_cor = cor(GPP, bare),
            rpi_bare_cor = cor(rpi, bare),
            pot_bare_cor = cor(pot_gpp, bare),
            gpp_herb_cor = cor(GPP, herb_tot),
            rpi_herb_cor = cor(rpi, herb_tot),
            pot_herb_cor = cor(pot_gpp, herb_tot),
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
  pivot_wider(id_cols = c("cell", "ObservationDate_GMT", "herb", "GPP", "rpi"), names_from = "group", values_from = "raindrop_n") %>%
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

cover_by_doy <- landpks_grid %>%
  mutate(yday = yday(ObservationDate_GMT)) %>%
  ggplot(aes(x = yday, y = total_n)) +
  geom_point(size = 0.5) +
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cc")) +
  facet_wrap(~group) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.placement = "outside",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

## 6. Indices as predictive tools ----

landpks_model <- landpks_grid %>%
  group_by(cell, hYear, ObservationDate_GMT, days_in_hyear) %>%
  mutate(x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2]) %>%
  summarise(x = mean(x), y = mean(y), cell = mean(cell),
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

bare_gam_all <- gam(bare / 100 ~ s(rpi, k = 6, bs = "tp") + 
                      s(GPP, k = 6, bs = "tp") +
                      s(pot_gpp, k = 6, bs = "tp") + 
                      s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.005)
)

bare_gam_rpi2 <- gam(bare / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                     data = landpks_model,
                     family = betar(link = "logit", eps = 0.001)
)

bare_gam_gpp2 <- gam(bare / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                     data = landpks_model,
                     family = betar(link = "logit", eps = 0.001)
)

bare_gam_pot_gpp2 <- gam(bare / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                         data = landpks_model,
                         family = betar(link = "logit", eps = 0.001)
)

# Total herbaceous cover

herb_gam_raw <- gam(herb / 100 ~ s(days_in_hyear, k = 12, bs = "cc") + 
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(eps = 0.001))

herb_gam_rpi <- gam(herb / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

herb_gam_gpp <- gam(herb / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

herb_gam_pot_gpp <- gam(herb / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc") +
                     te(x, y, k = 6, bs = "tp"),
                   data = landpks_model,
                   family = betar(link = "logit", eps = 0.001)
)

herb_gam_all <- gam(herb / 100 ~ s(rpi, k = 6, bs = "tp") + 
                      s(GPP, k = 6, bs = "tp") +
                      s(pot_gpp, k = 6, bs = "tp") + 
                      s(days_in_hyear, k = 12, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.005)
)

herb_gam_rpi2 <- gam(herb / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

herb_gam_gpp2 <- gam(herb / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                    data = landpks_model,
                    family = betar(link = "logit", eps = 0.001)
)

herb_gam_pot_gpp2 <- gam(herb / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                        data = landpks_model,
                        family = betar(link = "logit", eps = 0.001)
)

# Perennial vs. annual balance

landpks_per_frac <- landpks_model %>%
  filter(woody <= 50 & herb >= 20)

per_ann_gam_raw <- gam(per_frac / 100 ~ s(days_in_hyear, bs = "cc") + 
                      te(x, y, bs = "tp"),
                    data = landpks_per_frac,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_rpi <- gam(per_frac / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_per_frac,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_gpp <- gam(per_frac / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                      te(x, y, k = 6, bs = "tp"),
                    data = landpks_per_frac,
                    family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_pot_gpp <- gam(per_frac / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, bs = "cc") +
                         te(x, y, k = 6, bs = "tp"),
                       data = landpks_per_frac,
                       family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_rpi2 <- gam(per_frac / 100 ~ s(rpi, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                     data = landpks_model,
                     family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_gpp2 <- gam(per_frac / 100 ~ s(GPP, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                     data = landpks_model,
                     family = betar(link = "logit", eps = 0.001)
)

per_ann_gam_pot_gpp2 <- gam(per_frac / 100 ~ s(pot_gpp, k = 6, bs = "tp") + s(days_in_hyear, k = 12, bs = "cc"),
                         data = landpks_model,
                         family = betar(link = "logit", eps = 0.001)
)


# # Perennial fraction total
# 
# per_gam_raw <- gam(per_cover / 100 ~ s(days_in_hyear, bs = "cc") + 
#                          te(x, y, bs = "tp"),
#                        data = landpks_model,
#                        family = betar(link = "logit", eps = 0.001)
# )
# 
# per_gam_rpi <- gam(per_cover / 100 ~ rpi + s(days_in_hyear, bs = "cc") +
#                          te(x, y, bs = "tp"),
#                        data = landpks_model,
#                        family = betar(link = "logit", eps = 0.001)
# )
# 
# per_gam_gpp <- gam(per_cover / 100 ~ GPP + s(days_in_hyear, bs = "cc") +
#                          te(x, y, bs = "tp"),
#                        data = landpks_model,
#                        family = betar(link = "logit", eps = 0.001)
# )

## Plot function for a model

plot_gam <- function(model, exclude = "te(x,y)") {
  
  label_lookup <- list(
    "s(GPP)" = expression(GPP~(g~C~m^-2~yr^-1)),
    "s(pot_gpp)" = expression(Potential~GPP~(g~C~m^-2~yr^-1)),
    "s(rpi)" = "RPI",
    "s(days_in_hyear)" = "Day in hydrological year",
    "te(x,y)" = "Lat-lon"
  )
  
  vars <- smooth_label(model)
  
  vars_exclude <- vars[!(vars %in% exclude)]
  
  plots <- list()
  
  for (i in seq_along(vars_exclude)) {
    label <- label_lookup[[vars_exclude[i]]]
    
    p <- draw(
      model,
      ci_col = "steelblue",
      smooth_col = "blue",
      select = vars_exclude[i],
      caption = FALSE
    ) + labs(x = label, title = "") +
      coord_cartesian(ylim = c(-1.5, 1.5)) +
      theme_bw()
    
    plots[[i]] <- p
  }
  
  plots %>%
    patchwork::wrap_plots(axes = "collect")
}

bare_gam_plots <- list(bare_gam_gpp2, bare_gam_pot_gpp2, bare_gam_rpi2) %>%
  map(plot_gam)
herb_gam_plots <- list(herb_gam_gpp2, herb_gam_pot_gpp2, herb_gam_rpi2) %>%
  map(plot_gam)
per_ann_gam_plots <- list(per_ann_gam_gpp2, per_ann_gam_pot_gpp2, per_ann_gam_rpi2) %>%
  map(plot_gam)

bare_gam_fig <- ggpubr::ggarrange(plotlist = bare_gam_plots, nrow = 3)
herb_gam_fig <- ggpubr::ggarrange(plotlist = herb_gam_plots, nrow = 3)
per_ann_gam_fig <- ggpubr::ggarrange(plotlist = per_ann_gam_plots, nrow = 3)

ggsave("results/figures/bare_gam_fig.png", bare_gam_fig,
       height = 20, width = 16, units = "cm", dpi = 250)
ggsave("results/figures/herb_gam_fig.png", herb_gam_fig,
       height = 20, width = 16, units = "cm", dpi = 250)
ggsave("results/figures/per_ann_gam_fig.png", per_ann_gam_fig,
       height = 20, width = 16, units = "cm", dpi = 250)


### 7. Trend consistency

# Filter LandPKS data to those with multiple years sampled in the same location

landpks_multitemporal <- landpks_model %>%
  group_by(cell) %>%
  filter(length(unique(hYear)) > 1)

landpks_multitemporal_plot <- ggplot(landpks_multitemporal, aes(x = rpi, y = herb)) +
  geom_path(aes(colour = as.factor(cell)), show.legend = FALSE) +
  # geom_point(aes(fill = hYear)) +
  theme_bw()

landpks_trend <- landpks_multitemporal %>%
  group_by(cell) %>%
  mutate(centroid = find_centroid(geometry),
         dist_to_centroid = st_distance(geometry, centroid, by_element = TRUE)) %>%
  filter(dist_to_centroid <= set_units(100, "metres")) %>%
  filter(length(unique(hYear)) > 1) %>%
  summarise(n = length(unique(hYear)),
            start = first(hYear), end = last(hYear), length = end - start + 1,
            centroid = first(centroid),
            max_distance = max(dist_to_centroid),
            bare_trend = (lm(bare ~ hYear) %>% coef())[2],
            herb_trend = (lm(herb ~ hYear) %>% coef())[2],
            gpp_trend = (lm(GPP ~ hYear) %>% coef())[2],
            rpi_trend = (lm(rpi ~ hYear) %>% coef())[2],
            bare_trend_sign = ifelse(bare_trend > 0, "positive", "negative"),
            herb_trend_sign = ifelse(herb_trend > 0, "positive", "negative"),
            gpp_trend_sign = ifelse(gpp_trend > 0, "positive", "negative"),
            rpi_trend_sign = ifelse(rpi_trend > 0, "positive", "negative")) %>%
  filter(length >= 3)



# GPP vs perennial

ggplot(landpks_trend, aes(x = gpp_trend, y = herb_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits) +
  labs(x = "GPP trend", y = "Herbaceous cover trend")

# GPP vs bare ground
ggplot(landpks_trend, aes(x = gpp_trend, y = bare_trend)) +
  geom_point(aes(size = length, colour = n)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme_bw() + 
  scale_x_continuous(limits = symmetric_limits) + 
  scale_y_continuous(limits = symmetric_limits) +
  labs(x = "GPP trend", y = "Bare ground trend")

# RPI vs perennial
ggplot(landpks_trend, aes(x = rpi_trend, y = herb_trend)) +
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



# Path diagrams

# Group and join such that only plots within 100 m of each other in the same cell
# are retained (with the plot with the largest number of nearby plots in the
# same cell retained as the base plot for that cell)

# Filter to one row per sample and calculate number of samples in each cell
landpks_pre_clustered <- landpks_rpi_gpp %>%
  group_by(sample_id) %>%
  slice_head(n = 1) %>%
  select(Name, sample_id, cell, ObservationDate_GMT, hYear, days_in_hyear, geometry) %>%
  group_by(cell) %>%
  mutate(n_plots = n())

# Spatial join with a 100m buffered version of all points, keeping those within the
# same cell and those with at least two points.
landpks_clustered <- st_join(landpks_pre_clustered, landpks_pre_clustered %>% st_buffer(dist = set_units(100, "m")), join = st_within) %>%
  filter(cell.x == cell.y) %>%
  st_drop_geometry()

landpks_clustered_base <- landpks_clustered %>%
  group_by(sample_id.x, cell.x) %>%
  summarise(match_list = list(sample_id.y), n = length(unique(sample_id.y))) %>%
  filter(n > 1) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 1)

landpks_cell_groups <- landpks_clustered_base %>%
  group_by(cell.x) %>%
  unnest(match_list) %>%
  left_join(landpks_rpi_gpp %>% nest(sample_id, cell), by = c("cell.x" = "cell", "match_list" = "sample_id"))

ggplot(landpks_multitemporal, aes(x = herb, y = GPP, colour = as.factor(cell))) +
  geom_path(show.legend = FALSE) +
  theme_bw() +
  labs(x = "% Herbaceous", y = expression(GPP~(g~C~m^-2~yr^-1)))

