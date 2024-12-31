## Evaluation against Tanzania People and Wildlife field monitoring dataset
## 1st October 2024
## Guy Lomax
## G.Lomax@exeter.ac.uk

# TPW data downloaded on 1st October 2024

## 1. Setup ----

source("scripts/load.R")

START <- as.Date("2018-09-01")
END <- as.Date("2023-08-31")

## 2. Load data ----

# TPW field sampling data
tpw_locations <- st_read("data/raw/vector/tpw/tpw_locations.shp")
tpw_data <- foreign::read.dbf("data/raw/vector/tpw/tpw_samplemetrics.dbf")

inv_species_lookup <- read_csv("data/raw/vector/tpw/invasive_plant_lookup.csv")

# RPI and covariate rasters

rpi_map <- rast("data/processed/raster/rpi_rast_sp.tif")

static_vars <- rast("data/raw/raster/covariateMaps/staticVars.tif")

## 3. Pre-process data and extract raster values ----

# Convert transect strikes into site and date-level % cover and count values
tpw_data_summary <- x<- tpw_data %>%
  mutate(across(all_of(inv_species_lookup$colname), ~ as.numeric(as.character(.x)))) %>%
  group_by(village_pl, parentglob, date) %>%
  filter(markbare >= 0 & markbasal >= 0) %>%
  filter(n() == 20) %>%
  filter(date >= START & date <= END) %>%
  summarise(bare = sum(markbare),
            basal = sum(markbasal),
            invasive = sum(sum(!!!syms(inv_species_lookup$colname))))

# Join with locations
tpw_data_joined <- inner_join(tpw_locations, tpw_data_summary, by = c("globalid" = "parentglob")) %>%
  filter(village_pl.x == village_pl.y)  # Remove one mismatched data entry

# Remove superfluous columns, add hydrological year and day in hyear values

tpw_data_hyear <- tpw_data_joined %>%
  select(village, village_pl.x, date, plotburned, plotgrazed, bare, basal, invasive) %>%
  mutate(hYear = interval(START, date) %/% years() + year(START),
         day_in_hYear = interval(START, date) %% years() %/% days(),
         month = month(date))


# # Aggregate by hYear
# tpw_data_hyear_agg <- tpw_data_hyear %>%
#   group_by(village_pl.x, hYear) %>%
#   summarise(
#     across(
#       .cols = c("bare", "basal", "invasive"),
#       .fns = list("mean" = mean, "max" = max, "min" = min)
#     ), 
#     n = n()
#   )

# Extract RPI raster values

data_years <- unique(tpw_data_hyear$hYear) %>% sort()

tpw_rpi_gpp <- map(data_years, raster_extract, sf = tpw_data_hyear, raster = rpi_map) %>%
  bind_rows() %>%
  drop_na()

# Merge values that lie in the same raster cell

# tpw_rpi_gpp_merged <- tpw_rpi_gpp %>%
#   group_by(cell, hYear) %>%
#   summarise(
#     across(
#       starts_with(c("bare", "basal", "invasive")),
#       mean
#     ),
#     rpi = mean(rpi),
#     GPP = mean(GPP),
#     n = n()
#   )

# ## 4. Plot key indicators against GPP and RPI
# 
# tpw_gpp_tpw_rpi_gpp_merged %>%
#   st_drop_geometry() %>%
#   pivot_longer(ends_with(c("mean", "max", "min"))) %>%
#   separate_wider_delim(name, delim = "_", names = c("var", "fun")) %>%
#   ggplot(aes(x = value, y = GPP)) +
#   geom_point(size = 0.5) +
#   facet_grid(cols = vars(fun), rows = vars(var)) +
#   theme_bw() +
#   labs(x = "Cover")
# 
# tpw_rpi_gpp_merged %>%
#   st_drop_geometry() %>%
#   pivot_longer(ends_with(c("mean", "max", "min"))) %>%
#   separate_wider_delim(name, delim = "_", names = c("var", "fun")) %>%
#   ggplot(aes(x = value, y = rpi)) +
#   geom_point(size = 0.5) +
#   facet_grid(cols = vars(fun), rows = vars(var)) +
#   theme_bw() +
#   labs(x = "Cover", y = "rpi")


# ## 5. Simple GAMs for bare ground considering seasonal cycle
# 
# tpw_model <- tpw_rpi_gpp %>%
#   mutate(x = st_coordinates(geometry)[,1],
#          y = st_coordinates(geometry)[,2]) %>%
#   st_drop_geometry()
# 
# bare_gam_raw <- bam(bare / 100 ~ s(day_in_hYear, bs = "cc", k = 12) +
#                       te(x, y, bs = "ts", k = 8),
#                     data = tpw_model,
#                     method = "REML",
#                     family = betar(eps = 0.001)
# )
# 
# bare_gam_rpi <- bam(bare / 100 ~ s(rpi, bs = "ts", k = 8) +
#                       s(day_in_hYear, bs = "cc", k = 12) +
#                       te(x, y, bs = "ts", k = 8),
#                     data = tpw_model,
#                     method = "REML",
#                     family = betar(eps = 0.001)
# )
# 
# bare_gam_gpp <- bam(bare / 100 ~ s(GPP, bs = "ts", k = 8) +
#                       s(day_in_hYear, bs = "cc", k = 12) +
#                       te(x, y, bs = "ts", k = 8),
#                     data = tpw_model,
#                     method = "REML",
#                     family = betar(eps = 0.001)
# )

## 6. Condense TPW data to min wet-season values

wet_season_start <- yday("2001-11-01")
wet_season_end <- yday("2002-05-01")

tpw_wet_season <- tpw_rpi_gpp %>%
  filter(yday(date) >= wet_season_start | yday(date) < wet_season_end) %>%
  mutate(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2]) %>%
  st_drop_geometry() %>%
  group_by(cell, village_pl.x, hYear) %>%
  summarise(rpi = mean(rpi), GPP = mean(GPP),
            x = median(x), y = median(y),
            mean_wet_bare = mean(bare),
            min_wet_bare = min(bare),
            min_bare_month = month[which.min(bare)]) %>%
  group_by(cell, hYear) %>%
  summarise(rpi = mean(rpi), GPP = mean(GPP),
            x = mean(x), y = mean(y),
            mean_wet_bare = mean(mean_wet_bare),
            min_wet_bare = mean(min_wet_bare))


mean_wet_gam_rpi <- gam(mean_wet_bare / 100 ~ s(rpi, bs = "ts", k = 8) + te(x,y, bs = "ts", k = 8),
                        data = tpw_wet_season,
                        family = betar(eps = 0.001),
                        method = "REML"
)

mean_wet_gam_gpp <- gam(mean_wet_bare / 100 ~ s(GPP, bs = "ts", k = 8) + te(x,y, bs = "ts", k = 8),
                        data = tpw_wet_season,
                        family = betar(eps = 0.001),
                        method = "REML"
)

min_wet_gam_rpi <- gam(min_wet_bare / 100 ~ s(rpi, bs = "ts", k = 8) + te(x,y, bs = "ts", k = 8),
                        data = tpw_wet_season,
                        family = betar(eps = 0.001),
                       method = "REML"
)

min_wet_gam_gpp <- gam(min_wet_bare / 100 ~ s(GPP, bs = "ts", k = 8) + te(x,y, bs = "ts", k = 8),
                        data = tpw_wet_season,
                        family = betar(eps = 0.001),
                       method = "REML"
)


mean_wet_gam_rpi2 <- gam(mean_wet_bare / 100 ~ s(rpi, bs = "ts", k = 8) + te(x,y, bs = "ts", k = 8),
                        data = tpw_wet_season,
                        family = "gaussian",
                        method = "REML"
)

## Test BRMS models

tpw_wet_season_one_adjusted <- tpw_wet_season %>%
  mutate(across(
    .cols = contains("bare"),
    .fns = ~ ifelse(.x == 100, 99.5, .x)
  ))

mean_wet_gam_rpi_brms <- brm(
  bf(
    mean_wet_bare / 100 ~ rpi + t2(x, y, bs = "ts"),
    phi ~ rpi,
    zi ~ rpi
    ),
  data = tpw_wet_season_one_adjusted,
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4, seed = 111,
  backend = "rstan",
  control = list(adapt_delta = 0.95),
  file = "data/processed/rds/wet_gam_rpi_brms"
)

min_wet_gam_gpp_brms <- brm(
  bf(
    min_wet_bare / 100 ~ rpi + t2(x, y, bs = "ts"),
    phi ~ rpi,
    zi ~ rpi
  ),
  data = tpw_wet_season_one_adjusted,
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4, seed = 111,
  backend = "rstan",
  control = list(adapt_delta = 0.95),
  file = "data/processed/rds/min_wet_gam_rpi_brms"
)

# Scale GPP values to enable MCMC model fit

tpw_wet_season_scaled <- tpw_wet_season_one_adjusted %>%
  mutate(GPP = GPP / 1000)

mean_wet_gam_gpp_brms <- brm(
  bf(
    mean_wet_bare / 100 ~ GPP + t2(x, y, bs = "ts"),
    phi ~ GPP,
    zi ~ GPP
  ),
  data = tpw_wet_season_scaled,
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4, seed = 111,
  backend = "rstan",
  control = list(adapt_delta = 0.95),
  file = "data/processed/rds/wet_gam_gpp_brms"
)

min_wet_gam_gpp_brms <- brm(
  bf(
    min_wet_bare / 100 ~ GPP + t2(x, y, bs = "ts"),
    phi ~ GPP,
    zi ~ GPP
  ),
  data = tpw_wet_season_scaled,
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4, seed = 111,
  backend = "rstan",
  control = list(adapt_delta = 0.95),
  file = "data/processed/rds/min_wet_gam_gpp_brms"
)
