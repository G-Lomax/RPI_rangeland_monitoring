## Evaluation of the impact of Anthrax outbreaks in wild herbivores
## 8th January 2025
## Guy Lomax
## G.Lomax@exeter.ac.uk

## 1. Setup ----

source("scripts/load.R")

## 2. Load data ----

# RPI and covariate rasters

rpi_map <- rast("data/processed/raster/rpi_rast_sp.tif")

# Tree cover raster

tc <- rast("data/raw/raster/covariateMaps/staticVars.tif") %>%
  subset("treeCover") %>%
  crop(rpi_map)

# Country polygons
kenya <- st_read("data/raw/vector/kenya_tanzania.geojson") %>%
  filter(ADM0_NAME == "Kenya")

# Protected areas

pa <- st_read("data/raw/vector/WDPA/wdpa.shp")

## 3. Data cleaning and preparation -----

# Select Kenyan terrestrial national parks existing in 2010 and drop unneeded cols

np <- pa %>%
  filter(DESIG_ENG %in% c("National Park", "National Reserve") & STATUS_YR <= 2010) %>%
  select(NAME, DESIG_ENG, geometry) %>%
  filter(NAME != "Tsavo Road and Railways") %>%  # Not a national reserve
  filter(!(NAME == "Marsabit" & DESIG_ENG == "National Reserve"))  # Duplicate of national park

nakuru <- filter(np, NAME == "Lake Nakuru")

# Prepare rasters

rpi_gpp_only <- subset(rpi_map, str_detect(names(rpi_map), "rpi|GPP"))
rpi_only <- subset(rpi_map, str_detect(names(rpi_map), "rpi"))

rpi_nakuru <- crop(rpi_only, nakuru, mask = TRUE)

tc_mask <- mask(tc <= 60, tc <= 60, maskvalues = c(0, NA))

tc_nakuru <- (tc > 60) %>%
  mask(tc > 60, maskvalues = c(0, NA)) %>%
  crop(nakuru, mask = TRUE)

# Create masked rpi raster excluding high tree cover areas

rpi_gpp_low_tc <- mask(rpi_gpp_only, tc_mask)

# Filter to those > 25% covered by RPI dataset

np_rpi <- extract(!is.na(rpi_map$rpi.2000), np, fun = "mean", weights = TRUE, bind = TRUE) %>%
  st_as_sf() %>%
  filter(rpi.2000 >= 0.25) %>%
  rename(frac_cover = rpi.2000) %>%
  mutate(area = st_area(geometry) %>% set_units("km2") %>% drop_units()) %>%
  mutate(Designation = ifelse(NAME == "Lake Nakuru", "Lake Nakuru NP", DESIG_ENG))

# National parks remaining in Kenya

pa_map <- tm_shape(kenya) +
  tm_polygons(fill = "grey95") +
  tm_shape(rpi_map$rpi.2000 %>% mask(kenya)) +
  tm_raster(col.scale = tm_scale_continuous(values = "wheat"),
            col.legend = tm_legend(show = FALSE)) +
  tm_shape(kenya) +
  tm_borders(lwd = 2) +
  tm_shape(np_rpi) +
  tm_borders(col = "Designation",
             col.scale = tm_scale_categorical(values = c("red", "darkgreen", "blue4")),
             col.legend = tm_legend(frame = FALSE, position = tm_pos_out("right", "top")),
             lwd = 3) +
  tm_layout(frame = FALSE)

tmap_save(pa_map, "results/figures/np_map_kenya.png",
          height = 16, width = 16, units = "cm", dpi = 250)

## 4. Extract GPP and RPI time series for all parks ----

np_ts <- rpi_gpp_low_tc %>%
  extract(np_rpi, fun = "mean", weights = TRUE, bind = TRUE, na.rm = TRUE) %>%
  st_as_sf() %>%
  tidy_annual_vars()

nakuru_ts <- filter(np_ts, NAME == "Lake Nakuru")

gpp_ts <- ggplot(np_ts, aes(x = year, y = GPP, group = NAME, colour = Designation)) +
  geom_line(alpha = 0.2) +
  geom_line(data = nakuru_ts, lwd = 0.8, colour = "red") +
  theme_bw() +
  geom_vline(xintercept = 2015) +
  ylim(0, 6000) +
  labs(x = "Year", y = expression(GPP~(g~C~m^-2~y^-1)))

rpi_ts <- ggplot(np_ts, aes(x = year, y = rpi, group = NAME, colour = Designation)) +
  geom_line(alpha = 0.2) +
  geom_line(data = nakuru_ts, lwd = 0.8, colour = "red") +
  theme_bw() +
  geom_vline(xintercept = 2015) +
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(x = "Year", y = "RPI")

ts_plot <- rpi_ts + gpp_ts + plot_layout(nrow = 2, axes = "collect", guides = "collect") &
  theme(legend.position='bottom')

ggsave("results/figures/np_ts_plot.png", ts_plot,
       width = 20, height = 16, units = "cm", dpi = 250)

## 4B: Visualisation of RPI change after 2015

rpi_map_nakuru <- tm_shape(rpi_nakuru[[13:20]]) +
  tm_raster(col.scale = tm_scale_continuous(values = "viridis", limits = c(0.2, 1), outliers.trunc = c(TRUE, TRUE)),
            col.legend = tm_legend(reverse = TRUE, show = FALSE)) +
  tm_shape(tc_nakuru) +
  tm_raster(col.scale = tm_scale_discrete(values = "darkgreen"), col.legend = tm_legend(labels = "Tree cover")) +
  tm_shape(nakuru) +
  tm_borders(lwd = 3) +
  tm_facets(nrow = 2, ncol = 4) +
  tm_check_fix()

rpi_map_nakuru
  
## 5a: Difference in differences analysis

did_ts <- np_ts %>%
  filter(year >= 2013 & year <= 2018) %>%
  mutate(nakuru = (NAME == "Lake Nakuru"),
         after = (year > 2015)) %>%
  group_by(after, nakuru) %>%
  mutate(did = after * nakuru)

did <- lm(rpi ~ NAME + after + did, data = did_ts)


# Now controlling for pre-event trends



## 5: Synthetic control ----

np_ts_augsynth <- np_ts %>%
  # filter(DESIG_ENG == "National Park") %>%
  mutate(treated = ifelse(NAME == "Lake Nakuru" & year > 2015, 1, 0))

rpi_synth <- augsynth::augsynth(rpi ~ treated, NAME, year, np_ts_augsynth, t_int = 2016, progfunc = "Ridge")

gpp_synth <- augsynth::augsynth(GPP ~ treated, NAME, year, np_ts_augsynth, t_int = 2016, progfunc = "Ridge")

# Plot actual and predicted RPI and GPP with CI

rpi_att <- summary(rpi_synth)$att
gpp_att <- summary(gpp_synth)$att

synth_preds <- np_ts_augsynth %>%
  filter(NAME == "Lake Nakuru") %>%
  st_drop_geometry() %>%
  select(-geometry) %>%
  left_join(rpi_att, by = c("year" = "Time")) %>%
  left_join(gpp_att, by = c("year" = "Time"), suffix = c(".rpi", ".gpp")) %>%
  rename(truth.rpi = rpi, truth.gpp = GPP) %>%
  pivot_longer(ends_with(c(".gpp", ".rpi"))) %>%
  separate_wider_delim(cols = "name", delim = ".", names = c("param", "metric")) %>%
  pivot_wider(names_from = "param", values_from = "value") %>%
  mutate(across(.cols = c("Estimate", "lower_bound", "upper_bound"), .fns = \(x) x + truth)) %>%
  pivot_longer(cols = c("truth", "Estimate"))

rpi_preds <- filter(synth_preds, metric == "rpi")
gpp_preds <- filter(synth_preds, metric == "gpp")

# synth_plots <- ggplot(synth_preds, aes(x = year, y = value)) +
#   geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "steelblue", alpha = 0.5) +
#   geom_line(aes(colour = "Name"), colour = "black", lwd = 1) +
#   geom_line(aes(y = Estimate), colour = "steelblue", lwd = 1) +
#   geom_vline(xintercept = 2016) +
#   facet_wrap(~metric, nrow = 2, scales = "free") +
#   ggh4x::facetted_pos_scales(y = list(
#     metric == "gpp" ~ scale_y_continuous(limits = c(0, 8000)),
#     metric == "rpi" ~ scale_y_continuous(limits = c(0, 1.2))
#     )) +
#   labs(x = "Year", y = "Value") +
#   theme_bw()

rpi_synth_plot <- ggplot(rpi_preds, aes(x = year, y = value, colour = name)) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "steelblue", colour = "transparent", alpha = 0.2) +
  geom_line(lwd = 1) +
  geom_vline(xintercept = 2016) +
  scale_colour_manual(
    values = c("steelblue", "black"),
    labels = c("Synthetic Lake Nakuru NP", "Lake Nakuru NP")
  ) +
  ylim(0, 1.4) +
  labs(x = "Year", y = "RPI", colour = "") +
  theme_bw()

gpp_synth_plot <- ggplot(gpp_preds, aes(x = year, y = value, colour = name)) +
  geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "steelblue", colour = "transparent", alpha = 0.2) +
  geom_line(lwd = 1) +
  geom_vline(xintercept = 2016) +
  scale_colour_manual(
    values = c("steelblue", "black"),
    labels = c("Synthetic Lake Nakuru NP", "Lake Nakuru NP")
  ) +
  ylim(0, 8200) +
  labs(x = "Year", y = expression(GPP~(g~C~m^-2~y^-1)), colour = "") +
  theme_bw()

synth_plots <- rpi_synth_plot + gpp_synth_plot +
  plot_layout(nrow = 2, axes = "collect", axis_titles = "collect", guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results/figures/synth_plots.png", synth_plots,
       width = 20, height = 16, units = "cm", dpi = 250)

# Plot effect sizes from both indices

rpi_synth_plot <- plot(rpi_synth) + coord_cartesian(ylim = c(-0.4, 0.4)) + labs(x = "Year", y = "Estimated effect\nsize (RPI)")
gpp_synth_plot <- plot(gpp_synth) + coord_cartesian(ylim = c(-3000, 3000)) + labs(x = "Year", y = "Estimated effect\nsize (GPP)")

synth_plot_stack <- rpi_synth_plot + gpp_synth_plot + plot_layout(axes = "collect", nrow = 2)

ggsave("results/figures/synth_effects_all.png", synth_plot_stack,
       width = 20, height = 16, units = "cm", dpi = 250)

# Unit weights

rpi_weights <- rpi_synth$weights %>%
  as.data.frame() %>%
  rename(rpi = V1) %>%
  mutate(NAME = rownames(.)) %>%
  tibble::remove_rownames()

gpp_weights <- gpp_synth$weights %>%
  as.data.frame() %>%
  rename(gpp = V1) %>%
  mutate(NAME = rownames(.)) %>%
  tibble::remove_rownames()

all_weights <- left_join(rpi_weights, gpp_weights)

weights_plot <- all_weights %>%
  pivot_longer(cols = c("rpi", "gpp"), names_to = "index", values_to = "weight") %>%
  ggplot(aes(x = weight, y = NAME, fill = index)) +
  geom_col(position = "dodge") +
  theme_bw() +
  xlim(-0.6, 0.6) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("steelblue2", "green3"), labels = c("GPP", "RPI")) +
  labs(x = "Unit weight", y = "Unit name", fill = "Index")

ggsave("results/figures/sc_weights.png", weights_plot,
       width = 16, height = 20, units = "cm", dpi = 250)

# Performance metrics

synth_preds <- np_ts_augsynth %>%
  filter(NAME == "Lake Nakuru") %>%
  mutate(rpi_pred = predict(rpi_synth),
         gpp_pred = predict(gpp_synth),
         rpi_error = abs(rpi_pred - rpi),
         gpp_error = abs(gpp_pred - GPP),
         naive_step_rpi = rpi - lag(rpi),
         naive_step_gpp = GPP - lag(GPP))

scaled_errors <- synth_preds %>%
  filter(year <= 2015) %>%
  summarise(rpi_mase = mean(abs(rpi_error)) / mean(abs(naive_step_rpi), na.rm = T),
            gpp_mase = mean(abs(gpp_error)) / mean(abs(naive_step_gpp), na.rm = T),
            rpi_rmsse = sqrt(mean(rpi_error ^ 2) / mean(naive_step_rpi ^ 2, na.rm = T)),
            gpp_rmsse = sqrt(mean(gpp_error ^ 2) / mean(naive_step_gpp ^ 2, na.rm = T)),
            rpi_mae_z = mean(abs(rpi_error)) / sd(rpi),
            gpp_mae_z = mean(abs(gpp_error)) / sd(GPP),
            rpi_mase2 = mean(abs(rpi_error)) / mean(abs(rpi - mean(rpi))),
            gpp_mase2 = mean(abs(gpp_error)) / mean(abs(GPP - mean(GPP))),
            rpi_mae_pc = mean(abs(rpi_error)) / mean(rpi),
            gpp_mae_pc = mean(abs(gpp_error)) / mean(GPP),
            rpi_smape = 2 * mean(abs(rpi_error) / abs(rpi + rpi_pred)),
            gpp_smape = 2 * mean(abs(gpp_error) / abs(GPP + gpp_pred))
  )

# Plot trajectories

ggplot(synth_preds, aes(x = year)) +
  geom_line(aes(y = rpi)) +
  geom_line(aes(y = rpi_pred), colour = "blue") +
  theme_bw()



















## MSCMT implementation

np_ts_mscmt <- np_ts_augsynth %>%
  select(-geometry, -DESIG_ENG, -Designation) %>%
  group_by(NAME) %>%
  mutate(id = as.numeric(cur_group_id()),
         year2 = year,
         rpi2 = rpi,
         gpp2 = GPP) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("year2", "rpi2", "gpp2", "NAME", "id"), names_from = year, values_from = c(GPP, rpi)) %>%
  as.data.frame() %>%
  listFromLong(unit.variable = "id",  time.variable = "year2", unit.names.variable = "NAME")

controls_id <- setdiff(unique(np_ts_augsynth$NAME), "Lake Nakuru")
times_rpi <- cbind("rpi2" = c(2000, 2022))
times_gpp <- cbind("gpp2" = c(2000, 2022))
times_year_rpi <- cbind(
  "rpi_2000" = c(2000, 2000),
  "rpi_2001" = c(2001,2001),
  "rpi_2002" = c(2002,2002),
  "rpi_2003" = c(2003,2003),
  "rpi_2004" = c(2004,2004),
  "rpi_2005" = c(2005,2005),
  "rpi_2006" = c(2006,2006),
  "rpi_2007" = c(2007,2007),
  "rpi_2008" = c(2008,2008),
  "rpi_2009" = c(2009,2009),
  "rpi_2010" = c(2010,2010),
  "rpi_2011" = c(2011,2011),
  "rpi_2012" = c(2012,2012),
  "rpi_2013" = c(2013, 2013),
  "rpi_2014" = c(2014, 2014),
  "rpi_2015" = c(2015, 2015)
)

times_year_gpp <- cbind(
  "GPP_2000" = c(2000, 2000),
  "GPP_2001" = c(2001,2001),
  "GPP_2002" = c(2002,2002),
  "GPP_2003" = c(2003,2003),
  "GPP_2004" = c(2004,2004),
  "GPP_2005" = c(2005,2005),
  "GPP_2006" = c(2006,2006),
  "GPP_2007" = c(2007,2007),
  "GPP_2008" = c(2008,2008),
  "GPP_2009" = c(2009,2009),
  "GPP_2010" = c(2010,2010),
  "GPP_2011" = c(2011,2011),
  "GPP_2012" = c(2012,2012),
  "GPP_2013" = c(2013, 2013),
  "GPP_2014" = c(2014, 2014),
  "GPP_2015" = c(2015, 2015)
)


rpi_mscmt <- mscmt(
  np_ts_mscmt,
  treatment.identifier = "Lake Nakuru",
  controls.identifier = setdiff(unique(np_ts_augsynth$NAME), "Lake Nakuru"),
  times.dep = times_rpi,
  times.pred = times_year_rpi
)

gpp_mscmt <- mscmt(
  np_ts_mscmt,
  treatment.identifier = "Lake Nakuru",
  controls.identifier = setdiff(unique(np_ts_augsynth$NAME), "Lake Nakuru"),
  times.dep = times_gpp,
  times.pred = times_year_gpp
)


data_synth <-  np_ts_augsynth %>%
  select(-geometry) %>%
  group_by(NAME) %>%
  mutate(id = as.numeric(cur_group_id())) %>%
  ungroup() %>%
  as.data.frame() %>%
  dataprep(
    dependent = "rpi2",
    unit.variable = "id", 
    time.variable = "year2",
    unit.names.variable = "NAME",
    treatment.identifier = "Lake Nakuru",
    controls.identifier = controls_id,
    time.predictors.prior = 2000:2015,
    time.optimize.ssr = 2000:2015)
