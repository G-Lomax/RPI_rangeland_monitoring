## This script fits and tunes a quantile regression forest model to rangeland
## productivity data for calculation of the Relative Productivity Index (RPI)

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 18th September 2024

## 1. Setup ----

source("scripts/load.R")

# Parallelisation
nc <- availableCores() - 8
plan(multicore, workers = nc)

# Parameters

YEARS <- 2000:2022
MIN_STREAM_ORDER <- 4
N_SAMPLE <- 10000
SEED <- 451

## 2. Load data ----

static_covariates <- rast("data/raw/raster/covariate_maps/staticVars.tif")

if(!file.exists("data/processed/raster/dynamic_covariates.tif")) {
  dynamic_covariates <- read_annual_rasters(YEARS, "data/raw/raster/covariate_maps") %>%
    mask_zeros(layer = paste0("GPP.", YEARS[1])) %>%
    terra::resample(
      static_covariates, 
      threads = TRUE, 
      filename = "data/processed/raster/dynamic_covariates.tif",
      overwrite = TRUE)  # Need to resample because Google Earth Engine...
} else {
  dynamic_covariates <- rast("data/processed/raster/dynamic_covariates.tif")
}

dist_to_river <- rast(paste0("data/processed/raster/dist_to_river/dist_to_river_", MIN_STREAM_ORDER, ".tif"))

names(dist_to_river) <- "dist_to_river"

twi <- rast("data/processed/raster/merit/merit_twi_fd8.tif")

all_vars <- c(static_covariates, dist_to_river, twi, dynamic_covariates)

## 3. Generate random sample of training points ----
set.seed(SEED)

var_sample <- spatSample(
  x = all_vars,
  size = N_SAMPLE,
  na.rm = TRUE,
  cells = TRUE,
  xy = TRUE,
  exp = 10)

var_sample_long <- var_sample %>%
  pivot_longer(cols = contains(".")) %>%
  separate_wider_delim(name, delim = ".", names = c("var", "year")) %>%
  mutate(year = as.numeric(year)) %>%
  pivot_wider(names_from = "var", values_from = "value")

message("Grid cell sample complete")
pushoverr::pushover("Grid cell sample complete")

## 4. Additional variable processing ----
# Add PPT variables
# Designate ECO_ID as factor
var_sample_final <- var_sample_long %>%
  group_by(cell) %>%
  mutate(
    pptMean = mean(precipitation),
    pptAnomaly = (precipitation - pptMean) / pptMean * 100,
    pptMeanDayAnomaly = pptMeanDay - mean(pptMeanDay),
    ECO_ID = factor(ECO_ID)
  ) %>%
  ungroup()

write_csv(var_sample_final, "data/processed/csv/var_sample_final.csv")

## 5. Set up quantile regression model using mlr3 and ranger ----

# Create task

vars_to_retain <- c(
  # ID vars
  "cell",
  "year",
  "x", "y",
  
  # Annual vars
  "GPP",
  "pptAnomaly",
  "pptIntensity",
  "pptMeanDayAnomaly",
  "ugi",
  "temperature_2m_mean",
  "potential_evaporation_sum",
  "GMT_0900_PAR",
  
  # Static vars
  "pptMean",
  "treeCover",
  "ECO_ID",
  "DEM",
  "slope",
  "merit_twi_fd8",
  "dist_to_river",
  "sand"
)

var_sample_subset <- select(var_sample_final, all_of(vars_to_retain))

var_sample_sf <- var_sample_subset %>%
  st_as_sf(crs = "EPSG:4326", coords = c("x", "y")) %>%
  st_transform("ESRI:54034")

task_gpp <- as_task_regr_st(
  x = var_sample_sf,
  id = "potential_gpp",
  target = "GPP",
  coords_as_features = FALSE
)

task_gpp$set_col_roles("cell", roles = "space")
task_gpp$set_col_roles("year", roles = "time")

# Create learner (ranger random forest) with initial tuning values

lrn_ranger_untuned <- lrn("regr.ranger", 
                          predict_type = "response",
                          num.trees = 1001,
                          mtry.ratio = 0.33,
                          min.node.size = 100,
                          sample.fraction = 0.5,
                          respect.unordered.factors = "order"
)

# Create resampling strategy

spt_cv_plan <- rsmp("sptcv_cstf", folds = 20)
sp_cv_plan <- rsmp("spcv_coords", folds = 20)
# disc_cv_plan <- rsmp("spcv_disc", folds = 20, radius = 100000, buffer = 25000)

## 6. Feature selection with initial tuning values ----

# Create forward feature selection with untuned model, minimising RMSE

perf_msr <- msr("regr.rmse")

fs_method <- fs("sequential", min_features = 1)

fs_term <- trm("stagnation_batch")

spt_feature_select <- fsi(
  task = task_gpp,
  learner = lrn_ranger_untuned,
  resampling = spt_cv_plan,
  measures = perf_msr,
  terminator = fs_term
)

sp_feature_select <- fsi(
  task = task_gpp,
  learner = lrn_ranger_untuned,
  resampling = sp_cv_plan,
  measure = perf_msr,
  terminator = fs_term
)

# Identify optimal feature set for each resampling strategy and store
# Time ~ 10-12 hours

set.seed(123)

tic()
progressr::with_progress(
  spt_feature_set <- fs_method$optimize(spt_feature_select)
)
toc()


write_rds(spt_feature_select, "data/processed/rds/feature_selector_spt.rds")
write_rds(spt_feature_set, "data/processed/rds/features_spt.rds")
rm(spt_feature_select, spt_feature_set)
gc()
pushoverr::pushover("Feature selection complete: Spatiotemporal CV")


set.seed(456)
tic()
progressr::with_progress(
  sp_feature_set <- fs_method$optimize(sp_feature_select)
)
toc()

write_rds(sp_feature_select, "data/processed/rds/feature_selector_sp.rds")
write_rds(sp_feature_set, "data/processed/rds/features_sp.rds")
rm(sp_feature_select, sp_feature_set)
gc()

pushoverr::pushover("Feature selection complete: Spatial block CV")


## 7. Hyperparameter tuning with final feature set ----

# Load features to retain and assign to new tasks
spt_feature_set <- read_rds("data/processed/rds/features_spt.rds")
sp_feature_set <- read_rds("data/processed/rds/features_sp.rds")

task_gpp_spt <- task_gpp$clone()
task_gpp_spt$select(unlist(spt_feature_set$features))

task_gpp_sp <- task_gpp$clone()
task_gpp_sp$select(unlist(sp_feature_set$features))

# Define new quantile regression learners for (auto)tuning

lrn_ranger_tuned_spt <- lrn(
  "regr.ranger", 
  predict_type = "response",
  quantreg = TRUE,
  keep.inbag = TRUE,
  importance = "permutation",
  num.trees = 1001,
  mtry.ratio = to_tune(p_dbl(0, 1)),
  min.node.size = to_tune(p_int(100, 1000)),
  sample.fraction = to_tune(p_dbl(0.1, 0.9))
)

lrn_ranger_tuned_sp <- lrn(
  "regr.ranger", 
  predict_type = "response",
  quantreg = TRUE,
  keep.inbag = TRUE,
  importance = "permutation",
  num.trees = 1001,
  mtry.ratio = to_tune(p_dbl(0, 1)),
  min.node.size = to_tune(p_int(100, 1000)),
  sample.fraction = to_tune(p_dbl(0.1, 0.9))
)

# Define auto-tuner objects for each model

random_tuner <- tnr("random_search")

perf_msr <- msr("regr.rmse")

at_spt <- auto_tuner(
  tuner = random_tuner,
  learner = lrn_ranger_tuned_spt,
  resampling = spt_cv_plan,
  measure = perf_msr,
  term_evals = 25
)

at_sp <- auto_tuner(
  tuner = random_tuner,
  learner = lrn_ranger_tuned_sp,
  resampling = sp_cv_plan,
  measure = perf_msr,
  term_evals = 25
)

# Tune and train models

# Spatiotemporal model

set.seed(789)

tic()
rf_tuned_spt <- at_spt$train(task_gpp_spt)
toc()

write_rds(rf_tuned_spt, "data/processed/rds/rf_tuned_spt.rds")

rm(rf_tuned_spt)
gc()

pushoverr::pushover("Tuning complete: Spatiotemporal CV model")

# Spatial model

set.seed(987)

tic()
rf_tuned_sp <-at_sp$train(task_gpp_sp)
toc()

write_rds(rf_tuned_sp, "data/processed/rds/rf_tuned_sp.rds")

rm(rf_tuned_sp)
gc()

pushoverr::pushover("Tuning complete: Spatial block CV")
