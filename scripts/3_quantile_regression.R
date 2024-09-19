## This script fits and tunes a quantile regression forest model to rangeland
## productivity data for calculation of the Relative Productivity Index (RPI)

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 18th September 2024

## 1. Setup ----

source("scripts/load.R")

# Parallelisation
nc <- availableCores() / 2
plan(multicore, workers = nc)

# Parameters

YEARS <- 2000:2022
MIN_STREAM_ORDER <- 4
N_SAMPLE <- 10000
SEED <- 451

## 2. Load data ----

static_covariates <- rast("data/raw/raster/covariate_maps/staticVars.tif")

dynamic_covariates <- read_annual_rasters(YEARS, "data/raw/raster/covariate_maps") %>%
  mask_zeros(layer = paste0("GPP.", YEARS[1])) %>%
  terra::resample(
    static_covariates, 
    threads = TRUE, 
    filename = "data/processed/raster/dynamic_covariates.tif",
    overwrite = TRUE)  # Need to resample because Google Earth Engine...

dist_to_river <- rast(paste0("data/processed/raster/dist_to_river/dist_to_river_", MIN_STREAM_ORDER, ".tif"))

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
  "temperature_2m",
  # "potential_evaporation_sum",
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

task_gpp <- as_task_regr_st(
  x = var_sample_subset,
  id = "potential_gpp",
  target = "GPP",
  coords_as_features = TRUE,
  coordinate_names = c("x", "y"),
  crs = "epsg:4326"
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

spt_cv_plan <- rsmp("sptcv_cstf", folds = 10)
sp_cv_plan <- rsmp("spcv_coords", folds = 10)

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


write_rds(spt_feature_select, "data/processed/rds/spt_feature_selector.rds")
write_rds(spt_feature_set, "data/processed/rds/spt_features.rds")
rm(spt_feature_select, spt_feature_set)
gc()
pushoverr::pushover("Feature selection complete: Spatiotemporal CV")


set.seed(456)
tic()
progressr::with_progress(
  sp_feature_set <- fs_method$optimize(sp_feature_select)
)
toc()

write_rds(sp_feature_select, "data/processed/rds/sp_feature_selector.rds")
write_rds(sp_feature_set, "data/processed/rds/sp_features.rds")
rm(sp_feature_select, sp_feature_set)
gc()

pushoverr::pushover("Feature selection complete: Spatial block CV")

