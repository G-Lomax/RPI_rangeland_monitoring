## This script predicts the quantile regression forest model(s) to the
## covariate stack in order to estimate annual RPI for each grid cell
## in the study area.

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 25th September 2024

## 1. Set up and load data ----

source("scripts/load.R")

YEARS <- 2000:2022
MIN_STREAM_ORDER <- 4
THREADS <- 16

# Load raster data
static_covariates <- rast("data/raw/raster/covariate_maps/staticVars.tif")
dynamic_covariates <- rast("data/processed/raster/dynamic_covariates.tif")
dist_to_river <- rast(paste0("data/processed/raster/dist_to_river/dist_to_river_", MIN_STREAM_ORDER, ".tif"))
twi <- rast("data/processed/raster/merit/merit_twi_fd8.tif")

names(dist_to_river) <- "dist_to_river"

# Load models

qrf_model_spt <- read_rds("data/processed/rds/rf_tuned_spt.rds")
qrf_model_sp <- read_rds("data/processed/rds/rf_tuned_sp.rds")

## 2. Prepare additional precipitation covariates ----

pptAnnual <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "precipitation"))

pptMean <- mean(pptAnnual)
names(pptMean) <- "pptMean"

pptAnomaly <- (pptAnnual - pptMean) / pptMean * 100
names(pptAnomaly) <- str_replace(names(pptAnnual), "precipitation", "pptAnomaly")

pptMeanDayAnnual <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "pptMeanDay"))

pptMeanDayAnomaly <- pptMeanDayAnnual - mean(pptMeanDayAnnual)
names(pptMeanDayAnomaly) <- str_replace(names(pptAnnual), "precipitation", "pptMeanDayAnomaly")

# Define final sets of static and dynamic covariates
dynamic_covariates_all <- dynamic_covariates %>%
  subset(!str_detect(names(dynamic_covariates), c("precipitation|pptMeanDay"))) %>%
  c(pptAnomaly, pptMeanDayAnomaly)

static_covariates_all <- c(static_covariates, pptMean, dist_to_river, twi)

## 3. Predict iteratively to raster layers
tic()
annual_predictions_sp <- YEARS %>%
  map(predict_to_annual_raster,
      static_vars = static_covariates_all,
      dynamic_vars = dynamic_covariates_all,
      model = qrf_model_sp$learner$model,
      threads = THREADS) %>%
  rast()

toc()

writeRaster(annual_predictions_sp, "data/processed/raster/rpi_rast_sp.tif")

pushoverr::pushover("RPI raster written to disk")