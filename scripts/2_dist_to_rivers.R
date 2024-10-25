## Distance from nearest river of specified stream order
## 18th September 2024
## Guy Lomax
## G.Lomax@exeter.ac.uk

source("scripts/load.R")

# Input parameters
MIN_ORDER <- 4
BUFFER_DIST <- 5 # km

# Load data

covariates <- rast("data/raw/raster/covariate_maps/staticVars.tif")

river_files <- Sys.glob("data/raw/vector/Lin_rivers/*.shp")
rivers <- lapply(river_files, st_read) %>%
  bind_rows()

# Filter to study area and to rivers greater than specified Strahler stream order

study_area_buffered <- covariates %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(dist = set_units(BUFFER_DIST, "km"))

rivers_filtered <- rivers %>%
  filter(strmOrder >= MIN_ORDER) %>%
  st_filter(study_area_buffered)

# Calculate distance to nearest river for each grid cell in study region

rivers_rast <- rasterize(rivers_filtered, covariates, field = 1)

rivers_rast_reproject <- project(rivers_rast, "ESRI:54034")

dist_to_river <- distance(rivers_rast_reproject)

dist_to_river_sa <- dist_to_river %>%
  project(covariates$DEM) %>%
  mask(covariates$DEM)

names(dist_to_river_sa) <- "dist_to_river"

writeRaster(dist_to_river_sa, paste0("data/processed/raster/dist_to_river_", MIN_ORDER, ".tif"),
            overwrite = TRUE)
