## Topographic Wetness Index from DEM
## 18th September 2024
## Guy Lomax
## G.Lomax@exeter.ac.uk

library(terra)
library(raster)
library(sf)
library(tmap)
library(whitebox)

tmap_mode("view")

wbt_init()

# Load data and join

merit_dem <- rast("data/raw/raster/merit/meritDem.tif")
covariates <- rast("data/raw/raster/covariateMaps/staticVars.tif")

# Reproject DEM to equal area projection

merit_dem_reproj <- project(merit_dem, "ESRI:54034")
writeRaster(merit_dem_reproj, "data/processed/raster/merit/meritDemReproj.tif")

# Save slope as separate layer for Whitebox processing

writeRaster(covariates$slope, "data/processed/raster/slope.tif")

# Calculate flow accumulation (Freeman D8) algorithm and mask to study area
wbt_fd8_flow_accumulation(
  dem = "data/processed/raster/merit/meritDemReproj.tif",
  output = "data/processed/raster/merit/merit_fa_fd8.tif",
  out_type = "specific contributing area"
)

# Mask to study area and save
flow_accum <- rast("data/processed/raster/merit/merit_fa_fd8.tif")
flow_accum_sa <- flow_accum %>%
  project(covariates$DEM) %>%
  crop(covariates$DEM, mask = TRUE)

writeRaster(flow_accum_sa, "data/processed/raster/merit/merit_fa_fd8_sa.tif")

# Calculate Topographic Wetness Index
wbt_wetness_index(
  sca = "data/processed/raster/merit/merit_fa_fd8_sa.tif",
  slope = "data/processed/raster/slope.tif",
  output = "data/processed/raster/merit/merit_twi_fd8.tif"
)

twi_fd8 <- rast("data/processed/raster/merit/merit_twi_fd8.tif")

plot(twi_fd8, xlim = c(30, 30.5), ylim = c(-5, -4.5))

# Visualise
tmap_options(raster.max.cells = 8e6)
tm_shape(twi_fd8) +
  tm_raster(
    col.scale = tm_scale_continuous(
      limits = c(0, 25),
      outliers.trunc = c(TRUE, TRUE),
      values = "viridis")
  )
