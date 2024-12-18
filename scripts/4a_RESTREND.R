## This script calculates RESTREND residuals for the updated GPP, precipitation
## and temperature datasets and outputs the results as a .tif file.

## Guy Lomax
## G.Lomax@exeter.ac.uk
## 18th December 2024

## 1. Setup ----

source("scripts/load.R")

nc <- availableCores() / 4

plan(multicore, workers = nc)

rpi_rast <- rast("data/processed/raster/rpi_rast_sp.tif")

gpp <- subset(rpi_rast, str_detect(names(rpi_rast), "GPP"))

dynamic_covariates <- rast("data/processed/raster/dynamic_covariates.tif") %>%
  crop(rpi_rast[[1]], mask = TRUE)

precipitation <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "precipitation"))
tMean <- subset(dynamic_covariates, str_detect(names(dynamic_covariates), "temperature_2m"))

##### RESTREND function

fit_restrend <- function(gpp_rast, ppt_rast, tMean_rast) {
  
  # Convert rasters to a single data frame
  message("Converting to data frames\n")
  
  gpp_values <- as.data.frame(gpp_rast, cells = TRUE, xy = TRUE, na.rm = T)
  ppt_values <- as.data.frame(ppt_rast, cells = TRUE, xy = FALSE, na.rm = T)
  tMean_values <- as.data.frame(tMean_rast, cells = TRUE, xy = FALSE, na.rm = T)
  
  data_df <- gpp_values %>%
    left_join(ppt_values, by = "cell") %>%
    left_join(tMean_values, by = "cell") %>%
    pivot_longer(cols = starts_with(c("GPP", "precipitation", "temperature_2m_mean"))) %>%
    separate_wider_delim(cols = "name", delim = ".", names = c("var", "year")) %>%
    pivot_wider(id_cols = c(cell, x, y, year), names_from = var, values_from = value) %>%
    drop_na()
  
  # Apply RESTREND function to each cell in data frame
  message("Building models and calculating residuals\n")
  
  apply_restrend <- function(df) {
    
    years <- sort(unique(df$year))
    
    model <- lm(GPP ~ precipitation + temperature_2m_mean, data = df)
    
    yint <- coef(model)[1]
    ppt_slope <- coef(model)[2]
    t_slope <- coef(model)[3]
    
    r_sq <- broom::glance(model)$r.squared
    ppt_p_value <- broom::tidy(model)$p.value[2]
    t_p_value <- broom::tidy(model)$p.value[3]
    
    resids <- residuals(model)
    
    output <- c(yint, ppt_slope, t_slope, r_sq, ppt_p_value, t_p_value, resids)
    
    if (length(resids) == length(years)) {
      names(output) <- c("yint", "ppt_slope", "t_slope",
                         "rsq", "ppt_p_value", "t_p_value", 
                         paste0("resid_", years))
    } else {
      output <- rep(NA, length(years))
      names(output) <- c("yint", "ppt_slope", "t_slope",
                         "rsq", "ppt_p_value", "t_p_value", 
                         paste0("resid_", years))
    }
    
    output
  }
  
  restrend_df <- data_df %>%
    group_by(cell, x, y) %>%
    nest() %>%
    ungroup() %>%
    mutate(restrend_list = future_map(data, apply_restrend)) %>%
    select(x,y,restrend_list) %>%
    unnest_wider(restrend_list)
  
  # Convert to raster
  message("Converting back to raster")
  restrend_rast <- rast(restrend_df, crs = crs(gpp_rast), extent = ext(gpp_rast))
  
  restrend_rast
}

tic()
restrend_rast <- fit_restrend(gpp_rast = gpp, ppt_rast = precipitation, tMean_rast = tMean)

writeRaster(restrend_rast, "data/processed/raster/restrend.tif", overwrite = TRUE)
x <- toc()

pushoverr::pushover(paste0("RESTREND fit complete. ", x$callback_msg))
