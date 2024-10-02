## Helper functions for predicting from models

#' @title Predict to annual raster
#' @description Predicts expected and potential GPP (mean and 0.9 quantile)
#' for a raster of covariates in a given year 
#' 
#' @usage predict_to_annual_raster(year, static_vars, dynamic_vars, model, threads)
#' 
#' @param year numeric. Target year for prediction.
#' @param static_vars SpatRaster of temporally constant predictors
#' @param dynamic_vars SpatRaster of annually varying covariates with bands
#' named as {covariate}.{year}, e.g., "pptAnomaly.2004"
#' @param model ranger model object to use for prediction
#' @param threads numeric. Number of threads to use for prediction

predict_to_annual_raster <- function(year, static_vars, dynamic_vars, model, threads = NULL) {
  
  # Filter covariates to required year
  dynamic_vars_year <- subset(
    dynamic_vars,
    str_detect(names(dynamic_vars), as.character(year))
  )
  
  names(dynamic_vars_year) <- names(dynamic_vars_year) %>%
    str_replace("\\.\\d{4}", "")
  
  # Convert covariate rasters to data.frames
  predictors_all <- c(dynamic_vars_year, static_vars)
  
  predictors_all_df <- as.data.frame(predictors_all, xy = TRUE, na.rm = TRUE)
  
  # Predict from model
  message("Predicting mean GPP: ", year)
  gpp_predictions <- predict(model, predictors_all_df, type = "response", num.threads = threads)
  message("Predicting potential GPP: ", year)
  pot_gpp_predictions <- predict(model, predictors_all_df, type = "quantiles", quantiles = 0.9, num.threads = threads)
  message("Prediction complete: ", year)
  
  predictions <- predictors_all_df %>%
    mutate(
     gpp_predicted = gpp_predictions$predictions,
     potential_gpp_predicted = pot_gpp_predictions$predictions,
     rpi = GPP / potential_gpp_predicted
  )
  
  # Calculate RPI and convert back to raster
  predictions_rast <- predictions %>%
    select(x, y, GPP, gpp_predicted, potential_gpp_predicted, rpi) %>%
    rename_with(
      .cols = c("GPP", "gpp_predicted", "potential_gpp_predicted", "rpi"),
      .fn = ~ paste0(.x, ".", year)
    ) %>%
    rast(crs = crs(predictors_all), digits = 5)
  
  
  predictions_rast
}
