## Helper functions for data reading and pre-processing

#' @title Rename rasters
#' @description Appends a specified character string to the names of raster bands 
#' 
#' @usage rename_raster(raster, string)
#' 
#' @param raster a SpatRaster
#' @param string character or numeric to be appended to band names

rename_raster <- function(raster, string) {
  
  names(raster) <- paste0(names(raster), ".", string)
  
  raster
  
}

#' @title Read annual rasters
#' @description Reads in annual covariate rasters efficiently
#' 
#' @usage read_annual_rasters(years, folder)
#' 
#' @param years numeric. A vector of years to search for in the annual covariate filenames
#' @param folder character. The folder path containing the raster files

read_annual_rasters <- function(years, folder) {
  filenames <- Sys.glob(paste0(folder, "/*", years, "*.tif"))
  
  rasters <- filenames %>%
    map(rast) %>%
    map2(years, rename_raster) %>%
    rast()
    
  rasters
}

#' @title Mask zeros
#' @description Mask all layers in a spatRaster for which a designated layer is zero
#' 
#' @usage mask_zeros(raster, layer)
#' 
#' @param raster spatRaster. 
#' @param layer character. Layer name to use for the zero mask

mask_zeros <- function(raster, layer) {
  
  mask_layer <- raster[[layer]] != 0
  
  raster_masked <- mask(raster, mask_layer, maskvalues = c(0, NA))
  
  raster_masked
  
}

#' @title Tidy annual variables
#' @description Converts a data.frame of annual covariate values obtained
#' from terra::extract or terra::as.data.frame to a tidy data.frame with one
#' row per location and year
#' 
#' @usage tidy_annual_vars(data, delim = ".")
#' 
#' @param data a wide-format data.frame with column names "{var}.{year}", where
#' the delimiter "." can be changed with the "delim" argument.
#' @param delim character. A unique string separating the variable component
#' of the column name from the year. It cannot be present elsewhere in the column
#' names.

tidy_annual_vars <- function(data, delim = ".") {
  
  tidy_data <- data %>%
    pivot_longer(cols = contains(delim)) %>%
    separate_wider_delim(cols = "name", delim = delim, names = c("var", "year")) %>%
    mutate(year = as.numeric(year)) %>%
    pivot_wider(names_from = "var", values_from = "value")
  
  tidy_data
}
