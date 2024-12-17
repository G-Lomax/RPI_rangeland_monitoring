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

#' @title Extract raster values for field data points
#' @description For field collected data points, extract the appropriate cell values
#' from a given SpatRaster
#' 
#' @usage raster_extract(geom, hYear = NULL, raster)
#' 
#' @param geom an sfc geometry column in an sf object containing data 
#' @param hydro_year numeric. The hydrological year for which to extract data.
#' If numeric, spatRaster band names will be filtered by that year string before
#' extraction. If NULL, all raster layers are extracted
#' @param raster spatRaster object to extract cell values from

raster_extract <- function(hydro_year, sf, raster) {
  
  sf_filter <- filter(sf, hYear == hydro_year)
  
  raster_hYear <- subset(raster, str_detect(names(raster), as.character(hydro_year)))
  names(raster_hYear) <- str_replace(names(raster_hYear), paste0("\\.", as.character(hydro_year)), "")
  
  extracted <- terra::extract(raster_hYear, vect(sf_filter), cells = TRUE, ID = FALSE)
  
  bind_cols(sf_filter, extracted)

}

#' @title Check for unrecorded bare soil
#' @description For landPKS data, create a new row recording a "bare" soil
#' reading if there is no perennial or annual grass cover present. Otherwise,
#' tree or shrub canopy presence prevents recording of bare soil
#' 
#' @usage check_bare(groups, covers = c("perennial", "annual"))
#' 
#' @param groups a column name in a grouped data.frame containing land cover records
#' @param covers character. Land covers considered non-bare


check_bare <- function(groups, covers = c("perennial", "annual")) {
  is_bare <- sum(groups %in% covers) == 0
  already_contains_bare <- "bare" %in% groups
  
  if (is_bare == TRUE & already_contains_bare == FALSE) {
    "bare"
  } else {
    NA
  }
}


#' @title Find centroid by group
#' @description Find the centroid of a group of points
#' 
#' @usage find_centroid(geometry)
#' 
#' @param geometry An sfc column in a data frame

find_centroid <- function(geometry) {
  x <- st_coordinates(geometry)[,1]
  y <- st_coordinates(geometry)[,2]
  
  centroid_coords <- c(mean(x), mean(y))
  
  centroid_point <- st_point(centroid_coords)
  
  st_sfc(centroid_point, crs = "EPSG:4326")
  
}
