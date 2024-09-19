load <- function() {
  # load all the packages you need for the project
  message("Loading packages...")
  suppressPackageStartupMessages(
    source("scripts/packages.R")
  )
  
  # read the functions from R directory
  message("Loading functions...")
  function_files <- list.files("R", pattern = ".R$", full.names = TRUE)
  walk(function_files, source)
  message("Complete")
  
  # set global options
  terra::terraOptions(
    memfrac = 0.8,
    memmax = NA,
    tempdir = "data/processed/tmp"
  )
  
  options(scipen = 999)
  
}

load()
