#' Download Pretrained Model
#'
#' This function downloads the pretrained model from a GitHub repository.
#'
#' @param model_names Character vector specifying for which satellites you want to download the models. Any of "TM04","TM05","LE07","LC08","LC09" or "ALL"
#' @param model_dir Character, path to folder where the model objects will be stored.
#' @param timeout Maximum time to use to download the files (á ~80MB). Defaults to 300 seconds. Increase if very slow network.
#' @param force Logical. If FALSE (default) downloads the model object only if no previous model file exists.
#' @examples
#' \dontrun{
#' # The full workflow
#' library(snowman)
#' library(terra)
#' 
#' # Set the number of cores
#' n_workers <- 4
#' 
#' # Replace with your own path where all data will be downloaded
#' base_landsat_path <- "C:/MyTemp/RS/"
#' 
#' site <- "Sierra_nevada" # Name of the AOI
#' aoi_point <- list(lon = -3.311665, lat = 37.053188) # center point of AOI
#' 
#' # Download Landsat-8 imagery
#' image_df <- extract_landsat(aoi = aoi_point,
#'                             site_name = site,
#'                             base_landsat_dir = base_landsat_path,
#'                             sats = "LC08",
#'                             workers = n_workers)
#' 
#' # Calculate other geospatial information for the classifier
#' calc_predictors(image_df, site_name = site, base_landsat_dir = base_landsat_path)
#' 
#' # Download pretrained Random forest classifier
#' download_model(model_names = "LC08",
#'                model_dir = base_landsat_path)
#' 
#' # Run the classification across imagery
#' lss <- classify_landsat(image_df, 
#'                         site_name = site, 
#'                         base_landsat_dir = base_landsat_path, 
#'                         model_dir = base_landsat_path, 
#'                         workers = n_workers)
#' 
#' # Calculate snow variables over the AOI based on the classified imagery
#' snow_vars <- calc_snow_variables(image_df, 
#'                                  site_name = site, 
#'                                  base_landsat_dir = base_landsat_path, 
#'                                  workers = n_workers)
#' 
#' # Plot one of the resulting layers
#' plot(snow_vars$scd, col = rev(topo.colors(100)), 
#'      main = "Snow cover duration in Sierra Nevada")
#' 
#' # Save the resulting snow maps as GeoTiffs
#' writeRaster(snow_vars, paste0(base_landsat_path, "/", site, "/", "snow_variables.tif"), 
#'             datatype = "FLT4S")
#' }
#' @export
download_model <- function(model_names = c("TM04","TM05","LE07","LC08","LC09"), model_dir, timeout = 300, force = FALSE) {
  # URL of the GitHub repository
  options(timeout = max(timeout, getOption("timeout")))
  repo_url <- "https://github.com/poniitty/snowman_models/raw/main/"
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_TM05.rds")))){
    if(any(c("tm04","tm05","all") %in% tolower(model_names))){
      print("Downloading TM04 & TM05 model...")
      # Download the model object for Landsat-4 & Landsat 5
      download.file(paste0(repo_url, "Landsat_models/latest_TM05.rds"), 
                    destfile = paste0(model_dir, "/latest_TM05.rds"))
    }
  }
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_LE07.rds")))){
    if(any(c("le07","all") %in% tolower(model_names))){
      print("Downloading LE07 model...")
      # Download the model object for Landsat-7
      download.file(paste0(repo_url, "Landsat_models/latest_LE07.rds"), 
                    destfile = paste0(model_dir, "/latest_LE07.rds"))
    }
  }
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_LC08.rds")))){
    if(any(c("lc09","lc08","all") %in% tolower(model_names))){
      print("Downloading LC08 & LC09 model...")
      # Download the model object for Landsat-8 & Landsat 9
      download.file(paste0(repo_url, "Landsat_models/latest_LC08.rds"), 
                    destfile = paste0(model_dir, "/latest_LC08.rds"))
    }
  }
}
