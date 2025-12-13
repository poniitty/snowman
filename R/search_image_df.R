#' Extract image info
#'
#' This function extracts metadata from previously downloaded Landsat imagery.
#'
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing. 
#' Setting workers > 1 (default) enables parallel computing across multiple nodes.
#' @return A tibble containing the metadata of the previously downloaded Landsat scenes.
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
#' @import dplyr sf terra lubridate stringr parallel tibble readr
search_image_df <- function(site_name, base_landsat_dir, workers = 1){
  
  area_landsat_dir <- paste0(base_landsat_dir,"/",site_name, "/imagery")
  
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  
  if(file.exists(file.path(area_landsat_dir, "lss_final.csv"))){
    lss <- readr::read_csv(file.path(area_landsat_dir, "lss_final.csv"), show_col_types = FALSE)
    if(nrow(lss) == length(tifs)){
      return(lss)
    }
  }
  
  lss <- tibble(area = site_name,
                file = tifs)
  lss$collection <- gsub("0","C",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
  lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
  lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
  lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),1,3))
  lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])),4,6))
  lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
  lss$time <- gsub(".tif","",unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
  
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    # Use parLapply on Windows
    cl <- makeCluster(workers)
    on.exit(stopCluster(cl))
    lccs <- parLapply(cl, lss$file, calc_coverages, image_dir = area_landsat_dir) %>%
      bind_rows
  } else {
    # Use mclapply on unix
    lccs <- mclapply(lss$file, calc_coverages, image_dir = area_landsat_dir, mc.cores = workers) %>%
      bind_rows
  }
  
  lss <- full_join(lss,
                   lccs, by = "file")
  
  lss <- lss %>% arrange(date)
  write_csv(lss, paste0(area_landsat_dir, "/lss_final.csv"))
  
  return(lss)
}
