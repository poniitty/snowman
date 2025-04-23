#' Extract image info
#'
#' This function extracts metadata from previously downloaded Landsat imagery.
#'
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' @param data_source A character string representing the data source (default is "rstac").
#' @return A tibble containing the metadata of the previously downloaded Landsat scenes.
#' @examples
#' \dontrun{
#' search_image_df(site_name = "ExampleSite")
#' }
#' @export
#' @import dplyr sf terra lubridate stringr parallel tibble readr
search_image_df <- function(site_name, base_landsat_dir, workers = 1){
  
  area_landsat_dir <- paste0(base_landsat_dir,"/",site_name, "/imagery")
  
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  
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
