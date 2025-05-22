#' Extracts Landsat Imagery
#'
#' This function extracts Landsat imagery that fulfills the specified arguments.
#' It processes the area of interest (AOI), transforms coordinates, and downloads
#' the relevant Landsat scenes.
#'
#' @param aoi An `sf` object or a list containing the area of interest.
#' @param site_name A character string representing the name of the site.
#' @param aoi_size A numeric value representing the size of the AOI in kilometers.
#' @param start_date A character string representing the start date for the imagery.
#' @param end_date A character string representing the end date for the imagery.
#' @param months A numeric vector representing the months to include (default is all months, so 1:12).
#' @param sats A character string representing the satellite to use (default is Landsat-8 only, so "LC08").
#' @param minclouds A numeric value representing the minimum cloud cover percentage (default is 50).
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing. 
#' Setting workers > 1 (default) enables parallel computing across multiple nodes.
#' @param data_source A character string representing the data source (currently only option is "rstac").
#' @param force Logical, if TRUE deletes all previously downloaded imagery.
#' @return A tibble containing the metadata of the downloaded Landsat scenes.
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
#' @import rstac dplyr sf terra lubridate stringr parallel tibble readr ggplot2
extract_landsat <- function(aoi,
                            site_name,
                            aoi_size = 2,
                            start_date = "2020-01-01",
                            end_date = "2024-12-01",
                            months = 1:12,
                            sats = "LC08",
                            minclouds = 50,
                            base_landsat_dir,
                            workers = 1,
                            data_source = "rstac",
                            force = FALSE) {
  
  # load(file = "data/utm_zones.rda")
  utmall <- utm_zones
  
  suppressWarnings({
    if (inherits(aoi, "sf")) {
      aoi_mid <- aoi %>%
        st_centroid() %>%
        st_transform(crs = 4326)
    } else if (is.list(aoi)) {
      aoi_mid <- as_tibble(aoi) %>%
        mutate(name = site_name) %>%
        st_as_sf(coords = c("lon", "lat"), crs = 4326)
      aoi <- as_tibble(aoi) %>%
        mutate(name = site_name) %>%
        st_as_sf(coords = c("lon", "lat"), crs = 4326)
    } else {
      stop("aoi must be an sf object or a list")
    }
  })
  
  # WGS84 UTM zones to set the correct projection
  utm <- utmall[aoi_mid,] # Which zone the study points falls in
  
  if (nrow(utm) == 0) {
    utm <- utmall[st_nearest_feature(aoi_mid, utmall),]
  }
  
  lat <- st_coordinates(aoi_mid)[,"Y"]
  utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0", utm$ZONE), utm$ZONE)
  epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
  
  # From point to polygon (20km x 20km)
  if (st_geometry_type(aoi) == "POLYGON") {
    aoi <- aoi %>% st_transform(crs = epsg) %>%
      mutate(name = site_name)
  } else {
    aoi <- aoi %>% st_transform(crs = epsg) %>%
      st_buffer(aoi_size * 1000) %>% st_bbox() %>%
      st_as_sfc() %>% st_as_sf() %>%
      mutate(name = site_name)
  }
  
  # Create sub-directory where imagery will be saved if it does not exist
  area_landsat_dir <- paste0(base_landsat_dir, "/", site_name)
  if (!dir.exists(area_landsat_dir)) {
    dir.create(area_landsat_dir)
  }
  area_landsat_dir <- paste0(area_landsat_dir,"/imagery")
  if (!dir.exists(area_landsat_dir)) {
    dir.create(area_landsat_dir)
  }
  
  # If opted, delete all previously downloaded images
  if(force == TRUE){
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$", full.names = TRUE)
    unlink(tifs)
    unlink(paste0(area_landsat_dir, "/lss.csv"))
    unlink(paste0(area_landsat_dir, "/lss_final.csv"))
  }
  
  # List existing files
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  tif_dates <- tibble(satid = unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][1])),
                      date = as.character(ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][4])))))
  if (nrow(tif_dates) == 0) {
    tif_dates <- tibble(satid = NA,
                        date = NA)
  }
  
  if (data_source == "rstac") {
    lss <- extract_landsat_stac(workers = workers, start_date = start_date, end_date = end_date, months = months,
                                aoi = aoi, site_name = site_name, epsg = epsg, excl_dates = tif_dates,
                                sats = sats, area_landsat_dir = area_landsat_dir, minclouds = minclouds)
  }
  
  if (is.null(lss)) {
    lss <- tibble(id = NULL,
                  DATE_ACQUIRED = NULL)
  } else {
    lss <- lss %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
  }
  if (file.exists(paste0(area_landsat_dir, "/lss.csv"))) {
    lss <- bind_rows(read_csv(paste0(area_landsat_dir, "/lss.csv"), show_col_types = FALSE) %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED)),
                     lss) %>%
      distinct()
  }
  
  if (nrow(lss) == 0) {
    print("Zero Landsat scenes found and downloaded!")
    return(NULL)
  } else {
    write_csv(lss, paste0(area_landsat_dir, "/lss.csv"))
    
    # List downloaded images
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
    
    lss <- full_join(tibble(file = tifs),
                     lss, by = join_by(file))
    
    lss$area <- site_name
    lss$collection <- gsub("0", "C", unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
    lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
    lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
    lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])), 1, 3))
    lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])), 4, 6))
    lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
    lss$time <- gsub(".tif", "", unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
    
    lss <- lss %>% arrange(date)
    
    # Check if all rasters are working. Remove corrupted files.
    
    os <- Sys.info()["sysname"]
    
    if (os == "Windows") {
      # Use future_lapply on Windows
      plan(multisession, workers = workers)
      img_remove <- unlist(lapply(lss$file, check_raster, image_dir = area_landsat_dir))
    } else {
      # Use mclapply on unix
      img_remove <- unlist(mclapply(lss$file, check_raster, image_dir = area_landsat_dir, mc.cores = workers))
    }
    
    
    if (length(img_remove) > 0) {
      unlink(paste0(area_landsat_dir, "/", img_remove))
      tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
      tif_dates <- tibble(satid = unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][2])),
                          date = as.character(ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][5])))))
      if (nrow(tif_dates) == 0) {
        tif_dates <- tibble(satid = NA,
                            date = NA)
      }
      if (data_source == "rstac") {
        lss <- extract_landsat_stac(workers = workers, start_date = start_date, end_date = end_date, months = months,
                                    aoi = aoi, site_name = site_name, epsg = epsg, excl_dates = tif_dates,
                                    sats = sats, area_landsat_dir = area_landsat_dir, minclouds = minclouds)
      }
      
      if (is.null(lss)) {
        lss <- tibble(id = NULL,
                      DATE_ACQUIRED = NULL)
      } else {
        lss <- lss %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
      }
      if (file.exists(paste0(area_landsat_dir, "/lss.csv"))) {
        lss <- bind_rows(read_csv(paste0(area_landsat_dir, "/lss.csv"), show_col_types = FALSE) %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED)),
                         lss) %>%
          distinct()
      }
      write_csv(lss, paste0(area_landsat_dir, "/lss.csv"))
    }
    
    # List downloaded images
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
    
    lss <- full_join(tibble(file = tifs),
                     lss, by = join_by(file))
    
    lss$area <- site_name
    lss$collection <- gsub("0", "C", unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][6])))
    lss$tier <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][7]))
    lss$satid <- unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][1]))
    lss$path <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])), 1, 3))
    lss$row <- as.numeric(substr(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][3])), 4, 6))
    lss$date <- ymd(unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][4])))
    lss$time <- gsub(".tif", "", unlist(lapply(lss$file, function(x) str_split(x, "_")[[1]][8])))
    
    lss <- lss %>% arrange(date) %>% filter(!is.na(file))
    
    # Check if all rasters are working. Remove corrupted files.
    if (os == "Windows") {
      img_remove <- unlist(lapply(lss$file, check_raster, image_dir = area_landsat_dir))
    } else {
      # Use mclapply on unix
      img_remove <- unlist(mclapply(lss$file, check_raster, image_dir = area_landsat_dir, mc.cores = workers))
    }
    
    if (length(img_remove) > 0) {
      print(paste0(length(img_remove), " raster(s) not functional. REMOVED!!"))
    }
    
    unlink(paste0(area_landsat_dir, "/", img_remove))
    
    lss <- lss %>%
      filter(!file %in% img_remove)
    
    lss_d <- lss %>% group_by(date, satid, path) %>% count %>% filter(n > 1) %>% ungroup()
    if (nrow(lss_d) > 0) {
      for (ii in seq_len(nrow(lss_d))) {
        lss_dd <- lss_d %>% slice(ii)
        
        lss_dd <- right_join(lss, lss_dd, by = join_by(satid, path, date))
        
        rs <- lapply(lss_dd$file, function(x) {
          r <- rast(paste0(area_landsat_dir, "/", x))
          r[r == 0] <- NA
          return(r)
        })
        
        rs <- sprc(rs)
        rs <- mosaic(rs) %>% round
        writeRaster(rs, paste0(area_landsat_dir, "/", lss_dd$file[[1]]),
                    overwrite = TRUE, datatype = "INT2U")
        unlink(paste0(area_landsat_dir, "/", lss_dd$file[[2:nrow(lss_dd)]]))
        lss <- lss %>% filter(!file %in% lss_dd$file[[2:nrow(lss_dd)]])
      }
    }
    
    # First round of selection
    # remove images with very limited clear coverage over the AOI
    if (os == "Windows") {
      cl <- makeCluster(workers)
      on.exit(stopCluster(cl))
      lccs <- parLapply(cl, lss$file, calc_coverages, image_dir = area_landsat_dir) %>%
        bind_rows
    } else {
      # Use mclapply on unix
      lccs <- mclapply(lss$file, calc_coverages, image_dir = area_landsat_dir, mc.cores = workers) %>%
        bind_rows
    }
    
    lss <- full_join(lss %>% select(-ends_with("proportion")), lccs, by = "file")
    
    lss %>% arrange(desc(fill_proportion)) %>%
      filter(clear_proportion < 0.1) %>% pull(file) -> img_remove
    
    unlink(paste0(area_landsat_dir, "/", img_remove))
    
    lss <- lss %>%
      filter(!file %in% img_remove)
    
    unlink(list.files(tempdir(), full.names = T))
    write_csv(lss, paste0(area_landsat_dir, "/lss_final.csv"))
    
    return(lss)
  }
}

# Internal Function to Create a Base URL with Microsoft Planetary Computer
make_vsicurl_url <- function(base_url) {
  paste0(
    "/vsicurl",
    "?pc_url_signing=yes",
    "&pc_collection=landsat-c2-l2",
    "&url=",
    base_url
  )
}

# Internal Function to Extract Landsat Imagery from STAC
extract_landsat_stac <- function(aoi,
                                 epsg,
                                 excl_dates,
                                 site_name,
                                 sats = c("LT04", "LT05", "LE07", "LC08", "LC09"),
                                 start_date,
                                 end_date,
                                 months = 1:12,
                                 minclouds,
                                 area_landsat_dir,
                                 workers) {
  
  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  
  sats2 <- sats %>%
    ifelse(. == "LT04", "landsat-4", .) %>%
    ifelse(. == "LT05", "landsat-5", .) %>%
    ifelse(. == "LE07", "landsat-7", .) %>%
    ifelse(. == "LC08", "landsat-8", .) %>%
    ifelse(. == "LC09", "landsat-9", .)
  
  it_obj <- s_obj %>%
    stac_search(collections = "landsat-c2-l2",
                bbox = st_bbox(aoi %>% st_transform(4326)),
                datetime = paste0(start_date, "/", end_date),
                limit = 1000) %>%
    get_request()
  
  if (length(it_obj$features) == 1000) {
    dr <- tibble(date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day")) %>%
      mutate(gr = cut_number(row_number(.), 20)) %>%
      group_by(gr) %>%
      group_split()
    
    it_obj <- s_obj %>%
      stac_search(collections = "landsat-c2-l2",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  datetime = paste0(min(dr[[1]]$date), "/", max(dr[[1]]$date)),
                  limit = 1000) %>%
      get_request()
    
    for (i in 2:length(dr)) {
      it_obj2 <- s_obj %>%
        stac_search(collections = "landsat-c2-l2",
                    bbox = st_bbox(aoi %>% st_transform(4326)),
                    datetime = paste0(min(dr[[i]]$date), "/", max(dr[[i]]$date)),
                    limit = 1000) %>%
        get_request()
      
      it_obj$features <- c(it_obj$features, it_obj2$features)
    }
  }
  
  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) { x$properties$`eo:cloud_cover` < minclouds }) %>%
    items_filter(filter_fn = function(x) { x$properties$platform %in% sats2 }) %>%
    assets_select(asset_names = c("coastal", "blue", "green", "red", "nir08", "swir16", "swir22", "lwir", "lwir11", "qa_pixel", "qa_radsat")) %>%
    items_filter(filter_fn = function(x) { month(ymd_hms(x$properties$datetime)) %in% months }) %>%
    items_filter(filter_fn = function(x) { !as_date(ymd_hms(x$properties$datetime)) %in% excl_dates$date })
  
  if (length(it_obj$features) > 0) {
    
    itst <- items_as_sf(it_obj) %>% st_make_valid()
    itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
    
    itst$area <- as.numeric(st_area(itst)) / (1000 * 1000)
    
    tt <- itst %>%
      mutate(date = as_date(ymd_hms(datetime))) %>%
      group_by(date, instruments, `landsat:wrs_path`, .add = TRUE) %>%
      mutate(n = n()) %>%
      arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()
    
    suppressMessages({
      suppressWarnings({
        tt <- lapply(tt, function(x) {
          if (nrow(x) > 1) {
            if (diff(x %>% pull(area)) == 0) {
              x <- x %>% slice(1)
            } else {
              stint <- st_intersection(x) %>% st_make_valid() %>% slice(2) %>% st_area %>% as.numeric()
              x <- x %>% filter(area > stint / (1000 * 1000) + 1)
            }
          }
          return(x)
        }) %>% bind_rows()
      })
    })
    
    itst <- itst %>%
      filter(`landsat:scene_id` %in% (tt %>% pull(`landsat:scene_id`)))
    
    it_obj <- it_obj %>%
      items_filter(filter_fn = function(x) { x$properties$`landsat:scene_id` %in% (itst %>% pull(`landsat:scene_id`)) })
    
    
    juuh <- process_features_in_parallel(it_obj, area_landsat_dir, aoi, workers)
    
    itst$platform <- itst$platform %>%
      ifelse(. == "landsat-4", "LT04", .) %>%
      ifelse(. == "landsat-5", "LT05", .) %>%
      ifelse(. == "landsat-7", "LE07", .) %>%
      ifelse(. == "landsat-8", "LC08", .) %>%
      ifelse(. == "landsat-9", "LC09", .)
    
    tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
    
    idsALL <- itst %>%
      mutate(date = gsub("-", "", as_date(ymd_hms(datetime)))) %>%
      mutate(date2 = gsub("-", "", as_date(ymd_hms(created)))) %>%
      mutate(time = paste0(gsub(":", "", substr(datetime, 12, 19)), "GMT.tif")) %>%
      mutate(file = paste0(platform, "_", `landsat:correction`, "_",
                           `landsat:wrs_path`, `landsat:wrs_row`, "_",
                           date, "_", date2, "_",
                           `landsat:collection_number`, "_", `landsat:collection_category`, "_",
                           time)) %>%
      select(file, date, `eo:cloud_cover`, `view:sun_elevation`) %>%
      rename(DATE_ACQUIRED = date,
             CLOUD_COVER = `eo:cloud_cover`,
             SUN_ELEVATION = `view:sun_elevation`) %>%
      st_drop_geometry() %>%
      filter(file %in% tifs) %>%
      mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>%
      rownames_to_column("id")
  } else {
    print("No new Landsat scenes downloaded!")
    idsALL <- NULL
  }
  
  return(idsALL)
}

# Internal Function to Extract Landsat Imagery from STAC
process_feature <- function(ft, area_landsat_dir, aoi, tempdir) {
  # ft <- it_obj$features[[1]]
  
  make_vsicurl_url <- function(base_url) {
    paste0(
      "/vsicurl",
      "?pc_url_signing=yes",
      "&pc_collection=landsat-c2-l2",
      "&url=",
      base_url
    )
  }
  
  nm <- paste0(paste(str_split(ft$id, "_")[[1]][1:4], collapse = "_"), "_",
               gsub("-", "", substr(ft$properties$created, 1, 10)), "_",
               paste(str_split(ft$id, "_")[[1]][5:6], collapse = "_"), "_",
               gsub(":", "", substr(ft$properties$datetime, 12, 19)), "GMT.tif")
  
  if (!file.exists(paste0(area_landsat_dir, "/", nm))) {
    full_url <- make_vsicurl_url(assets_url(ft) %>% sort)
    file_names <- gsub("TIF$", "tif", basename(full_url))
    
    juuh <- lapply(seq_len(length(full_url)), function(nr) {
      # nr <- 1
      e <- try({
        gdal_utils(
          "warp",
          source = full_url[[nr]],
          destination = paste0(tempdir, "/", file_names[[nr]]),
          options = c(
            "-t_srs", st_crs(aoi)$wkt,
            "-te", st_bbox(aoi),
            "-tr", c(30, 30)
          )
        )
      }, silent = TRUE)
      if (class(e)[[1]] == "try-error") {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
    
    err <- file_names[!unlist(juuh)]
    if (length(err) > 0) {
      ll <- lapply(err, function(xx) {
        r <- rast(paste0(tempdir, "/", file_names[unlist(juuh)][1]))
        r[] <- NA
        writeRaster(r, paste0(tempdir, "/", xx), datatype = "INT2U", overwrite = TRUE)
      })
    }
    
    r <- rast(c(paste0(tempdir, "/", file_names[grepl("_SR_", file_names)]),
                paste0(tempdir, "/", file_names[grepl("_ST_", file_names)]),
                paste0(tempdir, "/", file_names[grepl("_QA_", file_names)])))
    names(r) <- lapply(names(r), function(x) paste(str_split(x, "_")[[1]][8:9], collapse = "_")) %>% unlist
    
    writeRaster(r, paste0(area_landsat_dir, "/", nm), datatype = "INT2U", overwrite = TRUE)
    
    unlink(c(paste0(tempdir, "/", file_names)))
  }
}

# Internal Function to Extract Landsat Imagery from STAC
process_features_in_parallel <- function(it_obj, area_landsat_dir, aoi, workers) {
  # Check the operating system
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    # Use parLapply on Windows
    cl <- makeCluster(workers)
    on.exit(stopCluster(cl))
    
    # Load required packages on each worker node
    clusterEvalQ(cl, {
      library(sf)
      library(terra)
      library(stringr)
      library(rstac)
      library(httr)
    })
    
    juuh <- parLapply(cl, it_obj$features, process_feature, area_landsat_dir, aoi, tempdir())
  } else {
    # Use mclapply on unix
    juuh <- mclapply(it_obj$features, process_feature, area_landsat_dir, aoi, tempdir(), mc.cores = workers)
  }
  
  return(juuh)
}

# Internal Function to check if raster files work properly
check_raster <- function(image, image_dir){
  
  xx <- try(terra::rast(paste0(image_dir,"/",image)))
  
  if(class(xx)[1] == "try-error"){
    return(image)
  } else {
    
    xx <- try(terra::values(terra::rast(paste0(image_dir,"/",image))[[1]]))
    
    if(class(xx)[1] == "try-error"){
      return(image)
    } else {
      return(NULL)
    }
  }
}

# Internal Function to calculate cloud cover within the rasters
calc_coverages <- function(image, image_dir){
  # image <- lss$file[[5]]
  require(terra)
  
  rs <- rast(paste0(image_dir,"/",image))
  rsn <- names(rs)
  # plot(rs)
  rs[[4]][is.na(rs[[4]])] <- 0
  
  cmask <- rs[[1]]
  cmask[] <- unlist(lapply(as.numeric(values(rs[[rsn[grepl("_pixel",rsn, ignore.case = T)]]])), 
                           function(x) {as.numeric(paste(as.numeric(intToBits(x)[9:10]), collapse = ""))}))
  
  fill_prop <- mean(values(rs[[4]]) == 0)
  rs[[4]][cmask == 1] <- NA
  cloud_prop <- mean(is.na(values(rs[[4]])))
  
  df <- dplyr::tibble(file = image,
               fill_proportion = fill_prop,
               cloud_proportion = cloud_prop,
               clear_proportion = 1-(fill_prop+cloud_prop))
  
  return(df)
  
}
