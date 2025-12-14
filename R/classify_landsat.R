#' Classify Landsat Imagery
#'
#' This function classifies Landsat imagery using pretrained models. It processes the input image data frame,
#' applies the models, and generates classification rasters.
#'
#' @param image_df A tibble, the output from the `extract_landsat` or `search_image_df` functions.
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param model_dir A character string representing the path to the directory containing the pretrained models.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' @param force Logical, if TRUE classification will be run also for images with previously classified.
#' @return A list of file paths to the classified rasters.
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
#' @import rstac dplyr sf terra lubridate stringr parallel tibble solartime randomForestSRC
classify_landsat <- function(image_df, site_name, base_landsat_dir, model_dir, workers, force = FALSE) {
  # Define directories for predictors and classifications
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  class_landsat_dir <- paste0(base_landsat_dir, "/", site_name, "/classifications")
  
  # Create the classifications directory if it does not exist
  if (!dir.exists(class_landsat_dir)) {
    dir.create(class_landsat_dir)
  }
  
  # Load pretrained models
  if(any(c("LT05","LT04") %in% unique(image_df$satid))){
    mod5 <- readRDS(paste0(model_dir, "latest_TM05.rds"))
  } else {
    mod5 <- NULL
  }
  if(any(c("LE07") %in% unique(image_df$satid))){
    mod7 <- readRDS(paste0(model_dir, "latest_LE07.rds"))
  } else {
    mod7 <- NULL
  }
  if(any(c("LC08","LC09") %in% unique(image_df$satid))){
    mod8 <- readRDS(paste0(model_dir, "latest_LC08.rds"))
  } else {
    mod8 <- NULL
  }
  
  # --- 1. PORTABLE PARALLEL SETUP (Final, Stable Version) ---
  # Check if running inside RStudio/RStudio Server (where multicore is unstable)
  # We use parallelly::supportsMulticore() == FALSE on Linux/Darwin to detect RStudio Server.
  is_rstudio_env <- parallelly::supportsMulticore() == FALSE && Sys.info()['sysname'] %in% c("Linux", "Darwin")
  
  # Determine how many workers are currently running in the active plan.
  current_plan_workers <- future::nbrOfWorkers()
  
  # We only set a new plan if:
  # 1. The user requested parallelism (workers > 1) AND
  # 2. No parallel workers are currently active (current_plan_workers <= 1).
  if (workers > 1 && current_plan_workers <= 1) {
    
    # Store the original plan state (which is currently sequential/idle)
    # We save this so snowman_cleanup() can restore it later.
    oplan <- future::plan() 
    
    os_type <- Sys.info()['sysname']
    
    # Decision: Use multisession for stability in RStudio/Windows, multicore for speed elsewhere.
    if (os_type == "Windows" || is_rstudio_env) {
      # USE MULTISESSION (SOCKETS): Stable for RStudio Server and required for Windows.
      future::plan(future::multisession, workers = workers)
      message(paste("Setting parallel plan to 'multisession' (Sockets) with", workers, "workers for cross-platform stability."))
    } else if (os_type %in% c("Linux", "Darwin")) {
      # USE MULTICORE (FORKING): Fastest and best for pure Linux terminal/HPC batch environments.
      future::plan(future::multicore, workers = workers)
      message(paste("Setting parallel plan to 'multicore' (Forking) with", workers, "workers."))
    }
    
    # Note: We rely on the user calling snowman_cleanup() to restore 'oplan'.
    
  } else if (workers > 1 && current_plan_workers > 1) {
    # If a plan is already set with workers, reuse it.
    message(paste("Reusing existing parallel plan with", current_plan_workers, "workers."))
  } else {
    # Default to sequential if workers <= 1.
    future::plan(future::sequential)
    # message is optional here, but good for debugging if workers=1.
  }
  # --- END OF SETUP ---
  
  # Apply the classifying function to each image file in parallel
  lss <- future.apply::future_lapply(image_df$file, function(f) {
    classifying_function(
      imageid = f, 
      image_df = image_df,
      predictor_dir = predictor_dir, 
      class_landsat_dir = class_landsat_dir, 
      base_landsat_dir = base_landsat_dir, 
      site_name = site_name, 
      force = force,
      mod7 = mod7, mod5 = mod5, mod8 = mod8
    )},
    future.packages = c("dplyr", "sf", "terra", "lubridate",
                        "solartime", "randomForestSRC"),
    future.seed = TRUE,
    future.scheduling = 10
  )
  
  snowman_cleanup()
  
  return(unlist(lss))
}

# Internal Function to Classify a Single Landsat Image
classifying_function <- function(imageid, image_df, predictor_dir, class_landsat_dir,
                                 base_landsat_dir, site_name, force = FALSE,
                                 mod8, mod7, mod5) {
  # imageid <- image_df$file[[1]]
  # Check if the classified raster already exists
  if (!file.exists(paste0(class_landsat_dir, "/", imageid)) | force == TRUE) {
    # Get the satellite ID from the image metadata
    satid <- image_df %>% dplyr::filter(file == imageid) %>% dplyr::pull(satid)
    
    # Select the appropriate model variables based on the satellite ID
    if (satid %in% c("LT05", "LT04")) {
      mod_vars <- mod5$xvar.names
    } else if (satid == "LE07") {
      mod_vars <- mod7$xvar.names
    } else if (satid %in% c("LC08", "LC09")) {
      mod_vars <- mod8$xvar.names
    }
    
    # Read the Landsat image raster
    r <- terra::rast(paste0(base_landsat_dir, "/", site_name, "/imagery/", imageid))
    
    # Set the band names for Landsat 8 and 9 imagery
    if (satid %in% c("LC08", "LC09")) {
      names(r) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
    } else {
      rr <- r[[1]]
      rr[] <- 0
      r <- c(rr, r)
      names(r) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
    }
    
    # Mask pixels with QA value of 0
    r[r[["QA"]] == 0] <- NA
    
    # Read and project the DEM raster
    dem <- terra::rast(paste0(predictor_dir, "/ALOSDEM.tif"))
    
    # Add predictor rasters to the raster stack if they are used in the model
    if ("esalc" %in% mod_vars) {
      r[["esalc"]] <- terra::rast(paste0(predictor_dir, "/ESALC.tif"))
    }
    if ("esaw" %in% mod_vars) {
      r[["esaw"]] <- terra::rast(paste0(predictor_dir, "/ESAW.tif"))
    }
    if ("slope" %in% mod_vars) {
      r[["slope"]] <- terra::rast(paste0(predictor_dir, "/slope.tif")) / 100
    }
    
    # Check if the number of non-NA pixels is greater than 100
    if (sum(!is.na(terra::values(r[[1]]))) > 100) {
      # Calculate hillshade if it is used in the model
      if ("hillshade" %in% mod_vars) {
        lat <- st_coordinates(st_transform(st_centroid(st_as_sf(st_as_sfc(st_bbox(dem)))), 4326))[1, "Y"]
        lon <- st_coordinates(st_transform(st_centroid(st_as_sf(st_as_sfc(st_bbox(dem)))), 4326))[1, "X"]
        time <- ymd_hms(paste(str_split(gsub(".tif", "", imageid), "_")[[1]][4],
                              str_split(gsub(".tif", "", imageid), "_")[[1]][8]))
        
        sa <- computeSunPositionDoyHour(doy = yday(time), hour = hour(time)+(minute(time)/60), 
                                        latDeg = lat, longDeg = lon, timeZone = 0L)
        
        hill <- shade(terra::terrain(dem, "slope", unit = "radians"),
                      terra::terrain(dem, "aspect", unit = "radians"),
                      (sa[3] * (180 / pi)), (sa[4] * (180 / pi)))
        esaw <- terra::rast(paste0(predictor_dir, "/ESAW.tif"))
        hill[esaw > 0.2] <- median(terra::values(hill, mat = FALSE), na.rm = TRUE)
        
        if (sum(is.na(terra::values(hill, mat = FALSE))) > 0) {
          repeat {
            hill <- terra::focal(hill, 5, "mean", na.policy = "only", na.rm = TRUE)
            if (sum(is.na(terra::values(hill, mat = FALSE))) == 0) {
              break
            }
          }
        }
        
        r[["hillshade"]] <- hill
      }
      
      # Calculate median Landsat differences if they are used in the model
      idate <- as.character(lubridate::ymd(stringr::str_split(imageid, "_")[[1]][4]))
      mo <- lubridate::month(idate)
      e <- try({
        mr <- terra::rast(paste0(predictor_dir, "/medianlandsat_", mo, ".tif"))
      })
      if (class(e) == "try-error") {
        e <- try({
          mr <- terra::rast(paste0(predictor_dir, "/medianlandsat_", mo + 1, ".tif"))
        })
      }
      if (class(e) == "try-error") {
        e <- try({
          mr <- terra::rast(paste0(predictor_dir, "/medianlandsat_", mo - 1, ".tif"))
        })
      }
      
      mr_diff <- mr - r[[1:7]]
      names(mr_diff) <- paste0("mr_diff_", names(mr_diff))
      if (any(grepl("mr_diff_B5_std5", mod_vars))) {
        mr_diff[["mr_diff_B5_std5"]] <- terra::focal(r[["B5"]], w = 5, fun = "sd") - terra::focal(mr[["B5"]], w = 5, fun = "sd")
      }
      if (any(grepl("mr_diff_B7_std5", mod_vars))) {
        mr_diff[["mr_diff_B7_std5"]] <- terra::focal(r[["B7"]], w = 5, fun = "sd") - terra::focal(mr[["B7"]], w = 5, fun = "sd")
      }
      
      # Calculate median indices if they are used in the model
      mi <- terra::rast(paste0(predictor_dir, "/medianindices.tif"))
      
      # Calculate focal statistics if they are used in the model
      filter <- matrix(1, nrow = 5, ncol = 5)
      filter[ceiling(length(filter) / 2)] <- 0
      
      if (any(grepl("B2_tpi", mod_vars))) {
        r[["B2_tpi"]] <- r[["B2"]] - terra::focal(r[["B2"]], w = filter, fun = mean,
                                                  na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B5_tpi", mod_vars))) {
        r[["B5_tpi"]] <- r[["B5"]] - terra::focal(r[["B5"]], w = filter, fun = mean,
                                                  na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B7_tpi", mod_vars))) {
        r[["B7_tpi"]] <- r[["B7"]] - terra::focal(r[["B7"]], w = filter, fun = mean,
                                                  na.policy = "omit", na.rm = TRUE)
      }
      
      if (any(grepl("B2_min5", mod_vars))) {
        r[["B2_min5"]] <- terra::focal(r[["B2"]], w = 5, fun = "min", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B5_min5", mod_vars))) {
        r[["B5_min5"]] <- terra::focal(r[["B5"]], w = 5, fun = "min", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B7_min5", mod_vars))) {
        r[["B7_min5"]] <- terra::focal(r[["B7"]], w = 5, fun = "min", na.policy = "omit", na.rm = TRUE)
      }
      
      if (any(grepl("B2_max5", mod_vars))) {
        r[["B2_max5"]] <- terra::focal(r[["B2"]], w = 5, fun = "max", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B5_max5", mod_vars))) {
        r[["B5_max5"]] <- terra::focal(r[["B5"]], w = 5, fun = "max", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B7_max5", mod_vars))) {
        r[["B7_max5"]] <- terra::focal(r[["B7"]], w = 5, fun = "max", na.policy = "omit", na.rm = TRUE)
      }
      
      if (any(grepl("B2_std5", mod_vars))) {
        r[["B2_std5"]] <- terra::focal(r[["B2"]], w = 5, fun = "sd", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B5_std5", mod_vars))) {
        r[["B5_std5"]] <- terra::focal(r[["B5"]], w = 5, fun = "sd", na.policy = "omit", na.rm = TRUE)
      }
      if (any(grepl("B7_std5", mod_vars))) {
        r[["B7_std5"]] <- terra::focal(r[["B7"]], w = 5, fun = "sd", na.policy = "omit", na.rm = TRUE)
      }
      
      # Calculate previous classification probabilities if they are used in the model
      rrrr <- terra::rast(ncols = 180, nrows = 180, xmin = 0)
      fm <- terra::focalMat(rrrr, 5, "circle")
      fm[fm > 0] <- 1
      
      if ("cirrus_probs5" %in% mod_vars) {
        cirrus_probs5 <- r[[1]]
        cirrus_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_cirrus))
        r[["cirrus_probs5"]] <- terra::focal(cirrus_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("cloud_probs5" %in% mod_vars) {
        cloud_probs5 <- r[[1]]
        cloud_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_cloud))
        r[["cloud_probs5"]] <- terra::focal(cloud_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("cshadow_probs5" %in% mod_vars) {
        cshadow_probs5 <- r[[1]]
        cshadow_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_cshadow))
        r[["cshadow_probs5"]] <- terra::focal(cshadow_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("snow_probs5" %in% mod_vars) {
        snow_probs5 <- r[[1]]
        snow_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_snow))
        r[["snow_probs5"]] <- terra::focal(snow_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("water_probs5" %in% mod_vars) {
        water_probs5 <- r[[1]]
        water_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_water))
        r[["water_probs5"]] <- terra::focal(water_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("clear_probs5" %in% mod_vars) {
        clear_probs5 <- r[[1]]
        clear_probs5[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_clear))
        r[["clear_probs5"]] <- terra::focal(clear_probs5, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      if ("cloud_probs25" %in% mod_vars) {
        fm <- terra::focalMat(rrrr, 25, "circle")
        fm[fm > 0] <- 1
        
        cloud_probs25 <- r[[1]]
        cloud_probs25[] <- unlist(lapply(as.numeric(terra::values(r[["QA"]])), mask_cloud))
        r[["cloud_probs25"]] <- terra::focal(cloud_probs25, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
      }
      
      # Combine all rasters into a single raster stack
      r <- c(r, mr_diff, mi)
      r[["QA"]][is.na(r[["B1"]])] <- NA
      r[["QA"]][is.na(r[["B2"]])] <- NA
      r[["QA"]][is.na(r[["B3"]])] <- NA
      r[["QA"]][is.na(r[["B4"]])] <- NA
      r[["QA"]][is.na(r[["B5"]])] <- NA
      r[["QA"]][is.na(r[["B6"]])] <- NA
      r[["QA"]][is.na(r[["B7"]])] <- NA
      
      # Convert the raster stack to a data frame and filter out rows with NA values in the QA band
      d <- as.data.frame(r, na.rm = FALSE) %>%
        dplyr::filter(!is.na(QA)) %>%
        dplyr::mutate(dplyr::across(dplyr::ends_with("tpi"), ~ ifelse(is.na(.x), 0, .x)))
      
      if(nrow(d) == 0){
        
        rr <- r[["QA"]]
        rr[["class"]] <- NA
        rr[["land"]] <- NA
        rr[["water"]] <- NA
        rr[["snow"]] <- NA
        rr[["cloud"]] <- NA
        rr[["artif"]] <- NA
        
        rr <- rr[[-1]]
        # plot(rr)
        terra::writeRaster(rr, paste0(class_landsat_dir, "/", imageid),
                           datatype = "INT1U", overwrite = T)
        
        return(imageid)
        
      } else {
        layers <- c("B1","B2","B3","B4","B5","B6","B7",
                    names(d %>% dplyr::select(dplyr::ends_with("_min5"), dplyr::ends_with("_max5"), dplyr::ends_with("_std5")) %>% 
                            dplyr::select(-starts_with("B2"))))
        for(i in layers){
          # print(i)
          layers <- layers[-1]
          for(ii in layers){
            d[,paste(i, ii, sep = "_")] <- (d[,i]-d[,ii])/(d[,i]+d[,ii])
          }
        }
        
        d <- d %>% 
          dplyr::mutate(satid = ifelse(satid == "LT04", "LT05", satid))
        
        # Replace NA values with 0 and convert factor columns to numeric
        d <- d %>% 
          dplyr::mutate(dplyr::across(dplyr::everything(), ~ifelse(is.nan(.x), 0, .x)))
        
        std_names <- names(d)[grepl("_std", names(d))]
        d <- d %>% 
          dplyr::mutate(dplyr::across(dplyr::all_of(std_names), ~ifelse(is.na(.x), 0, .x)))
        
        # PREDICTIONS
        if(satid %in% c("LT05","LT04")){
          preds <- predict(mod5, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`5`, ~round(.x*100)))
        }
        if(satid == "LE07"){
          preds <- predict(mod7, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`4`, ~round(.x*100)))
        }
        if(satid %in% c("LC08","LC09")){
          preds <- predict(mod8, d, importance = "none", na.action = "na.omit")
          preds <- tibble(pred = preds$class) %>% 
            bind_cols(preds$predicted %>% as.data.frame()) %>% 
            mutate(across(`1`:`4`, ~round(.x*100)))
        }
        # d %>% select(-RADSAT) %>% filter(!complete.cases(.))
        rm(d)
        # ¨71§ 
        if(nrow(preds) == length(na.omit(values(r[["QA"]], mat = F)))){
          rr <- r[["QA"]]
          rr[["class"]] <- 0
          rr[["class"]][!is.na(r[["QA"]])] <- as.numeric(preds$pred)
          rr[["land"]] <- 0
          rr[["land"]][!is.na(r[["QA"]])] <- as.numeric(preds$`1`)
          rr[["water"]] <- 0
          rr[["water"]][!is.na(r[["QA"]])] <- as.numeric(preds$`2`)
          rr[["snow"]] <- 0
          rr[["snow"]][!is.na(r[["QA"]])] <- as.numeric(preds$`3`)
          rr[["cloud"]] <- 0
          rr[["cloud"]][!is.na(r[["QA"]])] <- as.numeric(preds$`4`)
          rr[["artif"]] <- 0
          if(satid %in% c("LT05","LT04")){
            rr[["artif"]][!is.na(r[["QA"]])] <- as.numeric(preds$`5`)
          }
          
          rr[is.na(r[["QA"]])] <- NA
          rr <- rr[[-1]]
          # plot(rr)
          # plotRGB(r, r=7, g=5, b=2, stretch = "lin")
          writeRaster(rr, paste0(class_landsat_dir, "/", imageid),
                      datatype = "INT1U", overwrite = T)
        } else {
          stop("Error: The number of predictions does not match the number of non-NA pixels in the QA band.")
        }
        
        return(imageid)
      }
    } else { return(NULL) }
  } else { return(NULL) }
}

# Internal helper function to calculate a mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
