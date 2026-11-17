#' Classify Landsat Imagery
#'
#' This function classifies Landsat imagery using pretrained models. It processes the input image data frame,
#' applies the models, and generates classification rasters.
#'
#' @param image_df A tibble, the output from the `extract_landsat` or `search_image_df` functions.
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param model_dir Optional directory with your own `latest_TM05.rds`, `latest_LE07.rds` and/or
#' `latest_LC08.rds` (e.g. downloaded with [download_model()] or trained yourself). NULL (default), or a
#' model missing from the directory: the models included in snowman ([snowman_model_file()]).
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
classify_landsat <- function(image_df, site_name, base_landsat_dir, model_dir = NULL, workers = 1, force = FALSE) {
  # Define directories for predictors and classifications
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  class_landsat_dir <- paste0(base_landsat_dir, "/", site_name, "/classifications")
  
  # Create the classifications directory if it does not exist
  if (!dir.exists(class_landsat_dir)) {
    dir.create(class_landsat_dir)
  }
  
  tifs <- list.files(class_landsat_dir, pattern = ".tif$")
  
  if(!force){
    image_df <- image_df %>% 
      dplyr::filter(!file %in% tifs)
  }
  
  if(nrow(image_df) == 0 & !force){
    message("All imagery classified previously!")
    return(NULL)
  }
  
  
  # Load pretrained models
  if(any(c("LT05","LT04") %in% unique(image_df$satid))){
    mod5 <- s2_model_check(readRDS(snowman_model_file("TM05", model_dir)))
  } else {
    mod5 <- NULL
  }
  if(any(c("LE07") %in% unique(image_df$satid))){
    mod7 <- s2_model_check(readRDS(snowman_model_file("LE07", model_dir)))
  } else {
    mod7 <- NULL
  }
  if(any(c("LC08","LC09") %in% unique(image_df$satid))){
    mod8 <- s2_model_check(readRDS(snowman_model_file("LC08", model_dir)))
  } else {
    mod8 <- NULL
  }
  
  if(nrow(image_df) <= 10){
    message(paste0("Only ", nrow(image_df), " images to classify. Forced to sequential job!"))
    
    lss <- lapply(image_df$file, function(f) {
      safe_classify(
        imageid = f, 
        image_df = image_df,
        predictor_dir = predictor_dir, 
        class_landsat_dir = class_landsat_dir, 
        base_landsat_dir = base_landsat_dir, 
        site_name = site_name, 
        force = force,
        paral = FALSE,
        mod7 = mod7, mod5 = mod5, mod8 = mod8
      )}
    )
    
    return(summarise_classifications(lss))
    
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
    safe_classify(
      imageid = f, 
      image_df = image_df,
      predictor_dir = predictor_dir, 
      class_landsat_dir = class_landsat_dir, 
      base_landsat_dir = base_landsat_dir, 
      site_name = site_name, 
      force = force,
      paral = TRUE,
      mod7 = mod7, mod5 = mod5, mod8 = mod8
    )},
    future.packages = c("dplyr", "sf", "terra", "lubridate",
                        "solartime", "randomForestSRC",
                        if (any(vapply(list(mod5, mod7, mod8), function(m) isTRUE(m$type == "lightgbm"), logical(1)))) "lightgbm"),
    future.seed = TRUE,
    future.scheduling = 10
  )
  
  return(summarise_classifications(lss))
}

# Internal wrapper: an error in one image does not stop the whole run.
# Returns list(file, ok, msg).
safe_classify <- function(imageid, ...) {
  tryCatch({
    out <- classifying_function(imageid = imageid, ...)
    list(file = imageid, ok = !is.null(out), msg = NA_character_)
  }, error = function(e) {
    list(file = imageid, ok = FALSE, msg = conditionMessage(e))
  })
}

# Internal: collect results, warn about failed images, return classified ids
summarise_classifications <- function(lss) {
  failed <- Filter(function(x) !x$ok && !is.na(x$msg), lss)
  if (length(failed) > 0) {
    msgs <- vapply(failed, function(x) paste0(x$file, ": ", x$msg), character(1))
    warning(length(failed), " image(s) could not be classified:\n",
            paste(utils::head(msgs, 20), collapse = "\n"),
            if (length(msgs) > 20) "\n..." else "",
            call. = FALSE)
  }
  ok <- vapply(lss, function(x) isTRUE(x$ok), logical(1))
  unlist(lapply(lss[ok], `[[`, "file"))
}

# Internal Function to Classify a Single Landsat Image
classifying_function <- function(imageid, image_df, predictor_dir, class_landsat_dir,
                                 base_landsat_dir, site_name, force = FALSE,
                                 mod8, mod7, mod5, paral = TRUE) {
  # imageid <- image_df$file[[1]]
  # Check if the classified raster already exists
  if (!file.exists(paste0(class_landsat_dir, "/", imageid)) | force == TRUE) {
    
    # When distribute images across several cores, 
    # turn internal parallelization off
    if(paral){
      terra::terraOptions(threads = 1) 
      options(rf.cores = 1, mc.cores = 1)
    }
    
    # Get the satellite ID from the image metadata
    satid <- image_df %>% dplyr::filter(file == imageid) %>% dplyr::pull(satid)
    
    # Select the appropriate model variables based on the satellite ID
    if (satid %in% c("LT05", "LT04")) {
      mod_vars <- mod5$vars
    } else if (satid == "LE07") {
      mod_vars <- mod7$vars
    } else if (satid %in% c("LC08", "LC09")) {
      mod_vars <- mod8$vars
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

    # Only the predictors the model uses are calculated (each layer below is skipped otherwise)
    uses <- function(pattern) any(grepl(pattern, mod_vars, fixed = TRUE))

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
        dem <- terra::rast(paste0(predictor_dir, "/ALOSDEM.tif"))
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

      # Median Landsat differences, if used in the model
      extra <- list()
      if (uses("mr_diff")) {
        idate <- as.character(lubridate::ymd(stringr::str_split(imageid, "_")[[1]][4]))
        mo <- lubridate::month(idate)
        # Use the median image of the same month, or of the nearest available
        # month (circular distance, so December is next to January)
        mr_files <- list.files(predictor_dir, pattern = "^medianlandsat_[0-9]+\\.tif$")
        if (length(mr_files) == 0) {
          stop("No medianlandsat_*.tif files found. Run calc_predictors() first.")
        }
        avail_mo <- as.numeric(gsub("medianlandsat_|\\.tif", "", mr_files))
        mo_dist <- pmin(abs(avail_mo - mo), 12 - abs(avail_mo - mo))
        mo_use <- avail_mo[which.min(mo_dist)]
        mr <- terra::rast(paste0(predictor_dir, "/medianlandsat_", mo_use, ".tif"))

        # Median image minus the image, band by band in layer order, exactly as in the training data
        # of the models: for Landsat 4-7 the training images had no B1 layer, so the 7 median bands
        # (B1-B7) were paired with the image layers B2, B3, B4, B5, B6, B7, B10. The Landsat 4-7
        # models rely on these shifted differences (e.g. LE07 mr_diff_B7 = median B7 - B10).
        # Only the differences the model uses (directly) are calculated.
        img_layers <- if (satid %in% c("LC08", "LC09")) paste0("B", 1:7) else c("B2", "B3", "B4", "B5", "B6", "B7", "B10")
        md_names <- paste0("mr_diff_", c("B1", "B2", "B3", "B4", "B5", "B6", "B7"))
        for (k in which(md_names %in% mod_vars)) {
          x <- mr[[k]] - r[[img_layers[k]]]
          names(x) <- md_names[k]
          extra[[md_names[k]]] <- x
        }
        if (uses("mr_diff_B5_std5")) {
          x <- terra::focal(r[["B5"]], w = 5, fun = "sd") - terra::focal(mr[["B5"]], w = 5, fun = "sd")
          names(x) <- "mr_diff_B5_std5"
          extra[["mr_diff_B5_std5"]] <- x
        }
        if (uses("mr_diff_B7_std5")) {
          x <- terra::focal(r[["B7"]], w = 5, fun = "sd") - terra::focal(mr[["B7"]], w = 5, fun = "sd")
          names(x) <- "mr_diff_B7_std5"
          extra[["mr_diff_B7_std5"]] <- x
        }
      }

      # Median indices, if used in the model
      mi_names <- intersect(c("ndvi", "ndsi", "kbri", "swm", "bitm"), mod_vars)
      if (length(mi_names)) {
        mi <- terra::rast(paste0(predictor_dir, "/medianindices.tif"))
        extra[["mi"]] <- mi[[mi_names]]
      }

      # Focal statistics, if used in the model (directly, or in a pair such as B5_min5_B7_min5)
      filter <- matrix(1, nrow = 5, ncol = 5)
      filter[ceiling(length(filter) / 2)] <- 0
      vars_no_md <- gsub("mr_diff_B[57]_std5", "", mod_vars)   # (mr_diff_B5_std5 is not B5_std5)
      in_vars <- function(nm) any(grepl(paste0("(^|_)", nm, "(_|$)"), vars_no_md))
      for (b in c("B2", "B5", "B7")) {
        if (in_vars(paste0(b, "_tpi"))) {
          r[[paste0(b, "_tpi")]] <- r[[b]] - terra::focal(r[[b]], w = filter, fun = mean,
                                                          na.policy = "omit", na.rm = TRUE)
        }
      }
      for (st in c("min", "max", "sd")) {
        for (b in c("B2", "B5", "B7")) {
          nm <- paste0(b, "_", if (st == "sd") "std" else st, "5")
          if (in_vars(nm)) {
            r[[nm]] <- terra::focal(r[[b]], w = 5, fun = st, na.policy = "omit", na.rm = TRUE)
          }
        }
      }

      # Proportions of QA classes in the surroundings, if used in the model. The QA bits are read with
      # vectorised bit operations (the same values as the earlier per-pixel intToBits(), which took
      # ~1 s per layer and scene; pixels without data count as 0), and the cloud bit only once for the
      # 5- and 25-pixel proportions.
      probs_used <- intersect(c(paste0(c("cirrus", "cloud", "cshadow", "snow", "clear", "water"), "_probs5"),
                                "cloud_probs25"), mod_vars)
      if (length(probs_used)) {
        qa_vals <- as.integer(terra::values(r[["QA"]], mat = FALSE))
        qa_bit <- function(bit) {
          b <- bitwAnd(bitwShiftR(qa_vals, bit - 1L), 1L)
          b[is.na(b)] <- 0L
          x <- r[[1]]
          terra::values(x) <- b
          x
        }
        fm5 <- terra::focalMat(terra::rast(ncols = 180, nrows = 180, xmin = 0), 5, "circle")
        fm5[fm5 > 0] <- 1
        qa_bits <- c(cirrus = 3L, cloud = 4L, cshadow = 5L, snow = 6L, clear = 7L, water = 8L)
        cloud_bit <- NULL
        for (k in names(qa_bits)) {
          nm <- paste0(k, "_probs5")
          if (nm %in% mod_vars) {
            x <- qa_bit(qa_bits[[k]])
            if (k == "cloud") cloud_bit <- x
            r[[nm]] <- terra::focal(x, w = fm5, fun = "mean", na.policy = "omit", na.rm = TRUE)
          }
        }
        if ("cloud_probs25" %in% mod_vars) {
          fm25 <- terra::focalMat(terra::rast(ncols = 180, nrows = 180, xmin = 0), 25, "circle")
          fm25[fm25 > 0] <- 1
          if (is.null(cloud_bit)) cloud_bit <- qa_bit(qa_bits[["cloud"]])
          r[["cloud_probs25"]] <- terra::focal(cloud_bit, w = fm25, fun = "mean", na.policy = "omit", na.rm = TRUE)
        }
        rm(qa_vals, cloud_bit)
      }

      # Combine all rasters into a single raster stack
      if (length(extra)) r <- c(r, terra::rast(unname(extra)))
      r[["QA"]][is.na(r[["B1"]])] <- NA
      r[["QA"]][is.na(r[["B2"]])] <- NA
      r[["QA"]][is.na(r[["B3"]])] <- NA
      r[["QA"]][is.na(r[["B4"]])] <- NA
      r[["QA"]][is.na(r[["B5"]])] <- NA
      r[["QA"]][is.na(r[["B6"]])] <- NA
      r[["QA"]][is.na(r[["B7"]])] <- NA

      # The pixel table: QA, the bands and the calculated predictors only
      keep <- intersect(names(r), unique(c("QA", paste0("B", 1:7), mod_vars, pair_layers_landsat)))
      d <- as.data.frame(r[[keep]], na.rm = FALSE, cells = TRUE)
      d <- d[!is.na(d$QA), , drop = FALSE]
      for (v in grep("tpi$", names(d), value = TRUE)) d[[v]][is.na(d[[v]])] <- 0

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
        # Normalised differences: only the pairs the model uses. The pairs are named in the order of
        # pair_layers_landsat (the order the models were trained with), e.g. B5_B7, B5_min5_B7_min5.
        pr <- landsat_pairs(names(d))
        pr <- pr[pr$name %in% mod_vars, , drop = FALSE]
        for (k in seq_len(nrow(pr))) {
          d[[pr$name[k]]] <- (d[[pr$a[k]]] - d[[pr$b[k]]]) / (d[[pr$a[k]]] + d[[pr$b[k]]])
        }

        # NaN -> 0 and missing sd -> 0 (in the model's predictors), as in the training data
        for (v in intersect(mod_vars, names(d))) {
          if (is.numeric(d[[v]])) {
            d[[v]][is.nan(d[[v]])] <- 0
            if (grepl("_std", v)) d[[v]][is.na(d[[v]])] <- 0
          }
        }

        # PREDICTIONS
        # (random forest or LightGBM models: see s2_model(); predicted with s2_predict_probs())
        if (satid %in% c("LT05", "LT04")) {
          mod <- mod5
        } else if (satid == "LE07") {
          mod <- mod7
        } else {
          mod <- mod8
        }
        mod <- s2_model_load(mod, threads = if (paral) 1 else NULL)
        
        # Predict only pixels where all model predictors are available.
        # Pixels with any missing predictor are left NA in the output instead
        # of stopping the whole classification.
        pred_vars <- base::intersect(mod$vars, names(d))
        missing_vars <- setdiff(mod$vars, names(d))
        if (length(missing_vars) > 0) {
          stop("Predictors missing from the data: ", paste(missing_vars, collapse = ", "))
        }
        ok_rows <- stats::complete.cases(d[, pred_vars, drop = FALSE])
        d <- d[ok_rows, , drop = FALSE]
        
        out_names <- c("class", "land", "water", "snow", "cloud", "artif")
        vals <- matrix(NA_real_, nrow = terra::ncell(r), ncol = length(out_names))
        
        if (nrow(d) > 0) {
          prob <- s2_predict_probs(mod, d)
          if (nrow(prob) != nrow(d)) {
            stop("The number of predictions does not match the number of predicted pixels.")
          }
          # columns by class label (1 land, 2 water, 3 snow, 4 cloud, 5 artif); 0 for classes the
          # model does not have
          vals[d$cell, 1] <- as.numeric(colnames(prob)[max.col(prob, ties.method = "first")])
          vals[d$cell, 2:6] <- 0
          lab_col <- c("1" = 2, "2" = 3, "3" = 4, "4" = 5, "5" = 6)
          for (lab in base::intersect(colnames(prob), names(lab_col))) {
            vals[d$cell, lab_col[[lab]]] <- round(prob[, lab] * 100)
          }
        }
        rm(d)
        
        rr <- terra::rast(r[["QA"]], nlyrs = length(out_names), names = out_names)
        terra::values(rr) <- vals
        terra::writeRaster(rr, paste0(class_landsat_dir, "/", imageid),
                           datatype = "INT1U", overwrite = TRUE)
        
        return(imageid)
      }
    } else { return(NULL) }
  } else { return(NULL) }
}

# Internal: the layers of the Landsat band pairs (normalised differences), in the order the models were
# trained with: a pair of layers a and b (a before b) is named "a_b"
pair_layers_landsat <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7",
                         "B5_min5", "B7_min5", "B5_max5", "B7_max5", "B5_std5", "B7_std5",
                         "mr_diff_B5_std5", "mr_diff_B7_std5")
landsat_pairs <- function(nms) {
  layers <- pair_layers_landsat[pair_layers_landsat %in% nms]
  if (length(layers) < 2) return(data.frame(a = character(0), b = character(0), name = character(0)))
  cb <- utils::combn(layers, 2)
  data.frame(a = cb[1, ], b = cb[2, ], name = paste(cb[1, ], cb[2, ], sep = "_"), stringsAsFactors = FALSE)
}

# Internal helper function to calculate a mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
