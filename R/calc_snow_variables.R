#' Calculate Snow Variables
#'
#' This function calculates snow variables from classified Landsat imagery. It processes the input image data frame,
#' calculates snow proportions, and generates snow melting day (SCD) and trend rasters.
#'
#' @param image_df A tibble, the output from the `extract_landsat` or `search_image_df` functions.
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing. 
#' Setting workers > 1 (default) enables parallel computing across multiple nodes.
#' @param chunk_cells Approximate number of pixels processed at a time. Lower values
#' reduce memory use (memory scales with chunk_cells x number of images).
#' @return A SpatRaster object containing the snow variables.
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
#' @import dplyr sf terra lubridate stringr parallel tibble mgcv purrr tidyr zoo
calc_snow_variables <- function(image_df, site_name, base_landsat_dir, workers = 1,
                                chunk_cells = 20000) {
  # Define directories for classifications and predictors
  class_landsat_dir <- paste0(base_landsat_dir, "/", site_name, "/classifications")
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  
  # List files in the classifications directory
  class_tifs <- list.files(class_landsat_dir, pattern = "GMT.tif$")
  
  # Filter the image data frame to include only the classified images
  image_df <- image_df %>%
    filter(file %in% class_tifs)
  
  if (nrow(image_df) == 0) {
    message("No classified images found! Snow variables not calculated.")
    return(NULL)
  }
  
  ########################################################################
  # Calculate class stats
  
  # Initialize new columns for snow proportions
  image_df <- image_df %>%
    mutate(
      cloud_proportion_own = as.numeric(NA),
      clear_proportion_own = as.numeric(NA),
      land_proportion_own = as.numeric(NA),
      snow_proportion_own = as.numeric(NA)
    )
  
  # Calculate snow proportions for each image
  for (imageid in image_df %>% pull(file)) {
    r <- rast(paste0(class_landsat_dir, "/", imageid))
    vr2 <- values(r[[1]], mat = FALSE)
    
    image_df[image_df$file == imageid, "cloud_proportion_own"] <- round(mean(vr2 == 4, na.rm = TRUE), 3)
    image_df[image_df$file == imageid, "clear_proportion_own"] <- round(mean(vr2 < 4, na.rm = TRUE), 3)
    image_df[image_df$file == imageid, "land_proportion_own"] <- round(mean(vr2 == 1, na.rm = TRUE), 3)
    image_df[image_df$file == imageid, "snow_proportion_own"] <- round(mean(vr2 == 3, na.rm = TRUE), 3)
  }
  
  # Filter out mostly cloudy images
  image_df <- image_df %>%
    filter(cloud_proportion_own < 0.80)
  
  if (nrow(image_df) < 50) {
    print("Less than 50 individual images! Snow variables not calculated.")
    return(NULL)
  }
  
  # Create a data frame with the middle day of the year for each month
  middoys <- tibble(date = seq.Date(as.Date("2022-01-01"), as.Date("2022-12-31"), 1)) %>%
    mutate(month = month(date),
           doy = yday(date)) %>%
    group_by(month) %>%
    summarise(ndays = n(),
              middoy = mean(doy))
  
  # Count the number of images per month
  mcs <- image_df %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
    count
  
  # Create a data frame with extra days if needed
  if (nrow(mcs) != 12) {
    extradf <- middoys %>%
      filter(!month %in% mcs$month) %>%
      mutate(doy = round(middoy),
             snow = 1,
             prop = 80) %>%
      select(snow, prop, doy, month) %>%
      sample_n(size = mean(mcs$n) * nrow(.), replace = TRUE)
  } else {
    extradf <- NULL
  }
  
  ########################################################################
  # Write the per-image snow / no-snow probabilities to temporary files,
  # so that all images do not need to be held in memory at once
  tmp_dir <- file.path(base_landsat_dir, site_name, "tmp_snowvars")
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  
  prob_files <- read_classifications(image_df, class_landsat_dir, tmp_dir)
  
  st <- rast(prob_files)
  nosnow_st <- st[[seq(1, nlyr(st), by = 2)]]
  snow_st <- st[[seq(2, nlyr(st), by = 2)]]
  rtemp <- rast(st[[1]]) # empty template with the image geometry
  
  img_dates <- ymd(str_split_i(image_df$file, "_", 4))
  doys <- yday(img_dates)
  years <- year(img_dates)
  
  # Parallel backend
  os <- Sys.info()["sysname"]
  cl <- NULL
  if (os == "Windows" && workers > 1) {
    cl <- makeCluster(workers)
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
      library(mgcv)
    })
  }
  
  ########################################################################
  # Calculate snow melting day etc. in chunks of raster rows
  rows_per_chunk <- max(1, floor(chunk_cells / ncol(rtemp)))
  row_starts <- seq(1, nrow(rtemp), by = rows_per_chunk)
  res_list <- vector("list", length(row_starts))
  
  for (i in seq_along(row_starts)) {
    nr <- min(rows_per_chunk, nrow(rtemp) - row_starts[i] + 1)
    v_nosnow <- terra::values(nosnow_st, row = row_starts[i], nrows = nr, mat = TRUE)
    v_snow <- terra::values(snow_st, row = row_starts[i], nrows = nr, mat = TRUE)
    first_cell <- terra::cellFromRowCol(rtemp, row_starts[i], 1)
    cells <- first_cell + seq_len(nrow(v_nosnow)) - 1
    
    jobs <- lapply(seq_along(cells), function(j) {
      list(cell = cells[j], nosnow_prop = v_nosnow[j, ], snow_prop = v_snow[j, ])
    })
    rm(v_nosnow, v_snow)
    
    if (!is.null(cl)) {
      res <- parLapply(cl, jobs, cal_scd, doys = doys, years = years, extradf = extradf)
    } else if (workers > 1) {
      res <- mclapply(jobs, cal_scd, doys = doys, years = years, extradf = extradf,
                      mc.cores = workers)
    } else {
      res <- lapply(jobs, cal_scd, doys = doys, years = years, extradf = extradf)
    }
    
    failed <- vapply(res, function(x) !is.data.frame(x), logical(1))
    if (any(failed)) {
      stop("Snow variable calculation failed for some pixels: ",
           paste(unique(unlist(lapply(res[failed], as.character))), collapse = "; "))
    }
    res_list[[i]] <- bind_rows(res)
    rm(jobs, res)
  }
  
  results <- bind_rows(res_list) %>%
    arrange(cell) %>%
    mutate(melt = ifelse(melt > 365, melt - 365, melt),
           ns = ifelse(ns > 365, ns - 365, ns))
  
  gc()
  
  # Combine the results into a SpatRaster object
  rs <- lapply(names(results)[-1], function(x) {
    r1 <- rtemp
    r1 <- terra::setValues(r1, as.numeric(results[[x]]))
    names(r1) <- x
    return(r1)
  })
  
  rs <- rast(rs)
  
  # Add ESA Land Cover, DEM and NDVI
  esalc <- rast(paste0(predictor_dir, "/ESALC.tif")) / 10
  esalc <- aggregate(esalc, 3, getmode)
  rs[["esalc"]] <- project(esalc, rs, method = "near")
  
  dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
  rs[["elevation"]] <- project(dem, rs)
  
  mi <- rast(paste0(predictor_dir, "/medianindices.tif"))
  rs[["ndvi"]] <- project(mi$ndvi, rs)
  
  # Bring the result into memory before the temporary files are deleted
  rs <- rs * 1
  
  return(rs)
}

# Internal Function to Read Classifications.
# Writes a two-layer file (nosnow, snow probabilities) per image into tmp_dir
# and returns the file paths in the same order as imagedf$file.
read_classifications <- function(imagedf, basedir, tmp_dir) {
  vapply(seq_along(imagedf$file), function(i) {
    x <- imagedf$file[[i]]
    r <- rast(paste0(basedir, "/", x))
    
    cloudm <- r[["class"]]
    cloudm[r[["class"]] > 3] <- NA
    cloudm <- focal(cloudm, 3, min, na.rm = FALSE, fillvalue = 1, expand = FALSE)
    snow <- r[["snow"]]
    nosnow <- r[["land"]] + r[["water"]]
    snow <- mask(snow, cloudm)
    nosnow <- mask(nosnow, cloudm)
    snow[snow < 20] <- NA
    nosnow[nosnow < 20] <- NA
    
    rs <- c(nosnow, snow)
    names(rs) <- c("nosnow", "snow")
    
    out <- file.path(tmp_dir, sprintf("probs_%05d.tif", i))
    writeRaster(rs, out, datatype = "INT1U", overwrite = TRUE)
    return(out)
  }, character(1))
}

# Internal Function to Calculate Snow Melting Day and Trend for one pixel.
# x: list(cell, nosnow_prop, snow_prop), the last two are vectors over images
cal_scd <- function(x, doys, years, extradf) {
  ok_n <- !is.na(x$nosnow_prop)
  ok_s <- !is.na(x$snow_prop)
  
  d_clean <- data.frame(
    snow = c(rep(0, sum(ok_n)), rep(1, sum(ok_s))),
    prop = c(x$nosnow_prop[ok_n], x$snow_prop[ok_s]),
    doy = c(doys[ok_n], doys[ok_s]),
    year = c(years[ok_n], years[ok_s])
  )
  
  results <- list(
    cell = x$cell, nobs = length(unique(d_clean$doy + (d_clean$year * 1000))),
    nyears = length(unique(d_clean$year)), gamr2 = NA, 
    max_snow_doy = NA, min_snow_doy = NA, max_snow_prop = NA, min_snow_prop = NA,
    scd_raw = NA, snow_days = NA, scd = NA, melt = NA, ns = NA
  )
  
  if (results$nobs > 20) {
    
    d_clean <- d_clean[, c("snow", "prop", "doy")]
    if (!is.null(extradf)) {
      d_clean <- rbind(as.data.frame(extradf[, c("snow", "prop", "doy")]), d_clean)
    }
    
    results$scd_raw <- round(weighted.mean(d_clean$snow, d_clean$prop) * 365, 1)
    
    if ((results$scd_raw / 365) > 0.02 & (results$scd_raw / 365) < 0.98 & !is.nan(results$scd_raw)) {
      
      # Class probabilities (%) are used as observation weights on a 0-1
      # scale. With binomial family the weights are interpreted as numbers of
      # trials, so the raw 20-100 values would inflate the information content.
      d_clean$w <- d_clean$prop / 100
      
      e <- try({
        suppressWarnings({ # non-integer successes warning from fractional weights
          gammod <- gam(snow ~ s(doy, bs = "cc", k = 5),
                        knots = list(doy = c(1, 365)),
                        data = d_clean, family = "binomial", weights = w, method = "REML",
                        control = list(epsilon = 1e-04,  # Lower precision (default is 1e-07)
                                       maxit = 50))
        })
      }, silent = TRUE)
      
      if (!inherits(e, "try-error")) {
        
        pred <- predict(gammod, data.frame(doy = 1:365), type = "response", se.fit = FALSE)
        
        results$gamr2 <- 1 - (gammod$deviance / gammod$null.deviance)
        results$scd <- round(sum(pred), 1)
        results$snow_days <- round(sum(pred >= 0.5),1)
        
        padded_pred <- c(pred[336:365], pred, pred[1:30])
        mv_mean <- stats::filter(padded_pred, rep(1/30, 30), sides = 2)
        mv_mean <- as.numeric(mv_mean[31:395]) # Extract the 1:365 part
        
        results$max_snow_doy <- which.max(mv_mean)
        results$min_snow_doy <- which.min(mv_mean)
        results$max_snow_prop <- round(max(mv_mean), 3)
        results$min_snow_prop <- round(min(mv_mean), 3)
        
        
        if (results$min_snow_prop < 0.5 & results$max_snow_prop > 0.5) {
          
          p_expanded <- c(pred, pred)
          idx_max <- results$max_snow_doy
          idx_min <- results$min_snow_doy
          
          # Find melt and snow onset with a safety check
          search_melt <- p_expanded[idx_max:730] < 0.5
          m_idx <- match(TRUE, search_melt)
          if(!is.na(m_idx)) results$melt <- m_idx + idx_max - 1
          
          search_ns <- p_expanded[idx_min:730] > 0.5
          n_idx <- match(TRUE, search_ns)
          if(!is.na(n_idx)) results$ns <- n_idx + idx_min - 1
          
        }
      }
    }
  }
  
  return(as.data.frame(results))
}
