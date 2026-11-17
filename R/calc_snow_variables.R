#' Calculate Snow Variables
#'
#' This function calculates snow variables from classified Landsat or Sentinel-2 imagery. It processes the
#' input image data frame, calculates snow proportions, and generates snow melting day (SCD) and trend rasters.
#'
#' The snow variables come from one sensor: Landsat ([classify_landsat()], 30 m) or Sentinel-2
#' ([classify_sentinel2()], 10 m). The two are not combined. The output has the resolution of the
#' classifications. Steps that depend on the pixel size use ground distances: clear pixels next to a
#' cloud (within 30 m) are left out, and the land cover is aggregated to the output resolution.
#' Sentinel-2 has 9 times as many pixels per area as Landsat and more images, so the calculation
#' takes correspondingly longer.
#'
#' @param image_df A tibble, the output from [extract_landsat()], [extract_sentinel2()] or
#' [search_image_df()], with the scenes of one sensor (see `sensor`).
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory of the imagery: the
#' Landsat base directory, or the Sentinel-2 one (`base_sentinel2_dir`) for Sentinel-2.
#' @param sensor `"landsat"` or `"sentinel2"`: the scenes used. NULL (default): taken from the file
#' names in `image_df`, which must then all be of one sensor.
#' @param method How the seasonal snow curve of each pixel is fitted (the model is the same):
#' \itemize{
#'   \item `"fast"` (default): all pixels of a chunk are fitted together with matrix operations
#'   (penalised IRLS), each pixel with its own data and its own smoothing parameter, chosen by the same
#'   REML criterion as mgcv over a grid of values (`lambdas`). The knots of the cyclic spline are
#'   placed per pixel as gam() places them (at quantiles of the pixel's own observation days), rounded
#'   to `knot_step` days so that pixels with similar observation days share a basis (see `knots`).
#'   Several times faster than `"gam"`; in tests across snow regimes and pixel-specific data gaps, melt
#'   and onset dates agree with `"gam"` within 0-1 days (median) at the same accuracy.
#'   \item `"gam"`: one `mgcv::gam()` per pixel, as in earlier versions.
#' }
#' @param missing_months Months without any image get pseudo-observations of snow (weight 0.8, the
#' number of images of an average month): `"low_sun"` (default) only months whose noon sun is below
#' 15 degrees at the landscape's latitude (polar night and low winter sun, when no surface reflectance
#' products are made); a month missing for other reasons (e.g. clouds of a wet season) is bridged by the
#' fitted curve instead of being taken as snowy. `"all"`: every month without images (earlier
#' versions); `"none"`: no pseudo-observations.
#' @param lambdas Grid of smoothing parameters for `method = "fast"`.
#' @param knots `method = "fast"`: `"pixel"` (default) places the spline knots of each pixel at
#' quantiles of its own observation days, as gam() does, so that pixels with gaps in other seasons than
#' the rest of the landscape (e.g. persistent clouds or shade) keep their flexibility where they have
#' data; `"landscape"` uses knots at quantiles of all image dates for every pixel (faster, but about 1
#' day less accurate in melt dates for pixels with such gaps).
#' @param knot_step Rounding of the pixel knots in days (default 10; 0 = exact, as gam()).
#' @param min_images Minimum number of images (after `max_cloud`) for the snow variables (default 50).
#' @param max_cloud Images classified as cloud in this proportion of the AOI or more are left out
#' (default 0.8).
#' @param min_prob Class probabilities (%) below this are not used as observations (default 20): a
#' pixel's snow observation of an image is used if its snow probability is at least `min_prob`, and
#' likewise its no-snow observation (land + water probability).
#' @param min_obs The snow curve is fitted for pixels with more than `min_obs` observation days over
#' all years (default 20).
#' @param threshold Snow probability of the fitted curve that defines snow-covered days
#' (`snow_days`), melt-out (`melt`) and snow onset (`ns`) (default 0.5).
#' @param window Half-width in days of the windows around melt-out and onset in which the observations
#' are counted (`melt_nobs`, `ns_nobs`; default 30).
#' @param workers A numeric value representing the number of workers for parallel processing. 
#' Setting workers > 1 (default) enables parallel computing across multiple nodes.
#' @param chunk_cells Approximate number of pixels processed at a time. Lower values
#' reduce memory use (memory scales with chunk_cells x number of images). NULL (default): about
#' 40 million / number of images (e.g. 50000 pixels for 800 images), at most 200000.
#' @return A SpatRaster object containing the snow variables, multi-year averages over all images
#' (days of year): `nobs` (observation days), `nyears`, `gamr2` (deviance explained), `max_snow_doy`,
#' `min_snow_doy`, `max_snow_prop`, `min_snow_prop` (30-day moving mean of the fitted curve),
#' `scd_raw` (weighted snow share of the observations x 365), `snow_days` (days with fitted
#' probability >= `threshold`), `scd` (sum of the fitted daily probabilities), `melt`, `ns`; their
#' uncertainty: `scd_se`, `melt_se`, `ns_se` (standard errors in days, from the Bayesian covariance of
#' the fitted curve, by the delta method; capped at 182.5 days, e.g. where the curve is nearly flat at
#' the threshold). They are conservative: in simulations the true average
#' melt-out day was within +-1.96 standard errors in 98-100% of the pixels, mainly because the
#' probability weights (< 1) count an image as less than one observation; use them as a relative
#' measure of reliability (they grow with fewer, noisier or badly timed observations) rather than as
#' calibrated confidence intervals. The observation density at the critical times:
#' `melt_nobs`, `ns_nobs` (observations within +-`window` days of melt-out and onset, over all
#' years); and `esalc`, `elevation`, `ndvi`. The settings (package version, sensor, method, image
#' period, thresholds) are stored with the result as metadata tags: see `terra::metags(result)`.
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
                                chunk_cells = NULL, sensor = NULL, method = c("fast", "gam"),
                                missing_months = c("low_sun", "all", "none"),
                                lambdas = 10^seq(-2, 6, by = 0.5), knots = c("pixel", "landscape"),
                                knot_step = 10, min_images = 50, max_cloud = 0.8, min_prob = 20,
                                min_obs = 20, threshold = 0.5, window = 30) {
  method <- match.arg(method)
  knots <- match.arg(knots)
  missing_months <- match.arg(missing_months)
  # Define directories for classifications and predictors
  class_landsat_dir <- paste0(base_landsat_dir, "/", site_name, "/classifications")
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  
  # List files in the classifications directory
  class_tifs <- list.files(class_landsat_dir, pattern = "GMT.tif$")
  
  # One sensor only: the one chosen, or the one of all scenes in image_df
  sens <- image_sensor(image_df$file)
  if (is.null(sensor)) {
    u <- unique(stats::na.omit(sens))
    if (length(u) > 1) {
      stop("image_df has both Landsat and Sentinel-2 scenes. Snow variables are calculated from one ",
           "sensor: choose it with sensor = \"landsat\" or sensor = \"sentinel2\".", call. = FALSE)
    }
    sensor <- if (length(u)) u else "landsat"
  } else {
    sensor <- match.arg(sensor, c("landsat", "sentinel2"))
    image_df <- image_df[sens %in% sensor, , drop = FALSE]
  }

  # Filter the image data frame to include only the classified images
  image_df <- image_df %>%
    filter(file %in% class_tifs)
  message("Snow variables from ", nrow(image_df), " classified ",
          c(landsat = "Landsat", sentinel2 = "Sentinel-2")[[sensor]], " images")
  
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
    filter(cloud_proportion_own < max_cloud)
  
  if (nrow(image_df) < min_images) {
    message("Less than ", min_images, " images with less than ", 100 * max_cloud,
            "% cloud! Snow variables not calculated.")
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
  
  # Pseudo-observations of snow for months without images (see missing_months): as many per month as
  # an average month has images (earlier versions drew them at random, with this expected number)
  miss <- setdiff(1:12, mcs$month)
  if (length(miss) && missing_months == "low_sun") {
    r0 <- rast(paste0(class_landsat_dir, "/", image_df$file[1]))
    lat <- sf::st_coordinates(sf::st_transform(sf::st_centroid(sf::st_as_sfc(sf::st_bbox(r0))), 4326))[1, "Y"]
    noon <- 90 - abs(lat - 23.44 * sin(2 * pi * (284 + middoys$middoy[miss]) / 365))
    if (any(noon >= 15)) {
      message("Month(s) ", paste(miss[noon >= 15], collapse = ", "), " without images are bridged by the ",
              "fitted curve (noon sun above 15 degrees: not assumed snowy)")
    }
    miss <- miss[noon < 15]
  }
  if (length(miss) && missing_months != "none") {
    extradf <- middoys %>%
      filter(month %in% miss) %>%
      mutate(doy = round(middoy), snow = 1, prop = 80) %>%
      select(snow, prop, doy, month)
    extradf <- extradf[rep(seq_len(nrow(extradf)), each = max(1, round(mean(mcs$n)))), , drop = FALSE]
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
  
  # Clear pixels within 30 m of a cloud are left out: a 3 x 3 window at 30 m (Landsat), 7 x 7 at 10 m
  # (Sentinel-2)
  px <- terra::res(rast(paste0(class_landsat_dir, "/", image_df$file[1])))[1]
  cloud_window <- 2 * max(1, round(30 / px)) + 1
  prob_files <- read_classifications(image_df, class_landsat_dir, tmp_dir, cloud_window, min_prob)
  
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
  
  # Fast method: the spline basis, its penalty and the 30-day moving mean, the same for all pixels
  basis <- if (method == "fast") snow_curve_basis(c(doys, extradf$doy)) else NULL
  if (is.null(chunk_cells)) chunk_cells <- min(200000, max(5000, floor(4e7 / length(doys))))
  
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
    
    if (method == "fast") {
      # blocks of pixels, one per worker; each pixel is fitted with its own data only
      nb <- max(1, min(workers * 2, length(cells) %/% 200 + 1))
      blk <- split(seq_along(cells), ceiling(seq_along(cells) * nb / length(cells)))
      jobs <- lapply(blk, function(ix) list(cells = cells[ix], snow = v_snow[ix, , drop = FALSE],
                                            nosnow = v_nosnow[ix, , drop = FALSE]))
      obs_chunk <- !is.na(v_snow) | !is.na(v_nosnow)
      rm(v_nosnow, v_snow)
      res <- if (!is.null(cl)) {
        parLapply(cl, jobs, snow_vars_job, doys = doys, years = years, extradf = extradf, basis = basis, lambdas = lambdas,
                  knots = knots, knot_step = knot_step, min_obs = min_obs, threshold = threshold)
      } else if (workers > 1) {
        mclapply(jobs, snow_vars_job, doys = doys, years = years, extradf = extradf, basis = basis,
                 lambdas = lambdas, knots = knots, knot_step = knot_step, min_obs = min_obs,
                 threshold = threshold, mc.cores = workers)
      } else {
        lapply(jobs, snow_vars_job, doys = doys, years = years, extradf = extradf, basis = basis, lambdas = lambdas,
               knots = knots, knot_step = knot_step, min_obs = min_obs, threshold = threshold)
      }
      failed <- vapply(res, function(x) !is.data.frame(x), logical(1))
      if (any(failed)) {
        stop("Snow variable calculation failed: ", paste(unique(unlist(lapply(res[failed], as.character))), collapse = "; "))
      }
      res <- bind_rows(res)
      res$melt_nobs <- snow_obs_near(obs_chunk, doys, res$melt, window)
      res$ns_nobs <- snow_obs_near(obs_chunk, doys, res$ns, window)
      res_list[[i]] <- res
      rm(jobs, res, obs_chunk)
      next
    }
    
    jobs <- lapply(seq_along(cells), function(j) {
      list(cell = cells[j], nosnow_prop = v_nosnow[j, ], snow_prop = v_snow[j, ])
    })
    obs_chunk <- !is.na(v_snow) | !is.na(v_nosnow)
    rm(v_nosnow, v_snow)
    
    if (!is.null(cl)) {
      res <- parLapply(cl, jobs, cal_scd, doys = doys, years = years, extradf = extradf,
                       min_obs = min_obs, threshold = threshold)
    } else if (workers > 1) {
      res <- mclapply(jobs, cal_scd, doys = doys, years = years, extradf = extradf,
                      min_obs = min_obs, threshold = threshold, mc.cores = workers)
    } else {
      res <- lapply(jobs, cal_scd, doys = doys, years = years, extradf = extradf,
                    min_obs = min_obs, threshold = threshold)
    }
    
    failed <- vapply(res, function(x) !is.data.frame(x), logical(1))
    if (any(failed)) {
      stop("Snow variable calculation failed for some pixels: ",
           paste(unique(unlist(lapply(res[failed], as.character))), collapse = "; "))
    }
    res <- bind_rows(res)
    res$melt_nobs <- snow_obs_near(obs_chunk, doys, res$melt, window)
    res$ns_nobs <- snow_obs_near(obs_chunk, doys, res$ns, window)
    res_list[[i]] <- res
    rm(jobs, res, obs_chunk)
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
  # ESALC.tif is on the grid of the images (Landsat: already the most common class of 3 x 3 pixels of
  # 10 m; Sentinel-2: 10 m). Aggregated only if it is finer than the output (earlier versions
  # aggregated the Landsat layer a second time, to 90 m)
  esalc <- rast(paste0(predictor_dir, "/ESALC.tif")) / 10
  fact <- round(terra::res(rs)[1] / terra::res(esalc)[1])
  if (fact > 1) esalc <- aggregate(esalc, fact, getmode)
  rs[["esalc"]] <- project(esalc, rs, method = "near")
  
  dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
  rs[["elevation"]] <- project(dem, rs)
  
  mi <- rast(paste0(predictor_dir, "/medianindices.tif"))
  rs[["ndvi"]] <- project(mi$ndvi, rs)
  
  # Bring the result into memory before the temporary files are deleted
  rs <- rs * 1
  # the settings, as metadata tags of the result
  tags <- c(snowman_version = as.character(utils::packageVersion("snowman")), sensor = sensor,
            method = method, knots = if (method == "fast") knots else "gam",
            images = nrow(image_df), first_image = as.character(min(img_dates)),
            last_image = as.character(max(img_dates)), min_images = min_images, max_cloud = max_cloud,
            min_prob = min_prob, min_obs = min_obs, threshold = threshold, window = window,
            missing_months = missing_months,
            pseudo_snow_months = if (is.null(extradf)) "" else paste(unique(extradf$month), collapse = ","))
  try(terra::metags(rs) <- tags, silent = TRUE)
  
  return(rs)
}

# Internal Function to Read Classifications.
# Writes a two-layer file (nosnow, snow probabilities) per image into tmp_dir
# and returns the file paths in the same order as imagedf$file.
read_classifications <- function(imagedf, basedir, tmp_dir, cloud_window = 3, min_prob = 20) {
  vapply(seq_along(imagedf$file), function(i) {
    x <- imagedf$file[[i]]
    r <- rast(paste0(basedir, "/", x))
    
    cloudm <- r[["class"]]
    cloudm[r[["class"]] > 3] <- NA
    cloudm <- focal(cloudm, cloud_window, min, na.rm = FALSE, fillvalue = 1, expand = FALSE)
    snow <- r[["snow"]]
    nosnow <- r[["land"]] + r[["water"]]
    snow <- mask(snow, cloudm)
    nosnow <- mask(nosnow, cloudm)
    snow[snow < min_prob] <- NA
    nosnow[nosnow < min_prob] <- NA
    
    rs <- c(nosnow, snow)
    names(rs) <- c("nosnow", "snow")
    
    out <- file.path(tmp_dir, sprintf("probs_%05d.tif", i))
    writeRaster(rs, out, datatype = "INT1U", overwrite = TRUE)
    return(out)
  }, character(1))
}

# Internal Function to Calculate Snow Melting Day and Trend for one pixel.
# x: list(cell, nosnow_prop, snow_prop), the last two are vectors over images
cal_scd <- function(x, doys, years, extradf, min_obs = 20, threshold = 0.5) {
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
    scd_raw = NA, snow_days = NA, scd = NA, melt = NA, ns = NA,
    scd_se = NA, melt_se = NA, ns_se = NA
  )
  
  if (results$nobs > min_obs) {
    
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
        results$snow_days <- round(sum(pred >= threshold),1)
        
        padded_pred <- c(pred[336:365], pred, pred[1:30])
        mv_mean <- stats::filter(padded_pred, rep(1/30, 30), sides = 2)
        mv_mean <- as.numeric(mv_mean[31:395]) # Extract the 1:365 part
        
        results$max_snow_doy <- which.max(mv_mean)
        results$min_snow_doy <- which.min(mv_mean)
        results$max_snow_prop <- round(max(mv_mean), 3)
        results$min_snow_prop <- round(min(mv_mean), 3)
        
        
        if (results$min_snow_prop < threshold & results$max_snow_prop > threshold) {
          
          p_expanded <- c(pred, pred)
          idx_max <- results$max_snow_doy
          idx_min <- results$min_snow_doy
          
          # Find melt and snow onset with a safety check
          search_melt <- p_expanded[idx_max:730] < threshold
          m_idx <- match(TRUE, search_melt)
          if(!is.na(m_idx)) results$melt <- m_idx + idx_max - 1
          
          search_ns <- p_expanded[idx_min:730] > threshold
          n_idx <- match(TRUE, search_ns)
          if(!is.na(n_idx)) results$ns <- n_idx + idx_min - 1
          
        }
        # uncertainty (Bayesian covariance of the fit, as for the fast method)
        Xp <- predict(gammod, data.frame(doy = 1:365), type = "lpmatrix")
        np <- ncol(Xp)
        se <- snow_curve_se(Xp, matrix(stats::coef(gammod), 1), array(gammod$Vp, c(1, np, np)),
                            results$melt, results$ns, matrix(pred, ncol = 1))
        results$scd_se <- se$scd_se; results$melt_se <- se$melt_se; results$ns_se <- se$ns_se
      }
    }
  }
  
  return(as.data.frame(results))
}

# ==========================================================================
# Fast method: the snow curve of cal_scd() for many pixels at once
# ==========================================================================
# The model of cal_scd(): binomial, cyclic cubic regression spline of the day of year (k = 5, from day 1
# to 365), observation weights = class probability / 100, and the pseudo-observations of missing months.
# The snow and no-snow rows of one image (weights s and n) have the same likelihood as one binomial
# observation y = s / (s + n) with weight s + n, so each image is one column. The smoothing parameter of
# each pixel is chosen by the Laplace-approximate REML criterion that mgcv uses for a known scale,
#   V(lambda) = D / 2 + lambda b'Sb / 2 + log|X'WX + lambda S| / 2 - rank(S) log(lambda) / 2,
# over a grid of lambda values (fits warm-started along the grid). Pixels are independent: every
# operation is per pixel (row); only the basis (from the image dates) is shared.

# Internal: basis (365 x 4), penalty and 30-day moving-mean matrix. Knots at quantiles of the days of the
# observations, as mgcv's cc smooth places them when only the end knots 1 and 365 are given.
snow_curve_basis <- function(doys) {
  snow_curve_basis_knots(mgcv::place.knots(c(1, 365, doys), 5))
}

# Internal: solve A x = b for many small symmetric positive definite systems at once (A: m x p x p,
# b: m x p), with the log-determinants of A (Gaussian elimination; the pivots of an SPD matrix are > 0)
solve_spd_many <- function(A, b) {
  p <- dim(A)[2]
  ld <- numeric(nrow(b))
  for (k in seq_len(p)) {
    piv <- A[, k, k]
    ld <- ld + log(pmax(piv, 1e-300))
    if (k < p) for (i in (k + 1):p) {
      f <- A[, i, k] / piv
      A[, i, k:p] <- A[, i, k:p, drop = FALSE] - f * A[, k, k:p, drop = FALSE]
      b[, i] <- b[, i] - f * b[, k]
    }
  }
  x <- matrix(0, nrow(b), p)
  for (k in p:1) {
    s <- b[, k]
    if (k < p) for (j in (k + 1):p) s <- s - A[, k, j] * x[, j]
    x[, k] <- s / A[, k, k]
  }
  list(x = x, logdet = ld)
}

# Internal: penalised IRLS for many pixels (W, Y: m x n weights and snow proportions; X: n x p basis rows
# of the observations), for each lambda of the grid; returns the coefficients at each pixel's best
# lambda and whether that fit converged
fit_snow_curves <- function(W, Y, X, S, rankS, lambdas, maxit = 30, tol = 1e-6) {
  m <- nrow(W); p <- ncol(X)
  XX <- do.call(cbind, lapply(seq_len(p), function(a) X[, a] * X))   # n x p^2
  sat <- rowSums(W * (ifelse(Y > 0, Y * log(Y), 0) + ifelse(Y < 1, (1 - Y) * log(1 - Y), 0)))
  beta <- matrix(0, m, p)
  beta[, 1] <- stats::qlogis(pmin(pmax(rowSums(W * Y) / rowSums(W), 0.01), 0.99))
  best_v <- rep(Inf, m); best_b <- beta; best_conv <- rep(FALSE, m); best_A <- array(0, c(m, p, p))
  for (lam in lambdas) {
    conv <- rep(FALSE, m)
    for (it in seq_len(maxit)) {
      eta <- pmin(pmax(beta %*% t(X), -30), 30)
      mu <- stats::plogis(eta)
      v <- pmax(mu * (1 - mu), 1e-10)
      ww <- W * v
      z <- eta + (Y - mu) / v
      A <- array(ww %*% XX, c(m, p, p)) + rep(lam * S, each = m)
      r <- solve_spd_many(A, (ww * z) %*% X)
      conv <- rowSums(abs(r$x - beta)) < tol * p
      beta <- r$x
      if (all(conv)) break
    }
    eta <- pmin(pmax(beta %*% t(X), -30), 30)
    mu <- stats::plogis(eta)
    dev <- -2 * (rowSums(W * (Y * log(pmax(mu, 1e-12)) + (1 - Y) * log(pmax(1 - mu, 1e-12)))) - sat)
    V <- dev / 2 + lam * rowSums((beta %*% S) * beta) / 2 + r$logdet / 2 - rankS * log(lam) / 2
    better <- is.finite(V) & V < best_v
    best_v[better] <- V[better]; best_b[better, ] <- beta[better, ]; best_conv[better] <- conv[better]
    best_A[better, , ] <- A[better, , , drop = FALSE]
  }
  # Bayesian covariance of the coefficients, (X'WX + lambda S)^-1, as mgcv's Vp for a known scale
  list(beta = best_b, conv = best_conv, Vb = inv_spd_many(best_A))
}

# Internal: one block of pixels (list of cells, snow, nosnow) for the parallel workers
snow_vars_job <- function(j, doys, years, extradf, basis, lambdas, knots = "pixel", knot_step = 10,
                          min_obs = 20, threshold = 0.5) {
  snow_vars_fast(j$cells, j$snow, j$nosnow, doys, years, extradf, basis, lambdas, knots, knot_step,
                 min_obs, threshold)
}

# Internal: the snow variables of cal_scd() for the pixels `cells` (snow, nosnow: pixels x images,
# class probabilities in %, NA = not used); the same columns and rules as cal_scd().
# knots = "pixel": as gam() in cal_scd(), the knots of each pixel are at quantiles of its own observation
# days (and the pseudo-observation days), here rounded to knot_step days so that pixels with similar
# observation days share a basis and are fitted together. knots = "landscape": basis$X for all pixels.
snow_vars_fast <- function(cells, snow, nosnow, doys, years, extradf, basis, lambdas,
                           knots = "pixel", knot_step = 10, min_obs = 20, threshold = 0.5) {
  npix <- nrow(snow)
  ws <- snow / 100; ws[is.na(ws)] <- 0
  wn <- nosnow / 100; wn[is.na(wn)] <- 0
  obs <- ws > 0 | wn > 0
  # number of observation days (an image date counts once) and years
  key <- doys + years * 1000
  if (anyDuplicated(key)) {
    uk <- unique(key)
    nobs <- rowSums((obs %*% outer(key, uk, "==")) > 0)
  } else nobs <- rowSums(obs)
  uy <- unique(years)
  nyears <- rowSums(matrix(vapply(uy, function(y) rowSums(obs[, years == y, drop = FALSE]) > 0,
                                  logical(npix)), nrow = npix))
  d <- doys
  if (!is.null(extradf) && nrow(extradf)) {
    ws <- cbind(matrix(extradf$prop / 100, npix, nrow(extradf), byrow = TRUE), ws)
    wn <- cbind(matrix(0, npix, nrow(extradf)), wn)
    d <- c(extradf$doy, d)
  }
  wt <- ws + wn
  y <- ifelse(wt > 0, ws / wt, 0)
  scd_raw <- round(rowSums(ws) / rowSums(wt) * 365, 1)
  out <- data.frame(cell = cells, nobs = nobs, nyears = nyears, gamr2 = NA_real_,
                    max_snow_doy = NA_real_, min_snow_doy = NA_real_, max_snow_prop = NA_real_,
                    min_snow_prop = NA_real_, scd_raw = ifelse(nobs > min_obs, scd_raw, NA_real_),
                    snow_days = NA_real_, scd = NA_real_, melt = NA_real_, ns = NA_real_,
                    scd_se = NA_real_, melt_se = NA_real_, ns_se = NA_real_)
  fit <- which(nobs > min_obs & !is.nan(scd_raw) & scd_raw / 365 > 0.02 & scd_raw / 365 < 0.98)
  if (!length(fit)) return(out)
  if (knots == "pixel") {
    kn <- snow_pixel_knots(wt[fit, , drop = FALSE] > 0, d, knot_step)
    gkey <- paste(kn[, 2], kn[, 3], kn[, 4], sep = "_")
  } else {
    gkey <- rep("landscape", length(fit))
  }
  for (g in unique(gkey)) {
    ix <- fit[gkey == g]
    B <- if (g == "landscape") basis else snow_basis_cached(kn[match(g, gkey), ], basis$M)
    r <- snow_fit_derive(wt[ix, , drop = FALSE], y[ix, , drop = FALSE], B$X[d, , drop = FALSE], B, lambdas,
                         threshold)
    out[ix, names(r)] <- r
  }
  out
}

# Internal: knots of each pixel as mgcv::place.knots(c(1, 365, observation days), 5) places them (at
# quantiles of the unique days), interior knots rounded to `step` days. used: pixels x observations
# (TRUE if the observation is used), d: day of each observation. Returns pixels x 5.
snow_pixel_knots <- function(used, d, step = 10) {
  dm <- matrix(0, length(d), 365); dm[cbind(seq_along(d), d)] <- 1
  dayobs <- (used %*% dm) > 0                       # pixels x 365: days with observations
  dayobs[, c(1, 365)] <- TRUE
  C <- dayobs %*% upper.tri(diag(365), diag = TRUE)   # cumulative number of days
  n <- C[, 365]
  kth <- function(k) max.col(C >= k, ties.method = "first")   # the day of the k-th unique day
  kn <- matrix(0, nrow(used), 5); kn[, 1] <- 1; kn[, 5] <- 365
  delta <- (n - 1) / 4
  for (j in 1:3) {
    pos <- delta * j + 1; lbi <- floor(pos); frac <- pos - lbi
    kn[, j + 1] <- kth(lbi) * (1 - frac) + kth(pmin(lbi + 1, n)) * frac
  }
  if (step > 0) kn[, 2:4] <- round(kn[, 2:4] / step) * step
  kn[, 2] <- pmin(pmax(kn[, 2], 2), 362)            # strictly increasing within (1, 365)
  kn[, 3] <- pmin(pmax(kn[, 3], kn[, 2] + 1), 363)
  kn[, 4] <- pmin(pmax(kn[, 4], kn[, 3] + 1), 364)
  kn
}

# Internal: snow_curve_basis_knots() kept per knot set within a process (pixels of later chunks often
# share the knots)
snow_basis_env <- new.env(parent = emptyenv())
snow_basis_cached <- function(kn, M) {
  k <- paste(kn, collapse = "_")
  if (is.null(snow_basis_env[[k]])) {
    if (length(ls(snow_basis_env)) > 5000) rm(list = ls(snow_basis_env), envir = snow_basis_env)
    b <- snow_curve_basis_knots(kn, M)
    assign(k, b, envir = snow_basis_env)
  }
  snow_basis_env[[k]]
}

# Internal: fit the pixels W, Y (one basis X) and derive the variables of cal_scd()
snow_fit_derive <- function(W, Y, X, basis, lambdas, threshold = 0.5) {
  f <- fit_snow_curves(W, Y, X, basis$S, basis$rankS, lambdas)
  ok <- f$conv                                    # (not converged: NA, as a failed gam())
  b <- f$beta
  m <- nrow(W)
  # deviance explained
  eta <- pmin(pmax(b %*% t(X), -30), 30); mu <- stats::plogis(eta)
  ll <- function(mu) rowSums(W * (Y * log(pmax(mu, 1e-12)) + (1 - Y) * log(pmax(1 - mu, 1e-12))))
  sat <- rowSums(W * (ifelse(Y > 0, Y * log(Y), 0) + ifelse(Y < 1, (1 - Y) * log(1 - Y), 0)))
  mu0 <- rowSums(W * Y) / rowSums(W)
  gamr2 <- 1 - (sat - ll(mu)) / (sat - ll(matrix(mu0, m, ncol(W))))
  # the curve over the year and the variables derived from it, as in cal_scd()
  pred <- stats::plogis(basis$X %*% t(b))          # 365 x m
  mv <- basis$M %*% pred
  maxd <- max.col(t(mv), ties.method = "first")
  mind <- max.col(t(-mv), ties.method = "first")
  maxp <- round(mv[cbind(maxd, seq_len(m))], 3)
  minp <- round(mv[cbind(mind, seq_len(m))], 3)
  first_from <- function(L, start) {               # first TRUE at or after row `start` (of 730)
    L[row(L) < rep(start, each = nrow(L))] <- FALSE
    r <- max.col(t(L) * 1, ties.method = "first")
    r[colSums(L) == 0] <- NA
    r
  }
  P2 <- rbind(pred, pred)
  has <- minp < threshold & maxp > threshold
  melt <- ifelse(ok & has, first_from(P2 < threshold, maxd), NA)
  ns <- ifelse(ok & has, first_from(P2 > threshold, mind), NA)
  se <- snow_curve_se(basis$X, b, f$Vb, melt, ns, pred)
  data.frame(gamr2 = ifelse(ok, gamr2, NA), max_snow_doy = ifelse(ok, maxd, NA),
             min_snow_doy = ifelse(ok, mind, NA), max_snow_prop = ifelse(ok, maxp, NA),
             min_snow_prop = ifelse(ok, minp, NA), snow_days = ifelse(ok, colSums(pred >= threshold), NA),
             scd = ifelse(ok, round(colSums(pred), 1), NA), melt = melt, ns = ns,
             scd_se = ifelse(ok, se$scd_se, NA), melt_se = se$melt_se, ns_se = se$ns_se)
}

# Internal: the basis of snow_curve_basis() for given knots (5 values from 1 to 365)
snow_curve_basis_knots <- function(kn, M = snow_moving_mean_matrix()) {
  G <- mgcv::gam(y ~ s(doy, bs = "cc", k = 5), knots = list(doy = kn),
                 data = data.frame(y = rep(0:1, length.out = 365), doy = 1:365), family = stats::binomial,
                 fit = FALSE)
  p <- ncol(G$X)
  S <- matrix(0, p, p)
  S[G$off[1]:p, G$off[1]:p] <- G$S[[1]]
  list(X = unname(G$X), S = S, rankS = qr(S)$rank, M = M, knots = kn)
}
snow_moving_mean_matrix <- function() {
  vapply(1:365, function(j) {
    e <- numeric(365); e[j] <- 1
    pp <- c(e[336:365], e, e[1:30])
    as.numeric(stats::filter(pp, rep(1 / 30, 30), sides = 2))[31:395]
  }, numeric(365))
}

# Internal: inverse of many small SPD matrices (m x p x p)
inv_spd_many <- function(A) {
  m <- dim(A)[1]; p <- dim(A)[2]
  V <- array(0, c(m, p, p))
  for (k in seq_len(p)) {
    e <- matrix(0, m, p); e[, k] <- 1
    V[, , k] <- solve_spd_many(A, e)$x
  }
  V
}

# Internal: standard errors of the snow cover duration and of the melt-out and onset days of fitted
# snow curves, by the delta method from the Bayesian covariance Vb of the coefficients.
# X365: 365 x p basis; b: m x p coefficients; Vb: m x p x p; melt, ns: days (1-730, or NA);
# pred: 365 x m fitted probabilities. scd = sum_d plogis(eta_d), so its gradient is
# sum_d p_d (1 - p_d) x_d; a crossing day t of the threshold has se(t) = se(eta(t)) / |eta'(t)|
# (capped at 182.5 days where the curve is nearly flat at the threshold; scd_se likewise, e.g. for
# gam() fits that diverge when snow and no-snow observations are perfectly separated).
snow_curve_se <- function(X365, b, Vb, melt, ns, pred) {
  m <- nrow(b); p <- ncol(b)
  qf <- function(G) {                               # sqrt(g' Vb g) for each pixel (row of G)
    s <- 0
    for (a in seq_len(p)) for (c in seq_len(p)) s <- s + G[, a] * Vb[, a, c] * G[, c]
    sqrt(pmax(s, 0))
  }
  G <- t(crossprod(X365, pred * (1 - pred)))        # m x p
  scd_se <- qf(G)
  day_se <- function(t) {
    out <- rep(NA_real_, m)
    k <- !is.na(t)
    if (!any(k)) return(out)
    t0 <- ifelse(k, ((t - 1) %% 365) + 1, 1)
    tp <- t0 %% 365 + 1; tm <- (t0 - 2) %% 365 + 1
    slope <- rowSums((X365[tp, , drop = FALSE] - X365[tm, , drop = FALSE]) * b) / 2
    se_eta <- qf(X365[t0, , drop = FALSE])
    out[k] <- pmin(se_eta[k] / pmax(abs(slope[k]), 1e-12), 182.5)
    out
  }
  data.frame(scd_se = round(pmin(scd_se, 182.5), 1), melt_se = round(day_se(melt), 1), ns_se = round(day_se(ns), 1))
}

# Internal: number of observations (pixel x image, over all years) within +-window days (cyclic) of a
# day of each pixel (1-730, or NA)
snow_obs_near <- function(obs, doys, day, window = 30) {
  out <- rep(NA_real_, nrow(obs))
  k <- which(!is.na(day))
  if (!length(k)) return(out)
  t <- ((day[k] - 1) %% 365) + 1
  D <- abs(outer(t, doys, "-"))
  D <- pmin(D, 365 - D)
  out[k] <- rowSums(obs[k, , drop = FALSE] & D <= window)
  out
}
