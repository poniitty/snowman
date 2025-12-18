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
calc_snow_variables <- function(image_df, site_name, base_landsat_dir, workers = 1) {
  # Define directories for classifications and predictors
  class_landsat_dir <- paste0(base_landsat_dir, "/", site_name, "/classifications")
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  
  # List files in the classifications directory
  class_tifs <- list.files(class_landsat_dir, pattern = "GMT.tif$")
  
  # Filter the image data frame to include only the classified images
  image_df <- image_df %>%
    filter(file %in% class_tifs)
  
  ########################################################################
  # Calculate class stats
  
  if (nrow(image_df) > 0) {
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
  }
  
  # Filter out mostly cloudy images
  image_df <- image_df %>%
    filter(cloud_proportion_own < 0.80)
  
  if (nrow(image_df) >= 50) {
    # # Fit a generalized additive model (GAM) to predict snow proportion
    # image_df %>%
    #   filter(clear_proportion_own > 0.8 & fill_proportion < 0.2) %>%
    #   mutate(doy = yday(date)) %>%
    #   gam(snow_proportion_own ~ s(doy, bs = "cc"), data = ., knots = list(doy = c(1, 365)), family = "binomial") -> gmod
    # 
    # # Create a prediction data frame
    # preds <- tibble(doy = 1:365)
    # preds$preds <- predict(gmod, preds, type = "response")
    # 
    # # Determine the day of the year with the maximum snow proportion
    # maxsnowdoy <- preds$doy[which.max(preds$preds)]
    # maxsnowdoy <- ifelse(maxsnowdoy > 250, 30, maxsnowdoy)
    # maxsnowdoy <- ifelse(maxsnowdoy < (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min) + 14,
    #                      (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min) + 14, maxsnowdoy)
    # minsnowdoy <- preds$doy[which.min(preds$preds)]
    
    # Read the classified rasters
    rs <- read_classifications(image_df, class_landsat_dir)
    
    # Raster list to data frame
    mm <- rs %>%
      map(nosnow_values) %>%
      bind_rows
    
    rtemp <- rs[[1]][[1]]
    
    # Add day of the year and year to the data frame
    mm$doy <- rep(yday(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                  each = ncell(rs[[1]]))
    mm$year <- rep(year(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                   each = ncell(rs[[1]]))
    mm$week <- rep(week(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                   each = ncell(rs[[1]]))
    mm$month <- rep(month(ymd(unlist(lapply(names(rs), function(x) str_split(x, "_")[[1]][4])))),
                    each = ncell(rs[[1]]))
    
    mm <- mm %>%
      mutate(week = ifelse(week > 52, 52, week))
    
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
    
    # Clean up and return the results
    rm(rs)
    gc()
    
    # Calculate snow melting day and trend
    os <- Sys.info()["sysname"]
    
    options(show.error.messages = FALSE)
    if (os == "Windows") {
      # Use parLapply on Windows
      cl <- makeCluster(workers)
      on.exit(stopCluster(cl))
      
      # Load required packages on each worker node
      clusterEvalQ(cl, {
        library(mgcv)
      })
      
      results <- parLapply(cl, mm %>%
                          mutate(cell2 = cell) %>%
                          tidyr::nest(data = -cell2) %>%
                          pull(data),
                        cal_scd, extradf = extradf)
    } else {
      # Use mclapply on unix
      results <- mclapply(mm %>%
                            mutate(cell2 = cell) %>%
                            nest(data = -cell2) %>%
                            pull(data),
                          cal_scd, mc.cores = workers, extradf = extradf)
    }
    options(show.error.messages = TRUE)
    
    results <- bind_rows(results) %>%
      mutate(melt = ifelse(melt > 365, melt - 365, melt),
             ns = ifelse(ns > 365, ns - 365, ns))
    
    gc()
    
    # Combine the results into a SpatRaster object
    rs <- lapply(names(results)[-1], function(x) {
      rtemp[] <- results %>% pull(x)
      names(rtemp) <- x
      return(rtemp)
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
    
  } else {
    print("Less than 50 individual images! Snow variables not calculated.")
    rs <- NULL
  }
  
  return(rs)
}

# Internal Function to Read Classifications
read_classifications <- function(imagedf, basedir) {
  rs <- lapply(imagedf$file, function(x) {
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
    
    return(rs)
  })
  
  names(rs) <- imagedf$file
  return(rs)
}

# Internal Function to Calculate Snow Melting Day and Trend
cal_scd <- function(d_all, extradf) {
  # d_all <- mm %>% mutate(cell2 = cell) %>% nest(data = -cell2) %>% slice(1) %>% pull(data)
  # d_all <- d_all[[1]]
  
  d_clean <- rbind(setNames(d_all[,c("nosnow","nosnow_prop","doy","year")],c("snow","prop","doy","year")),
                   setNames(d_all[,c("snow","snow_prop","doy","year")],c("snow","prop","doy","year")))
  d_clean <- d_clean[!(is.na(d_clean$snow) | is.na(d_clean$prop)),]
  
  results <- list(
    cell = d_all$cell[1], nobs = length(unique(d_clean$doy + (d_clean$year * 1000))),
    nyears = length(unique(d_clean$year)), gamr2 = NA, 
    max_snow_doy = NA, min_snow_doy = NA, max_snow_prop = NA, min_snow_prop = NA,
    scd_raw = NA, snow_days = NA, scd = NA, melt = NA, ns = NA
  )
  
  if (results$nobs > 20) {
    
    if (!is.null(extradf)) {
      d_clean <- rbind(extradf[, c("snow", "prop", "doy")], d_clean[,-4])
    }
    
    results$scd_raw <- round(weighted.mean(d_clean$snow, d_clean$prop) * 365, 1)
    
    if ((results$scd_raw / 365) > 0.02 & (results$scd_raw / 365) < 0.98 & !is.nan(results$scd_raw)) {
      
      e <- try({
        gammod <- gam(snow ~ s(doy, bs = "cc", k = 5),
                      knots = list(doy = c(1, 365)),
                      data = d_clean, family = "binomial", weights = prop, method = "REML",
                      control = list(epsilon = 1e-04,  # Lower precision (default is 1e-07)
                                     maxit = 50))
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

# Internal Function to Prepare No-Snow Values
# x <- rs[[1]]
nosnow_values <- function(x) {
  r <- x
  x[["nosnow"]][!is.na(x[["nosnow"]])] <- 0
  x[["snow"]][!is.na(x[["snow"]])] <- 1
  v_mat <- values(c(x, r), mat = TRUE) 
  
  tibble::tibble(
    cell = 1:nrow(v_mat), # Assuming all rasters have same extent/cell count
    nosnow = v_mat[, 1], # nosnow from original
    snow = v_mat[, 2],    # snow from original
    nosnow_prop = v_mat[, 3], # nosnow from original
    snow_prop = v_mat[, 4]    # snow from original
  )
}
