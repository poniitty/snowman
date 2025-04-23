#' Calculate Snow Variables
#'
#' This function calculates snow variables from classified Landsat imagery. It processes the input image data frame,
#' calculates snow proportions, and generates snow melting day (SCD) and trend rasters.
#'
#' @param image_df A tibble, the output from the `extract_landsat` or `search_image_df` functions.
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' @return A SpatRaster object containing the snow variables.
#' @examples
#' calc_snow_variables(image_df, site_name = "ExampleSite", base_landsat_dir = "path/to/landsat")
#' @export
#' @import dplyr sf terra lubridate stringr parallel tibble mgcv purrr tidyr zoo
calc_snow_variables <- function(image_df, site_name, base_landsat_dir, workers) {
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
  
  # Filter out fully cloudy images
  image_df <- image_df %>%
    filter(cloud_proportion_own < 0.80)
  
  if (nrow(image_df) >= 50) {
    # Fit a generalized additive model (GAM) to predict snow proportion
    image_df %>%
      filter(clear_proportion_own > 0.8 & fill_proportion < 0.2) %>%
      mutate(doy = yday(date)) %>%
      gam(snow_proportion_own ~ s(doy, bs = "cc"), data = ., knots = list(doy = c(1, 365)), family = "binomial") -> gmod
    
    # Create a prediction data frame
    preds <- tibble(doy = 1:365)
    preds$preds <- predict(gmod, preds, type = "response")
    
    # Determine the day of the year with the maximum snow proportion
    maxsnowdoy <- preds$doy[which.max(preds$preds)]
    maxsnowdoy <- ifelse(maxsnowdoy > 250, 30, maxsnowdoy)
    maxsnowdoy <- ifelse(maxsnowdoy < (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min) + 14,
                         (image_df %>% mutate(doy = yday(date)) %>% pull(doy) %>% min) + 14, maxsnowdoy)
    minsnowdoy <- preds$doy[which.min(preds$preds)]
    
    # Read the classified rasters
    rs <- read_classifications(image_df, class_landsat_dir)
    
    # Raster list to data frame
    mm <- rs %>%
      map(nosnow_values) %>%
      reduce(bind_rows)
    
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
        library(dplyr)
        library(tidyr)
        library(mgcv)
        library(zoo)
      })
      
      results <- parLapply(cl, mm %>%
                          mutate(cell2 = cell) %>%
                          nest(data = -cell2) %>%
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
  }
  
  return(rs)
}

#' Internal Function to Read Classifications
#'
#' This function reads classified rasters and prepares them for further processing.
#'
#' @param imagedf A tibble containing the metadata of the classified images.
#' @param basedir A character string representing the base directory for the classified images.
#' @return A list of SpatRaster objects.
#' @keywords internal
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

#' Internal Function to Calculate Snow Melting Day and Trend
#'
#' This function calculates the snow melting day (SCD) and trend for a given cell.
#'
#' @param d_all A data frame containing snow and no-snow values for a given cell.
#' @param extradf An optional data frame containing extra days for prediction.
#' @return A data frame containing the calculated snow variables.
#' @keywords internal
cal_scd <- function(d_all, extradf) {
  # d_all <- mm %>% mutate(cell2 = cell) %>% nest(data = -cell2) %>% slice(1) %>% pull(data)
  # d_all <- d_all[[1]]
  cellid <- d_all$cell[1]
  
  d_all <- bind_rows(d_all %>% select(nosnow, nosnow_prop, doy, year, month, week) %>% drop_na() %>%
                       rename(snow = nosnow, prop = nosnow_prop),
                     d_all %>% select(snow, snow_prop, doy, year, month, week) %>% drop_na() %>%
                       rename(prop = snow_prop))
  
  results <- tibble(cell = cellid,
                    nobs = d_all %>% select(doy, year) %>% distinct() %>% nrow,
                    nyears = length(unique(d_all$year)),
                    gamr2 = as.numeric(NA),
                    max_snow_doy = as.numeric(NA), min_snow_doy = as.numeric(NA),
                    max_snow_prop = as.numeric(NA), min_snow_prop = as.numeric(NA),
                    scd_raw = as.numeric(NA), snow_days = as.numeric(NA),
                    scd = as.numeric(NA), melt = as.numeric(NA), ns = as.numeric(NA))
  
  
  
  if (nrow(d_all) > 20) {
    
    if (!is.null(extradf)) {
      d_all <- bind_rows(extradf, d_all)
    }
    
    results$scd_raw <- round(weighted.mean(d_all$snow, d_all$prop) * 365, 1)
    
    if ((results$scd_raw / 365) > 0.01 | (results$scd_raw / 365) < 0.99 | !is.nan(results$scd_raw)) {
      pred_grid <- expand.grid(doy = 1:365)
      
      e <- try({
        d_all %>%
          mutate(year = factor(year)) %>%
          gam(snow ~ s(doy, bs = "cc", k = 5),
              knots = list(doy = c(1, 365)),
              data = ., family = "binomial", weights = prop, method = "REML") -> gammod
      }, silent = TRUE)
      
      if (!inherits(e, "try-error")) {
        gamsum <- summary(gammod)
        pred_grid$pred <- predict(gammod, pred_grid, type = "response", se.fit = FALSE)
        # pred_grid$pred <- preds$fit
        # pred_grid$ci <- preds$se.fit * 1.96
        
        results$gamr2 <- gamsum$r.sq
        results$scd <- round(sum(pred_grid$pred), 1)
        results$snow_days <- round(sum(pred_grid$pred >= 0.5),1)
        
        pred_grid <- bind_rows(pred_grid, pred_grid) %>%
          mutate(doy = row_number())
        
        pred_grid <- pred_grid %>%
          mutate(mv_mean = rollmean(pred, k = 30, align = "center", na.pad = TRUE))
        
        max_doy_orig <- which.max(pred_grid$mv_mean)
        max_doy_orig <- ifelse(max_doy_orig > 365, max_doy_orig - 365, max_doy_orig)
        
        min_doy_orig <- which.min(pred_grid$mv_mean)
        min_doy_orig <- ifelse(min_doy_orig > 365, min_doy_orig - 365, min_doy_orig)
        
        results$max_snow_doy <- max_doy_orig
        results$min_snow_doy <- min_doy_orig
        results$max_snow_prop <- round(max(pred_grid$mv_mean, na.rm = TRUE), 3)
        results$min_snow_prop <- round(min(pred_grid$mv_mean, na.rm = TRUE), 3)
        
        if (min(pred_grid$pred) < 0.5 & max(pred_grid$pred) > 0.5) {
          pred_grid <- pred_grid %>%
            slice(max_doy_orig:nrow(.)) %>%
            mutate(doy = row_number())
          
          max_doy <- which.max(pred_grid$mv_mean)
          min_doy <- which.min(pred_grid$mv_mean)
          
          results$melt <- which.max(pred_grid$pred[max_doy:730] < 0.5) + max_doy - 1 + max_doy_orig - 1
          results$ns <- which.max(pred_grid$pred[min_doy:730] > 0.5) + min_doy - 1 + max_doy_orig - 1
        }
      }
    }
  }
  
  return(results)
}

#' Internal Function to Prepare No-Snow Values
#'
#' This function prepares no-snow values for further processing.
#'
#' @param x A SpatRaster object.
#' @return A data frame containing no-snow and snow values.
#' @keywords internal
nosnow_values <- function(x) {
  r <- x
  x[["nosnow"]][!is.na(x[["nosnow"]])] <- 0
  x[["snow"]][!is.na(x[["snow"]])] <- 1
  m <- values(c(x, r), dataframe = TRUE) %>%
    set_names(c("nosnow", "snow", "nosnow.1", "snow.1")) %>%
    rownames_to_column("cell") %>%
    mutate(cell = as.integer(cell)) %>%
    rename(nosnow_prop = nosnow.1,
           snow_prop = snow.1)
  return(m)
}
