#' Calculate Predictors for Landsat Imagery
#'
#' This function calculates various predictors from Landsat imagery and other data sources.
#' It processes the input image data frame, downloads necessary data, and computes several
#' predictors such as DEM variables, median Landsat images, and median indices.
#'
#' @param image_df A tibble, the output from the `extract_landsat` or `search_image_df` functions.
#' @param site_name A character string representing the name of the site.
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param data_source A character string representing the data source (default is "rstac").
#' @return Raster files written to disk.
#' @examples
#' \dontrun{
#' calc_predictors(image_df, site_name = "ExampleSite", base_landsat_dir = "path/to/landsat")
#' }
#' @export
#' @import rstac dplyr sf terra lubridate stringr parallel tibble readr
calc_predictors <- function(image_df, site_name, base_landsat_dir, data_source = "rstac") {
  # Create the predictor directory if it does not exist
  predictor_dir <- paste0(base_landsat_dir, "/", site_name, "/predictors")
  if (!dir.exists(predictor_dir)) {
    dir.create(predictor_dir)
  }
  
  # Read the first raster file to get the EPSG code
  r <- rast(paste0(base_landsat_dir, "/", site_name, "/imagery/", image_df$file[1]))
  epsg <- st_crs(r)$epsg
  
  # Create an area of interest (AOI) buffer
  aoi <- st_bbox(r) %>%
    st_as_sfc() %>%
    st_buffer(1000) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_as_sf()
  
  if (data_source == "rstac") {
    # ALOS DEM from STAC
    s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
    it_obj <- s_obj %>%
      stac_search(collections = "alos-dem",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  limit = 1000) %>%
      get_request()
    
    juuh <- lapply(it_obj$features, function(ft) {
      full_url <- make_vsicurl_url_dem(assets_url(ft) %>% sort)
      full_url <- full_url[endsWith(full_url, "_DSM.tif")]
      file_names <- gsub("TIF$", "tif", basename(full_url))
      
      juuh <- lapply(seq_len(length(full_url)), function(nr) {
        e <- try({
          gdal_utils(
            "warp",
            source = full_url[[nr]],
            destination = paste0(predictor_dir, "/", file_names[[nr]]),
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
    })
    
    dems <- list.files(predictor_dir, pattern = "_DSM.tif", full.names = TRUE)
    dems <- lapply(dems, function(x) {
      dem <- rast(x)
      dem[dem == 0] <- NA
      return(dem)
    })
    
    dem <- sprc(dems)
    dem <- mosaic(dem)
    dem[dem < 0] <- 0
    dem[is.na(dem)] <- 0
    
    writeRaster(dem, paste0(predictor_dir, "/ALOSDEM.tif"), overwrite = TRUE)
    unlink(list.files(predictor_dir, pattern = "_DSM.tif", full.names = TRUE))
    
    # ESA WorldCover from STAC
    it_obj <- s_obj %>%
      stac_search(collections = "esa-worldcover",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  datetime = "2021-01-01/2021-12-31",
                  limit = 1000) %>%
      get_request()
    
    juuh <- lapply(it_obj$features, function(ft) {
      full_url <- make_vsicurl_url_esa(assets_url(ft) %>% sort)
      full_url <- full_url[endsWith(full_url, "_Map.tif")]
      file_names <- gsub("TIF$", "tif", basename(full_url))
      
      juuh <- lapply(seq_len(length(full_url)), function(nr) {
        suppressMessages({
          suppressWarnings({
            e <- try({
              gdal_utils(
                "warp",
                source = full_url[[nr]],
                destination = paste0(predictor_dir, "/", file_names[[nr]]),
                options = c(
                  "-t_srs", st_crs(aoi)$wkt,
                  "-te", st_bbox(aoi),
                  "-tr", c(10, 10)
                )
              )
            }, silent = TRUE)
          })
        })
        
        if (class(e)[[1]] == "try-error") {
          return(FALSE)
        } else {
          return(TRUE)
        }
      })
    })
    
    esas <- list.files(predictor_dir, pattern = "_Map.tif", full.names = TRUE)
    esas <- lapply(esas, function(x) {
      esa <- rast(x)
      esa[esa == 0] <- NA
      return(esa)
    })
    
    esa <- sprc(esas)
    esa <- mosaic(esa, fun = "max")
    esa[esa < 0] <- 0
    esa[is.na(esa)] <- 0
    
    suppressMessages({
      suppressWarnings({
        writeRaster(esa, paste0(predictor_dir, "/ESALC.tif"), overwrite = TRUE, datatype = "INT1U")
      })
    })
    unlink(list.files(predictor_dir, pattern = "_Map.tif", full.names = TRUE))
    
  }
  
  # DEM Variables
  dem <- rast(paste0(predictor_dir, "/ALOSDEM.tif"))
  
  # SLOPE
  slp <- terrain(dem, "slope", unit = "degrees")
  names(slp) <- "slope"
  writeRaster(round(slp * 100), paste0(predictor_dir, "/slope.tif"),
              filetype = "GTiff", overwrite = TRUE, datatype = "INT2U")
  
  # Median Landsat Images
  if (nrow(image_df) > 0) {
    imagedf2 <- image_df %>%
      filter(cloud_proportion < 0.8,
             tier == "T1") %>%
      mutate(month = month(date)) %>%
      relocate(month, .after = date)
    
    if (nrow(imagedf2) < 50) {
      imagedf2 <- image_df %>%
        filter(cloud_proportion < 0.8) %>%
        mutate(month = month(date)) %>%
        relocate(month, .after = date)
    }
    
    for (mo in unique(imagedf2$month) %>% sort()) {
      imagedf3 <- imagedf2 %>%
        filter(month == mo)
      
      if (nrow(imagedf3) > 30) {
        imagedf3 <- imagedf3 %>%
          arrange(desc(clear_proportion)) %>%
          slice_head(n = 30)
      }
      
      rr <- lapply(imagedf3$file, function(x) {
        sat_id <- imagedf3 %>% filter(file == x) %>% pull(satid)
        r1 <- rast(paste0(base_landsat_dir, "/", site_name, "/imagery/", x))
        
        if (sat_id %in% c("LC08", "LC09")) {
          names(r1) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
        } else {
          rr <- r1[[1]]
          rr[] <- NA
          r1 <- c(rr, r1)
          names(r1) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
        }
        
        cmask <- r1[[1]]
        cmask[] <- unlist(lapply(as.numeric(values(r1[["QA"]])), mask_cloud))
        
        r1 <- r1[[1:7]]
        r1[cmask == 1] <- NA
        r1[r1 == 0] <- NA
        
        return(r1)
      })
      
      rr <- sprc(rr)
      rr <- terra::mosaic(rr, fun = "median")
      
      for (ilayer in names(rr)) {
        if (mean(is.na(values(rr[[ilayer]], mat = FALSE))) == 1) {
          rr[[ilayer]][] <- 0
        }
      }
      
      if (sum(is.na(values(rr, mat = FALSE))) > 0) {
        rr <- focal(rr, 5, "mean", na.rm = TRUE, na.policy = "only")
      }
      
      writeRaster(round(rr), paste0(predictor_dir, "/medianlandsat_", mo, ".tif"),
                  datatype = "INT2U", overwrite = TRUE)
    }
    
    r2 <- rast(paste0(predictor_dir, "/medianlandsat_",
                      unique(imagedf2$month) %>% sort(), ".tif"))
    
    all <- rast()
    for (ilayer in unique(names(r2))) {
      r3 <- r2[[which(names(r2) == ilayer)]]
      r4 <- approximate(r3, rule = 2, NArule = 2)
      all <- c(all, r4)
    }
    
    nmo <- unique(imagedf2$month) %>% sort() %>% length()
    for (mo in unique(imagedf2$month) %>% sort()) {
      wmo <- which(unique(imagedf2$month) %>% sort() == mo)
      rr <- all[[seq(wmo, nmo * 7, nmo)]]
      
      for (ilayer in names(rr)) {
        if (mean(is.na(values(rr[[ilayer]], mat = FALSE))) == 1) {
          rr[[ilayer]][] <- 0
        }
      }
      
      if (sum(is.na(values(rr, mat = FALSE))) > 0) {
        repeat {
          rr <- focal(rr, 5, "mean", na.policy = "only", na.rm = TRUE)
          if (sum(is.na(values(rr, mat = FALSE))) == 0) {
            break
          }
        }
      }
      
      writeRaster(round(rr), paste0(predictor_dir, "/medianlandsat_", mo, ".tif"),
                  datatype = "INT2U", overwrite = TRUE)
    }
  }
  
  # Median Indices
  imagedf2 <- image_df %>%
    filter(fill_proportion < 0.8,
           cloud_proportion < 0.8,
           tier == "T1") %>%
    mutate(month = month(date)) %>%
    relocate(month, .after = date) %>%
    group_by(month) %>%
    arrange(desc(clear_proportion)) %>%
    slice_head(n = 3) %>%
    ungroup()
  
  if (nrow(imagedf2) < 15) {
    imagedf2 <- image_df %>%
      filter(fill_proportion < 0.8,
             cloud_proportion < 0.8) %>%
      mutate(month = month(date)) %>%
      relocate(month, .after = date) %>%
      group_by(month) %>%
      arrange(desc(clear_proportion)) %>%
      slice_head(n = 3) %>%
      ungroup()
  }
  
  rr <- lapply(imagedf2$file, function(x) {
    sat_id <- imagedf2 %>% filter(file == x) %>% pull(satid)
    r <- rast(paste0(base_landsat_dir, "/", site_name, "/imagery/", x))
    
    if (sat_id %in% c("LC08", "LC09")) {
      names(r) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
    } else {
      rr <- r[[1]]
      rr[] <- NA
      r <- c(rr, r)
      names(r) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "QA", "RADSAT")
    }
    r[r[["QA"]] == 0] <- NA
    r[[1:7]] <- r[[1:7]] * 2.75e-05 - 0.2
    
    cmask <- r[[1]]
    cmask[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_cloud))
    smask <- r[[1]]
    smask[] <- unlist(lapply(as.numeric(values(r[["QA"]])), mask_snow))
    
    ndvi <- (r[["B5"]] - r[["B4"]]) / (r[["B5"]] + r[["B4"]])
    ndvi[smask == 1] <- NA
    names(ndvi) <- "ndvi"
    
    ndsi <- (r[["B3"]] - r[["B6"]]) / (r[["B3"]] + r[["B6"]])
    names(ndsi) <- "ndsi"
    
    kbri <- (r[["B6"]] - r[["B5"]]) / (20 * sqrt((r[["B6"]] + r[["B5"]])))
    kbri[smask == 1] <- NA
    names(kbri) <- "kbri"
    
    swm <- (r[["B2"]] + r[["B3"]]) / (r[["B5"]] + r[["B6"]])
    swm[smask == 1] <- NA
    names(swm) <- "swm"
    
    bitm <- (((r[["B2"]] ** 2.0) + (r[["B3"]] ** 2.0) + (r[["B4"]] ** 2.0)) / 3.0) ** 0.5
    bitm[smask == 1] <- NA
    names(bitm) <- "bitm"
    
    r1 <- c(ndvi, ndsi, kbri, swm, bitm)
    
    return(r1)
  })
  
  ndvi <- lapply(rr, function(x) {
    x <- x[["ndvi"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  ndvi <- sprc(ndvi)
  ndvi <- terra::mosaic(ndvi, fun = "mean")
  
  if (sum(is.na(values(ndvi, mat = FALSE))) > 0) {
    repeat {
      ndvi <- focal(ndvi, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (sum(is.na(values(ndvi, mat = FALSE))) == 0) {
        break
      }
    }
  }
  
  kbri <- lapply(rr, function(x) {
    x <- x[["kbri"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  kbri <- sprc(kbri)
  kbri <- terra::mosaic(kbri, fun = "mean")
  if (sum(is.na(values(kbri, mat = FALSE))) > 0) {
    repeat {
      kbri <- focal(kbri, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (sum(is.na(values(kbri, mat = FALSE))) == 0) {
        break
      }
    }
  }
  
  ndsi <- lapply(rr, function(x) {
    x <- x[["ndsi"]]
    x[x > 1] <- 1
    x[x < -1] <- -1
    return(x)
  })
  ndsi <- sprc(ndsi)
  ndsi <- terra::mosaic(ndsi, fun = "mean")
  if (sum(is.na(values(ndsi, mat = FALSE))) > 0) {
    repeat {
      ndsi <- focal(ndsi, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (sum(is.na(values(ndsi, mat = FALSE))) == 0) {
        break
      }
    }
  }
  
  swm <- lapply(rr, function(x) {
    x <- x[["swm"]]
    x[x > 10] <- NA
    x[x < -10] <- -NA
    return(x)
  })
  swm <- sprc(swm)
  swm <- terra::mosaic(swm, fun = "mean")
  if (sum(is.na(values(swm, mat = FALSE))) > 0) {
    repeat {
      swm <- focal(swm, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (sum(is.na(values(swm, mat = FALSE))) == 0) {
        break
      }
    }
  }
  
  bitm <- lapply(rr, function(x) {
    x <- x[["bitm"]]
    return(x)
  })
  bitm <- sprc(bitm)
  bitm <- terra::mosaic(bitm, fun = "median")
  if (sum(is.na(values(bitm, mat = FALSE))) > 0) {
    repeat {
      bitm <- focal(bitm, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (sum(is.na(values(bitm, mat = FALSE))) == 0) {
        break
      }
    }
  }
  
  r1 <- c(ndvi, ndsi, kbri, swm, bitm)
  names(r1) <- c("ndvi", "ndsi", "kbri", "swm", "bitm")
  
  writeRaster(round(r1, 4), paste0(predictor_dir, "/medianindices.tif"),
              overwrite = TRUE)
  
  unlink(list.files(tempdir(), full.names = TRUE))
}


# Internal Function to Create a Base URL with Microsoft Planetary Computer
make_vsicurl_url_dem <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=alos-dem",
    "&url=",
    base_url
  )
}

# Internal Function to Create a Base URL with Microsoft Planetary Computer
make_vsicurl_url_esa <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=esa-worldcover",
    "&url=",
    base_url
  )
}

# Internal Functions to extract pixel information from the Landsat Quality Assesment (QA) layer
mask_fill <- function(x) {as.numeric(intToBits(x)[1])}
mask_dilcloud <- function(x) {as.numeric(intToBits(x)[2])}
mask_cloud <- function(x) {as.numeric(intToBits(x)[4])}
mask_cshadow <- function(x) {as.numeric(intToBits(x)[5])}
mask_snow <- function(x) {as.numeric(intToBits(x)[6])}
mask_clear <- function(x) {as.numeric(intToBits(x)[7])}
mask_water <- function(x) {as.numeric(intToBits(x)[8])}
conf_cloud <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[9:10]), collapse = ""))}
conf_cshadow <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[11:12]), collapse = ""))}
conf_snow <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[13:14]), collapse = ""))}
conf_cirrus <- function(x) {as.numeric(paste(as.numeric(intToBits(x)[15:16]), collapse = ""))}