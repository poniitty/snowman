#' Extracts Landsat Imagery
#'
#' This function extracts Landsat imagery that fulfills the specified arguments.
#' It processes the area of interest (AOI), transforms coordinates, and downloads
#' the relevant Landsat scenes.
#'
#' @param aoi An `sf` object or a list containing the area of interest.
#' @param site_name A character string representing the name of the site.
#' @param aoi_size A numeric value representing the size of the AOI in kilometers. Ignored if aoi is a POLYGON.
#' @param start_date A character string representing the start date for the imagery.
#' @param end_date A character string representing the end date for the imagery.
#' @param months A numeric vector representing the months to include (default is all months, so 1:12).
#' @param sats A character string representing the satellite to use (default is Landsat-8 only, so "LC08").
#' @param minclouds A numeric value representing the minimum cloud cover percentage (default is 50).
#' @param base_landsat_dir A character string representing the base directory for Landsat imagery.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' Setting workers > 1 (default) enables parallel computing.
#' @param data_source A character string representing the data source (currently only option is "rstac").
#' @param force Logical, if TRUE deletes all previously downloaded imagery.
#' @return A tibble containing the metadata of the downloaded Landsat scenes.
#' @export
#' @import rstac dplyr sf terra lubridate stringr parallel tibble readr ggplot2 future future.apply lwgeom
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
  
  # --- 1. PORTABLE PARALLEL SETUP ---
  # We only set a plan if workers > 1 AND the current plan is idle (sequential).
  # This allows workers to persist across multiple rounds in a loop.
  os_type <- Sys.info()['sysname']
  
  if (workers > 1 && inherits(future::plan(), "sequential")) {
    if (os_type %in% c("Linux", "Darwin")) {
      future::plan(future::multicore, workers = workers)
    } else {
      # Multisession for Windows
      future::plan(future::multisession, workers = workers)
    }
  }
  
  # --- 2. INPUT VALIDATION & GEOMETRY FIX ---
  utmall <- utm_zones
  
  suppressWarnings({
    sf::sf_use_s2(FALSE)
    
    if (inherits(aoi, "sf")) {
      aoi <- sf::st_make_valid(aoi) 
      
      if(st_geometry_type(aoi) == "POINT"){
        aoi_mid <- aoi %>% st_centroid() %>% st_transform(crs = 4326)
      } else if(st_geometry_type(aoi) == "POLYGON"){
        aoi_mid <- aoi %>% st_centroid() %>% st_transform(crs = 4326)
        aoi <- aoi %>% st_transform(crs = 4326)
      } else {
        stop("ERROR: AOI sf object must be POINT or POLYGON")
      }
    } else if (is.list(aoi)) {
      aoi_mid <- as_tibble(aoi) %>% mutate(name = site_name) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
      aoi <- aoi_mid
    } else {
      stop("ERROR: AOI must be an sf object or a list")
    }
  })
  
  # Coordinate Reference System logic
  utm <- utmall[aoi_mid,] 
  if (nrow(utm) == 0) utm <- utmall[st_nearest_feature(aoi_mid, utmall),]
  
  lat <- st_coordinates(aoi_mid)[,"Y"]
  utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0", utm$ZONE), utm$ZONE)
  epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
  
  if (st_geometry_type(aoi) == "POLYGON") {
    aoi <- aoi %>% st_transform(crs = epsg) %>% mutate(name = site_name)
  } else {
    aoi <- aoi %>% st_transform(crs = epsg) %>% st_buffer(aoi_size * 1000) %>% 
      st_bbox() %>% st_as_sfc() %>% st_as_sf() %>% mutate(name = site_name)
  }
  
  # Guardrail for large AOIs
  if(as.numeric(st_area(aoi)) > 1000e6) stop("ERROR: AOI > $1000km^2$. Process in blocks.")
  
  # --- 3. DIRECTORY & TERRA SETUP ---
  area_landsat_dir <- file.path(base_landsat_dir, site_name, "imagery")
  if (!dir.exists(area_landsat_dir)) dir.create(area_landsat_dir, recursive = TRUE)
  
  my_temp_dir <- file.path(base_landsat_dir, "temp_terra")
  if(!dir.exists(my_temp_dir)) dir.create(my_temp_dir, recursive = TRUE)
  terra::terraOptions(tempdir = my_temp_dir)
  
  if(force){
    unlink(list.files(area_landsat_dir, pattern = "GMT.tif$", full.names = TRUE))
    unlink(file.path(area_landsat_dir, c("lss.csv", "lss_final.csv")))
  }
  
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  tif_dates <- if(length(tifs) > 0) {
    tibble(satid = sapply(str_split(tifs, "_"), `[`, 1),
           date = as.character(ymd(sapply(str_split(tifs, "_"), `[`, 4))))
  } else {
    tibble(satid = NA, date = NA)
  }
  
  # --- 4. EXTRACTION & PARALLEL PROCESSING ---
  if (data_source == "rstac") {
    lss <- extract_landsat_stac(workers = workers, start_date = start_date, end_date = end_date, months = months,
                                aoi = aoi, site_name = site_name, epsg = epsg, excl_dates = tif_dates,
                                sats = sats, area_landsat_dir = area_landsat_dir, minclouds = minclouds)
  }
  gc()
  
  if (is.null(lss)) return(NULL)
  
  lss <- lss %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
  if (file.exists(file.path(area_landsat_dir, "lss.csv"))) {
    existing <- read_csv(file.path(area_landsat_dir, "lss.csv"), show_col_types = FALSE) %>%
      mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
    lss <- bind_rows(existing, lss) %>% distinct()
  }
  
  write_csv(lss, file.path(area_landsat_dir, "lss.csv"))
  
  # --- 5. CLEANUP & QUALITY CONTROL ---
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  lss <- full_join(tibble(file = tifs), lss, by = "file") %>% filter(!is.na(file))
  
  # Use future_lapply for distributed file checking
  img_remove <- future.apply::future_lapply(
    lss$file, check_raster, image_dir = area_landsat_dir,
    future.seed = TRUE, future.packages = "terra", future.scheduling = 5
  ) %>% unlist()
  
  if (length(img_remove) > 0) {
    unlink(file.path(area_landsat_dir, img_remove))
    lss <- lss %>% filter(!file %in% img_remove)
  }
  
  # --- 6. COVERAGE CALCULATION ---
  lccs <- future.apply::future_lapply(
    lss$file, calc_coverages, image_dir = area_landsat_dir,
    future.seed = TRUE, future.packages = c("terra", "dplyr"), future.scheduling = 5
  ) %>% bind_rows()
  
  lss <- full_join(lss %>% select(-ends_with("proportion")), lccs, by = "file")
  
  # Final filtering: remove low coverage scenes
  img_remove_final <- lss %>% filter(clear_proportion < 0.1) %>% pull(file)
  unlink(file.path(area_landsat_dir, img_remove_final))
  lss <- lss %>% filter(!file %in% img_remove_final)
  
  write_csv(lss, file.path(area_landsat_dir, "lss_final.csv"))
  return(lss)
}

# Internal Function to Create a Base URL
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
extract_landsat_stac <- function(aoi, epsg, excl_dates, site_name,
                                 sats = c("LT04", "LT05", "LE07", "LC08", "LC09"),
                                 start_date, end_date, months = 1:12,
                                 minclouds, area_landsat_dir, workers) {
  
  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  
  # Mapping satellites
  sat_map <- c("LT04"="landsat-4", "LT05"="landsat-5", "LE07"="landsat-7", "LC08"="landsat-8", "LC09"="landsat-9")
  sats2 <- sat_map[sats]
  
  # Search logic
  it_obj <- s_obj %>%
    stac_search(collections = "landsat-c2-l2",
                bbox = st_bbox(aoi %>% st_transform(4326)),
                datetime = paste0(start_date, "/", end_date),
                limit = 1000) %>%
    get_request()
  
  # Logic to handle > 1000 items pagination manually if needed
  if (length(it_obj$features) == 1000) {
    # (Existing pagination logic kept as is)
    dr <- tibble(date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day")) %>%
      mutate(gr = cut_number(row_number(.), 20)) %>%
      group_by(gr) %>% group_split()
    
    it_obj <- s_obj %>%
      stac_search(collections = "landsat-c2-l2",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  datetime = paste0(min(dr[[1]]$date), "/", max(dr[[1]]$date)),
                  limit = 1000) %>% get_request()
    
    for (i in 2:length(dr)) {
      it_obj2 <- s_obj %>%
        stac_search(collections = "landsat-c2-l2",
                    bbox = st_bbox(aoi %>% st_transform(4326)),
                    datetime = paste0(min(dr[[i]]$date), "/", max(dr[[i]]$date)),
                    limit = 1000) %>% get_request()
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
    
    if(nrow(itst) > 0){
      itst$area <- as.numeric(st_area(itst)) / (1000 * 1000)
      
      # Duplicate handling
      tt <- itst %>%
        mutate(date = as_date(ymd_hms(datetime))) %>%
        group_by(date, instruments, `landsat:wrs_path`, .add = TRUE) %>%
        mutate(n = n()) %>%
        arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()
      
      suppressMessages({
        suppressWarnings({
          sf::sf_use_s2(FALSE) # GEOMETRY FIX
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
          sf::sf_use_s2(TRUE)
        })
      })
      
      itst <- itst %>%
        filter(`landsat:scene_id` %in% (tt %>% pull(`landsat:scene_id`)))
      
      it_obj <- it_obj %>%
        items_filter(filter_fn = function(x) { x$properties$`landsat:scene_id` %in% (itst %>% pull(`landsat:scene_id`)) })
      
      # --- CALL PROCESSING IN PARALLEL (FUTURE) ---
      juuh <- process_features_in_parallel(it_obj, area_landsat_dir, aoi, workers, epsg = epsg)
      
      # Renaming platforms for compatibility
      itst$platform <- case_when(
        itst$platform == "landsat-4" ~ "LT04",
        itst$platform == "landsat-5" ~ "LT05",
        itst$platform == "landsat-7" ~ "LE07",
        itst$platform == "landsat-8" ~ "LC08",
        itst$platform == "landsat-9" ~ "LC09",
        TRUE ~ itst$platform
      )
      
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
  } else {
    print("No new Landsat scenes downloaded!")
    idsALL <- NULL
  }
  
  return(idsALL)
}

#' Internal Function: Parallel Feature Processing with Grid Alignment
process_features_in_parallel <- function(it_obj, area_landsat_dir, aoi, workers, epsg) {
  
  # 1. Create a Grid Template from one valid feature to ensure alignment
  # We pick the first feature that matches the target epsg
  idx <- which(sapply(it_obj$features, function(x) x$properties$`proj:epsg`) == epsg)
  if(length(idx) == 0) idx <- 1 else idx <- idx[1]
  
  item_temp <- it_obj$features[[idx]]
  tmp_url <- make_vsicurl_url(rstac::assets_url(item_temp)[[1]])
  
  # Create numerical template for workers to avoid passing full SpatRaster objects
  r_temp <- terra::crop(terra::rast(tmp_url), aoi %>% sf::st_transform(crs = epsg))
  exte_vals <- as.numeric(sf::st_bbox(r_temp))[c(1,3,2,4)]
  reso_val <- terra::res(r_temp)[1]
  
  # 2. Parallel Processing with future_lapply
  # Note: workers is ignored here because plan is set in extract_landsat
  future.apply::future_lapply(
    it_obj$features, 
    process_feature, 
    area_landsat_dir = area_landsat_dir, 
    aoi = aoi,
    exte = exte_vals,
    reso = reso_val,
    future.packages = c("sf", "terra", "stringr", "rstac"),
    future.seed = TRUE,
    future.scheduling = 10,
    future.globals = c("make_vsicurl_url")
  )
}

#' Internal Function: Single Feature Extraction with Retry Logic
process_feature <- function(ft, area_landsat_dir, aoi, exte, reso) {
  
  make_vsicurl_url <- function(base_url) {
    paste0("/vsicurl", "?pc_url_signing=yes", "&pc_collection=landsat-c2-l2", "&url=", base_url)
  }
  
  nm <- paste0(paste(str_split(ft$id, "_")[[1]][1:4], collapse = "_"), "_",
               gsub("-", "", substr(ft$properties$created, 1, 10)), "_",
               paste(str_split(ft$id, "_")[[1]][5:6], collapse = "_"), "_",
               gsub(":", "", substr(ft$properties$datetime, 12, 19)), "GMT.tif")
  
  output_file <- file.path(area_landsat_dir, nm)
  if (file.exists(output_file)) return(paste0("EXISTS:", nm))
  
  full_url <- make_vsicurl_url(rstac::assets_url(ft))
  max_retries <- 3
  
  for (iii in 1:max_retries) {
    e <- try({
      # Attempt to load bands, handling vsicurl errors by dropping bad bands
      ee <- try(r_remote <- terra::rast(full_url), silent = TRUE)
      
      if(inherits(ee, "try-error")){
        # Regex to find failing URL
        raw_err <- as.character(ee)
        extracted_url <- gsub("\\n","", stringr::str_extract(raw_err, "(https://.+?)\\n"))
        r_remote <- terra::rast(full_url[!full_url %in% make_vsicurl_url(extracted_url)])
      }
      
      # Alignment Step: Use Numeric Template to force exact Grid
      template <- terra::rast(extent = terra::ext(exte), resolution = reso, crs = sf::st_crs(aoi)$wkt)
      
      r_cropped <- terra::crop(r_remote, sf::st_transform(aoi, sf::st_crs(r_remote)))
      r_final <- terra::project(r_cropped, template, method = "near")
      
      # Metadata formatting
      names(r_final) <- sapply(str_split(names(r_final), "_"), function(x) paste(x[8:9], collapse = "_"))
      terra::writeRaster(r_final, output_file, datatype = "INT2U", overwrite = TRUE, NAflag = 0)
    }, silent = TRUE)
    
    if (!inherits(e, "try-error")) return(paste0("SUCCESS:", nm))
    if (iii == max_retries) return(paste0("ERROR:", nm, " - ", as.character(e)))
    Sys.sleep(iii * 2) 
  }
}

#' Quality Control: Reset Workers
#' Call this function after your loop rounds are completely finished.
#' @export
snowman_cleanup <- function() {
  if (!inherits(future::plan(), "sequential")) {
    future::plan(future::sequential)
    message("Parallel workers shut down and session reset to sequential.")
  }
}

# Internal Function to check if raster files work properly
check_raster <- function(image, image_dir){
  fpath <- file.path(image_dir, image)
  
  # 1. Attempt to open header
  r <- try(terra::rast(fpath), silent = TRUE)
  if(inherits(r, "try-error")) return(image)
  
  # 2. Force a read of the min/max values 
  # This verifies that the actual data blocks/tiles are readable
  chk <- try(terra::minmax(r), silent = TRUE)
  
  if(inherits(chk, "try-error")) {
    return(image)
  } else {
    return(NULL)
  }
}

# Internal Function to calculate cloud cover within the rasters
calc_coverages <- function(image, image_dir){
  require(terra)
  require(dplyr)
  
  fpath <- file.path(image_dir, image)
  
  res_tibble <- tryCatch({
    rs <- terra::rast(fpath)
    
    # Identify QA band
    qa_band_name <- names(rs)[grepl("_pixel", names(rs), ignore.case = TRUE)]
    if(length(qa_band_name) == 0) stop("QA band not found")
    
    # 1. Fill Proportion (Using C++ internal freq)
    # We check band 4 (often the first SR band in your stack) for 0 values
    fill_tab <- terra::freq(rs[[4]], value = 0)
    fill_count <- if(nrow(fill_tab) > 0) fill_tab[1, "count"] else 0
    n_pixels <- terra::ncell(rs)
    fill_prop_val <- fill_count / n_pixels
    
    # 2. Cloud Proportion (Optimized Bitwise)
    qa_vals <- terra::values(rs[[qa_band_name]], mat = FALSE)
    
    # Extract unique values to minimize bitwise computation
    u_qa <- unique(qa_vals)
    u_qa <- u_qa[!is.na(u_qa)]
    
    # Landsat 8-9 QA_PIXEL bit logic (High Confidence Cloud: Bits 8 & 9 set)
    # In R intToBits, bit 8 is index 9, bit 9 is index 10.
    cloud_values <- u_qa[sapply(u_qa, function(x) {
      bits <- intToBits(x)
      # Check high confidence cloud (bits 8 and 9 both 1)
      return(as.integer(bits[9]) == 1 && as.integer(bits[10]) == 1)
    })]
    
    # Count pixels that match cloud-flagged values
    cloud_count <- sum(qa_vals %in% cloud_values, na.rm = TRUE)
    cloud_prop_val <- cloud_count / n_pixels
    
    # Clean up large objects before returning
    rm(qa_vals, rs); gc()
    
    dplyr::tibble(
      file = image,
      fill_proportion = fill_prop_val,
      cloud_proportion = cloud_prop_val,
      clear_proportion = 1 - (fill_prop_val + cloud_prop_val)
    )
    
  }, error = function(e) {
    # Return NA so the parallel loop doesn't fail
    return(dplyr::tibble(file = image, 
                         fill_proportion = NA, cloud_proportion = NA, clear_proportion = NA))
  })
  
  return(res_tibble)
}
