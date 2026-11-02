#' Extracts Landsat Imagery
#'
#' This function extracts Landsat imagery that fulfills the specified arguments.
#' It processes the area of interest (AOI), transforms coordinates, and downloads
#' the relevant Landsat scenes.
#'
#' Downloads are robust to temporary throttling by the Planetary Computer: GDAL
#' retries failed HTTP requests (`http_retries`), each scene is attempted up to
#' `max_retries` times with growing waits, and scenes that still fail are retried
#' in up to `retry_rounds` further rounds. Every scene's outcome is appended to
#' `imagery/download_log.csv`; scenes that could not be downloaded are reported
#' in a warning and attempted again when the function is run again.
#'
#' The function can be re-run to update a site (e.g. with a later `end_date`): it
#' downloads only scenes that are not yet known. Known scenes are images on disk,
#' scenes listed in `imagery/lss_final.csv` (also if their images have since been
#' deleted, e.g. after classification) and scenes rejected earlier as too cloudy
#' (`imagery/lss_rejected.csv`). Only the newly downloaded images are checked and
#' added to `lss_final.csv`.
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
#' @param force Logical, if TRUE deletes all previously downloaded imagery and the scene tables.
#' @param max_retries Attempts per scene within one round (default 5), with waits of about 5, 10, 20, 40 s.
#' @param retry_rounds Further rounds for scenes that failed in all their attempts (default 2).
#' @param retry_wait Pause in seconds before retry round k is `k * retry_wait` (default 60).
#' @param http_retries,http_retry_delay GDAL's own retries of failed HTTP requests and the initial
#' delay in seconds (GDAL_HTTP_MAX_RETRY, GDAL_HTTP_RETRY_DELAY).
#' @return A tibble containing the metadata of all usable Landsat scenes of the site (the content
#' of `lss_final.csv`), or NULL if there are none. Scenes that could not be downloaded in this run
#' are in `attr(result, "failed")`.
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
                            force = FALSE,
                            max_retries = 5,
                            retry_rounds = 2,
                            retry_wait = 60,
                            http_retries = 10,
                            http_retry_delay = 5) {

  # --- 1. PORTABLE PARALLEL SETUP (Final, Stable Version) ---
  # Check if running inside RStudio/RStudio Server (where multicore is unstable)
  # We use parallelly::supportsMulticore() == FALSE on Linux/Darwin to detect RStudio Server.
  is_rstudio_env <- parallelly::supportsMulticore() == FALSE && Sys.info()['sysname'] %in% c("Linux", "Darwin")

  # Determine how many workers are currently running in the active plan.
  current_plan_workers <- future::nbrOfWorkers()

  if (workers > 1 && current_plan_workers <= 1) {
    oplan <- future::plan()
    os_type <- Sys.info()['sysname']
    if (os_type == "Windows" || is_rstudio_env) {
      future::plan(future::multisession, workers = workers)
      message(paste("Setting parallel plan to 'multisession' (Sockets) with", workers, "workers for cross-platform stability."))
    } else if (os_type %in% c("Linux", "Darwin")) {
      future::plan(future::multicore, workers = workers)
      message(paste("Setting parallel plan to 'multicore' (Forking) with", workers, "workers."))
    }
    # Note: We rely on the user calling snowman_cleanup() to restore 'oplan'.
  } else if (workers > 1 && current_plan_workers > 1) {
    message(paste("Reusing existing parallel plan with", current_plan_workers, "workers."))
  } else {
    future::plan(future::sequential)
  }
  # --- END OF SETUP ---

  # Restore the user's global sf and terra settings when the function exits
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  old_terra_tmp <- tryCatch(terra::terraOptions(print = FALSE)$tempdir, error = function(e) NULL)
  if (!is.null(old_terra_tmp)) {
    on.exit(terra::terraOptions(tempdir = old_terra_tmp), add = TRUE)
  }

  # --- 2. INPUT VALIDATION & GEOMETRY FIX ---
  utmall <- utm_zones

  suppressWarnings({
    suppressMessages({
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
  })

  # Coordinate Reference System logic
  suppressMessages({
    utm <- utmall[aoi_mid,]
  })
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
  if(as.numeric(st_area(aoi)) > 1000e6) stop("ERROR: AOI > $1000km^2. Process in blocks.")

  # --- 3. DIRECTORY & TERRA SETUP ---
  area_landsat_dir <- file.path(base_landsat_dir, site_name, "imagery")
  if (!dir.exists(area_landsat_dir)) dir.create(area_landsat_dir, recursive = TRUE)

  my_temp_dir <- file.path(base_landsat_dir, "temp_terra")
  if(!dir.exists(my_temp_dir)) dir.create(my_temp_dir, recursive = TRUE)
  terra::terraOptions(tempdir = my_temp_dir)

  if(force){
    unlink(list.files(area_landsat_dir, pattern = "GMT.tif$", full.names = TRUE))
    unlink(file.path(area_landsat_dir, c("lss.csv", "lss_final.csv", "lss_rejected.csv")))
  }
  unlink(list.files(area_landsat_dir, pattern = "\\.part$", full.names = TRUE)) # leftovers of interrupted runs

  # Scenes that need no download (see the function description), and the scene
  # table of earlier runs, to which new scenes are added
  known <- known_groups(area_landsat_dir, "lss_final.csv", "lss_rejected.csv", sensor = "landsat")
  final_file <- file.path(area_landsat_dir, "lss_final.csv")
  prev_final <- if (file.exists(final_file)) {
    readr::read_csv(final_file, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else NULL

  # --- 4. SEARCH & DOWNLOAD ---
  if (data_source == "rstac") {
    dl <- extract_landsat_stac(workers = workers, start_date = start_date, end_date = end_date, months = months,
                               aoi = aoi, site_name = site_name, epsg = epsg, excl_groups = known,
                               sats = sats, area_landsat_dir = area_landsat_dir, minclouds = minclouds,
                               max_retries = max_retries, retry_rounds = retry_rounds, retry_wait = retry_wait,
                               http_retries = http_retries, http_retry_delay = http_retry_delay)
  }
  gc()
  write_download_log(dl$log, area_landsat_dir)
  report_downloads(dl$log, length(known), area_landsat_dir)
  failed <- if (is.null(dl$log)) NULL else dl$log %>% filter(status %in% c("ERROR", "BROKEN"))

  # lss.csv: catalogue metadata of every downloaded scene (also of those later
  # rejected), accumulated over runs
  meta_cols <- c("id", "file", "DATE_ACQUIRED", "CLOUD_COVER", "SUN_ELEVATION")
  meta <- if (!is.null(dl$meta)) dl$meta %>% mutate(across(everything(), as.character)) else NULL
  if (file.exists(file.path(area_landsat_dir, "lss.csv"))) {
    old_meta <- readr::read_csv(file.path(area_landsat_dir, "lss.csv"), show_col_types = FALSE,
                                col_types = readr::cols(.default = "c"))
    meta <- bind_rows(old_meta[, intersect(meta_cols, names(old_meta))], meta)
  }
  if (!is.null(meta) && nrow(meta) > 0) {
    meta <- meta %>% distinct(file, .keep_all = TRUE)
    write_csv(meta, file.path(area_landsat_dir, "lss.csv"))
  }

  # --- 5. NEW IMAGES: on disk but not yet in the final scene table ---
  tifs <- list.files(area_landsat_dir, pattern = "GMT.tif$")
  new_files <- setdiff(tifs, prev_final$file)
  if (length(new_files) == 0) {
    if (is.null(prev_final)) {
      print("Zero Landsat scenes found or downloaded!")
      return(NULL)
    }
    out <- readr::type_convert(prev_final, col_types = readr::cols())
    attr(out, "failed") <- failed
    return(out)
  }

  lss <- tibble(file = new_files)
  if (!is.null(meta)) lss <- left_join(lss, meta[, intersect(meta_cols, names(meta))], by = "file")
  lss <- parse_landsat_files(lss, site_name)

  # --- 6. CHECK RASTERS (FUTURE APPLY) ---
  # Unreadable images are deleted; they are not recorded as rejected, so they are
  # downloaded again on the next run.
  img_remove <- future.apply::future_lapply(
    lss$file,
    check_raster,
    image_dir = area_landsat_dir,
    future.seed = TRUE,
    future.packages = c("terra"),
    future.scheduling = 5
  ) %>% unlist()
  gc()
  if (length(img_remove) > 0) {
    message(length(img_remove), " unreadable image(s) deleted; they are downloaded again on the next run")
    unlink(file.path(area_landsat_dir, img_remove))
    lss <- lss %>% filter(!file %in% img_remove)
  }

  # --- 7. MOSAICKING (Sequential due to I/O intensity) ---
  lss_d <- lss %>% group_by(date, satid, path) %>% count %>% filter(n > 1) %>% ungroup()
  if (nrow(lss_d) > 0) {
    for (ii in seq_len(nrow(lss_d))) {
      lss_dd <- lss_d %>% slice(ii)
      lss_dd <- right_join(lss, lss_dd, by = join_by(satid, path, date))

      rs <- lapply(lss_dd$file, function(x) {
        r <- rast(file.path(area_landsat_dir, x))
        r[r == 0] <- NA
        return(r)
      })

      rs <- sprc(rs)
      rs <- mosaic(rs) %>% round()
      writeRaster(rs, file.path(area_landsat_dir, lss_dd$file[[1]]),
                  overwrite = TRUE, datatype = "INT2U")
      unlink(file.path(area_landsat_dir, lss_dd$file[2:nrow(lss_dd)]))
      lss <- lss %>% filter(!file %in% lss_dd$file[2:nrow(lss_dd)])
    }
  }

  # --- 8. CALC COVERAGES (FUTURE APPLY) ---
  lccs <- future.apply::future_lapply(
    lss$file,
    calc_coverages,
    image_dir = area_landsat_dir,
    future.seed = TRUE,
    future.packages = c("terra", "dplyr"),
    future.scheduling = 5
  ) %>% bind_rows()
  gc()

  lss <- full_join(lss %>% select(-ends_with("proportion")), lccs, by = "file")

  # Mostly cloudy or empty images are deleted and recorded, so they are not downloaded again
  rejected <- lss %>% filter(!is.na(clear_proportion), clear_proportion < 0.1)
  if (nrow(rejected) > 0) {
    unlink(file.path(area_landsat_dir, rejected$file))
    rej_file <- file.path(area_landsat_dir, "lss_rejected.csv")
    write_csv(rejected %>% select(file, date, satid, fill_proportion, cloud_proportion, clear_proportion) %>%
                mutate(across(everything(), as.character)),
              rej_file, append = file.exists(rej_file))
  }
  lss <- lss %>% filter(!file %in% rejected$file)

  # --- 9. FINAL SCENE TABLE: earlier scenes + the new ones ---
  final <- bind_rows(prev_final, lss %>% mutate(across(everything(), as.character))) %>%
    distinct(file, .keep_all = TRUE)
  final <- readr::type_convert(final, col_types = readr::cols()) %>% arrange(date)
  write_csv(final, final_file)
  message(nrow(lss), " new scene(s) added; ", nrow(final), " scenes in ", final_file)

  attr(final, "failed") <- failed
  return(final)
}

# Internal Function: scene metadata from the image file names
# (<satid>_<correction>_<pathrow>_<date>_<processing date>_<collection>_<tier>_<time>GMT.tif)
parse_landsat_files <- function(lss, site_name) {
  lss$area <- site_name
  lss_parts <- str_split(lss$file, "_")
  lss$collection <- gsub("0", "C", sapply(lss_parts, `[`, 6))
  lss$tier <- sapply(lss_parts, `[`, 7)
  lss$satid <- sapply(lss_parts, `[`, 1)
  lss$path <- as.numeric(substr(sapply(lss_parts, `[`, 3), 1, 3))
  lss$row <- as.numeric(substr(sapply(lss_parts, `[`, 3), 4, 6))
  lss$date <- ymd(sapply(lss_parts, `[`, 4))
  lss$time <- gsub(".tif", "", sapply(lss_parts, `[`, 8))
  lss %>% arrange(date)
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

# Internal Function to Extract Landsat Imagery from STAC.
# Returns list(meta = metadata of the downloaded scenes, log = outcome per scene).
extract_landsat_stac <- function(aoi, epsg, excl_groups, site_name,
                                 sats = c("LT04", "LT05", "LE07", "LC08", "LC09"),
                                 start_date, end_date, months = 1:12,
                                 minclouds, area_landsat_dir, workers,
                                 max_retries = 5, retry_rounds = 2, retry_wait = 60,
                                 http_retries = 10, http_retry_delay = 5) {

  # Mapping satellites
  sat_map <- c("LT04"="landsat-4", "LT05"="landsat-5", "LE07"="landsat-7", "LC08"="landsat-8", "LC09"="landsat-9")
  sats2 <- sat_map[sats]
  short <- setNames(names(sat_map), sat_map)

  # Search (all pages, retried if the API is busy)
  it_obj <- stac_search_all("landsat-c2-l2", aoi, start_date, end_date)

  # scene group of a catalogue item, in the same form as group_from_file()
  feat_group <- function(x) {
    paste(short[[x$properties$platform]],
          sprintf("%03d", as.integer(x$properties$`landsat:wrs_path`)),
          gsub("-", "", substr(x$properties$datetime, 1, 10)), sep = "_")
  }

  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) { x$properties$`eo:cloud_cover` < minclouds }) %>%
    items_filter(filter_fn = function(x) { x$properties$platform %in% sats2 }) %>%
    assets_select(asset_names = c("coastal", "blue", "green", "red", "nir08", "swir16", "swir22", "lwir", "lwir11", "qa_pixel", "qa_radsat")) %>%
    items_filter(filter_fn = function(x) { month(ymd_hms(x$properties$datetime)) %in% months }) %>%
    items_filter(filter_fn = function(x) { !feat_group(x) %in% excl_groups })

  if (length(it_obj$features) == 0) {
    print("No new Landsat scenes to download!")
    return(list(meta = NULL, log = NULL))
  }

  suppressMessages({
    suppressWarnings({
      itst <- items_as_sf(it_obj) %>% st_make_valid()
      itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
    })
  })

  if (nrow(itst) == 0) {
    print("No new Landsat scenes to download!")
    return(list(meta = NULL, log = NULL))
  }

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
          if (length(unique(x$area)) == 1) {
            x <- x %>% slice(1)
          } else if (nrow(x) > 2) {
            # More than two scenes: keep all, they are merged later in
            # the same-date mosaicking step
            x <- x
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

  # all remaining items can be dropped by the duplicate handling (e.g. a re-run in which the only
  # unknown items are partial duplicates of downloaded scenes)
  if (length(it_obj$features) == 0) {
    print("No new Landsat scenes to download!")
    return(list(meta = NULL, log = NULL))
  }
  message(length(it_obj$features), " Landsat scene(s) to download")

  # --- DOWNLOAD IN PARALLEL (FUTURE), WITH RETRIES ---
  res <- process_features_in_parallel(it_obj, area_landsat_dir, aoi, workers, epsg = epsg,
                                      max_retries = max_retries, retry_rounds = retry_rounds,
                                      retry_wait = retry_wait, http_retries = http_retries,
                                      http_retry_delay = http_retry_delay)

  # Renaming platforms for compatibility
  itst$platform <- unname(short[itst$platform])

  ok_files <- res$file[res$status %in% c("SUCCESS", "EXISTS")]
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
    filter(file %in% ok_files) %>%
    mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>%
    rownames_to_column("id")

  return(list(meta = idsALL, log = res))
}

#' Internal Function: Parallel Feature Processing with Grid Alignment
process_features_in_parallel <- function(it_obj, area_landsat_dir, aoi, workers, epsg,
                                         max_retries = 5, retry_rounds = 2, retry_wait = 60,
                                         http_retries = 10, http_retry_delay = 5) {

  # 1. Grid template: the extent comes from the AOI itself, and only the
  # pixel alignment (origin) is taken from a scene in the target UTM zone.
  # (Earlier the extent was the AOI cropped to one scene, so if that scene
  # covered only part of the AOI, e.g. a neighbouring WRS row, every image
  # got the truncated extent.)
  Sys.setenv(GDAL_HTTP_MAX_RETRY = http_retries, GDAL_HTTP_RETRY_DELAY = http_retry_delay)
  reso_val <- 30
  grid_origin <- c(15, 15) # Landsat C2 UTM grid: pixel corners at 15 m offsets
  idx <- which(vapply(it_obj$features, item_epsg, numeric(1)) == epsg)
  if (length(idx) > 0) {
    tmp_url <- make_vsicurl_url(rstac::assets_url(it_obj$features[[idx[1]]])[[1]])
    r_scene <- try(terra::rast(tmp_url), silent = TRUE)
    if (!inherits(r_scene, "try-error")) {
      grid_origin <- terra::origin(r_scene)
      reso_val <- terra::res(r_scene)[1]
    }
  }
  exte_vals <- aoi_grid_extent(aoi, epsg, reso_val, grid_origin)

  # 2. Parallel downloads with retry rounds (workers: the plan set in extract_landsat)
  download_rounds(it_obj$features, process_feature,
                  list(area_landsat_dir = area_landsat_dir, aoi = aoi, exte = exte_vals, reso = reso_val,
                       max_retries = max_retries, http_retries = http_retries,
                       http_retry_delay = http_retry_delay),
                  retry_rounds = retry_rounds, retry_wait = retry_wait)
}

item_epsg <- function(x) {
  p <- x$properties
  if (!is.null(p$`proj:epsg`)) return(as.numeric(p$`proj:epsg`))
  if (!is.null(p$`proj:code`)) return(as.numeric(sub("^EPSG:", "", p$`proj:code`)))
  NA_real_
}

# Internal Function: AOI bounding box in the target CRS, snapped outwards to
# the satellite pixel grid. Returns c(xmin, xmax, ymin, ymax).
aoi_grid_extent <- function(aoi, epsg, reso, origin) {
  bb <- as.numeric(sf::st_bbox(sf::st_transform(aoi, crs = epsg)))
  snap_down <- function(v, o) floor((v - o) / reso) * reso + o
  snap_up <- function(v, o) ceiling((v - o) / reso) * reso + o
  c(snap_down(bb[1], origin[1]), snap_up(bb[3], origin[1]),
    snap_down(bb[2], origin[2]), snap_up(bb[4], origin[2]))
}


#' Internal Function: Single Feature Extraction with Retry Logic
#'
#' All bands must be read: a band that cannot be opened fails the attempt (earlier
#' versions replaced it with an empty band, which made the image unusable for
#' classification). The image is written under a temporary name and renamed only
#' after it has been read back, so an interrupted run never leaves a partial image
#' that would count as downloaded. Returns a one-row tibble
#' (file, status = SUCCESS / EXISTS / ERROR / BROKEN, attempts, message).
#' @noRd
process_feature <- function(ft, area_landsat_dir, aoi, exte, reso,
                            max_retries = 5, http_retries = 10, http_retry_delay = 5) {

  # self-contained for parallel workers
  make_vsicurl_url <- function(base_url) {
    paste0("/vsicurl", "?pc_url_signing=yes", "&pc_collection=landsat-c2-l2", "&url=", base_url)
  }
  # GDAL retries throttled / failed HTTP requests itself before an attempt fails
  Sys.setenv(GDAL_HTTP_MAX_RETRY = http_retries, GDAL_HTTP_RETRY_DELAY = http_retry_delay)

  nm <- paste0(paste(str_split(ft$id, "_")[[1]][1:4], collapse = "_"), "_",
               gsub("-", "", substr(ft$properties$created, 1, 10)), "_",
               paste(str_split(ft$id, "_")[[1]][5:6], collapse = "_"), "_",
               gsub(":", "", substr(ft$properties$datetime, 12, 19)), "GMT.tif")
  out <- function(status, attempts, msg = "") {
    dplyr::tibble(file = nm, status = status, attempts = as.integer(attempts), message = substr(msg, 1, 500))
  }

  output_file <- file.path(area_landsat_dir, nm)
  if (file.exists(output_file)) return(out("EXISTS", 0))
  tmp_file <- paste0(output_file, ".part")

  full_url <- make_vsicurl_url(rstac::assets_url(ft))
  err <- ""
  prev_bad <- NULL; same_bad <- 0

  for (iii in seq_len(max_retries)) {
    # Open every band separately. A band that fails to open while all the others open is
    # probably a damaged file on the server (rare on the Planetary Computer); when the same
    # band(s) fail in 3 attempts in a row, the scene is given up as BROKEN instead of being
    # retried for minutes. Throttling makes all or varying bands fail, which is retried.
    r_list <- lapply(full_url, function(u) try(terra::rast(u), silent = TRUE))
    bad <- which(vapply(r_list, inherits, logical(1), "try-error"))
    if (length(bad)) {
      err <- paste0("cannot open ", paste(sub(".*/", "", full_url[bad]), collapse = ", "), ": ",
                    trimws(as.character(r_list[[bad[1]]])))
      same_bad <- if (identical(bad, prev_bad)) same_bad + 1 else 1
      prev_bad <- bad
      if (length(bad) < length(full_url) && same_bad >= 3) return(out("BROKEN", iii, err))
      if (iii < max_retries) Sys.sleep(min(300, 5 * 2^(iii - 1)) * stats::runif(1, 0.75, 1.25))
      next
    }
    prev_bad <- NULL; same_bad <- 0
    e <- try({
      r_remote <- terra::rast(r_list)       # all Landsat bands share the 30 m grid

      # Alignment Step: Use Numeric Template to force exact Grid
      template <- terra::rast(extent = terra::ext(exte), resolution = reso, crs = sf::st_crs(aoi)$wkt)

      r_cropped <- terra::crop(r_remote, sf::st_transform(aoi, sf::st_crs(r_remote)))
      r_final <- terra::project(r_cropped, template, method = "near")

      # Metadata formatting
      names(r_final) <- sapply(str_split(names(r_final), "_"), function(x) paste(x[8:9], collapse = "_"))
      r_final <- r_final[[c(names(r_final)[grepl("SR_", names(r_final))] %>% sort,
                            names(r_final)[grepl("ST_", names(r_final))] %>% sort,
                            names(r_final)[grepl("QA_", names(r_final))] %>% sort)]]
      terra::writeRaster(r_final, tmp_file, filetype = "GTiff", datatype = "INT2U",
                         overwrite = TRUE, NAflag = 0)

      # read the written image back before it counts as downloaded
      chk <- terra::rast(tmp_file)
      if (terra::nlyr(chk) != terra::nlyr(r_final)) stop("written image has missing bands")
      invisible(terra::global(chk[[1]], "notNA"))
      if (!file.rename(tmp_file, output_file)) stop("could not rename ", tmp_file)
    }, silent = TRUE)

    if (!inherits(e, "try-error")) return(out("SUCCESS", iii))
    err <- trimws(as.character(e))
    unlink(tmp_file)
    if (iii < max_retries) Sys.sleep(min(300, 5 * 2^(iii - 1)) * stats::runif(1, 0.75, 1.25))
  }
  out("ERROR", max_retries, err)
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
  
  fpath <- file.path(image_dir, image)
  
  res_tibble <- tryCatch({
    rs <- terra::rast(fpath)
    
    # Identify QA band
    qa_band_name <- names(rs)[grepl("_pixel", names(rs), ignore.case = TRUE)]
    if(length(qa_band_name) == 0) stop("QA band not found")
    
    # 1. Fill Proportion
    # Images are written with NAflag = 0, so fill pixels are read back as NA
    # (not 0). Count the NAs in band 4 (a surface reflectance band).
    n_pixels <- terra::ncell(rs)
    fill_count <- sum(is.na(terra::values(rs[[4]], mat = FALSE)))
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
