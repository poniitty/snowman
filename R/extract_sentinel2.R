#' Extracts Sentinel-2 Imagery
#'
#' This function extracts Sentinel-2 imagery that fulfills the specified arguments.
#' It processes the area of interest (AOI), transforms coordinates, and downloads
#' the relevant Sentinel-2 L2A scenes from Microsoft Planetary Computer.
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
#' scenes listed in `imagery/s2s_final.csv` (also if their images have since been
#' deleted, e.g. after classification) and scenes rejected earlier as too cloudy
#' (`imagery/s2s_rejected.csv`). Only the newly downloaded images are checked and
#' added to `s2s_final.csv`.
#'
#' @param aoi An `sf` object or a list containing the area of interest.
#' @param site_name A character string representing the name of the site.
#' @param aoi_size A numeric value representing the size of the AOI in kilometers. Ignored if aoi is a POLYGON.
#' @param start_date A character string representing the start date for the imagery.
#' @param end_date A character string representing the end date for the imagery.
#' @param months A numeric vector representing the months to include (default is all months, so 1:12).
#' @param sats A character vector of satellites (default all, c("S2A", "S2B", "S2C"); Sentinel-2C replaced
#' Sentinel-2A in routine operations in 2025).
#' @param minclouds A numeric value representing the minimum cloud cover percentage (default is 50).
#' @param base_sentinel2_dir A character string representing the base directory for Sentinel-2 imagery.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' Setting workers > 1 (default) enables parallel computing.
#' @param data_source A character string representing the data source (currently only option is "rstac").
#' @param force Logical, if TRUE deletes all previously downloaded imagery and the scene tables.
#' @param max_retries Attempts per scene within one round (default 5), with waits of about 5, 10, 20, 40 s.
#' @param retry_rounds Further rounds for scenes that failed in all their attempts (default 2).
#' @param retry_wait Pause in seconds before retry round k is `k * retry_wait` (default 60).
#' @param http_retries,http_retry_delay GDAL's own retries of failed HTTP requests and the initial
#' delay in seconds (GDAL_HTTP_MAX_RETRY, GDAL_HTTP_RETRY_DELAY).
#' @param groups Optional character vector of scene groups `"<sat>_<tile>_<yyyymmdd>"` (e.g.
#' `"S2A_T33WXR_20170903"`): only these scenes are downloaded (still subject to `minclouds`, so use
#' `minclouds = 101` to get them regardless of cloud cover). The search period is then the range of
#' their dates. Used e.g. to download the scenes of labelled training data.
#' @param max_per_month Optional maximum number of scenes per calendar month (over all years of the
#' period). Scenes covering at least 95\% of the AOI are preferred, then those with the lowest
#' catalogue cloud cover. Scenes already known for the site count towards the maximum, so a re-run
#' does not add more. Used to limit the data volume, e.g. for the monthly median images of
#' [calc_predictors_s2()]. NULL (default) = no limit.
#' @param min_clear Images whose clear proportion (no fill, no cloud) is below this are deleted and
#' recorded as rejected (default 0.1). With 0, only images without any data are rejected
#' (e.g. to keep cloudy scenes of training data).
#' @param bands Which bands to download:
#' \itemize{
#'   \item `"model"` (default): only the bands the classifier needs (see [s2_model_requirements()]): those
#'   of its predictors and those snowman always needs (B02, B03, B04, B08, B11, SCL). With the current
#'   model: B02, B03, B04, B06, B08, B11, B12 and SCL, about a third less data to download and 40\% less
#'   disk space than all bands. If no model is found, all bands are downloaded.
#'   \item `"all"`: all 12 reflectance bands and SCL (e.g. to gather training data for a new model).
#'   \item a character vector of bands, e.g. `c("B02", "B08")`: these, plus the bands snowman always
#'   needs and the bands the model uses (added whatever the choice, with a message).
#' }
#' Images downloaded earlier keep their bands. If they lack bands the model needs, a warning tells
#' how many (they are then not classified by [classify_sentinel2()]; download them again with
#' `force = TRUE` or after deleting them).
#' @param model_file The classifier whose bands are checked (`.rds`). NULL (default): the model
#' included in snowman ([snowman_model_file()]).
#' @return A tibble containing the metadata of all usable Sentinel-2 scenes of the site (the content
#' of `s2s_final.csv`, including the processing baseline of each scene in `PROCESSING_BASELINE`, used
#' to harmonise the reflectance offset of baseline 04.00 and later), or NULL if there are none. Scenes that could not be downloaded in this run
#' are in `attr(result, "failed")`.
#' @export
#' @import rstac dplyr sf terra lubridate stringr parallel tibble readr ggplot2 future future.apply lwgeom
extract_sentinel2 <- function(aoi,
                            site_name,
                            aoi_size = 2,
                            start_date = "2020-01-01",
                            end_date = "2024-12-01",
                            months = 1:12,
                            sats = c("S2A", "S2B", "S2C"),
                            minclouds = 50,
                            base_sentinel2_dir,
                            workers = 1,
                            data_source = "rstac",
                            force = FALSE,
                            max_retries = 5,
                            retry_rounds = 2,
                            retry_wait = 60,
                            http_retries = 10,
                            http_retry_delay = 5,
                            groups = NULL,
                            max_per_month = NULL,
                            min_clear = 0.1,
                            bands = "model",
                            model_file = NULL) {

  # Bands to download: the user's choice, completed with the bands snowman and the model need
  bands <- s2_resolve_bands(bands, model_file)
  message("Sentinel-2 bands downloaded: ", paste(bands, collapse = ", "))

  # Only given scene groups: search over the period of their dates
  if (!is.null(groups)) {
    groups <- unique(as.character(groups))
    g_dates <- lubridate::ymd(vapply(strsplit(groups, "_"), `[`, character(1), 3))
    if (anyNA(g_dates)) stop("groups must be of the form <sat>_<tile>_<yyyymmdd>, e.g. S2A_T33WXR_20170903")
    start_date <- as.character(min(g_dates) - 1)
    end_date <- as.character(max(g_dates) + 1)
  }

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
  area_sentinel2_dir <- file.path(base_sentinel2_dir, site_name, "imagery")
  check_single_sensor(file.path(base_sentinel2_dir, site_name), "sentinel2")
  if (!dir.exists(area_sentinel2_dir)) dir.create(area_sentinel2_dir, recursive = TRUE)

  my_temp_dir <- file.path(base_sentinel2_dir, "temp_terra")
  if(!dir.exists(my_temp_dir)) dir.create(my_temp_dir, recursive = TRUE)
  terra::terraOptions(tempdir = my_temp_dir)

  if(force){
    unlink(list.files(area_sentinel2_dir, pattern = "GMT.tif$", full.names = TRUE))
    unlink(file.path(area_sentinel2_dir, c("s2s.csv", "s2s_final.csv", "s2s_rejected.csv")))
  }
  unlink(list.files(area_sentinel2_dir, pattern = "\\.part$", full.names = TRUE)) # leftovers of interrupted runs

  # Scenes that need no download (see the function description), and the scene
  # table of earlier runs, to which new scenes are added
  known <- known_groups(area_sentinel2_dir, "s2s_final.csv", c("s2s_rejected.csv", "s2s_merged.csv"),
                        sensor = "sentinel2")
  final_file <- file.path(area_sentinel2_dir, "s2s_final.csv")
  prev_final <- if (file.exists(final_file)) {
    readr::read_csv(final_file, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else NULL

  # Images of earlier runs without bands needed now (e.g. downloaded for an earlier model)
  if (!is.null(prev_final)) {
    on_disk0 <- intersect(prev_final$file, list.files(area_sentinel2_dir, pattern = "GMT.tif$"))
    if (length(on_disk0)) {
      lack <- vapply(s2_file_bands(on_disk0, area_sentinel2_dir), function(b) length(setdiff(bands, b)) > 0, logical(1))
      if (any(lack)) {
        warning(sum(lack), " of ", length(on_disk0), " images downloaded earlier lack some of the bands ",
                paste(bands, collapse = ", "), ". They are kept as they are; download them again ",
                "(force = TRUE, or delete them and their rows of s2s_final.csv) if the model needs these bands.",
                call. = FALSE)
      }
    }
  }

  # --- 4. SEARCH & DOWNLOAD ---
  if (data_source == "rstac") {
    dl <- extract_sentinel2_stac(workers = workers, start_date = start_date, end_date = end_date, months = months,
                               aoi = aoi, site_name = site_name, epsg = epsg, excl_groups = known,
                               sats = sats, area_sentinel2_dir = area_sentinel2_dir, minclouds = minclouds,
                               groups = groups, max_per_month = max_per_month,
                               max_retries = max_retries, retry_rounds = retry_rounds, retry_wait = retry_wait,
                               http_retries = http_retries, http_retry_delay = http_retry_delay,
                               bands = bands)
  }
  gc()
  write_download_log(dl$log, area_sentinel2_dir)
  report_downloads(dl$log, length(known), area_sentinel2_dir)
  failed <- if (is.null(dl$log)) NULL else dl$log %>% filter(status %in% c("ERROR", "BROKEN"))
  # tiles that do not reach the AOI are recorded as rejected, so they are not searched for again
  no_ov <- if (is.null(dl$log)) NULL else dl$log %>% filter(status == "NO_OVERLAP")
  if (!is.null(no_ov) && nrow(no_ov) > 0) {
    rej_file <- file.path(area_sentinel2_dir, "s2s_rejected.csv")
    write_csv(parse_s2_files(tibble(file = no_ov$file), site_name) %>%
                transmute(file, date, satid, fill_proportion = 1, cloud_proportion = NA, clear_proportion = 0) %>%
                mutate(across(everything(), as.character)),
              rej_file, append = file.exists(rej_file))
  }

  # s2s.csv: catalogue metadata of every downloaded scene (also of those later
  # rejected), accumulated over runs
  meta_cols <- c("id", "file", "DATE_ACQUIRED", "CLOUD_COVER", "SUN_ELEVATION", "PROCESSING_BASELINE")
  meta <- if (!is.null(dl$meta)) dl$meta %>% mutate(across(everything(), as.character)) else NULL
  if (file.exists(file.path(area_sentinel2_dir, "s2s.csv"))) {
    old_meta <- readr::read_csv(file.path(area_sentinel2_dir, "s2s.csv"), show_col_types = FALSE,
                                col_types = readr::cols(.default = "c"))
    meta <- bind_rows(old_meta[, intersect(meta_cols, names(old_meta))], meta)
  }
  if (!is.null(meta) && nrow(meta) > 0) {
    meta <- meta %>% distinct(file, .keep_all = TRUE)
    write_csv(meta, file.path(area_sentinel2_dir, "s2s.csv"))
  }

  # --- 5. NEW IMAGES: on disk but not yet in the final scene table ---
  tifs <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")
  new_files <- setdiff(tifs, prev_final$file)
  # Images of the same acquisition (several tiles) are mosaicked into one (step 7). Earlier images of
  # an acquisition with a new tile image, and earlier images that were never mosaicked (downloaded by
  # older versions), are processed again with them.
  if (!is.null(prev_final)) {
    on_disk <- intersect(prev_final$file, tifs)
    cand <- c(new_files, on_disk)
    acq_k <- s2_acquisition_key(cand, meta)
    redo <- intersect(on_disk, cand[acq_k %in% acq_k[duplicated(acq_k)]])
    if (length(redo) > 0) {
      new_files <- c(new_files, redo)
      prev_final <- prev_final %>% filter(!file %in% redo)
    }
  }
  if (length(new_files) == 0) {
    if (is.null(prev_final)) {
      print("Zero Sentinel-2 scenes found or downloaded!")
      return(NULL)
    }
    out <- readr::type_convert(prev_final, col_types = readr::cols())
    attr(out, "failed") <- failed
    return(out)
  }

  s2s <- tibble(file = new_files)
  if (!is.null(meta)) s2s <- left_join(s2s, meta[, intersect(meta_cols, names(meta))], by = "file")
  s2s <- parse_s2_files(s2s, site_name)

  # --- 6. CHECK RASTERS (FUTURE APPLY) ---
  # Unreadable images are deleted; they are not recorded as rejected, so they are
  # downloaded again on the next run.
  img_remove <- future.apply::future_lapply(
    s2s$file,
    check_raster,
    image_dir = area_sentinel2_dir,
    future.seed = TRUE,
    future.packages = c("terra"),
    future.scheduling = 5
  ) %>% unlist()
  gc()
  if (length(img_remove) > 0) {
    message(length(img_remove), " unreadable image(s) deleted; they are downloaded again on the next run")
    unlink(file.path(area_sentinel2_dir, img_remove))
    s2s <- s2s %>% filter(!file %in% img_remove)
  }

  # --- 7. MOSAICKING ---
  # The tiles of one acquisition (satellite, date and time; the AOI may straddle MGRS tiles, which
  # overlap by ~10 km) and copies of one tile are mosaicked into one image, so that an acquisition is
  # one observation. Each pixel comes from one input image only, all bands and the SCL together (the
  # image with the most valid pixels first, the others fill its gaps); values are never averaged.
  # Tiles with and without the +1000 offset of processing baseline 04.00 are not mixed. The image
  # keeps the name of the first input; the others are deleted and recorded in s2s_merged.csv, so that
  # they are not downloaded again.
  s2s$.acq <- s2_acquisition_key(s2s$file, meta)
  dup_acq <- unique(s2s$.acq[duplicated(s2s$.acq)])
  if (length(dup_acq) > 0) {
    merged_log <- list()
    for (k in dup_acq) {
      fs <- s2s$file[s2s$.acq == k]
      nvalid <- vapply(fs, function(x) {
        r <- s2_named(rast(file.path(area_sentinel2_dir, x)), x)[["B04"]]
        as.numeric(terra::global(!is.na(r) & r > 0, "sum", na.rm = TRUE)[1, 1])
      }, numeric(1))
      fs <- fs[order(-nvalid)]
      # only the bands all of them have (images of earlier runs may have more bands than new ones)
      cb <- Reduce(intersect, s2_file_bands(fs, area_sentinel2_dir))
      cb <- s2_band_key[s2_band_key %in% cb]
      rs <- lapply(fs, function(x) {
        r <- s2_named(rast(file.path(area_sentinel2_dir, x)), x)[[cb]]
        terra::mask(r, r[["B04"]], maskvalues = c(NA, 0))   # a pixel without data is empty in all bands
      })
      rs <- terra::mosaic(terra::sprc(rs), fun = "first")
      tmp <- file.path(area_sentinel2_dir, paste0(fs[1], ".part"))
      writeRaster(rs, tmp, filetype = "GTiff", overwrite = TRUE, datatype = "INT2U", NAflag = 0)
      file.rename(tmp, file.path(area_sentinel2_dir, fs[1]))
      unlink(file.path(area_sentinel2_dir, fs[-1]))
      merged_log[[k]] <- tibble(file = fs[-1], merged_into = fs[1])
      s2s <- s2s %>% filter(!file %in% fs[-1])
    }
    merged_log <- bind_rows(merged_log)
    mfile <- file.path(area_sentinel2_dir, "s2s_merged.csv")
    write_csv(merged_log, mfile, append = file.exists(mfile))
    message(nrow(merged_log), " tile image(s) mosaicked into ", length(dup_acq), " acquisition(s)")
  }
  s2s$.acq <- NULL

  # --- 8. CALC COVERAGES (FUTURE APPLY) ---
  s2cs <- future.apply::future_lapply(
    s2s$file,
    calc_coverages_s2,
    image_dir = area_sentinel2_dir,
    future.seed = TRUE,
    future.packages = c("terra", "dplyr"),
    future.scheduling = 5
  ) %>% bind_rows()
  gc()

  s2s <- full_join(s2s %>% select(-ends_with("proportion")), s2cs, by = "file")

  # Mostly cloudy or empty images are deleted and recorded, so they are not downloaded again
  rejected <- s2s %>% filter(is.na(clear_proportion) | fill_proportion >= 0.99 | clear_proportion < min_clear)
  if (nrow(rejected) > 0) {
    unlink(file.path(area_sentinel2_dir, rejected$file))
    rej_file <- file.path(area_sentinel2_dir, "s2s_rejected.csv")
    write_csv(rejected %>% select(file, date, satid, fill_proportion, cloud_proportion, clear_proportion) %>%
                mutate(across(everything(), as.character)),
              rej_file, append = file.exists(rej_file))
  }
  s2s <- s2s %>% filter(!file %in% rejected$file)

  # --- 9. FINAL SCENE TABLE: earlier scenes + the new ones ---
  final <- bind_rows(prev_final, s2s %>% mutate(across(everything(), as.character))) %>%
    distinct(file, .keep_all = TRUE)
  final <- readr::type_convert(final, col_types = readr::cols()) %>% arrange(date)
  write_csv(final, final_file)
  message(nrow(s2s), " new scene(s) added; ", nrow(final), " scenes in ", final_file)

  attr(final, "failed") <- failed
  return(final)
}

# Internal Function: scene metadata from the image file names
# (<satid>_MSIL2A_<tile>_<date>_<processing date>_<time>GMT.tif)
parse_s2_files <- function(s2s, site_name) {
  s2s$area  <- site_name
  s2s_parts <- str_split(s2s$file, "_")
  s2s$satid <- sapply(s2s_parts, `[`, 1)                                  # S2A / S2B / S2C
  s2s$level <- sapply(s2s_parts, `[`, 2)                                  # MSIL2A
  s2s$tile  <- sapply(s2s_parts, `[`, 3)                                  # MGRS tile, e.g. T33VWF
  s2s$date  <- ymd(sapply(s2s_parts, `[`, 4))                             # acquisition date
  s2s$time  <- gsub(".tif", "", sapply(s2s_parts, `[`, 6))                # UTC time string
  s2s %>% arrange(date)
}

# Internal Function to Create a Signed Sentinel-2 vsicurl URL
make_vsicurl_url_s2 <- function(base_url) {
  paste0(
    "/vsicurl",
    "?pc_url_signing=yes",
    "&pc_collection=sentinel-2-l2a",   # differs from landsat-c2-l2
    "&url=",
    base_url
  )
}

# Internal: short satellite code of a Planetary Computer platform name
s2_short <- function(platform) {
  unname(c("Sentinel-2A" = "S2A", "Sentinel-2B" = "S2B", "Sentinel-2C" = "S2C", "Sentinel-2D" = "S2D")[platform])
}

# Internal Function to Extract Sentinel-2 Imagery from STAC.
# Returns list(meta = metadata of the downloaded scenes, log = outcome per scene).
extract_sentinel2_stac <- function(aoi, epsg, excl_groups, site_name,
                                   sats = c("S2A", "S2B", "S2C"),
                                   start_date, end_date, months = 1:12,
                                   minclouds, area_sentinel2_dir, workers,
                                   groups = NULL, max_per_month = NULL,
                                   max_retries = 5, retry_rounds = 2, retry_wait = 60,
                                   http_retries = 10, http_retry_delay = 5, bands = s2_band_key) {

  # Search (all pages, retried if the API is busy)
  it_obj <- stac_search_all("sentinel-2-l2a", aoi, start_date, end_date)

  # scene group of a catalogue item <sat>_<tile>_<date> (the form of the `groups` argument), and its
  # key <sat>_<tile>_<date>_<time> in the same form as group_from_file() (the time separates two
  # passes of one satellite on the same day, which happen at high latitudes)
  # (s2:mgrs_tile is stored without the leading "T")
  feat_group <- function(x) {
    paste(s2_short(x$properties$platform), paste0("T", x$properties$`s2:mgrs_tile`),
          gsub("-", "", substr(x$properties$datetime, 1, 10)), sep = "_")
  }
  feat_key <- function(x) {
    paste(feat_group(x), gsub(":", "", substr(x$properties$datetime, 12, 19)), sep = "_")
  }

  # Filter by cloud cover, platform, months, and already-known scenes.
  # Band selection: the bands of `bands` (reflectance bands + SCL, the Scene Classification Layer = QA
  # equivalent). All bands are resampled to 10 m during download (see process_feature_s2).
  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) { cc <- x$properties$`eo:cloud_cover`; is.null(cc) || isTRUE(cc < minclouds) }) %>%
    items_filter(filter_fn = function(x) { isTRUE(s2_short(x$properties$platform) %in% sats) }) %>%
    assets_select(asset_names = bands) %>%
    items_filter(filter_fn = function(x) { month(ymd_hms(x$properties$datetime)) %in% months }) %>%
    items_filter(filter_fn = function(x) { !feat_key(x) %in% excl_groups })
  if (!is.null(groups)) {
    it_obj <- it_obj %>% items_filter(filter_fn = function(x) { feat_group(x) %in% groups })
  }

  if (length(it_obj$features) == 0) {
    print("No new Sentinel-2 scenes to download!")
    return(list(meta = NULL, log = NULL))
  }

  suppressMessages({
    suppressWarnings({
      itst <- items_as_sf(it_obj) %>% st_make_valid()
      # the item id identifies a scene item (every item has one; properties such as s2:granule_id,
      # created or view:sun_elevation are missing from some older items)
      itst$.item_id <- vapply(it_obj$features, function(x) as.character(x$id), character(1))
      for (nm in c("eo:cloud_cover", "view:sun_elevation")) if (!nm %in% names(itst)) itst[[nm]] <- NA_real_
      # the part of each footprint inside the AOI itself (cropping with the AOI's bounding box in
      # degrees kept tiles that only touch the box's corners; their download then failed)
      sf::sf_use_s2(FALSE)
      itst <- st_intersection(itst, st_geometry(aoi %>% st_transform(st_crs(itst))))
      itst <- itst[!st_is_empty(itst), ]
      sf::sf_use_s2(TRUE)
    })
  })
  itst <- itst[as.numeric(st_area(itst)) > 1e4, ]   # at least 0.01 km2 of the AOI

  if (nrow(itst) == 0) {
    print("No new Sentinel-2 scenes to download!")
    return(list(meta = NULL, log = NULL))
  }

  itst$area <- as.numeric(st_area(itst)) / (1000 * 1000)

  # Only the tiles needed to cover the AOI: an acquisition's tiles (MGRS tiles overlap, by much more
  # than 10 km near UTM zone edges) are added largest-first, by their catalogue footprints, as long as
  # each adds at least 0.1% of the AOI; once 99.9% is covered the rest are not downloaded. Ties go to
  # tiles in the AOI's UTM zone (no reprojection). Not for `groups` (those tiles are wanted as such).
  if (is.null(groups)) {
    # the AOI's area measured like the footprints' (in geographic coordinates), so that a covering
    # tile gives exactly 100%
    aoi_km2_all <- as.numeric(st_area(st_transform(aoi, st_crs(itst)))) / (1000 * 1000)
    ep <- setNames(vapply(it_obj$features, item_epsg, numeric(1)),
                   vapply(it_obj$features, function(x) as.character(x$id), character(1)))
    itst$.epsg_ok <- ep[itst$.item_id] %in% epsg
    acq_id <- paste(itst$platform, substr(itst$datetime, 1, 10), substr(itst$datetime, 12, 19))
    # Copies of the same tile and acquisition (the catalogue keeps reprocessed versions, e.g. of 2022
    # scenes reprocessed in 2024): only the newest processing is considered here; the older copies are
    # downloaded only if the check after the download finds the AOI incompletely covered.
    prop0 <- function(nm) if (nm %in% names(itst)) itst[[nm]] else NULL
    base_num <- suppressWarnings(as.numeric(s2_baseline_prop(prop0("s2:processing_baseline"),
                                                             prop0("s2:product_uri"), nrow(itst))))
    pdate <- suppressWarnings(as.numeric(s2_procdate(prop0("created"), prop0("s2:generation_time"),
                                                     prop0("s2:product_uri"), itst$datetime)))
    copy_key <- paste(acq_id, itst$`s2:mgrs_tile`)
    newest <- vapply(seq_len(nrow(itst)), function(i) {
      ix <- which(copy_key == copy_key[i])
      ix[order(-ifelse(is.na(base_num[ix]), -1, base_num[ix]), -ifelse(is.na(pdate[ix]), 0, pdate[ix]))][1] == i
    }, logical(1))
    n_old <- sum(!newest)
    if (n_old > 0) message(n_old, " older processing(s) of a tile left out (the newest version of each is used)")
    geo <- st_geometry(itst)
    suppressMessages(suppressWarnings({
      sf::sf_use_s2(FALSE)
      keep <- unlist(lapply(split(which(newest), acq_id[newest]), function(ix) {
        if (length(ix) == 1) return(ix)
        ix <- ix[order(-round(itst$area[ix] / aoi_km2_all, 3), !itst$.epsg_ok[ix])]
        chosen <- ix[1]
        cov <- itst$area[ix[1]]
        for (j in ix[-1]) {
          if (cov >= 0.999 * aoi_km2_all) break
          new_cov <- as.numeric(st_area(st_union(geo[c(chosen, j)]))) / (1000 * 1000)
          if (new_cov - cov > 0.001 * aoi_km2_all) { chosen <- c(chosen, j); cov <- new_cov }
        }
        chosen
      }))
      sf::sf_use_s2(TRUE)
    }))
    n_skip <- sum(newest) - length(keep)
    if (n_skip > 0) message(n_skip, " tile(s) not needed: the AOI is covered by the other tiles of the same acquisition")
    itst$.epsg_ok <- NULL
    itst$.newest <- newest
    # kept for the check after the download: a catalogue footprint may claim more than the tile has data for
    itst_skipped <- itst[-keep, ]
    it_skipped <- it_obj %>% items_filter(filter_fn = function(x) as.character(x$id) %in% itst_skipped$.item_id)
    itst <- itst[sort(keep), ]
  } else {
    itst_skipped <- itst[0, ]
    it_skipped <- NULL
  }

  # At most max_per_month acquisitions per calendar month (known ones included). The unit is an
  # acquisition (satellite + date + time), with all its tiles (they are mosaicked). Within a month,
  # acquisitions are taken in turn from the different years (so that the medians do not come from a
  # few years only); within a year, those covering most of the AOI and with the lowest catalogue
  # cloud cover first.
  if (!is.null(max_per_month)) {
    aoi_km2 <- as.numeric(st_area(st_transform(aoi, st_crs(itst)))) / (1000 * 1000)
    kp <- strsplit(as.character(excl_groups), "_")
    known_acq <- unique(vapply(kp, function(z) paste(z[1], z[3], if (length(z) >= 4) z[4] else ""), character(1)))
    known_dates <- suppressWarnings(ymd(vapply(strsplit(known_acq, " "), `[`, character(1), 2)))
    known_mo <- table(factor(month(known_dates), levels = 1:12))
    acq <- itst %>%
      st_drop_geometry() %>%
      mutate(.acq = paste(platform, substr(datetime, 1, 10), substr(datetime, 12, 19)),
             .known = paste(s2_short(platform), gsub("-", "", substr(datetime, 1, 10)),
                            gsub(":", "", substr(datetime, 12, 19)))) %>%
      # other tiles of acquisitions already downloaded (e.g. training scenes) are not new acquisitions
      filter(!.known %in% known_acq) %>%
      group_by(.acq) %>%
      summarise(.month = month(ymd_hms(datetime[1])), .year = year(ymd_hms(datetime[1])),
                .cover = min(1, sum(area) / aoi_km2), .cc = mean(`eo:cloud_cover`, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(.cc = ifelse(is.na(.cc), 100, .cc),
             .covclass = cut(.cover, c(-Inf, 0.5, 0.9, Inf), labels = FALSE)) %>%
      group_by(.month, .year) %>%
      arrange(desc(.covclass), .cc, .by_group = TRUE) %>%
      mutate(.rank = row_number(), .best = min(.cc)) %>%
      group_by(.month) %>%
      arrange(.rank, .best, .by_group = TRUE) %>%
      filter(row_number() <= max(0, max_per_month - known_mo[as.character(.month[1])])) %>%
      ungroup()
    itst <- itst %>%
      mutate(.acq = paste(platform, substr(datetime, 1, 10), substr(datetime, 12, 19))) %>%
      filter(.acq %in% acq$.acq) %>%
      select(-.acq)
    if (nrow(itst) == 0) {
      print("No new Sentinel-2 scenes to download (max_per_month reached)!")
      return(list(meta = NULL, log = NULL))
    }
  }

  # Duplicate handling: items of the same tile and acquisition (satellite, date and time), e.g. copies
  # of one scene with different processing baselines
  tt <- itst %>%
    mutate(date = as_date(ymd_hms(datetime)), .time = substr(datetime, 12, 19)) %>%
    group_by(date, .time, platform, `s2:mgrs_tile`, .add = TRUE) %>%
    mutate(n = n()) %>%
    arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()

  suppressMessages({
    suppressWarnings({
      sf::sf_use_s2(FALSE)
      tt <- lapply(tt, function(x) {
        if (nrow(x) > 1) {
          if (length(unique(x$area)) == 1) {
            x <- x %>% slice(1)
          } else if (nrow(x) > 2) {
            # More than two items: keep all, they are merged in the mosaicking step
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

  # the item id is the unique scene identifier (analogous to landsat:scene_id)
  itst <- itst %>%
    filter(.item_id %in% (tt %>% pull(.item_id)))

  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) {
      as.character(x$id) %in% (itst %>% pull(.item_id))
    })

  # all remaining items can be dropped by the duplicate handling (e.g. a re-run in which the only
  # unknown items are partial duplicates of downloaded scenes)
  if (length(it_obj$features) == 0) {
    print("No new Sentinel-2 scenes to download!")
    return(list(meta = NULL, log = NULL))
  }
  message(length(it_obj$features), " Sentinel-2 scene(s) to download")

  # --- DOWNLOAD IN PARALLEL (FUTURE), WITH RETRIES ---
  res <- process_features_in_parallel_s2(it_obj, area_sentinel2_dir, aoi, workers, epsg = epsg,
                                         max_retries = max_retries, retry_rounds = retry_rounds,
                                         retry_wait = retry_wait, http_retries = http_retries,
                                         http_retry_delay = http_retry_delay, bands = bands)

  # Check of the skipped tiles against the data: if the downloaded tiles of an acquisition leave
  # more than 0.5% of the AOI without valid pixels (a footprint can be larger than the tile's
  # data), the skipped tiles of that acquisition are downloaded as well (and mosaicked later).
  # (first the other tiles; older processings of a tile only if a gap remains after those)
  if (!is.null(it_skipped) && length(it_skipped$features) > 0) {
    acq_of <- function(f) vapply(strsplit(basename(f), "_"), function(z) paste(z[1], z[4], sub("GMT.tif$", "", z[6])), character(1))
    sk_ids <- vapply(it_skipped$features, function(x) as.character(x$id), character(1))
    sk_files <- vapply(it_skipped$features, s2_item_file, character(1))
    sk_newest <- itst_skipped$.newest[match(sk_ids, itst_skipped$.item_id)]
    for (stage in c(TRUE, FALSE)) {
      ok_now <- res$file[res$status %in% c("SUCCESS", "EXISTS")]
      on_disk <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")   # incl. earlier runs' images
      cand <- which(sk_newest %in% stage & !sk_files %in% res$file & acq_of(sk_files) %in% acq_of(ok_now))
      gap_acq <- unique(Filter(function(a) s2_valid_fraction(on_disk[acq_of(on_disk) == a], area_sentinel2_dir, aoi) < 0.995,
                               unique(acq_of(sk_files[cand]))))
      need <- cand[acq_of(sk_files[cand]) %in% gap_acq]
      if (!length(need)) next
      message(length(need), " skipped ", if (stage) "tile(s)" else "older processing(s)",
              " downloaded after all: the other tiles left part of the AOI without data")
      it_need <- it_skipped
      it_need$features <- it_skipped$features[need]
      res2 <- process_features_in_parallel_s2(it_need, area_sentinel2_dir, aoi, workers, epsg = epsg,
                                              max_retries = max_retries, retry_rounds = retry_rounds,
                                              retry_wait = retry_wait, http_retries = http_retries,
                                              http_retry_delay = http_retry_delay, bands = bands)
      res <- bind_rows(res, res2)
      it_obj$features <- c(it_obj$features, it_need$features)
      itst <- bind_rows(itst, itst_skipped %>% filter(.item_id %in% sk_ids[need]))
    }
  }

  itst$platform <- s2_short(itst$platform)

  # File name as written by process_feature_s2:
  # S2A_MSIL2A_T33VWF_20210101_20210101_100031GMT.tif
  ok_files <- res$file[res$status %in% c("SUCCESS", "EXISTS")]
  # catalogue properties that not all items have (e.g. older scenes have no "created")
  prop <- function(nm) if (nm %in% names(itst)) itst[[nm]] else NULL
  itst$.baseline <- s2_baseline_prop(prop("s2:processing_baseline"), prop("s2:product_uri"), nrow(itst))
  # file names exactly as process_feature_s2() writes them
  item_files <- setNames(vapply(it_obj$features, s2_item_file, character(1)),
                         vapply(it_obj$features, function(x) as.character(x$id), character(1)))
  idsALL <- itst %>%
    mutate(date  = gsub("-", "", as_date(ymd_hms(datetime)))) %>%
    mutate(file  = unname(item_files[.item_id])) %>%
    mutate(PROCESSING_BASELINE = .baseline) %>%
    select(file, date, `eo:cloud_cover`, `view:sun_elevation`, PROCESSING_BASELINE) %>%
    rename(DATE_ACQUIRED = date,
           CLOUD_COVER   = `eo:cloud_cover`,
           SUN_ELEVATION = `view:sun_elevation`) %>%
    st_drop_geometry() %>%
    filter(file %in% ok_files) %>%
    mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>%
    rownames_to_column("id")

  return(list(meta = idsALL, log = res))
}

#' Internal Function: Parallel Feature Processing (Sentinel-2) with Grid Alignment
process_features_in_parallel_s2 <- function(it_obj, area_sentinel2_dir, aoi, workers, epsg,
                                            max_retries = 5, retry_rounds = 2, retry_wait = 60,
                                            http_retries = 10, http_retry_delay = 5, bands = s2_band_key) {

  # Build grid template from a feature whose native EPSG matches the target.
  # Use B04 (red, 10 m) as the template band. The extent comes from the AOI
  # itself (not from one tile, which may cover only part of the AOI); only the
  # pixel alignment is taken from a tile.
  Sys.setenv(GDAL_HTTP_MAX_RETRY = http_retries, GDAL_HTTP_RETRY_DELAY = http_retry_delay)
  reso_val    <- 10
  grid_origin <- c(0, 0) # Sentinel-2 UTM grid: pixel corners at 10 m multiples
  idx <- which(vapply(it_obj$features, item_epsg, numeric(1)) == epsg)
  if (length(idx) > 0) {
    # (rstac::assets_url() returns the URLs without asset names, so B04 is taken from the item itself)
    tmp_url <- make_vsicurl_url_s2(it_obj$features[[idx[1]]]$assets$B04$href)
    r_scene <- try(terra::rast(tmp_url), silent = TRUE)
    if (!inherits(r_scene, "try-error")) {
      grid_origin <- terra::origin(r_scene)
      reso_val    <- terra::res(r_scene)[1]
    }
  }
  exte_vals <- aoi_grid_extent(aoi, epsg, reso_val, grid_origin)

  download_rounds(it_obj$features, process_feature_s2,
                  list(area_sentinel2_dir = area_sentinel2_dir, aoi = aoi, exte = exte_vals, reso = reso_val,
                       max_retries = max_retries, http_retries = http_retries,
                       http_retry_delay = http_retry_delay, band_key = bands),
                  retry_rounds = retry_rounds, retry_wait = retry_wait)
}

#' Internal Function: Single Sentinel-2 Feature Extraction with Retry Logic
#'
#' All bands of `band_key` must be read: a band that cannot be opened fails the attempt (earlier
#' versions replaced it with an empty band). The image is written under a
#' temporary name and renamed only after it has been read back. Returns a one-row
#' tibble (file, status = SUCCESS / EXISTS / ERROR / BROKEN, attempts, message).
#' @noRd
process_feature_s2 <- function(ft, area_sentinel2_dir, aoi, exte, reso,
                               max_retries = 5, http_retries = 10, http_retry_delay = 5,
                               band_key = c("B01", "B02", "B03", "B04", "B05", "B06",
                                            "B07", "B08", "B8A", "B09", "B11", "B12", "SCL")) {

  # Redeclare locally so the function is self-contained in workers
  make_vsicurl_url_s2 <- function(base_url) {
    paste0("/vsicurl", "?pc_url_signing=yes", "&pc_collection=sentinel-2-l2a", "&url=", base_url)
  }
  Sys.setenv(GDAL_HTTP_MAX_RETRY = http_retries, GDAL_HTTP_RETRY_DELAY = http_retry_delay)
  # terra judges free memory from the whole node, not the job's limit: with many workers each
  # holding whole 10 m scenes in memory, R sessions were killed. Larger rasters go via temporary files.
  terra::terraOptions(memmax = 1)
  # GDAL's block cache defaults to 5% of the node's memory (tens of GB on an HPC node) in every
  # worker, and fills up while the scenes are read
  if (!nzchar(Sys.getenv("GDAL_CACHEMAX"))) Sys.setenv(GDAL_CACHEMAX = "256")
  terra::setGDALconfig("GDAL_CACHEMAX", Sys.getenv("GDAL_CACHEMAX"))

  # Build output file name ---------------------------------------------------
  # Format: S2A_MSIL2A_T33VWF_20210101_20210101_100031GMT.tif
  nm <- s2_item_file(ft)
  out <- function(status, attempts, msg = "") {
    dplyr::tibble(file = nm, status = status, attempts = as.integer(attempts), message = substr(msg, 1, 500))
  }
  output_file <- file.path(area_sentinel2_dir, nm)
  if (file.exists(output_file)) return(out("EXISTS", 0))
  tmp_file <- paste0(output_file, ".part")
  # --------------------------------------------------------------------------

  # band files in band_key order (taken by asset name; rstac::assets_url() returns unnamed URLs)
  hrefs <- vapply(band_key, function(b) { h <- ft$assets[[b]]$href; if (is.null(h)) NA_character_ else h }, character(1))
  if (anyNA(hrefs)) return(out("BROKEN", 0, paste("bands missing in the catalogue:", paste(band_key[is.na(hrefs)], collapse = ", "))))
  all_urls <- ifelse(grepl("^https?://", hrefs), make_vsicurl_url_s2(hrefs), hrefs)   # (local paths: tests)
  vrt <- paste0(tmp_file, ".vrt")
  # target grid: the AOI's UTM grid at 10 m
  t_srs <- if (!is.na(sf::st_crs(aoi)$epsg)) paste0("EPSG:", sf::st_crs(aoi)$epsg) else sf::st_crs(aoi)$wkt
  cut0 <- terra::as.polygons(terra::ext(exte), crs = sf::st_crs(aoi)$wkt)
  err <- ""
  prev_bad <- NULL; same_bad <- 0

  for (iii in seq_len(max_retries)) {
    # The bands (10, 20 and 60 m) are combined into one virtual raster (VRT, at 10 m) and GDAL's
    # warp reads the AOI window from it straight onto the 10 m UTM grid of the AOI (nearest
    # neighbour) and writes the image. This replaced reading the window of each band with terra and
    # reprojecting it with terra::project(), which used most of the time per scene and slowed down
    # sharply with parallel workers (test on Roihu: 8 s -> 59 s per scene with 1 -> 8 workers, while
    # GDAL's warp took 12 -> 21 s for open, read, warp and write together).
    e <- try({
      unlink(c(vrt, tmp_file))
      sf::gdal_utils("buildvrt", source = all_urls, destination = vrt, quiet = TRUE,
                     options = c("-separate", "-resolution", "highest", "-srcnodata", "0", "-vrtnodata", "0"))
      x <- readLines(vrt)
      # a band file that could not be opened is left out of the VRT (with a warning)
      present <- vapply(hrefs, function(h) any(grepl(basename(h), x, fixed = TRUE)), logical(1))
      bad <- which(!present)
      if (length(bad)) {
        same_bad <- if (identical(bad, prev_bad)) same_bad + 1 else 1
        prev_bad <- bad
        stop("cannot open ", paste(band_key[bad], collapse = ", "))
      }
      prev_bad <- NULL; same_bad <- 0
      # band names (GDAL band descriptions, which terra reads as layer names)
      bi <- grep("<VRTRasterBand", x)
      for (k in rev(seq_along(bi))) x <- append(x, paste0("    <Description>", band_key[k], "</Description>"), after = bi[k])
      writeLines(x, vrt)
      # a tile whose data do not reach the AOI: nothing to download (not an error, not retried)
      rv <- terra::rast(vrt)
      if (is.null(terra::intersect(terra::ext(rv), terra::ext(terra::project(cut0, terra::crs(rv)))))) stop("NO_OVERLAP")
      sf::gdal_utils("warp", vrt, tmp_file, quiet = TRUE,
                     options = c("-t_srs", t_srs, "-te", exte[1], exte[3], exte[2], exte[4],
                                 "-tr", reso, reso, "-r", "near", "-ot", "UInt16", "-dstnodata", "0",
                                 "-of", "GTiff", "-co", "COMPRESS=LZW", "-co", "TILED=YES", "-overwrite"))
      chk <- terra::rast(tmp_file)
      if (terra::nlyr(chk) != length(band_key)) stop("written image has missing bands")
      if (!identical(names(chk), band_key)) stop("unexpected band names in the written image")
      invisible(terra::global(chk[["B04"]], "notNA"))
      if (!file.rename(tmp_file, output_file)) stop("could not rename ", tmp_file)
    }, silent = TRUE)
    unlink(vrt)
    if (!inherits(e, "try-error")) return(out("SUCCESS", iii))
    err <- trimws(as.character(e))
    unlink(tmp_file)
    if (grepl("NO_OVERLAP", err)) return(out("NO_OVERLAP", iii, "the tile does not reach the AOI"))
    # the same band(s) failing to open in 3 attempts in a row while the others open: probably a
    # damaged file on the server; given up as BROKEN instead of being retried for minutes
    if (length(prev_bad) && length(prev_bad) < length(band_key) && same_bad >= 3) return(out("BROKEN", iii, err))
    if (iii < max_retries) Sys.sleep(min(300, 5 * 2^(iii - 1)) * stats::runif(1, 0.75, 1.25))
  }
  out("ERROR", max_retries, err)
}

# Internal Function to calculate cloud and fill coverage within Sentinel-2 rasters.
# Uses the Scene Classification Layer (SCL) instead of Landsat's QA_PIXEL bit logic.
# SCL class values relevant here:
#   0  = No data (fill)
#   3  = Cloud shadow
#   8  = Cloud medium probability
#   9  = Cloud high probability
#   10 = Thin cirrus
calc_coverages_s2 <- function(image, image_dir) {
  terra::terraOptions(memmax = 1)

  fpath <- file.path(image_dir, image)

  res_tibble <- tryCatch({
    rs <- s2_named(terra::rast(fpath), image)

    # Identify SCL band (analogous to QA_PIXEL in Landsat)
    scl_band_name <- names(rs)[grepl("^SCL$", names(rs), ignore.case = TRUE)]
    if (length(scl_band_name) == 0) stop("SCL band not found")

    n_pixels <- terra::ncell(rs)

    # 1. Fill Proportion: pixels where B04 == 0 indicate no-data fill areas.
    #    (Mirrors the Landsat approach of checking band 4 for zero values.)
    #    Images are written with NAflag = 0, so fill pixels are read back as NA.
    fill_count    <- sum(is.na(terra::values(rs[["B04"]], mat = FALSE)))
    fill_prop_val <- fill_count / n_pixels

    # 2. Cloud Proportion via SCL integer classes.
    #    Values 3 (shadow), 8 (medium cloud), 9 (high cloud), 10 (cirrus)
    #    are treated as cloud-contaminated — equivalent to Landsat's high-confidence
    #    cloud bits (bits 8 & 9 of QA_PIXEL).
    scl_vals       <- terra::values(rs[[scl_band_name]], mat = FALSE)
    cloud_pixels   <- sum(scl_vals %in% c(3L, 8L, 9L, 10L), na.rm = TRUE)
    cloud_prop_val <- cloud_pixels / n_pixels

    rm(scl_vals, rs); gc()

    dplyr::tibble(
      file             = image,
      fill_proportion  = fill_prop_val,
      cloud_proportion = cloud_prop_val,
      clear_proportion = 1 - (fill_prop_val + cloud_prop_val)
    )

  }, error = function(e) {
    dplyr::tibble(file = image, fill_proportion = NA, cloud_proportion = NA, clear_proportion = NA)
  })

  return(res_tibble)
}

# Internal: processing date (yyyymmdd) of Sentinel-2 catalogue items, for the file names. The item's
# "created" date when it has one; older items lack it: then the product generation time, the date at
# the end of the product name (S2A_MSIL2A_20170903T105031_N0205_R051_T33WXR_20170903T105029.SAFE) or,
# last, the acquisition date. Vectorised over items.
s2_procdate <- function(created, generation_time, product_uri, datetime) {
  n <- length(datetime)
  val <- function(x) {
    if (is.null(x) || length(x) == 0) return(rep(NA_character_, n))
    x <- as.character(unlist(lapply(x, function(z) if (is.null(z) || length(z) == 0) NA else z[1])))
    x[!nzchar(x)] <- NA
    x
  }
  ymd8 <- function(x) gsub("-", "", substr(x, 1, 10))
  out <- ymd8(val(created))
  gen <- ymd8(val(generation_time))
  pu <- stringr::str_match(val(product_uri), "_([0-9]{8})T[0-9]{6}(\\.SAFE)?$")[, 2]
  acq <- ymd8(val(datetime))
  out[is.na(out)] <- gen[is.na(out)]
  out[is.na(out)] <- pu[is.na(out)]
  out[is.na(out)] <- acq[is.na(out)]
  out
}

# Internal: processing baseline ("05.09") from the item property, or from the product name (N0509)
s2_baseline_prop <- function(baseline, product_uri, n) {
  val <- function(x) {
    if (is.null(x) || length(x) == 0) return(rep(NA_character_, n))
    as.character(unlist(lapply(x, function(z) if (is.null(z) || length(z) == 0) NA else z[1])))
  }
  b <- val(baseline)
  m <- stringr::str_match(val(product_uri), "_N([0-9]{2})([0-9]{2})_")
  from_uri <- ifelse(is.na(m[, 2]), NA_character_, paste0(m[, 2], ".", m[, 3]))
  b[is.na(b) | !nzchar(b)] <- from_uri[is.na(b) | !nzchar(b)]
  b
}

# Internal: output file name of a catalogue item:
# <sat>_MSIL2A_<tile>_<acquisition date>_<processing date>_<time>GMT.tif
s2_item_file <- function(ft) {
  p <- ft$properties
  platform_short <- unname(c("Sentinel-2A" = "S2A", "Sentinel-2B" = "S2B", "Sentinel-2C" = "S2C",
                             "Sentinel-2D" = "S2D")[p$platform])
  tile  <- paste0("T", p$`s2:mgrs_tile`)
  date  <- gsub("-", "", substr(p$datetime, 1, 10))
  date2 <- s2_procdate(p$created, p$`s2:generation_time`, p$`s2:product_uri`, p$datetime)
  time  <- paste0(gsub(":", "", substr(p$datetime, 12, 19)), "GMT.tif")
  paste0(platform_short, "_MSIL2A_", tile, "_", date, "_", date2, "_", time)
}

# Internal: acquisition key of image files: satellite, date, time and whether the scene has the
# +1000 reflectance offset (images with the same key are mosaicked). meta: the s2s.csv table (for the
# processing baseline), or NULL.
s2_acquisition_key <- function(files, meta = NULL) {
  if (!length(files)) return(character(0))
  p <- strsplit(basename(files), "_")
  base <- if (!is.null(meta) && "PROCESSING_BASELINE" %in% names(meta)) {
    meta$PROCESSING_BASELINE[match(files, meta$file)]
  } else rep(NA, length(files))
  off <- vapply(seq_along(files), function(i) s2_needs_offset(files[i], base[i]), logical(1))
  vapply(seq_along(p), function(i) {
    z <- p[[i]]
    paste(z[1], z[4], if (length(z) >= 6) sub("GMT.tif$", "", z[6]) else "", off[i], sep = "_")
  }, character(1))
}

# Internal: proportion of the AOI with valid data (B04) in at least one of the given images
s2_valid_fraction <- function(files, image_dir, aoi) {
  if (!length(files)) return(0)
  terra::terraOptions(memmax = 1)
  r <- terra::rast(lapply(files, function(f) s2_named(terra::rast(file.path(image_dir, f)), f)[["B04"]]))
  v <- if (terra::nlyr(r) > 1) terra::app(!is.na(r), "max") else !is.na(r)
  m <- terra::rasterize(terra::vect(sf::st_transform(aoi, terra::crs(v))), v)
  as.numeric(terra::global(terra::mask(v, m), "mean", na.rm = TRUE)[1, 1])
}
