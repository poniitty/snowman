#' Extracts Sentinel-2 Imagery
#'
#' This function extracts Sentinel-2 L2A imagery that fulfills the specified arguments.
#' It processes the area of interest (AOI), transforms coordinates, and downloads
#' the relevant Sentinel-2 scenes from Microsoft Planetary Computer.
#'
#' @param aoi An `sf` object or a list containing the area of interest.
#' @param site_name A character string representing the name of the site.
#' @param aoi_size A numeric value representing the size of the AOI in kilometers. Ignored if aoi is a POLYGON.
#' @param start_date A character string representing the start date for the imagery.
#' @param end_date A character string representing the end date for the imagery.
#' @param months A numeric vector representing the months to include (default is all months, so 1:12).
#' @param sats A character vector representing the satellite(s) to use (default is both, so c("S2A", "S2B")).
#' @param minclouds A numeric value representing the minimum cloud cover percentage (default is 50).
#' @param base_sentinel2_dir A character string representing the base directory for Sentinel-2 imagery.
#' @param workers A numeric value representing the number of workers for parallel processing.
#' Setting workers > 1 (default) enables parallel computing.
#' @param data_source A character string representing the data source (currently only option is "rstac").
#' @param force Logical, if TRUE deletes all previously downloaded imagery.
#' @return A tibble containing the metadata of the downloaded Sentinel-2 scenes.
#' @export
#' @import rstac dplyr sf terra lubridate stringr parallel tibble readr ggplot2 future future.apply lwgeom
extract_sentinel2 <- function(aoi,
                              site_name,
                              aoi_size = 2,
                              start_date = "2020-01-01",
                              end_date = "2024-12-01",
                              months = 1:12,
                              sats = c("S2A", "S2B"),
                              minclouds = 50,
                              base_sentinel2_dir,
                              workers = 1,
                              data_source = "rstac",
                              force = FALSE) {

  # --- 1. PORTABLE PARALLEL SETUP ---
  is_rstudio_env <- parallelly::supportsMulticore() == FALSE && Sys.info()['sysname'] %in% c("Linux", "Darwin")
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

  } else if (workers > 1 && current_plan_workers > 1) {
    message(paste("Reusing existing parallel plan with", current_plan_workers, "workers."))
  } else {
    future::plan(future::sequential)
  }
  # --- END OF SETUP ---

  # --- 2. INPUT VALIDATION & GEOMETRY FIX ---
  utmall <- utm_zones

  suppressWarnings({
    suppressMessages({
      sf::sf_use_s2(FALSE)

      if (inherits(aoi, "sf")) {
        aoi <- sf::st_make_valid(aoi)

        if (st_geometry_type(aoi) == "POINT") {
          aoi_mid <- aoi %>% st_centroid() %>% st_transform(crs = 4326)
        } else if (st_geometry_type(aoi) == "POLYGON") {
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
  if (as.numeric(st_area(aoi)) > 1000e6) stop("ERROR: AOI > 1000km^2. Process in blocks.")

  # --- 3. DIRECTORY & TERRA SETUP ---
  area_sentinel2_dir <- file.path(base_sentinel2_dir, site_name, "imagery")
  if (!dir.exists(area_sentinel2_dir)) dir.create(area_sentinel2_dir, recursive = TRUE)

  my_temp_dir <- file.path(base_sentinel2_dir, "temp_terra")
  if (!dir.exists(my_temp_dir)) dir.create(my_temp_dir, recursive = TRUE)
  terra::terraOptions(tempdir = my_temp_dir)

  if (force) {
    unlink(list.files(area_sentinel2_dir, pattern = "GMT.tif$", full.names = TRUE))
    unlink(file.path(area_sentinel2_dir, c("s2s.csv", "s2s_final.csv")))
  }

  tifs <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")
  # File name format: S2A_MSIL2A_T33VWF_20210101_20210101_100031GMT.tif
  # Position [4] = acquisition date, matching the Landsat convention
  tif_dates <- if (length(tifs) > 0) {
    tibble(satid = sapply(str_split(tifs, "_"), `[`, 1),
           date  = as.character(ymd(sapply(str_split(tifs, "_"), `[`, 4))))
  } else {
    tibble(satid = NA, date = NA)
  }

  # --- 4. EXTRACTION & PARALLEL PROCESSING ---
  if (data_source == "rstac") {
    s2s <- extract_sentinel2_stac(workers = workers, start_date = start_date, end_date = end_date,
                                  months = months, aoi = aoi, site_name = site_name, epsg = epsg,
                                  excl_dates = tif_dates, sats = sats,
                                  area_sentinel2_dir = area_sentinel2_dir, minclouds = minclouds)
  }
  gc()

  if (is.null(s2s)) {
    s2s <- tibble(id = NULL, DATE_ACQUIRED = NULL)
  } else {
    s2s <- s2s %>% mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
  }

  if (file.exists(file.path(area_sentinel2_dir, "s2s.csv"))) {
    existing_s2s <- read_csv(file.path(area_sentinel2_dir, "s2s.csv"), show_col_types = FALSE) %>%
      mutate(id = as.character(id), DATE_ACQUIRED = as.character(DATE_ACQUIRED))
    s2s <- bind_rows(existing_s2s, s2s) %>% distinct()
  }

  if (nrow(s2s) == 0) {
    print("Zero Sentinel-2 scenes found or downloaded!")
    return(NULL)
  } else {
    write_csv(s2s, file.path(area_sentinel2_dir, "s2s.csv"))

    # Refresh file list
    tifs <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")
    s2s <- full_join(tibble(file = tifs), s2s, by = join_by(file))

    # Metadata parsing from file name parts
    s2s$area  <- site_name
    s2s_parts <- str_split(s2s$file, "_")
    s2s$satid <- sapply(s2s_parts, `[`, 1)                                  # S2A / S2B
    s2s$level <- sapply(s2s_parts, `[`, 2)                                  # MSIL2A
    s2s$tile  <- sapply(s2s_parts, `[`, 3)                                  # MGRS tile, e.g. T33VWF
    s2s$date  <- ymd(sapply(s2s_parts, `[`, 4))                             # acquisition date
    s2s$time  <- gsub(".tif", "", sapply(s2s_parts, `[`, 6))                # UTC time string

    s2s <- s2s %>% arrange(date)

    # --- 5. CHECK RASTERS (FUTURE APPLY) ---
    img_remove <- future.apply::future_lapply(
      s2s$file,
      check_raster,                     # reuse the generic check_raster from the package
      image_dir = area_sentinel2_dir,
      future.seed = TRUE,
      future.packages = c("terra"),
      future.scheduling = 5
    ) %>% unlist()

    gc()

    if (length(img_remove) > 0) {
      unlink(file.path(area_sentinel2_dir, img_remove))
      s2s <- s2s %>% filter(!file %in% img_remove)
      write_csv(s2s, file.path(area_sentinel2_dir, "s2s.csv"))
    }

    tifs <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")
    s2s  <- full_join(tibble(file = tifs), s2s, by = join_by(file))
    s2s  <- s2s %>% arrange(date) %>% filter(!is.na(file))

    # --- 6. MOSAICKING ---
    # Each Sentinel-2 MGRS tile is geometrically unique, so mosaicking normally
    # only arises when the AOI straddles two tiles acquired on the same orbit pass.
    # Group by date + satellite + MGRS tile (instead of Landsat path/row).
    s2s_d <- s2s %>% group_by(date, satid, tile) %>% count() %>% filter(n > 1) %>% ungroup()
    if (nrow(s2s_d) > 0) {
      for (ii in seq_len(nrow(s2s_d))) {
        s2s_dd <- s2s_d %>% slice(ii)
        s2s_dd <- right_join(s2s, s2s_dd, by = join_by(satid, tile, date))

        rs <- lapply(s2s_dd$file, function(x) {
          r <- rast(file.path(area_sentinel2_dir, x))
          r[r == 0] <- NA
          return(r)
        })

        rs <- sprc(rs)
        rs <- mosaic(rs) %>% round()
        writeRaster(rs, file.path(area_sentinel2_dir, s2s_dd$file[[1]]),
                    overwrite = TRUE, datatype = "INT2U")
        unlink(file.path(area_sentinel2_dir, s2s_dd$file[[2:nrow(s2s_dd)]]))
        s2s <- s2s %>% filter(!file %in% s2s_dd$file[[2:nrow(s2s_dd)]])
      }
    }

    # --- 7. CALC COVERAGES (FUTURE APPLY) ---
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

    img_remove <- s2s %>% arrange(desc(fill_proportion)) %>%
      filter(clear_proportion < 0.1) %>% pull(file)

    unlink(file.path(area_sentinel2_dir, img_remove))
    s2s <- s2s %>% filter(!file %in% img_remove)

    write_csv(s2s, file.path(area_sentinel2_dir, "s2s_final.csv"))

    return(s2s)
  }
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

# Internal Function to Extract Sentinel-2 Imagery from STAC
extract_sentinel2_stac <- function(aoi, epsg, excl_dates, site_name,
                                   sats = c("S2A", "S2B"),
                                   start_date, end_date, months = 1:12,
                                   minclouds, area_sentinel2_dir, workers) {

  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

  # Planetary Computer stores platform as "Sentinel-2A" / "Sentinel-2B"
  sat_map <- c("S2A" = "Sentinel-2A", "S2B" = "Sentinel-2B")
  sats2   <- sat_map[sats]

  # Initial STAC search
  it_obj <- s_obj %>%
    stac_search(collections = "sentinel-2-l2a",
                bbox     = st_bbox(aoi %>% st_transform(4326)),
                datetime = paste0(start_date, "/", end_date),
                limit    = 1000) %>%
    get_request()

  # Handle pagination for > 1000 items (identical strategy as Landsat)
  if (length(it_obj$features) == 1000) {
    dr <- tibble(date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day")) %>%
      mutate(gr = cut_number(row_number(.), 20)) %>%
      group_by(gr) %>% group_split()

    it_obj <- s_obj %>%
      stac_search(collections = "sentinel-2-l2a",
                  bbox     = st_bbox(aoi %>% st_transform(4326)),
                  datetime = paste0(min(dr[[1]]$date), "/", max(dr[[1]]$date)),
                  limit    = 1000) %>% get_request()

    for (i in 2:length(dr)) {
      it_obj2 <- s_obj %>%
        stac_search(collections = "sentinel-2-l2a",
                    bbox     = st_bbox(aoi %>% st_transform(4326)),
                    datetime = paste0(min(dr[[i]]$date), "/", max(dr[[i]]$date)),
                    limit    = 1000) %>% get_request()
      it_obj$features <- c(it_obj$features, it_obj2$features)
    }
  }

  # Filter by cloud cover, platform, months, and already-downloaded dates.
  # Band selection: 12 reflectance bands + SCL (Scene Classification Layer = QA equivalent).
  # All bands are resampled to 10 m during download (see process_feature_s2).
  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) { x$properties$`eo:cloud_cover` < minclouds }) %>%
    items_filter(filter_fn = function(x) { x$properties$platform %in% sats2 }) %>%
    assets_select(asset_names = c("B01", "B02", "B03", "B04", "B05", "B06",
                                  "B07", "B08", "B8A", "B09", "B11", "B12", "SCL")) %>%
    items_filter(filter_fn = function(x) { month(ymd_hms(x$properties$datetime)) %in% months }) %>%
    items_filter(filter_fn = function(x) { !as_date(ymd_hms(x$properties$datetime)) %in% excl_dates$date })

  if (length(it_obj$features) > 0) {

    suppressMessages({
      suppressWarnings({
        itst <- items_as_sf(it_obj) %>% st_make_valid()
        itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
      })
    })

    if (nrow(itst) > 0) {
      itst$area <- as.numeric(st_area(itst)) / (1000 * 1000)

      # Duplicate handling: group by date + platform + MGRS tile.
      # (Replaces Landsat's grouping by date + instruments + wrs_path.)
      tt <- itst %>%
        mutate(date = as_date(ymd_hms(datetime))) %>%
        group_by(date, platform, `s2:mgrs_tile`, .add = TRUE) %>%
        mutate(n = n()) %>%
        arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()

      suppressMessages({
        suppressWarnings({
          sf::sf_use_s2(FALSE)
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

      # s2:granule_id is the unique scene identifier (analogous to landsat:scene_id)
      itst <- itst %>%
        filter(`s2:granule_id` %in% (tt %>% pull(`s2:granule_id`)))

      it_obj <- it_obj %>%
        items_filter(filter_fn = function(x) {
          x$properties$`s2:granule_id` %in% (itst %>% pull(`s2:granule_id`))
        })

      # --- CALL PROCESSING IN PARALLEL ---
      juuh <- process_features_in_parallel_s2(it_obj, area_sentinel2_dir, aoi, workers, epsg = epsg)

      # Recode platform to short satellite codes
      itst$platform <- case_when(
        itst$platform == "Sentinel-2A" ~ "S2A",
        itst$platform == "Sentinel-2B" ~ "S2B",
        TRUE ~ itst$platform
      )

      tifs <- list.files(area_sentinel2_dir, pattern = "GMT.tif$")

      # Build file name that matches what process_feature_s2 writes:
      # S2A_MSIL2A_T33VWF_20210101_20210101_100031GMT.tif
      # Note: s2:mgrs_tile on Planetary Computer is stored WITHOUT the leading "T"
      # (e.g. "33VWF"), so we prepend "T" to match the standard tile notation.
      idsALL <- itst %>%
        mutate(date  = gsub("-", "", as_date(ymd_hms(datetime)))) %>%
        mutate(date2 = gsub("-", "", as_date(ymd_hms(created)))) %>%
        mutate(time  = paste0(gsub(":", "", substr(datetime, 12, 19)), "GMT.tif")) %>%
        mutate(file  = paste0(platform, "_MSIL2A_",
                              paste0("T", `s2:mgrs_tile`), "_",
                              date, "_", date2, "_",
                              time)) %>%
        select(file, date, `eo:cloud_cover`, `view:sun_elevation`) %>%
        rename(DATE_ACQUIRED = date,
               CLOUD_COVER   = `eo:cloud_cover`,
               SUN_ELEVATION = `view:sun_elevation`) %>%
        st_drop_geometry() %>%
        filter(file %in% tifs) %>%
        mutate(DATE_ACQUIRED = ymd(DATE_ACQUIRED)) %>%
        rownames_to_column("id")

    } else {
      print("No new Sentinel-2 scenes downloaded!")
      idsALL <- NULL
    }
  } else {
    print("No new Sentinel-2 scenes downloaded!")
    idsALL <- NULL
  }

  return(idsALL)
}

#' Internal Function: Parallel Feature Processing (Sentinel-2) with Grid Alignment
process_features_in_parallel_s2 <- function(it_obj, area_sentinel2_dir, aoi, workers, epsg) {

  # Build grid template from a feature whose native EPSG matches the target.
  # Use B04 (red, 10 m) as the template band — it is always present in L2A
  # and defines the finest native resolution we want for the output stack.
  idx <- which(sapply(it_obj$features, function(x) x$properties$`proj:epsg`) == epsg)
  if (length(idx) == 0) idx <- 1 else idx <- idx[1]

  item_temp <- it_obj$features[[idx]]
  tmp_url   <- make_vsicurl_url_s2(rstac::assets_url(item_temp)[["B04"]])

  r_temp    <- terra::crop(terra::rast(tmp_url), aoi %>% sf::st_transform(crs = epsg))
  exte_vals <- as.numeric(sf::st_bbox(r_temp))[c(1, 3, 2, 4)]
  reso_val  <- terra::res(r_temp)[1]   # 10 m; coarser bands resampled to match

  future.apply::future_lapply(
    it_obj$features,
    process_feature_s2,
    area_sentinel2_dir = area_sentinel2_dir,
    aoi  = aoi,
    exte = exte_vals,
    reso = reso_val,
    future.packages   = c("sf", "terra", "stringr", "rstac"),
    future.seed       = TRUE,
    future.scheduling = 10,
    future.globals    = c("make_vsicurl_url_s2")
  )
}

#' Internal Function: Single Sentinel-2 Feature Extraction with Retry Logic
process_feature_s2 <- function(ft, area_sentinel2_dir, aoi, exte, reso) {

  # Redeclare locally so the function is self-contained in workers
  make_vsicurl_url_s2 <- function(base_url) {
    paste0("/vsicurl", "?pc_url_signing=yes", "&pc_collection=sentinel-2-l2a", "&url=", base_url)
  }

  # Build output file name ---------------------------------------------------
  # Format: S2A_MSIL2A_T33VWF_20210101_20210101_100031GMT.tif
  platform_short <- ifelse(ft$properties$platform == "Sentinel-2A", "S2A", "S2B")
  tile  <- paste0("T", ft$properties$`s2:mgrs_tile`)              # prepend "T" (e.g. T33VWF)
  date  <- gsub("-", "", substr(ft$properties$datetime, 1, 10))   # acquisition date
  date2 <- gsub("-", "", substr(ft$properties$created,  1, 10))   # processing date
  time  <- paste0(gsub(":", "", substr(ft$properties$datetime, 12, 19)), "GMT.tif")

  nm <- paste0(platform_short, "_MSIL2A_", tile, "_", date, "_", date2, "_", time)
  output_file <- file.path(area_sentinel2_dir, nm)
  if (file.exists(output_file)) return(paste0("EXISTS:", nm))
  # --------------------------------------------------------------------------

  all_urls   <- make_vsicurl_url_s2(rstac::assets_url(ft))
  max_retries <- 3

  for (iii in 1:max_retries) {
    e <- try({

      ee <- try(r_remote <- terra::rast(all_urls), silent = TRUE)

      if (inherits(ee, "try-error")) {
        # Identify and drop the failing band URL, pad with NA later
        raw_err       <- as.character(ee)
        extracted_url <- gsub("\\n", "", stringr::str_extract(raw_err, "(https://.+?)\\n"))
        r_remote      <- terra::rast(all_urls[!all_urls %in% make_vsicurl_url_s2(extracted_url)])
      } else {
        extracted_url <- NULL
      }

      # Grid alignment: project all bands to the 10 m UTM template.
      # This uniformly resamples 20 m bands (B05-B07, B8A, B11, B12, SCL) and
      # 60 m bands (B01, B09) to 10 m using nearest-neighbour, matching B02-B04/B08.
      template  <- terra::rast(extent = terra::ext(exte), resolution = reso, crs = sf::st_crs(aoi)$wkt)
      r_cropped <- terra::crop(r_remote, sf::st_transform(aoi, sf::st_crs(r_remote)))
      r_final   <- terra::project(r_cropped, template, method = "near")

      # Pad any band that failed to download with a blank NA layer
      if (length(extracted_url) > 0) {
        for (i_miss in extracted_url) {
          rtemp        <- r_final[[1]]
          rtemp[]      <- NA
          names(rtemp) <- stringr::str_split_i(basename(i_miss), "\\.", 1)
          r_final      <- c(r_final, rtemp)
        }
      }

      # Standardise band names to clean identifiers (B01, B02, … B12, B8A, SCL).
      # rstac asset names are already clean; this regex guards against
      # provider-side changes that may prefix extra path components.
      band_key  <- c("B01", "B02", "B03", "B04", "B05", "B06",
                     "B07", "B08", "B8A", "B09", "B11", "B12", "SCL")
      cur_names <- names(r_final)
      clean_names <- sapply(cur_names, function(n) {
        hit <- band_key[sapply(band_key, function(k) grepl(k, n, ignore.case = FALSE))]
        if (length(hit) == 1) hit else n
      })
      names(r_final) <- clean_names

      # Order: reflectance bands (sorted) then SCL — mirrors the SR_/ST_/QA_ ordering
      # used for Landsat, keeping the QA-equivalent layer last in the stack.
      refl_order <- c("B01", "B02", "B03", "B04", "B05", "B06",
                      "B07", "B08", "B8A", "B09", "B11", "B12")
      refl_bands <- refl_order[refl_order %in% names(r_final)]
      qa_bands   <- names(r_final)[names(r_final) == "SCL"]
      r_final    <- r_final[[c(refl_bands, qa_bands)]]

      terra::writeRaster(r_final, output_file, datatype = "INT2U", overwrite = TRUE, NAflag = 0)

    }, silent = TRUE)

    if (!inherits(e, "try-error")) return(paste0("SUCCESS:", nm))
    if (iii == max_retries) return(paste0("ERROR:", nm, " - ", as.character(e)))
    Sys.sleep(iii * 2)
  }
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
  require(terra)
  require(dplyr)

  fpath <- file.path(image_dir, image)

  res_tibble <- tryCatch({
    rs <- terra::rast(fpath)

    # Identify SCL band (analogous to QA_PIXEL in Landsat)
    scl_band_name <- names(rs)[grepl("^SCL$", names(rs), ignore.case = TRUE)]
    if (length(scl_band_name) == 0) stop("SCL band not found")

    n_pixels <- terra::ncell(rs)

    # 1. Fill Proportion: pixels where B04 == 0 indicate no-data fill areas.
    #    (Mirrors the Landsat approach of checking band 4 for zero values.)
    fill_tab      <- terra::freq(rs[["B04"]], value = 0)
    fill_count    <- if (nrow(fill_tab) > 0) fill_tab[1, "count"] else 0
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
