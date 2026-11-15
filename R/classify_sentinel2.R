# ==========================================================================
# Sentinel-2 classification: predictors, the feature stack shared by model
# training and classification, and the classification itself.
#
# The workflow mirrors the Landsat one:
#   extract_sentinel2() -> calc_predictors_s2() -> classify_sentinel2() -> calc_snow_variables()
# The classifications have the same layers as the Landsat ones (class, land, water, snow,
# cloud, artif), at 10 m, so calc_snow_variables() works on them unchanged.
#
# The features are built by s2_feature_stack() (rasters) and s2_prep_predictors()
# (band ratios etc. on the pixel table). Training data are extracted with the same two
# functions, so training and classification use identical features.
# ==========================================================================

# Layers of the Sentinel-2 images written by extract_sentinel2()
s2_band_key <- c("B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B11", "B12", "SCL")
s2_refl <- s2_band_key[1:12]
# Bands with focal statistics (Landsat: B2 blue, B5 NIR, B7 SWIR2)
s2_focal_bands <- c("B02", "B08", "B12")
# ESA WorldCover classes (fixed factor levels, so that training and prediction agree)
s2_esa_levels <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100)
# Processing baseline 04.00 (from 25 Jan 2022) adds 1000 to the reflectance values
s2_offset_date <- as.Date("2022-01-25")
# SCL classes behind the scene-classification probabilities (cf. the Landsat QA bits)
s2_scl_classes <- list(cloud = c(8, 9), cirrus = 10, cshadow = 3, snow = 11, water = 6,
                       clear = c(2, 4, 5, 6, 11), dark = 2, unclass = 7)

#' Read a Sentinel-2 image, harmonised to the reflectance scale of processing baselines before 04.00
#'
#' Reads an image written by [extract_sentinel2()] and returns its 12 reflectance bands and the SCL
#' layer. Pixels without data (SCL 0 or missing) are NA. For scenes of processing baseline 04.00 and
#' later, 1000 is subtracted from the reflectance values (minimum 1), so that all scenes share the
#' scale of the earlier baselines (reflectance x 10000). The baseline is taken from `baseline` (the
#' `PROCESSING_BASELINE` column of the scene table) when known, otherwise from the acquisition date
#' (baseline 04.00 was introduced on 25 January 2022), as recommended by the Planetary Computer.
#'
#' @param file File name of the image (e.g. `S2A_MSIL2A_T33WXR_20170903_20210512_105031GMT.tif`).
#' @param image_dir Directory of the image.
#' @param baseline Processing baseline of the scene (e.g. "05.09" or 5.09), or NA.
#' @return A SpatRaster with layers B01 ... B12, B8A and SCL.
#' @export
read_sentinel2 <- function(file, image_dir, baseline = NA) {
  r <- terra::rast(file.path(image_dir, file))
  if (!all(s2_band_key %in% names(r))) {
    if (terra::nlyr(r) != length(s2_band_key)) stop("Unexpected layers in ", file)
    names(r) <- s2_band_key
  }
  r <- r[[s2_band_key]]
  scl <- r[["SCL"]]
  refl <- r[[s2_refl]]
  if (s2_needs_offset(file, baseline)) {
    refl <- terra::clamp(refl - 1000, lower = 1, values = TRUE)
  }
  r <- c(refl, scl)
  r <- terra::mask(r, scl, maskvalues = c(NA, 0))
  names(r) <- s2_band_key
  r
}

# Internal: does the scene carry the +1000 offset of processing baseline 04.00?
s2_needs_offset <- function(file, baseline = NA) {
  b <- suppressWarnings(as.numeric(as.character(baseline)))
  if (length(b) == 1 && !is.na(b)) return(b >= 4)
  d <- lubridate::ymd(strsplit(basename(file), "_")[[1]][4])
  isTRUE(d >= s2_offset_date)
}

# Internal: processing baseline of a scene from the scene table (NA if not recorded)
s2_baseline <- function(image_df, file) {
  if (!"PROCESSING_BASELINE" %in% names(image_df)) return(NA)
  b <- image_df$PROCESSING_BASELINE[match(file, image_df$file)]
  if (length(b) == 0) NA else b
}

#' Calculate Predictors for the Sentinel-2 Classifier
#'
#' The Sentinel-2 counterpart of [calc_predictors()]. Writes into `<base_sentinel2_dir>/<site_name>/predictors`:
#' \itemize{
#'   \item `ALOSDEM30.tif` (ALOS DEM, 30 m, for the hillshade of each image), `ALOSDEM.tif` (DEM at 10 m),
#'   `slope.tif` (degrees x 100, 10 m)
#'   \item `ESALC.tif` (ESA WorldCover 2021, native 10 m) and `ESAW.tif` (proportion of water in the
#'   surrounding 3 x 3 pixels, i.e. 30 m, as in the Landsat version)
#'   \item `medians2_<month>.tif`: monthly median images (12 reflectance bands; SCL classes
#'   `median_mask` masked) of the up to `max_per_month` clearest scenes of each month, completed over
#'   months and space as in the Landsat version
#'   \item `medianindices.tif`: median NDVI, NDSI, KBRI, SWM and BITM of the up to `n_indices`
#'   clearest scenes per month. As in the Landsat version, clouds are not masked (cloud masks often
#'   take late snow beds for clouds); only snow is masked for NDVI, KBRI, SWM and BITM. The median
#'   (Landsat: mean, except BITM) limits the effect of the clouds left in the images.
#' }
#' All rasters are on the 10 m grid of the images.
#'
#' @param image_df The scene table, the output of [extract_sentinel2()] (or `imagery/s2s_final.csv`).
#' @param site_name Name of the site.
#' @param base_sentinel2_dir Base directory of the Sentinel-2 data.
#' @param max_per_month Maximum number of scenes per month for the median images (default 15).
#' @param n_indices Maximum number of scenes per month for the median indices (default 5, as in the training data of the Landsat models).
#' @param median_mask SCL classes masked in the median images (default cloud shadow 3, clouds 8 and 9,
#' cirrus 10; the Landsat version masks the QA cloud bit). NULL masks nothing.
#' @param index_mask SCL classes masked in the median indices (default NULL: none, as for Landsat).
#' @param force Logical, if TRUE all predictors are calculated again; otherwise the DEM/land cover part
#' and the median part are skipped when their files exist.
#' @param memmax Maximum memory (GB) terra may use before it works through temporary files (default 4;
#' 1 GB per parallel worker). terra judges the free memory from the whole node, not from the job's
#' memory limit, so without a limit it can exceed the memory of a Slurm job or RStudio session and
#' crash R.
#' @param workers Number of parallel workers for the per-scene and per-band steps (default 1). Plan
#' about 2 GB of memory per worker.
#' @return Invisibly, the predictor directory.
#' @export
calc_predictors_s2 <- function(image_df, site_name, base_sentinel2_dir, max_per_month = 15,
                               n_indices = 5, median_mask = c(3, 8, 9, 10), index_mask = NULL,
                               force = FALSE, memmax = 4, workers = 1) {
  old_mem <- s2_terra_memory(memmax)
  on.exit(terra::terraOptions(memmax = old_mem), add = TRUE)
  # parallel map over scenes / bands (own workers; the caller's plan is restored afterwards)
  pmap <- function(X, FUN, ...) {   # FUN: a package function (no closure over this frame)
    if (workers <= 1 || length(X) <= 1) return(lapply(X, FUN, ...))
    future.apply::future_lapply(X, FUN, ..., future.seed = TRUE, future.packages = c("terra", "snowman"),
                                future.globals = FALSE)
  }
  if (workers > 1) {
    old_plan <- future::plan()
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  image_dir <- file.path(base_sentinel2_dir, site_name, "imagery")
  predictor_dir <- file.path(base_sentinel2_dir, site_name, "predictors")
  if (!dir.exists(predictor_dir)) dir.create(predictor_dir, recursive = TRUE)

  image_df <- image_df %>%
    dplyr::filter(file %in% list.files(image_dir, pattern = "GMT.tif$")) %>%
    dplyr::mutate(date = lubridate::as_date(date))
  if (nrow(image_df) == 0) stop("No Sentinel-2 images of ", site_name, " in ", image_dir)

  tmpl <- terra::rast(terra::rast(file.path(image_dir, image_df$file[1]))[[1]])

  # ---- DEM, slope, ESA WorldCover --------------------------------------------------------------
  static_files <- file.path(predictor_dir, c("ALOSDEM30.tif", "ALOSDEM.tif", "slope.tif", "ESALC.tif", "ESAW.tif"))
  if (force || !all(file.exists(static_files))) {
    s2_static_predictors(tmpl, predictor_dir)
  }

  # ---- monthly median images -------------------------------------------------------------------
  med_df <- image_df %>%
    dplyr::filter(cloud_proportion < 0.8, fill_proportion < 0.8) %>%
    dplyr::mutate(month = lubridate::month(date))
  if (nrow(med_df) == 0) stop("No scenes with less than 80% cloud and fill for the median images")
  months <- sort(unique(med_df$month))
  med_files <- file.path(predictor_dir, paste0("medians2_", months, ".tif"))

  if (force || !all(file.exists(med_files))) {
    unlink(list.files(predictor_dir, pattern = "^medians2_[0-9]+\\.tif$", full.names = TRUE))
    tmp_dir <- file.path(predictor_dir, "tmp_medians")
    dir.create(tmp_dir, showWarnings = FALSE)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

    meds <- list()
    for (mo in months) {
      mdf <- med_df %>%
        dplyr::filter(month == mo) %>%
        dplyr::arrange(dplyr::desc(clear_proportion)) %>%
        dplyr::slice_head(n = max_per_month)
      # masked, harmonised reflectance of each image to a temporary file (one pass per image)
      fs <- file.path(tmp_dir, sprintf("m%02d_%03d.tif", mo, seq_len(nrow(mdf))))
      bl <- vapply(mdf$file, function(f) as.character(s2_baseline(mdf, f)), character(1))
      invisible(pmap(seq_len(nrow(mdf)), s2_task_reflectance, files = mdf$file, bl = bl,
                     image_dir = image_dir, mask = median_mask, out = fs))
      message("median image of month ", mo, " (", length(fs), " scenes)")
      bfs <- file.path(tmp_dir, sprintf("b%02d_%s.tif", mo, s2_refl))
      invisible(pmap(seq_along(s2_refl), s2_task_band_median, fs = fs, bands = s2_refl, out = bfs))
      med <- terra::rast(bfs)
      names(med) <- s2_refl
      med <- med * 1   # into memory / its own file before the inputs are removed
      unlink(c(fs, bfs))
      med <- s2_fill_na(med)   # all-NA layers -> 0, small gaps from neighbours
      f <- file.path(tmp_dir, paste0("med_", mo, ".tif"))
      terra::writeRaster(round(med), f, datatype = "INT2U", overwrite = TRUE)
      meds[[as.character(mo)]] <- f
    }

    # Gaps that remain are filled from the neighbouring months (linear interpolation over the months,
    # as in the Landsat version; the first and last month are carried to the ends), then from the
    # neighbouring pixels. Only pixels with a gap are interpolated (vectorised; terra::approximate()
    # on the whole stack took minutes per band).
    has_na <- vapply(meds, function(f) any(is.na(terra::values(terra::rast(f), mat = FALSE))), logical(1))
    if (any(has_na) && length(months) > 1) {
      for (b in s2_refl) {
        st <- terra::rast(lapply(meds, function(f) terra::rast(f)[[b]]))
        M <- terra::values(st, mat = TRUE)
        M <- s2_interp_months(M, months)
        for (k in seq_along(months)) {
          x <- terra::rast(st[[k]]); terra::values(x) <- M[, k]
          terra::writeRaster(x, file.path(tmp_dir, sprintf("fill_%02d_%s.tif", months[k], b)),
                             datatype = "INT2U", overwrite = TRUE)
        }
        rm(M, st); gc(verbose = FALSE)
      }
      message("gaps filled from neighbouring months")
    }
    for (k in seq_along(months)) {
      if (has_na[k] && length(months) > 1) {
        mk <- terra::rast(file.path(tmp_dir, sprintf("fill_%02d_%s.tif", months[k], s2_refl)))
      } else {
        mk <- terra::rast(meds[[k]])
      }
      names(mk) <- s2_refl
      mk <- s2_fill_na(mk, complete = TRUE)
      terra::writeRaster(round(mk), file.path(predictor_dir, paste0("medians2_", months[k], ".tif")),
                         datatype = "INT2U", overwrite = TRUE)
      rm(mk); gc(verbose = FALSE)
    }
  }

  # ---- median indices --------------------------------------------------------------------------
  mi_file <- file.path(predictor_dir, "medianindices.tif")
  if (force || !file.exists(mi_file)) {
    idf <- med_df %>%
      dplyr::group_by(month) %>%
      dplyr::arrange(dplyr::desc(clear_proportion)) %>%
      dplyr::slice_head(n = n_indices) %>%
      dplyr::ungroup()
    tmp_dir_i <- file.path(predictor_dir, "tmp_indices")
    dir.create(tmp_dir_i, showWarnings = FALSE)
    on.exit(unlink(tmp_dir_i, recursive = TRUE), add = TRUE)
    message("median indices from ", nrow(idf), " scenes")
    ifs <- file.path(tmp_dir_i, sprintf("i_%03d.tif", seq_len(nrow(idf))))   # on disk, not in memory
    bl <- vapply(idf$file, function(f) as.character(s2_baseline(idf, f)), character(1))
    invisible(pmap(seq_len(nrow(idf)), s2_task_indices, files = idf$file, bl = bl,
                   image_dir = image_dir, mask = index_mask, out = ifs))
    rr <- lapply(ifs, terra::rast)
    pick <- function(nm, fun) {
      st <- terra::rast(lapply(rr, function(x) x[[nm]]))
      s2_fill_na(terra::app(st, fun, na.rm = TRUE), complete = TRUE)
    }
    mi <- c(pick("ndvi", "median"), pick("ndsi", "median"), pick("kbri", "median"),
            pick("swm", "median"), pick("bitm", "median"))
    names(mi) <- c("ndvi", "ndsi", "kbri", "swm", "bitm")
    terra::writeRaster(round(mi, 4), mi_file, overwrite = TRUE)
  }
  invisible(predictor_dir)
}

# Internal: limit terra's memory use (GB); returns the previous setting. terra estimates the free
# memory from the whole node, which on a shared HPC node is far more than the job may use.
s2_terra_memory <- function(memmax) {
  old <- terra::terraOptions(print = FALSE)$memmax
  if (is.null(old) || length(old) == 0) old <- -1
  if (!is.null(memmax) && is.finite(memmax) && memmax > 0) terra::terraOptions(memmax = memmax)
  old
}

# Internal: DEM, slope and ESA WorldCover on the 10 m image grid
s2_static_predictors <- function(tmpl, predictor_dir) {
  aoi <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(tmpl)))
  aoi30 <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(sf::st_buffer(aoi, 100))))
  s_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

  # ALOS DEM: one virtual mosaic, warped in one step to 30 m (as in calc_predictors())
  it_obj <- stac_with_retry(function() {
    s_obj %>%
      rstac::stac_search(collections = "alos-dem", bbox = as.numeric(sf::st_bbox(sf::st_transform(aoi30, 4326))),
                         limit = 1000) %>%
      rstac::get_request()
  }, label = "STAC search alos-dem")
  dem_urls <- unlist(lapply(it_obj$features, function(ft) {
    u <- make_vsicurl_url_dem(sort(rstac::assets_url(ft)))
    u[endsWith(u, "_DSM.tif")]
  }))
  if (length(dem_urls) == 0) stop("No ALOS DEM tiles found")
  vrt_file <- file.path(predictor_dir, "alos_dsm.vrt")
  if (file.exists(vrt_file)) unlink(vrt_file)
  sf::gdal_utils("buildvrt", source = dem_urls, destination = vrt_file, options = c("-resolution", "highest"))
  dem_file <- file.path(predictor_dir, "alos_dsm_warped.tif")
  sf::gdal_utils("warp", source = vrt_file, destination = dem_file,
                 options = c("-t_srs", terra::crs(tmpl), "-te", sf::st_bbox(aoi30), "-tr", c(30, 30),
                             "-r", "bilinear", "-dstnodata", "-9999", "-overwrite"))
  dem30 <- terra::rast(dem_file) * 1
  dem30[dem30 < -500] <- NA
  for (i in 1:20) {
    if (!any(is.na(terra::values(dem30, mat = FALSE)))) break
    dem30 <- terra::focal(dem30, 3, "mean", na.policy = "only", na.rm = TRUE)
  }
  dem30[dem30 < 0] <- 0
  dem30[is.na(dem30)] <- 0
  names(dem30) <- "dem"
  unlink(c(vrt_file, dem_file))
  terra::writeRaster(dem30, file.path(predictor_dir, "ALOSDEM30.tif"), overwrite = TRUE)

  # ESA WorldCover 2021, warped tile by tile onto the image grid (10 m), mosaicked with max
  it_obj <- stac_with_retry(function() {
    s_obj %>%
      rstac::stac_search(collections = "esa-worldcover", bbox = as.numeric(sf::st_bbox(sf::st_transform(aoi30, 4326))),
                         datetime = stac_datetime("2021-01-01", "2021-12-31"), limit = 1000) %>%
      rstac::get_request()
  }, label = "STAC search esa-worldcover")
  esa_urls <- unlist(lapply(it_obj$features, function(ft) {
    u <- make_vsicurl_url_esa(sort(rstac::assets_url(ft)))
    u[endsWith(u, "_Map.tif")]
  }))
  e <- terra::ext(tmpl)
  esas <- lapply(seq_along(esa_urls), function(k) {
    f <- file.path(predictor_dir, paste0("esa_", k, "_Map.tif"))
    ok <- try(sf::gdal_utils("warp", source = esa_urls[k], destination = f,
                             options = c("-t_srs", terra::crs(tmpl),
                                         "-te", e$xmin, e$ymin, e$xmax, e$ymax,
                                         "-tr", terra::res(tmpl), "-r", "near", "-overwrite")),
              silent = TRUE)
    if (inherits(ok, "try-error")) return(NULL)
    x <- terra::rast(f) * 1
    x[x == 0] <- NA
    x
  })
  esas <- Filter(Negate(is.null), esas)
  if (length(esas) == 0) stop("No ESA WorldCover tiles could be read")
  esa <- if (length(esas) > 1) terra::mosaic(terra::sprc(esas), fun = "max") else esas[[1]]
  esa <- terra::resample(esa, tmpl, method = "near")   # exact grid (no-op when already aligned)
  esa[esa < 0] <- 0
  esa[is.na(esa)] <- 0
  names(esa) <- "esalc"
  unlink(list.files(predictor_dir, pattern = "_Map.tif$", full.names = TRUE))
  terra::writeRaster(esa, file.path(predictor_dir, "ESALC.tif"), overwrite = TRUE, datatype = "INT1U")

  esaw <- terra::focal(esa == 80, 3, "mean", na.rm = TRUE)   # 30 m water proportion, as for Landsat
  names(esaw) <- "esaw"
  terra::writeRaster(round(esaw, 3), file.path(predictor_dir, "ESAW.tif"), overwrite = TRUE)

  # Slope from the 30 m DEM (no facets from an interpolated DEM), interpolated to 10 m
  slp <- terra::project(terra::terrain(dem30, "slope", unit = "degrees"), tmpl, method = "bilinear")
  slp[esaw > 0.2] <- 0
  slp[is.na(slp)] <- 0
  names(slp) <- "slope"
  terra::writeRaster(round(slp * 100), file.path(predictor_dir, "slope.tif"), overwrite = TRUE, datatype = "INT2U")

  dem10 <- terra::project(dem30, tmpl, method = "bilinear")
  names(dem10) <- "dem"
  terra::writeRaster(dem10, file.path(predictor_dir, "ALOSDEM.tif"), overwrite = TRUE)
  invisible(TRUE)
}

# Internal: fill NA cells of each layer from the neighbouring cells. All-NA layers become 0.
# complete = TRUE fills until no NA is left (cells still empty after 100 rounds get the layer median).
s2_fill_na <- function(x, complete = FALSE) {
  layers <- lapply(seq_len(terra::nlyr(x)), function(k) {
    y <- x[[k]]
    v <- terra::values(y, mat = FALSE)
    if (all(is.na(v))) { y[] <- 0; return(y) }
    if (!any(is.na(v))) return(y)
    for (i in seq_len(if (complete) 100 else 1)) {
      y <- terra::focal(y, 5, "mean", na.policy = "only", na.rm = TRUE)
      if (!any(is.na(terra::values(y, mat = FALSE)))) break
    }
    if (complete) y[is.na(y)] <- stats::median(v, na.rm = TRUE)
    y
  })
  out <- terra::rast(layers)
  names(out) <- names(x)
  out
}

# Internal: the 12 reflectance bands of a scene, harmonised (baseline 04.00 offset removed, minimum 1)
# and with no-data pixels and the SCL classes `mask_classes` set to NA, written to `out`. One pass over
# the image (terra::app on blocks).
s2_masked_reflectance <- function(file, image_dir, baseline, mask_classes, out) {
  s2_worker_limits()
  r <- terra::rast(file.path(image_dir, file))
  if (!all(s2_band_key %in% names(r))) names(r) <- s2_band_key
  off <- if (s2_needs_offset(file, baseline)) 1000 else 0
  mc <- mask_classes
  terra::app(r[[s2_band_key]], fun = function(m) {
    scl <- m[, 13]
    v <- m[, 1:12, drop = FALSE]
    if (off > 0) v[] <- pmax(v - off, 1)
    v[is.na(scl) | scl == 0 | scl %in% mc, ] <- NA
    v
  }, filename = out, overwrite = TRUE, wopt = list(datatype = "INT2U", names = s2_refl))
  invisible(out)
}

# Internal: NDVI, NDSI, KBRI, SWM and BITM of a scene (reflectance 0-1, harmonised; snow (SCL 11)
# masked except for NDSI; SCL classes `mask_classes` masked), written to `out`. One pass over the image.
s2_scene_indices <- function(file, image_dir, baseline, mask_classes, out) {
  s2_worker_limits()
  r <- terra::rast(file.path(image_dir, file))
  if (!all(s2_band_key %in% names(r))) names(r) <- s2_band_key
  off <- if (s2_needs_offset(file, baseline)) 1000 else 0
  mc <- mask_classes
  terra::app(r[[c("B02", "B03", "B04", "B08", "B11", "SCL")]], fun = function(m) {
    scl <- m[, 6]
    v <- m[, 1:5, drop = FALSE]
    if (off > 0) v[] <- pmax(v - off, 1)
    v <- v / 10000
    v[is.na(scl) | scl == 0 | scl %in% mc, ] <- NA
    snow <- scl %in% 11
    b2 <- v[, 1]; b3 <- v[, 2]; b4 <- v[, 3]; b8 <- v[, 4]; b11 <- v[, 5]
    clamp1 <- function(x) pmin(pmax(x, -1), 1)
    ndvi <- clamp1((b8 - b4) / (b8 + b4)); ndvi[snow] <- NA
    ndsi <- clamp1((b3 - b11) / (b3 + b11))
    kbri <- clamp1((b11 - b8) / (20 * sqrt(b11 + b8))); kbri[snow] <- NA
    swm <- (b2 + b3) / (b8 + b11); swm[snow | abs(swm) > 10] <- NA
    bitm <- sqrt((b2^2 + b3^2 + b4^2) / 3); bitm[snow] <- NA
    cbind(ndvi, ndsi, kbri, swm, bitm)
  }, filename = out, overwrite = TRUE, wopt = list(names = c("ndvi", "ndsi", "kbri", "swm", "bitm")))
  invisible(out)
}

# Internal: tasks of calc_predictors_s2() for the parallel workers
s2_task_reflectance <- function(k, files, bl, image_dir, mask, out) {
  s2_masked_reflectance(files[k], image_dir, bl[k], mask, out[k])
}
s2_task_indices <- function(k, files, bl, image_dir, mask, out) {
  s2_scene_indices(files[k], image_dir, bl[k], mask, out[k])
}
s2_task_band_median <- function(i, fs, bands, out) {
  s2_worker_limits()
  st <- terra::rast(lapply(fs, function(f) terra::rast(f)[[bands[i]]]))
  terra::app(st, "median", na.rm = TRUE, filename = out[i], overwrite = TRUE)
  out[i]
}

# Internal: memory limits of terra and GDAL in a (worker) process, unless set already
s2_worker_limits <- function(memmax = 1) {
  if (terra::terraOptions(print = FALSE)$memmax <= 0) terra::terraOptions(memmax = memmax)
  if (!nzchar(Sys.getenv("GDAL_CACHEMAX"))) Sys.setenv(GDAL_CACHEMAX = "256")
  invisible(NULL)
}

# Internal: linear interpolation of the gaps (NA) in each row of M over the times tt (columns); values
# before the first / after the last observation are carried to the ends (approx(rule = 2)).
s2_interp_months <- function(M, tt) {
  rows <- which(rowSums(is.na(M)) > 0 & rowSums(!is.na(M)) > 0)
  if (!length(rows)) return(M)
  X <- M[rows, , drop = FALSE]
  n <- nrow(X); nt <- ncol(X)
  prev_i <- matrix(NA_integer_, n, nt); next_i <- matrix(NA_integer_, n, nt)
  last <- rep(NA_integer_, n)
  for (j in seq_len(nt)) { ok <- !is.na(X[, j]); last[ok] <- j; prev_i[, j] <- last }
  last <- rep(NA_integer_, n)
  for (j in rev(seq_len(nt))) { ok <- !is.na(X[, j]); last[ok] <- j; next_i[, j] <- last }
  na <- which(is.na(X))
  r <- (na - 1) %% n + 1
  cc <- (na - 1) %/% n + 1
  p <- prev_i[na]; q <- next_i[na]
  vp <- X[cbind(r, ifelse(is.na(p), 1L, p))]; vq <- X[cbind(r, ifelse(is.na(q), 1L, q))]
  w <- (tt[cc] - tt[ifelse(is.na(p), 1L, p)]) / (tt[ifelse(is.na(q), 1L, q)] - tt[ifelse(is.na(p), 1L, p)])
  X[na] <- ifelse(is.na(p), vq, ifelse(is.na(q), vp, vp + (vq - vp) * w))
  M[rows, ] <- X
  M
}

#' Sentinel-2 Feature Stack for the Classifier
#'
#' Builds the raster features of one Sentinel-2 image, the counterparts of the Landsat classifier's
#' predictors. Used by [classify_sentinel2()] and for extracting training data, so that both use the
#' same features. Requires the predictors of [calc_predictors_s2()].
#'
#' Layers (Landsat counterparts in brackets):
#' \itemize{
#'   \item B01 ... B12, B8A (harmonised reflectance x 10000) and SCL (QA; not a model predictor)
#'   \item `B02_`, `B08_`, `B12_` + `tpi`, `min5`, `max5`, `std5`: 5 x 5 pixel statistics (B2, B5, B7)
#'   \item `mr_diff_<band>`: monthly median image minus the image; `mr_diff_B08_std5`, `mr_diff_B12_std5`
#'   \item `slope`, `hillshade` (sun position at the acquisition time), `esalc`, `esaw`
#'   \item `<class>_probs5`: proportion of SCL cloud (8, 9), cirrus (10), cloud shadow (3), snow (11),
#'   water (6), clear (2, 4, 5, 6, 11), dark (2) and unclassified (7) pixels within 150 m, and
#'   `cloud_probs25` within ~750 m, the same ground distances as the Landsat QA-bit proportions (5 and 25
#'   pixels of 30 m). They are calculated at 30 m and 90 m and interpolated to 10 m.
#'   \item `ndvi`, `ndsi`, `kbri`, `swm`, `bitm`: median indices
#' }
#' The normalised differences between bands are added to the pixel table by [s2_prep_predictors()].
#'
#' @param file File name of the image.
#' @param image_dir Directory of the image.
#' @param predictor_dir Directory of the predictors of [calc_predictors_s2()].
#' @param baseline Processing baseline of the scene, or NA (see [read_sentinel2()]).
#' @param vars Optional names of the model predictors: only the layers needed for them are calculated.
#' NULL (default) calculates all.
#' @return A SpatRaster.
#' @export
s2_feature_stack <- function(file, image_dir, predictor_dir, baseline = NA, vars = NULL) {
  need <- function(...) {
    is.null(vars) || any(vapply(c(...), function(v) any(grepl(v, vars, fixed = TRUE)), logical(1)))
  }
  r <- read_sentinel2(file, image_dir, baseline)
  tmpl <- terra::rast(r[[1]])
  out <- list(r)
  pfile <- function(x) {
    p <- file.path(predictor_dir, x)
    if (!file.exists(p)) stop(x, " not found in ", predictor_dir, ". Run calc_predictors_s2() first.")
    y <- terra::rast(p)
    if (!terra::compareGeom(y, tmpl, stopOnError = FALSE)) {
      y <- terra::resample(y, tmpl, method = if (x %in% c("ESALC.tif")) "near" else "bilinear")
    }
    y
  }

  # 5 x 5 statistics (the same as for Landsat)
  ring <- matrix(1, 5, 5)
  ring[13] <- 0
  for (b in s2_focal_bands) {
    if (need(paste0(b, "_tpi"))) {
      x <- r[[b]] - terra::focal(r[[b]], w = ring, fun = mean, na.policy = "omit", na.rm = TRUE)
      names(x) <- paste0(b, "_tpi"); out <- c(out, x)
    }
    for (st in c("min", "max", "sd")) {
      nm <- paste0(b, "_", if (st == "sd") "std" else st, "5")
      if (need(nm)) {
        x <- terra::focal(r[[b]], w = 5, fun = st, na.policy = "omit", na.rm = TRUE)
        names(x) <- nm; out <- c(out, x)
      }
    }
  }

  # Difference to the median image of the same (or nearest) month
  mo <- lubridate::month(lubridate::ymd(strsplit(basename(file), "_")[[1]][4]))
  mr_files <- list.files(predictor_dir, pattern = "^medians2_[0-9]+\\.tif$")
  if (length(mr_files) == 0) stop("No medians2_*.tif files in ", predictor_dir, ". Run calc_predictors_s2() first.")
  avail <- as.numeric(gsub("medians2_|\\.tif", "", mr_files))
  mo_use <- avail[which.min(pmin(abs(avail - mo), 12 - abs(avail - mo)))]
  mr <- pfile(paste0("medians2_", mo_use, ".tif"))
  names(mr) <- s2_refl
  if (need("mr_diff")) {
    md <- mr - r[[s2_refl]]
    names(md) <- paste0("mr_diff_", s2_refl)
    out <- c(out, md)
    for (b in c("B08", "B12")) {
      nm <- paste0("mr_diff_", b, "_std5")
      if (need(nm)) {
        x <- terra::focal(r[[b]], w = 5, fun = "sd") - terra::focal(mr[[b]], w = 5, fun = "sd")
        names(x) <- nm; out <- c(out, x)
      }
    }
  }

  # Terrain and land cover
  esaw <- pfile("ESAW.tif"); names(esaw) <- "esaw"
  if (need("esaw")) out <- c(out, esaw)
  if (need("esalc")) { x <- pfile("ESALC.tif"); names(x) <- "esalc"; out <- c(out, x) }
  if (need("slope")) { x <- pfile("slope.tif") / 100; names(x) <- "slope"; out <- c(out, x) }
  if (need("hillshade")) out <- c(out, s2_hillshade(file, predictor_dir, tmpl, esaw))

  # Scene-classification proportions in the surroundings
  scl <- r[["SCL"]]
  agg_prop <- function(cls, fact, radius) {
    ind <- terra::ifel(is.na(scl), NA, scl %in% cls)
    a <- terra::aggregate(ind, fact, "mean", na.rm = TRUE)
    fm <- terra::focalMat(terra::rast(ncols = 180, nrows = 180, xmin = 0), radius, "circle")
    fm[fm > 0] <- 1
    a <- terra::focal(a, w = fm, fun = "mean", na.policy = "omit", na.rm = TRUE)
    terra::resample(a, tmpl, method = "bilinear")
  }
  for (k in names(s2_scl_classes)) {
    nm <- paste0(k, "_probs5")
    if (need(nm)) { x <- agg_prop(s2_scl_classes[[k]], 3, 5); names(x) <- nm; out <- c(out, x) }
  }
  if (need("cloud_probs25")) { x <- agg_prop(s2_scl_classes$cloud, 9, 8); names(x) <- "cloud_probs25"; out <- c(out, x) }

  # Median indices
  if (need("ndvi", "ndsi", "kbri", "swm", "bitm")) out <- c(out, pfile("medianindices.tif"))

  terra::rast(out)
}

# Internal: hillshade with the sun position at the acquisition time, from the 30 m DEM
s2_hillshade <- function(file, predictor_dir, tmpl, esaw) {
  dem <- terra::rast(file.path(predictor_dir, "ALOSDEM30.tif"))
  p <- strsplit(gsub("GMT.tif$", "", basename(file)), "_")[[1]]
  time <- lubridate::ymd_hms(paste(p[4], p[6]))
  cen <- sf::st_coordinates(sf::st_transform(sf::st_centroid(sf::st_as_sfc(sf::st_bbox(tmpl))), 4326))
  sa <- solartime::computeSunPositionDoyHour(doy = lubridate::yday(time),
                                             hour = lubridate::hour(time) + lubridate::minute(time) / 60,
                                             latDeg = cen[1, "Y"], longDeg = cen[1, "X"], timeZone = 0L)
  hill <- terra::shade(terra::terrain(dem, "slope", unit = "radians"),
                       terra::terrain(dem, "aspect", unit = "radians"),
                       sa[3] * (180 / pi), sa[4] * (180 / pi))
  hill <- terra::project(hill, tmpl, method = "bilinear")
  hill[esaw > 0.2] <- stats::median(terra::values(hill, mat = FALSE), na.rm = TRUE)
  hill <- s2_fill_na(hill, complete = TRUE)
  names(hill) <- "hillshade"
  hill
}

#' Prepare the Sentinel-2 Pixel Table for the Classifier
#'
#' Adds the normalised differences between all pairs of the reflectance bands and the B08/B12 focal
#' statistics (as in the Landsat classifier), replaces NaN with 0 and missing focal statistics with 0,
#' and makes `esalc` a factor with fixed levels. Used both for model training and in
#' [classify_sentinel2()].
#'
#' @param d A data frame of the layers of [s2_feature_stack()] (pixels or training points).
#' @param vars Optional names of the model predictors: only the band pairs among them are added.
#' @return The data frame with the added columns.
#' @export
s2_prep_predictors <- function(d, vars = NULL) {
  d <- as.data.frame(d)
  pr <- s2_pairs(names(d))
  if (!is.null(vars)) pr <- pr[pr$name %in% vars, , drop = FALSE]
  for (k in seq_len(nrow(pr))) {
    a <- d[[pr$a[k]]]
    b <- d[[pr$b[k]]]
    d[[pr$name[k]]] <- (a - b) / (a + b)
  }
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(x) { x[is.nan(x)] <- 0; x })
  fz <- grepl("_tpi$|_std5", names(d)) & num
  d[fz] <- lapply(d[fz], function(x) { x[is.na(x)] <- 0; x })
  if ("esalc" %in% names(d)) d$esalc <- factor(d$esalc, levels = s2_esa_levels)
  d
}

# Internal: band pairs of the normalised differences (the Landsat classifier uses all bands and the
# focal min/max/sd of B5 and B7, not of B2)
s2_pairs <- function(nms) {
  layers <- c(s2_refl, paste0(rep(c("B08", "B12"), each = 3), c("_min5", "_max5", "_std5")))
  layers <- layers[layers %in% nms]
  if (length(layers) < 2) return(data.frame(a = character(0), b = character(0), name = character(0)))
  cb <- utils::combn(layers, 2)
  data.frame(a = cb[1, ], b = cb[2, ], name = paste(cb[1, ], cb[2, ], sep = "_"), stringsAsFactors = FALSE)
}
s2_pair_names <- function(nms) s2_pairs(nms)$name

#' Classify Sentinel-2 Imagery
#'
#' The Sentinel-2 counterpart of [classify_landsat()]: classifies each image with a random forest model
#' trained on the features of [s2_feature_stack()] and [s2_prep_predictors()]. The output is written to
#' `<base_sentinel2_dir>/<site_name>/classifications/<image file>` with the same layers as the Landsat
#' classifications, at 10 m: `class` (1 land, 2 water, 3 snow, 4 cloud, 5 artif), and the probabilities
#' (%) `land`, `water`, `snow`, `cloud`, `artif` (0 when the model has no such class). The results can be
#' passed to [calc_snow_variables()] with `base_landsat_dir = base_sentinel2_dir`.
#'
#' @param image_df The scene table, output of [extract_sentinel2()].
#' @param site_name Name of the site.
#' @param base_sentinel2_dir Base directory of the Sentinel-2 data.
#' @param model_file Path of the model (`.rds`), e.g. `file.path(model_dir, "latest_S2.rds")`: a model
#' made with [s2_model()] (random forest or LightGBM), or a plain `randomForestSRC` model.
#' @param workers Number of images classified in parallel.
#' @param force Logical, if TRUE also images classified earlier are classified again.
#' @param chunk_cells Number of pixels predicted at a time (limits the memory use; default 250000).
#' @param memmax Maximum memory (GB) terra may use per worker before it works through temporary files
#' (default 2; see [calc_predictors_s2()]). Plan roughly `workers * (memmax + 2)` GB of memory.
#' @return The file names of the classified images.
#' @export
classify_sentinel2 <- function(image_df, site_name, base_sentinel2_dir, model_file, workers = 1,
                               force = FALSE, chunk_cells = 250000, memmax = 2) {
  image_dir <- file.path(base_sentinel2_dir, site_name, "imagery")
  predictor_dir <- file.path(base_sentinel2_dir, site_name, "predictors")
  class_dir <- file.path(base_sentinel2_dir, site_name, "classifications")
  if (!dir.exists(class_dir)) dir.create(class_dir, recursive = TRUE)

  image_df <- image_df %>% dplyr::filter(file %in% list.files(image_dir, pattern = "GMT.tif$"))
  if (!force) image_df <- image_df %>% dplyr::filter(!file %in% list.files(class_dir, pattern = ".tif$"))
  if (nrow(image_df) == 0) {
    message("All imagery classified previously!")
    return(NULL)
  }
  mod <- s2_model_check(readRDS(model_file))

  run_one <- function(f, paral) {
    safe_classify_s2(f, image_df = image_df, image_dir = image_dir, predictor_dir = predictor_dir,
                     class_dir = class_dir, mod = mod, chunk_cells = chunk_cells, paral = paral,
                     memmax = memmax)
  }
  if (workers <= 1 || nrow(image_df) <= 2) {
    lss <- lapply(image_df$file, run_one, paral = FALSE)
  } else {
    set_parallel_plan(workers)
    lss <- future.apply::future_lapply(image_df$file, run_one, paral = TRUE,
                                       future.packages = c("dplyr", "sf", "terra", "lubridate",
                                                           "solartime", "randomForestSRC",
                                                           if (mod$type == "lightgbm") "lightgbm"),
                                       future.seed = TRUE, future.scheduling = 2)
  }
  summarise_classifications(lss)
}

# Internal: an error in one image does not stop the run
safe_classify_s2 <- function(f, ...) {
  tryCatch({
    out <- classify_s2_image(f, ...)
    list(file = f, ok = !is.null(out), msg = NA_character_)
  }, error = function(e) list(file = f, ok = FALSE, msg = conditionMessage(e)))
}

# Internal: classify one image in chunks of rows
classify_s2_image <- function(f, image_df, image_dir, predictor_dir, class_dir, mod, chunk_cells = 250000,
                              paral = TRUE, memmax = 2) {
  old_mem <- s2_terra_memory(memmax)
  on.exit(terra::terraOptions(memmax = old_mem), add = TRUE)
  if (paral) {
    terra::terraOptions(threads = 1)
    options(rf.cores = 1, mc.cores = 1)
  }
  mod <- s2_model_load(mod, threads = if (paral) 1 else NULL)
  vars <- mod$vars
  st <- s2_feature_stack(f, image_dir, predictor_dir, baseline = s2_baseline(image_df, f), vars = vars)
  if (sum(!is.na(terra::values(st[["B04"]], mat = FALSE))) <= 100) return(NULL)
  missing <- setdiff(vars, c(names(st), s2_pair_names(names(st))))
  if (length(missing) > 0) stop("Predictors missing: ", paste(missing, collapse = ", "))

  out_names <- c("class", "land", "water", "snow", "cloud", "artif")
  lab_col <- c("1" = 2, "2" = 3, "3" = 4, "4" = 5, "5" = 6)   # class label -> output column
  lab_col <- lab_col[base::intersect(names(lab_col), mod$classes)]
  vals <- matrix(NA_real_, nrow = terra::ncell(st), ncol = length(out_names))
  rows_per_chunk <- max(1, floor(chunk_cells / terra::ncol(st)))
  for (r0 in seq(1, terra::nrow(st), by = rows_per_chunk)) {
    nr <- min(rows_per_chunk, terra::nrow(st) - r0 + 1)
    d <- as.data.frame(terra::values(st, row = r0, nrows = nr, mat = TRUE))
    d$cell <- terra::cellFromRowCol(st, r0, 1) + seq_len(nrow(d)) - 1
    d <- d[!is.na(d$B04) & !is.na(d$SCL), , drop = FALSE]
    if (nrow(d) == 0) next
    d <- s2_prep_predictors(d, vars)
    d <- d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
    if (nrow(d) == 0) next
    prob <- s2_predict_probs(mod, d)
    if (nrow(prob) != nrow(d)) stop("The number of predictions does not match the number of pixels.")
    vals[d$cell, 1] <- as.numeric(colnames(prob)[max.col(prob, ties.method = "first")])
    vals[d$cell, 2:6] <- 0
    for (lab in names(lab_col)) vals[d$cell, lab_col[[lab]]] <- round(prob[, lab] * 100)
  }
  rr <- terra::rast(st[[1]], nlyrs = length(out_names), names = out_names)
  terra::values(rr) <- vals
  terra::writeRaster(rr, file.path(class_dir, f), datatype = "INT1U", overwrite = TRUE)
  f
}

#' Sentinel-2 Classifier Model Object
#'
#' Wraps a fitted random forest (`randomForestSRC::rfsrc`) or LightGBM booster into the model object
#' that [classify_sentinel2()] reads, so that both kinds of model are predicted with the same code
#' ([s2_predict_probs()]), in model training (cross-validation) and in classification.
#'
#' @param model A fitted `rfsrc` classification forest, or a `lgb.Booster` (multiclass, labels
#' 0, 1, ... in the order of `classes`).
#' @param vars Names of the predictors (for LightGBM: the columns of [s2_model_matrix()], in order).
#' Taken from the forest for `rfsrc`.
#' @param classes Class labels ("1" land, "2" water, "3" snow, "4" cloud, "5" artif) in the order of the
#' model's classes. Taken from the forest for `rfsrc`.
#' @param info Optional list of extra information stored with the model (e.g. CV results).
#' @return A list of class `snowman_s2_model`: `type` ("rfsrc" or "lightgbm"), `vars`, `classes`,
#' `model` (the forest, or the LightGBM model as text) and `info`.
#' @export
s2_model <- function(model, vars = NULL, classes = NULL, info = NULL) {
  if (inherits(model, "rfsrc")) {
    out <- list(type = "rfsrc", vars = model$xvar.names, classes = levels(model$yvar), model = model)
  } else if (inherits(model, "lgb.Booster")) {
    if (is.null(vars) || is.null(classes)) stop("vars and classes are needed for a LightGBM model")
    out <- list(type = "lightgbm", vars = vars, classes = as.character(classes),
                model = model$save_model_to_string())
  } else stop("model must be a randomForestSRC forest or a LightGBM booster")
  out$info <- info
  class(out) <- "snowman_s2_model"
  out
}

# Internal: accept also a plain rfsrc forest
s2_model_check <- function(mod) {
  if (inherits(mod, "snowman_s2_model")) return(mod)
  if (inherits(mod, "rfsrc")) return(s2_model(mod))
  stop("Unknown model object: use s2_model()")
}

# Internal: LightGBM booster from its text (once per process); threads for prediction
s2_model_load <- function(mod, threads = NULL) {
  mod <- s2_model_check(mod)
  if (mod$type == "lightgbm" && is.null(mod$booster)) {
    if (!requireNamespace("lightgbm", quietly = TRUE)) stop("Package lightgbm is needed for this model")
    mod$booster <- lightgbm::lgb.load(model_str = mod$model)
  }
  mod$threads <- threads
  mod
}

#' Predictor Matrix for LightGBM
#'
#' The predictors `vars` of a prepared pixel table ([s2_prep_predictors()]) as a numeric matrix, with
#' `esalc` as integer category codes (0 = class 10, ..., 10 = class 100).
#'
#' @param d Prepared pixel table.
#' @param vars Predictor names.
#' @return A numeric matrix.
#' @export
s2_model_matrix <- function(d, vars) {
  m <- as.data.frame(d)[, vars, drop = FALSE]
  if ("esalc" %in% vars) m$esalc <- as.integer(factor(m$esalc, levels = s2_esa_levels)) - 1L
  as.matrix(m)
}

#' Class Probabilities of a Sentinel-2 Classifier
#'
#' @param model A model of [s2_model()] (or a plain `rfsrc` forest).
#' @param d Prepared pixel table ([s2_prep_predictors()]) with the model's predictors and no missing values.
#' @return A matrix of class probabilities (0-1), one column per class, named by the class labels.
#' @export
s2_predict_probs <- function(model, d) {
  model <- s2_model_load(model, threads = model$threads)
  if (model$type == "rfsrc") {
    p <- predict(model$model, as.data.frame(d)[, model$vars, drop = FALSE], importance = "none",
                 na.action = "na.omit")
    pr <- as.matrix(p$predicted)
    colnames(pr) <- levels(p$class)
  } else {
    prm <- if (is.null(model$threads)) list() else list(num_threads = model$threads)
    pr <- predict(model$booster, s2_model_matrix(d, model$vars), params = prm)
    if (!is.matrix(pr)) pr <- matrix(pr, ncol = length(model$classes), byrow = TRUE)
    colnames(pr) <- model$classes
  }
  pr
}

# Internal: parallel plan as in classify_landsat() (multisession in RStudio and on Windows,
# multicore elsewhere; an existing plan with workers is reused)
set_parallel_plan <- function(workers) {
  if (future::nbrOfWorkers() > 1) return(invisible(NULL))
  is_rstudio_env <- !parallelly::supportsMulticore() && Sys.info()["sysname"] %in% c("Linux", "Darwin")
  if (Sys.info()["sysname"] == "Windows" || is_rstudio_env) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::multicore, workers = workers)
  }
  invisible(NULL)
}
