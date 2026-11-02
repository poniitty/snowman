# ==========================================================================
# Internal helpers shared by extract_landsat() and extract_sentinel2():
# STAC searches with retries, bookkeeping of known scenes, download rounds with
# retries, and the download log.
#
# Scenes are identified by a "scene group": satellite + orbit unit (Landsat WRS
# path / Sentinel-2 MGRS tile) + acquisition date. Same-date images of one group
# are mosaicked into one file, so a group corresponds to one output image.
# ==========================================================================

# Internal: wait before a retry, growing exponentially, with jitter so that many
# parallel workers do not retry at the same moment.
backoff_wait <- function(attempt, base = 5, cap = 300) {
  Sys.sleep(min(cap, base * 2^(attempt - 1)) * stats::runif(1, 0.75, 1.25))
}

# Internal: run a STAC request (a function without arguments) with retries.
# Rate-limited or failed requests are retried with growing waits; the last error
# is raised, so a failed search stops the function instead of returning nothing.
stac_with_retry <- function(request, max_retries = 6, label = "STAC request") {
  for (i in seq_len(max_retries)) {
    res <- try(request(), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    if (i == max_retries) stop(label, " failed after ", max_retries, " attempts: ", as.character(res))
    message(label, " failed (attempt ", i, "), retrying: ", trimws(as.character(res)))
    backoff_wait(i, base = 10)
  }
}

# Internal: STAC search over the whole period, following the API's paging
# (no limit of 1000 items).
stac_search_all <- function(collection, aoi, start_date, end_date, max_retries = 6) {
  s_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  bb <- sf::st_bbox(sf::st_transform(aoi, 4326))
  stac_with_retry(function() {
    s_obj %>%
      rstac::stac_search(collections = collection, bbox = as.numeric(bb),
                         datetime = stac_datetime(start_date, end_date), limit = 1000) %>%
      rstac::get_request() %>%
      rstac::items_fetch(progress = FALSE)
  }, max_retries = max_retries, label = paste("STAC search", collection))
}

# Internal: scene group keys from output file names
# (<sat>_<x>_<path/row or tile>_<date>_...): satellite, first three characters of
# field 3 for Landsat (WRS path) or the whole field for Sentinel-2 (tile), date.
group_from_file <- function(files, sensor = c("landsat", "sentinel2")) {
  sensor <- match.arg(sensor)
  if (!length(files)) return(character(0))
  p <- strsplit(basename(files), "_")
  vapply(p, function(z) {
    if (length(z) < 4) return(NA_character_)
    unit <- if (sensor == "landsat") substr(z[3], 1, 3) else z[3]
    paste(z[1], unit, z[4], sep = "_")
  }, character(1))
}

# Internal: scene groups that need no download: images on disk, scenes already in
# the final scene table (their images may have been deleted after
# classification), and scenes rejected earlier as too cloudy. Failed downloads
# are not in any of these, so they are attempted again on the next run.
known_groups <- function(image_dir, final_csv, rejected_csv, sensor) {
  f <- list.files(image_dir, pattern = "GMT.tif$")
  for (csv in c(final_csv, rejected_csv)) {
    p <- file.path(image_dir, csv)
    if (file.exists(p)) f <- c(f, readr::read_csv(p, show_col_types = FALSE, col_types = readr::cols(.default = "c"))$file)
  }
  unique(stats::na.omit(group_from_file(f, sensor)))
}

# Internal: download features in parallel, then retry the failed ones in further
# rounds after a pause (throttling is usually temporary). FUN returns a one-row
# tibble (file, status, attempts, message) per feature. Scenes with a damaged band
# on the server (status BROKEN) are not retried in further rounds.
download_rounds <- function(features, FUN, fun_args, retry_rounds = 2, retry_wait = 60,
                            packages = c("sf", "terra", "stringr", "rstac")) {
  empty <- dplyr::tibble(file = character(0), status = character(0), attempts = integer(0),
                         message = character(0), round = integer(0))
  if (!length(features)) return(empty)
  out <- list()
  todo <- seq_along(features)
  for (round in 0:retry_rounds) {
    if (!length(todo)) break
    if (round > 0) {
      message(length(todo), " scene(s) failed; waiting ", retry_wait * round,
              " s before retry round ", round, " of ", retry_rounds)
      Sys.sleep(retry_wait * round)
    }
    res <- do.call(future.apply::future_lapply,
                   c(list(X = features[todo], FUN = FUN), fun_args,
                     list(future.packages = packages, future.seed = TRUE,
                          future.scheduling = if (length(todo) > 20) 10 else 1)))
    res <- dplyr::bind_rows(res)
    res$round <- round
    out[[length(out) + 1]] <- res
    todo <- todo[res$status == "ERROR"]
  }
  res <- dplyr::bind_rows(out)
  # one row per scene: its last attempt
  res %>% dplyr::group_by(file) %>% dplyr::slice_tail(n = 1) %>% dplyr::ungroup()
}

# Internal: append the outcome of this run's downloads to <image_dir>/download_log.csv
write_download_log <- function(res, image_dir) {
  if (is.null(res) || !nrow(res)) return(invisible(NULL))
  res$time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  f <- file.path(image_dir, "download_log.csv")
  readr::write_csv(res[, c("time", "file", "status", "attempts", "round", "message")], f, append = file.exists(f))
  invisible(NULL)
}

# Internal: summary message (and warning for failures) after the downloads
report_downloads <- function(res, n_known, image_dir) {
  if (is.null(res)) res <- dplyr::tibble(status = character(0))
  n_ok <- sum(res$status %in% c("SUCCESS", "EXISTS")); n_err <- sum(res$status == "ERROR")
  n_broken <- sum(res$status == "BROKEN")
  message("Scenes: ", n_known, " known from earlier runs, ", n_ok, " downloaded, ", n_err, " failed",
          if (n_broken) paste0(", ", n_broken, " with a damaged band on the server") else "")
  if (n_err + n_broken > 0)
    warning(n_err + n_broken, " scene(s) could not be downloaded (see ", file.path(image_dir, "download_log.csv"),
            "). They are attempted again when the function is run again.", call. = FALSE)
}
