# snowman (development version)

## Downloads (`extract_landsat()`, `extract_sentinel2()`)

* Robust to throttling by the Planetary Computer: GDAL retries failed HTTP
  requests (`http_retries`, `http_retry_delay`), each scene is attempted up to
  `max_retries` times with growing waits (with jitter, so parallel workers do not
  retry in step), and scenes that still fail get up to `retry_rounds` further
  rounds after a pause (`retry_wait`). STAC searches are retried as well, and
  follow the API's paging (no 1000-item limit).
* No silent losses: the outcome of every scene is appended to
  `imagery/download_log.csv`, scenes that could not be downloaded are reported in
  a warning and in `attr(result, "failed")`, and they are attempted again on the
  next run.
* A band that cannot be read now fails the scene (and it is retried) instead of
  being replaced by an empty band. When the same band(s) fail in three attempts in
  a row while all other bands open, the file is taken to be damaged on the server:
  the scene is logged as BROKEN after about 20 s instead of being retried for
  minutes (it is attempted again on the next run).
* Images are written under a temporary name and renamed only after they have
  been read back, so an interrupted run leaves no partial images.
* Incremental updates: a re-run downloads only scenes that are not yet known
  (images on disk, scenes in `lss_final.csv` / `s2s_final.csv` even if their
  images were deleted after classification, and scenes rejected as cloudy,
  recorded in `lss_rejected.csv` / `s2s_rejected.csv`). Only the new images are
  checked, mosaicked and added to the final scene table. Previously, scenes were
  skipped by acquisition date only (a second satellite on the same date was never
  downloaded) and cloudy scenes were downloaded again on every run.
* Sentinel-2: Sentinel-2C is supported (and in the default `sats`); bands are
  opened and resampled one by one (10, 20 and 60 m bands); fixed the mosaicking
  of more than two same-day images and the duplicate handling of more than two
  overlapping items.
