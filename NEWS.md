# snowman (development version)

## Landsat classification: features now match the models' training data

* `classify_landsat()`: for Landsat 4-7, the median-image differences (`mr_diff_*`) are now
  calculated as in the models' training data. There, the images had no B1 layer, so the
  median bands B1-B7 were subtracted from the image layers B2-B7 and B10 in layer order (e.g.
  `mr_diff_B7` = median B7 - thermal B10, about -30000; `mr_diff_B1` = median B1 - B2). The
  classifier instead used aligned bands with an empty B1 (`mr_diff_B1` = median B1, all
  values beyond the models' split points). The TM05 model uses `mr_diff_B1`, `_B5`, `_B6`
  and the LE07 model `mr_diff_B1`, `_B5`, `_B6`, `_B7`. Landsat 8-9 were not affected.
  Landsat 4-7 images classified earlier should be classified again (`force = TRUE`).
* `calc_predictors()`: the median indices use the 5 clearest scenes per month (before: 3)
  and the monthly median images the 50 clearest (before: 30), as in the models' training data.

## Sentinel-2 classification

* New `calc_predictors_s2()`, `classify_sentinel2()`: the Sentinel-2 counterparts of
  `calc_predictors()` and `classify_landsat()`. The classifications are at 10 m and
  have the same layers as the Landsat ones (`class`, `land`, `water`, `snow`,
  `cloud`, `artif`), so `calc_snow_variables()` works on them unchanged
  (`base_landsat_dir = base_sentinel2_dir`). One model for Sentinel-2A, 2B and 2C;
  `download_model("S2")` fetches `latest_S2.rds` once it is published.
* `classify_sentinel2()` takes random forest or LightGBM models: `s2_model()` wraps a fitted model,
  `s2_predict_probs()` predicts both kinds with the same code in training and classification
  (`s2_model_matrix()` gives LightGBM's predictor matrix). Only the features the model uses are
  calculated, and pixels are predicted in chunks. `lightgbm` is a suggested package.
* New `s2_feature_stack()` and `s2_prep_predictors()` build the classifier's
  features (shared by model training and classification), and `read_sentinel2()`
  reads an image harmonised across the processing baseline 04.00 change (+1000
  reflectance offset since 25 January 2022).
* `extract_sentinel2()`: new arguments `groups` (download only given scenes, e.g.
  those of training data), `max_per_month` (at most this many scenes per calendar
  month, preferring scenes covering the AOI and low cloud cover) and `min_clear`
  (rejection threshold, previously fixed at 0.1). The processing baseline of each
  scene is recorded in the scene tables (`PROCESSING_BASELINE`; from the product name when the
  catalogue item has no baseline property).
* `extract_sentinel2()` no longer fails on catalogue items without a `created` date (older scenes,
  e.g. 2017): the processing date in the file name then comes from the product generation time or
  the product name.

## Downloads (`extract_landsat()`, `extract_sentinel2()`)

* Images now cover the whole AOI grid: each scene is cut with the snapped AOI grid
  extent (plus 2 pixels) instead of the AOI itself. Before, `crop()` snapped the cut
  to the nearest pixel edge, so the outermost row or column was left empty on each
  side where the AOI edge lay more than half a pixel inside the grid edge; the
  cloud-mask step of `calc_snow_variables()` widened this to 2 pixels without snow
  variables. Images downloaded earlier keep the empty edges.
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
