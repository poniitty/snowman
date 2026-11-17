# snowman (development version)

## Parallel classification on systems with OpenMP thread binding

* `classify_landsat()` and `classify_sentinel2()`: the parallel workers now get one OpenMP thread
  each and no OpenMP thread binding. Where binding is switched on (`OMP_PROC_BIND` / `OMP_PLACES`, as
  on CSC Roihu), every worker bound its main thread to the same first core of the job, so all
  workers shared one core (test with 2 workers: 15.8 s before, 7.0 s now, 10.1 s sequential).
  Workers that were started earlier (a reused plan, forked processes) release the binding themselves.

## Snow variables from Sentinel-2; one sensor per base directory

* `calc_snow_variables()` works on Sentinel-2 classifications (10 m) as on Landsat ones. New argument
  `sensor` ("landsat" / "sentinel2"): the snow variables come from one sensor; NULL (default) takes it
  from the file names and stops if `image_df` has both. Steps tied to the pixel size use ground
  distances: clear pixels within 30 m of a cloud are left out (3 x 3 pixels for Landsat as before,
  7 x 7 for Sentinel-2), and the land cover is aggregated only if it is finer than the output.
  For Landsat, the `esalc` layer is now the 30 m class of ESALC.tif; earlier versions aggregated it
  a second time (most common class of 90 x 90 m). The snow variables themselves are unchanged.
* `search_image_df()` recognises Sentinel-2 images and returns `s2s_final.csv`, completed with images
  missing from it.
* Landsat and Sentinel-2 data of a site must be in separate base directories (they use the same
  predictor and table file names on different grids). `extract_landsat()`, `extract_sentinel2()`,
  `calc_predictors()` and `calc_predictors_s2()` stop with a message if the site folder holds data of
  the other sensor.

## Models included in the package; only what the model needs is downloaded and calculated

* The latest classifier models can be shipped with the package in `inst/extdata/models/`
  (`latest_TM05.rds`, `latest_LE07.rds`, `latest_LC08.rds`, `latest_S2.rds`). New
  `snowman_model_file()` finds them. `classify_landsat(model_dir = NULL)` and
  `classify_sentinel2(model_file = NULL)` use them by default; a `model_dir` / `model_file` given by
  the user takes precedence. lightgbm moved from Suggests to Imports (the models are LightGBM models).
* `extract_sentinel2()` has a new argument `bands`. The default `"model"` downloads only the bands
  the Sentinel-2 model needs (its predictors' bands plus B02, B03, B04, B08, B11 and SCL, which
  snowman itself needs); with the reduced model B02, B03, B04, B06, B08, B11, B12 and SCL: about a
  third less to download and 40% less disk space. `"all"` downloads all 13 layers (e.g. for training
  data); a vector of bands is completed with the bands snowman and the model need, whatever the
  choice. Images of one acquisition with different band sets (earlier and new downloads) are
  mosaicked from their common bands. A warning tells if images downloaded earlier lack bands needed
  now. New `s2_model_requirements()` lists what a model needs.
* `calc_predictors_s2()` has new arguments `median_images` and `model_file`. The monthly median
  images, its slowest step, are calculated only if the model uses them (default `NULL`), or always
  with `TRUE`; `FALSE` is overridden (with a message) when the model uses them. The median images
  have the reflectance bands present in the images.
* `read_sentinel2()`, `s2_feature_stack()` and `classify_sentinel2()` work with images of any band
  subset; `s2_feature_stack()` reads the median images only if the model uses them.
  `classify_sentinel2()` skips (with a warning) images that lack bands the model uses.
* `classify_landsat()` calculates only the predictors the model uses: the median-image differences
  and median indices only when used (before always), only the band pairs in the model (before all
  ~190 pairs of every pixel), and the pixel table only of the needed layers. The QA bits are read
  with vectorised bit operations instead of a per-pixel R function (about 1 s per QA layer and
  scene), and the cloud bit once for the 5- and 25-pixel proportions. Identical results, 30-45%
  faster per scene in tests.

## Landsat classification: features now match the models' training data

* `classify_landsat()` also takes LightGBM models (made with `s2_model()`, e.g. by
  fit_GBM_model_landsat.R); random forest models (plain `rfsrc` forests) work as before. Both are
  predicted with `s2_predict_probs()`; the output layers are unchanged.

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
* `extract_sentinel2()`: the tiles of one acquisition (satellite, date and time) are mosaicked into
  one image when the AOI straddles MGRS tiles, so an acquisition is one observation (before, each
  tile was a separate image and the ~10 km tile overlap was observed twice). Copies of one tile are
  mosaicked the same way. Each pixel comes from one input image with all its bands and SCL (the
  image with most valid pixels first); values are no longer averaged (averaging also mixed SCL class
  codes). Scenes with and without the baseline 04.00 offset are not mixed. Merged tile images are
  recorded in `s2s_merged.csv`. Images of earlier runs that were not mosaicked are mosaicked on the
  next run.
* `extract_sentinel2()` is several times faster: each scene is read with GDAL's warp from a virtual
  raster of its 13 band files, straight onto the AOI's 10 m grid, instead of reading each band with
  terra and reprojecting it with `terra::project()`. On Roihu (4 cores), reprojection with terra took
  8 s per scene alone and 59 s with 8 parallel workers; GDAL's warp took 12 s and 21 s for the whole
  scene (open, read, warp, write). Pixel values are the same (nearest neighbour; >99.8% identical,
  the rest differ only where a pixel centre falls on a boundary between source pixels).
* `calc_predictors_s2()` is faster: each scene is masked and harmonised in one pass (median images
  and indices), the gaps of the monthly median images are interpolated over the months only where
  there are gaps (vectorised; `terra::approximate()` took ~2 min per band at 10 m, now seconds), and
  the per-scene and per-band steps run in parallel with the new `workers` argument. The results are
  identical; the time axis of the interpolation is now the month number (before: the order of the
  available months; the same when no month is missing).
* `extract_sentinel2()` downloads only the tiles of an acquisition that are needed to cover the AOI
  (by the catalogue footprints; MGRS tiles overlap widely near UTM zone edges, so one tile often
  covers a site that touches two). After the download the decision is checked against the data:
  if the downloaded tiles leave more than 0.5% of the AOI without valid pixels, the skipped tiles of
  that acquisition are downloaded too. Of several processings of the same tile and acquisition (the
  catalogue keeps reprocessed versions, e.g. 2022-2023 scenes reprocessed in 2024), only the newest
  is downloaded; older ones only if a gap remains after the other tiles.
* `extract_sentinel2()`: two passes of one satellite on the same day (high latitudes) are separate
  scenes (the scene key includes the time).
* `extract_sentinel2()`: footprints are intersected with the AOI itself, not its bounding box in
  degrees; tiles that do not reach the AOI (status NO_OVERLAP) are recorded and not retried (before,
  they failed after five attempts in three rounds, delaying every run).
* `extract_sentinel2(max_per_month = )` counts acquisitions (with all their tiles) and takes them in
  turn from the different years, so that the scenes of the median images are spread over the years.
* `calc_predictors_s2()`, `classify_sentinel2()`: new `memmax` limits terra's memory use. terra
  estimates the free memory from the whole node, which on a shared HPC node is far more than a job
  may use; R could crash when computing the median images. Intermediates are kept on disk. The
  download workers of `extract_sentinel2()` are limited to 1 GB each in the same way.
* `extract_sentinel2()` works with catalogue items that lack optional properties (`created`,
  `view:sun_elevation`, `s2:granule_id`, `eo:cloud_cover`; common for older scenes): items are
  identified by their item id, missing values are left NA, and file names are built by one function
  for both the download and the scene table.
* `extract_sentinel2()`: fixed "subscript out of bounds" when taking the pixel grid from a scene's
  B04 band (`rstac::assets_url()` returns unnamed URLs).

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
