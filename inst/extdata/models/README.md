# Classifier models included in snowman

Put the latest models here with these exact names:

| File               | Sensor                        | Read by                                   |
|--------------------|-------------------------------|-------------------------------------------|
| `latest_TM05.rds`  | Landsat 4-5 TM                | `classify_landsat()`                      |
| `latest_LE07.rds`  | Landsat 7 ETM+                | `classify_landsat()`                      |
| `latest_LC08.rds`  | Landsat 8-9 OLI               | `classify_landsat()`                      |
| `latest_S2.rds`    | Sentinel-2 (S2A, S2B, S2C)    | `classify_sentinel2()`, `extract_sentinel2()`, `calc_predictors_s2()` |

Each file is a `snowman::s2_model()` object (LightGBM or random forest), saved with
`saveRDS(model, file, compress = TRUE)` by the model-fitting scripts.

When the package is installed, this folder becomes `<library>/snowman/extdata/models/`, and the
functions find the models with `snowman_model_file("LC08")` etc. A `model_dir` (classify_landsat) or
`model_file` (Sentinel-2 functions) given by the user takes precedence.

`extract_sentinel2(bands = "model")` and `calc_predictors_s2()` read `latest_S2.rds` to find out which
bands and median images the model needs, so replacing the Sentinel-2 model can change what is
downloaded and calculated.
