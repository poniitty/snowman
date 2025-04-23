#' Download Pretrained Model
#'
#' This function downloads the pretrained model from a GitHub repository.
#'
#' @param model_name Character vector specifying for which satellites you want to download the models. Any of "TM04","TM05","LE07","LC08","LC09" or "ALL"
#' @param model_dir Character, path to folder where the model objects will be stored.
#' @param timeout Maximum time to use to download the files (á ~80MB). Defaults to 300 seconds. Increase if very slow network.
#' @param force Logical. If FALSE (default) downloads the model object only if no previous model file exists.
#' @examples
#' \dontrun{
#' download_model(model_name = "landsat_latest", model_dir = "path/to/your/model")
#' }
#' @export
download_model <- function(model_names = c("TM04","TM05","LE07","LC08","LC09"), model_dir, timeout = 300, force = FALSE) {
  # URL of the GitHub repository
  options(timeout = max(timeout, getOption("timeout")))
  repo_url <- "https://github.com/poniitty/snowman_models/raw/main/"
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_TM05.rds")))){
    if(any(c("tm04","tm05","all") %in% tolower(model_names))){
      print("Downloading TM04 & TM05 model...")
      # Download the model object for Landsat-4 & Landsat 5
      download.file(paste0(repo_url, "Landsat_models/latest_TM05.rds"), 
                    destfile = paste0(model_dir, "/latest_TM05.rds"))
    }
  }
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_LE07.rds")))){
    if(any(c("le07","all") %in% tolower(model_names))){
      print("Downloading LE07 model...")
      # Download the model object for Landsat-7
      download.file(paste0(repo_url, "Landsat_models/latest_LE07.rds"), 
                    destfile = paste0(model_dir, "/latest_LE07.rds"))
    }
  }
  
  if(force == TRUE | (!file.exists(paste0(model_dir, "/latest_LC08.rds")))){
    if(any(c("lc09","lc08","all") %in% tolower(model_names))){
      print("Downloading LC08 & LC09 model...")
      # Download the model object for Landsat-8 & Landsat 9
      download.file(paste0(repo_url, "Landsat_models/latest_LC08.rds"), 
                    destfile = paste0(model_dir, "/latest_LC08.rds"))
    }
  }
}
