#' Example GeoMx DSP Dataset
#'
#' A small subset of data from the Nanostring Spatial Organ Atlas Dataset
#'
#' @format ## `example_ds`
#' A list of length 4:
#' \describe{
#'   \item{samples}{tibble of 2 rows and 36 columns, each row representing a region of interest}
#'   \item{meta}{tibble of 2 rows and 3 columns}
#'    \describe{
#'      \item{probe}{list of probes}
#'      \item{gene}{list of genes that correspond to the probes}
#'      \item{bg}{list correlating to background noise}
#'    }
#'   \item{counts}{tibble of 2 rows and 3 columns where each row contains the probes counts for each region of interest}
#'   \item{x}{tibble of 2 rows and 3 columns where each row contains the probe counts for each region of interest}
#'   ...
#' }
#' @source <https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/>
"example_ds"
