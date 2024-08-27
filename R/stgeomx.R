#
# Spatial Transcriptomics using NanoString GeoMx DSP
#

#'
#' initialize dataset
#'
#' @param samples tibble with sample information
#' @param meta tibble with gene metadata
#' @param counts tibble with gene count data
#' @param x tibble with normalized gene count data
#' @returns dataset
#' @export
init <- function(samples, meta, counts, x) {
  ds <- list(
    samples = samples,
    meta = meta,
    counts = counts,
    x = x
  )
  class(ds) <- "stgeomx"
  return(ds)
}

#'
#' copy dataset
#'
#' @param ds stgeomx dataset
#' @returns dataset
#' @export
copy <- function(ds) {
  x <- list(
    samples = ds$samples,
    meta = ds$meta,
    counts = ds$counts,
    x = ds$x
  )
  class(x) <- "stgeomx"
  return(x)
}
