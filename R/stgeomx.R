#
# Spatial Transcriptomics using NanoString GeoMx DSP
#

#'
#' initialize dataset
#'
#' @param samples tibble with sample information
#' @param meta tibble with gene metadata
#' @param counts tibble with gene count data
#' @param x tibble with gene count data
#' @returns dataset
#' @export
#'
#' @examples
#' samples <- data.frame(
#'   SlideName = c("brain_001", "brain_002"),
#'   ScanLabel = c("Brain 1", "Brain 2"),
#'   aoi = c("s01_r001_fullroi", "s02_r001_fullroi")
#' )
#' meta <- data.frame(
#'   probe = c("example01", "NegProbe-WTX01"),
#'   gene = c("example", "NegProbe-WTX"),
#'   bg = c(FALSE, TRUE)
#' )
#' counts <- data.frame(
#'   s01_r001_fullroi <- c(100, 32),
#'   s02_r001_fullroi <- c(95, 24)
#' )
#' x <- data.frame(
#'   s01_r001_fullroi <- c(100, 32),
#'   s02_r001_fullroi <- c(95, 24)
#' )
#' new_ds <- init(samples, meta, counts, x)
#'
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
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' duplicate <- copy(example_ds)
#'
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
