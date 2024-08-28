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


#'
#' Merge two datasets
#'
#' This function should only be used to merge datasets that are generated
#' from the same probe kit.
#'
#' This function can join datasets with two options: union and intersect.
#' Intersect will keep columns that are common to both datasets,
#' whereas union will keep all columns.
#'
#' @param a dataset produced by st_geomx_read
#' @param b dataset produced by st_geomx_read
#' @returns dataset
#' @export
#'
#' @examples
#' # example code
#' data(example_ds, package= "stgeomx")
#' merge(example_ds, example_ds)
#'
merge <- function(a, b) {
  # merge sample tables
  sample_cols <- intersect(colnames(a$samples), colnames(b$samples))
  s <- bind_rows(select(a$samples, all_of(sample_cols)),
                 select(b$samples, all_of(sample_cols)))
  # merge count data
  a_counts <- bind_cols(select(a$meta, probe, gene, bg), a$counts)
  b_counts <- bind_cols(select(b$meta, probe, gene, bg), b$counts)
  # counts <- inner_join(a_counts, b_counts, by=c("probe", "gene", "bg"),
  #                      relationship="one-to-one")

  # Madison edit: added suffix to avoid addition of ".y" to the end of column names
  counts <- inner_join(a_counts, b_counts, by=c("probe", "gene", "bg"), suffix = c("", ""),
                       relationship="one-to-one")

  # break into metadata and count data
  meta <- select(counts, -all_of(s$aoi))
  counts <- select(counts, all_of(s$aoi))

  # make new dataset object
  ds <- stgeomx::init(samples=s, meta=meta, counts=counts, x=counts)
  return(ds)
}
