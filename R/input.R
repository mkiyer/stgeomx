#
# Spatial Transcriptomics using NanoString GeoMx DSP
#

#'
#' Standardize names using Slide, ROI, Segment, and AOI columns of sample
#' annotation table. Partition count table into metadata and counts.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import rlang
#'
#' @param samples tibble of sample information
#' @param counts tibble of count data
#' @param slide_key_col sample column name containing scan
#' @param roi_key_col sample column name containing roi
#' @param segment_key_col sample column name containing segment
#' @param aoi_key_col sample column name corresponding to columns of count sheet
#' @param probe_key_col count sheet column name unique probe identifier
#' @param gene_key_col count sheet column name unique gene identifier
#' @param negprobe_regex regular expression corresponding to negative probes
#' @returns ds
prepare_input <- function(samples, counts,
                          slide_key_col,
                          roi_key_col,
                          segment_key_col,
                          aoi_key_col,
                          probe_key_col,
                          gene_key_col,
                          negprobe_regex) {
  # ensure all samples found in counts
  inds <- match(pull(samples, {{ aoi_key_col }}), colnames(counts))
  stopifnot(all(!is.na(inds)))

  # break counts into metadata and count data
  meta <- counts %>%
    select("probe" = {{ probe_key_col }}, "gene" = {{ gene_key_col }}) %>%
    mutate(bg = stringr::str_detect(.data$gene, negprobe_regex), .after="gene")

  # select valid samples
  counts <- select(counts, all_of(inds))
  # ensure samples and count tables are aligned
  stopifnot(all(pull(samples, {{ aoi_key_col }}) == colnames(counts)))

  # create new identifiers for slide, roi, segment, aoi
  slide_keys <- pull(samples, {{ slide_key_col }})
  samples <- samples %>% mutate(
    slide = sprintf("s%02d", match(slide_keys, unique(slide_keys))),
    roi = sprintf("r%s", .data[[roi_key_col]]),
    segment = tolower(gsub(" ", "", .data[[segment_key_col]], fixed = TRUE)),
    aoi = sprintf("%s_%s_%s", .data$slide, .data$roi, .data$segment)
  )
  stopifnot(n_distinct(samples$aoi) == nrow(samples))

  # relabel columns of count matrix with new aoi id
  colnames(counts) <- samples$aoi
  stopifnot(all(samples$aoi == colnames(counts)))

  # create dataset
  ds <- stgeomx::init(samples, meta, counts, counts)
  return(ds)
}


#'
#' Read Nanostring GeoMx DSP data from tab-delimited text files
#'
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @import dplyr
#' @import rlang
#'
#' @param sample_tsv_file tab-delimited text file with sample annotations
#' @param count_tsv_file tab-delimited text file with count matrix/metadata
#' @param slide_key_col sample column name containing scan
#' @param roi_key_col sample column name containing roi
#' @param segment_key_col sample column name containing segment
#' @param aoi_key_col sample column name corresponding to columns of count sheet
#' @param probe_key_col count sheet column name unique probe identifier
#' @param gene_key_col count sheet column name unique gene identifier
#' @param negprobe_regex regular expression corresponding to negative probes
#' @returns dataset
#' @export
read_tsv <- function(sample_tsv_file,
                     count_tsv_file,
                     slide_key_col = "ScanLabel",
                     roi_key_col = "ROILabel",
                     segment_key_col = "SegmentLabel",
                     aoi_key_col = "SegmentDisplayName",
                     probe_key_col = "ProbeDisplayName",
                     gene_key_col = "TargetName",
                     negprobe_regex = "^NegProbe") {
  # read samples (remove duplicates based on aoi key)
  samples <- readr::read_tsv(sample_tsv_file,
                             na = c("", "NA", "NaN"),
                             show_col_types=FALSE) %>%
    distinct(.data[[aoi_key_col]], .keep_all=TRUE)

  # read count data
  counts <- readr::read_tsv(count_tsv_file,
                            show_col_types=FALSE,
                            name_repair="minimal")
  # prepare dataset
  ds <- prepare_input(samples, counts, slide_key_col, roi_key_col,
                      segment_key_col, aoi_key_col, probe_key_col,
                      gene_key_col, negprobe_regex)
  return(ds)
}


#'
#' Read Nanostring GeoMx DSP data from an Excel (XLSX) file
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param xlsx_file Excel file containing input data
#' @param sample_sheet Excel sheet name containing sample information
#' @param count_sheet Excel sheet name containing counts
#' @param slide_key_col sample sheet column name containing scan
#' @param roi_key_col sample sheetcolumn name containing roi
#' @param segment_key_col sample sheet column name containing segment
#' @param aoi_key_col sample sheet column name corresponding to columns of count sheet
#' @param probe_key_col count sheet column name unique probe identifier
#' @param gene_key_col count sheet column name unique gene identifier
#' @param negprobe_regex regex corresponding to negative probes
#' @returns list containing samples, metadata, and counts
#' @export
read_xlsx <- function(xlsx_file,
                      sample_sheet = "SegmentProperties",
                      count_sheet = "BioProbeCountMatrix",
                      slide_key_col = "ScanLabel",
                      roi_key_col = "ROILabel",
                      segment_key_col = "SegmentLabel",
                      aoi_key_col = "SegmentDisplayName",
                      probe_key_col = "ProbeDisplayName",
                      gene_key_col = "TargetName",
                      negprobe_regex = "^NegProbe") {
  # read samples (remove duplicates based on aoi key)
  samples <- readxl::read_excel(xlsx_file, sheet=sample_sheet,
                                na = c("", "NaN")) %>%
    distinct(.data[[aoi_key_col]], .keep_all=TRUE)

  # read counts
  counts <- readxl::read_excel(xlsx_file, sheet=count_sheet,
                               .name_repair = "minimal")

  # prepare dataset
  ds <- prepare_input(samples, counts, slide_key_col, roi_key_col,
                      segment_key_col, aoi_key_col, probe_key_col,
                      gene_key_col, negprobe_regex)
  return(ds)
}
