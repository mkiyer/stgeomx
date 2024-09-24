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

  # ensure no duplicate probes
  if (sum(duplicated(meta$probe)) > 0) {
    stop("duplicate probe found")
  }

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
#'
#' @examples
#' sample_tsv <- system.file("extdata", "simple_geomx_sample.tsv", package = "stgeomx")
#' count_tsv <- system.file("extdata", "simple_geomx_count.tsv", package = "stgeomx")
#' read_tsv(sample_tsv, count_tsv)
#'
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

  # Madison edit: remove columns with blank names
  col_names <- colnames(counts)
  valid_cols <- col_names != ""
  counts <- counts[, valid_cols]

  # Madison edit: if count is blank, set to 0
  counts <- counts %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate_all(~ifelse(. == "", 0, .))

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
#'
#' @examples
#' xlsx_file <- system.file("extdata", "simple_geomx.xlsx", package = "stgeomx")
#' read_xlsx(xlsx_file)
#'
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

  # Madison edit: if count is blank, set to 0
  counts <- counts %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate_all(~ifelse(. == "", 0, .))

  # prepare dataset
  ds <- prepare_input(samples, counts, slide_key_col, roi_key_col,
                      segment_key_col, aoi_key_col, probe_key_col,
                      gene_key_col, negprobe_regex)
  return(ds)
}




#'
#' Read Nanostring probe definition file (.pkc)
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param pkc_file probe definition file (.pkc)
#' @returns tibble of probe information
#' @export
read_pkc <- function(pkc_file) {
  pkc <- yaml::read_yaml(pkc_file)
  probes <- tibble("target" = pkc$Targets)
  probes <- probes %>% tidyr::unnest_wider("target") %>%
    dplyr::rename("TargetName" = "DisplayName") %>%
    tidyr::unnest_longer("Probes") %>%
    tidyr::unnest_wider("Probes") %>%
    dplyr::rename("ProbeDisplayName" = "DisplayName") %>%
    dplyr::mutate(bg = startsWith(.data$CodeClass, "Negative"))
  return(probes)
}


#'
#' Read Nanostring digital count conversion (DCC) files
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param dcc_tar_file tarball of DCC files
#' @param dcc_dir temporary directory to store DCC files
#' @returns list of tibbles for 'samples' and 'counts'
#' @export
read_dcc <- function(dcc_tar_file, dcc_dir) {
  # cleanup and overwrite files in 'dcc_dir' if exists
  if (file.exists(dcc_dir)) {
    unlink(dcc_dir, recursive = TRUE)
  }
  # extract dcc tar file
  utils::untar(dcc_tar_file, exdir=dcc_dir)
  # read dcc files
  dcc_files <- list.files(dcc_dir, full.names=TRUE)
  dcc_data <- sapply(dcc_files, GeomxTools::readDccFile, simplify=FALSE)
  dcc_data <- tibble(dcc = dcc_data) %>%
    tidyr::unnest_wider("dcc", simplify=FALSE)

  samples <- dcc_data %>%
    dplyr::select("Header", "Scan_Attributes", "NGS_Processing_Attributes") %>%
    tidyr::unnest(cols = c("Header", "Scan_Attributes", "NGS_Processing_Attributes"))

  counts <- dcc_data %>%
    dplyr::select("Scan_Attributes", "Code_Summary") %>%
    tidyr::hoist("Scan_Attributes", SampleID="SampleID") %>%
    dplyr::select(-"Scan_Attributes") %>%
    tidyr::unnest("Code_Summary") %>%
    tidyr::pivot_wider(names_from="SampleID", values_from="Count", values_fill=1)

  # cleanup
  if (file.exists(dcc_dir)) {
    unlink(dcc_dir, recursive = TRUE)
  }
  return(list(samples=samples,
              counts=counts))
}


#'
#' Read Nanostring dataset from DCC files, probe config file (PKC).
#' and additional sample annotations
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param annot_xlsx_file sample annotations (.xlsx)
#' @param pkc_file probe definition file (.pkc)
#' @param dcc_tar_file tarball of DCC files
#' @param temp_dir temporary directory to store DCC files
#'
#' @param aoi_key_col sample column name corresponding to columns of count sheet
#' @param slide_key_col sample column name containing scan
#' @param roi_key_col sample column name containing roi
#' @param segment_key_col sample column name containing segment
#' @param rts_key_col probe metadata column name for probe RTS ID
#' @param probe_key_col probe metadata column name unique probe identifier
#' @param gene_key_col probe metadata column name unique gene identifier
#' @param negprobe_regex regular expression corresponding to negative probes
#'
#' @returns ds
#' @export
read_pkc_dcc <- function(annot_xlsx_file, pkc_file, dcc_tar_file, temp_dir,
                         slide_key_col = "SlideName",
                         roi_key_col = "ROILabel",
                         segment_key_col = "SegmentLabel",
                         aoi_key_col = "SampleID") {
  # constants (from PKC and DCC files)
  rts_key_col = "RTS_ID"
  probe_key_col = "ProbeDisplayName"
  gene_key_col = "TargetName"

  # read data
  annot <- readxl::read_xlsx(annot_xlsx_file)
  probes <- stgeomx::read_pkc(pkc_file)
  dcc_data <- stgeomx::read_dcc(dcc_tar_file, temp_dir)
  # join sample annotations with sample metadata from dcc files
  samples <- dplyr::inner_join(annot, dcc_data$samples, by=aoi_key_col)
  # join probe metadata with counts
  counts <- dplyr::left_join(probes, dcc_data$counts, by=rts_key_col)

  # ensure all samples found in counts
  inds <- match(pull(samples, {{ aoi_key_col }}), colnames(counts))
  stopifnot(all(!is.na(inds)))

  # break counts into metadata and count data
  meta <- select(counts, "rts"= {{ rts_key_col }},
                 "probe"= {{ probe_key_col }},
                 "gene"= {{ gene_key_col }},
                 "bg")
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
