# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}


geomean <- function(x) {
  return(exp(mean(log(x))))
}


calc_negprobe_stats <- function(meta, counts, negprobe_lod_quantile) {
  proc_auc <- function(responses, predictors) {
    # calculate AUC of signal vs background (negative probe)
    myroc <- pROC::roc(responses, predictors, levels=c(0, 1), direction="<")
    myauc <- as.numeric(myroc$auc)
    return(myauc)
  }
  negprobe_auc <- dplyr::bind_cols(dplyr::select(meta, probe_type), counts) %>%
    dplyr::summarize(dplyr::across(-probe_type, ~ proc_auc(probe_type, .x))) %>%
    unlist()

  negcounts <- counts[meta$probe_type == 0,]
  negprobe_counts <- colSums(negcounts)
  negprobe_geomean <- apply(negcounts, 2, geomean)
  negprobe_lod <- apply(negcounts, 2, quantile, negprobe_lod_quantile)
  return(bind_cols(negprobe_auc = negprobe_auc,
                   negprobe_counts=negprobe_counts,
                   negprobe_geomean=negprobe_geomean,
                   negprobe_lod=negprobe_lod))
}

calc_frac_expr <- function(meta, counts, bg_lod) {
  x <- sweep(counts, 2, bg_lod, FUN="-")
  x <- apply(x, 2, pmax, 0)
  x <- x > 0
  gene_frac_expr <- rowSums(x) / ncol(x)
  y <- x[meta$probe_type == 1,]
  sample_frac_expr <- colSums(y) / nrow(y)
  return(list(s=sample_frac_expr, g=gene_frac_expr))
}


#' @export
st_geomx_read <- function(xlsx_file,
                          sample_sheet = "SegmentProperties",
                          count_sheet = "BioProbeCountMatrix",
                          key = "SegmentDisplayName") {
  # read input data
  s <- read_excel(xlsx_file, sheet=sample_sheet)
  counts <- read_excel(xlsx_file, sheet=count_sheet)
  return(list(s=s, counts=counts))
}


#' @export
st_geomx_process_input <- function(s, meta, counts, negprobe_lod_quantile) {
  # background signal statistics
  x <- calc_negprobe_stats(meta, counts, negprobe_lod_quantile)
  s <- bind_cols(s, x)

  # probes expressed over background
  x <- calc_frac_expr(meta, counts, s$negprobe_lod)
  s$frac_expr <- x$s
  meta$frac_expr <- x$g

  # upper quartile per sample
  s$q3 <- bind_cols(select(meta, probe_type), counts) %>%
    filter(probe_type == 1) %>%
    summarize(across(-probe_type, quantile, 0.75, names = FALSE)) %>%
    unlist()

  # signal to noise ratios
  s$num_counts <- colSums(counts)
  s$snr <- s$num_counts / s$negprobe_counts
  s$q3snr <- s$q3 / s$negprobe_geomean

  list(s=s, meta=meta, counts=counts)
}
