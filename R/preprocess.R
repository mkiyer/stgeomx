# Spatial Transcriptomics using NanoString GeoMx DSP
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


#'
#' Geometric mean of a vector
#'
#' @param x vector of numeric values
#' @returns numeric
geomean <- function(x) {
  return(exp(mean(log(x))))
}


#'
#' Compute quantiles of values in vector 'x' along
#' compute ecdf of background distribution and use it to
#' interpolate quantiles of genes in vector
#'
#' @param x vector of numeric values
#' @param bg vector of TRUE/FALSE corresponding to x
#' @param bw.adjust bandwidth adjustment multiplier for 'density' function
#' @return vector of quantiles corresponding to array 'x'
compute_quantiles_kde <- function(x, bg, bw.adjust=2) {
  # setup
  n <- 2^13
  xlog <- log2(x)
  # select bandwidth
  bw <- stats::bw.nrd0(xlog[bg]) * bw.adjust
  # smooth ecdf function for neg probes
  dbg <- stats::density(xlog[bg], bw=bw, n=n)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)
  # find bg quantiles
  q <- stats::approx(dbg$x, dbgcdf, xout=xlog, yleft=0, yright=1, ties="ordered")$y
  return(q)
}


#'
#' Compute probe quantiles relative to background probe ecdf
#'
#' @param ds dataset
#' @returns matrix of quantiles corresponding to genes (rows) and aois (cols)
#
calc_bg_quantiles <- function(ds) {
  # calculate probe quantiles relative to negative probes
  bg_quants <- bind_cols(select(ds$meta, bg), ds$counts) %>%
    reframe(across(-bg, ~ compute_quantiles_kde(.x, bg)))
  return(bg_quants)
}


#'
#' Remove AOIs with no valid counts
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list produced by st_geomx_merge_probes
#' @param min_count integer count threshold for removing aoi
#' @returns list dataset
remove_empty_aois <- function(ds, min_count = 1) {
  # initial filter of "empty" samples with essential zero counts
  count_max = apply(ds$counts, 2, max)
  keep <- count_max > min_count
  fds <- list(
    samples = filter(ds$samples, keep),
    meta = ds$meta,
    counts = ds$counts[, keep]
  )
  return(fds)
}


#'
#' Calculate QC metrics
#' adds qc columns to the samples tibble
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds dataset
#' @param bg_lod_quantile number (0.0-1.0)
#' @returns samples tibble with added qc columns
calc_qc_metrics <- function(samples, meta, counts, bg_lod_quantile) {
  # unpack dataset for readability
  bg <- meta$bg
  x <- counts[!bg, ]
  xbg <- counts[bg, ]

  # add background metrics
  samples <- samples %>% mutate(
    bg_counts = colSums(xbg),
    bg_cpm = 1e6 * colSums(xbg) / colSums(counts),
    bg_geomean = apply(xbg, 2, geomean),
    bg_median = apply(xbg, 2, stats::median),
    bg_q3 = apply(xbg, 2, stats::quantile, 0.75),
    bg_q95 = apply(xbg, 2, stats::quantile, 0.95),
    bg_max = apply(xbg, 2, max),
    bg_lod = apply(xbg, 2, stats::quantile, bg_lod_quantile)
  )
  # add gene metrics
  samples <- samples %>% mutate(
    num_counts = colSums(x),
    complexity = apply(x, 2, function(x) length(unique(x))),
    count_mad = apply(x, 2, stats::mad),
    count_iqr = apply(x, 2, stats::IQR),
    count_geomean = apply(x, 2, geomean),
    count_median = apply(x, 2, stats::median),
    count_q3 = apply(x, 2, stats::quantile, 0.75),
    count_q95 = apply(x, 2, stats::quantile, 0.95),
    count_max = apply(x, 2, max)
  )
  # gene probe vs background probe metrics (signal-to-noise)
  samples <- samples %>% mutate(
    snr_geomean = .data$count_geomean / .data$bg_geomean,
    snr_median = .data$count_median / .data$bg_median,
    snr_q3 = .data$count_q3 / .data$bg_q3,
    snr_q95 = .data$count_q95 / .data$bg_q95
  )
  return(samples)
}


#'
#' Compute frac expressed genes above threshold
#'
#' @param counts tibble counts
#' @param bg boolean vector of TRUE/FALSE corresponding to count rows
#' @param bg_lod vector of limit of detection thresholds per sample
#' @returns list
calc_frac_expr <- function(counts, bg, bg_lod) {
  x <- sweep(counts, 2, bg_lod, FUN="-")
  x <- apply(x, 2, pmax, 0)
  x <- x > 0
  gene_frac_expr <- rowSums(x) / ncol(x)
  y <- x[bg == FALSE,]
  sample_frac_expr <- colSums(y) / nrow(y)
  return(list(s=sample_frac_expr, g=gene_frac_expr))
}


#'
#' Calculate gene vs background AUC (snAUC) for each AOI
#'
#' @param meta tibble count metadata
#' @param counts tibble counts
#' @returns vector of AUC values for each AOI
calc_snauc <- function(meta, counts) {
  # calculate AUC of gene vs background (negative probe)
  calc_roi_auc <- function(responses, predictors) {
    myroc <- pROC::roc(responses, predictors, levels=c(TRUE, FALSE), direction="<")
    myauc <- as.numeric(myroc$auc)
    return(myauc)
  }
  bg_auc <- bind_cols(select(meta, bg), counts) %>%
    summarize(across(-bg, ~ calc_roi_auc(bg, .x))) %>%
    unlist()
  return(bg_auc)
}


#'
#' Preprocess input data
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds dataset produced by st_geomx_read
#' @param bg_lod_quantile number (0.0-1.0)
#' @returns dataset
#' @export
preprocess <- function(ds, bg_lod_quantile=0.9) {
  # first remove empty/control AOIs
  ds <- remove_empty_aois(ds)

  # unpack dataset
  samples <- ds$samples
  meta <- ds$meta
  bg <- meta$bg
  counts <- ds$counts

  # if dataset is 1-based, convert to 0-based for consistency
  if (min(ds$counts) > 0) {
    counts <- counts - 1
  }

  # compute basic qc metrics (returns sample table)
  samples <- calc_qc_metrics(samples, meta, counts, bg_lod_quantile)

  # separate bg and gene counts
  x <- counts[!bg, ]
  xbg <- counts[bg, ]

  # probes expressed over detection limit
  x <- calc_frac_expr(counts, bg, samples$bg_lod)
  samples$frac_expr <- x$s
  meta$frac_expr <- x$g

  # calculate AUC of gene vs background (negative probe)
  samples$bg_auc <- calc_snauc(meta, counts)

  # KS test of gene vs background
  bg_ks_dist <- NULL
  bg_ks_pval <- NULL
  for (i in 1:ncol(counts)) {
    x <- unlist(counts[!bg, i])
    xbg <- unlist(counts[bg, i])
    res <- suppressWarnings(stats::ks.test(x, xbg, alternative="less"))
    bg_ks_dist[i] <- res$statistic
    bg_ks_pval[i] <- res$p.value
  }
  samples <- mutate(samples,
    bg_ks_dist = bg_ks_dist,
    bg_ks_pval = bg_ks_pval
  )

  # filtering defaults
  samples$keep <- TRUE
  meta$keep <- TRUE

  return(list(samples=samples,
              meta=meta,
              counts=counts))
}

#'
#' Merge multiple probes per gene
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import rlang
#'
#' @param ds list produced by st_geomx_process
#' @returns list
#' @export
merge_probes <- function(ds) {
  # convenience function for geomean
  geomean_fast <- function(x) {
    ifelse(length(x) == 1, x, exp(mean(log(x))))
  }

  # do not merge the negative probes, instead just move the probe name
  # to the gene name
  ybg <- bind_cols(ds$meta, ds$counts) %>%
    filter(.data$bg) %>%
    rowwise() %>%
    summarize(
      gene = .data$probe,
      num_probes = 1,
      probes = .data$probe,
      bg = TRUE,
      frac_expr = .data$frac_expr,
      across(colnames(ds$counts), geomean_fast)
    )
  bgmeta <- ybg %>% select(-colnames(ds$counts))
  bgcounts <- ybg %>% select(colnames(ds$counts))

  # merge gene probes
  y <- bind_cols(ds$meta, ds$counts) %>%
    filter(!.data$bg) %>%
    group_by(.data$gene) %>%
    summarize(
      num_probes = n(),
      probes = stringr::str_flatten(.data$probe, collapse=","),
      bg = FALSE,
      frac_expr = mean(.data$frac_expr),
      across(colnames(ds$counts), geomean_fast)
    )
  meta <- y %>% select(-colnames(ds$counts))
  counts <- y %>% select(colnames(ds$counts))

  ds <- list(samples=ds$samples,
             meta=bind_rows(meta, bgmeta),
             counts=bind_rows(counts, bgcounts))
  return(ds)
}



#'
#' Apply user-defined thresholds to AOIs in preparation of filtering
#'
#' This function does not perform filtering. It sets the 'keep' column in the
#' sample and gene metadata allowing users to generate plots of different
#' thresholds.
#'
#' Sequencing depth and signal-to-noise ratio can vary dramatically across
#' experiments. The default parameter values will not perform any filtering.
#' We recommend using the plotting functions to aid in choosing
#' appropriate parameter values.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list produced by st_geomx_merge_probes
#' @param min_counts integer minimum counts per aoi
#' @param min_auc number minimum snAUC per aoi
#' @param aoi_min_frac_expr number minimum frac of genes expressed per aoi
#' @param gene_min_frac_expr number minimum frac of samples expressing gene
#' @returns list dataset
#' @export
st_geomx_set_thresholds <- function(ds,
                                    min_counts = 0,
                                    min_auc = 0.5,
                                    aoi_min_frac_expr = 0.0,
                                    gene_min_frac_expr = 0.0) {
  # apply filtering parameters to samples
  ds$samples <- ds$samples %>% mutate(
    keep = ((.data$num_counts > min_counts) &
              (.data$bg_auc > min_auc) &
              (.data$frac_expr > aoi_min_frac_expr))
  )
  # apply filtering parameters to genes
  ds$meta <- ds$meta %>% mutate(
    keep = .data$frac_expr > gene_min_frac_expr,
    expressed = case_when(.data$bg ~ "bg",
                          !(.data$keep) ~ "no",
                          TRUE ~ "yes")
  )
  return(ds)
}


#'
#' Filter AOIs based on thresholds set by st_geomx_filter_thresholds
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import rlang
#'
#' @param ds list dataset produced by st_geomx_filter_thresholds
#' @returns list dataset
#' @export
st_geomx_filter <- function(ds) {
  keep_genes <- ds$meta$keep & (!ds$meta$bg)
  fds <- list(
    samples = filter(ds$samples, .data$keep),
    meta = filter(ds$meta, keep_genes),
    counts = ds$counts[keep_genes, ds$samples$keep],
    bgmeta = filter(ds$meta, ds$meta$bg),
    bgcounts = ds$counts[ds$meta$bg, ds$samples$keep]
  )
  return(fds)
}

