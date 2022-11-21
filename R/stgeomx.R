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
#' Read Nanostring GeoMx DSP data
#'
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
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
#' @param negprobe_gene_name name of gene corresponding to negative probes
#' @returns list containing samples, metadata, and counts
#' @export
st_geomx_read <- function(xlsx_file,
                          sample_sheet = "SegmentProperties",
                          count_sheet = "BioProbeCountMatrix",
                          slide_key_col = "ScanLabel",
                          roi_key_col = "ROILabel",
                          segment_key_col = "SegmentLabel",
                          aoi_key_col = "SegmentDisplayName",
                          probe_key_col = "ProbeDisplayName",
                          gene_key_col = "TargetName",
                          negprobe_gene_name = "NegProbe-WTX") {
  # read samples (remove duplicates based on aoi key)
  samples <- read_excel(xlsx_file, sheet=sample_sheet) %>%
    distinct(.data[[aoi_key_col]], .keep_all=TRUE)

  # read counts
  counts <- read_excel(xlsx_file, sheet=count_sheet,
                       .name_repair = "minimal")
  # ensure all samples found in counts
  inds <- match(pull(samples, !!aoi_key_col), colnames(counts))
  stopifnot(all(!is.na(inds)))

  # break counts into metadata and count data
  meta <- counts %>%
    select(probe = !!probe_key_col, gene = !!gene_key_col) %>%
    mutate(bg = (gene == negprobe_gene_name), .after="gene")

  # select valid samples
  counts <- select(counts, all_of(inds))
  # ensure samples and count tables are aligned
  stopifnot(all(pull(samples, !!aoi_key_col) == colnames(counts)))

  # create new identifiers for slide, roi, segment, aoi
  slide_keys <- pull(samples, !!slide_key_col)
  samples$slide <- sprintf("s%02d", match(slide_keys, unique(slide_keys)))
  samples <- samples %>% mutate(
    roi = sprintf("r%s", .data[[roi_key_col]]),
    segment = tolower(gsub(" ", "", .data[[segment_key_col]], fixed = TRUE)),
    aoi = sprintf("%s_%s_%s", slide, roi, segment)
  )
  stopifnot(n_distinct(samples$aoi) == nrow(samples))

  # relabel columns of count matrix with new aoi id
  colnames(counts) <- samples$aoi
  stopifnot(all(samples$aoi == colnames(counts)))

  ds <- list(samples=samples, meta=meta, counts=counts)
  return(ds)
}


#'
#' Geometric mean of a vector
#'
#' @param x vector of numeric values
geomean <- function(x) {
  return(exp(mean(log(x))))
}


calc_ks_stats <- function(x, bg, iter=10) {
  # KS test of gene vs background
  ks_dist <- NULL
  ks_pval <- NULL
  bg.length <- sum(bg)
  for (i in 1:ncol(x)) {
    print(i)
    xi <- unlist(x[!bg, i])
    xibg <- unlist(x[bg, i])

    d <- NULL
    p <- NULL
    for (j in 1:iter) {
      xj <- sample(xi, size=bg.length)
      res <- suppressWarnings(stats::ks.test(xj, xibg, alternative="less"))
      d[j] <- res$statistic
      p[j] <- res$p.value
    }
    ks_dist[i] <- median(d[j])
    ks_pval[i] <- median(p[j])
  }
  return(list(ks_dist=ks_dist, ks_pval=ks_pval))
}

#'
#' Compute qc metrics
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param meta tibble of count metadata
#' @param counts tibble of counts
#' @param bg_lod_quantile number range (0.0-1.0)
calc_bg_stats <- function(meta, counts, bg_lod_quantile) {
  # calculate AUC of gene vs background (negative probe)
  proc_auc <- function(responses, predictors) {
    myroc <- pROC::roc(responses, predictors, levels=c(TRUE, FALSE), direction="<")
    myauc <- as.numeric(myroc$auc)
    return(myauc)
  }
  bg_auc <- bind_cols(select(meta, bg), counts) %>%
    summarize(across(-bg, ~ proc_auc(bg, .x))) %>%
    unlist()

  # KS test of gene vs background
  # ks_stats <- calc_ks_stats(counts, meta$bg)
  ks_dist <- NULL
  ks_pval <- NULL
  bg <- meta$bg
  for (i in 1:ncol(counts)) {
    x <- unlist(counts[!bg, i])
    xbg <- unlist(counts[bg, i])
    res <- suppressWarnings(stats::ks.test(x, xbg, alternative="less"))
    ks_dist[i] <- res$statistic
    ks_pval[i] <- res$p.value
  }

  x <- counts[meta$bg, ]
  y <- bind_cols(
    bg_ks_dist = ks_dist,
    bg_ks_pval = ks_pval,
    bg_auc = bg_auc,
    bg_counts = colSums(x),
    bg_cpm = 1e6 * colSums(x) / colSums(counts),
    bg_geomean = apply(x, 2, geomean),
    bg_median = apply(x, 2, stats::median),
    bg_q3 = apply(x, 2, stats::quantile, 0.75),
    bg_max = apply(x, 2, max),
    bg_lod = apply(x, 2, stats::quantile, bg_lod_quantile)
  )
  return(y)
}

#'
#' Compute frac expressed genes above threshold
#'
#' @param meta tibble count metadata
#' @param counts tibble counts
#' @param bg_lod vector of limit of detection thresholds per sample
#' @returns list
calc_frac_expr <- function(meta, counts, bg_lod) {
  x <- sweep(counts, 2, bg_lod, FUN="-")
  x <- apply(x, 2, pmax, 0)
  x <- x > 0
  gene_frac_expr <- rowSums(x) / ncol(x)
  y <- x[meta$bg == FALSE,]
  sample_frac_expr <- colSums(y) / nrow(y)
  return(list(s=sample_frac_expr, g=gene_frac_expr))
}

#'
#' Process input data
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list produced by st_geomx_read
#' @param bg_lod_quantile number (0.0-1.0)
#' @returns list
#' @export
st_geomx_process_input <- function(ds, bg_lod_quantile) {
  # unpack dataset
  samples <- ds$samples
  meta <- ds$meta
  counts <- ds$counts

  # background noise statistics
  x <- calc_bg_stats(meta, counts, bg_lod_quantile)
  samples <- bind_cols(samples, x)

  # probes expressed over background
  x <- calc_frac_expr(meta, counts, samples$bg_lod)
  samples$frac_expr <- x$s
  meta$frac_expr <- x$g

  # sample metrics
  x <- counts[!meta$bg, ]
  samples <- bind_cols(samples,
    num_counts = colSums(x),
    count_mad = apply(x, 2, stats::mad),
    count_iqr = apply(x, 2, stats::IQR),
    count_geomean = apply(x, 2, geomean),
    count_median = apply(x, 2, stats::median),
    count_q3 = apply(x, 2, stats::quantile, 0.75),
    count_max = apply(x, 2, max)
  ) %>% mutate(
    snr_geomean = count_geomean / bg_geomean,
    snr_median = count_median / bg_median,
    snr_q3 = count_q3 / bg_q3
  )
  return(list(samples=samples, meta=meta, counts=counts))
}

#'
#' Merge multiple probes per gene
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list produced by st_geomx_process
#' @returns list
#' @export
st_geomx_merge_probes <- function(ds) {
  # convenience function for geomean
  geomean_fast <- function(x) {
    ifelse(length(x) == 1, x, exp(mean(log(x))))
  }

  # do not merge the negative probes, instead just move the probe name
  # to the gene name
  ybg <- bind_cols(ds$meta, ds$counts) %>%
    filter(bg) %>%
    rowwise() %>%
    summarize(
      gene = probe,
      num_probes = 1,
      probes = probe,
      bg = TRUE,
      frac_expr = frac_expr,
      across(colnames(ds$counts), geomean_fast)
    )
  bgmeta <- ybg %>% select(-colnames(ds$counts))
  bgcounts <- ybg %>% select(colnames(ds$counts))

  # merge gene probes
  y <- bind_cols(ds$meta, ds$counts) %>%
    filter(!bg) %>%
    group_by(gene) %>%
    summarize(
      num_probes = n(),
      probes = str_flatten(probe, collapse=","),
      bg = FALSE,
      frac_expr = mean(frac_expr),
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
#' Remove AOIs with no valid counts
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list produced by st_geomx_merge_probes
#' @returns list dataset
#' @export
st_geomx_rm_empty_aois <- function(ds) {
  # initial filter of "empty" samples with essential zero counts
  keep <- ds$samples$count_max > 1
  fds <- list(
    samples = filter(ds$samples, keep),
    meta = ds$meta,
    counts = ds$counts[, keep]
  )
  return(fds)
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
#' @returns list dataset
#' @export
st_geomx_set_thresholds <- function(ds,
                                    min_counts = 0,
                                    min_auc = 0.5,
                                    aoi_min_frac_expr = 0.0,
                                    gene_min_frac_expr = 0.0) {
  # apply filtering parameters to samples
  ds$samples <- ds$samples %>% mutate(
    keep = (num_counts > min_counts) & (bg_auc > min_auc) & (frac_expr > aoi_min_frac_expr)
  )
  # apply filtering parameters to genes
  ds$meta <- ds$meta %>% mutate(
    keep = frac_expr > gene_min_frac_expr,
    expressed = case_when(bg ~ "bg", !keep ~ "no", TRUE ~ "yes")
  )
  return(ds)
}


#'
#' Filter AOIs based on thresholds set by st_geomx_filter_thresholds
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds list dataset produced by st_geomx_filter_thresholds
#' @returns list dataset
#' @export
st_geomx_filter <- function(ds) {
  fds <- list(
    samples = filter(ds$samples, keep),
    meta = filter(ds$meta, keep),
    counts = ds$counts[ds$meta$keep, ds$samples$keep]
  )
  return(fds)
}


#'
#' convert to counts per million
#'
#' @param x tibble/dataframe/matrix of counts
#' @returns matrix counts-per-million
cpm <- function(x) {
  sweep(x, 2, colSums(x) / 1e6, FUN="/")
}

#'
#' convert from cpm to counts
#'
#' @param x tibble/dataframe/matrix of cpm
#' @returns matrix counts
cpm2counts <- function(x, tot_counts) {
  sweep(x, 2, tot_counts / 1e6, FUN="*")
}


#'
#' subtract background noise from each sample
#'
#' @param x count matrix
#' @param offset vector of bg noise levels per sample
#' @param min.count minimum count threshold
#' @returns matrix of background-subtracted counts
#'
background_subtract <- function(x, offset=0, min.count=1) {
  x <- sweep(x, 2, offset, FUN="-")
  x <- apply(x, 2, pmax, min.count)
  x
}


bgcorrect_norm <- function(x, bg, log=TRUE) {
  if (log) {
    x <- log2(x)
  }
  a <- sd(x[bg])^2 / sd(x[!bg])^2
  c <- a * mean(x[!bg]) - mean(x[bg])
  y <- (1-a)*x + c
  if (log) {
    y <- 2^y
  }
  return(y)
}


bgcorrect_kde <- function(x, bg, bw.adjust=1) {
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  n <- 10000
  epsbg <- 1e-15
  eps <- 1e-10

  bw <- bw.nrd0(x[!bg]) * bw.adjust
  bwbg <- bw.nrd0(x[bg]) * bw.adjust

  d <- density(x[!bg], bw=bw, from=xmin, to=xmax, n=n)
  dbg <- density(x[bg], bw=bwbg, from=xmin, to=xmax, n=n)

  # LR
  lr <- (dbg$y + epsbg) / (d$y + eps)
  # LR cannot be >1
  lr <- pmin(lr, 1)
  # LR to the left of median of noise distribution equals 1
  lr[(d$x <= median(x[bg]))] <- 1

  # scale LR to count values of x
  scaled_lr <- xmin + (xmax - xmin) * cumsum(lr) / n
  # linear interpolation to estimate noise component for each value of x
  noise <- approx(d$x, scaled_lr, xout=x)$y
  # noise can't be greater than original vector
  noise <- pmin(noise, x)
  # perform background subtraction
  y <- 2^x - 2^noise + 1

  return(list(
    y = y,
    dx = d$x,
    dy = d$y,
    dbg = dbg$y,
    lr = lr,
    scaled_lr = scaled_lr,
    bw = bw,
    bwbg = bwbg,
    eps = eps,
    epsbg = epsbg
  ))
}


bgcorrect_wfpr <- function(x, bg, bw.adjust=2, w=2) {
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  n <- 10000
  eps <- 1e-10

  bwbg <- bw.nrd0(x[bg]) * bw.adjust
  dbg <- density(x[bg], bw=bwbg, from=xmin, to=xmax, n=n)
  dfpr <- 1 - (cumsum(dbg$y + eps) / sum(dbg$y + eps))
  dwfpr <- pmin(1, w * fpr)
  scaled_wfpr <- xmin + (xmax - xmin) * cumsum(dwfpr) / n

  noise <- approx(dbg$x, scaled_wfpr, xout=x)$y
  # noise can't be greater than original vector
  noise <- pmin(noise, x)
  # perform background subtraction
  #y <- 2^(x - noise)
  y <- 2^x - 2^noise + 1

  return(list(
    y = y,
    dx = dbg$x,
    dy = dbg$y,
    dwfpr = dwfpr,
    scaled_wfpr = scaled_wfpr,
    bwbg = bwbg,
    eps = eps
  ))
}


#'
#' background correction
#'
#' @param ds list dataset
#' @returns matrix of background corrected counts-per-million
#' @export
st_geomx_bgcorrect <- function(ds, method=c("norm", "kde", "wfpr"), bw.adjust=2, w=2) {
  bg <- ds$meta$bg
  x <- ds$counts

  apply_bgcorrect_norm <- function(x) {
    return(bgcorrect_norm(x, bg))
  }
  apply_bgcorrect_kde <- function(x) {
    return(bgcorrect_kde(x, bg, bw.adjust=bw.adjust)$y)
  }
  apply_bgcorrect_wfpr <- function(x) {
    return(bgcorrect_wfpr(x, bg, bw.adjust=bw.adjust, w=w)$y)
  }

  if (method == "norm") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_norm)
  } else if (method == "kde") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_kde)
  } else if (method == "wfpr") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_wfpr)
  } else {
    stop("method not found")
  }

  x <- apply(x, 2, pmax, 1)
  return(x)
}


#'
#' quantile normalization
#'
#' @param x matrix of numeric values
#' @returns matrix of quantile normalized counts
#'
normalize_quantile <- function(x) {
  # use row sums as a tiebreaker for genes with equal counts
  row_rank <- rank(rowSums(x), ties.method="random")
  row_ord <- apply(x, 2, order, row_rank)
  # rank matrix with ties broken by row sums
  xrank <- apply(row_ord, 2, order)

  # now standard quantile norm (geomean)
  xsort <- apply(x, 2, sort)
  xgeomean <- apply(xsort, 1, geomean)

  index_to_value <- function(my_index, my_value){
    return(my_value[my_index])
  }

  xnorm <- apply(xrank, 2, index_to_value, xgeomean)
  return(xnorm)
}


#'
#' edgeR package
#' https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
#' Cite Mark Robinson
#' Scale factors as in Anders et al (2010)
#' Mark Robinson
#' Created 16 Aug 2010
#'
#' @param data matrix of counts
calc_norm_factors_rle <- function(data) {
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}


#'
#' Normalization
#'
#' background subtraction and scale/quantile normalization
#'
#' @param ds list dataset
#' @param x matrix of transformed counts
#' @param method string normalization method
#' @returns matrix of normalized counts
#' @export
st_geomx_normalize <- function(ds, x, method=c("qn", "rle", "cpm")) {

  if (method == "qn") {
    x <- normalize_quantile(cpm(x))
  } else if (method == "rle") {
    sf <- colSums(x) / 1e6
    nf <- calc_norm_factors_rle(x) / sf
    nf <- nf / geomean(nf)
    x <- sweep(x, 2, sf * nf, FUN="/")
  } else if (method == "cpm") {
    x <- cpm(x)
  } else {
    stop("method not found")
  }

  rownames(x) <- ds$meta$gene
  return(x)
}


