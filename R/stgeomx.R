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
#' Standardize names using Slide, ROI, Segment, and AOI columns of sample
#' annotation table. Partition count table into metadata and counts.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import rlang
#'
#' @param slide_key_col sample column name containing scan
#' @param roi_key_col sample column name containing roi
#' @param segment_key_col sample column name containing segment
#' @param aoi_key_col sample column name corresponding to columns of count sheet
#' @param probe_key_col count sheet column name unique probe identifier
#' @param gene_key_col count sheet column name unique gene identifier
#' @param negprobe_gene_name name of gene corresponding to negative probes
#' @returns dataset
prepare_input <- function(samples, counts,
                          slide_key_col,
                          roi_key_col,
                          segment_key_col,
                          aoi_key_col,
                          probe_key_col,
                          gene_key_col,
                          negprobe_gene_name) {
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

  ds <- list(samples=samples,
             meta=meta,
             counts=counts)
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
#' @param negprobe_gene_name name of gene corresponding to negative probes
#' @returns dataset
#' @export
st_geomx_read_tsv <- function(sample_tsv_file,
                              count_tsv_file,
                              slide_key_col = "ScanLabel",
                              roi_key_col = "ROILabel",
                              segment_key_col = "SegmentLabel",
                              aoi_key_col = "SegmentDisplayName",
                              probe_key_col = "ProbeDisplayName",
                              gene_key_col = "TargetName",
                              negprobe_gene_name = "NegProbe-WTX") {
  # read samples (remove duplicates based on aoi key)
  samples <- read_tsv(sample_tsv_file,
                      na = c("", "NA", "NaN"),
                      show_col_types=FALSE) %>%
    distinct(.data[[aoi_key_col]], .keep_all=TRUE)
  # read count data
  counts <- read_tsv(count_tsv_file,
                     show_col_types=FALSE,
                     name_repair="minimal")
  # prepare dataset
  ds <- prepare_input(samples, counts, slide_key_col, roi_key_col,
                      segment_key_col, aoi_key_col, probe_key_col,
                      gene_key_col, negprobe_gene_name)
  return(ds)
}


#'
#' Read Nanostring GeoMx DSP data from an Excel (XLSX) file
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
st_geomx_read_xlsx <- function(xlsx_file,
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
  samples <- read_excel(xlsx_file, sheet=sample_sheet,
                        na = c("", "NaN")) %>%
    distinct(.data[[aoi_key_col]], .keep_all=TRUE)

  # read counts
  counts <- read_excel(xlsx_file, sheet=count_sheet,
                       .name_repair = "minimal")

  # prepare dataset
  ds <- prepare_input(samples, counts, slide_key_col, roi_key_col,
                      segment_key_col, aoi_key_col, probe_key_col,
                      gene_key_col, negprobe_gene_name)
  return(ds)
}


#'
#' Geometric mean of a vector
#'
#' @param x vector of numeric values
#' @returns numeric
geomean <- function(x) {
  return(exp(mean(log(x))))
}

#'
#' Geometric mean of a vector
#'
#' @param x vector of numeric values
#' @param w vector of numeric weights with length equal to x
#' @returns numeric
weighted_geomean <- function(x, w) {
  return(exp(sum(w * log(x))/sum(w)))
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
#' Preprocess input data
#'
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @param ds dataset produced by st_geomx_read
#' @param bg_lod_quantile number (0.0-1.0)
#' @returns dataset
#' @export
st_geomx_preprocess <- function(ds, bg_lod_quantile) {
  # unpack dataset
  samples <- ds$samples
  meta <- ds$meta
  counts <- ds$counts
  bg <- meta$bg

  # sample metrics
  x <- counts[!bg, ]
  xbg <- counts[bg, ]
  samples <- bind_cols(
    samples,
    # background metrics
    bg_counts = colSums(xbg),
    bg_cpm = 1e6 * colSums(xbg) / colSums(counts),
    bg_geomean = apply(xbg, 2, geomean),
    bg_median = apply(xbg, 2, stats::median),
    bg_q3 = apply(xbg, 2, stats::quantile, 0.75),
    bg_q95 = apply(xbg, 2, stats::quantile, 0.95),
    bg_max = apply(xbg, 2, max),
    bg_lod = apply(xbg, 2, stats::quantile, bg_lod_quantile),
    # gene metrics
    num_counts = colSums(x),
    complexity = apply(x, 2, function(x) length(unique(x))),
    count_mad = apply(x, 2, stats::mad),
    count_iqr = apply(x, 2, stats::IQR),
    count_geomean = apply(x, 2, geomean),
    count_median = apply(x, 2, stats::median),
    count_q3 = apply(x, 2, stats::quantile, 0.75),
    count_q95 = apply(x, 2, stats::quantile, 0.95),
    count_max = apply(x, 2, max)
  ) %>% mutate(
    snr_geomean = count_geomean / bg_geomean,
    snr_median = count_median / bg_median,
    snr_q3 = count_q3 / bg_q3,
    snr_q95 = count_q95 / bg_q95
  )

  # probes expressed over detection limit
  x <- calc_frac_expr(counts, bg, samples$bg_lod)
  samples$frac_expr <- x$s
  meta$frac_expr <- x$g

  # calculate AUC of gene vs background (negative probe)
  calc_auc <- function(responses, predictors) {
    myroc <- pROC::roc(responses, predictors, levels=c(TRUE, FALSE), direction="<")
    myauc <- as.numeric(myroc$auc)
    return(myauc)
  }
  bg_auc <- bind_cols(select(meta, bg), counts) %>%
    summarize(across(-bg, ~ calc_auc(bg, .x))) %>%
    unlist()
  samples$bg_auc <- bg_auc

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

  # calculate probe quantiles relative to negative probes
  bg_quants <- bind_cols(select(meta, bg), counts) %>%
    reframe(across(-bg, ~ compute_quantiles_kde(.x, bg)))

  return(list(samples=samples,
              meta=meta,
              counts=counts))
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
  keep_genes <- ds$meta$keep & (!ds$meta$bg)
  fds <- list(
    samples = filter(ds$samples, keep),
    meta = filter(ds$meta, keep_genes),
    counts = ds$counts[keep_genes, ds$samples$keep],
    bgmeta = filter(ds$meta, ds$meta$bg),
    bgcounts = ds$counts[ds$meta$bg, ds$samples$keep]
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
#' cluster gene probes into two groups: pure noise and signal mixed w noise
#'
#' @param x tibble/dataframe/matrix of count data
#' @param bg boolean vector indicating bg probe
cluster_noise <- function(x, bg) {
  xgene <- sort(x[!bg])
  xgene.cumsum <- cumsum(xgene)
  xgene.sum <- xgene.cumsum[length(xgene)]

  bg_center <- mean(x[bg])
  bg_dists <- cumsum(abs(xgene - bg_center))

  hi_centers <- rep(0, length(xgene))
  hi_dists <- rep(0, length(xgene))
  for (i in 1:(length(xgene)-1)) {
    hi_centers[i] <- (xgene.sum - xgene.cumsum[i])/(length(xgene) - i)
    hi_dists[i] <- sum(abs(xgene[(i+1):length(xgene)] - hi_centers[i]))
  }
  tot_dists <- bg_dists + hi_dists
  cutoff_index <- which.min(tot_dists)
  cutoff <- xgene[cutoff_index]
  hi_center <- hi_centers[cutoff_index]
  noise <- (x < cutoff)
  fpr <- sum(x[bg] > cutoff)/length(x[bg])

  # test average instead of sum
  bg_dists2 <- bg_dists / (1:length(xgene))
  hi_dists2 <- hi_dists / (length(xgene) - 1:length(xgene))
  #tot_dists2 <- bg_dists2 + hi_dists2

  return(list(noise=noise,
              x=xgene,
              bg_center=bg_center,
              hi_center=hi_center,
              bg_dist=bg_dists,
              hi_dist=hi_dists,
              bg_dist2=bg_dists2,
              hi_dist2=hi_dists2,
              fpr=fpr,
              cutoff=cutoff))
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

#'
#' background correction
#'
#' assumes normally distributed background (negative probes)
#' and gene probes
#'
#' @param x matrix of counts
#' @param bg vector of TRUE/FALSE corresponding to rows of x
#' @returns matrix of background-subtracted counts
#'
bgcorrect_norm <- function(x, bg) {
  a <- sd(x[bg])^2 / sd(x[!bg])^2
  c <- a * mean(x[!bg]) - mean(x[bg])
  y <- (1-a)*x + c
  return(y)
}


#'
#' background correction using quantile-quantile approach
#'
#' @param x matrix of counts
#' @param bg vector of TRUE/FALSE corresponding to rows of x
#' @returns list where element 'y' contains corrected counts
#'
bgcorrect_qq <- function(x, bg, bw.adjust=1, bg.quant=0.5) {
  # setup
  n <- 2^13
  xlog <- log2(x)

  # select bandwidth
  bwbg <- bw.nrd0(xlog[bg]) * bw.adjust
  bw <- bw.nrd0(xlog[!bg]) * bw.adjust
  bw <- max(bw, bwbg)

  # smooth ecdf functions for gene probes and neg probes
  d <- density(xlog[!bg], bw=bw, n=n)
  dcdf <- cumsum(d$y) / sum(d$y)

  dbg <- density(xlog[bg], bw=bw, n=n)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  q <- approx(d$x, dcdf, xout=xlog, ties="ordered")$y
  q <- pmax(bg.quant, q)

  log2noise <- approx(dbgcdf, dbg$x, xout=q, ties="ordered")$y
  noise <- 2^log2noise

  # ensure noise is not greater than (x - 1)
  noise <- pmin(noise, x - 1)
  y <- x - noise
  # ensure ymin == 1
  y <- y - (min(y) - 1)
  # ensure isotonic behavior
  ir <- isoreg(x, y)
  yf <- ir$yf[order(ir$ord)]

  res <- list(
    y = y,
    yf = yf,
    noise = noise,
    dx = d$x,
    dy = d$y,
    dcdf = dcdf,
    dbgx = dbg$x,
    dbgy = dbg$y,
    dbgcdf = dbgcdf,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  )
}


bgcorrect_kdeppv <- function(x, bg, bw.adjust=2, bg.quant=0.5) {
  # setup
  eps <- 1e-10
  n <- 2^13
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  expressed <- !bg & (x > quantile(x[bg], bg.quant))

  # select bandwidth
  bw.gene <- bw.nrd0(x[expressed])
  bw.bg <- bw.nrd0(x[bg])
  bw <- max(bw.gene, bw.bg) * bw.adjust

  d <- density(x[expressed], bw=bw, from=xmin, to=xmax, n=n)
  dbg <- density(x[bg], bw=bw, from=xmin, to=xmax, n=n)

  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)



  # convert density to ecdf
  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  # calc false omission rate (FOR) as FN / PN
  dppv <- (1 - dcdf + eps) / (1 - dcdf + 1 - dbgcdf + eps)
  dppv <- (dppv - min(dppv)) / (max(dppv) - min(dppv))

  # ensure isotonic behavior
  dppv.ir <- isoreg(d$x, dppv)$yf
  # interpolate
  ppv <- approx(d$x, dppv.ir, xout=x, ties="ordered")$y
  # transform
  y <- (2^x) * ppv
  y <- pmax(1, y)

  return(list(
    y = y,
    dx = d$x,
    dcdf = dcdf,
    dbgcdf = dbgcdf,
    dp = dppv.ir,
    bw = bw,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  ))
}


bgcorrect_kdefor <- function(x, bg, bw.adjust=2, bg.quant=0.5) {
  # setup
  eps <- 1e-10
  n <- 2^13
  x <- log2(x)
  xmin <- min(x)
  xmax <- max(x)
  expressed <- !bg & (x > quantile(x[bg], bg.quant))

  # select bandwidth
  bw.gene <- bw.nrd0(x[expressed])
  bw.bg <- bw.nrd0(x[bg])
  bw <- max(bw.gene, bw.bg) * bw.adjust

  # use kde to produce a smooth density estimate in log space
  d <- density(x[expressed], bw=bw, from=xmin, to=xmax, n=n)
  dbg <- density(x[bg], bw=bw, from=xmin, to=xmax, n=n)

  # convert density to ecdf
  dcdf <- cumsum(d$y) / sum(d$y)
  dbgcdf <- cumsum(dbg$y) / sum(dbg$y)

  # calc false omission rate (FOR) as FN / PN
  dp <- dcdf / (dbgcdf + dcdf + eps)
  dp <- (dp - min(dp)) / (max(dp) - min(dp))

  # ensure isotonic behavior
  dp.ir <- isoreg(d$x, dp)$yf

  # interpolate
  p <- approx(d$x, dp.ir, xout=x, ties="ordered")$y
  #f <- splinefun(d$x, dp.ir, method="natural", ties="ordered")
  #p <- f(x)

  # transform
  y <- (2^x) * p
  y <- pmax(1, y)

  return(list(
    y = y,
    dx = d$x,
    dcdf = dcdf,
    dbgcdf = dbgcdf,
    dp = dp.ir,
    bw = bw,
    bw.adjust = bw.adjust,
    bg.quant = bg.quant
  ))
}


#'
#' background correction
#'
#' @param ds list dataset
#' @returns matrix of background corrected counts-per-million
#' @export
st_geomx_bgcorrect <- function(ds, method=c("qq", "norm", "bgsub", "kdeppv", "kdefor", "none"),
                               bw.adjust=1, bg.quant=0.5) {
  bg <- ds$meta$bg
  x <- ds$counts

  apply_bgcorrect_qq <- function(x) {
    y <- bgcorrect_qq(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }
  apply_bgcorrect_norm <- function(x) {
    y <- bgcorrect_norm(x, bg)
    y <- pmax(y, 1)
    return(y)
  }
  apply_bgcorrect_bgsub <- function(x) {
    xsub <- 2^quantile(log2(x)[bg], bg.quant)
    y <- pmax(1, x - xsub)
    return(y)
  }
  apply_bgcorrect_kdeppv <- function(x) {
    y <- bgcorrect_kdeppv(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }
  apply_bgcorrect_kdefor <- function(x) {
    y <- bgcorrect_kdefor(x, bg, bw.adjust=bw.adjust, bg.quant=bg.quant)$y
    return(y)
  }

  if (method == "none") {
    return(x)
  } else if (method == "qq") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_qq)
  } else if (method == "norm") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_norm)
  } else if (method == "bgsub") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_bgsub)
  } else if (method == "kdeppv") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_kdeppv)
  } else if (method == "kdefor") {
    x <- apply(x, MARGIN=2, FUN=apply_bgcorrect_kdefor)
  } else {
    stop("method not found")
  }
  return(x)
}


#'
#' tie-breaking quantile normalization
#'
#' @import dplyr
#'
#' @param x matrix of raw counts
#' @returns matrix of normalized relative counts
#'
# tbqnorm <- function(x) {
#   # convert to relative counts (counts-per-million)
#   xcpm <- apply(x-1, 2, function(x) { 1e6 * x / sum(x) })
#   # rank genes across samples by total relative expression
#   global_ranks <- rank(rowSums(xcpm), ties.method="random")
#   # use correlation to compute sample-sample similarity matrix
#   xcor <- cor(xcpm, method="spearman")
#   # exclude samples with negative correlation
#   xcor[xcor < 0] <- NA
#   # ordered list of nearest neighbors for each sample
#   nn_list <- apply(xcor, 2, order, na.last=NA, decreasing=TRUE)
#   # starting ranks (with ties present)
#   xrank_init <- as.data.frame(apply(x, 2, rank, ties.method="min"))
#   #xrank_init <- reframe(as_tibble(x), across(everything(), min_rank))
#   # reordering procedure to break ties for each sample
#   xrank <- matrix(nrow=nrow(x), ncol=ncol(x))
#   for(i in 1:ncol(x)) {
#     nn_ordered_cols <- xrank_init[, nn_list[,i]]
#     nn_ordered_cols <- bind_cols(nn_ordered_cols, tiebrk_ranks=global_ranks)
#     xrank[,i] <- do.call(order, unname(nn_ordered_cols))
#   }
#   # convert ordering to ranks
#   xrank <- apply(xrank, 2, order)
#
#   # compute quantiles
#   xquantiles <- apply(xcpm, 2, sort)
#   xquantiles <- apply(xquantiles, 1, geomean)
#   # renormalize quantiles such that columns sum to one million (e.g. cpm)
#   xquantiles <- 1e6 * (xquantiles / sum(xquantiles))
#   # apply quantile normalization procedure using new ranks
#   xnorm <- apply(xrank, 2, function(x) { xquantiles[x] })
# }

# tbqnorm <- function(x,
#                     final_tiebreak=c("random", "global"),
#                     noise_weights=NULL,
#                     noise_alpha=1,
#                     smooth_alpha=0) {
#   # parse arguments
#   final_tiebreak <- match.arg(final_tiebreak)
#
#   # if weights not specified, make weights equal
#   if (is.null(noise_weights)) {
#     noise_weights <- rep(1, ncol(x))
#   }
#   noise_weights <- noise_weights ^ noise_alpha
#
#   # convert to relative counts (counts-per-million)
#   xcpm <- apply(x-1, 2, function(x) { 1e6 * x / sum(x) })
#
#   # compute sample similarities
#   #xquantiles <- apply(xcpm, 2, sort)
#   #xsimilarity <- as.matrix(stats::dist(t(xquantiles), method="manhattan"))
#   #xsimilarity <- 1 - (xsimilarity / (maxdist * 1e6))
#
#   # rank genes across samples by total relative expression
#   global_ranks <- rank(rowSums(xcpm), ties.method="random")
#
#   # use correlation to compute sample-sample similarity matrix
#   xcor <- cor(xcpm, method="spearman")
#   # exclude samples with negative correlation
#   xcor[xcor < 0] <- NA
#   # ordered list of nearest neighbors for each sample
#   nn_list <- apply(xcor, 2, order, na.last=NA, decreasing=TRUE)
#
#   # starting ranks (with ties present)
#   xrank_init <- as.data.frame(apply(x, 2, rank, ties.method="min"))
#   #xrank_init <- reframe(as_tibble(x), across(everything(), min_rank))
#
#   # reordering procedure to break ties for each sample
#   xrank <- matrix(nrow=nrow(x), ncol=ncol(x))
#   for(i in 1:ncol(x)) {
#     # if unable to break ties using nearest neighbors will defer to
#     # either global ranks or random ranks
#     if (final_tiebreak == "global") {
#       tiebrk_ranks <- global_ranks
#     } else {
#       tiebrk_ranks <- rank(xrank_init[,i], ties.method="random")
#     }
#     nn_ordered_cols <- xrank_init[, nn_list[[i]]]
#     nn_ordered_cols <- bind_cols(nn_ordered_cols, tiebrk_ranks=tiebrk_ranks)
#     xrank[,i] <- do.call(order, unname(nn_ordered_cols))
#   }
#   # convert ordering to ranks
#   xrank <- apply(xrank, 2, order)
#
#   # compute weighted quantiles
#   xquantiles <- apply(xcpm, 2, sort)
#   xquantiles <- apply(xquantiles, 1, geomean)
#   #xquantiles <- apply(xquantiles, 1, weighted_geomean, noise_weights)
#   # renormalize quantiles such that columns sum to one million (e.g. cpm)
#   xquantiles <- 1e6 * (xquantiles / sum(xquantiles))
#
#   # apply quantile normalization procedure using new ranks
#   xnorm <- apply(xrank, 2, function(x) { xquantiles[x] })
#   return(xnorm)
# }


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
#' edgeR package
#' TMM between two libraries
#' Mark Robinson
#' Copied from https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
#'
calc_norm_factors_tmm <- function(x, logratioTrim=0.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {

  num_counts <- colSums(x)
  f75 <- pmax(1, apply(x, 2, quantile, 0.75)) / num_counts
  nsamples <- ncol(x)

  # choose reference column to normalize against
  if(median(f75) < 1e-20) {
    refColumn <- which.max(colSums(sqrt(x)))
  } else {
    refColumn <- which.min(abs(f75-mean(f75)))
  }
  ref <- unlist(x[, refColumn])
  nR <- num_counts[refColumn]

  calc_factor_tmm <- function(i) {
    obs <- unlist(x[, i])
    nO <- num_counts[i]

    logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
    absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
    v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance

    #	remove infinite values, cutoff based on A
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]

    if(max(abs(logR)) < 1e-6) return(1)

    #	taken from the original mean() function
    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS

    #	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    #	a fix from leonardo ivan almonacid cardenas, since rank() can return
    #	non-integer values when there are a lot of ties
    keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

    if(doWeighting) {
      f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
    } else {
      f <- mean(logR[keep], na.rm=TRUE)
    }

    #	Results will be missing if the two libraries share no features with positive counts
    #	In this case, return unity
    if(is.na(f)) f <- 0
    2^f
  }

  nf <- rep_len(NA_real_, nsamples)
  for(i in 1:nsamples) {
    obs <- x[, i]
    nO <- num_counts[i]
    nf[i] <- calc_factor_tmm(i)
  }
  return(nf)
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
st_geomx_normalize <- function(ds, x, method=c("qn", "rle", "cpm", "q3", "tmm", "none")) {
  # check method
  method <- match.arg(method)

  if (method == "none") {
    # do nothing
    x <- x
  } else if (method == "qn") {
    x <- normalize_quantile(cpm(x))
  } else if (method == "cpm") {
    x <- cpm(x)
  } else {
    # both Q3, RLE, and TMM produce normalization factors
    sf <- colSums(x) / 1e6

    if (method == "rle") {
      nf <- calc_norm_factors_rle(x) / sf

    } else if (method == "q3") {
      nf <- pmax(1, apply(x, 2, quantile, 0.75)) / sf
    } else if (method == "tmm") {
      nf <- calc_norm_factors_tmm(x)
    }
    # norm factors should multiply to 1
    nf <- nf/exp(mean(log(nf)))
    x <- sweep(x, 2, sf * nf, FUN="/")
  }

  rownames(x) <- ds$meta$gene
  return(x)
}

