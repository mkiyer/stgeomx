#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


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
#' @param tot_counts vector of total counts per sample
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
  apply(data, 2, function(u) stats::median((u/gm)[gm > 0]))
}

#'
#' edgeR package
#' TMM between two libraries
#' Mark Robinson
#' Copied from https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R
#'
#' @param x count matrix
#' @param logratioTrim see edgeR documentation
#' @param sumTrim see edgeR documentation
#' @param doWeighting see edgeR documentation
#' @param Acutoff see edgeR documentation
calc_norm_factors_tmm <- function(x, logratioTrim=0.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {

  num_counts <- colSums(x)
  f75 <- pmax(1, apply(x, 2, stats::quantile, 0.75)) / num_counts
  nsamples <- ncol(x)

  # choose reference column to normalize against
  if(stats::median(f75) < 1e-20) {
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
#' @param x matrix of transformed counts
#' @param method string normalization method
#' @returns matrix of normalized counts
#' @export
st_geomx_normalize <- function(x, method=c("qn", "rle", "cpm", "q3", "tmm", "none")) {
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
      nf <- pmax(1, apply(x, 2, stats::quantile, 0.75)) / sf
    } else if (method == "tmm") {
      nf <- calc_norm_factors_tmm(x)
    }
    # norm factors should multiply to 1
    nf <- nf/exp(mean(log(nf)))
    x <- sweep(x, 2, sf * nf, FUN="/")
  }

  return(x)
}
