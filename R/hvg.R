#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


#'
#' Model and rank highly variable genes in dataset
#'
#' Models mean-variance relationship using loess. Ratio of actual variance
#' versus modeled variance is used to rank genes
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import rlang
#'
#' @param x matrix/data frame of normalized log2 transformed gene expression,
#' where gene identifiers are present in the row names of the matrix
#' @param loess.span span parameter for loess model
#' @param min_var_quantile [0.0-1.0] used as lower bound for predicted
#' variance by loess model (in rare cases the predicted variance by loess
#' can be less than zero or very close to zero, resulting in high variance
#' ratios)
#' @returns tibble object with per-gene parameters
#' @export
#'
get_hvg_tbl <- function(x, loess.span=0.25, min_var_quantile=0.01) {
  tbl <- dplyr::bind_cols(
    gene = rownames(x),
    u = rowMeans(x),
    v = apply(x, 1, function(x) { stats::var(x) })
  )
  lovar <- stats::loess(v ~ u, data=tbl, span=loess.span)
  # stabilize fitted variance (loess can produce negative values
  # and values very close to zero), setting it to a low positive number
  tbl$vloess <- stats::predict(lovar)
  tbl$vminbound <- stats::quantile(tbl$v[tbl$v > 0], min_var_quantile)
  tbl$vloess_bounded <- pmax(tbl$vloess, tbl$vminbound)
  # compute variance ratios
  tbl <- dplyr::mutate(tbl,
    vratio = .data$v / .data$vloess_bounded,
    rank = dplyr::min_rank(-.data$vratio)
  )
  return(tbl)
}


#'
#' plot of gene mean versus variance with line showing predicted (expected)
#' variance by loess model. dashed line shows lower bound to predicted
# 'variance
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @param x matrix/data frame (normalized log2 gene expression)
#' @param nhvg number of top highly variable genes to color differently
#' @param nlabel number of top highly variable genes to label
#' @param loess.span span parameter for loess model
#' @param min_var_quantile [0.0-1.0] variance quantile to use as lower bound
#' @returns ggplot object
#' @export
#'
plot_hvg_model <- function(x, nhvg=500, nlabel=100,
                           loess.span=0.25,
                           min_var_quantile=0.01) {
  tbl <- get_hvg_tbl(x, loess.span, min_var_quantile)
  label_data <- dplyr::filter(tbl, rank <= nlabel)
  top_data <- dplyr::filter(tbl, rank <= nhvg)
  bot_data <- dplyr::filter(tbl, rank > nhvg)
  p <- ggplot(tbl, aes(.data$u, .data$v)) +
    geom_point(data=bot_data, alpha=0.5, size=1, color="#aaaaaa") +
    geom_point(data=top_data, alpha=0.8, size=2, color="#ffaaaa") +
    geom_line(aes(.data$u, .data$vloess), color="blue") +
    geom_line(aes(.data$u, .data$vminbound), color="red", linetype="dashed") +
    geom_text_repel(data=label_data, aes(label=.data$gene), color="black", size=3, max.overlaps=Inf) +
    theme_minimal() +
    labs(x="mean", y="variance", title="HVG Model")
  return(p)
}
