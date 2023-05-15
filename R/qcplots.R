#
# Spatial Transcriptomics using NanoString GeoMx DSP
#

#'
#' Set of plots to guide aoi filtering
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds dataset
#' @param min_counts number counts
#' @param min_auc number minimum area under the signal-noise ROC curve (snauc)
#' @param min_frac_expr number threshold range (0-1.0)
#' @returns ggplot object
#' @export
st_geomx_plot_aoi_filter <- function(ds,
                                     min_counts=0,
                                     min_auc=0.5,
                                     min_frac_expr=0) {
  # apply cutoffs
  s <- ds$samples %>% mutate(
    keep = ((.data$num_counts > min_counts) &
              (.data$bg_auc > min_auc) &
              (.data$frac_expr > min_frac_expr))
  )

  keep_aois <- table(s$keep)
  keep_pct <- 100 * keep_aois["TRUE"] / nrow(s)
  subtitle <- sprintf("Retained %d/%d (%.2f%%) AOIs", keep_aois["TRUE"], nrow(s), keep_pct)

  p1 <- ggplot(s, aes(x=.data$num_counts, y=.data$frac_expr, color=.data$bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_frac_expr, linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts", y="Frac detectable genes",
         title="AOI QC Metrics",
         subtitle=subtitle) +
    theme_minimal()

  p2 <- ggplot(ds$samples, aes(x=.data$num_counts, y=.data$bg_auc, color=.data$bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_auc, linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts", y="snAUC") +
    theme_minimal()

  p <- p1 + p2 + patchwork::plot_layout(guides="collect")
  return(p)
}


#'
#' Set of plots to guide gene filtering
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param min_frac_expr number threshold range (0-1.0)
#' @returns ggplot object
#' @export
st_geomx_plot_gene_filter <- function(ds, min_frac_expr=0) {
  # apply cutoffs
  meta <- ds$meta %>% mutate(
    keep = .data$frac_expr > min_frac_expr,
    expressed = case_when(.data$bg ~ "bg",
                          !.data$keep ~ "no",
                          TRUE ~ "yes")
  )

  # ROC curve analysis for expressed vs bg
  myroc <- pROC::roc(meta$bg, meta$frac_expr, levels=c(TRUE, FALSE), direction="<")
  mystats <- pROC::coords(myroc, x=min_frac_expr, input="threshold", ret="all", transpose=FALSE)

  label_auc <- paste0("AUC = ", round(myroc$auc, 3))
  label_fpr <- paste0("FPR = ", round(mystats$fpr, 3))
  label_stats <- sprintf("Threshold=%.2f Expr=%d Filtered=%d Recall=%.2f",
                         min_frac_expr, mystats$tp, mystats$fn, mystats$recall)

  # density plot showing bg, undetectable (filtered), and detected (expressed)
  p1 <- ggplot(meta, aes(x=.data$frac_expr, y=.data$expressed, fill=.data$expressed)) +
    ggridges::geom_density_ridges(alpha=0.5, scale=5, alpha=0.5) +
    #ggridges::geom_density_ridges(stat="binline", bins=50, scale=5, alpha=0.5) +
    geom_vline(xintercept = min_frac_expr, linetype="dashed", color="red") +
    scale_fill_manual(values = pals::cols25()) +
    ggridges::theme_ridges() +
    labs(x="Frac detectable above BG", y="Expressed", fill="Expressed",
         title="Gene detection distribution",
         subtitle=label_stats)

  # threshold vs false positive rate
  x <- pROC::coords(myroc, x="all", input="threshold", ret=c("threshold", "fpr", transpose = FALSE))
  p2 <- ggplot(x, aes(x=.data$threshold, y=.data$fpr)) +
    geom_line(linewidth=1, alpha=0.5, color="blue") +
    geom_vline(xintercept = min_frac_expr, linetype="dashed", color="red") +
    geom_hline(yintercept = mystats$fpr, linetype="dashed", color="red") +
    #annotate("text", label=label_fpr, x=min_frac_expr, y=mystats$fpr, hjust=0, vjust=-0.25) +
    scale_x_continuous(breaks = seq(0, 1, by=0.1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust=0.5, hjust=1, size=6)) +
    labs(x="Frac Detectable", y="FPR", title=label_fpr)

  p3 <- pROC::ggroc(myroc) +
    geom_vline(xintercept = mystats$specificity, linetype="dashed", color="red") +
    geom_hline(yintercept = mystats$sensitivity, linetype="dashed", color="red") +
    theme_minimal() +
    labs(title = label_auc)

  p <- p1 / (p2 + p3)
  return(p)
}


#'
#' boxplot showing distribution of count values
#'
#' @import ggplot2
#'
#' @param ds list dataset
#' @param x tibble expr matrix
#' @param ymax.quantile number between 0 and 1
#' @returns ggplot object
#' @export
st_geomx_plot_dist_boxplot <- function(ds, x, ymax.quantile=0.999) {
  s <- select(ds$samples, "aoi", "slide", "keep")
  y <- bind_cols(ds$meta, x)
  y <- y %>%
    tidyr::pivot_longer(s$aoi, names_to="aoi", values_to="value") %>%
    dplyr::inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi")) %>%
    dplyr::arrange(.data$slide)
  ymax <- stats::quantile(y$value, ymax.quantile)

  p <- ggplot(y, aes(x=.data$aoi, y=.data$value+1, color=.data$slide, fill=.data$slide)) +
    geom_boxplot(middle=NA, outlier.shape=NA, alpha=0.5) +
    stat_summary(fun=mean, geom="point") +
    #coord_cartesian(ylim=c(1, ymax)) +
    scale_y_continuous(trans="log10") +
    scale_color_manual(values = pals::cols25()) +
    scale_fill_manual(values = pals::cols25()) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  return(p)
}


#'
#' qc plots to assess relationship of bg versus gene
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds dataset
#' @param x count matrix
#' @param quantiles vector of quantiles [0.0-1.0]
#' @returns ggplot object
#' @export
st_geomx_plot_bg_dotplot <- function(ds, x, quantiles=c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95)) {
  # convert quantiles to string names
  quant_names <- paste0("q", round(100 * quantiles))

  # setup
  s <- select(ds$samples, "aoi", "slide", "keep", "bg_auc")
  bg <- ds$meta$bg

  # compute quantiles for bg and gene
  qbg <- summarise(x[bg,],
                   across(everything(), ~ stats::quantile(.x, quantiles)))
  qbg <- bind_cols(qbg, bg="bg", q=quant_names)
  #rownames(qbg) <- quant_names

  qgene <- summarise(x[!bg,],
                     across(everything(), ~ stats::quantile(.x, quantiles)))
  qgene <- bind_cols(qgene, bg="gene", q=quant_names)
  #rownames(qgene) <- quant_names

  # combine bg and gene quantiles
  q <- bind_rows(qbg, qgene) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  gg_dotplot_gene_vs_bg_quantiles <- function(color_by) {
    p <- ggplot(q,
                aes(x = .data$value,
                    y = stats::reorder(.data$aoi, .data$bg_auc),
                    shape = q,
                    color = {{color_by}})) +
      geom_point(alpha=0.5) +
      scale_x_log10() +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()) +
      labs(x="value", y="aoi") +
      facet_wrap(~ bg, nrow=2)
    return(p)
  }

  p1 <- gg_dotplot_gene_vs_bg_quantiles("q")
#    scale_color_manual(values = pals::cols25())
  p2 <- gg_dotplot_gene_vs_bg_quantiles("keep")
#    scale_color_manual(values = pals::cols25())
  p3 <- gg_dotplot_gene_vs_bg_quantiles("bg_auc")
    scale_color_viridis_c()
  p <- p1 + p2 + p3
  return(p)
}

#'
#' qc plots to assess relationship of bg versus gene
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param x tibble expr matrix
#' @param quantiles vector of quantiles [0.0-1.0]
#' @returns ggplot object
#' @export
st_geomx_plot_bg_scatter <- function(ds, x, quantiles=c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95)) {
  # convert quantiles to string names
  quant_names <- paste0("q", round(100 * quantiles))

  # setup
  s <- select(ds$samples, "aoi", "slide", "keep", "bg_auc")
  bg <- ds$meta$bg

  # compute quantiles for bg and gene
  qbg <- summarise(cpm(ds$counts)[bg,], across(everything(), ~ stats::quantile(.x, quantiles)))
  qbg <- bind_cols(qbg, bg="bg", q=quant_names)
  rownames(qbg) <- quant_names

  qgene <- summarise(x[!bg,], across(everything(), ~ stats::quantile(.x, quantiles)))
  qgene <- bind_cols(qgene, bg="gene", q=quant_names)
  rownames(qgene) <- quant_names

  # pivot to create scatter plots
  q <- bind_rows(qbg, qgene) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    tidyr::pivot_wider(names_from="bg", values_from="value", names_prefix="") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  gg_scatter_gene_vs_bg_quantiles <- function(color_by) {
    p <- ggplot(q, aes(x=bg, y=.data$gene, color={{color_by}})) +
      geom_point(alpha=0.6) +
      geom_abline(slope = 1, intercept = 0, color="black", linetype="dashed") +
      geom_smooth(data=filter(q, q == "q50"), method = 'lm', formula = y ~ x, se = FALSE, color="black", linetype="dashed") +
      scale_x_log10() +
      scale_y_log10() +
      theme_minimal()
    return(p)
  }

  p1 <- gg_scatter_gene_vs_bg_quantiles("q") +
    scale_color_manual(values = pals::cols25())
  p2 <- gg_scatter_gene_vs_bg_quantiles("keep") +
    scale_color_manual(values = pals::cols25())
  p3 <- gg_scatter_gene_vs_bg_quantiles("bg_auc") +
    scale_color_viridis_c()
  p <- p1 + p2 + p3
  return(p)
}


#'
#' qc plots to assess the effect of background subtraction
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param x tibble expr matrix
#' @returns ggplot object
#' @export
st_geomx_plot_bgcorrect <- function(ds, x) {
  s <- ds$samples
  m <- ds$meta

  s <- select(s, "aoi", "slide", "keep")
  y <- bind_cols(m, x)
  y <- y %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="count") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi")) %>%
    group_by(.data$aoi, .data$bg) %>%
    summarise(keep_aoi = first(.data$keep_aoi),
              slide = first(.data$slide),
              gmcpm = geomean(count))

  p1 <- ggplot(y, aes(x=.data$gmcpm,
                      y=stats::reorder(factor(.data$aoi), .data$gmcpm),
                      color=factor(.data$bg))) +
    geom_point(size=1, alpha=0.5) +
    scale_x_log10() +
    scale_color_manual(values = pals::cols25()) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(color = "BG", x="CPM", y="AOI") +
    facet_wrap(~ .data$keep_aoi)

  y <- y %>%
    tidyr::pivot_wider(names_from="bg", values_from="gmcpm", names_prefix="bg_")
  ykeep <- filter(y, .data$keep_aoi)
  bg_gene_cor <- stats::cor.test(ykeep$bg_TRUE, ykeep$bg_FALSE)$estimate
  label_cor <- sprintf("Correlation r = %.2f", bg_gene_cor)
  p2 <- ggplot(y, aes(x=.data$bg_TRUE, y=.data$bg_FALSE, color=.data$keep_aoi)) +
    geom_point(alpha=0.6) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color="black", linetype="dashed") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = pals::cols25()) +
    theme_minimal() +
    labs(color="Keep AOI",
         x="bg cpm", y="gene cpm",
         subtitle=label_cor) +
    facet_wrap(~ .data$keep_aoi)
  return(p1 / p2)
}


#'
#' qc plots to compare count densities
#'
#' @import ggplot2
#'
#' @param ds list dataset
#' @param x tibble expr matrix
#' @returns ggplot object
#' @export
st_geomx_plot_bgcorrect_density <- function(ds, x) {

  s <- select(ds$samples, "aoi", "slide", "keep")
  x <- bind_cols(ds$meta, x)
  x <- x %>%
    tidyr::pivot_longer(s$aoi, names_to="aoi", values_to="value") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))
  colnames(x)

  p <- ggplot(x, aes(x=.data$value, y=stats::reorder(factor(.data$aoi), .data$value), fill=factor(.data$bg))) +
    ggridges::stat_density_ridges(color="#ffffff00", alpha=0.4, scale=20,
                                  rel_min_height=0.001, quantile_lines=FALSE) +
    scale_x_continuous(trans="log10") +
    scale_fill_manual(values = pals::cols25()) +
    ggridges::theme_ridges() +
    theme(axis.text.y = element_blank()) +
    labs(x="Value", y="AOI", fill="BG") +
    facet_wrap(~ .data$keep_aoi, scales="free_y", ncol=1)
  return(p)
}

