#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


#'
#' QC plots for background correction
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds dataset
#' @param aoi name of aoi
#' @param bw.adjust integer adjustment factor for bandwidth function
#' @param bg.quant background quantile to subtract to the counts
#' @param bg.lod.quant limit of detection quantile
#' @returns ggplot object
#' @export
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_aoi_qc(example_ds, "s01_r001_fullroi")
#'
plot_aoi_qc <- function(ds, aoi,
                        bw.adjust=2,
                        bg.quant=0.5,
                        bg.lod.quant=0.9) {
  x <- pull(ds$counts, aoi)
  bg <- ds$meta$bg

  # ROC curve analysis for gene vs bg
  myroc <- pROC::roc(bg, x, levels=c(TRUE, FALSE), direction="<")
  label_auc <- paste0("AUC = ", round(myroc$auc, 3))
  p1 <- pROC::ggroc(myroc) +
    theme_bw() +
    labs(title=aoi, subtitle=label_auc)

  # density plot of magnitude of gene vs bg counts (snr)
  bg.lod <- stats::quantile(x[bg], bg.lod.quant)
  snr <- log2(x+1) - log2(bg.lod+1)
  frac_expr <- sum(snr > 0) / length(snr)
  label_frac_expr <- paste0("Frac > LOD = ", round(frac_expr, 3))
  d <- stats::density(snr)
  probs <- c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
  quantiles <- stats::quantile(snr, prob=probs)
  bins <- findInterval(d$x, quantiles)+1
  tbl <- tibble(x=d$x, y=d$y, q=factor(c(0, probs)[bins]))
  p2 <- ggplot(tbl, aes(.data$x, .data$y, fill=.data$q)) +
    geom_line() +
    geom_ribbon(aes(ymin=0, ymax=.data$y)) +
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    scale_fill_viridis_d() +
    scale_x_continuous(breaks = scales::breaks_extended(n=10)) +
    theme_bw() +
    labs(x="log2(SNR)", y="Density", subtitle=label_frac_expr)

  # background correction
  res <- bgcorrect_qq(x, bg, bw.adjust, bg.quant)
  tbl <- bind_rows(
    bind_cols(name="raw", x=x, value=x, bg=bg),
    bind_cols(name="noise", x=x, value=res$noise, bg=bg),
    bind_cols(name="corrected", x=x, value=res$y, bg=bg)
  )
  p3 <- ggplot(tbl, aes(x=.data$x+1, y=.data$value+1, color=.data$name)) +
    geom_line(alpha=0.5) +
    geom_point(data=filter(tbl, bg==TRUE), alpha=0.5) +
    geom_vline(xintercept=res$bg.loq, color="red", linetype="dashed", alpha=0.5) +
    scale_color_manual(values = pals::cols25()) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(x="Raw counts + 1", y="Counts + 1", color="BG Correction")

  p4 <- ggplot(tbl, aes(x=.data$value+1, y=.data$name, fill=.data$bg)) +
    ggridges::stat_density_ridges(alpha=0.5, scale=2, rel_min_height=0.001) +
    scale_fill_manual(values = pals::cols25()) +
    scale_x_log10() +
    theme_bw() +
    labs(x="Counts + 1", y="BG Correction", fill="Background")

  return((p1 + p2) / p3 / p4)
}


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
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_aoi_filter(example_ds)
#'
plot_aoi_filter <- function(ds,
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

  p1 <- ggplot(s, aes(x=.data$num_counts+1, y=.data$frac_expr, color=.data$bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_frac_expr, linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts+1", y="Frac detectable genes", color="BG CPM",
         title="AOI QC Metrics",
         subtitle=subtitle) +
    theme_minimal()

  p2 <- ggplot(ds$samples, aes(x=.data$num_counts+1, y=.data$bg_auc, color=.data$bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_auc, linetype = "dashed", color = "red") +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts+1", y="snAUC", color="BG CPM") +
    theme_minimal()

  p3 <- ggplot(ds$samples, aes(x=.data$bg_auc, y=.data$frac_expr, color=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_auc, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_frac_expr, linetype = "dashed", color = "red") +
    scale_color_manual(values = pals::cols25()) +
    scale_x_log10() +
    labs(x="snAUC", y="Frac detectable genes", color="QC Filter") +
    theme_minimal()

  p <- p1 + p2 + p3 + patchwork::plot_layout(guides="collect")
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
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_gene_filter(example_ds)
#'
plot_gene_filter <- function(ds, min_frac_expr=0) {
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
    ggridges::geom_density_ridges(scale=5, alpha=0.5) +
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
#' qc plot to assess relationship of bg versus gene
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param quantiles vector of quantiles [0.0-1.0]
#' @returns ggplot object
#' @export
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_bgcorrect1(example_ds)
#'
plot_bgcorrect1 <- function(ds, quantiles=c(0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)) {
  # labels for plot legend
  quant_names <- paste0("q", round(100 * quantiles))
  s <- select(ds$samples, "aoi", "slide", "keep")
  x <- ds$x

  # bg cpm before bgcorrect
  bg <- ds$meta$bg
  xbg <- filter(cpm(ds$counts), bg)
  qbg <- xbg %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="bg", quant=quant_names)

  # gene cpm before bgcorrect
  xgene <- filter(ds$counts, !bg)
  xgene <- as_tibble(cpm(xgene))
  qgene <- xgene %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="gene", quant=quant_names)

  # cpm after bgcorrect
  y <- as_tibble(cpm(x[!bg, ]))
  qy <- y %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="y", quant=quant_names)

  # combine bg and gene quantiles
  q <- bind_rows(qbg, qgene, qy) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    tidyr::pivot_wider(names_from="name", values_from="value", names_prefix="") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  p1 <- ggplot(q, aes(x=.data$bg+1, y=.data$gene+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()

  p2 <- ggplot(q, aes(x=.data$bg+1, y=.data$y+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()
  return(p1 + p2)
}


#'
#' plot total background noise (x axis) versus gene quantiles (y axis)
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param quantiles vector of quantiles [0.0-1.0]
#' @returns ggplot object
#' @export
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_bgcorrect2(example_ds)
#'
plot_bgcorrect2 <- function(ds, quantiles=c(0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)) {
  s <- select(ds$samples, "aoi", "slide", "keep")
  x <- ds$x
  bg <- ds$meta$bg

  # compute background level per aoi
  s$bg_cpm <- 1e6 * colSums(ds$counts[bg,]) / colSums(ds$counts)
  # labels for plot legend
  quant_names <- paste0("q", round(100 * quantiles))

  # gene cpm before bgcorrect
  xgene <- filter(ds$counts, !bg)
  xgene <- as_tibble(cpm(xgene))
  qgene <- xgene %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="gene", quant=quant_names)

  # cpm after bgcorrect
  y <- as_tibble(cpm(x[!bg, ]))
  qy <- y %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="y", quant=quant_names)

  # combine bg and gene quantiles
  q <- bind_rows(qgene, qy) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    tidyr::pivot_wider(names_from="name", values_from="value", names_prefix="") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  p1 <- ggplot(q, aes(x=.data$bg_cpm, y=.data$gene+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()

  p2 <- ggplot(q, aes(x=.data$bg_cpm, y=.data$y+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()
  return(p1 + p2)
}



#'
#' plot total background noise (x axis) versus gene quantiles (y axis)
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param ds list dataset
#' @param quantiles vector of quantiles [0.0-1.0]
#' @returns ggplot object
#' @export
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_bgcorrect3(example_ds)
#'
plot_bgcorrect3 <- function(ds, quantiles=c(0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99)) {
  s <- select(ds$samples, "aoi", "keep")
  x <- ds$x
  bg <- ds$meta$bg

  # compute background level per aoi
  s$bg_cpm <- 1e6 * colSums(ds$counts[bg,]) / colSums(ds$counts)
  # labels for plot legend
  quant_names <- paste0("q", round(100 * quantiles))

  # raw gene counts
  xgene <- filter(ds$counts, !bg)
  qgene <- xgene %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="gene", quant=quant_names)

  # normalized counts
  y <- as_tibble(x[!bg, ])
  qy <- y %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(name="y", quant=quant_names)

  # combine bg and gene quantiles
  q <- bind_rows(qgene, qy) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    tidyr::pivot_wider(names_from="name", values_from="value", names_prefix="") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  p1 <- ggplot(q, aes(x=.data$bg_cpm, y=.data$gene+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()

  p2 <- ggplot(q, aes(x=.data$bg_cpm, y=.data$y+1, color=.data$quant, shape=.data$keep)) +
    geom_point(alpha=0.6) +
    geom_smooth(data=filter(q, .data$keep), linetype='dashed', method = 'lm',
                formula = y ~ x, se=FALSE, alpha=0.6) +
    scale_color_viridis_d(option="plasma") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal()
  return(p1 + p2)
}



#'
#' boxplot showing distribution of count values
#'
#' @import ggplot2
#'
#' @param ds list dataset
#' @param bg.lod.quantile number between 0 and 1
#' @returns ggplot object
#' @export
#'
#' @examples
#' data(example_ds, package = "stgeomx")
#' example_ds <- preprocess(example_ds)
#' plot_expr_dist(example_ds)
#'
plot_expr_dist <- function(ds, bg.lod.quantile=0.9) {
  s <- select(ds$samples, "aoi", "slide", "keep")
  x <- ds$x
  bg <- ds$meta$bg

  # compute background level per aoi
  s$bg_lod <- apply(x[bg,], 2, stats::quantile, bg.lod.quantile)

  # quantiles
  quantiles <- c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99)
  quant_names <- paste0("q", round(100 * quantiles))

  q <- as_tibble(x) %>%
    filter(!bg) %>%
    reframe(across(everything(), ~ stats::quantile(.x, quantiles))) %>%
    mutate(quant=quant_names) %>%
    tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="value") %>%
    tidyr::pivot_wider(names_from="quant", values_from="value", names_prefix="") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  data1 <- q %>%
    rowwise() %>%
    arrange(.data$keep, .data$q50) %>%
    mutate(x=factor(.data$aoi, .data$aoi)) %>%
    ungroup()

  data2 <- data1 %>%
    select("x", "slide", "keep", "q50", "q90", "q95", "q99") %>%
    tidyr::pivot_longer(cols = c("q50", "q90", "q95", "q99"),
                        names_to="quant", values_to="value")

  p <- ggplot() +
    geom_errorbar(data=data1, aes(x=.data$x, ymin=.data$q25+1, ymax=.data$q75+1, color=.data$keep), alpha=0.75) +
    #geom_segment(data=data1, aes(x=x, xend=x, y=q25, yend=q75), color="grey", alpha=0.75) +
    geom_point(data=data1, aes(x=.data$x, y=.data$bg_lod+1), color="black", shape=4, alpha=0.75) +
    geom_point(data=data2, aes(x=.data$x, y=.data$value+1, fill=.data$quant), shape=21, alpha=0.75) +
    scale_color_manual(values=pals::cols25()) +
    scale_fill_viridis_d() +
    scale_y_log10() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    facet_grid(~ slide, scales="free_x") +
    labs(x = "AOI", y="Gene Expression")

  return(p)
}



# st_geomx_plot_bgcorrect <- function(ds, x) {
#   s <- ds$samples
#   m <- ds$meta
#
#   s <- select(s, "aoi", "slide", "keep")
#   y <- bind_cols(m, x)
#   y <- y %>%
#     tidyr::pivot_longer(colnames(x), names_to="aoi", values_to="count") %>%
#     inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi")) %>%
#     group_by(.data$aoi, .data$bg) %>%
#     summarise(keep_aoi = first(.data$keep_aoi),
#               slide = first(.data$slide),
#               gmcpm = geomean(count))
#
#   p1 <- ggplot(y, aes(x=.data$gmcpm,
#                       y=stats::reorder(factor(.data$aoi), .data$gmcpm),
#                       color=factor(.data$bg))) +
#     geom_point(size=1, alpha=0.5) +
#     scale_x_log10() +
#     scale_color_manual(values = pals::cols25()) +
#     theme_minimal() +
#     theme(axis.text.y = element_blank(),
#           panel.grid.major.y = element_blank(),
#           panel.grid.minor.y = element_blank()) +
#     labs(color = "BG", x="CPM", y="AOI") +
#     facet_wrap(~ .data$keep_aoi)
#
#   y <- y %>%
#     tidyr::pivot_wider(names_from="bg", values_from="gmcpm", names_prefix="bg_")
#   ykeep <- filter(y, .data$keep_aoi)
#   bg_gene_cor <- stats::cor.test(ykeep$bg_TRUE, ykeep$bg_FALSE)$estimate
#   label_cor <- sprintf("Correlation r = %.2f", bg_gene_cor)
#   p2 <- ggplot(y, aes(x=.data$bg_TRUE, y=.data$bg_FALSE, color=.data$keep_aoi)) +
#     geom_point(alpha=0.6) +
#     geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color="black", linetype="dashed") +
#     scale_x_log10() +
#     scale_y_log10() +
#     scale_color_manual(values = pals::cols25()) +
#     theme_minimal() +
#     labs(color="Keep AOI",
#          x="bg cpm", y="gene cpm",
#          subtitle=label_cor) +
#     facet_wrap(~ .data$keep_aoi)
#   return(p1 / p2)
# }



# st_geomx_plot_bgcorrect_density <- function(ds, x) {
#
#   s <- select(ds$samples, "aoi", "slide", "keep")
#   x <- bind_cols(ds$meta, x)
#   x <- x %>%
#     tidyr::pivot_longer(s$aoi, names_to="aoi", values_to="value") %>%
#     inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))
#   colnames(x)
#
#   p <- ggplot(x, aes(x=.data$value, y=stats::reorder(factor(.data$aoi), .data$value), fill=factor(.data$bg))) +
#     ggridges::stat_density_ridges(color="#ffffff00", alpha=0.4, scale=20,
#                                   rel_min_height=0.001, quantile_lines=FALSE) +
#     scale_x_continuous(trans="log10") +
#     scale_fill_manual(values = pals::cols25()) +
#     ggridges::theme_ridges() +
#     theme(axis.text.y = element_blank()) +
#     labs(x="Value", y="AOI", fill="BG") +
#     facet_wrap(~ .data$keep_aoi, scales="free_y", ncol=1)
#   return(p)
# }

