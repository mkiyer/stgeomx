#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


#'
#' Set of plots to guide aoi filtering
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param list dataset
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
    keep = (num_counts > min_counts) & (bg_auc > min_auc) & (frac_expr > min_frac_expr)
  )

  keep_aois <- table(s$keep)
  keep_pct <- 100 * keep_aois["TRUE"] / nrow(s)
  subtitle <- sprintf("Retained %d/%d (%.2f%%) AOIs", keep_aois["TRUE"], nrow(s), keep_pct)

  p1 <- ggplot(s, aes(x=num_counts, y=frac_expr, color=bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_frac_expr, linetype = "dashed", color = "red") +
#    scale_color_manual(values = pals::cols25()) +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts", y="Frac detectable genes",
         title="AOI QC Metrics",
         subtitle=subtitle) +
    theme_minimal()

  p2 <- ggplot(ds$samples, aes(x=num_counts, y=bg_auc, color=bg_cpm)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_auc, linetype = "dashed", color = "red") +
#    scale_color_manual(values = pals::cols25()) +
    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts", y="snAUC") +
    theme_minimal()

  p <- p1 + p2 + plot_layout(guides="collect")
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
    keep = frac_expr > min_frac_expr,
    expressed = case_when(bg ~ "bg", !keep ~ "no", TRUE ~ "yes")
  )

  # ROC curve analysis for expressed vs bg
  myroc <- pROC::roc(meta$bg, meta$frac_expr, levels=c(TRUE, FALSE), direction="<")
  mystats <- pROC::coords(myroc, x=min_frac_expr, input="threshold", ret="all", transpose=FALSE)

  label_auc <- paste0("AUC = ", round(myroc$auc, 3))
  label_fpr <- paste0("FPR = ", round(mystats$fpr, 3))
  label_stats <- sprintf("Threshold=%.2f Expr=%d Filtered=%d Recall=%.2f",
                         min_frac_expr, mystats$tp, mystats$fn, mystats$recall)

  # density plot showing bg, undetectable (filtered), and detected (expressed)
  p1 <- ggplot(meta, aes(x=frac_expr, y=expressed, fill=expressed)) +
    #geom_density_ridges(alpha=0.7, scale=5) +
    ggridges::geom_density_ridges(stat="binline", bins=50, scale=5, alpha=0.5) +
    geom_vline(xintercept = min_frac_expr, linetype="dashed", color="red") +
    scale_fill_manual(values = pals::cols25()) +
    ggridges::theme_ridges() +
    labs(x="Frac detectable above BG", y="Expressed", fill="Expressed",
         title="Gene detection distribution",
         subtitle=label_stats)

  # threshold vs false positive rate
  x <- pROC::coords(myroc, x="all", input="threshold", ret=c("threshold", "fpr", transpose = FALSE))
  p2 <- ggplot(x, aes(x=threshold, y=fpr)) +
    geom_line(size=1, alpha=0.5, color="blue") +
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


#' compare median gene to bg cpm
#'
#' @param x list dataset
#' @returns ggplot object
dotplot_aoi_bg_vs_signal <- function(x, cpm_quantile = 0.75) {
  y <- x %>%
    group_by(aoi, bg) %>%
    summarise(keep_aoi = first(keep_aoi),
              slide = first(slide),
              qcpm = stats::quantile(count, cpm_quantile),
              gmcpm = geomean(count))

  p <- ggplot(y, aes(x=qcpm,
                     y=reorder(factor(aoi), qcpm),
                     color=factor(bg),
                     shape=factor(keep_aoi))) +
    geom_point(size=1, alpha=0.5) +
    scale_color_manual(values = pals::cols25()) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(color = "BG", shape = "Keep AOI", x="Median CPM", y="AOI",
         title=paste0("BG vs Signal CPM Quantile = ", cpm_quantile))
  return(p)
}

#'
#' scatter plot signal vs bg
#'
#' @param x list dataset
#' @returns ggplot object
scatter_aoi_bg_vs_signal <- function(x, cpm_quantile = 0.75) {
  y <- x %>%
    group_by(aoi, bg) %>%
    summarise(keep_aoi = first(keep_aoi),
              slide = first(slide),
              qcpm = stats::quantile(count, cpm_quantile),
              gmcpm = geomean(count)) %>%
    pivot_wider(names_from="bg", values_from=c("gmcpm", "qcpm"))

  bg_gene_cor <- cor.test(y$qcpm_TRUE, y$qcpm_FALSE)$estimate
  label_cor <- sprintf("Correlation r = %.2f", bg_gene_cor)
  p <- ggplot(y, aes(x=qcpm_TRUE, y=qcpm_FALSE, color=keep_aoi)) +
    geom_point(alpha=0.6) +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color="black", linetype="dashed") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = pals::cols25()) +
    theme_minimal() +
    labs(color="Keep AOI",
         x="CPM background", y="CPM genes",
         title=paste0("BG vs Gene CPM quantile = ", cpm_quantile),
         subtitle=label_cor)
  return(p)
}


#'
#' qc plots to assess the effect of background subtraction
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param list dataset
#' @returns ggplot object
#' @export
st_geomx_plot_bgsub <- function(ds, cpm_quantile) {
  # sample annotations to use in plots
  s <- select(ds$samples, aoi, slide, keep)

  # compute cpm (before background subtract), bind to meta, join to annot
  x <- cpm(ds$counts)
  x <- bind_cols(ds$meta, x)
  x <- x %>%
    pivot_longer(colnames(ds$counts), names_to="aoi", values_to="count") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))
  dotplot1 <- dotplot_aoi_bg_vs_signal(x)
  scatter1 <- scatter_aoi_bg_vs_signal(x, cpm_quantile)

  # background subtraction, bind meta, join sample annotation
  x <- background_subtract(ds$counts, ds$samples$bg_geomean)
  x <- cpm(x)
  x <- bind_cols(ds$meta, x)
  x <- x %>%
    pivot_longer(colnames(ds$counts), names_to="aoi", values_to="count") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))
  dotplot2 <- dotplot_aoi_bg_vs_signal(x)
  scatter2 <- scatter_aoi_bg_vs_signal(x, cpm_quantile)

  p <- (dotplot1 | scatter1) / (dotplot2 | scatter2)
  return(p)
}


#'
#' qc plots to assess the effect of background subtraction
#'
#' @import ggplot2
#'
#' @param list dataset
#' @param aes_y string
#' @param aes_fill string
#' @returns ggplot object
#' @export
st_geomx_plot_count_distribution <- function(ds, aes_y, aes_fill, ridge_scale=5) {
  s <- select(ds$samples, aoi, slide, keep, {{ aes_y }}, {{ aes_fill }})
  x <- bind_cols(ds$meta, cpm(ds$counts))
  x <- x %>%
    pivot_longer(colnames(ds$counts), names_to="aoi", values_to="count") %>%
    inner_join(s, by=c("aoi"="aoi"), suffix=c("_gene", "_aoi"))

  # ridge plot showing variability in background noise
  p <- ggplot(x, aes(x=count, y=reorder({{ aes_y }}, count), fill={{ aes_fill }})) +
    stat_density_ridges(color="#00000066", alpha=0.4, scale=ridge_scale, rel_min_height=0.001, quantile_lines=TRUE) +
    scale_x_log10() +
    scale_fill_manual(values = pals::cols25()) +
    xlab("Count") +
    ylab("AOIs") +
    labs(fill = "Slide") +
    theme_ridges() +
    theme(axis.text.y = element_blank())
    #facet_grid(slide ~ ., scales="free_y")
  return(p)
}

