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
                                     min_snauc=0.5,
                                     min_frac_expr=0) {
  keep_aois <- table(ds$samples$keep)
  keep_pct <- 100*keep_aois["TRUE"] / nrow(ds$samples)
  subtitle <- sprintf("Retained %d/%d (%.2f%%) AOIs", keep_aois["TRUE"], nrow(ds$samples), keep_pct)

  p1 <- ggplot(ds$samples, aes(x=num_counts, y=frac_expr, color=keep)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_frac_expr, linetype = "dashed", color = "red") +
#    scale_color_viridis_c() +
    scale_x_log10() +
    labs(x="Counts", y="Frac detectable genes",
         title="AOI QC Metrics",
         subtitle=subtitle) +
    theme_minimal()

  p2 <- ggplot(ds$samples, aes(x=num_counts, y=bg_auc, color=keep)) +
    geom_point(alpha=0.6) +
    geom_vline(xintercept = min_counts, linetype = "dashed", color = "red") +
    geom_hline(yintercept = min_snauc, linetype = "dashed", color = "red") +
#    scale_color_viridis_c() +
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
#' @param list dataset
#' @param min_frac_expr number threshold range (0-1.0)
#' @returns ggplot object
#' @export
st_geomx_plot_gene_filter <- function(ds, min_frac_expr=0) {
  # apply cutoff
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


