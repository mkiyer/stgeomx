#
# Spatial Transcriptomics using NanoString GeoMx DSP
#


# de_format_results <- function(res, analysis, method,
#                               padj_cutoff, log2fc_cutoff) {
#   mutate(res,
#          de = case_when(padj > padj_cutoff ~ "no",
#                         log2fc < -log2fc_cutoff ~ "dn",
#                         log2fc > log2fc_cutoff ~ "up",
#                         TRUE ~ "no"),
#          analysis = analysis,
#          method = method)
# }

# run_limma_trend <- function(y, design, contrasts) {
#   fit <- lmFit(y, design)
#   fit <- contrasts.fit(fit, contrasts)
#   fit <- eBayes(fit, trend=TRUE)
#   return(fit)
# }

# run_limma_trend <- function(y, design, coef, padj_cutoff, log2fc_cutoff,
#                             analysis, method) {
#   fit <- lmFit(y, design)
#   fit <- eBayes(fit, trend=TRUE)
#   res <- topTable(fit, coef=coef, number=Inf, sort.by="none")
#   res <- limma_rename_result_cols(res)
#   res <- de_format_results(res, analysis=analysis, method=method,
#                            padj_cutoff = padj_cutoff, log2fc_cutoff = log2fc_cutoff)
#   return(res)
# }

# process_limma <- function(fit, contrasts, method,
#                           padj_cutoff, log2fc_cutoff) {
#   x <- NULL
#   for (coef in colnames(contrasts)) {
#     res <- topTable(fit, coef=coef, number=Inf, sort.by="none")
#     res <- limma_rename_result_cols(res)
#     res <- de_format_results(res,
#                              analysis=coef,
#                              method=method,
#                              padj_cutoff=padj_cutoff,
#                              log2fc_cutoff=log2fc_cutoff)
#     x <- bind_rows(x, res)
#   }
#   return(x)
# }


run_limma_trend <- function(y, design, contrasts, padj_cutoff, log2fc_cutoff) {
  fit <- lmFit(y, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit, trend=TRUE)

  x <- NULL
  for (coef in colnames(contrasts)) {
    res <- topTable(fit, coef=coef, number=Inf, sort.by="none")
    res <- as_tibble(res, rownames="gene")
    res <- res %>%
      select(gene,
             log2fc=logFC,
             avgexpr=AveExpr,
             pval=P.Value,
             padj=adj.P.Val) %>%
      mutate(de = case_when(padj > padj_cutoff ~ "no",
                            log2fc < -log2fc_cutoff ~ "dn",
                            log2fc > log2fc_cutoff ~ "up",
                            TRUE ~ "no"),
             contrast = coef)
    x <- bind_rows(x, res)
  }
  return(x)
}



plot_de_volcano <- function(res) {
  # res <- res %>%
  #   mutate(de = case_when(padj > 0.05 ~ "no",
  #                         log2fc < 0 ~ "dn",
  #                         log2fc > 0 ~ "up"))
  p <- ggplot(res, aes(x=log2fc, y=-log10(padj), color=de, size=de)) +
    geom_point(alpha=0.7) +
    #    geom_text_repel(data=subset(res, de != "no"), color="black", size=3, aes(label=gene), max.overlaps=Inf) +
    geom_text_repel(data=subset(res, de != "no"), color="black", size=3, aes(label=gene)) +
    scale_color_manual(values=c("no"="grey", "dn"="blue", "up"="red")) +
    scale_size_manual(values=c("no"=1, "dn"=2, "up"=2)) +
    theme_minimal() +
    theme(legend.position="bottom") +
    theme(axis.line = element_line(color = "black"))
  labs(x="log2fc", y="-log10(padj)")
  return(p)
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
