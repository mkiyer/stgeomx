#
# Spatial Transcriptomics using NanoString GeoMx DSP
#

# get_de_ranks <- function(de, a) {
#   ranks <- de %>%
#     filter(analysis == a) %>%
#     select(gene, log2fc, padj) %>%
#     mutate(rank = log2fc)
#   ranks = sort(setNames(ranks$rank, ranks$gene), decreasing = TRUE)
#   return(ranks)
# }

# run_batch_fgsea <- function(my_analyses, my_de, my_gs, my_prefix, my_plot_dir, my_padj_cutoff=0.05) {
#   gsea <- NULL
#   for (a in my_analyses) {
#     print(a)
#     ranks <- get_de_ranks(my_de, a)
#     res <- fgsea(pathways = my_gs, stats = ranks, minSize = 10, eps = 0, nPermSimple = 10000, nproc=1)
#     res <- res %>%
#       replace_na(list(pval=1, padj=1, log2err=0, NES=1)) %>%
#       mutate(analysis = a)
#
#     gsea <- bind_rows(gsea, res)
#     # for (i in 1:nrow(res)) {
#     #   x <- res[i,]
#     #   if (x$padj >= my_padj_cutoff) { next }
#     #   print(x$pathway)
#     #   p <- plot_gsea_enrichment(x, ranks, my_gs)
#     #   f <- file.path(my_plot_dir, paste0(my_prefix, "_", a, "_gs_", x$pathway, ".pdf"))
#     #   ggsave(f, plot=p, device="pdf", width=5, height=3)
#     # }
#   }
#   return(gsea)
# }
#
#
# get_de_ranks <- function(de, a, method=c("lfc", "p", "lfcp", "q")) {
#   x <- de %>%
#     filter(analysis == a) %>%
#     mutate(
#       lfc = log2fc,
#       p = sign(log2fc) * -log10(pval),
#       lfcp = log2fc * -log10(pval),
#       q = percent_rank(lfcp) - 0.5
#     )
#   ranks <- sort(setNames(pull(x, method), x$gene), decreasing = TRUE)
#   return(ranks)
# }
#
# run_batch_fgsea <- function(gs, de, analyses, rank_method,
#                             fgsea.nperm=10000) {
#   gsea <- NULL
#   for (a in analyses) {
#     print(a)
#     ranks <- get_de_ranks(de, a, rank_method)
# #     res <- fgseaSimple(pathways=gs,
# #                        stats=ranks,
# #                        nperm=fgsea.nperm)
# # #                       nproc=1)
#     res <- fgsea(pathways=gs,
#                  stats=ranks,
#                  eps=0,
#                  nPermSimple=fgsea.nperm)
# #                 nproc=1)
#     res <- res %>%
#       replace_na(list(pval=1, padj=1, log2err=0, NES=1)) %>%
#       mutate(analysis = a)
#     gsea <- bind_rows(gsea, res)
#   }
#   return(gsea)
# }
#
# plot_gsea_enrichment <- function(x, ranks, gs) {
#   txt_stat <- paste0("ES=", round(x$ES,2), " NES=", round(x$NES, 2), " padj=", format(x$padj, scientific=TRUE, digits=3))
#   txt_title <- paste0("ranks: ", a, " gs: ", x$pathway)
#   p <- plotEnrichment(gs[[x$pathway]], ranks) +
#     annotate(geom="text", x=150, y=0.1, label=txt_stat, hjust=0) +
#     labs(title = txt_title,
#          xlab = "Rank",
#          ylab = "Enrichment Score") +
#     theme(plot.title = element_text(size=8))
#   return(p)
# }
#
#
#
# gmt_to_tbl <- function(gmt) {
#   gs_to_tbl <- function(gs_name, gs_genes) {
#     bind_cols(gs_name = rep(gs_name, length(gs_genes)),
#               gene_symbol = gs_genes)
#   }
#   tbl <- NULL
#   for (i in 1:length(gmt)) {
#     tbl <- bind_rows(tbl, gs_to_tbl(names(gmt)[i], gmt[[i]]))
#   }
#   return(tbl)
# }
#
#
# plot_gsea_volcano <- function(x, padj_cutoff = 0.01, max.overlaps = Inf) {
#   p <- ggplot(x, aes(x=ES, y=-log10(pval), color=analysis)) +
#     geom_point(data=subset(x, padj > padj_cutoff), size=1, alpha=0.4) +
#     geom_point(data=subset(x, padj <= padj_cutoff), size=2, alpha=0.8) +
#     geom_text_repel(data=subset(x, padj <= padj_cutoff), color="black", size=3, aes(label=pathway), max.overlaps=max.overlaps) +
#     #ylab("-log10(adjusted p-value)") +
#     #xlab("NES") +
#     theme_bw() +
#     theme(legend.position="bottom") +
#     theme(axis.line = element_line(color = "black")) +
#     coord_flip()
#   return(p)
# }
#
# plot_gsea_barplot <- function(x, padj_cutoff = 0.01) {
#   x <- mutate(x, sig = ifelse(padj < padj_cutoff, ifelse(NES < 0, "dn", "up"), "no"),
#               neglog10padj = ifelse(padj < padj_cutoff, -log10(padj), NA))
#   p <- ggplot(x, aes(x=reorder(pathway, NES), y=NES, fill=neglog10padj)) +
#     geom_col() +
#     scale_fill_viridis_c(limits=c(0, 5)) +
#     #scale_fill_gradient(low="blue", high="cyan", limits=c(-5, 5), na.value="#aaaaaa") +
#     coord_flip() +
#     labs(x="Analysis", y="Normalized Enrichment Score") +
#     theme_bw()
#   return(p)
# }
#
# gsea_leading_edge_matrix <- function(x, my_padj_cutoff) {
#   # filter for significant results
#   sig_pathways <- x %>%
#     filter(padj < my_padj_cutoff) %>% select(pathway) %>% distinct() %>% pull()
#   x <- x %>%
#     filter(pathway %in% sig_pathways) %>%
#     arrange(desc(NES))
#
#   # get leading edge genes
#   leading_edge_genes <- summarise(x, result = unlist(leadingEdge, use.names=FALSE)) %>% distinct() %>% pull()
#
#   # matrix genes (rows) vs pathways (columns)
#   m <- matrix(data = 0, nrow=length(leading_edge_genes), ncol=nrow(x))
#   m <- as.data.frame(m)
#   rownames(m) <- leading_edge_genes
#   colnames(m) <- x$pathway
#   for (i in 1:nrow(x)) {
#     p <- x$pathway[i]
#     nes <- x$NES[i]
#     padj <- x$padj[i]
#     for (j in x$leadingEdge[i]) {
#       m[j, p] <- sign(nes)
#     }
#   }
#   # filter genes only found in a single pathway
#   m <- m[abs(rowSums(m)) >= 2, ]
#   # filter empty columns
#   m <- m[, abs(colSums(m)) > 1]
#   # order by number of pathways sharing the gene
#   m <- m[order(rowSums(m), decreasing=TRUE),]
#   m <- m[,order(colSums(m), decreasing=TRUE)]
#   return(m)
# }

