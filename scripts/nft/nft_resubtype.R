library(parallel)
library(SingleCellExperiment)
library(scFlow)
library(scFlowData)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
library(gplots)
library(viridis)
library(cowplot)

ensembl_fp <- system.file("extdata","ensembl_mappings.tsv",package="scFlowData")
ctd_fp <- system.file("extdata","ctd",package="scFlowData")
sce <- read_sce("~/SCADSplice/scrnaseq/nft/results/typed_sce", read_metadata = FALSE)

#__________________________________________________
# EN neruo

celltype <- c("EN-L2", "EN-L2-3", "EN-L3-4", "EN-L4-5", "EN-L3-5",
              "EN-L4-6", "EN-L5-6","EN-L6")
EN_sce <- sce[, sce$cluster_celltype == celltype]

EN_sce <- reduce_dims_sce(EN_sce,
                          pca_dims = 10)

set.seed(321)
EN_sce <- cluster_sce(EN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.01,
                      k = 50
)

plot_reduced_dim(EN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


set.seed(321)
EN_sce <- map_celltypes_sce(
  EN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

plot_reduced_dim(EN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                 highlight_feature = NA, label_clusters = F, size = 2,
                 alpha = 1)

EN_sce <- annotate_celltype_metrics(
  EN_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest","diagnosis","sex", "age", "PMI", "RIN"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(EN_sce)

write_celltype_mappings(EN_sce, folder_path = getwd())
write_sce(
  sce,
  "~/SCADSplice/scrnaseq/nft/results/EN_typed_sce",
  write_metadata = T
)

write_celltype_mappings(EN_sce, folder_path = "~")

png("~/SCADSplice/scrnaseq/nft/results/plots/heatmap_clusters.png",
    width = 160, height = 330, units = "mm", res = 600)


png("~/SCADSplice/scrnaseq/nft/results/EN_plots/umap_celltype.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                 highlight_feature = NA, label_clusters = F, size = 0.8,
                 alpha = 1) + scale_colour_manual(values = brewer.pal(n = 7, name = "Paired")) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()


png("~/SCADSplice/scrnaseq/nft/results/EN_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.8) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c(brewer.pal(n = 12, name = "Set3")))
dev.off()


png("~/SCADSplice/scrnaseq/nft/results/EN_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

markers <- c("SNAP25","SLC17A7", "GAD1", "NRGN", "CHN1", "SYT1",
             "SST", "LAMP5", "PVALB", "VIP", "RORB", "COL5A2",
             "THEMIS", "GABRG1", "VGF", "ENC1", "NTNG2", "CXCL14",
             "KIT","PTGDS")

plotlist <- list()
for (i in markers) {
  plotlist[[i]] <- plot_reduced_dim_gene(
    EN_sce,
    reduced_dim = "UMAP_Liger",
    gene = i,
    size = 0.1,
    alpha = 10,
    palette = c("grey80", "#440154FF")
  )
}

png("~/SCADSplice/scrnaseq/nft/results/EN_plots/umap_genes2.png",
    width = 170, height = 212, units = "mm", res = 600)
plot_grid(ncol = 4, plotlist = plotlist)
dev.off()


dt <- EN_sce@metadata[["markers"]][["cluster_celltype"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 7, nrow = 31)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nft/results/EN_plots/heatmap_celltype.png",
    width = 140, height = 330, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

cmp <- read.delim("~/SCADSplice/scrnaseq/nft/results/celltype_mappings.tsv", header = T)
tot <- data.frame(cmp$cluster_celltype, cmp$clusters)
tot <- tot[order(tot$cmp.cluster_celltype),]
j <- 1
v <- "z"
for (i in tot$cmp.cluster_celltype) {
  tot$col[which(tot$cmp.cluster_celltype %in% i)] <-
    rainbow(15)[j]
  if (i != v) {
    j <- j+1
  }
  v = i
}
tot <- tot[order(tot$cmp.clusters),]
dt <- EN_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 11, nrow = 44)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1),
          col_side_colors = tot$cmp.cluster_celltype)
ggheatmap(mat,
          ColSideColors = tot$col)

png("~/SCADSplice/scrnaseq/nft/results/EN_plots/heatmap_clusters.png",
    width = 160, height = 330, units = "mm", res = 600)
ggheatmap(log2(mat+1),
          col_side_colors = tot$cmp.cluster_celltype)
dev.off()

png("~/SCADSplice/scrnaseq/nft/results/EN_plots/heatmap_clusters2.png",
    width = 170, height = 360, units = "mm", res = 600)
heatmap.2(log2(mat+1), col=viridis(1024), ColSideColors = tot$col,
          trace = "none", density.info = "none")
dev.off()


#__________________________________________________
# IN neruo

celltype <- c("IN-VIP", "IN-SST", "IN-LAMP5", "IN-PVALB")
IN_sce <- sce[, sce$cluster_celltype == celltype]

IN_sce <- reduce_dims_sce(IN_sce,
                          pca_dims = 10)

set.seed(321)
IN_sce <- cluster_sce(IN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.01,
                      k = 50
)

plot_reduced_dim(IN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


set.seed(321)
IN_sce <- map_celltypes_sce(
  IN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

plot_reduced_dim(IN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                 highlight_feature = NA, label_clusters = F, size = 2,
                 alpha = 1)

IN_sce <- annotate_celltype_metrics(
  IN_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest","diagnosis","sex", "age", "PMI", "RIN"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(IN_sce)


write_sce(
  sce,
  "~/SCADSplice/scrnaseq/nft/results/IN_typed_sce",
  write_metadata = T
)

write_celltype_mappings(IN_sce, folder_path = "~")

png("~/SCADSplice/scrnaseq/nft/results/IN_plots/umap_celltype.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                 highlight_feature = NA, label_clusters = F, size = 0.8,
                 alpha = 1) + scale_colour_manual(values = brewer.pal(n = 7, name = "Paired")) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()


png("~/SCADSplice/scrnaseq/nft/results/IN_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.8) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c(brewer.pal(n = 12, name = "Set3")))
dev.off()


png("~/SCADSplice/scrnaseq/nft/results/IN_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

markers <- c("SNAP25","SLC17A7", "GAD1", "NRGN", "CHN1", "SYT1",
             "SST", "LAMP5", "PVALB", "VIP", "RORB", "COL5A2",
             "THEMIS", "GABRG1", "VGF", "ENC1", "NTNG2", "CXCL14",
             "KIT","PTGDS")

plotlist <- list()
for (i in markers) {
  plotlist[[i]] <- plot_reduced_dim_gene(
    IN_sce,
    reduced_dim = "UMAP_Liger",
    gene = i,
    size = 0.1,
    alpha = 10,
    palette = c("grey80", "#440154FF")
  )
}

png("~/SCADSplice/scrnaseq/nft/results/IN_plots/umap_genes2.png",
    width = 170, height = 212, units = "mm", res = 600)
plot_grid(ncol = 4, plotlist = plotlist)
dev.off()


dt <- IN_sce@metadata[["markers"]][["cluster_celltype"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 4, nrow = 20)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nft/results/IN_plots/heatmap_celltype.png",
    width = 140, height = 330, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

cmp <- read.delim("~/SCADSplice/scrnaseq/nft/results/celltype_mappings.tsv", header = T)
tot <- data.frame(cmp$cluster_celltype, cmp$clusters)
tot <- tot[order(tot$cmp.cluster_celltype),]
j <- 1
v <- "z"
for (i in tot$cmp.cluster_celltype) {
  tot$col[which(tot$cmp.cluster_celltype %in% i)] <-
    rainbow(15)[j]
  if (i != v) {
    j <- j+1
  }
  v = i
}
tot <- tot[order(tot$cmp.clusters),]
dt <- IN_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 5, nrow = 24)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1),
          col_side_colors = tot$cmp.cluster_celltype)
ggheatmap(mat,
          ColSideColors = tot$col)

png("~/SCADSplice/scrnaseq/nft/results/IN_plots/heatmap_clusters.png",
    width = 160, height = 330, units = "mm", res = 600)
ggheatmap(log2(mat+1),
          col_side_colors = tot$cmp.cluster_celltype)
dev.off()

png("~/SCADSplice/scrnaseq/nft/results/IN_plots/heatmap_clusters2.png",
    width = 170, height = 360, units = "mm", res = 600)
heatmap.2(log2(mat+1), col=viridis(1024), ColSideColors = tot$col,
          trace = "none", density.info = "none")
dev.off()

write_sce(
  IN_sce,
  "~/SCADSplice/scrnaseq/nft/results/IN_typed_IN_sce",
  write_metadata = TRUE
)
