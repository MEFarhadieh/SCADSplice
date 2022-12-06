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
sce <- read_sce("~/SCADSplice/scrnaseq/nge/results/sce_merged", read_metadata = T)

# Integration with LIGER

sce <- integrate_sce(
  sce,
  method             = "Liger",
  unique_id_var      = "manifest",
  take_gene_union    = FALSE,
  remove_missing     = TRUE,
  num_genes          = 3000,
  combine            = "union",
  keep_unique        = FALSE,
  capitalize         = FALSE,
  use_cols           = TRUE,
  k                  = 30,
  lambda             = 5.0,
  thresh             = 0.0001,
  max_iters          = 100,
  nrep               = 1,
  rand_seed          = 1,
  knn_k              = 20,
  k2                 = 500,
  prune_thresh       = 0.2,
  min_cells          = 2,
  quantiles          = 50,
  nstart             = 10,
  resolution         = 1,
  dims_use           = NULL,
  dist_use           = "CR",
  center             = FALSE,
  small_clust_thresh = 0
)

#  Save SCE 

write_sce(
  sce,
  file.path("~/SCADSplice/scrnaseq/nge/results/sce_integrated"),
  write_metadata = TRUE
)

sce@colData@listData[["total_features_by_counts"]] <-
  as.double(sce@colData@listData[["total_features_by_counts"]])

# reducu dimention 

sce <- reduce_dims_sce(sce,
                       pca_dims = 20)

# visualize clusters 

png("~/SCADSplice/scrnaseq/nge/results/plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 0.3, size = 0.1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "RIN", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

set.seed(321)
sce <- cluster_sce(sce,
                   cluster_method = "leiden",
                   reduction_method = "UMAP_Liger",
                   pca_dims = 20,
                   res = 0.001,
                   k = 100
)

write_sce(
  sce,
  "~/SCADSplice/scrnaseq/nge/results/sce_clustered"
)

sce <- annotate_integrated_sce(
  sce,
  categorical_covariates = c("manifest","diagnosis", "age", "sex", "PMI", "APOE"),
  input_reduced_dim = "UMAP"
)

report_integrated_sce(sce)



png("~/SCADSplice/scrnaseq/nge/results/plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.8) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c(brewer.pal(n = 12, name = "Set3"),
                                 brewer.pal(n = 12, name = "Set3"),
                                 brewer.pal(n = 12, name = "Set3")))
dev.off()

# automatic celltyping

set.seed(123)
sce <- map_celltypes_sce(
  sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.8,
                      alpha = 1)
p <- p + scale_colour_manual(values = rainbow(15)) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())

png("~/SCADSplice/scrnaseq/nge/results/plots/umap_celltype.png",
    width = 170, height = 130, units = "mm", res = 600)
p
dev.off()

# visualize celltypes

p <- plot_reduced_dim(sce, feature_dim = "cluster_celltype2", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.8,
                      alpha = 1)
lightp <- c("#9bd0b7", "#67CFDF", "#BD0000", "#d5b6d5", "#f6b4bf", "#f6aa90", "#f8df81")
p <- p + scale_colour_manual(values = lightp) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())

p

manifest <- sce@colData@listData[["manifest"]]
df <- data.frame(manifest)
df$UMAP1 = SingleCellExperiment::reducedDim(sce, "UMAP3D_Liger")[, 1]
df$UMAP2 = SingleCellExperiment::reducedDim(sce, "UMAP3D_Liger")[, 2]
df$UMAP3 = SingleCellExperiment::reducedDim(sce, "UMAP3D_Liger")[, 3]
df$cluster = sce@colData@listData[["clusters"]]
df$celltype = sce@colData@listData[["cluster_celltype"]]

plotly::plot_ly(df,
                x = df$UMAP1,
                y = df$UMAP2,
                z = df$UMAP3,
                type = "scatter3d",
                mode = "markers",
                color = df$cluster,
                size = 0.1

)

sce <- annotate_celltype_metrics(
  sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest","diagnosis", "age", "sex", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(sce)

write_celltype_mappings(sce, folder_path = "~/SCADSplice/scrnaseq/nge/metadata/")

sce@metadata$markers$cluster_celltype
sce@metadata$markers$clusters$marker_plot

# assment markers per cluster

markers <- c("AQP4", "GJA1", "ALDH1L1", "SLC1A3", "PDGFRA", "VCAN",
             "SOX10", "PLP1", "MBP", "MOG", "SNAP25", "NRGN", "SLC17A7",
             "GAD1", "C1QB", "TYROBP", "P2RY12", "CLDN5", "NOSTRIN", "PECAM1")

plotlist <- list()
for (i in markers) {
  plotlist[[i]] <- plot_reduced_dim_gene(
    sce,
    reduced_dim = "UMAP_Liger",
    gene = i,
    size = 0.1,
    alpha = 10,
    palette = c("grey80", "#440154FF")
  )
}

png("~/SCADSplice/scrnaseq/nge/results/plots/umap_genes2.png",
    width = 170, height = 212, units = "mm", res = 600)
plot_grid(ncol = 4, plotlist = plotlist)
dev.off()

# custom celltyping

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/metadata/celltype_mappings.tsv")
sce <- map_custom_celltypes(
  sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(sce, feature_dim = "cluster_celltype2", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.8,
                      alpha = 1)
lightp <- c("#9bd0b7", "#67CFDF", "#BD0000", "#d5b6d5", "#f6b4bf", "#f6aa90", "#f8df81")
p <- p + scale_colour_manual(values = lightp) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())

png("~/SCADSplice/scrnaseq/nge/results/plots/umap_celltype_modify.png",
    width = 170, height = 130, units = "mm", res = 600)
p
dev.off()

# heatmaps

dt <- sce@metadata[["markers"]][["cluster_celltype"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 7, nrow = 35)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/plots/heatmap_celltype.png",
    width = 100, height = 240, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

cmp <- read.delim("~/celltype_mappings.tsv", header = T)
tot <- data.frame(cmp$cluster_celltype, cmp$clusters)
tot <- tot[order(tot$cmp.cluster_celltype),]
j <- 1
v <- "z"
for (i in tot$cmp.cluster_celltype) {
  tot$col[which(tot$cmp.cluster_celltype %in% i)] <-
    lightp[j]
  if (i != v) {
    j <- j+1
  }
  v = i
}
tot <- tot[order(tot$cmp.clusters),]
rownames(tot) <- tot$cmp.clusters
dt <- sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 33, nrow = 117)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1),
          col_side_colors = tot$cmp.cluster_celltype)
ggheatmap(mat,
          ColSideColors = tot$col)

png("~/SCADSplice/scrnaseq/nge/results/plots/heatmap_clusters.png",
    width = 220, height = 420, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/plots/heatmap_clusters2.png",
    width = 170, height = 360, units = "mm", res = 600)
heatmap.2(log2(mat+1), col=viridis(1024), ColSideColors = tot$col,
          trace = "none", density.info = "none")
dev.off()

# Save SCE

write_sce(
  sce,
  "~/SCADSplice/scrnaseq/nge/results/typed_sce",
  write_metadata = TRUE
)

# build annotation file for Spliz and ReadZs

barcode <- sce@colData@listData[["barcode"]]
type <- sce@colData@listData[["cluster_celltype"]]
subtype <- sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/nge_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
