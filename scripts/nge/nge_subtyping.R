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

# allocate 2 cores for annotate_celltype_metrics()
trace(annotate_celltype_metrics, edit = T)
sce <- read_sce("~/SCADSplice/scrnaseq/nge/results/typed_sce")

##__________________________________________________
# IN neruo

IN_sce <- sce[, sce$cluster_celltype == "IN"]

IN_sce <- reduce_dims_sce(IN_sce,
                          input_reduced_dim = "Liger",
                          reduction_methods = "UMAP",
                          pca_dims = 10)

plot_reduced_dim(IN_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(IN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(IN_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(IN_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
IN_sce <- cluster_sce(IN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(IN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  IN_sce,
  reduced_dim = "UMAP_Liger",
  gene = "SLC17A7",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
IN_sce <- map_celltypes_sce(
  IN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(IN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.2,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_INtypes.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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


celltype <- c("IN-VIP", "IN-SST", "IN-LAMP5", "IN-PVALB")
IN_sce <- IN_sce[, IN_sce$cluster_celltype == celltype]

set.seed(321)
IN_sce <- cluster_sce(IN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
IN_sce <- map_celltypes_sce(
  IN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(IN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_INtypes_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggsci::default_locuszoom") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

IN_sce <- annotate_celltype_metrics(
  IN_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(IN_sce)

dt <- IN_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 5, nrow = 24)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/heatmap_clusters.png",
    width = 70, height = 130, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/IN_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(IN_sce, feature_dim = "RIN", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_celltype_mappings(IN_sce, folder_path = getwd())
write_sce(
  IN_sce,
  "~/SCADSplice/scrnaseq/nge/results/IN_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  IN_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/model_CTR.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  IN_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/IN_plots/model_ADM.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- IN_sce@colData@listData[["barcode"]]
type <- IN_sce@colData@listData[["cluster_celltype"]]
subtype <- IN_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/IN_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#__________________________________________________
# Endo

Endo_sce <- sce[, sce$cluster_celltype == "Endo"]

Endo_sce <- reduce_dims_sce(Endo_sce,
                          input_reduced_dim = "Liger",
                          reduction_methods = "UMAP",
                          pca_dims = 10)

plot_reduced_dim(Endo_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(Endo_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(Endo_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(Endo_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
Endo_sce <- cluster_sce(Endo_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(Endo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  Endo_sce,
  reduced_dim = "UMAP_Liger",
  gene = "SLC17A7",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
Endo_sce <- map_celltypes_sce(
  Endo_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Endo_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Endo_sce <- map_custom_celltypes(
  Endo_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Endo_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 0.2,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_INtypes.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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


celltype <- c("IN-VIP", "IN-SST", "IN-LAMP5", "IN-PVALB")
Endo_sce <- Endo_sce[, Endo_sce$cluster_celltype == celltype]

set.seed(321)
Endo_sce <- cluster_sce(Endo_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
Endo_sce <- map_celltypes_sce(
  Endo_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(Endo_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_INtypes_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggsci::default_locuszoom") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

Endo_sce <- annotate_celltype_metrics(
  Endo_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(Endo_sce)

dt <- Endo_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 3, nrow = 15)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/heatmap_clusters.png",
    width = 60, height = 80, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Endo_sce, feature_dim = "RIN", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  Endo_sce,
  "~/SCADSplice/scrnaseq/nge/results/Endo_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  Endo_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/model_CTR.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  Endo_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Endo_plots/model_ADM.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- Endo_sce@colData@listData[["barcode"]]
type <- Endo_sce@colData@listData[["cluster_celltype"]]
subtype <- Endo_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/Endo_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

##__________________________________________________
# Astro

Astro_sce <- sce[, sce$cluster_celltype == "Astro"]

Astro_sce <- reduce_dims_sce(Astro_sce,
                          input_reduced_dim = "Liger",
                          reduction_methods = "UMAP",
                          pca_dims = 10)

plot_reduced_dim(Astro_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(Astro_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(Astro_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(Astro_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
Astro_sce <- cluster_sce(Astro_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(Astro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  Astro_sce,
  reduced_dim = "UMAP_Liger",
  gene = "SLC17A7",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
Astro_sce <- map_celltypes_sce(
  Astro_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Astro_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Astro_sce <- map_custom_celltypes(
  Astro_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Astro_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_INtypes.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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




set.seed(321)
Astro_sce <- cluster_sce(Astro_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
Astro_sce <- map_celltypes_sce(
  Astro_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(Astro_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_types_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggsci::default_locuszoom") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

Astro_sce <- annotate_celltype_metrics(
  Astro_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(Astro_sce)

dt <- Astro_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 3, nrow = 15)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/heatmap_clusters.png",
    width = 60, height = 80, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Astro_sce, feature_dim = "RIN", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  Astro_sce,
  "~/SCADSplice/scrnaseq/nge/results/Astro_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  Astro_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/model_CTR.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  Astro_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Astro_plots/model_ADM.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- Astro_sce@colData@listData[["barcode"]]
type <- Astro_sce@colData@listData[["cluster_celltype"]]
subtype <- Astro_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/Astro_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#trace(annotate_celltype_metrics, edit = T)
sce <- read_sce("~/SCADSplice/scrnaseq/nge/results/sn_typed_sce")

trace(annotate_celltype_metrics, edit = T)
EN_sce <- read_sce("~/SCADSplice/scrnaseq/nge/results/sn_typed_sce")

#__________________________________________________
# EN

EN_sce <- EN_sce[, EN_sce$cluster_celltype == "EN"]

plot_reduced_dim(EN_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(EN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(EN_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(EN_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
EN_sce <- cluster_sce(EN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(EN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  EN_sce,
  reduced_dim = "UMAP_Liger",
  gene = "GAD2",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
EN_sce <- map_celltypes_sce(
  EN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

p <- plot_reduced_dim(EN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_INtypes.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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


EN_sce <- reduce_dims_sce(EN_sce,
                          input_reduced_dim = "Liger",
                          reduction_methods = "UMAP",
                          pca_dims = 10)

set.seed(321)
EN_sce <- cluster_sce(EN_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
EN_sce <- map_celltypes_sce(
  EN_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(EN_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

EN_sce <- map_custom_celltypes(
  EN_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(EN_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_types_AA.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggthemes::manyeys") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

EN_sce <- annotate_celltype_metrics(
  EN_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(EN_sce)

dt <- EN_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 19, nrow = 70)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/heatmap_clusters.png",
    width = 170, height = 300, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/EN_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(EN_sce, feature_dim = "RIN", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  EN_sce,
  "~/SCADSplice/scrnaseq/nge/results/EN_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  EN_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/model_CTR.png",
    width = 450, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  EN_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/EN_plots/model_ADM.png",
    width = 450, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- EN_sce@colData@listData[["barcode"]]
type <- EN_sce@colData@listData[["cluster_celltype"]]
subtype <- EN_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/EN_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##__________________________________________________
# Micro

Micro_sce <- sce[, sce$cluster_celltype == "Micro"]


plot_reduced_dim(Micro_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(Micro_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(Micro_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(Micro_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
Micro_sce <- cluster_sce(Micro_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(Micro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  Micro_sce,
  reduced_dim = "UMAP_Liger",
  gene = "GAD1",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
Micro_sce <- map_celltypes_sce(
  Micro_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Micro_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Micro_sce <- map_custom_celltypes(
  Micro_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Micro_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_types.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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

Micro_sce <- reduce_dims_sce(Micro_sce,
                             input_reduced_dim = "Liger",
                             reduction_methods = "UMAP",
                             pca_dims = 10)


set.seed(321)
Micro_sce <- cluster_sce(Micro_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
Micro_sce <- map_celltypes_sce(
  Micro_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Micro_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Micro_sce <- map_custom_celltypes(
  Micro_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Micro_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_types_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggsci::default_locuszoom") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

Micro_sce <- annotate_celltype_metrics(
  Micro_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(Micro_sce)

dt <- Micro_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 4, nrow = 13)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/heatmap_clusters.png",
    width = 60, height = 80, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Micro_sce, feature_dim = "APOE", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  Micro_sce,
  "~/SCADSplice/scrnaseq/nge/results/Micro_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  Micro_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/model_CTR.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  Micro_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Micro_plots/model_ADM.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- Micro_sce@colData@listData[["barcode"]]
type <- Micro_sce@colData@listData[["cluster_celltype"]]
subtype <- Micro_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/Micro_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

##__________________________________________________
# OPC

OPC_sce <- sce[, sce$cluster_celltype == "OPC"]


plot_reduced_dim(OPC_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(OPC_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(OPC_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(OPC_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
OPC_sce <- cluster_sce(OPC_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(OPC_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  OPC_sce,
  reduced_dim = "UMAP_Liger",
  gene = "GAD1",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
OPC_sce <- map_celltypes_sce(
  OPC_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(OPC_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

OPC_sce <- map_custom_celltypes(
  OPC_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(OPC_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_types.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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

OPC_sce <- reduce_dims_sce(OPC_sce,
                             input_reduced_dim = "Liger",
                             reduction_methods = "UMAP",
                             pca_dims = 10)


set.seed(321)
OPC_sce <- cluster_sce(OPC_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
OPC_sce <- map_celltypes_sce(
  OPC_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(OPC_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

OPC_sce <- map_custom_celltypes(
  OPC_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(OPC_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_types_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggsci::default_locuszoom") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

OPC_sce <- annotate_celltype_metrics(
  OPC_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(OPC_sce)

dt <- OPC_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 3, nrow = 15)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/heatmap_clusters.png",
    width = 60, height = 80, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(OPC_sce, feature_dim = "APOE", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  OPC_sce,
  "~/SCADSplice/scrnaseq/nge/results/OPC_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  OPC_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/model_CTR.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  OPC_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/OPC_plots/model_ADM.png",
    width = 170, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- OPC_sce@colData@listData[["barcode"]]
type <- OPC_sce@colData@listData[["cluster_celltype"]]
subtype <- OPC_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/OPC_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

##__________________________________________________
# Oligo

Oligo_sce <- sce[, sce$cluster_celltype == "Oligo"]


plot_reduced_dim(Oligo_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger", alpha = 2, size = 1)
plot_reduced_dim(Oligo_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger", alpha = 1, size = 1)
plot_reduced_dim(Oligo_sce, feature_dim = "age", reduced_dim = "UMAP_Liger", alpha = 2, size = 0.3)
plot_reduced_dim(Oligo_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",  alpha = 2, size = 0.3)

set.seed(321)
Oligo_sce <- cluster_sce(Oligo_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

plot_reduced_dim(Oligo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger", label_clusters = T,
                 alpha = 2, size = 1)


plot_reduced_dim_gene(
  Oligo_sce,
  reduced_dim = "UMAP_Liger",
  gene = "GAD1",
  size = 1,
  alpha = 1,
  palette = c("grey80", "#440154FF")
)

set.seed(321)
Oligo_sce <- map_celltypes_sce(
  Oligo_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Oligo_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Oligo_sce <- map_custom_celltypes(
  Oligo_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Oligo_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_types.png",
    width = 170, height = 130, units = "mm", res = 600)
p  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_clusters.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 0.2) +
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

Oligo_sce <- reduce_dims_sce(Oligo_sce,
                             input_reduced_dim = "Liger",
                             reduction_methods = "UMAP",
                             pca_dims = 10)


set.seed(321)
Oligo_sce <- cluster_sce(Oligo_sce,
                      cluster_method = "leiden",
                      reduction_method = "UMAP_Liger",
                      pca_dims = 10,
                      res = 0.001,
                      k = 50
)

set.seed(321)
Oligo_sce <- map_celltypes_sce(
  Oligo_sce,
  ctd_folder = ctd_fp,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

write_celltype_mappings(Oligo_sce, folder_path = "~")

celltype_mappings <- read_celltype_mappings("~/SCADSplice/scrnaseq/nge/results/celltype_mappings.tsv")

Oligo_sce <- map_custom_celltypes(
  Oligo_sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

p <- plot_reduced_dim(Oligo_sce, feature_dim = "cluster_celltype", reduced_dim = "UMAP_Liger",
                      highlight_feature = NA, label_clusters = F, size = 1,
                      alpha = 1)

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_types_A.png",
    width = 170, height = 130, units = "mm", res = 600)
p  + scale_color_paletteer_d("ggthemes::manyeys") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_clusters_A.png",
    width = 170, height = 170, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "clusters", reduced_dim = "UMAP_Liger",
                 label_clusters = T,
                 alpha = 2, size = 1) +
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

Oligo_sce <- annotate_celltype_metrics(
  Oligo_sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest", "diagnosis", "age", "PMI", "APOE"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(Oligo_sce)

dt <- Oligo_sce@metadata[["markers"]][["clusters"]][["marker_plot"]][["data"]]
dm <- setNames(data.frame(matrix(ncol = 13, nrow = 40)), c(as.character(unique(dt$Group))))
rownames(dm) <- unique(dt$Gene)
for (i in c(as.character(unique(dt$Group)))) {
  dm[i] <- dt[dt$Group == i,3]
}

mat <- as.matrix(dm)
ggheatmap(log2(mat+1))
ggheatmap(mat)

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/heatmap_clusters.png",
    width = 110, height = 190, units = "mm", res = 600)
ggheatmap(log2(mat+1))
dev.off()

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_manifest.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "manifest", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_diagnosis.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "diagnosis", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FF0000", "#FFC400","#00BBFF"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_age.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "age", reduced_dim = "UMAP_Liger",
                 alpha = 0.5, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_sex.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "sex", reduced_dim = "UMAP_Liger",
                 alpha = 1, size = 1) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values = c("#FFAFC7", "#73D7EE"))
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_PMI.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "PMI", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()
png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/umap_RIN.png",
    width = 170, height = 130, units = "mm", res = 600)
plot_reduced_dim(Oligo_sce, feature_dim = "APOE", reduced_dim = "UMAP_Liger",
                 alpha = 2, size = 0.5) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(axis.line = element_line(color="black", size = 0.2),
        line = element_line(colour = "black"),
        text = element_text(color = "black"),
        title = element_text(color = "black")) +
  theme(legend.title=element_blank())
dev.off()

write_sce(
  Oligo_sce,
  "~/SCADSplice/scrnaseq/nge/results/Oligo_typed_sce",
  write_metadata = T
)

results <- model_celltype_freqs(
  Oligo_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "CTR",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/model_CTR.png",
    width = 390, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

results <- model_celltype_freqs(
  Oligo_sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "ADM",
  var_order = c("CTR", "ADM", "ADH")
)

report_celltype_model(results)

png("~/SCADSplice/scrnaseq/nge/results/Oligo_plots/model_ADM.png",
    width = 390, height = 100, units = "mm", res = 600)
results[["dirichlet_plot"]]
dev.off()

barcode <- Oligo_sce@colData@listData[["barcode"]]
type <- Oligo_sce@colData@listData[["cluster_celltype"]]
subtype <- Oligo_sce@colData@listData[["clusters"]]
celanot <- data.frame(barcode, type, subtype)

write.table(celanot, "~/SCADSplice/scrnaseq/nge/metadata/Oligo_cells_annotation.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
