### Following codes were used to achive DGE

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(DTUrtle)
library(paletteer)

### DGE in Astrocyte
# import counts files

list.files("~/SCADSplice/rnaseq/results/astro/results/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/astro/results/salmon/*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

txi <- tximport(files, 
                type="salmon", 
                tx2gene=tx2gene[,c("transcript_id", "gene_name")], 
                countsFromAbundance = "lengthScaledTPM")


sampletype <- factor(c(rep("AAD",4), rep("ACt", 4)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

all(colnames(txi$counts) == rownames(metadata))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta, 
                                design = ~ sampletype)

# Normalize counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="~/SCADSplice/rnaseq/results/dge/astro_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, 
         annotation = meta)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Specify contrast for comparison of interest
contrast <- c("sampletype", "AAD", "ACt")

# Output results of Wald test for contrast
res <- DESeq2::results(dds, 
               contrast = contrast, 
               alpha = 0.05)

# Shrink the log2 fold changes to be more accurate
res <- lfcShrink(dds, 
                 coef = "sampletype_ACt_vs_AAD", 
                 type = "apeglm")

# Set thresholds
padj.cutoff <- 0.05

# Turn the results object into a tibble 
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res <- filter(res_tbl, 
                  padj < padj.cutoff)

write.table(sig_res, file="~/astro_DGE.tsv", sep="\t", quote=F, col.names=NA)

# Visualizing

# heatmap
norm_sig <- normalized_counts[which(rownames(normalized_counts) %in% sig_res$gene),] 
 
heat_colors <- paletteer_c("grDevices::Viridis", 6)
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Volcano plot
res_vol <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)
res_vol$de <- "Not sig"
res_vol$de[res_vol$padj <=0.05 & res_vol$log2FoldChange >= 1] <- "Up"
res_vol$de[res_vol$padj < 0.05 & res_vol$log2FoldChange <= -1] <- "Down"
res_vol$de <- factor(res_vol$de, levels = c("Up", "Down", "Not sig"))
n_label <- 5
n_up <- sum(res_vol$de == "Up")
n_down <- sum(res_vol$de == "Down")
res_vol$label <- NA
if (n_up > 0) {
  top_up <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = -padj)
  top_up_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = abs(log2FoldChange))
  top_up <- rbind(top_up, top_up_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_up$gene] <- "Yes"
}
if (n_down > 0) {
  top_down <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = -padj)
  top_down_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = abs(log2FoldChange))
  top_down <- rbind(top_down, top_down_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_down$gene] <- "Yes"
}
ggplot(res_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = de)) +
  ggtitle("astrocyte AD vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::scale_colour_manual(name = NULL,
                               aesthetics = c("colour", "fill"),
                               values = c("#DC0000FF", "#3C5488FF", "grey"),
                               label = c("Up-regulated", "Down-regulated"),
                               breaks = c("Up", "Down")) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
  ggrepel::geom_text_repel(
    data = res_vol,
    ggplot2::aes(log2FoldChange, y = -log10(padj), label = ifelse(label == "Yes", as.character(.data[["gene"]]), "")),
    max.iter = 1000, size = 4, na.rm = TRUE) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black", size = 16),
    axis.title = ggplot2::element_text(color = "black", size = 18),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    legend.text = ggplot2::element_text(size = 12, colour = "black"),
    legend.position = "top",
    legend.direction="horizontal",
    plot.margin = ggplot2::margin(c(1, 1, 1, 1), unit = "cm")
  )

### DGE in Endothelial
# import counts files

list.files("~/SCADSplice/rnaseq/results/endo/results/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/endo/results/salmon/*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

txi <- tximport(files, 
                type="salmon", 
                tx2gene=tx2gene[,c("transcript_id", "gene_name")], 
                countsFromAbundance = "lengthScaledTPM")


sampletype <- factor(c(rep("EAD",4), rep("ECt", 4)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

all(colnames(txi$counts) == rownames(metadata))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta, 
                                design = ~ sampletype)

# Normalize counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="~/SCADSplice/rnaseq/results/dge/endo_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, 
         annotation = meta)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Specify contrast for comparison of interest
contrast <- c("sampletype", "EAD", "ECt")

# Output results of Wald test for contrast
res <- DESeq2::results(dds, 
               contrast = contrast, 
               alpha = 0.05)

# Shrink the log2 fold changes to be more accurate
res <- lfcShrink(dds, 
                 coef = "sampletype_ECt_vs_EAD", 
                 type = "apeglm")

# Set thresholds
padj.cutoff <- 0.05

# Turn the results object into a tibble 
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res <- filter(res_tbl, 
                  padj < padj.cutoff)

write.table(sig_res, file="~/endo_DGE.tsv", sep="\t", quote=F, col.names=NA)

# Visualizing

# heatmap
norm_sig <- normalized_counts[which(rownames(normalized_counts) %in% sig_res$gene),] 
 
heat_colors <- paletteer_c("grDevices::Viridis", 6)
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Volcano plot
res_vol <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)
res_vol$de <- "Not sig"
res_vol$de[res_vol$padj <=0.05 & res_vol$log2FoldChange >= 1] <- "Up"
res_vol$de[res_vol$padj < 0.05 & res_vol$log2FoldChange <= -1] <- "Down"
res_vol$de <- factor(res_vol$de, levels = c("Up", "Down", "Not sig"))
n_label <- 5
n_up <- sum(res_vol$de == "Up")
n_down <- sum(res_vol$de == "Down")
res_vol$label <- NA
if (n_up > 0) {
  top_up <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = -padj)
  top_up_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = abs(log2FoldChange))
  top_up <- rbind(top_up, top_up_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_up$gene] <- "Yes"
}
if (n_down > 0) {
  top_down <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = -padj)
  top_down_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = abs(log2FoldChange))
  top_down <- rbind(top_down, top_down_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_down$gene] <- "Yes"
}
ggplot(res_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = de)) +
  ggtitle("endocyte AD vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::scale_colour_manual(name = NULL,
                               aesthetics = c("colour", "fill"),
                               values = c("#DC0000FF", "#3C5488FF", "grey"),
                               label = c("Up-regulated", "Down-regulated"),
                               breaks = c("Up", "Down")) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
  ggrepel::geom_text_repel(
    data = res_vol,
    ggplot2::aes(log2FoldChange, y = -log10(padj), label = ifelse(label == "Yes", as.character(.data[["gene"]]), "")),
    max.iter = 1000, size = 4, na.rm = TRUE) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black", size = 16),
    axis.title = ggplot2::element_text(color = "black", size = 18),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    legend.text = ggplot2::element_text(size = 12, colour = "black"),
    legend.position = "top",
    legend.direction="horizontal",
    plot.margin = ggplot2::margin(c(1, 1, 1, 1), unit = "cm")
  )

### DGE in Microglia
# import counts files

list.files("~/SCADSplice/rnaseq/results/micro/results/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/micro/results/salmon/*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

txi <- tximport(files, 
                type="salmon", 
                tx2gene=tx2gene[,c("transcript_id", "gene_name")], 
                countsFromAbundance = "lengthScaledTPM")


sampletype <- factor(c(rep("MAD",4), rep("MCt", 4)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

all(colnames(txi$counts) == rownames(metadata))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta, 
                                design = ~ sampletype)

# Normalize counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="~/SCADSplice/rnaseq/results/dge/micro_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, 
         annotation = meta)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Specify contrast for comparison of interest
contrast <- c("sampletype", "MAD", "MCt")

# Output results of Wald test for contrast
res <- DESeq2::results(dds, 
               contrast = contrast, 
               alpha = 0.05)

# Shrink the log2 fold changes to be more accurate
res <- lfcShrink(dds, 
                 coef = "sampletype_MCt_vs_MAD", 
                 type = "apeglm")

# Set thresholds
padj.cutoff <- 0.05

# Turn the results object into a tibble 
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res <- filter(res_tbl, 
                  padj < padj.cutoff)

write.table(sig_res, file="~/micro_DGE.tsv", sep="\t", quote=F, col.names=NA)

# Visualizing

# heatmap
norm_sig <- normalized_counts[which(rownames(normalized_counts) %in% sig_res$gene),] 
 
heat_colors <- paletteer_c("grDevices::Viridis", 6)
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Volcano plot
res_vol <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)
res_vol$de <- "Not sig"
res_vol$de[res_vol$padj <=0.05 & res_vol$log2FoldChange >= 1] <- "Up"
res_vol$de[res_vol$padj < 0.05 & res_vol$log2FoldChange <= -1] <- "Down"
res_vol$de <- factor(res_vol$de, levels = c("Up", "Down", "Not sig"))
n_label <- 5
n_up <- sum(res_vol$de == "Up")
n_down <- sum(res_vol$de == "Down")
res_vol$label <- NA
if (n_up > 0) {
  top_up <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = -padj)
  top_up_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = abs(log2FoldChange))
  top_up <- rbind(top_up, top_up_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_up$gene] <- "Yes"
}
if (n_down > 0) {
  top_down <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = -padj)
  top_down_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = abs(log2FoldChange))
  top_down <- rbind(top_down, top_down_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_down$gene] <- "Yes"
}
ggplot(res_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = de)) +
  ggtitle("microcyte AD vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::scale_colour_manual(name = NULL,
                               aesthetics = c("colour", "fill"),
                               values = c("#DC0000FF", "#3C5488FF", "grey"),
                               label = c("Up-regulated", "Down-regulated"),
                               breaks = c("Up", "Down")) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
  ggrepel::geom_text_repel(
    data = res_vol,
    ggplot2::aes(log2FoldChange, y = -log10(padj), label = ifelse(label == "Yes", as.character(.data[["gene"]]), "")),
    max.iter = 1000, size = 4, na.rm = TRUE) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black", size = 16),
    axis.title = ggplot2::element_text(color = "black", size = 18),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    legend.text = ggplot2::element_text(size = 12, colour = "black"),
    legend.position = "top",
    legend.direction="horizontal",
    plot.margin = ggplot2::margin(c(1, 1, 1, 1), unit = "cm")
  )

### DGE in Neuron
# import counts files

list.files("~/SCADSplice/rnaseq/results/neuro/results/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/neuro/results/salmon/*/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

txi <- tximport(files, 
                type="salmon", 
                tx2gene=tx2gene[,c("transcript_id", "gene_name")], 
                countsFromAbundance = "lengthScaledTPM")


sampletype <- factor(c(rep("NAD",4), rep("NCt", 4)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

all(colnames(txi$counts) == rownames(metadata))

# Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta, 
                                design = ~ sampletype)

# Normalize counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="~/SCADSplice/rnaseq/results/dge/neuro_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA 
plotPCA(rld, intgroup="sampletype")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, 
         annotation = meta)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Specify contrast for comparison of interest
contrast <- c("sampletype", "NAD", "NCt")

# Output results of Wald test for contrast
res <- DESeq2::results(dds, 
               contrast = contrast, 
               alpha = 0.05)

# Shrink the log2 fold changes to be more accurate
res <- lfcShrink(dds, 
                 coef = "sampletype_NCt_vs_NAD", 
                 type = "apeglm")

# Set thresholds
padj.cutoff <- 0.05

# Turn the results object into a tibble 
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res <- filter(res_tbl, 
                  padj < padj.cutoff)

write.table(sig_res, file="~/neuro_DGE.tsv", sep="\t", quote=F, col.names=NA)

# Visualizing

# heatmap
norm_sig <- normalized_counts[which(rownames(normalized_counts) %in% sig_res$gene),] 
 
heat_colors <- paletteer_c("grDevices::Viridis", 6)
pheatmap(norm_sig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Volcano plot
res_vol <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 1)
res_vol$de <- "Not sig"
res_vol$de[res_vol$padj <=0.05 & res_vol$log2FoldChange >= 1] <- "Up"
res_vol$de[res_vol$padj < 0.05 & res_vol$log2FoldChange <= -1] <- "Down"
res_vol$de <- factor(res_vol$de, levels = c("Up", "Down", "Not sig"))
n_label <- 5
n_up <- sum(res_vol$de == "Up")
n_down <- sum(res_vol$de == "Down")
res_vol$label <- NA
if (n_up > 0) {
  top_up <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = -padj)
  top_up_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Up") %>%
    dplyr::top_n(min(n_up, n_label), wt = abs(log2FoldChange))
  top_up <- rbind(top_up, top_up_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_up$gene] <- "Yes"
}
if (n_down > 0) {
  top_down <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = -padj)
  top_down_log2FoldChange <- res_vol %>%
    dplyr::filter(de == "Down") %>%
    dplyr::top_n(min(n_down, n_label), wt = abs(log2FoldChange))
  top_down <- rbind(top_down, top_down_log2FoldChange)
  res_vol$label[res_vol$gene %in% top_down$gene] <- "Yes"
}
ggplot(res_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = de)) +
  ggtitle("neurocyte AD vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  ggplot2::geom_vline(xintercept = c(-log2(2), log2(2)),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = 2, size = 0.5, alpha = 0.5) +
  ggplot2::scale_colour_manual(name = NULL,
                               aesthetics = c("colour", "fill"),
                               values = c("#DC0000FF", "#3C5488FF", "grey"),
                               label = c("Up-regulated", "Down-regulated"),
                               breaks = c("Up", "Down")) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3))) +
  ggrepel::geom_text_repel(
    data = res_vol,
    ggplot2::aes(log2FoldChange, y = -log10(padj), label = ifelse(label == "Yes", as.character(.data[["gene"]]), "")),
    max.iter = 1000, size = 4, na.rm = TRUE) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black", size = 16),
    axis.title = ggplot2::element_text(color = "black", size = 18),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    legend.text = ggplot2::element_text(size = 12, colour = "black"),
    legend.position = "top",
    legend.direction="horizontal",
    plot.margin = ggplot2::margin(c(1, 1, 1, 1), unit = "cm")
  )
