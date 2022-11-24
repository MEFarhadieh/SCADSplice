setwd("~/SCADSplice/rnaseq/results/dtu")
library(DTUrtle)
biocpar <- BiocParallel::MulticoreParam(24) # depends on your machine cores

### Following codes were used to identify DTU in Astrocyte
# Importing and processing GTF annotation

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

# Reading in quantification data

list.files("~/SCADSplice/rnaseq/results/astro/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/astro/salmon/*.quant/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

cts <- import_counts(files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])
rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

# Sample metadata

pd <- data.frame("id"=colnames(cts), 
                 "group"=c(rep("AAD",6), rep("ACt",6)), 
                 stringsAsFactors = FALSE)

# DTU analysis

dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                       cond_col = "group", cond_levels = c("AAD", "ACt"), filtering_strategy = "bulk", 
                       BPPARAM = biocpar)

dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

# DTU table creation

dturtle <- create_dtu_table(
  dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
  add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))

# Gene proportion barplot

dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)


# Gene proportion heatmap

dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = TRUE,
                                    treeheight_col=20,
                                    savepath = "images", 
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

# Transcript overview

dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                 gtf = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf", 
                                 genome = 'hg38', 
                                 one_to_one = TRUE,
                                 savepath = "images", 
                                 add_to_table = "transcript_view",
                                 BPPARAM = biocpar)

# Visualize DTU table

column_formatter_list <- list(
  "gene_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "number_tx" = formattable::color_tile('white', "lightblue"),
  "number_significant_tx" = formattable::color_tile('white', "lightblue"),
  "max(AAD-ACt)" = table_percentage_bar('lightgreen', "#FF9999", digits=2),
  "tx_expr_in_max" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2))

plot_dtu_table(dturtle = dturtle,
               savepath = "~/SCADSplice/rnaseq/results/dtu/astro_DTU_DGE.html", 
               column_formatters = column_formatter_list)

write.table(dturtle$sig_gene,
            "~/SCADSplice/rnaseq/results/dtu/astro_genelist_dtu.txt",
            col.names = F, row.names = F, sep = "\t",
            quote = F)

### Following codes were used to identify DTU in Endothelial
# Importing and processing GTF annotation

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

# Reading in quantification data

list.files("~/SCADSplice/rnaseq/results/endo/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/endo/salmon/*.quant/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

cts <- import_counts(files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])
rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

# Sample metadata

pd <- data.frame("id"=colnames(cts), 
                 "group"=c(rep("EAD",6), rep("ECt",6)), 
                 stringsAsFactors = FALSE)

# DTU analysis

dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                       cond_col = "group", cond_levels = c("EAD", "ECt"), filtering_strategy = "bulk", 
                       BPPARAM = biocpar)

dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

# DTU table creation

dturtle <- create_dtu_table(
  dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
  add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))

# Gene proportion barplot

dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)


# Gene proportion heatmap

dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = TRUE,
                                    treeheight_col=20,
                                    savepath = "images", 
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

# Transcript overview

dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                 gtf = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf", 
                                 genome = 'hg38', 
                                 one_to_one = TRUE,
                                 savepath = "images", 
                                 add_to_table = "transcript_view",
                                 BPPARAM = biocpar)

# Visualize DTU table

column_formatter_list <- list(
  "gene_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "number_tx" = formattable::color_tile('white', "lightblue"),
  "number_significant_tx" = formattable::color_tile('white', "lightblue"),
  "max(EAD-ECt)" = table_percentage_bar('lightgreen', "#FF9999", digits=2),
  "tx_expr_in_max" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2))


plot_dtu_table(dturtle = dturtle,
               savepath = "~/SCADSplice/rnaseq/results/dtu/endo_DTU_DGE.html", 
               column_formatters = column_formatter_list)

write.table(dturtle$sig_gene,
            "~/SCADSplice/rnaseq/results/dtu/endo_genelist_dtu.txt",
            col.names = F, row.names = F, sep = "\t",
            quote = F)

### Following codes were used to identify DTU in Microglia
# Importing and processing GTF annotation

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

# Reading in quantification data

list.files("~/SCADSplice/rnaseq/results/micro/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/micro/salmon/*.quant/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

cts <- import_counts(files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])
rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

# Sample metadata

pd <- data.frame("id"=colnames(cts), 
                 "group"=c(rep("MAD",6), rep("MCt",6)), 
                 stringsAsFactors = FALSE)

# DTU analysis

dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                       cond_col = "group", cond_levels = c("MAD", "MCt"), filtering_strategy = "bulk", 
                       BPPARAM = biocpar)

dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

# DTU table creation

dturtle <- create_dtu_table(
  dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
  add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))

# Gene proportion barplot

dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)


# Gene proportion heatmap

dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = TRUE,
                                    treeheight_col=20,
                                    savepath = "images", 
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

# Transcript overview

dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                 gtf = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf", 
                                 genome = 'hg38', 
                                 one_to_one = TRUE,
                                 savepath = "images", 
                                 add_to_table = "transcript_view",
                                 BPPARAM = biocpar)

# Visualize DTU table

column_formatter_list <- list(
  "gene_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "number_tx" = formattable::color_tile('white', "lightblue"),
  "number_significant_tx" = formattable::color_tile('white', "lightblue"),
  "max(MAD-MCt)" = table_percentage_bar('lightgreen', "#FF9999", digits=2),
  "tx_expr_in_max" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2))


plot_dtu_table(dturtle = dturtle,
               savepath = "~/SCADSplice/rnaseq/results/dtu/micro_DTU_DGE.html", 
               column_formatters = column_formatter_list)

write.table(dturtle$sig_gene,
            "~/SCADSplice/rnaseq/results/dtu/micro_genelist_dtu.txt",
            col.names = F, row.names = F, sep = "\t",
            quote = F)

### Following codes were used to identify DTU in Neuron
# Importing and processing GTF annotation

tx2gene <- import_gtf(
  gtf_file = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
tx2gene$gene_name <- one_to_one_mapping(name = tx2gene$gene_name, id = tx2gene$gene_id)
tx2gene$transcript_name <- one_to_one_mapping(name = tx2gene$transcript_name, id = tx2gene$transcript_id)

tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_name", "gene_name"))

# Reading in quantification data

list.files("~/SCADSplice/rnaseq/results/neuro/salmon/")
files <- Sys.glob("~/SCADSplice/rnaseq/results/neuro/salmon/*.quant/quant.sf")
names(files) <- gsub(".*/","",gsub("/quant.sf","",files))

cts <- import_counts(files, type = "salmon", tx2gene=tx2gene[,c("transcript_id", "gene_name")])
rownames(cts) <- tx2gene$transcript_name[match(rownames(cts), tx2gene$transcript_id)]

# Sample metadata

pd <- data.frame("id"=colnames(cts), 
                 "group"=c(rep("MAD",6), rep("MCt",6)), 
                 stringsAsFactors = FALSE)

# DTU analysis

dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                       cond_col = "group", cond_levels = c("MAD", "MCt"), filtering_strategy = "bulk", 
                       BPPARAM = biocpar)

dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

# DTU table creation

dturtle <- create_dtu_table(
  dturtle = dturtle, add_gene_metadata = list("chromosome"="seqnames"), 
  add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))

# Gene proportion barplot

dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)


# Gene proportion heatmap

dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = TRUE,
                                    treeheight_col=20,
                                    savepath = "images", 
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

# Transcript overview

dturtle <- plot_transcripts_view(dturtle = dturtle, 
                                 gtf = "~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf", 
                                 genome = 'hg38', 
                                 one_to_one = TRUE,
                                 savepath = "images", 
                                 add_to_table = "transcript_view",
                                 BPPARAM = biocpar)

# Visualize DTU table

column_formatter_list <- list(
  "gene_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 6),
  "number_tx" = formattable::color_tile('white', "lightblue"),
  "number_significant_tx" = formattable::color_tile('white', "lightblue"),
  "max(MAD-MCt)" = table_percentage_bar('lightgreen', "#FF9999", digits=2),
  "tx_expr_in_max" = table_percentage_bar('white', "lightblue", color_break = 0, digits=2))


plot_dtu_table(dturtle = dturtle,
               savepath = "~/SCADSplice/rnaseq/results/dtu/neuro_DTU_DGE.html", 
               column_formatters = column_formatter_list)

write.table(dturtle$sig_gene,
            "~/SCADSplice/rnaseq/results/dtu/neuro_genelist_dtu.txt",
            col.names = F, row.names = F, sep = "\t",
            quote = F)
