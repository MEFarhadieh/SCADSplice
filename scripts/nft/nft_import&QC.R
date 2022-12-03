library(SingleCellExperiment)
library(scFlow)
library(scFlowData)
library(DropletUtils)
library(ggplot2)
library(dplyr)

#import data and set variables

outputDir <- "~/SCADSplice/scrnaseq/nft/results/GSE129308"
outdir <-paste0(outputDir,"/qc_results")
dir.create(outdir, showWarnings = FALSE)
manifest_fp <- "~/SCADSplice/scrnaseq/nftmetadata/ss_manifest.tsv"
samplesheet_fp <- "~/SCADSplice/scrnaseq/nftmetadata//ss_samplesheet.tsv"
ensembl_fp <- system.file("extdata","ensembl_mappings.tsv",package="scFlowData")
ctd_fp <- system.file("extdata","ctd",package="scFlowData")
sembl <- read.table(ensembl_fp, header = TRUE)

j <- 1
manifest <- read.delim(manifest_fp)
samplesheet <- read.delim(samplesheet_fp)
dir_list <- manifest$filepath
dir_list

# create SingleCellExoeiment objects and QCs with annotaitons
# we did not run emptydrop() due to low counts had been filterd before

for (i in dir_list) {
  mat <- Read10X_h5(i, use.names = TRUE, unique.features = TRUE)
  gene.symb <- mat@Dimnames[[1]]
  ges <- data.frame(gene.symb)
  ges$embl <- sembl$ensembl_gene_id[
    match(ges$gene.symb, sembl$external_gene_name)
    ]
  ges$embl <- as.character(ges$embl)
  for (k in ges[is.na(ges)]){
    ges[is.na(ges)] <- ensembl$ensembl_gene_id[which(ensembl$ensembl_gene_id %in% ges[is.na(ges)][k])]
  }
  mat@Dimnames[[1]] <- ges$embl
  var_classes <- c(
    age = "factor",
    PMI = "factor",
    RIN = "factor"
  )
  metadata <- read_metadata(
    unique_key = manifest$key[j],
    key_colname = "manifest",
    samplesheet_path = samplesheet_fp,
    col_classes = var_classes
  )
  sce <- generate_sce(mat, metadata)

  sce <- annotate_sce(
    sce,
    min_library_size = 500,
    max_library_size = "adaptive",
    min_features = 200,
    max_features = "adaptive",
    max_mito = 0.1,
    min_ribo = 0,
    max_ribo = 1,
    min_cells = 3,
    drop_mito = TRUE,
    drop_ribo = FALSE,
    ensembl_mapping_file = ensembl_fp
  )
  sce <- filter_sce(sce)
  sce <- find_singlets(sce, "doubletfinder", pK = 0.005,
                       vars_to_regress_out = c("nCount_RNA", "pc_mito"),
                       num.cores = 1)
  sce <- filter_sce(sce)

  dir_report <- file.path(outdir, "qc_report")
  dir.create(dir_report, showWarnings = FALSE)
  report_qc_sce(
    sce = sce,
    report_folder_path = file.path(dir_report),
    report_file = paste0(manifest$key[j],"_qc_report")
  )
  dir_sce <- file.path(outdir, "sce_individual")
  dir.create(dir_sce, showWarnings = FALSE)
  write_sce(sce = sce,
            folder_path = file.path(dir_sce,
                                    paste("sce", basename(i), sep = "_")),
            overwrite = TRUE)
  j <- j+1
}
