library(SingleCellExperiment)
library(scFlow)
library(scFlowData)

ensembl_fp <- system.file("extdata","ensembl_mappings.tsv",package="scFlowData")
ctd_fp <- system.file("extdata","ctd",package="scFlowData")
manifest <- read.delim(manifest_fp)

# read all sce individuals together
dir_sce <- "~/nft_qc_results/sce_individual"

sce_path <- dir(
  path = dir_sce,
  pattern = "sce_",
  full.names = TRUE
)

sce_pathlist <- list()

for (i in sce_path) {
  sce_pathlist[[i]] <- i
}

sce_list <- lapply(sce_pathlist, read_sce)

# extract barcodes of each cells as a whitelist file for umi_tools
for (i in sce_path) {
  write.table(
    sce_list[[i]]@colData@rownames,
    file = paste0("~/metadata/nft_whitelist/",
                  sce_list[[i]]@colData@listData[["manifest"]][1],
                  "_whitelist.txt"),
    quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
  )
}

# merge of post-QC samples
sce <- merge_sce(
  sce_list,
  ensembl_mapping_file = ensembl_fp
)

# generate merged report
sce <- annotate_merged_sce(
  sce,
  plot_vars = c("total_features_by_counts","total_counts","pc_mito","pc_ribo"),
  unique_id_var = "manifest",
  facet_vars = "diagnosis",
  outlier_vars = c("total_features_by_counts", "total_counts")
)

report_merged_sce(sce)

write_sce(
  sce,
  "~/nft_merge/sce_merged"
)
