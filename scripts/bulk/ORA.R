### ORA for significant DTU or DGE or both

library(readr)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(WebGestaltR)
library(paletteer)
library(cowplot)

# unify gene list with ensembl_gene_id for DTU

load("~/SCADSplice/rnaseq/results/dtu/astro_dtu.RData")
raw_dtu <- dturtle[["meta_table_gene"]]$gene_id.1
sig_dtu <- read_tsv("~/SCADSplice/rnaseq/results/dtu/astro_genelist_dtu.txt")

annotations_ahb <- read_csv("~/SCADSplice/rnaseq/references/annotations_ahb.csv")

raw <- data.frame(unique(c(raw_dtu$gene_ID))
colnames(raw) <- "genes"

sig <- data.frame(unique(c(sig_dtu$gene_ID))
colnames(sig) <- "genes"

raw <- dplyr::left_join(raw, annotations_ahb, 
                        by=c("genes" ="gene_name"))
raw <- na.omit(raw)
all_gene <- as.character(unique(raw$gene_id))

sig <- dplyr::left_join(sig, annotations_ahb, 
                        by=c("genes" ="gene_name"))
sig <- na.omit(sig)
sig_gene <- as.character(c(sig$gene_id, ase$ensembl_gene_id))


# ORA with ego & WebGestaltR

ego_bp <- enrichGO(gene = sig_gene, 
                universe = all_gene,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego_cc <- enrichGO(gene = sig_gene, 
                universe = all_gene,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "CC", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
                
ego_mf <- enrichGO(gene = sig_gene, 
                universe = all_gene,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "MF", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

ego_all <- enrichGO(gene = sig_gene, 
                universe = all_gene,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = c("BP", "CC", "MF") 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
                                                
res <- WebGestaltR::WebGestaltR(
  enrichMethod = "ORA",
  organism = "hsapiens",
  enrichDatabase = c(
    "geneontology_Biological_Process",
    "geneontology_Cellular_Component",
    "geneontology_Molecular_Function",
    "pathway_KEGG",
    "pathway_Panther",
    "pathway_Reactome",
    "pathway_Wikipathway"),
  interestGene = sig_gene,
  interestGeneType = "ensembl_gene_id",
  referenceGene = all_gene,
  referenceGeneType = "ensembl_gene_id",
  isOutput = FALSE,
  projectName = NULL)

write.table(ego_bp@result,
            "~/SCADSplice/rnaseq/results/ora/astro_ora.txt",
            sep = "\t", quote = F)  

# visualization

dotplot(ego_all, showCategory = 20)  
enrichplot::cnetplot(ego_bp, circular = T, colorEdge = T,showCategory = 20)
d <- setNames(dge_sig$logFC, dge_sig$ensembl_gene_id)
enrichplot::cnetplot(ego_bp, circular = F, colorEdge = T,showCategory = 20,
                     node_label="gene",foldChange=d, color_category="#BC7A3C") +
  scale_color_gradient2(name='associated data', low='#26456E', high='#9C0824')
enrichplot::cnetplot(ego_bp, circular = F, colorEdge = F,showCategory = 20)

go_table <- read.delim("~/SCADSplice/rnaseq/results/ora/astro_ora_modified.txt")
go_genes <- data.frame(go_table$genes, go_table$source)
colnames(go_genes) <- c("ID", "source")
go_terms <- unique(go_table$term)
crd <- chord_dat(data = go_table, genes = go_genes, process = go_terms)

GOChord(crd, space = 0.02, gene.order = 'source', gene.space = 0.25, gene.size = 5) +
        scale_color_gradient(low='black', high='grey50')
