### Following codes were used to identify statistics of ASE
### after performing ASE in shell

library(maser)
library(rtracklayer)

### Astrocyte ASE evaluation 
#Importing splicing events

path <-("~/SCADSplice/rnaseq/results/astro/splicing/output/")
AD_raw <- maser(path, c("AAD", "ACt"), ftype = "JC")
AD_raw

# Filtering events

AD_filt <- filterByCoverage(AD_raw, avg_reads = 5)
AD_filt

# Top events by FDR < 0.05 & deltaPSI > 0.1

AD_top <- topEvents(AD_filt, fdr = 0.05, deltaPSI = 0.1)
AD_top

save(AD, AD_filt, AD_top,file="~/SCADSplice/rnaseq/results/astro/bulk_ASE_astro.Rdata")

# Gene specific events for example RUFY1

AD_RUFY1 <- geneEvents(AD_filt, geneS = "RUFY1", fdr = 0.05, deltaPSI = 0.1)
AD_RUFY1

maser::display(AD_RUFY1, "SE")

plotGenePSI(AD_RUFY1, type = "SE", show_replicates = T)

# Global splicing

volcano(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
volcano(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)


pca(AD_top)
boxplot_PSI_levels(AD_top, type = "RI") # type = c(A3SS,A5SS,MXE,RI,SE)
splicingDistribution(AD_top)

# visualization splicing events

gtf_path <- ("~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
ens_gtf <- rtracklayer::import.gff(gtf_path)

Atr_events <- geneEvents(AD_filt, geneS = "PDE8A", fdr = 0.05,
                         deltaPSI = 0.1)
maser::display(Atr_events, "MXE")
plotTranscripts(Atr_events, type = "MXE", event_id = 1181,
                gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)

AD_top@A3SS_events$ASE <- "A3SS"
AD_top@A5SS_events$ASE <- "A5SS"
AD_top@A5SS_events$ASE <- "MXE"
AD_top@RI_events$ASE <- "RI"
AD_top@SE_events$ASE <- "SE"

ase_sig <- rbind(merge(AD_top@A3SS_stats, AD_top@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@A5SS_stats, AD_top@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@MXE_stats, AD_top@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@RI_stats, AD_top@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@SE_stats, AD_top@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             

write.table(ase_sig,
            "~/astro_sig_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)

AD_raw@A3SS_events$ASE <- "A3SS"
AD_raw@A5SS_events$ASE <- "A5SS"
AD_raw@A5SS_events$ASE <- "MXE"
AD_raw@RI_events$ASE <- "RI"
AD_raw@SE_events$ASE <- "SE"

ase_raw <- rbind(merge(AD_raw@A3SS_stats, AD_raw@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@A5SS_stats, AD_raw@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@MXE_stats, AD_raw@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@RI_stats, AD_raw@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@SE_stats, AD_raw@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             
write.table(ase_raw,
            "~/astro_raw_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)
            
### Endothelial ASE evaluation
#Importing splicing events

path <-("~/SCADSplice/rnaseq/results/endo/splicing/output/")
AD_raw <- maser(path, c("EAD", "ECt"), ftype = "JC")
AD_raw

# Filtering events

AD_filt <- filterByCoverage(AD_raw, avg_reads = 5)
AD_filt

# Top events by FDR < 0.05 & deltaPSI > 0.1

AD_top <- topEvents(AD_filt, fdr = 0.05, deltaPSI = 0.1)
AD_top

save(AD, AD_filt, AD_top,file="~/SCADSplice/rnaseq/results/endo/bulk_ASE_endo.Rdata")

# Gene specific events for example PGS1

AD_PGS1 <- geneEvents(AD_filt, geneS = "PGS1", fdr = 0.05, deltaPSI = 0.1)
AD_PGS1

maser::display(AD_PGS1, "SE")

plotGenePSI(AD_PGS1, type = "SE", show_replicates = T)

# Global splicing

volcano(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
volcano(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)


pca(AD_top)
boxplot_PSI_levels(AD_top, type = "RI") # type = c(A3SS,A5SS,MXE,RI,SE)
splicingDistribution(AD_top)

# visualization splicing events

gtf_path <- ("~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
ens_gtf <- rtracklayer::import.gff(gtf_path)

Atr_events <- geneEvents(AD_filt, geneS = "PGS1", fdr = 0.05,
                         deltaPSI = 0.1)
maser::display(Atr_events, "SE")
plotTranscripts(Atr_events, type = "SE", event_id = 1327,
                gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)

AD_top@A3SS_events$ASE <- "A3SS"
AD_top@A5SS_events$ASE <- "A5SS"
AD_top@A5SS_events$ASE <- "MXE"
AD_top@RI_events$ASE <- "RI"
AD_top@SE_events$ASE <- "SE"

ase_sig <- rbind(merge(AD_top@A3SS_stats, AD_top@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@A5SS_stats, AD_top@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@MXE_stats, AD_top@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@RI_stats, AD_top@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@SE_stats, AD_top@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             

write.table(ase_sig,
            "~/endo_sig_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)

AD_raw@A3SS_events$ASE <- "A3SS"
AD_raw@A5SS_events$ASE <- "A5SS"
AD_raw@A5SS_events$ASE <- "MXE"
AD_raw@RI_events$ASE <- "RI"
AD_raw@SE_events$ASE <- "SE"

ase_raw <- rbind(merge(AD_raw@A3SS_stats, AD_raw@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@A5SS_stats, AD_raw@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@MXE_stats, AD_raw@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@RI_stats, AD_raw@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@SE_stats, AD_raw@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             
write.table(ase_raw,
            "~/endo_raw_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)
            
### Microglia ASE evaluation
#Importing splicing events

path <-("~/SCADSplice/rnaseq/results/micro/splicing/output/")
AD_raw <- maser(path, c("MAD", "MCt"), ftype = "JC")
AD_raw

# Filtering events

AD_filt <- filterByCoverage(AD_raw, avg_reads = 5)
AD_filt

# Top events by FDR < 0.05 & deltaPSI > 0.1

AD_top <- topEvents(AD_filt, fdr = 0.05, deltaPSI = 0.1)
AD_top

save(AD, AD_filt, AD_top,file="~/SCADSplice/rnaseq/results/micro/bulk_ASE_micro.Rdata")

# Gene specific events for example DROSHA

AD_DROSHA <- geneEvents(AD_filt, geneS = "DROSHA", fdr = 0.05, deltaPSI = 0.1)
AD_DROSHA

maser::display(AD_DROSHA, "SE")

plotGenePSI(AD_DROSHA, type = "SE", show_replicates = T)

# Global splicing

volcano(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
volcano(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)


pca(AD_top)
boxplot_PSI_levels(AD_top, type = "RI") # type = c(A3SS,A5SS,MXE,RI,SE)
splicingDistribution(AD_top)

# visualization splicing events

gtf_path <- ("~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
ens_gtf <- rtracklayer::import.gff(gtf_path)

Atr_events <- geneEvents(AD_filt, geneS = "DROSHA", fdr = 0.05,
                         deltaPSI = 0.1)
maser::display(Atr_events, "SE")
plotTranscripts(Atr_events, type = "SE", event_id = 539,
                gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)

AD_top@A3SS_events$ASE <- "A3SS"
AD_top@A5SS_events$ASE <- "A5SS"
AD_top@A5SS_events$ASE <- "MXE"
AD_top@RI_events$ASE <- "RI"
AD_top@SE_events$ASE <- "SE"

ase_sig <- rbind(merge(AD_top@A3SS_stats, AD_top@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@A5SS_stats, AD_top@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@MXE_stats, AD_top@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@RI_stats, AD_top@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@SE_stats, AD_top@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             

write.table(ase_sig,
            "~/micro_sig_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)

AD_raw@A3SS_events$ASE <- "A3SS"
AD_raw@A5SS_events$ASE <- "A5SS"
AD_raw@A5SS_events$ASE <- "MXE"
AD_raw@RI_events$ASE <- "RI"
AD_raw@SE_events$ASE <- "SE"

ase_raw <- rbind(merge(AD_raw@A3SS_stats, AD_raw@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@A5SS_stats, AD_raw@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@MXE_stats, AD_raw@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@RI_stats, AD_raw@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@SE_stats, AD_raw@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             
write.table(ase_raw,
            "~/micro_raw_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)
            
### Neuron ASE evaluation
#Importing splicing events

path <-("~/SCADSplice/rnaseq/results/neuro/splicing/output/")
AD_raw <- maser(path, c("NAD", "NCt"), ftype = "JC")
AD_raw

# Filtering events

AD_filt <- filterByCoverage(AD_raw, avg_reads = 5)
AD_filt

# Top events by FDR < 0.05 & deltaPSI > 0.1

AD_top <- topEvents(AD_filt, fdr = 0.05, deltaPSI = 0.1)
AD_top

save(AD, AD_filt, AD_top,file="~/SCADSplice/rnaseq/results/neuro/bulk_ASE_neuro.Rdata")

# Gene specific events for example CLPX

AD_CLPX <- geneEvents(AD_filt, geneS = "CLPX", fdr = 0.05, deltaPSI = 0.1)
AD_CLPX

maser::display(AD_CLPX, "SE")

plotGenePSI(AD_CLPX, type = "SE", show_replicates = T)

# Global splicing

volcano(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_filt, fdr = 0.05, deltaPSI = 0.1, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
dotplot(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)
volcano(AD_top, type = "SE") # type = c(A3SS,A5SS,MXE,RI,SE)


pca(AD_top)
boxplot_PSI_levels(AD_top, type = "RI") # type = c(A3SS,A5SS,MXE,RI,SE)
splicingDistribution(AD_top)

# visualization splicing events

gtf_path <- ("~/SCADSplice/rnaseq/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf")
ens_gtf <- rtracklayer::import.gff(gtf_path)

Atr_events <- geneEvents(AD_filt, geneS = "CLPX", fdr = 0.05,
                         deltaPSI = 0.1)
maser::display(Atr_events, "SE")
plotTranscripts(Atr_events, type = "SE", event_id = 26087,
                gtf = ens_gtf, zoom = FALSE, show_PSI = TRUE)

AD_top@A3SS_events$ASE <- "A3SS"
AD_top@A5SS_events$ASE <- "A5SS"
AD_top@A5SS_events$ASE <- "MXE"
AD_top@RI_events$ASE <- "RI"
AD_top@SE_events$ASE <- "SE"

ase_sig <- rbind(merge(AD_top@A3SS_stats, AD_top@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@A5SS_stats, AD_top@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@MXE_stats, AD_top@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@RI_stats, AD_top@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_top@SE_stats, AD_top@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             

write.table(ase_sig,
            "~/neuro_sig_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)

AD_raw@A3SS_events$ASE <- "A3SS"
AD_raw@A5SS_events$ASE <- "A5SS"
AD_raw@A5SS_events$ASE <- "MXE"
AD_raw@RI_events$ASE <- "RI"
AD_raw@SE_events$ASE <- "SE"

ase_raw <- rbind(merge(AD_raw@A3SS_stats, AD_raw@A3SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@A5SS_stats, AD_raw@A5SS_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@MXE_stats, AD_raw@MXE_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@RI_stats, AD_raw@RI_events, by.x = "ID",by.y = "ID", fill=-99999),
             merge(AD_raw@SE_stats, AD_raw@SE_events, by.x = "ID",by.y = "ID", fill=-99999))
             
write.table(ase_raw,
            "~/neuro_raw_dtu.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)
            
