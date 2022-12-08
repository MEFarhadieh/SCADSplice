/*
 * -------------------------------------------------
 *  nf-core/readzs Nextflow config file
 * -------------------------------------------------
 *  config options for conda environment on a local machine.
 */

 process {
   executor = 'local'
   memory = { 1.GB * task.attempt }
   time = { 1.h * task.attempt }
   errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
   maxRetries = 3
 }

 params {
     input = "~/SCADSplice/scrnaseq/nge/metadata/nge_samplesheet.csv"
     useChannels = true
     chr_lengths = "~/SCADSplice/references/chrLength.txt"
     libType = "10X"
     isSICILIAN = true
     isCellranger = false
     runName = "nge_ReadZs"
     metadata = "~/SCADSplice/scrnaseq/nge/metadata/nge_cells_annotation.tsv"
     ontologyCols = "'group, type'"
     annotation_bed = "~/SCADSplice/references/hg38_genes.bed"
     gff = "~/SCADSplice/references/gencode.v39.annotation.gff3"
     max_memory                 = '126.GB'
     max_cpus                   = 24
     max_time                   = '360.h'
 }


params.outdir = "~/SCADSplice/scrnaseq/nge/results/${params.runName}"
params.tracedir = "~/SCADSplice/scrnaseq/nge/results/${params.runName}/pipeline_info"
