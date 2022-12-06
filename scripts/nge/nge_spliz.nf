/*
 * -------------------------------------------------
 *  nf-core/spliz Nextflow config file
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
  dataname = "nge_all"
  input_file = "~/SCADSplice/scrnaseq/nge/results/nge_prespliz.tsv"
  SICILIAN = true
  pin_S = 0.1
  pin_z = 0.0
  bounds = 5
  light = false
  svd_type = "normdonor"
  n_perms = 100
  grouping_level_1 = "type"
  grouping_level_2 = "group"
  libraryType = "10X"
  run_analysis = true
  max_memory = 126.GB
  max_cpus = 24
  max_time = 360.h
}

params.outdir = "~/SCADSplice/scrnaseq/nge/results/${params.dataname}"
params.tracedir = "~/SCADSplice/scrnaseq/nge/results/${params.dataname}/pipeline_info"
params.schema_ignore_params = "input,single_end,show_hidden_params,validate_params,igenomes_ignore,tracedir,igenomes_base,help,monochrome_logs,plaintext_email,max_multiqc_email_size,email_on_fail,email,multiqc_config,publish_dir_mode,genome,genomes"
