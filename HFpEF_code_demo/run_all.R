#!/usr/bin/env Rscript
# Terminal: Rscript run_all.R
# RStudio: open HFpEF_code_demo.Rproj, then source("run_all.R").
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
if (!file.exists(file.path(root, "config.R"))) stop("Open the project folder before sourcing run_all.R")
arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i+1L], "--")) stop("Missing value for ", flag)
  args[i+1L]
}
if ("--help" %in% args) {
  cat("Usage: Rscript run_all.R [--data DIRECTORY] [--out DIRECTORY]\n")
} else {
  source(file.path(root, "config.R"))
  for (file in sort(list.files(file.path(root,"R"), pattern="\\.R$", full.names=TRUE))) source(file)
  require_packages()
  data_dir <- normalizePath(arg("--data", file.path(root,"data","synthetic")), mustWork=TRUE)
  default_out <- file.path(root,"results","demo")
  if (dir.exists(default_out) && length(list.files(default_out,all.files=TRUE,no..=TRUE)))
    default_out <- file.path(root,"results",paste0("demo_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  out_dir <- arg("--out", default_out)
  assert(!dir.exists(out_dir) || length(list.files(out_dir, all.files=TRUE, no..=TRUE)) == 0,
         "Output directory is not empty. Choose a new --out directory to preserve earlier results.")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  out_dir <- normalizePath(out_dir)
  label <- if (file.exists(file.path(data_dir,"SYNTHETIC_DATA.txt"))) "SYNTHETIC DEMONSTRATION" else "USER-SUPPLIED DATA"
  set.seed(config$seed)
  message("Validating inputs (",label,")")
  input <- load_inputs(data_dir,config)
  message("1/7 Paired limma and exercise-response overlaps")
  results <- run_paired_limma(input,config,out_dir)
  run_overlap_ternary(results,config,out_dir,label)
  message("2/7 Tissue-associated exercise response")
  run_tissue_plots(input,results,config,out_dir,label)
  message("3/7 Unadjusted resting comparisons (reconstructed)")
  run_resting_comparisons(input,config,out_dir,label)
  message("4/7 PCA and PLS-DA")
  run_pca_plsda(input,config,out_dir,label)
  message("5/7 CAMERA pathways (reconstructed)")
  run_camera(input,config,out_dir,label)
  message("6/7 Table S14 medication sensitivity (reconstructed)")
  run_medication_sensitivity(input,results,config,out_dir,label)
  message("7/7 Record versions and input checksums")
  write_csv(versions_table(),file.path(out_dir,"software_versions.csv"))
  files <- list.files(data_dir,pattern="\\.(csv|tsv)$",full.names=TRUE)
  write_csv(data.frame(File=basename(files),MD5=unname(tools::md5sum(files))), file.path(out_dir,"input_checksums.csv"))
  writeLines(c("Run completed.",label,
    "Table S14 uses an explicit reconstructed medication contrast.",
    "See docs/METHODS_AND_PROVENANCE.md for reconstructed versus verified analytical components."), file.path(out_dir,"RUN_COMPLETED.txt"))
  message("Completed. Output folder: ",basename(out_dir))
}
