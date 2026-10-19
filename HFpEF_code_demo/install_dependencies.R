# Run manually once. Analysis scripts do not install packages or download data.
cran <- c("ggplot2","ggrepel","statmod","VennDiagram","BiocManager")
missing <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
bio <- c("limma","mixOmics")
missing <- bio[!vapply(bio,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) BiocManager::install(missing,ask=FALSE,update=FALSE)
if ("--enrichr" %in% commandArgs(trailingOnly=TRUE) && !requireNamespace("enrichR",quietly=TRUE))
  install.packages("enrichR",repos="https://cloud.r-project.org")
