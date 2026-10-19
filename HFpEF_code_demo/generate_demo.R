#!/usr/bin/env Rscript
# Fully synthetic generator. It never reads study data or external annotations.
generate_demo <- function(out_dir, seed = 20261019L) {
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("Destination is not empty. Use a new directory; no existing data will be overwritten.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  n <- 48L; ng <- 400L; na <- 480L
  person <- sprintf("SYN_P%03d",seq_len(n))
  assay <- sprintf("SYN_A%04d",seq_len(na))
  genes <- sprintf("SYN_G%04d",seq_len(ng))
  gene_index <- c(seq_len(ng),seq_len(na-ng))
  group <- rep(c("Control","higherEF","lowerEF"),each=n/3)
  sex <- rep(rep(c("Female","Male"),each=8),3)
  age <- round(runif(n,48,82),1); bmi <- round(runif(n,21,36),1)
  z_ci <- rnorm(n); z_la <- rnorm(n)
  pheno <- data.frame(PatientID=person,Phenotype=group,Age=age,Sex=sex,BMI=bmi,
    DeltaLVSVi=round(8+4*z_ci+rnorm(n,0,3),3),
    DeltaLVEF=round(3+rnorm(n,0,5),3),
    DeltaCI=round(.8+.3*z_ci,3),DeltaLAESVi=round(2+2*z_la,3),
    DeltaLungWater=round(5+rnorm(n,0,4),3),SixMinuteWalk=round(runif(n,250,550)),
    KCCQ=round(runif(n,35,95),1),check.names=FALSE)
  pheno$DeltaLungWater[c(19,36)] <- NA_real_
  pheno$SixMinuteWalk[22] <- NA_real_
  # Independent invented subject effects and assay signals, not fitted to the study.
  shared <- matrix(rnorm(ng*n,0,.30),ng,n)
  baseline_gene <- rnorm(ng,10,1.0) + shared
  baseline_gene[81:100,group!="Control"] <- baseline_gene[81:100,group!="Control"] + .65
  baseline <- baseline_gene[gene_index,] + matrix(rnorm(na*n,0,.06),na,n)
  change_gene <- matrix(rnorm(ng*n,0,.14),ng,n)
  change_gene[1:35,] <- change_gene[1:35,] + .55
  change_gene[36:60,] <- change_gene[36:60,] - .50
  change_gene[101:120,group=="higherEF"] <- change_gene[101:120,group=="higherEF"] + .65
  change_gene[121:140,group=="lowerEF"] <- change_gene[121:140,group=="lowerEF"] - .65
  change_gene[141:155,group=="Control"] <- change_gene[141:155,group=="Control"] + .65
  change_gene[161:180,group!="Control"] <- change_gene[161:180,group!="Control"] + .55
  change_gene[201:220,] <- change_gene[201:220,] + rep(.35*z_ci,each=20)
  change_gene[241:260,] <- change_gene[241:260,] - rep(.35*z_la,each=20)
  # Correlated responses within invented pathway blocks exercise CAMERA's correlation handling.
  for (block in seq_len(20)) {
    rows <- ((block-1)*20+1):(block*20)
    change_gene[rows,] <- change_gene[rows,] + rep(rnorm(n,0,.06),each=20)
  }
  change <- change_gene[gene_index,] + matrix(rnorm(na*n,0,.035),na,n)
  baseline[1,] <- 10 + .015*(age-65)+rnorm(n,0,.15)
  change[1,] <- rnorm(n,0,.12)
  rest <- 2^baseline; stress <- 2^(baseline+change)
  dimnames(rest) <- dimnames(stress) <- list(assay,person)
  samples <- data.frame(SampleID=c(paste0(person,"_REST"),paste0(person,"_STRESS")),
    PatientID=rep(person,2),Condition=rep(c("Rest","Stress"),each=n))
  combined <- cbind(rest,stress); colnames(combined) <- samples$SampleID
  features <- data.frame(SeqId=assay,GeneSymbol=genes[gene_index],
    TargetName=paste("Synthetic target",gene_index),UniProtID=sprintf("SYN_U%04d",gene_index),Role="protein")
  features$Role[1] <- "cystatin_covariate"
  # Include multiple UniProt mappings to exercise the corrected SeqId handling.
  multi <- c(12,27,81,126,205,246,303,407)
  alt <- sprintf("SYN_ALT%02d",seq_along(multi))
  features$UniProtID[multi] <- paste(features$UniProtID[multi],alt,sep="|")
  tissues <- c("skeletal muscle","heart muscle","liver","lung","lymphoid tissue","bone marrow","brain")
  make_tissue <- function(i) {
    a <- tissues[(i-1)%%length(tissues)+1]; b <- tissues[i%%length(tissues)+1]
    if (i%%11 == 0) paste0(a,": 30;",b,": 20")
    else if (i%%7 == 0) paste0(a,": 60")
    else paste0(a,": 80;",b,": 8")
  }
  atlas <- data.frame(Uniprot=c(sprintf("SYN_U%04d",seq_len(ng)),alt),
    RNA.tissue.specific.nTPM=vapply(seq_len(ng+length(alt)),make_tissue,character(1)),
    Reference="Fictional annotation for synthetic demonstration",check.names=FALSE)
  pathways <- do.call(rbind,lapply(seq_len(30),function(k) {
    selected <- if (k<=20) ((k-1)*20+1):(k*20) else ((k-21)*25+seq_len(35)-1)%%ng+1
    data.frame(Pathway=sprintf("DEMO_PATHWAY_%02d",k),GeneSymbol=genes[selected],
      Family=paste("Demo family",LETTERS[(k-1)%%6+1]),Level=ifelse(k%%2,2,3),
      Reference="Synthetic pathway membership v1; not Reactome")
  }))
  # A separate random stream keeps medication generation independent of
  # expression, clinical and reference data generation.
  set.seed(seed + 14L)
  drug_names <- c("ACEi","ARB","MRA","BetaBlocker","SGLT2i","LoopDiuretic")
  exposed <- c(6L,8L,5L,10L,4L,8L)
  for (j in seq_along(drug_names)) {
    pheno[[drug_names[j]]] <- unlist(lapply(seq_len(3),function(g)
      sample(c(rep(1L,exposed[j]),rep(0L,n/3-exposed[j])))))
  }
  write_matrix <- function(x,name) write.csv(data.frame(SeqId=rownames(x),x,check.names=FALSE),file.path(out_dir,name),row.names=FALSE)
  write_matrix(combined,"expressiondata_restandstress.csv")
  write_matrix(rest,"expressiondata-rest.csv"); write_matrix(stress,"expressiondata-stress.csv")
  write.csv(pheno,file.path(out_dir,"phenotypedata.csv"),row.names=FALSE,na="")
  write.csv(samples,file.path(out_dir,"sampledata.csv"),row.names=FALSE)
  write.csv(features,file.path(out_dir,"featuredata.csv"),row.names=FALSE)
  write.table(atlas,file.path(out_dir,"proteinatlas.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
  write.table(pathways,file.path(out_dir,"pathways.tsv"),sep="\t",quote=TRUE,row.names=FALSE)
  writeLines(c("SYNTHETIC DATA ONLY",paste("Seed:",seed),
    "48 fictional participants; 96 paired samples; 480 fictional assays mapping to 400 fictional genes.",
    "Created from independent random numbers and explicitly invented effects.",
    "No study matrix, patient values, patient identifiers, or real annotation files were read.",
    "Tissue annotations and pathway memberships are fictional; they have no biological meaning.",
    "Six binary medication fields are independently invented; 1=use, 0=no use. No medication effect was inserted."),
    file.path(out_dir,"SYNTHETIC_DATA.txt"))
  invisible(out_dir)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  script <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
  root <- dirname(normalizePath(script))
  i <- match("--out",args)
  if (!is.na(i) && i==length(args)) stop("Missing --out directory")
  destination <- if (is.na(i)) file.path(root,"data","synthetic") else args[i+1L]
  generate_demo(destination)
  message("Synthetic data generated.")
}
