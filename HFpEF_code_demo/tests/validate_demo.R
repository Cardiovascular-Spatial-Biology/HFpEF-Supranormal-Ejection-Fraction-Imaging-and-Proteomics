#!/usr/bin/env Rscript
# Scientific workflow checks using only the supplied fictional data.
file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- dirname(dirname(normalizePath(file_arg)))
source(file.path(root,"config.R"))
for (f in sort(list.files(file.path(root,"R"),pattern="\\.R$",full.names=TRUE))) source(f)
source(file.path(root,"generate_demo.R"))
require_packages()
temporary <- tempfile("hfpef_validation_"); dir.create(temporary)
checks <- character()
check <- function(condition, name) {
  assert(condition, paste("FAIL:",name)); checks <<- c(checks,name); message("PASS: ",name)
}
data <- file.path(root,"data","synthetic")
regenerated <- file.path(temporary,"regenerated")
generate_demo(regenerated)
names <- list.files(data,all.files=FALSE)
check(setequal(names,list.files(regenerated)) &&
      identical(unname(tools::md5sum(file.path(data,names))),unname(tools::md5sum(file.path(regenerated,names)))),
      "Fixed-seed generator reproduces every synthetic input byte")
input <- load_inputs(data,config)
check(nrow(input$x)==480 && ncol(input$x)==96 && nrow(input$pheno)==48,
      "Synthetic dimensions and all sample/feature mappings")
shuffled <- file.path(temporary,"shuffled"); dir.create(shuffled)
invisible(file.copy(list.files(data,full.names=TRUE),shuffled))
set.seed(37)
for (file in c("phenotypedata.csv","sampledata.csv","featuredata.csv")) {
  d <- read_table(file.path(shuffled,file))
  write_csv(d[sample(nrow(d)),,drop=FALSE],file.path(shuffled,file))
}
for (file in c("expressiondata_restandstress.csv","expressiondata-rest.csv","expressiondata-stress.csv")) {
  d <- read_table(file.path(shuffled,file))
  d <- d[sample(nrow(d)),c(1,sample(2:ncol(d))),drop=FALSE]
  write_csv(d,file.path(shuffled,file))
}
permuted <- load_inputs(shuffled,config)
check(isTRUE(all.equal(unname(input$rest),unname(permuted$rest[rownames(input$rest),colnames(input$rest)]),tolerance=1e-12)),
      "Paired values are invariant to all input row/column permutations")
original_out <- file.path(temporary,"original"); dir.create(original_out)
permuted_out <- file.path(temporary,"permuted"); dir.create(permuted_out)
a <- run_paired_limma(input,config,original_out)
b <- run_paired_limma(permuted,config,permuted_out)
for (group in config$groups) {
  x <- a[[group]]; y <- b[[group]][match(x$SeqId,b[[group]]$SeqId),]
  check(isTRUE(all.equal(unname(as.matrix(x[c("logFC","t","adj.P.Val")])),
                       unname(as.matrix(y[c("logFC","t","adj.P.Val")])),tolerance=1e-8)),
        paste("ID matching preserves paired limma estimates:",group))
}
bad <- read_table(file.path(shuffled,"phenotypedata.csv")); bad$PatientID[2] <- bad$PatientID[1]
write_csv(bad,file.path(shuffled,"phenotypedata.csv"))
error <- tryCatch({load_inputs(shuffled,config); ""},error=function(e) conditionMessage(e))
check(grepl("PatientID must be nonempty and unique",error,fixed=TRUE),"Duplicate participants fail explicitly")
expanded <- expand_uniprot(input$features)
check(nrow(expanded)==488 && all(expanded$SeqId %in% input$features$SeqId),
      "Multi-UniProt expansion preserves assay IDs")
selected <- select_gene_assays(input)
means <- rowMeans(input$rest[,input$pheno$PatientID[input$pheno$Phenotype!="Control"]])
check(nrow(selected)==400 && all(vapply(seq_len(nrow(selected)),function(i) {
  candidates <- input$features$SeqId[input$features$GeneSymbol==selected$GeneSymbol[i]]
  means[selected$SeqId[i]] == max(means[candidates])
},logical(1))),"CAMERA selects one highest-mean-resting assay per gene")
x <- t(log2(input$x+1))
groups <- factor(input$pheno[input$samples$PatientID,"Phenotype"])
fit <- mixOmics::plsda(x,groups,ncomp=2,scale=TRUE)
projection <- vapply(1:2,function(j) {
  score <- fit$variates$X[,j]
  sum(crossprod(score,fit$X)^2)/sum(score^2)/sum(fit$X^2)
},numeric(1))
check(isTRUE(all.equal(unname(projection),unname(fit$prop_expl_var$X),tolerance=1e-10)),
      "PLS-DA axis percentages agree with an independent X-projection calculation")
camera_file <- file.path(root,"results","demo","Fig4g_camera_all_results.csv")
if (file.exists(camera_file)) {
  camera <- read_table(camera_file)
  check(all(vapply(split(camera,camera$Variable),function(d)
    isTRUE(all.equal(d$FDR,stats::p.adjust(d$PValue,"BH"),tolerance=1e-12)),logical(1))),
    "CAMERA controls BH FDR separately within each variable")
}
# Fig 4d: independent small-example probability check and explicit universe.
background <- sprintf("TEST_G%02d",1:20)
membership <- data.frame(Pathway=c(rep("A",5),rep("B",4)),
  GeneSymbol=c(background[c(1,2,5,6,7)],background[10:13]))
ora <- local_overrepresentation(background[1:4],background,membership)
fisher <- stats::fisher.test(matrix(c(2,2,3,13),nrow=2,byrow=TRUE),alternative="greater")$p.value
check(abs(ora$PValue[ora$Term=="A"]-fisher)<1e-12,"Offline enrichment agrees with independent Fisher exact calculation")
check(ora$PValue[ora$Term=="B"]==1 &&
  isTRUE(all.equal(ora$FDR,stats::p.adjust(ora$PValue,"BH"))),"BH includes tested zero-overlap pathways")
empty <- local_overrepresentation(character(),background,membership)
check(all(empty$PValue==1) && all(empty$FDR==1),"Empty enrichment query gives no enrichment")
extra <- rbind(membership,data.frame(Pathway="A",GeneSymbol="UNMEASURED"))
check(isTRUE(all.equal(ora,local_overrepresentation(background[1:4],background,extra))),
  "Unmeasured reference genes do not inflate the background or set size")
gene_lists <- enrichment_gene_lists(input,a,config)
check(length(gene_lists$background)==400 && !anyDuplicated(gene_lists$background),
  "Enrichment counts unique measured genes, not repeated assays")
fixture <- data.frame(Term="TEST_SET",P.value=.001,Adjusted.P.value=.02,Overlap="2/5",Genes="TEST_G01;TEST_G02")
normalised <- normalise_enrichr(fixture)
check(normalised$FDR==.02 && normalised$Count==2 && normalised$PathwaySize==5,
  "Enrichr adapter retains supplied adjusted p-values and overlaps")
guard_cfg <- config; guard_cfg$enrichment_mode <- "enrichr_live"; guard_cfg$allow_enrichr_network <- FALSE
guard_out <- file.path(temporary,"network_guard"); dir.create(guard_out)
error <- tryCatch({run_pathway_enrichment(input,a,guard_cfg,guard_out,"SYNTHETIC DEMONSTRATION"); ""},error=function(e) conditionMessage(e))
check(grepl("allow-enrichr-network",error,fixed=TRUE),"Live Enrichr cannot run without explicit network selection")

# Table S14: preserve primary model and independently verify the estimand.
ids <- input$pheno$PatientID[input$pheno$Phenotype=="higherEF"]
baseline <- fit_medication_response(input,ids)
reference <- a$higherEF[match(baseline$response$SeqId,a$higherEF$SeqId),]
check(isTRUE(all.equal(unname(as.matrix(baseline$response[c("logFC","t","adj.P.Val")])),
  unname(as.matrix(reference[c("logFC","t","adj.P.Val")])),tolerance=1e-8)),
  "Medication baseline reproduces the existing paired primary model")
model <- fit_medication_response(input,ids,"BetaBlocker")
p <- input$pheno[ids,,drop=FALSE]; md <- rbind(p,p)
md$Condition <- factor(rep(c("Rest","Stress"),each=nrow(p)),levels=c("Rest","Stress"))
md$Drug <- md$BetaBlocker
design <- stats::model.matrix(~Condition*Drug+Age+Sex+BMI+BaselineCystatin,md)
y <- cbind(log2(input$rest[,ids]),log2(input$stress[,ids]))
fit <- limma::lmFit(y,design,block=md$PatientID,correlation=model$correlation)
contrast <- setNames(rep(0,ncol(design)),colnames(design))
contrast["ConditionStress"] <- 1
contrast["ConditionStress:Drug"] <- mean(p$BetaBlocker)
weighted <- limma::eBayes(limma::contrasts.fit(fit,contrast))
weighted <- limma::topTable(weighted,coef=1,number=Inf,sort.by="none")
check(max(abs(weighted$logFC-model$response$logFC))<1e-9 && max(abs(weighted$t-model$response$t))<1e-8,
  "Centered exercise effect equals the prevalence-weighted user/non-user contrast")
reverse <- input; reverse$pheno$BetaBlocker <- 1-reverse$pheno$BetaBlocker
reversed_model <- fit_medication_response(reverse,ids,"BetaBlocker")
check(max(abs(reversed_model$response$logFC-model$response$logFC))<1e-9 &&
      max(abs(reversed_model$interaction$logFC+model$interaction$logFC))<1e-9,
  "Reversing medication coding preserves average response and reverses interaction sign")
check(isTRUE(all.equal(model$interaction$adj.P.Val,stats::p.adjust(model$interaction$P.Value,"BH"),tolerance=1e-12)),
  "Interaction FDR uses all tested assays for the group and drug")
med_out <- file.path(temporary,"medication_edge"); dir.create(med_out)
edge <- input; edge$pheno$ACEi[edge$pheno$Phenotype=="higherEF"] <- 0
edge_cfg <- config; edge_cfg$medication_columns <- "ACEi"
edge_summary <- run_medication_sensitivity(edge,a,edge_cfg,med_out,"SYNTHETIC DEMONSTRATION")
check(grepl("not estimable",edge_summary$Status[edge_summary$Group=="higherEF"],fixed=TRUE),
  "A constant medication exposure is reported as not estimable, not as a null interaction")
message(length(checks)," checks passed.")
unlink(temporary,recursive=TRUE)
