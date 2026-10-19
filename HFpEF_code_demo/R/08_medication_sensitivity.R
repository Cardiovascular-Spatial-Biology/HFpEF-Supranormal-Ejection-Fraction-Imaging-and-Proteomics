# Table S14 reconstructed from the table caption and written methods.
# Centering drug exposure defines the exercise effect at observed drug prevalence.
fit_medication_response <- function(input,patient_ids,drug=NULL) {
  p <- input$pheno[patient_ids,,drop=FALSE]
  y <- cbind(log2(input$rest[,patient_ids,drop=FALSE]),log2(input$stress[,patient_ids,drop=FALSE]))
  metadata <- rbind(p,p)
  metadata$Condition <- factor(rep(c("Rest","Stress"),each=nrow(p)),levels=c("Rest","Stress"))
  prevalence <- if(is.null(drug)) NA_real_ else mean(p[[drug]])
  if(is.null(drug)) design <- stats::model.matrix(~Condition+Age+Sex+BMI+BaselineCystatin,metadata)
  else {
    metadata$DrugCentered <- metadata[[drug]]-prevalence
    design <- stats::model.matrix(~Condition*DrugCentered+Age+Sex+BMI+BaselineCystatin,metadata)
  }
  check_design(design,if(is.null(drug)) "Medication baseline" else paste("Medication",drug))
  correlation <- limma::duplicateCorrelation(y,design,block=metadata$PatientID)$consensus.correlation
  assert(is.finite(correlation),"Medication model has non-finite paired correlation")
  fit <- limma::eBayes(limma::lmFit(y,design,block=metadata$PatientID,correlation=correlation))
  table <- function(coef) {
    d <- limma::topTable(fit,coef=coef,number=Inf,sort.by="none",adjust.method="BH")
    cbind(SeqId=rownames(d),GeneSymbol=input$features[rownames(d),"GeneSymbol"],d)
  }
  list(response=table("ConditionStress"),
    interaction=if(is.null(drug)) NULL else table("ConditionStress:DrugCentered"),
    fit=fit,design=design,prevalence=prevalence,correlation=correlation)
}

run_medication_sensitivity <- function(input,primary_results,cfg,out_dir,label) {
  summaries <- list(); all_response <- list(); all_interaction <- list(); details <- list()
  for(group in c("lowerEF","higherEF")) {
    p <- input$pheno[input$pheno$Phenotype==group,,drop=FALSE]
    primary <- primary_results[[group]]
    original_sig <- primary$SeqId[primary$SignificantFDRandFC]
    for(drug in cfg$medication_columns) {
      complete <- !is.na(p[[drug]]); cohort <- p[complete,,drop=FALSE]
      n <- nrow(cohort); users <- sum(cohort[[drug]]==1)
      summary <- data.frame(Group=group,Drug=drug,DrugLabel=unname(cfg$medication_labels[drug]),
        Participants=n,MissingDrug=sum(!complete),Users=users,NonUsers=n-users,
        OriginalPrimaryHits=length(original_sig),BaselineSameCohortHits=NA_integer_,
        RetainedFDR=NA_integer_,RetainedFDRandFC=NA_integer_,RetainedSameDirection=NA_integer_,
        NewlySignificantFDRandFC=NA_integer_,InteractionFDRHits=NA_integer_,
        DrugPrevalence=if(n) users/n else NA_real_,Status="not estimated",stringsAsFactors=FALSE)
      if(users<2 || n-users<2) {
        summary$Status <- "not estimable: fewer than two participants in a medication stratum"
        summaries[[length(summaries)+1L]] <- summary
        next
      }
      model <- tryCatch(list(base=fit_medication_response(input,cohort$PatientID),
        adjusted=fit_medication_response(input,cohort$PatientID,drug)),error=function(e)e)
      if(inherits(model,"error")) {
        summary$Status <- paste("not estimable:",conditionMessage(model))
        summaries[[length(summaries)+1L]] <- summary
        next
      }
      base <- model$base$response; adjusted <- model$adjusted$response
      interaction <- model$adjusted$interaction
      hit <- !is.na(base$adj.P.Val) & base$adj.P.Val<cfg$fdr & abs(base$logFC)>cfg$log2fc
      after_fdr <- !is.na(adjusted$adj.P.Val) & adjusted$adj.P.Val<cfg$fdr
      after_both <- after_fdr & abs(adjusted$logFC)>cfg$log2fc
      same_direction <- sign(base$logFC)==sign(adjusted$logFC)
      summary$BaselineSameCohortHits <- sum(hit)
      summary$RetainedFDR <- sum(hit & after_fdr)
      summary$RetainedFDRandFC <- sum(hit & after_both)
      summary$RetainedSameDirection <- sum(hit & after_both & same_direction)
      summary$NewlySignificantFDRandFC <- sum(!hit & after_both)
      summary$InteractionFDRHits <- sum(interaction$adj.P.Val<cfg$fdr,na.rm=TRUE)
      summary$Status <- "estimated"
      adjusted$Group <- interaction$Group <- group
      adjusted$Drug <- interaction$Drug <- drug
      adjusted$Participants <- interaction$Participants <- n
      adjusted$DrugPrevalence <- model$adjusted$prevalence
      adjusted$BaselinePValue <- base$P.Value; adjusted$BaselineFDR <- base$adj.P.Val
      adjusted$BaselineLogFC <- base$logFC; adjusted$InBaselineSet <- hit
      adjusted$RetainedFDR <- hit & after_fdr
      adjusted$RetainedFDRandFC <- hit & after_both
      adjusted$RetainedSameDirection <- hit & after_both & same_direction
      all_response[[length(all_response)+1L]] <- adjusted
      all_interaction[[length(all_interaction)+1L]] <- interaction
      details[[length(details)+1L]] <- data.frame(Group=group,Drug=drug,
        Correlation=model$adjusted$correlation,DesignRank=qr(model$adjusted$design)$rank,
        ResidualDF=nrow(model$adjusted$design)-ncol(model$adjusted$design))
      summaries[[length(summaries)+1L]] <- summary
    }
  }
  summary <- do.call(rbind,summaries)
  write_csv(summary,file.path(out_dir,"TableS14_medication_summary.csv"))
  if(length(all_response)) write_csv(do.call(rbind,all_response),file.path(out_dir,"TableS14_adjusted_response_all_assays.csv"))
  if(length(all_interaction)) {
    interactions <- do.call(rbind,all_interaction)
    interactions$FDR_global_exploratory <- stats::p.adjust(interactions$P.Value,"BH")
    write_csv(interactions,file.path(out_dir,"TableS14_interactions_all_assays.csv"))
    write_csv(do.call(rbind,details),file.path(out_dir,"TableS14_model_diagnostics.csv"))
  }
  rows <- c("Primary model",unname(cfg$medication_labels[cfg$medication_columns]))
  table <- data.frame(DrugClass=rows,stringsAsFactors=FALSE)
  for(group in c("lowerEF","higherEF")) {
    z <- summary[summary$Group==group,,drop=FALSE]; z <- z[match(cfg$medication_columns,z$Drug),]
    cells <- ifelse(z$Status=="estimated",paste0(z$RetainedFDR," / ",z$BaselineSameCohortHits),"Not estimable")
    table[[group]] <- c(as.character(sum(primary_results[[group]]$SignificantFDRandFC)),cells)
  }
  write_csv(table,file.path(out_dir,"TableS14_display.csv"))
  writeLines(c("Table S14-style medication sensitivity (reconstructed)",plot_caption(label),
    "Primary set: FDR <0.05 and |log2FC| >log2(1.25). Main table retention: FDR <0.05 after adjustment.",
    "Adjustment: each drug separately as a main effect and drug-by-exercise interaction, plus Age, Sex, BMI and baseline cystatin; participant block.",
    "Drug exposure is centred at its observed prevalence in the analysed group, so the exercise coefficient is a prevalence-weighted average response.",
    "The original final contrast/centering was not available; this is an explicit reconstruction choice.",
    "Detailed outputs additionally report retention using both FDR and fold change, direction, new hits and interactions.",
    "Missing medication data: use matched complete cases for both baseline and adjusted models; denominator may differ from the full primary model.",
    "BH correction: all assays within each group/drug/contrast; pooled interaction BH is also exported as exploratory.",
    "Non-significant interactions do not prove absence of medication effects, particularly in small strata."),file.path(out_dir,"TableS14_caption.txt"))
  plotdata <- rbind(data.frame(Row=seq_along(rows),Column="Medication class",Value=table$DrugClass),
    data.frame(Row=seq_along(rows),Column="HF-EF 50-60%",Value=table$lowerEF),
    data.frame(Row=seq_along(rows),Column="HF-EF >60%",Value=table$higherEF))
  plotdata$Column <- factor(plotdata$Column,levels=c("Medication class","HF-EF 50-60%","HF-EF >60%"))
  plotdata$Row <- factor(plotdata$Row,levels=rev(seq_along(rows)))
  plot <- ggplot2::ggplot(plotdata,ggplot2::aes(Column,Row,label=Value))+
    ggplot2::geom_tile(fill="white",color="grey85")+ggplot2::geom_text(size=3)+
    ggplot2::scale_x_discrete(position="top",expand=c(0,0))+
    ggplot2::theme_minimal(base_size=11)+ggplot2::theme(panel.grid=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())+
    ggplot2::labs(title="Medication sensitivity: exercise-responsive assays retained",
      subtitle="Reconstructed; retained at FDR <0.05 / primary set on the same participants",
      x=NULL,y=NULL,caption=paste(plot_caption(label),"\nDrug-centred model; see TableS14_caption.txt for definitions"))
  save_plot(plot,file.path(out_dir,"TableS14_medication_sensitivity"),14,5.7)
  invisible(summary)
}
