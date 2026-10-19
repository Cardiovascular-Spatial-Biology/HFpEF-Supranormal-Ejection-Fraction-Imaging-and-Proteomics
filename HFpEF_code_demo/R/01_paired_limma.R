# Fig 4b analytical core: retain the original group-wise model and moderation.
run_paired_limma <- function(input, cfg, out_dir) {
  results <- list(); diagnostics <- list()
  for (group in cfg$groups) {
    p <- input$pheno[input$pheno$Phenotype == group, , drop = FALSE]
    ids <- p$PatientID
    y <- cbind(log2(input$rest[, ids, drop = FALSE]), log2(input$stress[, ids, drop = FALSE]))
    metadata <- rbind(p, p)
    metadata$Condition <- factor(rep(c("rest", "stress"), each = nrow(p)), levels = c("rest", "stress"))
    design <- stats::model.matrix(~ 0 + Condition + Age + Sex + BMI + BaselineCystatin, metadata)
    check_design(design, paste("Paired", group))
    contrast <- limma::makeContrasts(Conditionstress - Conditionrest, levels = design)
    corfit <- limma::duplicateCorrelation(y, design, block = metadata$PatientID)
    assert(is.finite(corfit$consensus.correlation), paste("Cannot estimate paired correlation for", group))
    fit <- limma::lmFit(y, design, block = metadata$PatientID, correlation = corfit$consensus.correlation)
    fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
    tt <- limma::topTable(fit, number = Inf, sort.by = "none", adjust.method = "BH")
    tt <- cbind(SeqId = rownames(tt), input$features[rownames(tt), c("GeneSymbol", "TargetName"), drop = FALSE], tt)
    tt$SignificantFDR <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < cfg$fdr
    tt$SignificantFDRandFC <- tt$SignificantFDR & abs(tt$logFC) > cfg$log2fc
    results[[group]] <- tt
    write_csv(tt, file.path(out_dir, paste0("paired_", group, ".csv")))
    diagnostics[[group]] <- data.frame(Group = group, Participants = nrow(p), Samples = ncol(y),
      Assays = nrow(y), DesignRank = qr(design)$rank, ResidualDF = nrow(design) - ncol(design),
      PairedCorrelation = corfit$consensus.correlation, FDR_hits = sum(tt$SignificantFDR),
      FDR_and_FC_hits = sum(tt$SignificantFDRandFC))
  }
  write_csv(do.call(rbind, diagnostics), file.path(out_dir, "paired_model_diagnostics.csv"))
  results
}
