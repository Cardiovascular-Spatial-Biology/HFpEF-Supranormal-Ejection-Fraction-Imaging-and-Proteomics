# 
run_resting_comparisons <- function(input, cfg, out_dir, label) {
  for (comparison in c("HFpEF_vs_Control", "higherEF_vs_lowerEF")) {
    p <- input$pheno
    if (comparison == "HFpEF_vs_Control") {
      p$Group <- factor(ifelse(p$Phenotype == "Control", "Control", "HFpEF"), levels = c("Control", "HFpEF"))
    } else {
      p <- p[p$Phenotype != "Control", , drop = FALSE]
      p$Group <- factor(p$Phenotype, levels = c("lowerEF", "higherEF"))
    }
    design <- stats::model.matrix(~ Group, p); check_design(design, comparison)
    fit <- limma::eBayes(limma::lmFit(log2(input$rest[, p$PatientID, drop = FALSE]), design))
    d <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none", adjust.method = "BH")
    d <- cbind(SeqId = rownames(d), GeneSymbol = input$features[rownames(d), "GeneSymbol"], d)
    d$Significant <- !is.na(d$adj.P.Val) & d$adj.P.Val < cfg$fdr
    d$neglogFDR <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    labels <- head(d[order(d$adj.P.Val), ], 8)
    p <- ggplot2::ggplot(d, ggplot2::aes(logFC, neglogFDR, color = Significant)) +
      ggplot2::geom_point(alpha = .7, size = 1.6) +
      ggplot2::geom_hline(yintercept = -log10(cfg$fdr), linetype = 2, color = "grey50") +
      ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = GeneSymbol), seed = cfg$seed, color = "black", size = 2.6) +
      ggplot2::scale_color_manual(values = c("FALSE" = "#bdbdbd", "TRUE" = "#2f7fa8")) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::labs(x = "log2 fold change", y = "-log10(BH-adjusted p-value)",
        title = paste("Rest:", gsub("_", " ", comparison)), subtitle = "Unadjusted limma; reconstructed implementation", caption = plot_caption(label))
    panel <- if (comparison == "HFpEF_vs_Control") "FigS4a" else "FigS4b"
    save_plot(p, file.path(out_dir, paste0(panel, "_", comparison)))
    write_csv(d, file.path(out_dir, paste0(panel, "_", comparison, ".csv")))
  }
}
