run_pca_plsda <- function(input, cfg, out_dir, label) {
  rest <- log2(input$rest + 1)
  rest <- rest[apply(rest, 1, stats::var) > 0, , drop = FALSE]
  fit <- stats::prcomp(t(rest), center = TRUE, scale. = TRUE)
  variance <- 100 * fit$sdev[1:2]^2 / sum(fit$sdev^2)
  d <- data.frame(PatientID = rownames(fit$x), Component1 = fit$x[,1], Component2 = fit$x[,2],
                  Phenotype = input$pheno[rownames(fit$x), "Phenotype"])
  p <- ggplot2::ggplot(d, ggplot2::aes(Component1, Component2, color = Phenotype)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::stat_ellipse(ggplot2::aes(group = Phenotype), type = "norm", level = .95, linetype = 2) +
    ggplot2::scale_color_manual(values = group_colors, labels = group_labels) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = "Resting PCA", x = sprintf("PC1 (%.1f%%)", variance[1]),
      y = sprintf("PC2 (%.1f%%)", variance[2]), caption = plot_caption(label))
  save_plot(p, file.path(out_dir, "FigS5a_rest_PCA"))
  write_csv(d, file.path(out_dir, "FigS5a_scores.csv"))
  var_results <- list(data.frame(Panel = "S5a", Model = "PCA", Component = 1:2, PercentXVariance = variance))
  # All assays, not a differential-expression-selected or top-variable subset.
  x <- t(log2(input$x + 1))
  assert(all(apply(x, 2, stats::var) > 0), "PLS-DA input has zero-variance assays")
  annotation <- input$samples[match(rownames(x), input$samples$SampleID), , drop = FALSE]
  annotation$Phenotype <- input$pheno[annotation$PatientID, "Phenotype"]
  annotation$Condition <- factor(annotation$Condition, levels = c("Rest", "Stress"))
  annotation$Phenotype <- factor(annotation$Phenotype, levels = c("higherEF", "lowerEF", "Control"))
  for (target in c("Condition", "Phenotype")) {
    fit <- mixOmics::plsda(x, annotation[[target]], ncomp = 2, scale = TRUE)
    variance <- 100 * fit$prop_expl_var$X[1:2]
    d <- cbind(annotation, Component1 = fit$variates$X[,1], Component2 = fit$variates$X[,2])
    if (target == "Condition") {
      p <- ggplot2::ggplot(d, ggplot2::aes(Component1, Component2, color = Condition, shape = Phenotype)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::stat_ellipse(ggplot2::aes(group = Condition, shape = NULL), type = "norm", level = .95, linetype = 2) +
        ggplot2::scale_color_manual(values = c(Rest = "#2b83ba", Stress = "#d7191c")) +
        ggplot2::scale_shape_manual(values = c(Control = 15, higherEF = 16, lowerEF = 17), labels = group_labels)
      panel <- "FigS5b"
    } else {
      p <- ggplot2::ggplot(d, ggplot2::aes(Component1, Component2, color = Phenotype, shape = Condition)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::stat_ellipse(ggplot2::aes(group = Phenotype, shape = NULL), type = "norm", level = .95, linetype = 2) +
        ggplot2::scale_color_manual(values = group_colors, labels = group_labels) +
        ggplot2::scale_shape_manual(values = c(Rest = 16, Stress = 17))
      panel <- "FigS5c"
    }
    p <- p + ggplot2::theme_classic(base_size = 11) +
      ggplot2::labs(title = paste("PLS-DA supervised by", tolower(target)),
        subtitle = "Exploratory in-sample display; not classification validation",
        x = sprintf("Component 1 (%.1f%% of X variance)", variance[1]),
        y = sprintf("Component 2 (%.1f%% of X variance)", variance[2]), caption = plot_caption(label))
    save_plot(p, file.path(out_dir, paste0(panel, "_PLSDA_", target)))
    write_csv(d, file.path(out_dir, paste0(panel, "_scores.csv")))
    var_results[[length(var_results)+1L]] <- data.frame(Panel = sub("Fig", "", panel),
      Model = paste("PLSDA", target), Component = 1:2, PercentXVariance = variance)
  }
  write_csv(do.call(rbind, var_results), file.path(out_dir, "FigS5_explained_variance.csv"))
}
