
run_overlap_ternary <- function(results, cfg, out_dir, label) {
  ids <- results$Control$SeqId
  merged <- results$Control[, c("SeqId", "GeneSymbol"), drop = FALSE]
  for (group in cfg$groups) {
    d <- results[[group]][match(ids, results[[group]]$SeqId), , drop = FALSE]
    assert(identical(d$SeqId, ids), "Top-table feature IDs do not match")
    merged[[paste0("logFC_", group)]] <- d$logFC
    merged[[paste0("FDR_", group)]] <- d$adj.P.Val
    merged[[paste0("sig_", group)]] <- d$SignificantFDRandFC
  }
  merged$Category <- category_from_flags(merged$sig_Control, merged$sig_higherEF, merged$sig_lowerEF)
  counts <- as.data.frame(table(factor(merged$Category, levels = names(category_colors))))
  names(counts) <- c("Category", "Assays")
  write_csv(counts, file.path(out_dir, "Fig4b_overlap_counts.csv"))
  sets <- lapply(results, function(d) d$SeqId[d$SignificantFDRandFC])
  venn <- VennDiagram::venn.diagram(sets[cfg$groups], filename = NULL,
    category.names = unname(group_labels[cfg$groups]), fill = unname(group_colors[cfg$groups]),
    alpha = 0.5, cex = 1.3, cat.cex = 1.05, margin = 0.12, main = "Exercise-responsive assays",
    sub = plot_caption(label), sub.cex = 0.65, disable.logging = TRUE)
  for (extension in c("pdf", "png")) {
    file <- file.path(out_dir, paste0("Fig4b_overlap.", extension))
    if (extension == "pdf") grDevices::pdf(file, width = 8, height = 7)
    else grDevices::png(file, width = 1280, height = 1120, res = 160)
    grid::grid.newpage(); grid::grid.draw(venn); grDevices::dev.off()
  }
  effects <- abs(as.matrix(merged[paste0("logFC_", cfg$groups)]))
  bary <- (effects + 0.01) / rowSums(effects + 0.01)
  colnames(bary) <- paste0("fraction_", cfg$groups)
  merged <- cbind(merged, bary)
  # Exact barycentric-to-Cartesian conversion avoids ggtern label-layer incompatibility.
  merged$x <- merged$fraction_lowerEF + 0.5 * merged$fraction_higherEF
  merged$y <- sqrt(3) / 2 * merged$fraction_higherEF
  merged$max_effect <- apply(effects, 1, max)
  merged$distance <- sqrt(rowSums((bary - 1/3)^2))
  merged$Category <- factor(merged$Category, levels = names(category_colors))
  merged <- merged[order(merged$Category), ]
  sig <- merged[merged$Category != "Not significant", ]
  labels <- sig[order(-sig$distance, sig$SeqId), ]
  labels <- head(labels[!duplicated(labels$GeneSymbol), ], 9)
  triangle <- data.frame(x = c(0, .5, 1, 0), y = c(0, sqrt(3)/2, 0, 0))
  p <- ggplot2::ggplot(merged, ggplot2::aes(x, y)) +
    ggplot2::geom_path(data = triangle, color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(fill = Category, size = max_effect), shape = 21, color = "grey35", stroke = .15, alpha = .8, show.legend = TRUE) +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = GeneSymbol), seed = cfg$seed,
      size = 2.5, max.overlaps = Inf, ylim = c(.035,.82), box.padding = .4, point.padding = .2,
      min.segment.length = 0, max.iter = 10000, max.time = 2, show.legend = FALSE) +
    ggplot2::annotate("text", x = c(0, .5, 1), y = c(-.035, sqrt(3)/2+.04, -.035),
      label = c("Control", "HF-EF >60%", "HF-EF 50-60%"), size = 3.3) +
    ggplot2::scale_fill_manual(values = category_colors, drop = FALSE,
      guide = ggplot2::guide_legend(override.aes = list(size = 2.8, alpha = 1))) +
    ggplot2::scale_size_continuous(range = c(.5, 4), name = "Maximum |log2FC|") +
    ggplot2::coord_equal(xlim = c(-.12,1.12), ylim = c(-.1,1), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = "Relative exercise response across groups", subtitle = "Position uses |log2FC| + 0.01 in each group", caption = plot_caption(label))
  save_plot(p, file.path(out_dir, "Fig4c_ternary"), 10, 7)
  write_csv(merged, file.path(out_dir, "Fig4c_ternary_data.csv"))
  invisible(merged)
}
