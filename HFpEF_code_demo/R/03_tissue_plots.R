# Retains the available script's rule on the LISTED tissue nTPM values.
# This summary field alone does not verify a comparison with all other tissues.
parse_tissue <- function(value, fold) {
  if (is.na(value) || !nzchar(value)) return(NA_character_)
  parts <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
  names <- trimws(sub(":.*$", "", parts))
  values <- suppressWarnings(as.numeric(trimws(sub("^.*:", "", parts))))
  assert(all(is.finite(values)) && all(values >= 0), "Malformed tissue nTPM field")
  order <- order(-values, names); names <- names[order]; values <- values[order]
  if (length(values) == 1L || (values[1] > 0 && values[1] >= fold * values[2])) names[1] else NA_character_
}

expand_uniprot <- function(features) {
  pieces <- lapply(seq_len(nrow(features)), function(i) {
    ids <- features$UniProtID[i]
    if (is.na(ids) || !nzchar(ids)) return(NULL)
    ids <- unique(trimws(strsplit(ids, "|", fixed = TRUE)[[1]]))
    # Preserve the actual assay ID, never replace it with a data-frame row number.
    data.frame(SeqId = features$SeqId[i], GeneSymbol = features$GeneSymbol[i], UniProtID = ids)
  })
  do.call(rbind, pieces)
}

run_tissue_plots <- function(input, results, cfg, out_dir, label) {
  atlas <- read_table(file.path(input$data_dir, "proteinatlas.tsv"), sep = "\t")
  check_columns(atlas, c("Uniprot", "RNA.tissue.specific.nTPM"), "proteinatlas.tsv")
  check_ids(atlas$Uniprot, "Atlas UniProt IDs")
  expanded <- expand_uniprot(input$features)
  assert(!is.null(expanded) && nrow(expanded) > 0, "No UniProt mappings supplied")
  index <- match(expanded$UniProtID, atlas$Uniprot)
  expanded$Tissue <- vapply(atlas$RNA.tissue.specific.nTPM[index], parse_tissue, character(1), fold = cfg$tissue_specificity_fold)
  write_csv(data.frame(ExpandedMappings = nrow(expanded), AtlasMatched = sum(!is.na(index)),
    AssayIDsPreserved = all(expanded$SeqId %in% input$features$SeqId)), file.path(out_dir, "tissue_mapping_diagnostics.csv"))
  mapped <- expanded[!is.na(expanded$Tissue) & expanded$Tissue %in% cfg$tissues, , drop = FALSE]
  assert(nrow(mapped) > 0, "No assays satisfy the tissue annotation rule")
  # Different UniProt entries for the same assay/tissue are one plotted assay.
  mapped <- mapped[!duplicated(mapped[c("SeqId", "Tissue")]), , drop = FALSE]
  write_csv(mapped, file.path(out_dir, "tissue_assay_mapping.csv"))
  for (group in c("lowerEF", "higherEF")) {
    d <- mapped
    current <- results[[group]][match(d$SeqId, results[[group]]$SeqId), , drop = FALSE]
    assert(!anyNA(current$SeqId), "A tissue assay is absent from limma output")
    d$logFC <- current$logFC; d$FDR <- current$adj.P.Val
    d$Direction <- ifelse(d$logFC > 0, "Increase", "Decrease")
    flags <- lapply(cfg$groups, function(g) {
      z <- results[[g]][match(d$SeqId, results[[g]]$SeqId), , drop = FALSE]
      !is.na(z$adj.P.Val) & z$adj.P.Val < cfg$fdr & sign(z$logFC) == sign(d$logFC)
    }); names(flags) <- cfg$groups
    d$Category <- category_from_flags(flags$Control, flags$higherEF, flags$lowerEF)
    # Non-significant assays in this panel stay grey, even if another group is significant.
    d$Category[is.na(d$FDR) | d$FDR >= cfg$fdr] <- "Not significant"
    d$Tissue <- factor(d$Tissue, levels = cfg$tissues)
    d <- d[order(d$Tissue, d$FDR, d$SeqId), ]
    d$x <- NA_real_
    for (i in seq_along(cfg$tissues)) {
      rows <- which(d$Tissue == cfg$tissues[i])
      d$x[rows] <- i - .5 + seq_along(rows) / (length(rows) + 1)
    }
    d$neglogFDR <- -log10(pmax(d$FDR, .Machine$double.xmin))
    labels <- d[order(d$FDR, d$SeqId), ]
    labels <- head(labels[!duplicated(labels[c("GeneSymbol", "Tissue")]), ], 12)
    p <- ggplot2::ggplot(d, ggplot2::aes(x, neglogFDR)) +
      ggplot2::geom_vline(xintercept = seq(1.5, 6.5, 1), color = "grey85", linetype = 3) +
      ggplot2::geom_hline(yintercept = -log10(cfg$fdr), color = "grey50", linetype = 2) +
      ggplot2::geom_point(ggplot2::aes(fill = Category, shape = Direction, size = abs(logFC)), color = "grey25", stroke = .2) +
      ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = GeneSymbol), seed = cfg$seed, size = 2.5, max.overlaps = Inf) +
      ggplot2::scale_fill_manual(values = category_colors, drop = FALSE,
        guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 3))) +
      ggplot2::scale_shape_manual(values = c(Increase = 24, Decrease = 25),
        guide = ggplot2::guide_legend(override.aes = list(fill = "white", size = 3))) +
      ggplot2::scale_size_continuous(range = c(.5,4), name = "|log2FC|") +
      ggplot2::scale_x_continuous(breaks = seq_along(cfg$tissues), labels = gsub(" ", "\n", cfg$tissues), limits = c(.5,7.5)) +
      ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::labs(x = "Tissue annotation", y = "-log10(BH-adjusted p-value)",
        title = paste("Exercise response:", group_labels[group]),
        subtitle = "Colours: FDR <0.05 and shared direction; no fold-change cutoff", caption = plot_caption(label))
    panel <- if (group == "lowerEF") "Fig4e" else "Fig4f"
    save_plot(p, file.path(out_dir, paste0(panel, "_tissue_", group)), 12, 6)
    write_csv(d, file.path(out_dir, paste0(panel, "_tissue_", group, ".csv")))
  }
}
