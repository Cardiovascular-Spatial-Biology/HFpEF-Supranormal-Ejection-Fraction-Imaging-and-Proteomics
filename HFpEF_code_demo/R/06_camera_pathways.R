# Reconstructed from the final written methods plus the available CAMERA core.
# Reference pathway membership/families must be supplied explicitly and versioned.
select_gene_assays <- function(input) {
  hf <- input$pheno$PatientID[input$pheno$Phenotype != "Control"]
  f <- input$features
  f$MeanRestRFU <- rowMeans(input$rest[, hf, drop = FALSE])
  # Ambiguous multi-gene annotations cannot be treated as a single gene.
  keep <- !is.na(f$GeneSymbol) & nzchar(f$GeneSymbol) & !grepl("[|;,[:space:]]", f$GeneSymbol)
  f <- f[keep, , drop = FALSE]
  f <- f[order(f$GeneSymbol, -f$MeanRestRFU, f$SeqId), , drop = FALSE]
  f[!duplicated(f$GeneSymbol), , drop = FALSE]
}

run_camera <- function(input, cfg, out_dir, label) {
  annotation <- read_table(file.path(input$data_dir, "pathways.tsv"), sep = "\t")
  check_columns(annotation, c("Pathway", "GeneSymbol", "Family", "Level", "Reference"), "pathways.tsv")
  assert(!anyNA(annotation), "Pathway annotations must not contain missing values")
  assert(all(nzchar(annotation$Reference)), "Each pathway needs reference/version provenance")
  for (p in split(annotation, annotation$Pathway)) {
    assert(length(unique(p$Family)) == 1 && length(unique(p$Level)) == 1,
           "Each pathway must have one explicit family and hierarchy level")
  }
  annotation <- annotation[annotation$Level %in% cfg$camera_levels, , drop = FALSE]
  selected <- select_gene_assays(input)
  assert(nrow(selected) > cfg$camera_set_min, "Too few unambiguous gene annotations for CAMERA")
  write_csv(selected[, c("SeqId", "GeneSymbol", "MeanRestRFU")], file.path(out_dir, "camera_selected_assays.csv"))
  p <- input$pheno[input$pheno$Phenotype != "Control", , drop = FALSE]
  delta <- log2(input$stress[selected$SeqId, p$PatientID, drop = FALSE]) -
           log2(input$rest[selected$SeqId, p$PatientID, drop = FALSE])
  rownames(delta) <- selected$GeneSymbol
  sets <- split(annotation$GeneSymbol, annotation$Pathway)
  sets <- lapply(sets, unique)
  metadata <- annotation[!duplicated(annotation$Pathway), c("Pathway", "Family", "Level", "Reference")]
  all_results <- list(); diagnostics <- list()
  for (v in cfg$camera_variables) {
    complete <- stats::complete.cases(p[, c(v, "Age", "Sex", "BMI")])
    d <- p[complete, , drop = FALSE]; d$x <- d[[v]]
    design <- stats::model.matrix(~ x + Age + Sex + BMI, d)
    check_design(design, paste("CAMERA", v))
    y <- delta[, complete, drop = FALSE]
    variable <- apply(y, 1, stats::var) > .Machine$double.eps
    y <- y[variable, , drop = FALSE]
    index <- lapply(sets, function(g) which(rownames(y) %in% g))
    index <- index[lengths(index) >= cfg$camera_set_min & lengths(index) <= cfg$camera_set_max]
    assert(length(index) > 0, paste("No eligible pathways for", v))
    assert(all(lengths(index) < nrow(y)), "A pathway cannot cover the entire measured gene universe")
    result <- limma::camera(y, index, design, contrast = "x", inter.gene.cor = cfg$camera_inter_gene_cor, sort = TRUE)
    result$Pathway <- rownames(result)
    result$Variable <- v; result$Participants <- nrow(d)
    result$Family <- metadata$Family[match(result$Pathway, metadata$Pathway)]
    result$Reference <- metadata$Reference[match(result$Pathway, metadata$Pathway)]
    result$CorrelationSetting <- if (is.na(cfg$camera_inter_gene_cor)) "estimated" else as.character(cfg$camera_inter_gene_cor)
    result$Rank <- seq_len(nrow(result))
    all_results[[v]] <- result
    diagnostics[[v]] <- data.frame(Variable = v, Participants = nrow(d), Genes = nrow(y),
      Pathways = length(index), FDR_hits = sum(result$FDR < cfg$fdr),
      ResidualDF = nrow(design) - ncol(design), CorrelationSetting = result$CorrelationSetting[1])
  }
  all <- do.call(rbind, all_results); rownames(all) <- NULL
  all$FDR_global_exploratory <- stats::p.adjust(all$PValue, method = "BH")
  write_csv(all, file.path(out_dir, "Fig4g_camera_all_results.csv"))
  write_csv(do.call(rbind, diagnostics), file.path(out_dir, "camera_model_diagnostics.csv"))
  chosen <- lapply(all_results, function(d) {
    if (any(d$FDR < cfg$fdr)) d[d$FDR < cfg$fdr, , drop = FALSE]
    else d[d$Rank <= cfg$camera_top_nominal & d$PValue < .05, , drop = FALSE]
  })
  chosen <- do.call(rbind, chosen); rownames(chosen) <- NULL
  write_csv(chosen, file.path(out_dir, "Fig4g_displayed_pathways.csv"))
  if (!nrow(chosen)) {
    p <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No pathways meet the prespecified display criteria") + ggplot2::theme_void() +
      ggplot2::labs(title = "CAMERA pathway associations", caption = plot_caption(label))
    save_plot(p, file.path(out_dir, "Fig4g_camera_map"), 11, 7)
    write_csv(data.frame(Variable = character(), Family = character(), Pathways = integer(), FDR_hits = integer(),
      ProportionUp = numeric(), Direction = character()), file.path(out_dir, "Fig4g_family_links.csv"))
    return(invisible(all))
  }
  pieces <- split(chosen, interaction(chosen$Variable, chosen$Family, drop = TRUE))
  links <- do.call(rbind, lapply(pieces, function(d) {
    up <- mean(d$Direction == "Up")
    data.frame(Variable = d$Variable[1], Family = d$Family[1], Pathways = nrow(d), FDR_hits = sum(d$FDR < cfg$fdr),
      ProportionUp = up, Direction = if (up >= .75) "positive" else if (up <= .25) "negative" else "mixed")
  }))
  rownames(links) <- NULL
  write_csv(links, file.path(out_dir, "Fig4g_family_links.csv"))
  draw_camera_map(links, cfg, out_dir, label)
  invisible(all)
}

draw_camera_map <- function(links, cfg, out_dir, label) {
  left <- cfg$camera_variables[cfg$camera_variables %in% links$Variable]
  right <- sort(unique(links$Family))
  totals_l <- vapply(left, function(v) sum(links$Pathways[links$Variable == v]), numeric(1))
  totals_r <- vapply(right, function(v) sum(links$Pathways[links$Family == v]), numeric(1))
  gap <- .025; height <- .73
  scale <- min((height - gap*(length(left)-1))/sum(totals_l), (height-gap*(length(right)-1))/sum(totals_r))
  make_nodes <- function(nms, totals) {
    top <- .89 - c(0, head(cumsum(totals*scale + gap), -1))
    data.frame(Name = nms, Top = top, Bottom = top - totals*scale)
  }
  l <- make_nodes(left, totals_l); r <- make_nodes(right, totals_r)
  links <- links[order(match(links$Variable,left), match(links$Family,right)), ]
  cursor_l <- setNames(l$Top,l$Name); cursor_r <- setNames(r$Top,r$Name)
  for (i in seq_len(nrow(links))) {
    a <- links$Variable[i]; b <- links$Family[i]; w <- links$Pathways[i]*scale
    links$la[i] <- cursor_l[a]; links$lb[i] <- cursor_l[a]-w; cursor_l[a] <- cursor_l[a]-w
    links$ra[i] <- cursor_r[b]; links$rb[i] <- cursor_r[b]-w; cursor_r[b] <- cursor_r[b]-w
  }
  draw <- function() {
    graphics::par(mar = c(0,0,0,0), family = "sans")
    graphics::plot.new(); graphics::plot.window(xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
    graphics::text(.5,.97,"Imaging / clinical variables and pathway families", font = 2, cex = 1.1)
    graphics::text(.5,.93,"Reconstructed CAMERA workflow - width counts pathways, not independent evidence", cex = .72)
    cols <- c(positive = "#b2182b", negative = "#2166ac", mixed = "#888888")
    for (i in seq_len(nrow(links))) {
      t <- seq(0,1,length.out = 80); u <- t*t*(3-2*t)
      x <- .29 + .40*t
      a <- links$la[i] + (links$ra[i]-links$la[i])*u
      b <- links$lb[i] + (links$rb[i]-links$lb[i])*u
      color <- cols[links$Direction[i]]
      graphics::polygon(c(x,rev(x)),c(a,rev(b)), border = if (links$FDR_hits[i]>0) color else NA,
        col = grDevices::adjustcolor(color, alpha.f = if (links$FDR_hits[i]>0) .75 else .2))
    }
    graphics::rect(.28,l$Bottom,.29,l$Top,col="#34657f",border=NA)
    graphics::rect(.69,r$Bottom,.70,r$Top,col="#555555",border=NA)
    graphics::text(.27,(l$Top+l$Bottom)/2,l$Name,adj=1,cex=.8)
    graphics::text(.71,(r$Top+r$Bottom)/2,r$Name,adj=0,cex=.75)
    graphics::legend("bottom", inset = .075, legend = c("Positive", "Negative", "Mixed"),
      fill = cols, horiz = TRUE, bty = "n", cex = .75)
    graphics::text(.5,.07,"Solid: at least one pathway FDR <0.05; pale: nominal p <0.05 only",cex=.72)
    graphics::text(.5,.035,plot_caption(label),cex=.68)
  }
  for (extension in c("pdf", "png")) {
    file <- file.path(out_dir, paste0("Fig4g_camera_map.",extension))
    if (extension == "pdf") grDevices::pdf(file,width=12,height=8)
    else grDevices::png(file,width=1920,height=1280,res=160)
    draw(); grDevices::dev.off()
  }
}
