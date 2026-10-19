# Fig 4d: offline fictional-set demonstration, cached Enrichr results, or
# explicitly requested live Enrichr with the measured-gene background.
# The local hypergeometric test is NOT described as an Enrichr reproduction.
enrichment_gene_lists <- function(input, results, cfg) {
  valid <- !is.na(input$features$GeneSymbol) & nzchar(input$features$GeneSymbol) &
    !grepl("[|;,[:space:]]",input$features$GeneSymbol)
  background <- sort(unique(input$features$GeneSymbol[valid]))
  queries <- lapply(results[cfg$groups],function(d)
    sort(unique(d$GeneSymbol[d$SignificantFDRandFC & d$GeneSymbol %in% background])))
  assert(length(background)>0,"No unambiguous measured genes for enrichment")
  list(background=background,queries=queries)
}

local_overrepresentation <- function(query, background, memberships) {
  background <- unique(background); query <- intersect(unique(query),background)
  sets <- lapply(split(memberships$GeneSymbol,memberships$Pathway),function(g) intersect(unique(g),background))
  sets <- sets[lengths(sets)>0]
  assert(length(sets)>0,"No pathway genes overlap the measured background")
  n <- length(query); universe <- length(background)
  d <- do.call(rbind,lapply(names(sets),function(term) {
    genes <- intersect(query,sets[[term]]); k <- length(genes); m <- length(sets[[term]])
    p <- if (k==0 || n==0) 1 else stats::phyper(k-1,m,universe-m,n,lower.tail=FALSE)
    data.frame(Term=term,PValue=p,Count=k,PathwaySize=m,QuerySize=n,BackgroundSize=universe,
      Genes=paste(sort(genes),collapse=";"),stringsAsFactors=FALSE)
  }))
  d$FDR <- stats::p.adjust(d$PValue,"BH")
  d
}

normalise_enrichr <- function(raw) {
  keys <- tolower(gsub("[^a-zA-Z0-9]","",names(raw)))
  get <- function(key) {
    i <- which(keys==key); assert(length(i)==1,paste("Missing or duplicate Enrichr field:",key)); raw[[i]]
  }
  terms <- as.character(get("term")); p <- as.numeric(get("pvalue")); fdr <- as.numeric(get("adjustedpvalue"))
  check_ids(terms,"Enrichr terms")
  assert(all(is.finite(p) & p>=0 & p<=1) && all(is.finite(fdr) & fdr>=0 & fdr<=1),"Invalid Enrichr p-values")
  genes <- if ("genes" %in% keys) as.character(get("genes")) else rep("",length(terms))
  if ("overlap" %in% keys) {
    overlap <- as.character(get("overlap"))
    assert(all(grepl("^[0-9]+/[0-9]+$",overlap)),"Enrichr Overlap must be count/set-size")
    count <- as.integer(sub("/.*","",overlap)); size <- as.integer(sub(".*/","",overlap))
  } else {
    assert("genes" %in% keys,"Enrichr export needs Overlap or Genes")
    count <- vapply(genes,function(x) if(is.na(x) || !nzchar(x)) 0L else length(unique(strsplit(x,";",fixed=TRUE)[[1]])),integer(1))
    size <- rep(NA_integer_,length(count))
  }
  data.frame(Term=terms,PValue=p,FDR=fdr,Count=count,PathwaySize=size,Genes=genes,stringsAsFactors=FALSE)
}

run_pathway_enrichment <- function(input, results, cfg, out_dir, label) {
  lists <- enrichment_gene_lists(input,results,cfg)
  paths <- file.path(out_dir,paste0("Fig4d_",c("background",cfg$groups),"_genes.txt"))
  writeLines(lists$background,paths[1])
  for(i in seq_along(cfg$groups)) writeLines(lists$queries[[cfg$groups[i]]],paths[i+1])
  checksums <- unname(tools::md5sum(paths))
  mode <- cfg$enrichment_mode
  assert(mode %in% c("demo_local","enrichr_exports","enrichr_live"),"Unknown enrichment_mode")
  manifest <- data.frame(Group=cfg$groups,Library=cfg$enrichment_library,
    BackgroundMD5=checksums[1],QueryMD5=checksums[-1],Mode=mode,
    BackgroundGenes=length(lists$background),QueryGenes=lengths(lists$queries),stringsAsFactors=FALSE)
  outputs <- list()
  if(mode=="demo_local") {
    assert(label=="SYNTHETIC DEMONSTRATION" && all(startsWith(lists$background,"SYN_")),
      "demo_local is restricted to fictional SYN_ genes. Use Enrichr exports or explicitly enabled live Enrichr for study data.")
    membership <- read_table(file.path(input$data_dir,"pathways.tsv"),sep="\t")
    check_columns(membership,c("Pathway","GeneSymbol"),"pathways.tsv")
    assert(all(startsWith(membership$Pathway,"DEMO_")),"Offline demonstration requires DEMO_ pathways")
    for(group in cfg$groups) outputs[[group]] <- local_overrepresentation(lists$queries[[group]],lists$background,membership)
    manifest$Library <- "Fictional pathways; local hypergeometric/BH demonstration"
    method <- "Offline fictional-set enrichment; not an Enrichr analysis"
  } else if(mode=="enrichr_exports") {
    folder <- file.path(input$data_dir,cfg$enrichment_export_dir)
    recorded <- read_table(file.path(folder,"manifest.csv"))
    check_columns(recorded,c("Group","Library","BackgroundMD5","QueryMD5"),"Enrichr manifest")
    check_ids(recorded$Group,"Enrichr manifest groups")
    assert(setequal(recorded$Group,cfg$groups),"Enrichr manifest must cover all three groups")
    recorded <- recorded[match(cfg$groups,recorded$Group),]
    assert(identical(recorded$Library,manifest$Library) && identical(recorded$BackgroundMD5,manifest$BackgroundMD5) &&
      identical(recorded$QueryMD5,manifest$QueryMD5),"Enrichr export provenance does not match current library/background/query lists")
    for(group in cfg$groups) outputs[[group]] <- normalise_enrichr(read_table(file.path(folder,paste0(group,".tsv")),sep="\t"))
    method <- "Enrichr exports with declared measured-gene background"
  } else {
    assert(isTRUE(cfg$allow_enrichr_network),"Live Enrichr sends gene lists to its server; add --allow-enrichr-network only when intended")
    assert(!any(startsWith(lists$background,"SYN_")),"Fictional gene names must not be sent to Enrichr")
    assert(requireNamespace("enrichR",quietly=TRUE),"Install enrichR >=3.4 for live custom-background analysis")
    assert("background" %in% names(formals(enrichR::enrichr)),"Installed enrichR lacks custom-background support")
    enrichR::setEnrichrSite("Enrichr")
    libraries <- enrichR::listEnrichrDbs()
    assert(cfg$enrichment_library %in% libraries$libraryName,"Requested Enrichr library is not available; do not substitute another release silently")
    cache <- file.path(out_dir,"enrichr_exports"); dir.create(cache,showWarnings=FALSE)
    for(group in cfg$groups) {
      query <- lists$queries[[group]]
      if(!length(query)) {
        raw <- data.frame(Term=character(),P.value=numeric(),Adjusted.P.value=numeric(),Overlap=character(),Genes=character())
      } else {
        response <- enrichR::enrichr(query,databases=cfg$enrichment_library,
          background=lists$background,include_overlap=TRUE)
        raw <- response[[cfg$enrichment_library]]
        assert(is.data.frame(raw),paste("Enrichr did not return a result table for",group))
      }
      write.table(raw,file.path(cache,paste0(group,".tsv")),sep="\t",row.names=FALSE,quote=TRUE)
      outputs[[group]] <- normalise_enrichr(raw)
    }
    manifest$RetrievedUTC <- format(Sys.time(),tz="UTC",usetz=TRUE)
    manifest$enrichRVersion <- as.character(utils::packageVersion("enrichR"))
    write_csv(manifest,file.path(cache,"manifest.csv"))
    method <- "Enrichr with explicit measured-gene background"
  }
  all <- do.call(rbind,lapply(cfg$groups,function(group) {
    d <- outputs[[group]]
    if(!"QuerySize" %in% names(d)) d$QuerySize <- rep(length(lists$queries[[group]]),nrow(d))
    if(!"BackgroundSize" %in% names(d)) d$BackgroundSize <- rep(length(lists$background),nrow(d))
    d$Group <- rep(group,nrow(d)); d$Method <- rep(method,nrow(d)); d
  })); rownames(all) <- NULL
  write_csv(all,file.path(out_dir,"Fig4d_enrichment_all_results.csv"))
  write_csv(manifest,file.path(out_dir,"Fig4d_enrichment_manifest.csv"))
  plot_enrichment(all,cfg,out_dir,label,method)
  invisible(all)
}

plot_enrichment <- function(all,cfg,out_dir,label,method) {
  # These broad-parent exclusions match the available grouped-dotplot source.
  excluded <- c("Immune System","Adaptive Immune System","Cytokine Signaling in Immune System",
    "Metabolism of Proteins","Metabolism of Lipids","Metabolism of Carbohydrates","Metabolism of RNA",
    "Disease","Signal Transduction","Hemostasis","Cell Cycle","Gene Expression (Transcription)",
    "Developmental Biology","Programmed Cell Death","Transport of Small Molecules","Vesicle-mediated Transport")
  d <- all[!all$Term %in% excluded,,drop=FALSE]
  terms <- unique(d$Term[!is.na(d$FDR) & d$FDR<cfg$fdr])
  if(!length(terms)) {
    p <- ggplot2::ggplot()+ggplot2::annotate("text",x=0,y=0,label="No pathways meet FDR <0.05")+
      ggplot2::theme_void()+ggplot2::labs(title="Pathway enrichment",subtitle=method,caption=plot_caption(label))
    save_plot(p,file.path(out_dir,"Fig4d_pathway_enrichment"),9,6)
    write_csv(d[FALSE,,drop=FALSE],file.path(out_dir,"Fig4d_displayed_pathways.csv"))
    return(invisible(NULL))
  }
  classify <- function(term) {
    z <- d[d$Term==term,,drop=FALSE]
    sig <- setNames(cfg$groups %in% z$Group[z$FDR<cfg$fdr],cfg$groups)
    if(all(sig)) "1. Common to all three"
    else if(sig["higherEF"] && sig["lowerEF"]) "2. HFpEF groups only"
    else if(sig["Control"]) "3. Control-containing (other)"
    else if(sig["higherEF"]) "4. Higher EF only"
    else "5. Lower EF only"
  }
  categories <- vapply(terms,classify,character(1))
  best <- vapply(terms,function(term) min(d$PValue[d$Term==term]),numeric(1))
  selection <- data.frame(Term=terms,Category=categories,BestP=best)
  selection <- do.call(rbind,lapply(split(selection,selection$Category),function(z)
    head(z[order(z$BestP,z$Term),],cfg$enrichment_top_per_category)))
  plotdata <- merge(expand.grid(Term=selection$Term,Group=cfg$groups,stringsAsFactors=FALSE),
    d,by=c("Term","Group"),all.x=TRUE,sort=FALSE)
  plotdata$Category <- selection$Category[match(plotdata$Term,selection$Term)]
  plotdata$Term <- factor(plotdata$Term,levels=rev(selection$Term))
  plotdata$Group <- factor(plotdata$Group,levels=cfg$groups)
  plotdata$Significant <- !is.na(plotdata$FDR) & plotdata$FDR<cfg$fdr
  p <- ggplot2::ggplot(plotdata,ggplot2::aes(Group,Term))+
    ggplot2::geom_point(data=plotdata[!plotdata$Significant,,drop=FALSE],color="grey80",size=1.4)+
    ggplot2::geom_point(data=plotdata[plotdata$Significant,,drop=FALSE],
      ggplot2::aes(size=Count,color=-log10(pmax(FDR,.Machine$double.xmin))))+
    ggplot2::scale_color_gradient(low="#fdcc8a",high="#b30000",name="-log10(FDR)")+
    ggplot2::scale_size_continuous(range=c(2,6),name="Overlapping genes")+
    ggplot2::scale_x_discrete(labels=group_labels)+
    ggplot2::facet_grid(Category~.,scales="free_y",space="free_y")+
    ggplot2::theme_bw(base_size=10)+ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
      strip.text.y=ggplot2::element_text(angle=0,size=8),axis.text.x=ggplot2::element_text(angle=25,hjust=1))+
    ggplot2::labs(title="Pathway enrichment across groups",subtitle=method,x=NULL,y=NULL,caption=plot_caption(label))
  save_plot(p,file.path(out_dir,"Fig4d_pathway_enrichment"),11,max(6,nrow(selection)*.35+2))
  write_csv(plotdata,file.path(out_dir,"Fig4d_displayed_pathways.csv"))
}
