require_packages <- function() {
  needed <- c("limma", "statmod", "ggplot2", "ggrepel", "mixOmics", "VennDiagram")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages using install_dependencies.R: ",
                            paste(missing, collapse = ", "), call. = FALSE)
}

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

check_columns <- function(x, columns, name) {
  absent <- setdiff(columns, names(x))
  assert(!length(absent), paste(name, "is missing:", paste(absent, collapse = ", ")))
  assert(!anyDuplicated(names(x)), paste(name, "has duplicate column names"))
}

check_ids <- function(ids, name) {
  assert(!anyNA(ids) && all(nzchar(trimws(ids))) && !anyDuplicated(ids),
         paste(name, "must be nonempty and unique"))
}

read_table <- function(path, sep = ",") {
  assert(file.exists(path), paste("Missing input:", basename(path)))
  read.table(path, header = TRUE, sep = sep, quote = '"', comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

read_expression <- function(path) {
  d <- read_table(path)
  assert(names(d)[1] == "SeqId", paste(basename(path), "must start with SeqId"))
  check_ids(d$SeqId, "Expression SeqId")
  check_ids(names(d)[-1], "Expression sample IDs")
  assert(all(vapply(d[-1], is.numeric, logical(1))), "All expression values must be numeric")
  x <- as.matrix(d[-1]); rownames(x) <- d$SeqId
  assert(all(is.finite(x)) && all(x > 0), "Expression must contain finite positive RFU without missing values")
  x
}

load_inputs <- function(data_dir, cfg) {
  x <- read_expression(file.path(data_dir, "expressiondata_restandstress.csv"))
  samples <- read_table(file.path(data_dir, "sampledata.csv"))
  pheno <- read_table(file.path(data_dir, "phenotypedata.csv"))
  features <- read_table(file.path(data_dir, "featuredata.csv"))
  check_columns(samples, c("SampleID", "PatientID", "Condition"), "sampledata.csv")
  check_columns(pheno, c("PatientID", "Phenotype", "Age", "Sex", "BMI", cfg$camera_variables, cfg$medication_columns), "phenotypedata.csv")
  check_columns(features, c("SeqId", "GeneSymbol", "TargetName", "UniProtID", "Role"), "featuredata.csv")
  check_ids(samples$SampleID, "SampleID"); check_ids(pheno$PatientID, "PatientID")
  check_ids(features$SeqId, "Feature SeqId")
  assert(setequal(samples$SampleID, colnames(x)), "Sample metadata and expression columns must match exactly")
  assert(setequal(features$SeqId, rownames(x)), "Feature metadata and expression rows must match exactly")
  assert(setequal(samples$PatientID, pheno$PatientID), "Participant IDs must agree between sample and phenotype tables")
  assert(!anyNA(samples$Condition) && all(samples$Condition %in% c("Rest", "Stress")), "Condition must be Rest or Stress")
  assert(all(table(factor(samples$PatientID, levels = pheno$PatientID),
                   factor(samples$Condition, levels = c("Rest", "Stress"))) == 1),
         "Each participant must have exactly one Rest and one Stress sample")
  assert(!anyNA(pheno$Phenotype) && setequal(pheno$Phenotype, cfg$groups), "Phenotype must contain all three configured groups")
  assert(!anyNA(pheno$Sex) && all(pheno$Sex %in% c("Female", "Male")), "Sex must be Female or Male")
  for (v in c("Age", "BMI")) assert(is.numeric(pheno[[v]]) && all(is.finite(pheno[[v]])), paste(v, "must be finite numeric"))
  for (v in cfg$camera_variables) assert(is.numeric(pheno[[v]]) && all(is.na(pheno[[v]]) | is.finite(pheno[[v]])), paste(v, "must be numeric, allowing NA"))
  for (v in cfg$medication_columns) assert(is.numeric(pheno[[v]]) && all(is.na(pheno[[v]]) | pheno[[v]] %in% c(0,1)), paste(v, "must be 0/1, allowing NA"))
  samples <- samples[match(colnames(x), samples$SampleID), , drop = FALSE]
  features <- features[match(rownames(x), features$SeqId), , drop = FALSE]
  rownames(features) <- features$SeqId
  rownames(pheno) <- pheno$PatientID
  pheno$Sex <- factor(pheno$Sex, levels = c("Female", "Male"))
  paired_columns <- function(condition) {
    sub <- samples[samples$Condition == condition, , drop = FALSE]
    sub$SampleID[match(pheno$PatientID, sub$PatientID)]
  }
  rest <- x[, paired_columns("Rest"), drop = FALSE]
  stress <- x[, paired_columns("Stress"), drop = FALSE]
  colnames(rest) <- colnames(stress) <- pheno$PatientID
  cystatin <- features$SeqId[!is.na(features$Role) & features$Role == "cystatin_covariate"]
  assert(length(cystatin) == 1, "Exactly one feature must have Role=cystatin_covariate")
  pheno$BaselineCystatin <- as.numeric(log2(rest[cystatin, ]))
  # Optional split matrices are redundant copies; verify them if present.
  for (condition in c("Rest", "Stress")) {
    file <- file.path(data_dir, paste0("expressiondata-", tolower(condition), ".csv"))
    if (file.exists(file)) {
      m <- read_expression(file)
      target <- if (condition == "Rest") rest else stress
      assert(setequal(rownames(m), rownames(target)) && setequal(colnames(m), colnames(target)),
             paste(basename(file), "has inconsistent IDs"))
      assert(isTRUE(all.equal(unname(m[rownames(target), colnames(target), drop = FALSE]), unname(target), tolerance = 1e-12)),
             paste(basename(file), "does not agree with the combined matrix"))
    }
  }
  list(x = x, rest = rest, stress = stress, samples = samples,
       pheno = pheno, features = features, data_dir = data_dir)
}

check_design <- function(design, label) {
  assert(all(is.finite(design)), paste(label, "design contains non-finite values"))
  assert(qr(design)$rank == ncol(design), paste(label, "design is rank deficient; check covariate variation"))
  assert(nrow(design) - ncol(design) >= 4, paste(label, "needs at least four residual degrees of freedom"))
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")

group_colors <- c(Control = "#808182", higherEF = "#9cd0c6", lowerEF = "#8a70b2")
group_labels <- c(Control = "Control", higherEF = "HF-EF >60%", lowerEF = "HF-EF 50-60%")
category_colors <- c("Not significant" = "#d6d6d6", "Control only" = "#808182",
  "higherEF only" = "#9cd0c6", "lowerEF only" = "#8a70b2", "Control & higherEF" = "#cc84b4",
  "Control & lowerEF" = "#ad669f", "HFpEF only" = "#c7ad32", "All three" = "#2375c8")

category_from_flags <- function(control, higher, lower) {
  code <- as.integer(control) + 2L * as.integer(higher) + 4L * as.integer(lower)
  c("Not significant", "Control only", "higherEF only", "Control & higherEF",
    "lowerEF only", "Control & lowerEF", "HFpEF only", "All three")[code + 1L]
}

save_plot <- function(plot, stem, width = 8, height = 5.5) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                  device = grDevices::pdf, useDingbats = FALSE, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 160, bg = "white")
}

plot_caption <- function(label) {
  if (label == "SYNTHETIC DEMONSTRATION") "SYNTHETIC DEMONSTRATION - fictional data; no biological inference"
  else "User-supplied data - release provenance must be verified separately"
}

versions_table <- function() {
  pkgs <- c("limma", "statmod", "ggplot2", "ggrepel", "mixOmics", "VennDiagram")
  data.frame(Software = c("R", pkgs), Version = c(as.character(getRversion()),
    vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))))
}
