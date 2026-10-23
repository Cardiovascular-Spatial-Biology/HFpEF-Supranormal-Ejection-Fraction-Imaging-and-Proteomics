# Explicit analysis settings. Please see docs/METHODS_AND_PROVENANCE.md before reuse.
config <- list(
  groups = c("Control", "higherEF", "lowerEF"),
  fdr = 0.05,
  log2fc = log2(1.25),
  tissue_specificity_fold = 5,
  tissues = c("skeletal muscle", "heart muscle", "liver", "lung",
              "lymphoid tissue", "bone marrow", "brain"),
  camera_variables = c("DeltaLVSVi", "DeltaLVEF", "DeltaCI", "DeltaLAESVi",
                       "DeltaLungWater", "SixMinuteWalk", "KCCQ"),
  camera_inter_gene_cor = 0.01, # Primary analysis setting; NA estimates correlations.
  camera_levels = c(2L, 3L),
  camera_set_min = 10L,
  camera_set_max = 500L,
  camera_top_nominal = 10L,
  medication_columns = c("ACEi", "ARB", "MRA", "BetaBlocker", "SGLT2i", "LoopDiuretic"),
  medication_labels = c(ACEi = "ACE inhibitor", ARB = "Angiotensin receptor blocker",
    MRA = "Mineralocorticoid receptor antagonist", BetaBlocker = "Beta blocker",
    SGLT2i = "SGLT2 inhibitor", LoopDiuretic = "Loop diuretic"),
  seed = 20261019L
)
