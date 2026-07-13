library(renv)
renv::restore()
library(tidyverse)
library(parallel)
library(data.table)
library(matrixStats) # Added for fast matrix operations

source("scripts/UTILITIES.R")
# file_path <- "data/external_input/depmap_26Q1_internal.h5" # !!!
file_path <- "data/external_input/depmap_oncref_manuscript.h5"

set.seed(23)

# -----
# 1. Load the PRISM data ----
# -----

LAUC <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumLog2AUCMatrix.csv") %>% 
  tibble::column_to_rownames("V1") %>% 
  as.matrix() 

# LFC <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumLog2ViabilityCollapsedMatrix.csv") %>% 
#  tibble::column_to_rownames("V1") %>% 
#  as.matrix() 

CompoundList <- data.table::fread("data/input_data/PRISMOncologyReferenceCompoundList.csv") %>% 
  dplyr::mutate(cn = paste0(SampleID, "::", CompoundPlate)) %>%
  dplyr::filter(Prioritized, cn %in% colnames(LAUC))

selected_compounds <- CompoundList %>% dplyr::pull(cn) %>% unique()


LAUC <- LAUC[, selected_compounds, drop = FALSE]

# LFC <- LFC[, stringr::word(colnames(LFC), 1, 2, sep = stringr::fixed("::")) %in% selected_compounds, drop = FALSE]


# -----
# 2. Target-recovery functions -----
# -----

# Compute the top 100 correlates with q < 0.1 and n > 250 along with annotations
bm.lauc <- target_recovery(Y = LAUC, file = file_path, compound_annotations = CompoundList)
readr::write_csv(bm.lauc, "data/biomarker_results/lauc_univariate_biomarkers.csv")

## Compute per lfc univariate biomarkers with the same constraints
# bm.lfc <- target_recovery(Y = LFC, file = file_path, compound_annotations = CompoundList)
# readr::write_csv(bm.lfc, "data/biomarker_results/lfc_univariate_biomarkers.csv")


# ----
# 3. Load DepMap data into a single matrix (Optimized) ----
# ----

datasets_to_load <- c(
  "CRISPR" = "CRISPR", 
  "RNAi" = "RNAi", 
  "Expression" = "EXP", 
  "CopyNumber" = "CN", 
  "Mutation" = "MUT", 
  "Lineage" = "LIN", 
  "Fusion" = "FUS"
)

target_rows <- rownames(LAUC)

# Loop over datasets, clean, align to LAUC, and impute medians natively as matrices
matrix_list <- lapply(names(datasets_to_load), function(ds_name) {
  prefix <- datasets_to_load[[ds_name]]
  mat <- read_dataset(file = file_path, dataset = ds_name)
  
  # Intersect with LAUC initially
  cl <- intersect(rownames(mat), target_rows)
  mat <- mat[cl, , drop = FALSE]
  
  # Keep columns with > 90% finite values
  mat <- mat[, colMeans(is.finite(mat)) > 0.9, drop = FALSE]
  colnames(mat) <- paste0(prefix, "_", colnames(mat))
  
  # Initialize an empty matrix perfectly aligned to LAUC rows
  aligned_mat <- matrix(NA_real_, nrow = length(target_rows), ncol = ncol(mat), 
                        dimnames = list(target_rows, colnames(mat)))
  
  # Fill existing rows
  common_rows <- intersect(target_rows, rownames(mat))
  aligned_mat[common_rows, ] <- mat[common_rows, ]
  
  # Fast median imputation for missing values
  col_meds <- matrixStats::colMedians(aligned_mat, na.rm = TRUE)
  na_idx <- which(is.na(aligned_mat), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    aligned_mat[na_idx] <- col_meds[na_idx[, 2]]
  }
  
  return(aligned_mat)
})

# Column-bind all matrices instantly (Replaces reshape2::melt -> join -> acast)
X <- do.call(cbind, matrix_list)

# Filter final matrix by variance
X <- X[, matrixStats::colVars(X, na.rm = TRUE) > 0.005, drop = FALSE]


# -----
# 4. RANDOM FOREST MODELS -----
# -----

RF.lauc <- biomarker_suite_rf_cv(
  X, LAUC,
  biomarker_file = file_path,
  CompoundList = CompoundList,
  bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 10, seed = 23
  )

saveRDS(RF.lauc, "data/biomarker_results/biomarkers.RDS")


# ----
# 5. BIOMARKER SUMMARY TABLES ----
# ----

# Auxiliary tables
BM.Summary.Table <- RF.lauc$model_performances %>%
  dplyr::filter(K > 0) %>%
  dplyr::group_by(model, cn, CompoundName) %>%
  dplyr::summarise(mse = mean(mse), r.sd = sd(r), r = mean(r), var.y = mean(var.y.test), .groups = "drop_last") %>%
  dplyr::mutate(r2 = 1 - mse / var.y) %>%
  dplyr::group_by(cn) %>%
  dplyr::mutate(
    r.m = max(r[!model %in% c("targets", "extended")], na.rm = TRUE),
    n.t = length(setdiff(model, c("targets", "extended"))),
    model.class = dplyr::case_when(
      model == "extended" ~ "Extended",
      model == "targets" ~ "Targets",
      r == r.m ~ "Best Single Target",
      TRUE ~ "Other Targets"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-r.m)

DF <- RF.lauc$predictions %>%
  dplyr::filter(type == "test") %>%
  dplyr::select(-y.hat.oob) %>%
  dplyr::distinct() %>%
  tidyr::drop_na() %>%
  dplyr::group_by(CompoundName, model) %>%
  dplyr::arrange(y.hat) %>%
  dplyr::mutate(
    n = dplyr::row_number(),
    N = dplyr::n(),
    p = var(y) * (1 / n + 1 / (N - n)),
    cs = cumsum(y),
    s = sum(y),
    m1 = cs / n,
    m2 = (s - cs) / (N - n),
    t = -(m1 - m2) / sqrt(p)
  ) %>%
  dplyr::select(-p, -s, -cs, -m1, -m2, -K) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()


DF <- BM.Summary.Table %>%
  dplyr::left_join(
    DF %>%
      dplyr::filter(is.finite(t)) %>%
      dplyr::group_by(CompoundName, model, N) %>%
      dplyr::summarize(
        t.mean = mean(t),
        t.peak = max(t),
        n.peak = min(n[t == t.peak]),
        .groups = "drop"
      ),
    by = c("CompoundName", "model")
  ) %>%
  dplyr::left_join(
    CompoundList %>% dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism),
    by = "CompoundName"
  )


# Scores table - Pivot directly without splitting the dataframe
Scores.Table <- BM.Summary.Table %>%
  dplyr::filter(model.class != "Other Targets") %>%
  dplyr::distinct(CompoundName, cn, n.t, r, r.sd, model.class) %>%
  tidyr::pivot_wider(
    names_from = model.class,
    values_from = c(r, r.sd),
    names_glue = "{.value}_{model.class}" # Formats as r_Targets, r.sd_Targets, etc.
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    PolypharmacologyScore = ifelse(n.t > 1, (r_Targets - `r_Best Single Target`) / sqrt((`r.sd_Targets`^2 + `r.sd_Best Single Target`^2) / 10), 0),
    ExcessPredictabilityScore = ifelse(PolypharmacologyScore > 0,
                                       (r_Extended - r_Targets) / sqrt((`r.sd_Targets`^2 + `r.sd_Extended`^2) / 10),
                                       (r_Extended - `r_Best Single Target`) / sqrt((`r.sd_Extended`^2 + `r.sd_Best Single Target`^2) / 10))
  ) %>%
  dplyr::mutate(
    PolypharmacologyScore = pmax(PolypharmacologyScore, 0),
    ExcessPredictabilityScore = pmax(ExcessPredictabilityScore, 0),
    Best.r = pmax(r_Extended, r_Targets, `r_Best Single Target`, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(CompoundName, cn, n.t, PolypharmacologyScore, ExcessPredictabilityScore, Best.r, r_Extended, r_Targets, `r_Best Single Target`)


Scores.Table <- DF %>%
  dplyr::filter(model %in% c("targets", "extended")) %>%
  dplyr::distinct(cn, CompoundName, model, t.mean, t.peak, N, n.peak) %>%
  dplyr::mutate(SelectivityScore = t.peak - pmax(t.mean, 0)) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(SelectivityScore, t.mean, t.peak, N, n.peak)) %>%
  dplyr::left_join(Scores.Table, by = c("cn", "CompoundName")) %>%
  dplyr::rename(
    OnTargetPolypharmacologyScore = PolypharmacologyScore,
    OffTargetPolypharmacologyScore = ExcessPredictabilityScore,
    n.targets = n.t
  ) %>%
  dplyr::select(cn, CompoundName, Best.r,
                OnTargetPolypharmacologyScore, OffTargetPolypharmacologyScore,
                SelectivityScore_extended, SelectivityScore_targets,
                r_Extended, r_Targets, `r_Best Single Target`,
                n.peak_extended, t.peak_extended, t.mean_extended, N_extended,
                n.peak_targets, t.peak_targets, t.mean_targets, N_targets,
                n.targets)


# Variable importances
Importance.Table <- RF.lauc$variable_importances %>%
  dplyr::filter(K > 0) %>%
  dplyr::group_by(cn, model, K) %>%
  dplyr::mutate(imp = imp / sum(imp)) %>%
  dplyr::group_by(cn, CompoundName, model, var) %>%
  dplyr::summarise(imp = sum(imp) / 10, .groups = "drop_last") %>%
  dplyr::arrange(desc(imp)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::ungroup()

# Predictability results
Predictability.Table <- RF.lauc$model_performances %>%
  dplyr::filter(K > 0) %>%
  dplyr::group_by(cn, CompoundName, model) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  dplyr::select(cn, CompoundName, model, mse, r2, r, var.y.test)


# ----
# 6. SAVE RESULTS ----
# ----

readr::write_csv(Predictability.Table, "data/biomarker_results/model_performances.csv")
readr::write_csv(Importance.Table, "data/biomarker_results/variable_importances.csv")
readr::write_csv(Scores.Table, "data/biomarker_results/model_scores.csv")



# ----
# 7. Biomarkers for TK/RTK Vignette ----
# ----

TKRTK_gene_symbols <- data.table::fread("data/external_input/TKRTK_gene_symbols.csv")



TK.RTK.CL <- CompoundList %>% 
  dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism) %>% 
  tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>% 
  dplyr::filter(GeneSymbolOfTargets %in% TKRTK_gene_symbols$GeneSymbol) %>% 
  dplyr::select(-GeneSymbolOfTargets) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(CompoundList) %>% 
  dplyr::mutate(GeneSymbolOfTargets = paste0(sort(unique(TKRTK_gene_symbols$GeneSymbol)), collapse = ";"))


TK.RTK.LAUC <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumLog2AUCMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 

TK.RTK.LAUC <- TK.RTK.LAUC[, TK.RTK.CL$cn]

TK.RTK.BM <- biomarker_suite_rf_cv(X = X, Y = TK.RTK.LAUC,
                                   biomarker_file = file_path,
                                   CompoundList = TK.RTK.CL,
                                   bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 10, seed = 23)


TK.RTK.BM %>% saveRDS("data/biomarker_results_for_TK_RTK_vignette/tk_rtk_biomarkers.RDS")


# TK.RTK.BM <- readRDS("data/biomarker_results_for_TK_RTK_vignette/tk_rtk_biomarkers.RDS")
# RF.lauc <- readRDS("data/biomarker_results/biomarkers.RDS")

# Variable importances 
TK.RTK.Importance.Table <- TK.RTK.BM$variable_importances %>% 
  dplyr::filter(K> 0) %>% 
  dplyr::group_by(cn, model, K) %>% 
  dplyr::mutate(imp = imp / sum(imp)) %>%  
  dplyr::group_by(cn, CompoundName, model, var) %>% 
  dplyr::summarise(imp = sum(imp)/10) %>% 
  dplyr::group_by(cn, CompoundName, model) %>%
  dplyr::arrange(desc(imp)) %>% 
  dplyr::mutate(rank = 1:n()) %>% 
  dplyr::ungroup() 

# Predictability results
TK.RTK.Predictability.Table <- TK.RTK.BM$model_performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(cn, CompoundName, model) %>% 
  dplyr::summarise_all(function(x) mean(x, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cn, CompoundName, model, mse, r2, r, var.y.test)


# Scores table
TK.RTK.performances <- RF.lauc$model_performances %>%
  dplyr::filter(CompoundName %in% TK.RTK.BM$model_performances$CompoundName,
                model == "targets") %>% 
  dplyr::mutate(model = "real_targets") %>% 
  dplyr::bind_rows(TK.RTK.BM$model_performances)


TK.RTK.BM.Summary.Table <- TK.RTK.performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(model, cn, CompoundName) %>%  
  dplyr::summarise(mse = mean(mse), r.sd = sd(r),  r = mean(r), var.y = mean(var.y.test)) %>% 
  dplyr::mutate(r2 = 1 - mse / var.y) %>% 
  dplyr::group_by(cn) %>%
  dplyr::mutate(r.m = max(r[!model %in% c("targets", "extended", "real_targets")]),
                n.t = length(setdiff(model, c("targets", "extended", "real_targets"))),
                model.class = ifelse(model == "extended", "Extended", 
                                     ifelse(model == "targets", "TK.RTK",
                                            ifelse(model == "real_targets", "Targets",
                                                   ifelse(r == r.m, "Best Single Target", "Other Targets"))))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-r.m)



TK.RTK.Scores.Table <- TK.RTK.BM.Summary.Table %>% 
  dplyr::filter(model.class != "Other Targets") %>% 
  dplyr::distinct(CompoundName, cn, n.t, r, r.sd, model.class) 


TK.RTK.Scores.Table <- TK.RTK.Scores.Table %>% 
  dplyr::filter(model.class == "Best Single Target") %>% 
  tidyr::pivot_wider(names_from = "model.class", values_from = c("r", "r.sd")) %>% 
  dplyr::full_join(TK.RTK.Scores.Table %>% 
                     dplyr::filter(model.class != "Best Single Target") %>% 
                     tidyr::pivot_wider(names_from = "model.class", values_from = c("r", "r.sd"))) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(OnTargetPolypharmacologyScore = ifelse(n.t > 1, (r_Targets - `r_Best Single Target`) / sqrt((r.sd_Targets^2 + `r.sd_Best Single Target`)/10)  , 0),
                RTKPolypharmacologyScore = (r_TK.RTK - `r_Targets`) / sqrt((r.sd_Targets^2 + `r.sd_TK.RTK`)/10)) %>% 
  dplyr::mutate(OnTargetPolypharmacologyScore = pmax(OnTargetPolypharmacologyScore, 0),
                RTKPolypharmacologyScore = pmax(RTKPolypharmacologyScore, 0), 
                Best.r = pmax(r_Extended, pmax(r_Targets, pmax(r_TK.RTK, `r_Best Single Target`)))) %>%
  dplyr::distinct(CompoundName, cn, Best.r, OnTargetPolypharmacologyScore, RTKPolypharmacologyScore) %>% 
  dplyr::ungroup()



TK.RTK.Predictability.Table %>%
  write_csv("data/biomarker_results_for_TK_RTK_vignette/tk_rtk_model_performances.csv")

TK.RTK.Importance.Table %>% 
  write_csv("data/biomarker_results_for_TK_RTK_vignette/tk_rtk_variable_importances.csv")

TK.RTK.Scores.Table %>% 
  write_csv("data/biomarker_results_for_TK_RTK_vignette/tk_rtk_model_scores.csv")