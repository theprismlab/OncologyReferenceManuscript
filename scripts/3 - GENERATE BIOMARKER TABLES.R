library(renv)
renv::restore()
library(magrittr)
library(tidyverse)
library(parallel)
source("scripts/UTILITIES.R")
file <- "data/external inputs/depmap_oncref_manuscript.h5"

set.seed(23)

# -----
# Load the PRISM data ----
# -----

CompoundList <- data.table::fread("data/input data/PRISMOncologyReferenceCompoundList.csv") %>% 
  dplyr::mutate(cn = paste0(SampleID, "::", CompoundPlate))

selected_compounds <- CompoundList %>% 
  dplyr::filter(Prioritized) %>% 
  .$cn %>% unique
  

LAUC <- data.table::fread("data/processed data/PRISMOncologyReferenceLog2AUCMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LAUC <- LAUC[, selected_compounds]

LFC <- data.table::fread("data/processed data/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LFC <- LFC[, word(colnames(LFC), 1,2, sep = fixed("::")) %in% selected_compounds]

CompoundList <- CompoundList %>% 
  dplyr::filter(Prioritized)



# -----
# Target-recovery functions -----
# -----

# This computes univariate biomarkers without any filtering
bm.lauc.complete <- univariate_biomarker_table(Y = LAUC, file = file, q_val_max = 1, rank.max = 1e6)

bm.lauc.complete %>% 
  write_csv("data/results/biomarker results/lauc_univariate_biomarkers_complete.csv")


# Compute the top 100 correlates with q < 0.1 and n > 250 along with annotations
bm.lauc <- target_recovery(Y = LAUC, file = file, compound_annotations = CompoundList) 

bm.lauc %>% 
  write_csv("data/results/biomarker results/lauc_univariate_biomarkers.csv")

# # Compoute per lfc univariate biomarkers with the same constraints
# bm.lfc <- target_recovery(Y = LFC, file = file, compound_annotations = CompoundList)
# 
# bm.lfc %>% 
#   write_csv("results/biomarker results/lfc_univariate_biomarkers.csv")



# ----
# Load DepMap data into a single matrix ----
# ----

X.XPR <- read_dataset(file = file, dataset = "CRISPR")
cl = intersect(rownames(X.XPR), rownames(LAUC));  X.XPR <- X.XPR[cl, ]
X.XPR <- X.XPR[ , colMeans(is.finite(X.XPR)) > 0.9]
colnames(X.XPR) %<>% paste0("CRISPR_", .)

X.RNAi <- read_dataset(file = file, dataset = "RNAi")
cl = intersect(rownames(X.RNAi), rownames(LAUC));  X.RNAi <- X.RNAi[cl, ]
X.RNAi <- X.RNAi[ , colMeans(is.finite(X.RNAi)) > 0.9]
colnames(X.RNAi) %<>% paste0("RNAi_", .)


X.EXP <- read_dataset(file = file, dataset = "Expression")
cl = intersect(rownames(X.EXP), rownames(LAUC));  X.EXP <- X.EXP[cl, ]
X.EXP <- X.EXP[ , colMeans(is.finite(X.EXP)) > 0.9]
colnames(X.EXP) %<>% paste0("EXP_", .)

X.CN <- read_dataset(file = file, dataset = "CopyNumber")
cl = intersect(rownames(X.CN), rownames(LAUC));  X.CN <- X.CN[cl, ]
X.CN <- X.CN[ , colMeans(is.finite(X.CN)) > 0.9]
colnames(X.CN) %<>% paste0("CN_", .)


X.MUT <-  read_dataset(file = file, dataset = "Mutation")
cl = intersect(rownames(X.MUT), rownames(LAUC));  X.MUT <- X.MUT[cl, ]
X.MUT <- X.MUT[ , colMeans(is.finite(X.MUT)) > 0.9]
colnames(X.MUT) %<>% paste0("MUT_", .)


X.LIN <-  read_dataset(file = file, dataset = "Lineage")
cl = intersect(rownames(X.LIN), rownames(LAUC));  X.LIN <- X.LIN[cl, ]
X.LIN <- X.LIN[ , colMeans(is.finite(X.LIN)) > 0.9]
colnames(X.LIN) %<>% paste0("LIN_", .)

X.FUS <-  read_dataset(file = file, dataset = "Fusion")
cl = intersect(rownames(X.FUS), rownames(LAUC));  X.FUS <- X.FUS[cl, ]
X.FUS <- X.FUS[ , colMeans(is.finite(X.FUS)) > 0.9]
colnames(X.FUS) %<>% paste0("FUS_", .)


X <- X.XPR %>%
  reshape2::melt() %>% 
  dplyr::bind_rows(reshape2::melt(X.RNAi)) %>%
  dplyr::bind_rows(reshape2::melt(X.EXP)) %>%
  dplyr::bind_rows(reshape2::melt(X.CN)) %>%
  dplyr::bind_rows(reshape2::melt(X.MUT)) %>%
  dplyr::bind_rows(reshape2::melt(X.FUS)) %>%
  dplyr::bind_rows(reshape2::melt(X.LIN)) %>%
  dplyr::filter(is.finite(value))

X <- tibble(Var1 = rownames(LAUC), dummy = 1) %>% 
  dplyr::left_join(tibble(Var2 = unique(X$Var2), dummy = 1)) %>% 
  dplyr::select(-dummy) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(X) %>% 
  dplyr::group_by(Var2) %>% 
  dplyr::mutate(value = ifelse(is.na(value), median(value, na.rm = T), value)) %>% 
  dplyr::ungroup() %>% 
  reshape2::acast(Var1 ~ Var2)


rm(X.CN, X.MUT, X.EXP, X.FUS, X.RNAi, X.XPR, X.LIN)

cl <- intersect(rownames(LAUC), rownames(X))
X <- X[cl, ]; X <- X[, apply(X, 2, var) > 0.005] 


# -----
# RANDOM FOREST MODELS -----
# -----

RF.lauc <- biomarker_suite_rf_cv(X, LAUC, biomarker_file = file, CompoundList = CompoundList,                                       
                                    bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 20, seed = 23)

RF.lauc %>%  saveRDS("data/results/biomarker results/biomarkers.RDS") 

# ----
# BIOMARKER SUMMARY TABLES ----
# ----


# auxiliary tables
BM.Summary.Table <- RF.lauc$model_performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(model, cn, CompoundName) %>%  
  dplyr::summarise(mse = mean(mse), r.sd = sd(r),  r = mean(r), var.y = mean(var.y.test)) %>% 
  dplyr::mutate(r2 = 1 - mse / var.y) %>% 
  dplyr::group_by(cn) %>%
  dplyr::mutate(r.m = max(r[!model %in% c("targets", "extended")]),
                n.t = length(setdiff(model, c("targets", "extended"))),
                model.class = ifelse(model == "extended", "Extended", 
                                     ifelse(model == "targets", "Targets", 
                                            ifelse(r == r.m, "Best Single Target", "Other Targets")))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-r.m)

DF <- RF.lauc$predictions %>% 
  dplyr::filter(type == "test") %>%
  dplyr::select(-y.hat.oob) %>% 
  dplyr::distinct() %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(CompoundName, model) %>% 
  dplyr::arrange(y.hat) %>% 
  dplyr::mutate(n = 1:n(), N = n(),
                p = var(y) * (1 / n + 1 / (N - n)), 
                cs = cumsum(y), 
                s = sum(y)) %>% 
  dplyr::mutate(m1 = cs / n, m2 = (s - cs) / (N - n),
                t = -(m1 - m2) / sqrt(p), 
                df = N -2 ) %>%
  dplyr::select(-p, -s, -cs, -m1, -m2, -df, -K) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct()


DF <- BM.Summary.Table %>% 
  dplyr::left_join(DF %>% 
                     dplyr::filter(is.finite(t)) %>% 
                     dplyr::group_by(CompoundName, model, N) %>% 
                     dplyr::summarize(t.mean = mean(t),
                                      t.peak = max(t),
                                      n.peak = min(n[t == t.peak])) %>% 
                     dplyr::ungroup()) %>% 
  dplyr::left_join(CompoundList %>% 
                     dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism))




# Scores table
Scores.Table <- BM.Summary.Table %>% 
  dplyr::filter(model.class != "Other Targets") %>% 
  dplyr::distinct(CompoundName, cn, n.t, r, r.sd, model.class) 



Scores.Table <- Scores.Table %>% 
  dplyr::filter(model.class == "Best Single Target") %>% 
  tidyr::pivot_wider(names_from = "model.class", values_from = c("r", "r.sd")) %>% 
  dplyr::full_join(Scores.Table %>% 
                     dplyr::filter(model.class != "Best Single Target") %>% 
                     tidyr::pivot_wider(names_from = "model.class", values_from = c("r", "r.sd"))) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(PolypharmacologyScore = ifelse(n.t > 1, (r_Targets - `r_Best Single Target`) / sqrt((r.sd_Targets^2 + `r.sd_Best Single Target`)/20)  , 0),
                ExcessPredictabilityScore = ifelse(PolypharmacologyScore > 0, 
                                                   (r_Extended - r_Targets) / sqrt((r.sd_Targets^2 + r.sd_Extended)/20) ,
                                                   (r_Extended - `r_Best Single Target`) / sqrt((r.sd_Extended^2 + `r.sd_Best Single Target`)/20))) %>% 
  dplyr::mutate(PolypharmacologyScore = pmax(PolypharmacologyScore, 0),
                ExcessPredictabilityScore = pmax(ExcessPredictabilityScore, 0),
                Best.r = pmax(r_Extended, pmax(r_Targets, `r_Best Single Target`))) %>%
  dplyr::distinct(CompoundName, cn, n.t, PolypharmacologyScore, ExcessPredictabilityScore, Best.r, r_Extended, r_Targets, `r_Best Single Target`) %>% 
  dplyr::ungroup()


Scores.Table <- DF %>% 
  dplyr::filter(model %in% c("targets", "extended")) %>% 
  dplyr::distinct(cn, CompoundName, model, t.mean, t.peak, N, n.peak) %>% 
  dplyr::mutate(SelectivityScore = t.peak - pmax(t.mean, 0)) %>% 
  tidyr::pivot_wider(names_from = model, values_from = c("SelectivityScore", "t.mean", "t.peak", "N", "n.peak")) %>% 
  dplyr::left_join(Scores.Table) %>% 
  dplyr::rename(OnTargetPolypharmacologyScore = PolypharmacologyScore,
                OffTargetPolypharmacologyScore = ExcessPredictabilityScore,
                n.targets = n.t) %>% 
  dplyr::select(cn, CompoundName, Best.r,
                OnTargetPolypharmacologyScore, OffTargetPolypharmacologyScore,
                SelectivityScore_extended, SelectivityScore_targets,
                r_Extended, r_Targets, `r_Best Single Target`,
                n.peak_extended, t.peak_extended, t.mean_extended, N_extended, 
                n.peak_targets, t.peak_targets, t.mean_targets, N_targets,
                n.targets)


# Variable importances 
Importance.Table <- RF.lauc$variable_importances %>% 
  dplyr::filter(K> 0) %>% 
  dplyr::group_by(cn, model, K) %>% 
  dplyr::mutate(imp = imp / sum(imp)) %>%  
  dplyr::group_by(cn, CompoundName, model, var) %>% 
  dplyr::summarise(imp = sum(imp)/20) %>% 
  dplyr::group_by(cn, CompoundName, model) %>%
  dplyr::arrange(desc(imp)) %>% 
  dplyr::mutate(rank = 1:n()) %>% 
  dplyr::ungroup() 

# Predictability results
Predictability.Table <- RF.lauc$model_performances %>% 
  dplyr::filter(K > 0) %>% 
  dplyr::group_by(cn, CompoundName, model) %>% 
  dplyr::summarise_all(function(x) mean(x, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cn, CompoundName, model, mse, r2, r, var.y.test)



Predictability.Table %>%
  write_csv("data/results/biomarker results/model_performances.csv")

Importance.Table %>% 
  write_csv("data/results/biomarker results/variable_importances.csv")

Scores.Table %>%
  write_csv("data/results/biomarker results/model_scores.csv")


# ----
# Biomarkers for TK/RTK Vignette ----
# ----

TKRTK_gene_symbols <- data.table::fread("data/external inputs/TKRTK_gene_symbols.csv")
TKRTK_gene_symbols$GeneSymbol
 
 
TK.RTK.CL <- CompoundList %>% 
   dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism) %>% 
   tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>% 
   dplyr::filter(GeneSymbolOfTargets %in% TKRTK_gene_symbols$GeneSymbol) %>% 
   dplyr::select(-GeneSymbolOfTargets) %>% 
   dplyr::distinct() %>% 
   dplyr::left_join(CompoundList) %>% 
   dplyr::mutate(GeneSymbolOfTargets = paste0(sort(unique(TKRTK_gene_symbols$GeneSymbol)), collapse = ";"))


TK.RTK.LAUC <- data.table::fread("data/processed data/PRISMOncologyReferenceLog2AUCMatrix.csv") %>% 
   column_to_rownames("V1") %>% 
   as.matrix() 

TK.RTK.LAUC <- TK.RTK.LAUC[, TK.RTK.CL$cn]

TK.RTK.BM <- biomarker_suite_rf_cv(X = X, Y = TK.RTK.LAUC,  biomarker_file = file, 
                       CompoundList = TK.RTK.CL,                                       
                       bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 20, seed = 23)


TK.RTK.BM %>% saveRDS("data/results/biomarker results for TK:RTK vignette/tk_rtk_biomarkers.RDS") 


# Variable importances 
TK.RTK.Importance.Table <- TK.RTK.BM$variable_importances %>% 
  dplyr::filter(K> 0) %>% 
  dplyr::group_by(cn, model, K) %>% 
  dplyr::mutate(imp = imp / sum(imp)) %>%  
  dplyr::group_by(cn, CompoundName, model, var) %>% 
  dplyr::summarise(imp = sum(imp)/20) %>% 
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
  dplyr::mutate(OnTargetPolypharmacologyScore = ifelse(n.t > 1, (r_Targets - `r_Best Single Target`) / sqrt((r.sd_Targets^2 + `r.sd_Best Single Target`)/20)  , 0),
                RTKPolypharmacologyScore = (r_TK.RTK - `r_Targets`) / sqrt((r.sd_Targets^2 + `r.sd_TK.RTK`)/20)) %>% 
  dplyr::mutate(OnTargetPolypharmacologyScore = pmax(OnTargetPolypharmacologyScore, 0),
                RTKPolypharmacologyScore = pmax(RTKPolypharmacologyScore, 0), 
                Best.r = pmax(r_Extended, pmax(r_Targets, pmax(r_TK.RTK, `r_Best Single Target`)))) %>%
  dplyr::distinct(CompoundName, cn, Best.r, OnTargetPolypharmacologyScore, RTKPolypharmacologyScore) %>% 
  dplyr::ungroup()



TK.RTK.Predictability.Table %>%
  write_csv("data/results/biomarker results for TK:RTK vignette/tk_rtk_model_performances.csv")

TK.RTK.Importance.Table %>% 
  write_csv("data/results/biomarker results for TK:RTK vignette/tk_rtk_variable_importances.csv")

TK.RTK.Scores.Table %>% 
  write_csv("data/results/biomarker results for TK:RTK vignette/tk_rtk_model_scores.csv")


