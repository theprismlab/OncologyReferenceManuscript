library(magrittr)
library(tidyverse)
library(parallel)


source("scripts/UTILITIES.R")
file <- "extra/depmap_25Q3_public.h5"


# -----
# Load the PRISM data ----
# -----

CompoundList <- data.table::fread("data/PRISMOncologyReferenceCompoundList.csv") %>% 
  dplyr::mutate(cn = paste0(SampleID, "::", CompoundPlate))

selected_compounds <- CompoundList %>% 
  dplyr::filter(Prioritized) %>% 
  .$cn %>% unique
  

LAUC <- data.table::fread("data/PRISMOncologyReferenceLog2AUCMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LAUC <- LAUC[, selected_compounds]

LFC <- data.table::fread("data/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LFC <- LFC[, word(colnames(LFC), 1,2, sep = fixed("::")) %in% selected_compounds]


LAUC0 <- data.table::fread("data/PRISMOncologyReferenceLog2AUCMatrix_no_reg.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LAUC0 <- LAUC0[, selected_compounds]

LFC0 <- data.table::fread("data/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix_no_reg.csv") %>% 
  column_to_rownames("V1") %>% 
  as.matrix() 
LFC0 <- LFC0[, word(colnames(LFC0), 1,2, sep = fixed("::")) %in% selected_compounds]




CompoundList <- CompoundList %>% 
  dplyr::filter(Prioritized)



# -----
# Target-recovery functions -----
# -----

# This computes univariate biomarkers without any filtering
bm.lauc.complete <- univariate_biomarker_table(Y = LAUC, file = file, q_val_max = 1, rank.max = 1e6)

bm.lauc.complete %>% 
  write_csv("results/lauc_univariate_biomarkers_complete.csv")


# Compute the top 100 correlates with q < 0.1 and n > 250 along with annotations
bm.lauc <- target_recovery(Y = LAUC, file = file, compound_annotations = CompoundList) 

bm.lauc %>% 
  write_csv("results/lauc_univariate_biomarkers.csv")

# Compoute per lfc univariate biomarkers with the same constraints
bm.lfc <- target_recovery(Y = LFC, file = file, compound_annotations = CompoundList)

bm.lfc %>% 
  write_csv("results/lfc_univariate_biomarkers.csv")



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

RF.lauc %>%  saveRDS("results/biomarkers.RDS") 


# ----
# Biomarkers for TK/RTK Vignette ----
# ----

TKRTK_gene_symbols <- data.table::fread("extra/TKRTK_gene_symbols.csv")
TKRTK_gene_symbols$GeneSymbol
 
 
TK.RTK.CL <- CompoundList %>% 
   dplyr::distinct(CompoundName, GeneSymbolOfTargets, TargetOrMechanism) %>% 
   tidyr::separate_rows(GeneSymbolOfTargets, sep = ";") %>% 
   dplyr::filter(GeneSymbolOfTargets %in% TKRTK_gene_symbols$GeneSymbol) %>% 
   dplyr::select(-GeneSymbolOfTargets) %>% 
   dplyr::distinct() %>% 
   dplyr::left_join(CompoundList) %>% 
   dplyr::mutate(GeneSymbolOfTargets = paste0(sort(unique(TKRTK_gene_symbols$GeneSymbol)), collapse = ";"))


TK.RTK.LAUC <- data.table::fread("data/PRISMOncologyReferenceLog2AUCMatrix.csv") %>% 
   column_to_rownames("V1") %>% 
   as.matrix() 

TK.RTK.LAUC <- TK.RTK.LAUC[, TK.RTK.CL$cn]

TK.RTK.BM <- biomarker_suite_rf_cv(X = X, Y = TK.RTK.LAUC,  biomarker_file = file, 
                       CompoundList = TK.RTK.CL,                                       
                       bm_th = 0.05, bm_R = 10, bm_R2 = 50, K = 20, seed = 23)


TK.RTK.BM %>% saveRDS("results/tk_rtk_biomarkers.RDS") 
