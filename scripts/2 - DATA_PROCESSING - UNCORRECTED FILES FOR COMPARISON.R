library(tidyverse)
library(magrittr)
library(scam)
library(parallel)
library(renv)
source("scripts/UTILITIES.R")
set.seed(23)
renv::restore()

#----
# LOAD THE LFC DATA
#----

inst_meta <- data.table::fread("data/input data/PRISMOncologyReferenceInstMeta.csv")
analyte_meta <- data.table::fread("data/input data/PRISMOncologyReferenceAnalyteMeta.csv")
LFC <- data.table::fread("data/processed data/PRISMOncologyReferenceLFC.csv")

# ----
# COLLAPSE REPLICATES 
# ----

LFC.collapsed <- LFC %>% 
  dplyr::filter(PASS,  is.finite(LFC_corrected)) %>% 
  dplyr::left_join(inst_meta) %>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(pert_type == "trt_cp", !is.na(depmap_id), pool_id != "CTLBC") %>%  
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, pool_id, depmap_id) %>% 
  dplyr::filter(n() > 1) %>% 
  dplyr::summarise(LFC = median(LFC_corrected, na.rm = T)) %>% 
  dplyr::ungroup()


# ----
# FLAG AND FILTER OUTLIER TREATMENT WELLS
# ----

dose_indices <- inst_meta %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose) %>% 
  dplyr::group_by(CompoundPlate, SampleID) %>% 
  dplyr::arrange(pert_dose) %>% 
  dplyr::mutate(dose_ix = 1:n()) %>%
  dplyr::ungroup()


monotonicity_flags <- dose_indices %>% 
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.prev = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix + 1)) %>%   
  dplyr::left_join(dose_indices %>% 
                     dplyr::rename(pert_dose.next = pert_dose) %>% 
                     dplyr::mutate(dose_ix = dose_ix - 1)) %>%   
  dplyr::select(-dose_ix) %>% 
  dplyr::left_join(LFC.collapsed) %>%  
  dplyr::left_join(LFC.collapsed %>% 
                     dplyr::rename(pert_dose.prev = pert_dose,  LFC.prev = LFC)) %>%  
  dplyr::left_join(LFC.collapsed %>% 
                     dplyr::rename(pert_dose.next = pert_dose,  LFC.next = LFC)) %>%  
  dplyr::mutate(flag.down = LFC < -2, flag.up =  LFC > -1,
                flag.down.prev = ifelse(is.na(pert_dose.prev), FALSE, LFC.prev < -2), 
                flag.up.prev =  ifelse(is.na(pert_dose.prev), TRUE, LFC.prev > -1),
                flag.down.next = ifelse(is.na(pert_dose.next), TRUE, LFC.next < -2),
                flag.up.next =  ifelse(is.na(pert_dose.next), FALSE, LFC.next > -1)) %>% 
  dplyr::mutate(flag1 = flag.down & flag.up.next & flag.up.prev,
                flag2 = flag.up & flag.down.next & flag.down.prev)


outlier.trt.pools <- monotonicity_flags %>% 
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, cellset, pool_id) %>% 
  dplyr::summarise(n.f1 = mean(flag1, na.rm = T),
                   n.f2 = mean(flag2, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(pmax(n.f1, n.f2) > 0.25) %>% 
  dplyr::distinct(CompoundPlate, SampleID, pert_dose, cellset, pool_id) %>% 
  dplyr::mutate(outlier = TRUE)



LFC.collapsed %<>% 
  dplyr::left_join(outlier.trt.pools) %>% 
  dplyr::mutate(outlier = !is.na(outlier))


LFC %<>%  
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(outlier.trt.pools) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, pool_id, cellset, screen, LFC, LFC_corrected, LFC_regressed, PASS, outlier) %>% 
  dplyr::mutate(outlier = !is.na(outlier))



# -------
# FITTING DOSE RESPONSE CURVES 
# ------

LFC.FILTERED <- LFC %>%
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::filter(is.finite(LFC_corrected), PASS, !outlier, pool_id != "CTLBC", pert_type == "trt_cp") %>%
  dplyr::group_by(CompoundPlate, SampleID, pert_dose, analyte_id, depmap_id, pool_id, cellset) %>%
  dplyr::filter(abs(pmin(2^LFC_corrected, 1.5) - median(pmin(2^LFC_corrected, 1.5))) < 0.75) %>%
  dplyr::group_by(CompoundPlate, SampleID, analyte_id, cellset) %>%
  dplyr::filter(length(unique(pert_dose)) > 4, n() > 12) %>%
  dplyr::ungroup()

CompoundPlates <- dplyr::distinct(LFC.FILTERED, SampleID, CompoundPlate) %>% 
  dplyr::count(CompoundPlate) %>% 
  dplyr::arrange(n)

f <- function(df){
  df %>% 
    dplyr::group_by(SampleID, CompoundPlate, analyte_id, depmap_id, pool_id, cellset, screen) %>%    
    dplyr::summarise(get_best_fit(pmin(2^LFC_corrected,1.5), pert_dose)) %>%
    dplyr::ungroup()
}



DRC <- list() 
time = Sys.time()
for(ix in 1:nrow(CompoundPlates)){
  
  print(CompoundPlates[ix, ])
  print(Sys.time() - time)
  time = Sys.time()
  
  DRC[[ix]] <- CompoundPlates[ix, ] %>% 
    dplyr::select(-n) %>% 
    dplyr::left_join(LFC.FILTERED) %>% 
    dplyr::group_split(SampleID) %>% 
    parallel::mclapply(f, mc.cores = parallel::detectCores() - 1) %>%
    dplyr::bind_rows()
}

DRC <- bind_rows(DRC)

DRC %<>% 
  dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(note = ifelse(is.na(note), "", note)) %>% 
  dplyr::filter(note != "REMOVE") %>% 
  dplyr::group_by(CompoundPlate, SampleID, depmap_id) %>% 
  dplyr::arrange(!is.finite(auc), note == "DEPRIORITIZE", desc(cellset)) %>% 
  dplyr::mutate(priority = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-note)


DRC %>% 
  write_csv("data/processed data without artifact correction (for benchmarking)/PRISMOncologyReferenceDoseResponseParameters_no_reg.csv")


# ----
# CALCULATE THE FITTED LFC VALUES
# ----

LFC.fitted <- inst_meta %>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit) %>%  
  dplyr::inner_join(DRC) %>% 
  dplyr::filter(successful_fit) %>% 
  dplyr::mutate(inflection = as.numeric(inflection)) %>% 
  # dplyr::left_join(analyte_meta) %>% 
  dplyr::mutate(LFC_fitted = log2(lower_limit + (upper_limit - lower_limit) / (1 + (pert_dose / inflection)^-slope))) %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, depmap_id, pool_id, LFC_fitted, priority)


LFC.collapsed %<>% 
  dplyr::full_join(LFC.fitted) 

LFC.collapsed %>%  
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, pool_id, depmap_id, LFC, LFC_fitted, outlier, priority) %>%
  write_csv("data/processed data without artifact correction (for benchmarking)/PRISMOncologyReferenceLFCCollapsed_no_reg.csv")


# -----
# WRITE DATA MATRICES FILES
# -----

DRC %>% 
  dplyr::filter(priority == 1, successful_fit) %>% 
  dplyr::mutate(cn = paste0(SampleID , "::", CompoundPlate)) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "log2_auc") %>% 
  write.csv("data/processed data without artifact correction (for benchmarking)/PRISMOncologyReferenceLog2AUCMatrix_no_reg.csv")


LFC.collapsed %>% 
  dplyr::distinct(screen, CompoundPlate, SampleID, pert_dose, pert_dose_unit, cellset, pool_id, depmap_id, LFC, LFC_fitted, outlier, priority) %>%
  dplyr::filter(!outlier, priority == 1) %>% 
  dplyr::mutate(cn = paste0(SampleID,"::", CompoundPlate, "::", pert_dose )) %>% 
  reshape2::acast(depmap_id ~ cn, value.var = "LFC_fitted") %>%  
  write.csv("data/processed data without artifact correction (for benchmarking)/PRISMOncologyReferenceLog2ViabilityCollapsedMatrix_no_reg.csv")




