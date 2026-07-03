library(tidyverse)

# ----
# 1. READ INPUT FILES 
# ----

compound_annotations <- data.table::fread("data/input_data/PRISMOncologyReferenceCompoundList.csv")
inst_meta <- data.table::fread("data/input_data/PRISMOncologyReferenceLumInstMeta.csv")
analyte_meta <- data.table::fread("data/input_data/PRISMOncologyReferenceLumAnalyteMeta.csv")

# Combine Corrected and Uncorrected DRC
DRC <- dplyr::bind_rows(
  data.table::fread("data/processed_data/PRISMOncologyReferenceLumDoseResponseParameters.csv") %>% 
    dplyr::mutate(response = "corrected"),
  data.table::fread("data/processed_data_without_artifact_correction/PRISMOncologyReferenceLumDoseResponseParameters_no_reg.csv") %>% 
    dplyr::mutate(response = "uncorrected")
)


# ----
# 2. DRC RELEASE FILES
# ----

# Cache the filtered DRC dataset to avoid redundant operations
DRC_filtered <- DRC %>% 
  dplyr::filter(response == "corrected", priority == 1, successful_fit) %>% 
  dplyr::inner_join(
    compound_annotations %>% dplyr::filter(Prioritized) %>% dplyr::distinct(SampleID, CompoundPlate),
    by = c("SampleID", "CompoundPlate")
  )

# Write Response Curves
DRC_filtered %>% 
  dplyr::rename(
    ModelID = depmap_id, 
    EC50 = inflection, 
    LowerAsymptote = lower_limit, 
    UpperAsymptote = upper_limit, 
    Slope = slope
  ) %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  readr::write_csv("data/release_files/PRISMOncologyReferenceLumResponseCurves.csv")

# Write Log2AUC Matrix
DRC_filtered %>% 
  dplyr::distinct(depmap_id, SampleID, log2_auc) %>% 
  tidyr::pivot_wider(names_from = SampleID, values_from = log2_auc) %>% 
  tibble::column_to_rownames("depmap_id") %>% 
  write.csv("data/release_files/PRISMOncologyReferenceLumLog2AUCMatrix.csv")


# ----
# 3. COLLAPSED LFC RELEASE FILES
# ----

# Combine Corrected and Uncorrected Collapsed LFC
LFC.collapsed <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumLFCCollapsed.csv") %>% 
  dplyr::filter(is.finite(LFC)) %>% 
  dplyr::mutate(dose_ix = round(log2(pert_dose) / log2(3))) %>% 
  dplyr::select(-screen) %>% 
  dplyr::distinct() %>% 
  dplyr::full_join(
    data.table::fread("data/processed_data_without_artifact_correction/PRISMOncologyReferenceLumLFCCollapsed_no_reg.csv") %>% 
      dplyr::filter(is.finite(LFC)) %>% 
      dplyr::rename(
        LFC_unccorected = LFC, 
        LFC_uncorrected_fitted = LFC_fitted, 
        outlier_uncorrected = outlier
      ) %>% 
      dplyr::mutate(dose_ix = round(log2(pert_dose) / log2(3))) %>% 
      dplyr::select(-priority, -pert_dose) %>% 
      dplyr::distinct(),
    by = c("CompoundPlate", "SampleID", "pert_dose_unit", "cellset", "pool_id", "depmap_id", "dose_ix") # Specifying join keys safely
  )



# Cache the filtered LFC base dataset
LFC_filtered <- LFC.collapsed %>%
  dplyr::inner_join(
    compound_annotations %>% dplyr::filter(Prioritized) %>% dplyr::distinct(SampleID, CompoundName, CompoundPlate),
    by = c("SampleID", "CompoundPlate")
  ) %>% 
  dplyr::filter(priority == 1, !outlier) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, LFC_fitted, CompoundName, screen, CompoundPlate, pool_id, cellset) %>% 
  dplyr::mutate(Label = paste0(CompoundName, " (", SampleID, ") @", pert_dose, " ", pert_dose_unit)) %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit)

# Write Collapsed Matrix
LFC_filtered %>% 
  dplyr::distinct(Label, depmap_id, LFC_fitted) %>% 
  tidyr::pivot_wider(names_from = Label, values_from = LFC_fitted) %>% 
  tibble::column_to_rownames("depmap_id") %>% 
  write.csv("data/release_files/PRISMOncologyReferenceLumLog2ViabilityCollapsedMatrix.csv")

# Write Collapsed Conditions
LFC_filtered %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit) %>% 
  readr::write_csv("data/release_files/PRISMOncologyReferenceLumLog2ViabilityCollapsedConditions.csv")


# ----
# 4. RAW LFC & VIABILITY PROCESSING
# ----

LFC <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumLFC.csv")  %>% 
  dplyr::select(-LFC) %>% 
  dplyr::rename(LFC = LFC_regressed, LFC_uncorrected = LFC_corrected) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, pool_id, screen, PASS, LFC, LFC_uncorrected)

LFC_ <- LFC_filtered %>%
  dplyr::distinct(screen, CompoundPlate, SampleID, Dose, DoseUnit, depmap_id, pool_id, cellset) %>%
  dplyr::rename(pert_dose = Dose, pert_dose_unit = DoseUnit) %>% 
  dplyr::left_join(inst_meta) %>%
  dplyr::left_join(analyte_meta) %>% 
  dplyr::distinct(prism_replicate, pert_well, analyte_id, cellset, screen, depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::left_join(LFC) %>% 
  dplyr::filter(is.finite(LFC), PASS) %>% 
  dplyr::distinct(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate, LFC) %>% 
  dplyr::group_by(depmap_id, SampleID, pert_dose, pert_dose_unit, CompoundPlate) %>% 
  dplyr::mutate(Replicate = dplyr::row_number()) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(Dose = pert_dose, DoseUnit = pert_dose_unit) 

LFC.conditions <- LFC_ %>% 
  dplyr::distinct(SampleID, Dose, DoseUnit, CompoundPlate, Replicate) %>% 
  dplyr::arrange(SampleID, CompoundPlate, Dose, DoseUnit, Replicate) %>% 
  dplyr::mutate(Label = dplyr::row_number() - 1) 

# Write Viability Matrix
LFC_ %>% 
  dplyr::left_join(LFC.conditions) %>% 
  dplyr::mutate(viability = pmin(2^LFC, 1.5)) %>% 
  dplyr::distinct(depmap_id, Label, viability) %>% 
  tidyr::pivot_wider(names_from = Label, values_from = viability) %>% 
  tibble::column_to_rownames("depmap_id") %>% 
  write.csv("data/release_files/PRISMOncologyReferenceLumViabilityMatrix.csv")

# Write Viability Conditions
LFC.conditions %>% 
  readr::write_csv("data/release_files/PRISMOncologyReferenceLumViabilityConditions.csv")


# ----
# 5. QC & CONFOUNDER MATRIX
# ----

QC <- data.table::fread("data/processed_data/PRISMOncologyReferenceLumQCTable.csv")

# Optimize Confounder generation (Bypasses reshape2::melt entirely)
analyte_filtered <- analyte_meta %>% 
  dplyr::filter(pool_id != "CTLBC", !is.na(depmap_id), !is.na(pool_id), is.na(note)) %>% 
  dplyr::mutate(ix = 1)

conf_pools_1 <- analyte_filtered %>% 
  dplyr::mutate(col = paste(pool_id, cellset, screen, sep = "_")) %>% 
  dplyr::distinct(depmap_id, col, ix) %>% 
  tidyr::pivot_wider(names_from = col, values_from = ix, values_fill = 0)

conf_pools_2 <- analyte_filtered %>% 
  dplyr::mutate(col = paste(cellset, screen, sep = "_")) %>% 
  dplyr::distinct(depmap_id, col, ix) %>% 
  tidyr::pivot_wider(names_from = col, values_from = ix, values_fill = 0)

conf_pools <- conf_pools_1 %>% 
  dplyr::full_join(conf_pools_2) %>% 
  tibble::column_to_rownames("depmap_id")

# Drop identical columns safely
conf_pools <- conf_pools[, !duplicated(t(conf_pools))]
conf_pools <- conf_pools %>% tibble::rownames_to_column("depmap_id")

# Aggregate QC metrics
conf_QC <- QC %>% 
  dplyr::inner_join(
    LFC.collapsed %>% dplyr::filter(priority == 1) %>% dplyr::distinct(screen, cellset, pool_id, depmap_id)
  ) %>% 
  dplyr::group_by(depmap_id) %>% 
  dplyr::summarise(dplyr::across(c(NC.median, NC.mad, PC.median, PC.mad, DR, SSMD), ~ median(.x[is.finite(.x)], na.rm = TRUE)))

# Merge and Write Confounder Matrix directly
conf_pools %>% 
  dplyr::full_join(conf_QC) %>% 
  tibble::column_to_rownames("depmap_id") %>% 
  write.csv("data/release_files/PRISMOncologyReferenceLumConfounderMatrix.csv")


# ----
# 6. WRITE FINAL BASE TABLES
# ----

readr::write_csv(QC, "data/release_files/PRISMOncologyReferenceLumQCTable.csv")
readr::write_csv(LFC, "data/release_files/PRISMOncologyReferenceLumLFC.csv")
readr::write_csv(LFC.collapsed, "data/release_files/PRISMOncologyReferenceLumLFCCollapsed.csv")
readr::write_csv(DRC, "data/release_files/PRISMOncologyReferenceLumDoseResponseParameters.csv")
