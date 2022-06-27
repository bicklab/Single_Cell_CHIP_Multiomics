source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")

tet2_patient_data = readRDS("input_data/CH21014_7754G_C_mutant.rds")
dnmt3a_patient_data = readRDS("input_data/CH21046_mutants.rds")

quality_control = function(data) {
  subset(data, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & nCount_RNA < 10000)
}

index = c("NK", "Treg", "CD14 Mono", "CD8 T cell", "CD4 T cell", "Eryth", "CD16 Mono", "DC",
          "CD8 T cell", "CD4 Naive", "CD8 Naive", "HSPC", "MAIT", "CD8 T cell", "CD4 T cell", "B",
          "B", "Platelet", "NK", "B", "DC", "CD4 CTL", "NK", "CD4 T cell", "Plasmablast",
          "NK", "dnT", "DC", "CD8 T cell")

values = c("NK", "CD4 T cell", "CD14 Mono", "CD8 T cell", "CD4 T cell", "Eryth", "CD16 Mono", "DC",
           "CD8 T cell", "CD4 T cell", "CD8 T cell", "HSPC", "MAIT", "CD8 T cell", "CD4 T cell", "B",
           "B", "Platelet", "NK", "B", "DC", "CD4 T cell", "NK", "CD4 T cell", "Plasmablast",
           "NK", "CD4 T cell", "DC", "CD8 T cell")

index_de = c("CD8 TCM", "CD8 TEM", "CD8 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Naive", "Treg", "CD14 Mono", "CD16 Mono",
             "B", "Eryth", "NK", "DC", "HSPC", "MAIT", "Platelet", "Plasmablast", "dnT", "gdT",
             "B naive", "NK_CD56bright", "B intermediate", "pDC", "B memory", "cDC2", "ILC", "ASDC", "CD4 Proliferating", "CD8 Proliferating",
             "NK Proliferating")
values_de = c("CD8 T cell", "CD8 T cell", "CD8 Naive", "CD4 T cell", "CD4 T cell", "CD4 CTL", "CD4 Naive", "Treg", "CD14 Mono", "CD16 Mono",
              "B", "Eryth", "NK", "DC", "HSPC", "MAIT", "Platelet", "Plasmablast", "dnT", "gdT",
              "B", "NK", "B", "DC", "B", "DC", "ILC", "DC", "CD4 Proliferating", "CD8 Proliferating",
              "NK Proliferating")

combine_cell_types = function(data) {
  data@meta.data$celltype.bh = values[match(data@meta.data$collapsed.celltype, index)]

  # make new column for de, pull from predicted celltype l2
  data@meta.data$celltype.de = values_de[match(data@meta.data$predicted.celltype.l2, index_de)]

  return(data)
}

tet2_patient_data = combine_cell_types(tet2_patient_data) %>% quality_control()
dnmt3a_patient_data = combine_cell_types(dnmt3a_patient_data) %>% quality_control()

tet2_patient_data@meta.data = tet2_patient_data@meta.data %>%
  rownames_to_column("UMI") %>%
  left_join(singlet_data@meta.data) %>%
  column_to_rownames("UMI")

Idents(tet2_patient_data) = "DF.classifications_0.25_0.09_4434"
tet2_patient_singlet_data = subset(tet2_patient_data, idents = "Singlet")
saveRDS(tet2_patient_singlet_data, "rds_objects/tet2_patient_singlet_data.rds")

dnmt3a_patient_data@meta.data = dnmt3a_patient_data@meta.data %>%
  rownames_to_column("UMI") %>%
  left_join(singlet_data@meta.data) %>%
  column_to_rownames("UMI")

Idents(dnmt3a_patient_data) = "DF.classifications_0.25_0.09_4434"
dnmt3a_patient_singlet_data = subset(dnmt3a_patient_data, idents = "Singlet")
saveRDS(dnmt3a_patient_singlet_data, "rds_objects/dnmt3a_patient_singlet_data.rds")

