# Prep object for CZI


# SAVE SIMPLIFIED SEURAT OBJECT --------------------------------------------------------

library(Seurat)
library(tidyverse)

diet_seurat_trimmed = readRDS("~/Alyssa/Single_Cell_CHIP/rds_objects/diet_seurat_trimmed.rds")
chip_pbmc_for_czi = diet_seurat_trimmed %>% subset(subset = STIM == "VEH")

chip_pbmc_for_czi$cell_type[chip_pbmc_for_czi$cell_type == "CD14 Monos"] = "CD14+ Monocytes"
chip_pbmc_for_czi$cell_type[chip_pbmc_for_czi$cell_type == "CD16 Monos"] = "CD16+ Monocytes"
chip_pbmc_for_czi$cell_type[chip_pbmc_for_czi$cell_type == "B"] = "B cells"
chip_pbmc_for_czi$cell_type[chip_pbmc_for_czi$cell_type == "NK"] = "NK cells"

chip_pbmc_for_czi@meta.data = chip_pbmc_for_czi@meta.data %>% dplyr::select(-STIM)
saveRDS(chip_pbmc_for_czi, "~/Alyssa/Single_Cell_CHIP/rds_objects/chip_pbmc_for_czi.rds")


# SAVE ETHNICITY METADATA --------------------------------------------------------------

library(tidyverse)

metadata = readxl::read_excel("~/Alyssa/Single_Cell_CHIP/input_data/clinical data for brett v5.xlsx")
chip_pbmc_for_czi = readRDS("rds_objects/chip_pbmc_for_czi.rds")

metadata_for_czi = metadata %>%
  filter(patient_id != "CH-21-045") %>%
  select(patient_id, race, ethnicity) %>%
  rename(CHIVEID = patient_id)

write_tsv(metadata_for_czi, "~/Alyssa/Single_Cell_CHIP/tables/metadata_for_czi.tsv")

# GENERATE H5AD FORM FOR BRIAN FOR HIS PRIVATE SERVER --------------------------------

library(sceasy)
library(reticulate)
setwd("/home/rstudio/Alyssa/Single_Cell_CHIP")
sceasy::convertFormat(chip_pbmc_for_czi, from="seurat", to="anndata",
                      outFile='rds_objects/chip_pbmc_for_czi.h5ad')


library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)

chip_pbmc_for_czi = readRDS("rds_objects/chip_pbmc_for_czi.rds")
seurat_with_umap = readRDS("rds_objects/filtered_annotated_seurat.rds")

seurat_with_umap = seurat_with_umap %>% subset(subset = STIM == "VEH")

seurat_with_umap$cell_type[seurat_with_umap$cell_type == "CD14 Monos"] = "CD14+ Monocytes"
seurat_with_umap$cell_type[seurat_with_umap$cell_type == "CD16 Monos"] = "CD16+ Monocytes"
seurat_with_umap$cell_type[seurat_with_umap$cell_type == "B"] = "B cells"
seurat_with_umap$cell_type[seurat_with_umap$cell_type == "NK"] = "NK cells"

seurat_with_umap@meta.data = seurat_with_umap@meta.data %>% dplyr::select(-STIM)

SaveH5Seurat(seurat_with_umap, filename = "rds_objects/seurat_with_umap.h5Seurat")
Convert("rds_objects/seurat_with_umap.h5Seurat", dest = "h5ad") 
system("gsutil -m cp rds_objects/seurat_with_umap.h5ad gs://fc-112f611a-aca5-42eb-9970-5050086b3e8d/Single_Cell_CHIP/rds_objects/")


system("gsutil cp gs://bicklab-main-storage/Users/Brian_Sharber/temp/pbmc3k.h5ad rds_objects/")
Convert("rds_objects/pbmc3k.h5ad", "rds_objects/pbmc3k.h5Seurat")
example_h5ad = LoadH5Seurat("rds_objects/pbmc3k.h5Seurat")

# seurat_with_umap actually only has UMAP in the metadata - we need it in the graphs for the cellxgene viewer
chip_pbmc_for_czi = readRDS("rds_objects/chip_pbmc_for_czi.rds")

chip_pbmc_with_umap = SCTransform(chip_pbmc_for_czi)
chip_pbmc_with_umap = RunPCA(chip_pbmc_for_czi)
chip_pbmc_with_umap = RunUMAP(chip_pbmc_for_czi)




