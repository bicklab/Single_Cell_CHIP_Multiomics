# Prep object for CZI

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
