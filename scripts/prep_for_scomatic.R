# prep for scomatic

library(tidyverse)

diet_seurat_trimmed = readRDS("~/Alyssa/Single_Cell_CHIP/rds_objects/diet_seurat_trimmed.rds")

# SComatic relies on cell barcodes from bam files, which must be unique (no doubles and no doublets)
good_cells = diet_seurat_trimmed@meta.data %>%
  rownames_to_column("Index") %>%
  dplyr::select(Index, HTO_maxID) %>% 
  dplyr::mutate(Index = gsub(".*_", "", Index)) %>%
  dplyr::group_by(Index) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::filter(Count == 1) %>%
  dplyr::select(Index)

celltype_barcode_map = diet_seurat_trimmed@meta.data %>%
  rownames_to_column("Index") %>%
  dplyr::select(Index, cell_type) %>% 
  dplyr::mutate(Index = gsub(".*_", "", Index)) %>%
  dplyr::inner_join(good_cells) %>%
  dplyr::rename("Cell_type" = cell_type) %>%
  dplyr::select(Index, Cell_type)

write_tsv(celltype_barcode_map, "~/Alyssa/Single_Cell_CHIP/tables/scomatic/chip_pbmc_celltype_barcode_map.tsv")


seurat_with_clones_PB_10.16.23 = readRDS("~/Alyssa/SC_CHIP_HIV_Diabetes/input_data/seurat_with_clones_PB_10.16.23.rds")

good_cells = seurat_with_clones_PB_10.16.23@meta.data %>%
  rownames_to_column("Index") %>%
  dplyr::select(Index, HTO_maxID) %>% 
  dplyr::mutate(Index = gsub(".*_", "", Index)) %>%
  dplyr::group_by(Index) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::filter(Count == 1) %>%
  dplyr::select(Index)

clone_map = seurat_with_clones_PB_10.16.23@meta.data %>%
  rownames_to_column("Index") %>%
  dplyr::mutate(Index = gsub(".*_", "", Index)) %>%
  dplyr::inner_join(good_cells) %>%
  dplyr::rename("Cell_type" = cell_type) %>%
  dplyr::select(Index, Cell_type, Clone)

write_tsv(clone_map, "~/Alyssa/Single_Cell_CHIP/tables/scomatic/clone_map.tsv")



