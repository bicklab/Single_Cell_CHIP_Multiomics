library(Seurat)
library(tidyverse)
library(org.Hs.eg.db)
library(SingleCellExperiment)


# Convert to loom format

setwd("~/Alyssa/Single_Cell_CHIP/")
diet_seurat_trimmed = readRDS("rds_objects/diet_seurat_trimmed.rds")

diet_seurat_trimmed_gene_names = mapIds(
  org.Hs.eg.db,
  keys = rownames(diet_seurat_trimmed),
  column = 'ENSEMBL',
  keytype = 'SYMBOL')

diet_seurat_trimmed = as.SingleCellExperiment(diet_seurat_trimmed)
rownames(diet_seurat_trimmed) = diet_seurat_trimmed_gene_names
diet_seurat_trimmed = diet_seurat_trimmed[!is.na(rownames(diet_seurat_trimmed)),]
rownames(diet_seurat_trimmed) = make.unique(rownames(diet_seurat_trimmed)) 

diet_seurat_trimmed_ensembl = as.Seurat(diet_seurat_trimmed)

SeuratDisk::SaveLoom(diet_seurat_trimmed_ensembl, "rds_objects/geneformer/diet_seurat_trimmed_ensembl.loom")

Idents(diet_seurat_trimmed_ensembl) = "CHIP"
split_pbmc_seurat = SplitObject(diet_seurat_trimmed_ensembl)
SeuratDisk::SaveLoom(split_pbmc_seurat$chip, "rds_objects/geneformer/chip_pbmc_ensembl.loom")
SeuratDisk::SaveLoom(split_pbmc_seurat$control, "rds_objects/geneformer/control_pbmc_ensembl.loom")


# Pull HSC data from CZI CELLxGENE
library(cellxgene.census)


census = open_soma(census_version = "2023-12-15")

# Open obs SOMADataFrame
cell_metadata <-  census$get("census_data")$get("homo_sapiens")$get("obs")

# Read as Arrow Table
cell_metadata <-  cell_metadata$read(
  value_filter = "cell_type %in% c('hematopoietic stem cell') & assay %in% c('10x 3\\' v3') & disease %in% c('normal')",
  column_names = c("assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"),
)

cell_metadata <-  cell_metadata$concat()
cell_metadata <-  as.data.frame(cell_metadata)

hsc_seurat = get_seurat(
  census = census,
  organism = "Homo sapiens",
  obs_value_filter = "cell_type %in% c('hematopoietic stem cell') & assay %in% c('10x 3\\' v3') & disease %in% c('normal') & tissue %in% c('blood', 'bone marrow')",
  obs_column_names = c("assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease")
)

rownames(seurat_obj_2@assays$RNA@counts) = seurat_obj_2@assays$RNA@meta.features$feature_name
rownames(seurat_obj_2@assays$RNA@data) = seurat_obj_2@assays$RNA@meta.features$feature_name

saveRDS(seurat_obj_2, "rds_objects/czi/endo_sample_2.rds")






