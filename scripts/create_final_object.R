#Brett Heimlich 

#Take the cell annotations identified from the harmony integrated assay and add them to the SCT integrated object. 

library(Seurat)
library(tidyverse)

seu.int <- readRDS("~/Brett/seurat_objects/old/Previous-versions/integrated_seurat.rds")

seu.name <- read_tsv("~/Brett/seurat_objects/integrated_seurat_metadata.tsv", col_names = T)
seu.name <- as.data.frame(seu.name)
head(seu.name)

cells.use <- seu.name$cell_ID
#subset based on cell ids in the seu.name df 
seu <- subset(seu.int, cells = cells.use)
seu 

#Add cell annotations and UMAP embeddings 
row.names(seu.name) <- seu.name$cell_ID
head(seu.name)
seu <- AddMetaData(object = seu, metadata = seu.name)
head(seu@meta.data)

#import final object 
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/filtered_annotated_seurat.rds")

#import PANN data 
pann <- read_tsv(file = "~/Brett/CH_Multiomics/seurat_objects/pANN_frame.tsv", col_names = T)
rownames(pann) <- pann$cell_ID
pann$cell_ID <- NULL
#add back PANN data
seu <- AddMetaData(object = seu, metadata = pann)
head(seu@meta.data)

#remove known doublet cells 
cells_to_remove = seu@meta.data %>%
  filter(LANE == "7079-3") %>%
  filter(pANN >= 0.28) %>%
  pull(cell_ID)

seu = seu[,!colnames(seu) %in% cells_to_remove]

saveRDS(seu, file = "~/Brett/CH_Multiomics/seurat_objects/final_annotated_seurat_correctedDR_5.11.23.rds")
