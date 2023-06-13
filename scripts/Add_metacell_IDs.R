#Add meta-cell IDs to seurat object 

import seurat object
library(Seurat)
library(tidyverse)

#import seurat object(s)
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/final_normscalevehonly.rds")

#import csvs w/IDs 
tetveh <- read_csv("~/Brett/CH_Multiomics/seurat_objects/meta_cells/metacell_IDs/tet2_veh_metacells.csv", col_names = T)
noneveh <- read_csv("~/Brett/CH_Multiomics/seurat_objects/meta_cells/metacell_IDs/dnmt3a_veh_metacells.csv", col_names = T)
dnmt3aveh <- read_csv("~/Brett/CH_Multiomics/seurat_objects/meta_cells/metacell_IDs/none_veh_metacells.csv", col_names = T)

#add unique identifier to each
tetveh$metacell <- paste0(tetveh$metacell, "T")
dnmt3aveh$metacell <- paste0(dnmt3aveh$metacell, "D")
noneveh$metacell <- paste0(noneveh$metacell, "N")

#combine IDs
meta_ids <- rbind(tetveh, noneveh, dnmt3aveh)
meta_ids <- as.data.frame(meta_ids)

#rename
colnames(meta_ids) <- c("cell_id", "metacellID")
head(meta_ids)

#Add cell annotations and UMAP embeddings 
row.names(meta_ids) <- meta_ids$cell_id
head(meta_ids)
seu <- AddMetaData(object = seu, metadata = meta_ids)
head(seu@meta.data)

#visualize 
table(seu@meta.data$metacellID, seu@meta.data$cell_type)

seu@meta.data %>% 
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = metacellID)) +
  facet_wrap(~cell_type) +
  geom_point(size = 0.5)



  #save 
saveRDS(seu, file = "~/Brett/seurat_objects/final_annotated_object.rds")