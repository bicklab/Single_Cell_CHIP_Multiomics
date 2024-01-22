library(Seurat)
library(tidyverse)

setwd("~/Alyssa/Single_Cell_CHIP/")
diet_seurat_trimmed = readRDS("rds_objects/diet_seurat_trimmed.rds")
data_7079 = subset(diet_seurat_trimmed, subset = ProjectID == 7079)

count_table = data_7079$nCount_RNA %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(cell = substr(cell, 5, 22))
colnames(count_table)[2] = "id"


filtered_count_table = count_table[!duplicated(count_table$cell) & !duplicated(count_table$cell, fromLast=TRUE),]
write_csv(filtered_count_table, "tables/filtered_count_table_7079.csv")
