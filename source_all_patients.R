source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
# 
# all_data = readRDS("input_data/complete_filtered_az_june22.rds")
# 
# index = c("NK", "Treg", "CD14 Mono", "CD8 T cell", "CD4 T cell", "Eryth", "CD16 Mono", "DC", "CD8 T cell", "CD4 Naive",
#           "CD8 Naive", "HSPC", "MAIT", "CD8 T cell", "CD4 T cell", "B", "B", "Platelet", "NK", "B",
#           "DC", "CD4 CTL", "NK", "CD4 T cell", "Plasmablast",  "NK", "dnT", "DC", "CD8 T cell")
# values = c("NK", "CD4 T cell", "CD14 Mono", "CD8 T cell", "CD4 T cell", "Eryth", "CD16 Mono", "DC", "CD8 T cell", "CD4 T cell",
#            "CD8 T cell", "HSPC", "MAIT", "CD8 T cell", "CD4 T cell", "B", "B", "Platelet", "NK", "B",
#            "DC", "CD4 T cell", "NK", "CD4 T cell", "Plasmablast", "NK", "CD4 T cell", "DC", "CD8 T cell")
# 
# index_de = c("CD8 TCM", "CD8 TEM", "CD8 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Naive", "Treg", "CD14 Mono", "CD16 Mono",
#              "B", "Eryth", "NK", "DC", "HSPC", "MAIT", "Platelet", "Plasmablast", "dnT", "gdT",
#              "B naive", "NK_CD56bright", "B intermediate", "pDC", "B memory", "cDC2", "ILC", "ASDC", "CD4 Proliferating", "CD8 Proliferating",
#              "NK Proliferating")
# values_de = c("CD8 T cell", "CD8 T cell", "CD8 Naive", "CD4 T cell", "CD4 T cell", "CD4 CTL", "CD4 Naive", "Treg", "CD14 Mono", "CD16 Mono",
#               "B", "Eryth", "NK", "DC", "HSPC", "MAIT", "Platelet", "Plasmablast", "dnT", "gdT",
#               "B", "NK", "B", "DC", "B", "DC", "ILC", "DC", "CD4 Proliferating", "CD8 Proliferating",
#               "NK Proliferating")
# 
# combine_cell_types = function(data) {
#   data@meta.data$celltype.bh = values[match(data@meta.data$collapsed.celltype, index)]
# 
#   # make new column for de, pull from predicted celltype l2
#   data@meta.data$celltype.de = values_de[match(data@meta.data$predicted.celltype.l2, index_de)]
# 
#   return(data)
# }
# 
# all_data = combine_cell_types(all_data)
# all_data$MUTATION.GROUP[all_data$MUTATION.GROUP == "none"] = "Control"
# 
# # this patient has very low cell counts
# all_data = subset(all_data, orig.ident != "CH-21-020")
# # this value was missing from the data
# all_data@meta.data$MUTATION[all_data@meta.data$MUTATION == "DNMT3A R882C"] = "DNMT3A R882C (13.3%)"
# 
# # add UMAP groupings to metadata
# all_data$UMAP_1 = all_data@reductions$proj.umap@cell.embeddings[,1]
# all_data$UMAP_2 = all_data@reductions$proj.umap@cell.embeddings[,2]
# 
# options('future.globals.maxSize' = 9000*1024^2)
# future(saveRDS(all_data, "rds_objects/cleaned_complete_data.rds"))

all_data = readRDS("rds_objects/cleaned_complete_data.rds")

