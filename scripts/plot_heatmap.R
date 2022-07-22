
# Heatmap plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

cell_types_of_interest = c("CD4 T cell", "CD14 Mono", "CD8 T cell", "NK", "B", "CD16 Mono", "DC", "Platelet")
num_cell_per_type = 500
num_top_genes = 10
Idents(singlet_data) = "celltype.de"

set.seed(123)
heatmap_data = subset(singlet_data, idents = cell_types_of_interest, downsample = 500)

pbmc.markers = FindAllMarkers(heatmap_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

levels(heatmap_data) = c("NK", "CD8 T cell", "CD14 Mono", "CD4 T cell",  "B", "CD16 Mono", "Platelet", "DC")

(heatmap = DoHeatmap(heatmap_data, features = top10$gene, raster = FALSE) + 
  guides(color = "none"))

filename = "figures/supplemental/cell_identity_heatmap"
save_plot(filename, heatmap, width = 7, height = 13)
