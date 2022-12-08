
# Density plots -------------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

# All patients plot ---------------------------------------------------------

output_dir = "figures/figure_1/"


Idents(singlet_data) = "celltype.de"
relevant_cell_type_data = subset(singlet_data, idents = c("CD4 T cell", "CD14 Mono", "CD8 T cell", "NK", "B", "CD16 Mono", "DC", "Platelet"))

set.seed(123)
royal = wes_palette("Royal1", 10, "continuous")[c(1,2,3,4,5,6,9,10)] %>% sample(8)


(all_data_plot = DimPlot(object = relevant_cell_type_data, reduction = "proj.umap", group.by = "celltype.de", label = FALSE, repel = T, raster = F, cols = royal, pt.size = 0.1) +
  NoLegend() +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size=6), ncol=1)) +
  theme(plot.title = element_blank()) +
  annotate(x=-16, xend=-16, y=-16, yend=-8, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate(x=-16, xend=-11, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate("text", label="UMAP 1", x=-17, -13, size = 4, angle = 90) +
  annotate("text", label="UMAP 2", x=-14, y=-18, size = 4))

filename = paste0(output_dir, "umap_all_data")
save_plot(filename, all_data_plot, height = 5, width = 5, png = TRUE)


# Density plot by genotype -----------------------------------------------


(overlaid_plot = ggplot(singlet_data@meta.data, aes(x = UMAP_1, y = UMAP_2)) +
   geom_point(fill = "#bdbdbd", color = "white", size = 1, shape = 21, stroke = 0.05) +
   stat_density2d(aes(alpha = ..level..), fill = wes_palette("Royal1", 20, "continuous")[7], bins = 10, geom = "polygon") +
   format_add_arrows +
   labs(fill = "Density") +
   facet_wrap(~MUTATION.GROUP))

filename = paste0(output_dir, "density_all_patients")
save_plot(filename, overlaid_plot, height = 5, width = 15, png = TRUE)
write_tsv(singlet_data@meta.data, "tables/figure_1/density_all_patients.tsv")

# Individual patient plot ---------------------------------------------------

output_dir = "figures/figure_3/"

tet2_patient_singlet_metadata = tet2_patient_singlet_data@meta.data
tet2_patient_singlet_metadata$Clone[tet2_patient_singlet_metadata$Clone == "Mutant"] = "TET2"

(tet2_patient_plot = ggplot(tet2_patient_singlet_metadata, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(fill = "#bdbdbd", color = "white", shape = 21, stroke = 0.05) +
    stat_density2d(aes(alpha = ..level..), fill = wes_palette("Royal1", 20, "continuous")[7], bins = 100, geom = "polygon") +
    format_add_arrows +
  labs(fill = "Density") +
  facet_wrap(~Clone))

filename = paste0(output_dir, "density_tet2_patient")
save_plot(filename, tet2_patient_plot, height = 4, width = 8, png = TRUE)
write_tsv(tet2_patient_singlet_metadata, "tables/figure_3/density_tet2_patient.tsv")

dnmt3a_patient_singlet_metadata = dnmt3a_patient_singlet_data@meta.data
dnmt3a_patient_singlet_metadata$Clone[dnmt3a_patient_singlet_metadata$Clone == "Mutant"] = "DNMT3A"
(dnmt3a_patient_plot = ggplot(dnmt3a_patient_singlet_metadata, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(fill = "#bdbdbd", color = "white", shape = 21, stroke = 0.05) +
    stat_density2d(aes(alpha = ..level..), fill = wes_palette("Royal1", 20, "continuous")[7], bins = 6, geom = "polygon") +
    format_add_arrows +
  labs(fill = "Density") +
  facet_wrap(~Clone))

filename = paste0(output_dir, "density_dnmt3a_patient")
save_plot(filename, dnmt3a_patient_plot, height = 4, width = 8, png = TRUE)
write_tsv(dnmt3a_patient_singlet_metadata, "tables/figure_3/density_dnmt3a_patient.tsv")



  
