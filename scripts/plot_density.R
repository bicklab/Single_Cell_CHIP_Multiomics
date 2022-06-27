
# Density plots -------------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

# All patients plot ---------------------------------------------------------

output_dir = "figures/figure_1/"

royal = wes_palette("Royal1", 22, "continuous")
(all_data_plot = DimPlot(object = singlet_data, reduction = "proj.umap", group.by = "collapsed.celltype", label = T, label.size = 5, repel = T, raster = F, cols = royal, pt.size = 0.1) +
  NoLegend() +
  ggtitle("Baseline UMAP") +
  guides(color = guide_legend(override.aes = list(size=6), ncol=1)) +
  NoAxes() +
  annotate(x=-16, xend=-16, y=-16, yend=-10, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate(x=-16, xend=-12, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))) +
  annotate("text", label="UMAP 1", x=-16.5, -14, size = 4, angle = 90) +
  annotate("text", label="UMAP 2", x=-15, y=-17, size = 4))
filename = paste0(output_dir, "umap_all_data")
save_plot(filename, all_data_plot, width = 8)


overlaid_data = singlet_data@meta.data %>% 
  group_by(MUTATION.GROUP) %>%
  mutate(MUTATION.LABEL = MUTATION.GROUP)

counts = overlaid_data %>%
  group_by(MUTATION.GROUP) %>%
  summarise(count = n()) %>%
  dplyr::select(count, MUTATION.GROUP)

counts_vector = deframe(counts)

overlaid_data$MUTATION.LABEL[overlaid_data$MUTATION.LABEL == counts_vector[1]] = paste0(counts_vector[1], " (n = ", names(counts_vector[1]), ")")
overlaid_data$MUTATION.LABEL[overlaid_data$MUTATION.LABEL == counts_vector[2]] = paste0(counts_vector[2], " (n = ", names(counts_vector[2]), ")")
overlaid_data$MUTATION.LABEL[overlaid_data$MUTATION.LABEL == counts_vector[3]] = paste0(counts_vector[3], " (n = ", names(counts_vector[3]), ")")


(overlaid_plot = ggplot(overlaid_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(color = wes_palette("Royal1", 20, "continuous")[15], size = 0.1) +
    stat_density2d(aes(fill = ..level.., alpha = ..level..), bins = 20, geom = "polygon") +
    format_add_arrows +
    labs(fill = "Density") +
    facet_wrap(~MUTATION.LABEL))

filename = paste0(output_dir, "density_all_patients")
save_plot(filename, overlaid_plot, width = 24)
write_tsv(overlaid_data, "tables/figure_1/density_all_patients.tsv")

# Individual patient plot ---------------------------------------------------

output_dir = "figures/figure_3/"

tet2_patient_singlet_metadata = tet2_patient_singlet_data@meta.data
tet2_patient_singlet_metadata$Clone[tet2_patient_singlet_metadata$Clone == "Mutant"] = "TET2"
# (tet2_patient_plot = ggplot(tet2_patient_singlet_metadata, aes(x = UMAP_1, y = UMAP_2)) +
#   geom_point(color = wes_palette("Royal1", 20, "continuous")[15]) +
#   stat_density2d(data = ~ subset(., Clone == "Wildtype"), aes(fill = ..level.., alpha = ..level..), bins = 40, geom = "polygon") +
#   format_add_arrows_2 +
#     labs(fill = "Wildtype density") +
#   new_scale_fill() +
#   stat_density2d(data = ~ subset(., Clone == "TET2"), aes(x = UMAP_1, y = UMAP_2, fill = ..level.., alpha = ..level..), bins = 40, geom = "polygon") +
#   format_add_arrows +
#     labs(fill = "Mutant density") +
#     ggtitle("TET2") +
#     theme(plot.title = element_text(size = 30)))

(tet2_patient_plot = ggplot(tet2_patient_singlet_metadata, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = wes_palette("Royal1", 20, "continuous")[15]) +
  stat_density2d(aes(fill = ..level.., alpha = ..level..), bins = 100, geom = "polygon") +
  format_add_arrows +
  labs(fill = "Density") +
  facet_wrap(~Clone))

filename = paste0(output_dir, "density_tet2_patient")
save_plot(filename, tet2_patient_plot, width = 12)
write_tsv(tet2_patient_singlet_metadata, "tables/figure_3/density_tet2_patient.tsv")

dnmt3a_patient_singlet_metadata = dnmt3a_patient_singlet_data@meta.data
dnmt3a_patient_singlet_metadata$Clone[dnmt3a_patient_singlet_metadata$Clone == "Mutant"] = "DNMT3A"
(dnmt3a_patient_plot = ggplot(dnmt3a_patient_singlet_metadata, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = wes_palette("Royal1", 20, "continuous")[15]) +
  stat_density2d(aes(fill = ..level.., alpha = ..level..), bins = 10, geom = "polygon") +
  format_add_arrows +
  labs(fill = "Density") +
  facet_wrap(~Clone))

filename = paste0(output_dir, "density_dnmt3a_patient")
save_plot(filename, dnmt3a_patient_plot, width = 12)
write_tsv(dnmt3a_patient_singlet_metadata, "tables/figure_3/density_dnmt3a_patient.tsv")



  
