
# Radar plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


summarize_data = function(data, grouping_category) {
  ignore_cells = c("Eryth", "HSPC", "Plasmablast", "MAIT")
  summary_stats = data@meta.data %>%
    filter(!(celltype.de %in% ignore_cells)) %>%
    group_by(across(all_of(grouping_category)), celltype.de) %>%
    summarize(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    mutate(celltype.de = factor(celltype.de, levels = c("CD4 T cell", "CD14 Mono", "CD8 T cell", "NK", 
                                                        "B", "CD16 Mono", 
                                                        "CD4 Naive", "CD8 Naive", "DC", "Platelet"))) %>%
    arrange(celltype.de) %>%
    dplyr::select(-count) %>%
    na.omit() %>%
    pivot_wider(names_from = celltype.de, values_from = proportion)
  return(summary_stats)
}

# All patients ------------------------------------------------------------

output_dir = "figures/figure_1/"

summary_stats = summarize_data(singlet_data, "MUTATION.GROUP")

radar_stats = summary_stats %>%
  mutate(MUTATION.GROUP = factor(MUTATION.GROUP, levels = c("Control", "TET2", "DNMT3A")))

(radar_plot = ggradar(radar_stats, 
                values.radar = c(0, 0.2,0.5), 
                grid.min = 0, grid.mid = 0.2, grid.max = 0.5,
                group.line.width = 1,
                group.point.size = 3,
                group.colours = wes_palette("Royal1", 20, "continuous")[c(1,7,18)],
                background.circle.colour = "white",
                gridline.mid.colour = "grey",
                legend.position = "bottom"))

filename = paste0(output_dir, "radar_myeloid_skew")
save_plot(filename, radar_plot)
write_tsv(radar_stats, "tables/figure_1/radar_myeloid_skew.tsv")


# TET2 patient ------------------------------------------------------

output_dir = "figures/figure_3/"

summary_stats_tet2 = summarize_data(tet2_patient_singlet_data, "Clone")
summary_stats_tet2[is.na(summary_stats_tet2)] = 0
summary_stats_tet2[3,] = summary_stats[1,c(1:6, 10:11)]

tet2_patient_radar_stats = summary_stats_tet2 %>%
  mutate(Clone = factor(Clone, levels = c("Wildtype", "Mutant", "Control")))

(tet2_patient_radar = ggradar(tet2_patient_radar_stats, 
                values.radar = c(0, 0.5, 1), 
                group.line.width = 1,
                group.point.size = 3,
                group.colours = wes_palette("Royal1", 20, "continuous")[c(1,7,18)],
                background.circle.colour = "white",
                gridline.mid.colour = "grey",
                legend.position = "bottom",
                plot.title = tet2_patient_singlet_data$MUTATION.GROUP[1]))

filename = paste0(output_dir, "radar_tet2_patient")
save_plot(filename, tet2_patient_radar)
write_tsv(tet2_patient_radar_stats, "tables/figure_3/radar_tet2_patient.tsv")


# DNMT3A patient ------------------------------------------------------

output_dir = "figures/figure_3/"

summary_stats_dnmt3a = summarize_data(dnmt3a_patient_singlet_data, "Clone")
summary_stats_dnmt3a[is.na(summary_stats_dnmt3a)] = 0
summary_stats_dnmt3a[3,] = summary_stats[1,c(1:7, 10:11)]

dnmt3a_patient_radar_stats = summary_stats_dnmt3a %>%
  mutate(Clone = factor(Clone, levels = c("Wildtype", "Mutant", "Control")))

(dnmt3a_patient_radar = ggradar(dnmt3a_patient_radar_stats, 
                              values.radar = c(0, 0.5, 1), 
                              group.line.width = 1,
                              group.point.size = 3,
                              group.colours = wes_palette("Royal1", 20, "continuous")[c(1,7,18)],
                              background.circle.colour = "white",
                              gridline.mid.colour = "grey",
                              legend.position = "bottom",
                              plot.title = dnmt3a_patient_singlet_data$MUTATION.GROUP[1]))

filename = paste0(output_dir, "radar_dnmt3a_patient")
save_plot(filename, dnmt3a_patient_radar)
write_tsv(dnmt3a_patient_radar_stats, "tables/figure_3/radar_dnmt3a_patient.tsv")
