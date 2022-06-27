
# Scatterplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

# set variables
cell_type = "CD14_Mono"
folder = "tables/supplemental/"

gene_label_data = read.xlsx("input_data/gene_listp_61322.xlsx", sheet = "CD14 Mono")

curate_data = function(exp_data, cell_type = "CD14_Mono") {
  if (cell_type == "CD14_Mono") {
    exp_data %>% 
      dplyr::mutate(Pathway = factor(case_when(Gene %in% gene_label_data$Inflammation ~ "Inflammation",
                                                  Gene %in% gene_label_data$Antigen.Presentation ~ "Antigen Presentation", 
                                                  Gene %in% gene_label_data$Cell.Adhesion ~ "Cell Adhesion",
                                                  Gene %in% gene_label_data$Monocyte.Activation ~ "Monocyte Activation",
                                                  TRUE ~ "Other" ),
                                        levels = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation", "Other"))) %>% 
      arrange(desc(Pathway)) %>%
      return()
  } else if (cell_type == "CD8_T_cell") {
    exp_data %>% 
      dplyr::mutate(Pathway = factor(case_when(Gene %in% gene_label_data_t_cell$T.Cell.Activation ~ "T Cell Activation",
                                                  Gene %in% gene_label_data_t_cell$Intracellular.Signaling ~ "Intracellular Signaling", 
                                                  Gene %in% gene_label_data_t_cell$Immunomodulation ~ "Immunomodulation",
                                                  TRUE ~ "Other" ),
                                        levels = c("T Cell Activation", "Intracellular Signaling", "Immunomodulation", "Other"))) %>% 
      arrange(desc(Pathway)) %>%
      return()
  }
}

for (stim_status in c("unstim", "stim")) {
  # read in log fc results from tet2 and dnmt3a
  tet2_input_filename = list.files(folder, pattern = paste0("^TET2_", cell_type, ".*_", stim_status))
  tet2_results = read_csv(paste0(folder, tet2_input_filename)) %>%
    dplyr::select(Gene, avg_log2FC) %>%
    dplyr::rename(TET2_log2FC = avg_log2FC)
  
  dnmt3a_input_filename = list.files(folder, pattern = paste0("^DNMT3A_", cell_type, ".*_", stim_status))
  dnmt3a_results = read_csv(paste0(folder, dnmt3a_input_filename)) %>%
    dplyr::select(Gene, avg_log2FC) %>%
    dplyr::rename(DNMT3A_log2FC = avg_log2FC)
  
  # combine based on gene
  combined_fc = inner_join(tet2_results, dnmt3a_results)
  combined_fc_labeled = curate_data(combined_fc)
  
  # use to make plot
  (scatter_plot = ggplot(combined_fc_labeled, aes(x = TET2_log2FC, y = DNMT3A_log2FC)) +
      geom_point(aes(color = Pathway)) +
      theme_bw() +
      scale_color_manual(values = colors) +
      geom_smooth(method = 'lm', color = "darkgrey") +
      theme(panel.border = element_blank(),
            panel.grid.minor = element_blank()) +
      ggtitle(paste(gsub("_", " ", cell_type), stim_status)))
  
  # save plot and data
  output_dir = "figures/supplemental/"
  filename = paste0(output_dir, "scatter_tet2_vs_dnmt3a_", cell_type, "_", stim_status)
  save_plot(filename, scatter_plot, height = 5, width = 6)
  write_tsv(combined_fc_labeled, paste0(gsub("figures", "tables", filename), ".tsv"))
}


