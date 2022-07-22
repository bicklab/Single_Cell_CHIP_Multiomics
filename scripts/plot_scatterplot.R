
# Scatterplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

# set variables
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

make_scatterplot = function(cell_type) {
  plots = list()
  for (stim_status in c("unstim", "stim")) {
    # read in log fc results from tet2 and dnmt3a
    tet2_input_filename = list.files(folder, pattern = paste0("^TET2_", cell_type, ".*_", stim_status, "_wilcox"))
    print(tet2_input_filename)
    tet2_results = read_csv(paste0(folder, tet2_input_filename)) %>%
      dplyr::select(Gene, avg_log2FC) %>%
      dplyr::rename(TET2_log2FC = avg_log2FC)
    
    dnmt3a_input_filename = list.files(folder, pattern = paste0("^DNMT3A_", cell_type, ".*_", stim_status, "_wilcox"))
    print(dnmt3a_input_filename)
    
    dnmt3a_results = read_csv(paste0(folder, dnmt3a_input_filename)) %>%
      dplyr::select(Gene, avg_log2FC) %>%
      dplyr::rename(DNMT3A_log2FC = avg_log2FC)
    
    # combine based on gene
    combined_fc = inner_join(tet2_results, dnmt3a_results)
    combined_fc_labeled = curate_data(combined_fc)
    
    # use to make plot
    (scatter_plot = ggplot(combined_fc_labeled, aes(x = TET2_log2FC, y = DNMT3A_log2FC)) +
        geom_point(aes(color = Pathway)) +
        xlim(-2, 2) +
        ylim(-2, 2) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        scale_color_manual(values = colors) +
        geom_smooth(method = 'lm', color = "darkgrey") +
        theme(panel.border = element_blank(),
              panel.grid.minor = element_blank()) +
        ggtitle(paste(gsub("_", " ", cell_type), stim_status)))
    
    # save plot and data
    output_dir = "figures/supplemental/"
    filename = paste0(output_dir, "scatter_tet2_vs_dnmt3a_", cell_type, "_", stim_status)
    #save_plot(filename, scatter_plot, height = 5, width = 6)
    
    combined_fc_labeled = combined_fc_labeled %>%
      mutate(log2_FC_diff_tet2_minus_dnmt3a = TET2_log2FC - DNMT3A_log2FC) %>%
      arrange(log2_FC_diff_tet2_minus_dnmt3a)
    
    interesting_results = rbind(head(combined_fc_labeled, 50), tail(combined_fc_labeled, 50))
    write_csv(interesting_results, paste0(gsub("figures", "tables", filename), ".csv"))
    plots[[length(plots)+1]] = plot(scatter_plot)  
  }
  return(plots)
}

cell_types = c("CD14_Mono", "CD8_T_cell", "CD4_T_cell")
cd14_scatter = make_scatterplot("CD14_Mono")
cd8_scatter = make_scatterplot("CD8_T_cell")
cd4_scatter = make_scatterplot("CD4_T_cell")


combined_scatter = ggarrange(plotlist = c(cd14_scatter, cd8_scatter, cd4_scatter), nrow = 3, ncol = 2)
save_plot(paste0(output_dir, "scatter_combined"), combined_scatter, height = 12, width = 10)

# Compile excel file ------------------------------------------------------

sample_genotypes = c("DNMT3A", "TET2")
folder = "tables/supplemental/"
output_dir = "tables/supplemental/"

get_cell_data = function(sample_genotype, cell_type) {
  input_filename = list.files(folder, pattern = paste0("scatter_tet2_vs_dnmt3a_", cell_type, "_unstim"))
  data = read_csv(paste0(folder, input_filename))
  return(data)
}
  
cd8_t_cell_data = get_cell_data(sample_genotype, "CD8_T_cell")
cd14_mono_data = get_cell_data(sample_genotype, "CD14_Mono")
cd4_t_cell_data = get_cell_data(sample_genotype, "CD4_T_cell")

dataframes = list("CD14_Mono" = cd14_mono_data, "CD4_T_cell" = cd4_t_cell_data, "CD8_T_cell" = cd8_t_cell_data)

filename = paste0(output_dir, "tet2_v_dnmt3a_unstim.xlsx")
write.xlsx(dataframes, file = filename)




# VAF and myeloid scatter ----------------------------------------------------

vaf_and_myeloid = read_tsv("tables/supplemental/vaf_and_myeloid.tsv")

vaf_myeloid_plot = ggplot(vaf_and_myeloid, aes(x = VAF, y = cd14_mono_prop)) +
  geom_smooth(method = "lm") +
  geom_point() +
  ylab("CD14 monocyte proportion") +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

save_plot("figures/supplemental/vaf_and_myeloid", vaf_myeloid_plot, height = 4, width = 4)
