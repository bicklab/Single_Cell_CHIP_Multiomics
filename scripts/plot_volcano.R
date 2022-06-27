
# Volcano plots -----------------------------------------------------------
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


# load gene classes
gene_df = read_csv("input_data/Gene_class.csv")
gene_label_data = read.xlsx("input_data/gene_listp_61322.xlsx", sheet = "CD14 Mono")
gene_label_data_t_cell = read.xlsx("input_data/gene_listp_61322.xlsx", sheet = "CD8 T")

# Helper functions -----------------------------------------------------------

curate_data = function(exp_data, cell_type = "CD14_Mono") {
  if (cell_type == "CD14_Mono") {
    exp_data %>% 
    dplyr::mutate(significant = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
    dplyr::mutate(gene_class = factor(case_when(significant & Gene %in% gene_label_data$Inflammation ~ "Inflammation",
                                                significant & Gene %in% gene_label_data$Antigen.Presentation ~ "Antigen Presentation", 
                                                significant & Gene %in% gene_label_data$Cell.Adhesion ~ "Cell Adhesion",
                                                significant & Gene %in% gene_label_data$Monocyte.Activation ~ "Monocyte Activation",
                                                TRUE ~ "Other" ),
                                      levels = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation", "Other"))) %>% 
    mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
    arrange(desc(gene_class)) %>%
    return()
  } else if (cell_type == "CD8_T_cell") {
    exp_data %>% 
      dplyr::mutate(significant = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
      dplyr::mutate(gene_class = factor(case_when(significant & Gene %in% gene_label_data_t_cell$T.Cell.Activation ~ "T Cell Activation",
                                                  significant & Gene %in% gene_label_data_t_cell$Intracellular.Signaling ~ "Intracellular Signaling", 
                                                  significant & Gene %in% gene_label_data_t_cell$Immunomodulation ~ "Immunomodulation",
                                                  TRUE ~ "Other" ),
                                        levels = c("T Cell Activation", "Intracellular Signaling", "Immunomodulation", "Other"))) %>% 
      mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
      arrange(desc(gene_class)) %>%
      return()
  }
}

get_lims = function(top, space_per_gene, n) {
  lower_limit = top - space_per_gene * n
  upper_limit = top
  return(c(lower_limit, upper_limit))
}

get_lims_t = function(bottom, space_per_gene, n) {
  lower_limit = bottom
  upper_limit = bottom + space_per_gene * n
  return(c(lower_limit, upper_limit))
}

gene_labels = function(data, direction, xlim, ylim) {
  sorted_data = data %>% 
    arrange(desc(`-log10_p_val`))
  if (direction == "pos") {
    xlims = c(xlim - 0.25, xlim + 0.25)
  } else {
    xlims = c(-xlim - 0.25, -xlim + 0.25)
  }
  labels = geom_label_repel(data = sorted_data, mapping = aes(label = Gene), size = 3, direction = "y", max.overlaps =  Inf,
                            segment.alpha = 0.25, force_pull = 0, xlim = xlims, ylim = ylim, show.legend = FALSE, box.padding = 0, 
                            segment.linetype = "solid", label.size = 0.08)
  return(list(labels))
}

make_volcano_plots_mono = function(cell_type, output_dir, table_dir, x_limit = 1.5, y_limit = 300, space_per_gene = 12.1, stim_status = "unstim") {
  sample_genotypes = c("TET2", "DNMT3A")
  plots = list()
  
  for (sample_genotype in sample_genotypes) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data, cell_type)
    
    # prep data subsets for labeling
    inflammation_data = subset(data, gene_class == "Inflammation")
    antigen_data = subset(data, gene_class == "Antigen Presentation")
    adhesion_data = subset(data, gene_class == "Cell Adhesion")
    activation_data = subset(data, gene_class == "Monocyte Activation")
    
    inflammation_data_pos = inflammation_data %>% filter(avg_log2FC > 0)
    antigen_data_pos = antigen_data %>% filter(avg_log2FC > 0)
    adhesion_data_pos = adhesion_data %>% filter(avg_log2FC > 0)
    activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
    
    inflammation_data_neg = inflammation_data %>% filter(avg_log2FC < 0)
    antigen_data_neg = antigen_data %>% filter(avg_log2FC < 0)
    adhesion_data_neg = adhesion_data %>% filter(avg_log2FC < 0)
    activation_data_neg = activation_data %>% filter(avg_log2FC < 0)

    #space_per_gene_pos = y_limit / (nrow(inflammation_data_pos) + nrow(antigen_data_pos) + nrow(adhesion_data_pos) + nrow(activation_data_pos))
    #space_per_gene_neg = y_limit / (nrow(inflammation_data_neg) + nrow(antigen_data_neg) + nrow(adhesion_data_neg) + nrow(activation_data_neg))
    
    legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano +
        xlim(-x_limit, x_limit) +
        ylim(0, y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
        
        gene_labels(data = inflammation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit, space_per_gene, nrow(inflammation_data_pos))) +
        gene_labels(data = antigen_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_pos)), space_per_gene, nrow(antigen_data_pos))) +
        gene_labels(data = adhesion_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_pos) + nrow(antigen_data_pos)), space_per_gene, nrow(adhesion_data_pos))) +
        gene_labels(data = activation_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_pos) + nrow(adhesion_data_pos) + nrow(antigen_data_pos)), space_per_gene, nrow(activation_data_pos))) +
        
        gene_labels(data = inflammation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit, space_per_gene, nrow(inflammation_data_neg))) +
        gene_labels(data = antigen_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_neg)), space_per_gene, nrow(antigen_data_neg))) +
        gene_labels(data = adhesion_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_neg) + nrow(antigen_data_neg)), space_per_gene, nrow(adhesion_data_neg))) +
        gene_labels(data = activation_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data_neg) + nrow(adhesion_data_neg) + nrow(antigen_data_neg)), space_per_gene, nrow(activation_data_neg))) +
        
        ggtitle(paste(gsub("_", " ", cell_type), sample_genotype, paste0(stim_status, "ulated"))) +
        theme(legend.position = ifelse(legend, "right", "none"), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
    
    filename = paste0(output_dir, "volcano_", sample_genotype, "_", cell_type, "_", stim_status)
    save_plot(filename, plot, legend, width = 9)
    write_tsv(data, paste0(table_dir, sample_genotype, "_", cell_type, "_", stim_status, ".tsv"))
    plots[[length(plots)+1]] = plot(plot)  
  }
  return(plots)
}

make_volcano_plots_t = function(cell_type, output_dir, table_dir, x_limit = 1.5, y_limit = 300, space_per_gene = 12.1, stim_status = "unstim", y_lower_bound = 100) {
  sample_genotypes = c("TET2", "DNMT3A")
  plots = list()
  
  for (sample_genotype in sample_genotypes) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data, cell_type)
    
    # prep data subsets for labeling
    activation_data = subset(data, gene_class == "T Cell Activation")
    intracellular_data = subset(data, gene_class == "Intracellular Signaling")
    immunomodulation_data = subset(data, gene_class == "Immunomodulation")
    
    activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
    intracellular_data_pos = intracellular_data %>% filter(avg_log2FC > 0)
    immunomodulation_data_pos = immunomodulation_data %>% filter(avg_log2FC > 0)
    
    activation_data_neg = activation_data %>% filter(avg_log2FC < 0)
    intracellular_data_neg = intracellular_data %>% filter(avg_log2FC < 0)
    immunomodulation_data_neg = immunomodulation_data %>% filter(avg_log2FC < 0)

    legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano_2 +
        xlim(-x_limit, x_limit) +
        ylim(0, y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
        
        gene_labels(data = activation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, nrow(activation_data_pos))) +
        gene_labels(data = intracellular_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * nrow(activation_data_pos), space_per_gene, nrow(intracellular_data_pos))) +
        gene_labels(data = immunomodulation_data_pos, direction = "pos",  xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (nrow(activation_data_pos) + nrow(intracellular_data_pos)), space_per_gene, nrow(immunomodulation_data_pos))) +
        
        gene_labels(data = activation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, nrow(activation_data_neg))) +
        gene_labels(data = intracellular_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * nrow(activation_data_neg), space_per_gene, nrow(intracellular_data_neg))) +
        gene_labels(data = immunomodulation_data_neg, direction = "neg",xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (nrow(activation_data_neg) + nrow(intracellular_data_neg)), space_per_gene, nrow(immunomodulation_data_neg))) +
        
        ggtitle(paste(gsub("_", " ", cell_type), sample_genotype, paste0(stim_status, "ulated"))) +
        theme(legend.position = ifelse(legend, "right", "none"), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
    
    filename = paste0(output_dir, "volcano_", sample_genotype, "_", cell_type, "_", stim_status)
    save_plot(filename, plot, legend, width = 9)
    write_tsv(data, paste0(table_dir, sample_genotype, "_", cell_type, "_", stim_status, ".tsv"))
    plots[[length(plots)+1]] = plot(plot)  
  }
  return(plots)
}


# All patients - unstimulated ---------------------------------------------

output_dir = "figures/figure_2/"

# for main figure
cd14_plots = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/")
cd8_plots = make_volcano_plots_t(cell_type = "CD8_T_cell", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/", y_limit = 250, x_limit = 1.75, space_per_gene = 9)

all_plots = c(cd14_plots, cd8_plots)
combined_plot = ggarrange(plotlist = all_plots,
                           nrow = 2, ncol = 2, widths = c(3, 4))

save_plot(paste0(output_dir, "volcano_combined"), combined_plot, height = 12.4, width = 16)

# # supplemental
# NK_plots = make_volcano_plots(cell_type = "NK", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", y_limit = 50)
# B_plots = make_volcano_plots(cell_type = "B", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", y_limit = 50)
# 
# all_plots = c(NK_plots, B_plots)
# combined_plot = ggarrange(plotlist = all_plots,
#                           nrow = 2, ncol = 2, widths = c(3, 4))
# 
# output_dir = "figures/supplemental/"
# save_plot(paste0(output_dir, "volcano_unstim_supplement_combined"), combined_plot, height = 28, width = 16)

# All patients - stimulated -----------------------------------------------

# for main figure
output_dir = "figures/figure_4/"

cd14_plots = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", x_limit = 1.75, y_limit = 240, space_per_gene = 10)
#cd8_plots = make_volcano_plots_t(cell_type = "CD8_T_cell", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", y_limit = 250, space_per_gene = 9, y_lower_bound = 160)

all_plots = c(cd14_plots)
combined_plot = ggarrange(plotlist = all_plots,
                          nrow = 1, ncol = 2, widths = c(3, 4))

save_plot(paste0(output_dir, "volcano_stim_combined"), combined_plot, height = 6.5, width = 16)

# 
# # supplemental
# output_dir = "figures/supplemental/"
# #cd16_plots = make_volcano_plots(cell_type = "CD16_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", x_limit = 2.5, y_limit = 10)
# #cd4_plots = make_volcano_plots(cell_type = "CD4_T_cell", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim")
# NK_plots = make_volcano_plots(cell_type = "NK", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", stim_status = "stim", y_limit = 50)
# B_plots = make_volcano_plots(cell_type = "B", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", stim_status = "stim", y_limit = 50)
# 
# all_plots = c(NK_plots, B_plots)
# combined_plot = ggarrange(plotlist = all_plots,
#                           nrow = 2, ncol = 2, widths = c(3, 4))
# 
# save_plot(paste0(output_dir, "volcano_stim_supplement_combined"), combined_plot, height = 28, width = 16)



# Stim vs unstim - all patients -------------------------------------------

make_stim_volcano_plots = function(cell_type, output_dir, table_dir, x_limit = 2.5, y_limit = 300) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("stim_v_unstim_", cell_type, "_wilcox.csv"))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data)
    data$`-log10_p_val`[data$`-log10_p_val` == Inf] = 292.2792
    
    # prep data subsets for labeling
    inflammation_data = subset(data, gene_class == "Inflammation")
    antigen_data = subset(data, gene_class == "Antigen Presentation")
    adhesion_data = subset(data, gene_class == "Cell Adhesion")
    activation_data = subset(data, gene_class == "Monocyte Activation")
    
    buffer = y_limit / 10
    space_per_gene = (y_limit - buffer) / (nrow(inflammation_data) + nrow(antigen_data) + nrow(adhesion_data) + nrow(activation_data))
    print(space_per_gene)
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano +
        xlim(-x_limit, x_limit) +
        ylim(0,y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
        gene_labels(data = inflammation_data, xlim = x_limit, ylim = get_lims(y_limit, space_per_gene, nrow(inflammation_data))) +
        gene_labels(data = antigen_data, xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * nrow(inflammation_data), space_per_gene, nrow(antigen_data))) +
        gene_labels(data = adhesion_data,  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data) + nrow(antigen_data)), space_per_gene, nrow(adhesion_data))) +
        gene_labels(data = activation_data,  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data) + nrow(antigen_data) + nrow(adhesion_data)), space_per_gene, nrow(activation_data))) +
        ggtitle(paste("Stim v unstim", gsub("_", " ", cell_type))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))

    filename = paste0(output_dir, "volcano_stim_v_unstim_", cell_type)
    save_plot(filename, plot, height = 12, width = 9)
    write_tsv(data, paste0(table_dir, cell_type, ".tsv"))
  return(plot)
}

output_dir = "figures/figure_4/"
cd14_plots = make_stim_volcano_plots(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/")
#cd8_plots = make_stim_volcano_plots(cell_type = "CD8_T_cell", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", x_limit = 3.5, y_limit = 500)


# TET2 vs DNMT3A ---------------------------------------------------------

make_tet2_dnmt3a_volcano_plots = function(cell_type, output_dir, table_dir, stim_status, x_limit = 2.5, y_limit = 300) {
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0("CHIP_comparison_", cell_type, "_", stim_status, "_wilcox.csv"))
  data = read_csv(paste0(folder, input_filename))
  data = curate_data(data)
  data$`-log10_p_val`[data$`-log10_p_val` == Inf] = 292.2792
  
  # prep data subsets for labeling
  inflammation_data = subset(data, gene_class == "Inflammation")
  antigen_data = subset(data, gene_class == "Antigen Presentation")
  adhesion_data = subset(data, gene_class == "Cell Adhesion")
  activation_data = subset(data, gene_class == "Monocyte Activation")
  
  buffer = y_limit / 10
  space_per_gene = (y_limit - buffer) / (nrow(inflammation_data) + nrow(antigen_data) + nrow(adhesion_data) + nrow(activation_data))
  print(space_per_gene)
  
  # render plot
  (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
      geom_point() +
      format_volcano +
      xlim(-x_limit, x_limit) +
      ylim(0,y_limit) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
      geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
      gene_labels(data = inflammation_data, xlim = x_limit, ylim = get_lims(y_limit, space_per_gene, nrow(inflammation_data))) +
      gene_labels(data = antigen_data, xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * nrow(inflammation_data), space_per_gene, nrow(antigen_data))) +
      gene_labels(data = adhesion_data,  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data) + nrow(antigen_data)), space_per_gene, nrow(adhesion_data))) +
      gene_labels(data = activation_data,  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * (nrow(inflammation_data) + nrow(antigen_data) + nrow(adhesion_data)), space_per_gene, nrow(activation_data))) +
      ggtitle(paste("TET2 v DNMT3A", gsub("_", " ", cell_type), paste0(stim_status, "ulated"))) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  
  filename = paste0(output_dir, "volcano_CHIP_comparison_", cell_type, "_", stim_status)
  save_plot(filename, plot, width = 9)
  write_tsv(data, paste0(table_dir, cell_type, ".tsv"))
  return(plot)
}

output_dir = "figures/figure_4/"
cd14_unstim_plots = make_tet2_dnmt3a_volcano_plots(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "unstim")
cd14_stim_plots = make_tet2_dnmt3a_volcano_plots(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim")

