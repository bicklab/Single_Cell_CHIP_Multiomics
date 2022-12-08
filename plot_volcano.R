
# Volcano plots -----------------------------------------------------------
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


# load gene classes
gene_label_data_mono = read.xlsx("input_data/gene_listp_81722.xlsx", sheet = "CD14 Mono")
gene_label_data_cd8t = read.xlsx("input_data/gene_listp_81722.xlsx", sheet = "CD8 T")
gene_label_data_cd4t = read.xlsx("input_data/gene_listp_81722.xlsx", sheet = "CD4T")

# Helper functions -----------------------------------------------------------

curate_data = function(exp_data, cell_type = "CD14_Mono") {
  if (cell_type == "CD14_Mono") {
    exp_data %>% 
      dplyr::mutate(significant = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
      dplyr::mutate(gene_class = factor(case_when(significant & Gene %in% gene_label_data_mono$Inflammation ~ "Inflammation",
                                                  significant & Gene %in% gene_label_data_mono$Antigen.Presentation ~ "Antigen presentation", 
                                                  significant & Gene %in% gene_label_data_mono$Cell.Adhesion ~ "Cell adhesion",
                                                  significant & Gene %in% gene_label_data_mono$Monocyte.Activation ~ "Monocyte activation",
                                                  TRUE ~ "Other" ),
                                        levels = c("Inflammation", "Antigen presentation", "Cell adhesion", "Monocyte activation", "Other"))) %>% 
      mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
      arrange(desc(gene_class)) %>%
      return()
  } else if (cell_type == "CD8_T_cell") {
    exp_data %>% 
      dplyr::mutate(significant = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
      dplyr::mutate(gene_class = factor(case_when(significant & Gene %in% gene_label_data_cd8t$T.Cell.Activation ~ "T cell activation",
                                                  significant & Gene %in% gene_label_data_cd8t$Intracellular.Signaling ~ "Intracellular signaling", 
                                                  significant & Gene %in% gene_label_data_cd8t$Immunomodulation ~ "Immunomodulation",
                                                  TRUE ~ "Other" ),
                                        levels = c("T cell activation", "Intracellular signaling", "Immunomodulation", "Other"))) %>% 
      mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
      arrange(desc(gene_class)) %>%
      return()
  } else if (cell_type == "CD4_T_cell") {
    exp_data %>% 
      dplyr::mutate(significant = abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) %>%
      dplyr::mutate(gene_class = factor(case_when(significant & Gene %in% gene_label_data_cd4t$`T.cell.activation./.TCR.signaling` ~ "T cell activation/TCR signaling",
                                                  significant & Gene %in% gene_label_data_cd4t$Intracellular.Signaling ~ "Intracellular signaling", 
                                                  significant & Gene %in% gene_label_data_cd4t$Immune.Recruitment ~ "Immune recruitment",
                                                  TRUE ~ "Other" ),
                                        levels = c("T cell activation/TCR signaling", "Intracellular signaling", "Immune recruitment", "Other"))) %>% 
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

gene_labels = function(data, direction, xlim, ylim, few = FALSE) {
  sorted_data = data %>% 
    arrange(desc(`-log10_p_val`))
  if (few) {
    sorted_data = head(sorted_data, 3)
  }
  # we want to highlight IL1B even if it's not top 3, so replace non-interesting gene label with IL1B
  if ("SNHG29" %in% sorted_data$Gene) {
    sorted_data[which(sorted_data$Gene == "SNHG29"),] = data[which(data$Gene == "IL1B"),]
  }
  if (direction == "pos") {
    xlims = c(xlim - 0.5, xlim + 0.5)
  } else {
    xlims = c(-xlim - 0.5, -xlim + 0.5)
  }
  labels = geom_label_repel(data = sorted_data, mapping = aes(label = Gene), size = 3, fill = NA, direction = "y", max.overlaps =  Inf, nudge_x = 0.05, 
                            segment.alpha = 0.25, force_pull = 0, xlim = xlims, ylim = ylim, show.legend = FALSE, box.padding = 0, 
                            segment.linetype = "solid", label.size = 0.08, segment.curvature = ifelse(direction == "neg", -1e-40, 1e-40))
  return(list(labels))
}

make_volcano_plots_mono = function(cell_type, output_dir, table_dir, x_limit = 1.5, y_limit = 315, space_per_gene = 12.1, stim_status = "unstim", few = FALSE, height = 7, width = 9) {
  sample_genotypes = c("TET2", "DNMT3A")
  plots = list()
  
  for (sample_genotype in sample_genotypes) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status, "_wilcox"))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data, cell_type)
    
    # prep data subsets for labeling
    inflammation_data = subset(data, gene_class == "Inflammation")
    antigen_data = subset(data, gene_class == "Antigen presentation")
    adhesion_data = subset(data, gene_class == "Cell adhesion")
    activation_data = subset(data, gene_class == "Monocyte activation")
    
    inflammation_data_pos = inflammation_data %>% filter(avg_log2FC > 0)
    antigen_data_pos = antigen_data %>% filter(avg_log2FC > 0)
    adhesion_data_pos = adhesion_data %>% filter(avg_log2FC > 0)
    activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
    
    inflammation_data_neg = inflammation_data %>% filter(avg_log2FC < 0)
    antigen_data_neg = antigen_data %>% filter(avg_log2FC < 0)
    adhesion_data_neg = adhesion_data %>% filter(avg_log2FC < 0)
    activation_data_neg = activation_data %>% filter(avg_log2FC < 0)
    
    legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
    buffer = 20
    
    check = function(subset) {
      if (few) {
        return(ifelse(nrow(subset) < 3, nrow(subset), 3))
      } else {
        return(nrow(subset))
      }
    }
    
    clean_cell_type = gsub("_", " ", cell_type)
    clean_cell_type = gsub("Mono", "Monocyte", clean_cell_type)
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano +
        xlim(-x_limit, x_limit) +
        ylim(0, y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash", alpha = 0.5) +
        
        gene_labels(data = inflammation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer, space_per_gene, check(inflammation_data_pos)), few) +
        gene_labels(data = antigen_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos)), space_per_gene, check(antigen_data_pos)), few) +
        gene_labels(data = adhesion_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos) + check(antigen_data_pos)), space_per_gene, check(adhesion_data_pos)), few) +
        gene_labels(data = activation_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos) + check(adhesion_data_pos) + check(antigen_data_pos)), space_per_gene, check(activation_data_pos)), few) +
        
        gene_labels(data = inflammation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit - buffer, space_per_gene, check(inflammation_data_neg)), few) +
        gene_labels(data = antigen_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg)), space_per_gene, check(antigen_data_neg)), few) +
        gene_labels(data = adhesion_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg) + check(antigen_data_neg)), space_per_gene, check(adhesion_data_neg)), few) +
        gene_labels(data = activation_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg) + check(adhesion_data_neg) + check(antigen_data_neg)), space_per_gene, check(activation_data_neg)), few) +
        
        ggtitle(paste(sample_genotype)) +
        theme(legend.position = ifelse(legend, "right", "none")))
    
    write_tsv(data, paste0(table_dir, sample_genotype, "_", cell_type, "_", stim_status, ".tsv"))
    plots[[length(plots)+1]] = plot(plot)  
  }
  return(plots)
}

make_volcano_plots_cd8t = function(cell_type, output_dir, table_dir, x_limit = 1.5, y_limit = 315, space_per_gene = 12.1, stim_status = "unstim", y_lower_bound = 100, few = FALSE, height = 7, width = 9) {
  sample_genotypes = c("TET2", "DNMT3A")
  plots = list()
  
  for (sample_genotype in sample_genotypes) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status, "_wilcox"))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data, cell_type)
    
    # prep data subsets for labeling
    activation_data = subset(data, gene_class == "T cell activation")
    intracellular_data = subset(data, gene_class == "Intracellular signaling")
    immunomodulation_data = subset(data, gene_class == "Immunomodulation")
    
    activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
    intracellular_data_pos = intracellular_data %>% filter(avg_log2FC > 0)
    immunomodulation_data_pos = immunomodulation_data %>% filter(avg_log2FC > 0)
    
    activation_data_neg = activation_data %>% filter(avg_log2FC < 0)
    intracellular_data_neg = intracellular_data %>% filter(avg_log2FC < 0)
    immunomodulation_data_neg = immunomodulation_data %>% filter(avg_log2FC < 0)
    
    legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
    
    check = function(subset) {
      if (few) {
        return(ifelse(nrow(subset) < 3, nrow(subset), 3))
      } else {
        return(nrow(subset))
      }
    }
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano_2 +
        xlim(-x_limit, x_limit) +
        ylim(0, y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash", alpha = 0.5) +
        
        gene_labels(data = activation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, check(activation_data_pos)), few) +
        gene_labels(data = intracellular_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * check(activation_data_pos), space_per_gene, check(intracellular_data_pos)), few) +
        gene_labels(data = immunomodulation_data_pos, direction = "pos",  xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (check(activation_data_pos) + check(intracellular_data_pos)), space_per_gene, check(immunomodulation_data_pos)), few) +
        
        gene_labels(data = activation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, check(activation_data_neg)), few) +
        gene_labels(data = intracellular_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * check(activation_data_neg), space_per_gene, check(intracellular_data_neg)), few) +
        gene_labels(data = immunomodulation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (check(activation_data_neg) + check(intracellular_data_neg)), space_per_gene, check(immunomodulation_data_neg)), few) +
        theme(legend.position = ifelse(legend, "right", "none")) +
        ggtitle(paste(sample_genotype, gsub("_", " ", cell_type))))
    
    write_tsv(data, paste0(table_dir, sample_genotype, "_", cell_type, "_", stim_status, ".tsv"))
    plots[[length(plots)+1]] = plot(plot)  
  }
  return(plots)
}

make_volcano_plots_cd4t = function(cell_type, output_dir, table_dir, x_limit = 1.5, y_limit = 315, space_per_gene = 12.1, stim_status = "unstim", y_lower_bound = 100, few = FALSE, height = 7, width = 9) {
  sample_genotypes = c("TET2", "DNMT3A")
  plots = list()
  
  for (sample_genotype in sample_genotypes) {
    folder = "tables/supplemental/"
    input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status, "_wilcox"))
    data = read_csv(paste0(folder, input_filename))
    data = curate_data(data, cell_type) %>%
      filter(!(Gene %in% c("HBB", "EPB41", "MRPS6")))
    
    # prep data subsets for labeling
    activation_data = subset(data, gene_class == "T cell activation/TCR signaling")
    intracellular_data = subset(data, gene_class == "Intracellular signaling")
    immune_data = subset(data, gene_class == "Immune recruitment")
    
    activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
    intracellular_data_pos = intracellular_data %>% filter(avg_log2FC > 0)
    immune_data_pos = immune_data %>% filter(avg_log2FC > 0)
    
    activation_data_neg = activation_data %>% filter(avg_log2FC < 0)
    intracellular_data_neg = intracellular_data %>% filter(avg_log2FC < 0)
    immune_data_neg = immune_data %>% filter(avg_log2FC < 0)
    
    legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
    
    check = function(subset) {
      if (few) {
        return(ifelse(nrow(subset) < 3, nrow(subset), 3))
      } else {
        return(nrow(subset))
      }
    }
    
    # render plot
    (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
        geom_point() +
        format_volcano_3 +
        xlim(-x_limit, x_limit) +
        ylim(0, y_limit) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash", alpha = 0.5) +
        
        gene_labels(data = activation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, check(activation_data_pos)), few) +
        gene_labels(data = intracellular_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * check(activation_data_pos), space_per_gene, check(intracellular_data_pos)), few) +
        gene_labels(data = immune_data_pos, direction = "pos",  xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (check(activation_data_pos) + check(intracellular_data_pos)), space_per_gene, check(immune_data_pos)), few) +
        
        gene_labels(data = activation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound, space_per_gene, check(activation_data_neg)), few) +
        gene_labels(data = intracellular_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * check(activation_data_neg), space_per_gene, check(intracellular_data_neg)), few) +
        gene_labels(data = immune_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims_t(y_lower_bound + space_per_gene * (check(activation_data_neg) + check(intracellular_data_neg)), space_per_gene, check(immune_data_neg)), few) +
        theme(legend.position = ifelse(legend, "right", "none")) +
        ggtitle(paste(sample_genotype, gsub("_", " ", cell_type))))
    
    write_tsv(data, paste0(table_dir, sample_genotype, "_", cell_type, "_", stim_status, ".tsv"))
    plots[[length(plots)+1]] = plot(plot)  
  }
  return(plots)
}

remove_legends = function(plots) {
  for (i in c(1:length(plots))) {
    plots[[i]] = plots[[i]] + guides(color = "none")
  }
  return(plots)
}

# All patients - unstimulated ---------------------------------------------

output_dir = "figures/supplemental/"

cd14_plots = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/", x_limit = 2)
cd8_plots = make_volcano_plots_cd8t(cell_type = "CD8_T_cell", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/", x_limit = 2.5, space_per_gene = 14)
cd4_plots = make_volcano_plots_cd4t(cell_type = "CD4_T_cell", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", x_limit = 2.5, space_per_gene = 14)

all_plots = c(cd14_plots, cd8_plots)
combined_plot = ggarrange(plotlist = all_plots, nrow = 2, ncol = 2, widths = c(3, 4))
save_plot(paste0(output_dir, "volcano_combined"), combined_plot, height = 13, width = 14) #13 and 14

t_cell_plots = ggarrange(plotlist = c(cd4_plots, cd8_plots), nrow = 2, ncol = 2, widths = c(3, 4))
save_plot(paste0(output_dir, "t_cell_volcano_combined"), t_cell_plots, height = 11.5, width = 14) #png looks nice at height = 10.9 and width = 14

output_dir = "figures/figure_2/"

cd14_plots_few = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/", x_limit = 2, space_per_gene = 20, few = TRUE, height = 4, width = 7)
cd8_plots_few = make_volcano_plots_cd8t(cell_type = "CD8_T_cell", output_dir = "figures/figure_2/", table_dir = "tables/figure_2/", x_limit = 2, space_per_gene = 20, few = TRUE, height = 4, width = 7)
all_plots_few = c(cd14_plots_few, cd8_plots_few) %>% remove_legends()

mono_plots = ggarrange(plotlist = c(cd14_plots_few), nrow = 1, ncol = 2, widths = c(3, 4))
save_plot(paste0(output_dir, "mono_volcano_combined_few"), mono_plots, height = 4.1, width = 11.75, png = TRUE)

combined_plot_few = ggarrange(plotlist = all_plots_few[c(1,3,2,4)], nrow = 2, ncol = 2)
save_plot(paste0(output_dir, "volcano_combined_few"), combined_plot_few, height = 8, width = 12) #8.5 and 13

legend_mono = cowplot::get_legend(cd14_plots_few[[2]]) 
legend_gg_plot = as_ggplot(legend_mono)
save_plot(paste0(output_dir, "volcano_legend_mono.pdf"), legend_gg_plot, height = 2, width = 2, transparent = TRUE)

legend_t = cowplot::get_legend(cd8_plots_few[[2]]) 
as_ggplot(legend_t)
ggsave(filename = paste0(output_dir, "volcano_legend_t.pdf"), height = 2, width = 2)


# All patients - stimulated -----------------------------------------------

output_dir = "figures/supplemental/"
cd14_plots = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", x_limit = 2, y_limit = 240, space_per_gene = 10)
all_plots = c(cd14_plots)
combined_plot = ggarrange(plotlist = all_plots, nrow = 1, ncol = 2, widths = c(3, 4))
save_plot(paste0(output_dir, "volcano_stim_combined"), combined_plot, height = 8.5, width = 14)

output_dir = "figures/figure_4/"
cd14_plots_few = make_volcano_plots_mono(cell_type = "CD14_Mono", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", x_limit = 2, y_limit = 240, space_per_gene = 15, few = TRUE)
all_plots_few = c(cd14_plots_few) %>% remove_legends()
combined_plot_few = ggarrange(plotlist = all_plots_few, nrow = 1, ncol = 2)
save_plot(paste0(output_dir, "volcano_stim_combined_few"), combined_plot_few, height = 4, width = 12)
as_ggplot(legend_mono)
ggsave(filename = paste0(output_dir, "volcano_legend_mono.pdf"), height = 2, width = 2)


# TET2 vs DNMT3A ---------------------------------------------------------

make_tet2_dnmt3a_volcano_plots = function(cell_type, output_dir, table_dir, stim_status, x_limit = 2.5, y_limit = 300, space_per_gene = 12.1, few = FALSE) {
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0("CHIP_comparison_", cell_type, "_", stim_status, "_wilcox.csv"))
  data = read_csv(paste0(folder, input_filename))
  data = curate_data(data)
  data$`-log10_p_val`[data$`-log10_p_val` == Inf] = 292.2792
  
  # prep data subsets for labeling
  inflammation_data = subset(data, gene_class == "Inflammation")
  antigen_data = subset(data, gene_class == "Antigen presentation")
  adhesion_data = subset(data, gene_class == "Cell adhesion")
  activation_data = subset(data, gene_class == "Monocyte activation")
  
  inflammation_data_pos = inflammation_data %>% filter(avg_log2FC > 0)
  antigen_data_pos = antigen_data %>% filter(avg_log2FC > 0)
  adhesion_data_pos = adhesion_data %>% filter(avg_log2FC > 0)
  activation_data_pos = activation_data %>% filter(avg_log2FC > 0)
  
  inflammation_data_neg = inflammation_data %>% filter(avg_log2FC < 0)
  antigen_data_neg = antigen_data %>% filter(avg_log2FC < 0)
  adhesion_data_neg = adhesion_data %>% filter(avg_log2FC < 0)
  activation_data_neg = activation_data %>% filter(avg_log2FC < 0)
  
  check = function(subset) {
    if (few) {
      return(ifelse(nrow(subset) < 3, nrow(subset), 3))
    } else {
      return(nrow(subset))
    }
  }
  
  buffer = y_limit / 10
  
  # render plot
  (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
      geom_point() +
      format_volcano +
      xlim(-x_limit, x_limit) +
      ylim(0,y_limit) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
      geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
      
      gene_labels(data = inflammation_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer, space_per_gene, check(inflammation_data_pos)), few) +
      gene_labels(data = antigen_data_pos, direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos)), space_per_gene, check(antigen_data_pos)), few) +
      gene_labels(data = adhesion_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos) + check(antigen_data_pos)), space_per_gene, check(adhesion_data_pos)), few) +
      gene_labels(data = activation_data_pos,  direction = "pos", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_pos) + check(adhesion_data_pos) + check(antigen_data_pos)), space_per_gene, check(activation_data_pos)), few) +
      
      gene_labels(data = inflammation_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit - buffer, space_per_gene, check(inflammation_data_neg)), few) +
      gene_labels(data = antigen_data_neg, direction = "neg", xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg)), space_per_gene, check(antigen_data_neg)), few) +
      gene_labels(data = adhesion_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg) + check(antigen_data_neg)), space_per_gene, check(adhesion_data_neg)), few) +
      gene_labels(data = activation_data_neg, direction = "neg",  xlim = x_limit, ylim = get_lims(y_limit - buffer - space_per_gene * (check(inflammation_data_neg) + check(adhesion_data_neg) + check(antigen_data_neg)), space_per_gene, check(activation_data_neg)), few) +
      
      ggtitle(paste("TET2 v DNMT3A", gsub("_", " ", gsub("Mono", "Monocyte", cell_type)), paste0(stim_status, "ulated"))) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  
  filename = paste0(output_dir, "volcano_CHIP_comparison_", cell_type, "_", stim_status)
  save_plot(filename, plot, height = 5, width = 7)
  write_tsv(data, paste0(table_dir, cell_type, "_TET2_vs_DNMT3A.tsv"))
  return(plot)
}

output_dir = "figures/supplemental/"
cd14_unstim_plots = make_tet2_dnmt3a_volcano_plots(cell_type = "CD14_Mono", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", stim_status = "unstim", space_per_gene = 15)
cd14_stim_plots = make_tet2_dnmt3a_volcano_plots(cell_type = "CD14_Mono", output_dir = "figures/supplemental/", table_dir = "tables/supplemental/", stim_status = "stim", space_per_gene = 15)


