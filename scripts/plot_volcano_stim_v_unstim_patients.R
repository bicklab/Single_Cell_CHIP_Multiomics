
# Volcano plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


# load gene classes
gene_df = read_csv("input_data/Gene_class.csv")
x_chr_genes = read_tsv("input_data/x_genes.csv")

# Helper functions -----------------------------------------------------------

curate_data = function(exp_data, il6 = FALSE) {
  if (il6) {
    exp_data %>% 
      dplyr::mutate(gene_class = factor(case_when(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$STAT3 ~ "IL6 Response",
                                                  abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$Inflammatory_Response ~ "Inflammatory Response",
                                                  abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$Cell_adhesion_molecules ~ "Cell Adhesion Molecules", 
                                                  abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$APC ~ "Antigen Processing and Presentation",
                                                  TRUE ~ "Other" ),
                                        levels = c("Antigen Processing and Presentation", "Cell Adhesion Molecules", "IL6 Response", "Inflammatory Response", "Other"))) %>% 
      mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
      arrange(desc(gene_class)) %>%
      return()
  } else {
    exp_data %>% 
      dplyr::mutate(gene_class = factor(case_when(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$Inflammatory_Response ~ "Inflammatory Response",
                                                  abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$Cell_adhesion_molecules ~ "Cell Adhesion Molecules", 
                                                  abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 & Gene %in% gene_df$APC ~ "Antigen Processing and Presentation",
                                                  TRUE ~ "Other" ),
                                        levels = c("Antigen Processing and Presentation", "Cell Adhesion Molecules", "Inflammatory Response", "Other"))) %>% 
      mutate(`-log10_p_val` = -log10(p_val_adj)) %>%
      arrange(desc(gene_class)) %>%
      return()
  }
}

gene_labels = function(data, xlim, ylim) {
  filtered_data = data %>% 
    arrange(desc(`-log10_p_val`))
  positive = filtered_data %>%
    filter(avg_log2FC > 0) %>%
    head(4)
  negative = filtered_data %>%
    filter(avg_log2FC < 0) %>%
    head(4)
  return(list(geom_label_repel(data = positive, aes(label = Gene), direction = "y", max.overlaps =  Inf,
                               xlim = c(xlim - 0.5, xlim), ylim = ylim, show.legend = FALSE, segment.linetype = "dashed", box.padding = 0.1),
              geom_label_repel(data = negative, aes(label = Gene), direction = "y", max.overlaps =  Inf,
                               xlim = c(-xlim, -xlim + 0.5), ylim = ylim, show.legend = FALSE, segment.linetype = "dashed", box.padding = 0.1)))
}

get_lims = function(top, space_per_gene) {
  return(c(top - space_per_gene * 5 - 10, top - 10))
}

prep_data = function(seurat_data) {
  DefaultAssay(seurat_data) = "RNA"
  normalized_data = NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
  all_genes = rownames(normalized_data)
  scaled_data = ScaleData(normalized_data, features = all_genes)
  return(scaled_data)
}

perform_stats = function(data, cell_type, mutation, stim_status) {
  Idents(data) = "celltype.de"
  cell_type_data = subset(data, idents = cell_type)
  
  Idents(cell_type_data) = "Clone"
  sample_names = rownames(cell_type_data$Clone)
  
  num_cells = table(cell_type_data@meta.data$Clone, cell_type_data@meta.data$celltype.de) %>% min()
  print(num_cells)
  if (num_cells < 3 | (cell_type_data@meta.data$Clone %>% unique() %>% length()) < 2) {
    print("Cell counts below 3.")
    return()
  }
  # cell_list = WhichCells(cell_type_data, idents = sample_names, downsample = num_cells, seed = 100)
  # data_subset = cell_type_data[,cell_list]
  
  method = "wilcox"
  response = FindMarkers(cell_type_data, 
                         ident.1 = "Mutant",
                         ident.2 = "Wildtype",
                         test.use = method,
                         assay = "RNA")
  
  clean_response = response %>%
    arrange(p_val_adj) %>%
    na.omit() %>%
    rownames_to_column("Gene") %>%
    filter(!(Gene %in% x_chr_genes$Gene))
  
  filename = paste0(output_dir, "patient_", mutation, "_", gsub(" ", "_", cell_type), "_", method, "_", stim_status, ".csv")
  write_csv(clean_response, file = filename)
}

make_stim_v_unstim_volcano_plots = function(cell_type, sample_genotype, output_dir, table_dir, x_limit = 1.5, y_limit = 300) {
  cell_type = gsub(" ", "_", cell_type)
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0("patient_", sample_genotype, "_", cell_type, ".*"))
  data = read_csv(paste0(folder, input_filename))
  data = curate_data(data)
  
  # prep data subsets for labeling
  inflammatory_data = subset(data, gene_class == "Inflammatory Response")
  antigen_data = subset(data, gene_class == "Antigen Processing and Presentation")
  adhesion_data = subset(data, gene_class == "Cell Adhesion Molecules")
  
  space_per_gene = (y_limit - y_limit / 10) / 20
  legend = ifelse(sample_genotype == "TET2", FALSE, TRUE)
  
  # render plot
  (plot = ggplot(data, aes(x = avg_log2FC, y = `-log10_p_val`, color = gene_class)) +
      geom_point() +
      format_volcano +
      xlim(-x_limit, x_limit) +
      ylim(0,y_limit) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "longdash") +
      geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
      gene_labels(data = inflammatory_data, xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * 10, space_per_gene)) +
      gene_labels(data = antigen_data,  xlim = x_limit, ylim = get_lims(y_limit - space_per_gene * 5, space_per_gene)) +
      gene_labels(data = adhesion_data, xlim = x_limit, ylim = get_lims(y_limit, space_per_gene)) +
      ggtitle(paste(gsub("_", " ", cell_type), sample_genotype, "patient", paste0(stim_status, "ulated"))) +
      theme(legend.position = ifelse(legend, "right", "none"), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()))
  
  filename = paste0(output_dir, "patient_", sample_genotype, "_", cell_type, )
  save_plot(filename, plot, legend, width = 9)
  write_tsv(data, paste0(table_dir, "patient_", sample_genotype, "_", cell_type, ".tsv"))
  
  return(plot)
}

# Run for each cell type -----------------------------------------------------------
# FUNCTION IS NOT YET FUNCTIONAL HAHA COMING BACK TO THIS LATER
for (cell_type in cell_types) {
  make_volcano_plots(cell_type = cell_type, sample_genotype = "TET2", output_dir = "figures/figure_3/", table_dir = "tables/figure_3/", stim_status = "unstim", y_limit = 20)
  make_volcano_plots(cell_type = cell_type, sample_genotype = "TET2", output_dir = "figures/figure_4/", table_dir = "tables/figure_4/", stim_status = "stim", y_limit = 20)
}


