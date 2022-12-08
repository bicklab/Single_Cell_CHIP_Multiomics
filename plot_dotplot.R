
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/figure_4/"
gene_label_data = read.xlsx("input_data/gene_listp_61322.xlsx")[,1:4]
individ_label_data = read.xlsx("input_data/individ_gene_list_fordots.xlsx")
total_label_data = rbind(gene_label_data, individ_label_data)

# pulling top genes for each cell type based on unstim differential expression between patients and controls (all patient data)
collect_interesting_genes = function(cell_type, patient_genes = TRUE, stim_status, all_genes = FALSE) {
  folder = "tables/supplemental/"
  
  tet2_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "TET2_", gsub(" ", "_", cell_type), ".*_", stim_status))
  dnmt3a_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "DNMT3A_", gsub(" ", "_", cell_type), ".*_", stim_status))
  tet2_significance_data = read_csv(paste0(folder, tet2_input_filename))
  dnmt3a_significance_data = read_csv(paste0(folder, dnmt3a_input_filename))
  
  
  index = c("Inflammation", "Antigen.Presentation", "Cell.Adhesion", "Monocyte.Activation")
  values = c("Inflammation", "Antigen presentation", "Cell adhesion", "Monocyte activation")
  
  pathway_genes = total_label_data %>%
    pivot_longer(cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
    na.omit() %>%
    full_join(tet2_significance_data, by = "Gene") %>%
    full_join(dnmt3a_significance_data, by = "Gene") %>%
    filter(!is.na(Pathway))
  
  if (!all_genes) {
    pathway_genes = pathway_genes %>% 
      filter(p_val_adj.x < 0.05 | p_val_adj.y < 0.05)
  }
  
  pathway_genes = pathway_genes %>%
    dplyr::select(Gene, Pathway)
  
  pathway_genes$Pathway = values[match(pathway_genes$Pathway, index)]
  pathway_genes = pathway_genes %>%
    mutate(Pathway = factor(Pathway, levels = c("Inflammation", "Antigen presentation", "Cell adhesion", "Monocyte activation"))) %>%
    na.omit() %>%
    unique()
  return(pathway_genes)
}

dotplot_stats = function(data, sample_genotype, cell_type, stim_status = NULL, split = "Clone", patient_genes = TRUE, all_genes = FALSE) {
  if (!is.null(stim_status)) {
    Idents(data) = "STIM"
    data = subset(data, idents = stim_status)
  }

  DefaultAssay(data) = "RNA"
  Idents(data) = "celltype.de"
  dotplot_data = subset(data, idents = cell_type)
  
  interesting_genes = collect_interesting_genes(cell_type, stim_status = ifelse(stim_status == "STIM", "stim", "unstim"), all_genes = all_genes)
  features = interesting_genes$Gene %>% unique()
  dotplot = DotPlot(dotplot_data, features = features, split.by = split, cluster.idents = TRUE)
  
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), sample_genotype, "_", gsub(" ", "_", cell_type), ".*_", ifelse(stim_status == "STIM", "stim", "unstim")))
  stats_data = read_csv(paste0(folder, input_filename))
  
  pathway_genes = full_join(interesting_genes, dotplot$data, by = c("Gene" = "features.plot")) %>%
    na.omit() %>% 
    arrange(id) %>% 
    mutate(log_avg_exp = log(avg.exp)) %>%
    mutate(id = factor(id, levels = c("CD14 Mono_Mutant", "CD14 Mono_Wildtype", 
                                      "CD4 T cell_Mutant", "CD4 T cell_Wildtype", 
                                      "CD8 T cell_Mutant", "CD8 T cell_Wildtype"))) %>%
    left_join(stats_data) %>%
    mutate(significant = factor(p_val_adj < 0.05, levels = c(FALSE, TRUE)))
  
  write_tsv(pathway_genes, paste0("tables/figure_4/pathway_gene_stats_", stim_status, "_", sample_genotype, "_", cell_type, ".tsv"))
  return(pathway_genes)
}

remove_legends = function(plot) {
  return(plot + guides(color = "none") + guides(fill = "none") + guides(size = "none") + guides(alpha = "none") + guides(shape = "none"))
}


make_dotplot = function(pathway_genes_tet2, pathway_genes_dnmt3a, stim_status, all_genes = FALSE) {
  
  pathway_genes = rbind(pathway_genes_tet2 %>% mutate(id = paste0("TET2_", id)), 
                        pathway_genes_dnmt3a %>% mutate(id = paste0("DNMT3A_", id)))
  
  pathway_genes = pathway_genes %>%
    mutate(`Adj p-value` = case_when(p_val_adj < 0.0001 ~ "< 0.0001",
                                        p_val_adj < 0.05 ~ "0.0001 - 0.05",
                                        TRUE ~ "Not significant")) %>%
    mutate(`Avg exp` = case_when(avg.exp < 50 ~ "< 50",
                                 avg.exp > 50 & avg.exp < 100 ~ "> 50 & < 100",
                                 avg.exp >= 100 ~ ">= 100")) %>%
    mutate(`Pct exp` = case_when(pct.exp < 10 ~ "< 10%",
                                 pct.exp >= 10 ~ ">= 10%")) %>%
    mutate(`Significant` = case_when(p_val_adj < 0.05 ~ "yes",
                                 pct.exp >= 0.05 ~ "no")) %>%
    mutate(`Avg log2 fold change` = case_when(avg_log2FC < 0 ~ "< 0",
                                              avg_log2FC > 0 ~ "> 0")) %>%
    arrange(p_val_adj) %>%
    mutate(Gene = factor(Gene, levels = unique(Gene))) %>%
    filter(id == "TET2_CD14 Mono_Wildtype" | id == "DNMT3A_CD14 Mono_Wildtype")
  
  pathway_genes$id[pathway_genes$id == "TET2_CD14 Mono_Wildtype"] = "TET2"
  pathway_genes$id[pathway_genes$id == "DNMT3A_CD14 Mono_Wildtype"] = "DNMT3A"

  pathway_genes_matched = pathway_genes %>%
    select(Gene, Pathway, id, avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp)
  
  shapes = c(21, 1)
  names(shapes) = c("yes", "no")
  
  alphas = c(1, 0.65)
  names(alphas) = c("> 0", "< 0")

  borders = c("white", "black")
  names(borders) = c("< 10%", ">= 10%")
  
  sizes = c(3, 5)
  names(alphas) = c("> 0", "< 0")
  
  (dotplot_compiled = ggplot(pathway_genes, aes(x = Gene, y = id)) +
      geom_point(aes(size = `Avg log2 fold change`, fill = Pathway, shape = Significant), color = "black") + # , stroke = 1.5 (border), shape = 21
      scale_alpha_manual(values = alphas) +
      scale_size_manual(values = sizes) +
      scale_fill_manual(values = dot_colors) +
      scale_color_manual(values = dot_colors) +
      scale_shape_manual(values = shapes) +
      facet_grid(~Pathway, scales = "free", space = "free") +
      theme_bw() +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text = element_blank(),
            strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            legend.key.size = unit(0, "lines")) +
      geom_hline(yintercept = seq(1.5, length(unique(pathway_genes$id))-0.5, 1), lwd = 0.5, colour = "grey"))
  print(dotplot_compiled)
  
  filename = paste0(output_dir, "celltype_dotplot_", stim_status, "_", gsub(" ", "_", cell_type), ifelse(all_genes, "_all_genes", ""))
  save_plot(filename, dotplot_compiled %>% remove_legends(), width = ifelse(all_genes, 8, 6), height = 1.5)
  
  legend = cowplot::get_legend(dotplot_compiled) 
  ggsave(paste0(filename, "_legend.pdf"), as_ggplot(legend), height = 3, width = 2)
  
  print(filename)
  write_tsv(pathway_genes, paste0(gsub("figures", "tables", filename), ".tsv"))
}

cell_types = "CD14 Mono"

for (cell_type in cell_types) {
  output_dir = "figures/figure_3/"
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "VEH"), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "VEH"), 
               stim_status = "VEH")
  output_dir = "figures/figure_4/"
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "STIM"), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "STIM"), 
               stim_status = "STIM")
  
  output_dir = "figures/supplemental/"
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "VEH", all_genes = TRUE), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "VEH", all_genes = TRUE), 
               stim_status = "VEH", all_genes = TRUE)
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "STIM", all_genes = TRUE), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "STIM", all_genes = TRUE), 
               stim_status = "STIM", all_genes = TRUE)
}
  
  
