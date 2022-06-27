
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/figure_4/"
gene_df = read_csv("input_data/Gene_class.csv")
gene_label_data = read.xlsx("input_data/gene_listp_61322.xlsx")

# pulling top genes for each cell type based on unstim differential expression between patients and controls (all patient data)
collect_interesting_genes = function(cell_type, patient_genes = TRUE, stim_status) {
  folder = "tables/supplemental/"
  
  tet2_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "TET2_", gsub(" ", "_", cell_type), ".*_", stim_status))
  dnmt3a_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "DNMT3A_", gsub(" ", "_", cell_type), ".*_", stim_status))
  tet2_significance_data = read_csv(paste0(folder, tet2_input_filename))
  dnmt3a_significance_data = read_csv(paste0(folder, dnmt3a_input_filename))
  
  
  index = c("Inflammation", "Antigen.Presentation", "Cell.Adhesion", "Monocyte.Activation")
  values = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation")
  
  pathway_genes = gene_label_data %>%
    pivot_longer(cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
    na.omit() %>%
    full_join(tet2_significance_data, by = "Gene") %>%
    full_join(dnmt3a_significance_data, by = "Gene") %>%
    filter(!is.na(Pathway)) %>%
    dplyr::select(Gene, Pathway)
  
  pathway_genes$Pathway = values[match(pathway_genes$Pathway, index)]
  pathway_genes = pathway_genes %>%
    mutate(Pathway = factor(Pathway, levels = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation"))) %>%
    na.omit() %>%
    unique()
  return(pathway_genes)
}

dotplot_stats = function(data, sample_genotype, cell_type, stim_status = NULL, split = "Clone", patient_genes = TRUE) {
  if (!is.null(stim_status)) {
    Idents(data) = "STIM"
    data = subset(data, idents = stim_status)
  }

  DefaultAssay(data) = "RNA"
  Idents(data) = "celltype.de"
  dotplot_data = subset(data, idents = cell_type)
  
  interesting_genes = collect_interesting_genes(cell_type, stim_status = ifelse(stim_status == "STIM", "stim", "unstim"))
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

make_dotplot = function(pathway_genes_tet2, pathway_genes_dnmt3a, stim_status) {
  
  pathway_genes = rbind(pathway_genes_tet2 %>% mutate(id = paste0("TET2_", id)), 
                        pathway_genes_dnmt3a %>% mutate(id = paste0("DNMT3A_", id)))
  
  pathway_genes = pathway_genes %>%
    mutate(`p-value category` = case_when(p_val_adj < 0.001 ~ 5,
                                        p_val_adj < 0.01 ~ 4,
                                        p_val_adj < 0.05 ~ 3,
                                        TRUE ~ 2))

  shapes = c(21, 22)
  names(shapes) = c(TRUE, FALSE)
  (dotplot_compiled = ggplot(pathway_genes, aes(x = Gene, y = id)) +
      geom_point(aes(size = p_val_adj, alpha = avg.exp.scaled, fill = Pathway), shape = 21) +
      scale_size_continuous(trans = "reverse", breaks = c(0.05), range = c(3, 7)) +
      scale_alpha(range = c(0.25, 1)) +
      scale_fill_manual(values = colors) +
      scale_shape_manual(values = shapes) +
      facet_grid(~Pathway, scales = "free", space = "free") +
      theme_bw() +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text = element_text(size = rel(1.2), face = "bold"),
            strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_hline(yintercept=seq(1.5, length(unique(pathway_genes$id))-0.5, 1), lwd=0.5, colour="grey") +
      ggtitle(paste(cell_type, stim_status)) +
      guides(fill = "none"))
  print(dotplot_compiled)
  
  filename = paste0(output_dir, "celltype_dotplot_", stim_status, "_", gsub(" ", "_", cell_type))
  save_plot(filename, dotplot_compiled, width = 16, height = 4)
  
  print(filename)
  write_tsv(pathway_genes, paste0(gsub("figures", "tables", filename), ".tsv"))
}

# sample_genotypes = c("DNMT3A", "TET2")
# cell_types = c("CD14 Mono", "CD16 Mono", "CD8 T cell", "CD4 T cell", "NK", "B")
cell_type = "CD14 Mono"

for (cell_type in cell_types) {
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "VEH"), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "VEH"), 
               stim_status = "VEH")
  make_dotplot(pathway_genes_tet2 = dotplot_stats(tet2_patient_singlet_data, "TET2", cell_type, stim_status = "STIM"), 
               pathway_genes_dnmt3a = dotplot_stats(dnmt3a_patient_singlet_data, "DNMT3A", cell_type, stim_status = "STIM"), 
               stim_status = "STIM")
}
  
  
