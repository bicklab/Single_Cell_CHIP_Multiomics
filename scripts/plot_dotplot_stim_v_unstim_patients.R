
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/figure_4/"
gene_df = read_csv("input_data/Gene_class.csv")

# pulling top genes for each cell type based on differential expression between stim and unstim (all patient data)
collect_interesting_genes = function(sample_genotype, cell_type, patient_genes = FALSE) {
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0("stim_v_unstim_", gsub(" ", "_", cell_type), "_wilcox.csv"))
  significance_data = read_csv(paste0(folder, input_filename))
  
  index = c("APC", "Cell_adhesion_molecules", "Inflammatory_Response", "STAT3")
  values = c("Antigen", "Adhesion", "Inflammatory", "IL6")
  
  # get top 8 genes based on p_val for each pathway
  pathway_genes = gene_df %>%
    pivot_longer(cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
    na.omit() %>%
    inner_join(significance_data) %>%
    group_by(Pathway) %>%
    arrange(p_val_adj) %>%
    head(30) %>%
    dplyr::select(Gene, Pathway)
  
  pathway_genes$Pathway = values[match(pathway_genes$Pathway, index)]
  pathway_genes = pathway_genes %>%
    mutate(Pathway = factor(Pathway, levels = c("Antigen", "Adhesion", "IL6"))) %>%
    na.omit() %>%
    unique()
  return(pathway_genes)
}

data = singlet_data

make_dotplot = function(data, sample_genotype, cell_type, patient_genes = FALSE) {
  
  interesting_genes = collect_interesting_genes(sample_genotype, cell_type, patient_genes)
  features = interesting_genes$Gene %>% unique()
  
  DefaultAssay(data) = "RNA"
  Idents(data) = "celltype.de"
  cell_type_data = subset(data, idents = cell_type)
  
  all_pathway_genes = c("id", "Gene", "avg.exp.scaled", "pct.exp", "Pathway")
  names(all_pathway_genes) = c("id", "Gene", "avg.exp.scaled", "pct.exp", "Pathway")
  
  patients = data@meta.data$orig.ident %>% unique()
  
  for (patient in patients[1:11]) {
    # Filter to patient
    Idents(cell_type_data) = "orig.ident"
    dotplot_data = subset(cell_type_data, idents = patient)
    
    split = "STIM"  
    dotplot = DotPlot(dotplot_data, features = features, split.by = split, cluster.idents = TRUE)
    
    pathway_genes = full_join(interesting_genes, dotplot$data, by = c("Gene" = "features.plot")) %>%
      na.omit() %>% 
      arrange(id) %>% 
      mutate(log_avg_exp = log(avg.exp)) %>%
      dplyr::select(id, Gene, avg.exp.scaled, pct.exp, Pathway)
    
    if (is.null(nrow(all_pathway_genes))) {
      all_pathway_genes = pathway_genes
    } else {
      all_pathway_genes = rbind(all_pathway_genes, pathway_genes)
    }
  }
  
  (dotplot_compiled = ggplot(all_pathway_genes, aes(x = Gene, y = id)) +
      scale_color_manual(values = wes_palette("Royal1", 20, "continuous")[c(1,7,18,20)]) +
      geom_point(aes(size = pct.exp, alpha = avg.exp.scaled, color = Pathway)) +
      facet_grid(~Pathway, scales = "free") +
      theme_bw() +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text = element_text(size = rel(1.2), face = "bold"),
            strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_hline(yintercept=seq(1.5, length(unique(all_pathway_genes$id))-0.5, 1), 
                 lwd=0.5, colour="grey") +
      ggtitle(paste(sample_genotype, cell_type, stim_status)))
  
  filename = paste0(output_dir, "all_stim_patients_dotplot_", sample_genotype, "_", gsub(" ", "_", cell_type))
  save_plot(filename, dotplot_compiled, width = 10)
}


make_dotplot(singlet_data, "DNMT3A", "CD14 Mono")
make_dotplot(singlet_data, "TET2", "CD14 Mono")


