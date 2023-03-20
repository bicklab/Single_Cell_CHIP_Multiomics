
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/supplemental/"
gene_df = read_csv("input_data/Gene_class.csv")

# pulling top genes for each cell type based on differential expression between stim and unstim (all patient data)
collect_interesting_genes = function(cell_type, patient_genes = FALSE) {
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
cell_type = "CD14 Mono"
umpa_data = read_tsv("tables/supplemental/umap_patient_level.tsv")

make_dotplot = function(data, cell_type, patient_genes = FALSE) {
  
  interesting_genes = collect_interesting_genes(cell_type, patient_genes)
  features = interesting_genes$Gene %>% unique()
  
  DefaultAssay(data) = "RNA"
  Idents(data) = "celltype.de"
  cell_type_data = subset(data, idents = cell_type)
  
  all_pathway_genes = c("id", "Gene", "avg.exp.scaled", "pct.exp", "Pathway", "orig.ident", "MUTATION.GROUP", "STIM")
  names(all_pathway_genes) = c("id", "Gene", "avg.exp.scaled", "pct.exp", "Pathway", "orig.ident", "MUTATION.GROUP", "STIM")
  
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
      dplyr::select(c("id", "Gene", "avg.exp.scaled", "pct.exp", "Pathway")) %>% 
      inner_join(cell_type_data@meta.data %>% dplyr::select(orig.ident, STIM, MUTATION.GROUP) %>% dplyr::mutate(id = paste(orig.ident, STIM, sep = "_"))) %>%
      unique()
    
    if (is.null(nrow(all_pathway_genes))) {
      all_pathway_genes = pathway_genes
    } else {
      all_pathway_genes = rbind(all_pathway_genes, pathway_genes)
    }
  }
  
  ordered_genes = all_pathway_genes %>%
    filter(grepl("STIM", id)) %>%
    group_by(Gene) %>%
    summarize(Avg_Exp = mean(avg.exp.scaled), Pct_Exp = mean(pct.exp)) %>%
    arrange(Avg_Exp, desc(Pct_Exp)) %>%
    mutate(Gene = factor(Gene, levels = Gene)) %>%
    pull(Gene)
  
  stim_dotplot_data = all_pathway_genes %>%
    mutate(`Mutation group` = factor(MUTATION.GROUP, levels = c("Control", "TET2", "DNMT3A"))) %>%
    mutate(Gene = factor(Gene, levels = c(ordered_genes))) %>%
    arrange(desc(`Mutation group`)) %>%
    mutate(id = factor(id, levels = id %>% unique())) %>%
    rename(`Percent expression` = pct.exp) %>%
    rename(`Average expression scaled` = avg.exp.scaled)
  
  stim_dotplot_data = inner_join(stim_dotplot_data, umap_data %>% select(orig.ident, vaf_identifier, STIM) %>% unique()) %>%
    rowwise() %>%
    mutate(identifier = paste(strsplit(vaf_identifier, ":")[[1]][1], STIM)) %>%
    arrange(desc(identifier)) %>%
    mutate(identifier = factor(identifier, levels = identifier %>% unique()))
  
  (dotplot_compiled = ggplot(stim_dotplot_data, aes(x = Gene, y = identifier)) +
      scale_color_manual(values = wes_palette("Royal1", 20, "continuous")[c(18,1,7)]) +
      geom_point(aes(size = `Percent expression`, alpha = `Average expression scaled`, color = `Mutation group`)) +
      theme_bw() +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text = element_text(size = rel(1.2), face = "bold"),
            strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_hline(yintercept=seq(1.5, length(unique(stim_dotplot_data$identifier))-0.5, 1), 
                 lwd=0.5, colour="grey") +
      ggtitle(paste(gsub("Mono", "Monocytes", cell_type))))
  
  filename = paste0(output_dir, "all_stim_patients_dotplot_", gsub(" ", "_", cell_type))
  write_tsv(stim_dotplot_data, "tables/supplemental/stim_dotplot_data.tsv")
  save_plot(filename, dotplot_compiled, height = 6, width = 8, png = TRUE)
}


make_dotplot(singlet_data, "CD14 Mono")


