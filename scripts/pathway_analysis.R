
# Pathway analysis ------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

output_dir = "figures/supplemental/"

# Load cluster profiler input
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

find_pathways = function(sample_genotype, cell_type, stim_status, use_cutoff = TRUE) {
  # load expression data from stats output
  folder = "tables/supplemental/"
  input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_", stim_status))
  exp_data = read_csv(paste0(folder, input_filename))

  # format gene list
  gene_list = exp_data %>%
    dplyr::select(Gene, avg_log2FC) %>%
    arrange(desc(avg_log2FC)) %>%
    deframe()
  
  # written with trycatch to avoid errors with empty lists
  tryCatch({
    # identify enriched pathways
    gse = gseGO(geneList = gene_list, keyType = "SYMBOL", OrgDb = organism, ont = "ALL", minGSSize = 2, maxGSSize = 800, pvalueCutoff = ifelse(use_cutoff, 0.05, 1))
    saveRDS(gse, paste0("rds_objects/gse_", sample_genotype, "_", cell_type, ".rds"))

    # save plots
    filename = paste0(output_dir, "go/dotplot_", sample_genotype, "_", cell_type, ifelse(use_cutoff, "", "_no_cutoff"), "_", stim_status)
    title = paste(sample_genotype, gsub("_", " ", cell_type))
    (plot = dotplot(gse %>% filter(NES > 0), split = ".sign") + 
        format_colors_desc + 
        ggtitle(title) + 
        facet_grid(.~.sign) +
        theme(strip.background = element_blank(),
              strip.text.x = element_blank()))
    save_plot(filename, plot, height = 6, width = 6)
  },
  error = function(cond) {
    message(cond)
    return("Error rendering plot.")
  })
}

sample_genotypes = c("DNMT3A", "TET2")
cell_types = c("CD14_Mono", "CD8_T_cell")

for (sample_genotype in sample_genotypes) {
  for (cell_type in cell_types) {
    find_pathways(sample_genotype, cell_type, stim_status = "unstim")
  }
}

