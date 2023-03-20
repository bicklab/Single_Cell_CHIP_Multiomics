
# Pathway analysis ------------------------------------------------------

source("~/Alyssa/Single_Cell_CHIP/scripts/load_libraries.R")
source("~/Alyssa/Single_Cell_CHIP/scripts/image_formatting.R")

output_dir = "figures/supplemental/"

# Load cluster profiler input
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
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

    # save plots
    filename = paste0(output_dir, "go/dotplot_", sample_genotype, "_", cell_type, ifelse(use_cutoff, "", "_no_cutoff"), "_", stim_status)
    title = paste(sample_genotype, gsub("_", " ", cell_type))
    (plot = dotplot(gse) + 
        format_colors_desc + 
        ggtitle(title))
    save_plot(filename, plot, height = 6, width = 6, png = TRUE)
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

# Pull out gene function from GSE object
# testing with gse_TET2_CD14_Mono

# BiocManager::install("GO.db")
# library(GO.db)
# GO = as.list(GOTERM)
# 
# go_map = data.frame(ID = NA, Description = NA)
# for (item in GO) {
#   go_map[nrow(go_map) + 1,] = c(item@GOID, item@Term)
# }
# go_map = na.omit(go_map)
# 
# gse_TET2_CD14_Mono = readRDS("rds_objects/gse_TET2_CD14_Mono.rds")
# gse_genesets = gse_TET2_CD14_Mono@geneSets %>%
#   unlist() %>%
#   enframe() %>%
#   inner_join(go_map, by = c("name" = "ID")) %>%
#   rename(Gene = "value")
# 
# BiocManager::install("biomaRt")
# library(biomaRt)
# 
# install.packages("readGAF")
# 
# hg19 = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#                 host = "https://grch37.ensembl.org", 
#                 path = "/biomart/martservice",
#                 dataset = "hsapiens_gene_ensembl")
# 
# ids_and_gene_names = getBM(attributes = c("hgnc_symbol", "go_id"),
#                filters = "go",
#                values = gse_genesets$name,
#                mart = hg19)
# 
# complete_table = inner_join(go_map, ids_and_gene_names, by = c("ID" = "go_id"))
# write_tsv(complete_table, "input_data/GO_reference.tsv")
# 
# input_filename = "tables/supplemental/TET2_CD14_Mono_unstim_wilcox.csv"
# input_data = read_csv(input_filename)
# 
# exp_data = input_data %>%
#   left_join(complete_table, by = c("Gene" = "hgnc_symbol"))  %>%
#   group_by(Gene) %>%
#   mutate(Description = paste(Description, collapse = ", ")) %>%
#   dplyr::select(-ID) %>%
#   unique()

