
# Test ALRA for imputation ------------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")

remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

Idents(singlet_data) = "celltype.de"
cd14_data = subset(singlet_data, idents = "CD14 Mono")
alra_cd14 = RunALRA(cd14_data)

# visualize original and imputed values
features.plot = c("CD3D", "MS4A1", "CD8A", "GZMK", "NCAM1", "FCGR3A")
DefaultAssay(cd14_data) = "RNA"
plot1 = FeaturePlot(cd14_data, features.plot, ncol = 2)

DefaultAssay(alra_cd14) = "alra"
plot2 = FeaturePlot(alra_cd14, features.plot, ncol = 2, cols = c("lightgrey", "red"))
CombinePlots(list(plot1, plot2), ncol = 1)

cell_types = c("CD14 Mono")
sample_genotypes = c("TET2", "DNMT3A")
x_chr_genes = read_csv("input_data/x_genes.csv")
output_dir = "tables/supplemental/"

cell_type = "CD14 Mono"
sample_genotype = "TET2"

run_stats = function(data, stim_status) {
  for (cell_type in cell_types) {
    for (sample_genotype in sample_genotypes) {
      
      # Identify case and control data
      if (sample_genotype == "CHIP") {
        case_idents = c("TET2", "DNMT3A")
      } else {
        case_idents = sample_genotype
      }
      
      control_idents = "Control"
      method = "wilcox"
      Idents(data) = "celltype.de"
      
      # Run statistics
      response = FindMarkers(data, 
                             subset.ident = cell_type,
                             group.by = "MUTATION.GROUP",
                             ident.1 = case_idents,
                             ident.2 = control_idents,
                             min.pct = 0,
                             min.diff.pct = -Inf,
                             logfc.threshold = 0.25,
                             test.use = method,
                             assay = "alra")
      
      # Arrange data and remove X chromosome genes
      clean_response = response %>%
        arrange(p_val_adj) %>%
        na.omit() %>%
        rownames_to_column("Gene") %>%
        filter(!(Gene %in% x_chr_genes$Gene))
      
      # Save output
      filename = paste0(output_dir, "ALRA_", sample_genotype, "_", gsub(" ", "_", cell_type), "_", stim_status, "_", method, ".csv")
      write_csv(clean_response, file = filename)
    }
  }
}

Idents(alra_cd14) = "STIM"
alra_cd14_stim = subset(alra_cd14, idents = "STIM")
alra_cd14_unstim = subset(alra_cd14, idents = "VEH")


plan(multisession)
options('future.globals.maxSize' = 8000*1024^2)
future(run_stats(alra_cd14_unstim, "unstim"))
future(run_stats(alra_cd14_stim, "stim"))

data = alra_cd14_stim
stim_status = "stim"
