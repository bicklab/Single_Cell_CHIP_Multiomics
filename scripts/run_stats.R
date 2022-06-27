
# Wilcoxon differential expression analysis  -------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")

singlet_data = readRDS("rds_objects/singlet_data.rds")
output_dir = "tables/supplemental/"

# Prep data ---------------------------------------------------------------

x_chr_genes = read_csv("input_data/x_genes.csv")
gene_df = read_csv("input_data/Gene_class.csv")


# Separate vehicle and stim data
Idents(singlet_data) = "STIM"
veh_data = subset(singlet_data, idents = "VEH")
stim_data = subset(singlet_data, idents = "STIM")

# Prep for stats
#ignore_cells = c("Eryth", "Plasmablast", "Platelet", "HSPC", "MAIT", "gdT", "DC", "CD4 Proliferating", "NK Proliferating", "ILC", "dnT", "CD8 Proliferating")
#cell_types = unique(stim_data@meta.data$celltype.de) %>% setdiff(ignore_cells)
cell_types = "CD8 T cell"
sample_genotypes = c("DNMT3A", "TET2") #CHIP

# Run stats ---------------------------------------------------------------
plan(multisession)
options('future.globals.maxSize' = 6000*1024^2)

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
                             logfc.threshold = -Inf,
                             test.use = method,
                             assay = "RNA")
      
      # Arrange data and remove X chromosome genes
      clean_response = response %>%
        arrange(p_val_adj) %>%
        na.omit() %>%
        rownames_to_column("Gene") %>%
        filter(!(Gene %in% x_chr_genes$Gene))
      
      # Save output
      filename = paste0(output_dir, sample_genotype, "_", gsub(" ", "_", cell_type), "_", stim_status, "_", method, ".csv")
      write_csv(clean_response, file = filename)
    }
  }
}

future(run_stats(veh_data, "unstim"))
future(run_stats(stim_data, "stim"))


# Compile excel file ------------------------------------------------------

sample_genotypes = c("DNMT3A", "TET2")
folder = "tables/supplemental/"

get_cell_data = function(sample_genotype, cell_type) {
  input_filename = list.files(folder, pattern = paste0("^", sample_genotype, "_", cell_type, ".*_unstim"))
  data = read_csv(paste0(folder, input_filename))
  return(data)
}

for (sample_genotype in sample_genotypes) {
  
  nk_data = get_cell_data(sample_genotype, "NK")
  cd8_t_cell_data = get_cell_data(sample_genotype, "CD8_T_cell")
  cd8_naive_data = get_cell_data(sample_genotype, "CD8_Naive")
  cd14_mono_data = get_cell_data(sample_genotype, "CD14_Mono")
  cd4_t_cell_data = get_cell_data(sample_genotype, "CD4_T_cell")
  cd4_naive_data = get_cell_data(sample_genotype, "CD4_Naive")
  cd4_ctl_data = get_cell_data(sample_genotype, "CD4_CTL")
  b_data = get_cell_data(sample_genotype, "B")
  cd16_mono_data = get_cell_data(sample_genotype, "CD16_Mono")
  treg_data = get_cell_data(sample_genotype, "Treg")

  dataframes = list("NK" = nk_data, "CD8_T_cell" = cd8_t_cell_data, "CD8_Naive" = cd8_naive_data,
                    "CD14_Mono" = cd14_mono_data, "CD4_T_cell" = cd4_t_cell_data, "CD4_Naive" = cd4_naive_data,
                    "CD4_CTL" = cd4_ctl_data,
                    "B" = b_data,
                    "CD16_Mono" = cd16_mono_data, "Treg" = treg_data,
                    "gene_sets" = gene_df)
  
  filename = paste0(output_dir, sample_genotype, "_unstim_wilcox_diff_exp.xlsx")
  write.xlsx(dataframes, file = filename)
}





# Stats for stim vs unstim ------------------------------------------------


cell_types = c("CD14 Mono", "CD16 Mono", "CD8 T cell", "CD4 T cell")
run_stim_stats = function(data) {
  for (cell_type in cell_types) {

    case_idents = "STIM"
    control_idents = "VEH"
    method = "wilcox"
    Idents(data) = "celltype.de"
    
    # Run statistics
    response = FindMarkers(data, 
                           subset.ident = cell_type,
                           group.by = "STIM",
                           ident.1 = case_idents,
                           ident.2 = control_idents,
                           min.pct = 0,
                           min.diff.pct = -Inf,
                           test.use = method,
                           assay = "RNA")
      
      # Arrange data and remove X chromosome genes
      clean_response = response %>%
        arrange(p_val_adj) %>%
        na.omit() %>%
        rownames_to_column("Gene") %>%
        filter(!(Gene %in% x_chr_genes$Gene))
      
      # Save output
      filename = paste0(output_dir, "stim_v_unstim_", gsub(" ", "_", cell_type), "_", method, ".csv")
      write_csv(clean_response, file = filename)
  }
}

cell_types = "CD14 Mono"
run_stim_stats(singlet_data)


# Stim vs. unstim stats at patient level -----------------------------------------

run_patient_level_stim_stats = function(data) {
  patients = data@meta.data$orig.ident %>% unique()
  for (cell_type in cell_types) {
      for (patient in patients) {
      
      # Filter to patient
      Idents(data) = "orig.ident"
      data_subset = subset(data, idents = patient)
      
      case_idents = "STIM"
      control_idents = "VEH"
      method = "wilcox"
      Idents(data_subset) = "celltype.de"
      
      # Run statistics
      response = FindMarkers(data_subset, 
                             subset.ident = cell_type,
                             group.by = "STIM",
                             ident.1 = case_idents,
                             ident.2 = control_idents,
                             min.pct = 0,
                             min.diff.pct = -Inf,
                             test.use = method,
                             assay = "RNA")
      # Arrange data and remove X chromosome genes
      clean_response = response %>%
        arrange(p_val_adj) %>%
        na.omit() %>%
        rownames_to_column("Gene") %>%
        filter(!(Gene %in% x_chr_genes$Gene))
      
      # Save output
      filename = paste0(output_dir, "stim_v_unstim_", gsub(" ", "_", cell_type), "_", patient, "_", method, ".csv")
      write_csv(clean_response, file = filename)
    }
  }
}

cell_types = "CD14 Mono"
run_patient_level_stim_stats(singlet_data)

# Stats for tet2 vs. dnmt3a stim and unstim ---------------------------------------------------------------

run_tet2_dnmt3a_stats = function(data, stim_status) {
  for (cell_type in "CD14 Mono") {
    case_idents = "TET2"
    control_idents = "DNMT3A"
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
                           test.use = method,
                           assay = "RNA")
    
    # Arrange data and remove X chromosome genes
    clean_response = response %>%
      arrange(p_val_adj) %>%
      na.omit() %>%
      rownames_to_column("Gene") %>%
      filter(!(Gene %in% x_chr_genes$Gene))
    
    # Save output
    filename = paste0(output_dir, "CHIP_comparison_", gsub(" ", "_", cell_type), "_", stim_status, "_", method, ".csv")
    write_csv(clean_response, file = filename)
  }
}

run_tet2_dnmt3a_stats(veh_data, "unstim")
run_tet2_dnmt3a_stats(stim_data, "stim")


