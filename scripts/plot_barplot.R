
# Barplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/figure_3/"
gene_df = read_csv("input_data/Gene_class.csv")
gene_label_data = read.xlsx("input_data/gene_listp_61322.xlsx")

folder = "tables/supplemental/"

patient_genes = TRUE
cell_type = "CD14 Mono"
stim_status = "unstim"

tet2_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "TET2_", gsub(" ", "_", cell_type), ".*_", stim_status))
dnmt3a_input_filename = list.files(folder, pattern = paste0(ifelse(patient_genes, "^patient_", "^"), "DNMT3A_", gsub(" ", "_", cell_type), ".*_", stim_status))
tet2_significance_data = read_csv(paste0(folder, tet2_input_filename))
dnmt3a_significance_data = read_csv(paste0(folder, dnmt3a_input_filename))


index = c("Inflammation", "Antigen.Presentation", "Cell.Adhesion", "Monocyte.Activation")
values = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation")

pathway_genes = gene_label_data %>%
  pivot_longer(cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
  na.omit()

pathway_genes$Pathway = values[match(pathway_genes$Pathway, index)]
pathway_genes = pathway_genes %>%
  mutate(Pathway = factor(Pathway, levels = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation"))) %>%
  na.omit() %>%
  unique()

# tet2_pathway_genes = gene_label_data %>%
#   pivot_longer(cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
#   na.omit() %>%
#   inner_join(tet2_significance_data, by = "Gene") %>%
#   mutate(Pathway = factor(Pathway, levels = c("Inflammation", "Antigen.Presentation", "Cell.Adhesion", "Monocyte.Activation"))) %>%
#   na.omit() %>%
#   mutate(Group = ifelse(avg_log2FC > 0, "Mutant", "Wildtype")) %>%
#   group_by(Pathway, Group) %>%
#   summarise(Fraction = sum(p_val_adj < 0.05) / n()) %>%
#   arrange(Pathway)
# 
# tet2_pathway_genes$Pathway = factor(values[match(tet2_pathway_genes$Pathway, index)], levels = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation"))
# tet2_pathway_genes[7,] = list("Antigen Presentation", "Wildtype", 0)
# tet2_pathway_genes[8,] = list("Monocyte Activation", "Wildtype", 0)

tet2_significance_data$Genotype = "TET2"
dnmt3a_significance_data$Genotype = "DNMT3A"

all_significance_data = rbind(tet2_significance_data, dnmt3a_significance_data) %>%
  mutate(Significant = p_val_adj < 0.05) %>%
  filter(Gene %in% pathway_genes$Gene)


(barplot_all = ggplot(all_significance_data, aes(x = Genotype)) +
  geom_bar(position = "stack", aes(fill = Significant)) +
  scale_fill_manual(values = wes_palette("Royal1", 3)[c(1,2)]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
        panel.border = element_blank()) +
  xlab(element_blank()) +
  ggtitle("Fraction pathway genes significant"))

filename = paste0(output_dir, "barplot_all_sig_fraction")
save_plot(filename, barplot_all)
write_tsv(all_significance_data, paste0(gsub("figures", "tables", filename), ".tsv"))

