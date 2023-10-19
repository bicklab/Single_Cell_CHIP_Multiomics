# Review CCUS vs control and CHIP vs control results.
setwd("Alyssa/Single_Cell_CHIP/")

# CD14 Monos

chip_cd14_tet2_results = read_tsv("tables/metacells_de_cell_types_filtered_ccus/CD14 Monos_TET2_CHIP_VEH_vs_none_Control_VEH.tsv")
ccus_cd14_tet2_results = read_tsv("tables/metacells_de_cell_types_filtered_ccus/CD14 Monos_TET2_CCUS_VEH_vs_none_Control_VEH.tsv")

# 752 sig genes
sig_chip = chip_cd14_tet2_results %>%
  filter(padj < 0.05)

# 1036 sig genes
sig_ccus = ccus_cd14_tet2_results %>%
  filter(padj < 0.05)

# 437 overlap
intersect(sig_chip$gene, sig_ccus$gene)


# 671
pos_chip = chip_cd14_tet2_results %>%
  filter(log2FoldChange > 0)
# 848
pos_ccus = ccus_cd14_tet2_results %>%
  filter(log2FoldChange > 0)

# 362 overlap
intersect(pos_chip$gene, pos_ccus$gene) %>% length()

# 614
neg_chip = chip_cd14_tet2_results %>%
  filter(log2FoldChange < 0)
# 915
neg_ccus = ccus_cd14_tet2_results %>%
  filter(log2FoldChange < 0)

# 367 intersect
intersect(neg_chip$gene, neg_ccus$gene) %>% length()


intersect(chip_cd14_tet2_results$gene, ccus_cd14_tet2_results$gene) %>% length()



# 479
sig_pos_chip = chip_cd14_tet2_results %>%
  filter(log2FoldChange > 0) %>%
  filter(padj < 0.05)
# 402
sig_pos_ccus = ccus_cd14_tet2_results %>%
  filter(log2FoldChange > 0) %>%
  filter(padj < 0.05)

# 142 overlap
intersect(sig_pos_chip$gene, sig_pos_ccus$gene) %>% length()



# 350
sig_neg_chip = chip_cd14_tet2_results %>%
  filter(log2FoldChange < 0) %>%
  filter(padj < 0.05)
# 557
sig_neg_ccus = ccus_cd14_tet2_results %>%
  filter(log2FoldChange < 0) %>%
  filter(padj < 0.05)

# 156 overlap
intersect(sig_neg_chip$gene, sig_neg_ccus$gene) %>% length()

