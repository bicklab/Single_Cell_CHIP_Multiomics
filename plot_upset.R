
# Upset plots ---------------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/figure_1/"

patient_level_data = read_tsv("tables/supplemental/umap_patient_level.tsv")

get_mutations = function(mutation) {
  mutations = str_split(mutation, ", ")
  mutation_vector = lapply(mutations, str_extract, pattern = "[^ ]+")
  return(mutation_vector)
}

get_vafs = function(mutation) {
  mutations = str_split(mutation, ", ")
  mutation_vector = lapply(mutations, str_extract, pattern = "\\([^%]+")
  clean_mutation_vector = lapply(mutation_vector, substring, 2)
  numeric_mutation_vector = lapply(clean_mutation_vector, as.numeric)
  return(clean_mutation_vector)
}

mutation_data = patient_level_data %>%
  dplyr::select(orig.ident, MUTATION, VAF) %>%
  unique() %>%
  rowwise() %>%
  mutate(Mutations = get_mutations(MUTATION)) %>%
  mutate(VAFs = get_vafs(MUTATION)) %>%
  mutate(Mutation_Group = paste(Mutations %>% unique() %>% sort(), collapse = " "))

write_tsv(mutation_data, "tables/supplemental/mutation_data.tsv")


mutation_list = mutation_data$Mutations %>% lapply(paste, sep = "&")

# list of sets, each set is a mutation
unique_mutations = mutation_data$Mutations %>% unlist() %>% unique()

tet2_data = c()
srsf2_data = c()
dnmt3a_data = c()
idh2_data = c()
tp53_data = c()
none_data = c()

tet2_vaf = c()
srsf2_vaf = c()
dnmt3a_vaf = c()
idh2_vaf = c()
tp53_vaf = c()
none_vaf = c()

full_data = c("Entry", "VAF", "Mutation_Group")
  
  
for (i in 1:length(mutation_data$Mutations)) {
  entry = mutation_data$Mutations[i] %>% unlist()
  vafs = mutation_data$VAFs[i] %>% unlist()
  mutation_group = mutation_data$Mutation_Group[i]
  
  for (j in 1:length(entry)) {
    sub_entry = entry[j]
    vaf = vafs[j]
    
    new_line = c(sub_entry, as.numeric(vaf), mutation_group)
    full_data = rbind(full_data, new_line)
    
    if ("TET2" == sub_entry) {
      tet2_data = c(tet2_data, i)
      tet2_vaf = c(tet2_vaf, vaf)
    }
    if ("SRSF2" == sub_entry) {
      srsf2_data = c(srsf2_data, i)
      srsf2_vaf = c(srsf2_vaf, vaf)
    }
    if ("DNMT3A" == sub_entry) {
      dnmt3a_data = c(dnmt3a_data, i)
      dnmt3a_vaf = c(dnmt3a_vaf, vaf)
    }
    if ("IDH2" == sub_entry) {
      idh2_data = c(idh2_data, i)
      idh2_vaf = c(idh2_vaf, vaf)
    }
    if ("TP53" == sub_entry) {
      tp53_data = c(tp53_data, i)
      tp53_vaf = c(tp53_vaf, vaf)
    }
    if ("none" == sub_entry) {
      none_data = c(none_data, i)
      none_vaf = c(none_vaf, "0")
    }
  }
}

list_mutation_data = list(
  "TET2" = tet2_data, "SRSF2" = srsf2_data, "DNMT3A" = dnmt3a_data, "IDH2" = idh2_data, "TP53" = tp53_data, "None" = none_data
)
list_vaf_data = list(
  "TET2" = tet2_vaf, "SRSF2" = srsf2_vaf, "DNMT3A" = dnmt3a_vaf, "IDH2" = idh2_vaf, "TP53" = tp53_vaf, "None" = none_vaf
)


(upset_plot = upset(fromList(list_mutation_data), nsets = 7, keep.order = T, mb.ratio = c(0.3, 0.7),
      sets = c("None", "TP53", "IDH2", "SRSF2", "TET2", "DNMT3A"), text.scale = 1.5,
      main.bar.color = wes_palette("Royal1", 8, "continuous")[c(1,4,3,7,2)], mainbar.y.label = element_blank(), mainbar.y.max = 11, sets.x.label = element_blank()))

filename = paste0(output_dir, "upset_mutations")
pdf(file = paste0(filename, ".pdf"), width = 4, height = 4)
upset_plot
dev.off()


matrix_mutation_data = plyr::ldply(list_mutation_data, rbind)
write_tsv(matrix_mutation_data, paste0(gsub("figures", "tables", filename), ".tsv"))


colnames(full_data) = full_data[1,]
full_data = full_data[-1,]
full_data[is.na(full_data)] = 0
full_data[full_data == "?ASXL1G645fs"] = "ASXL1"
full_data[full_data == "?ASXL1G645fs TET2"] = "TET2"

full_data = as.data.frame(full_data) %>%
  mutate(VAF = as.numeric(VAF)) %>%
  mutate(Entry = factor(Entry, levels = c("none", "TP53", "IDH2", "SRSF2", "TET2", "DNMT3A"))) %>%
  mutate(Mutation_Group = factor(Mutation_Group, levels = c("DNMT3A", "none", "TET2", "SRSF2 TET2", "DNMT3A IDH2 TP53"))) %>%
  na.omit()


(dot_vafs = ggplot(full_data %>% as.data.frame(), aes(x = Entry, y = VAF)) +
    geom_rect(aes(fill = Entry, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.05), show.legend = FALSE) +
    geom_point(size = 3, aes(color = Mutation_Group)) +
  theme_bw() +
  theme(text = element_text(size = 18), 
        axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.y = element_text(angle = -90, vjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = -90),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
    ylab("VAF") +
    scale_y_reverse() +
    scale_color_manual(values = wes_palette("Royal1", 8, "continuous")[c(1,4,3,7,8,2)]) +
    scale_fill_manual(values = c("#f0f0f0", "#ffffff", "#f0f0f0", "#ffffff", "#f0f0f0", "#ffffff", "#f0f0f0")) +
  facet_grid(~Entry, scales = "free"))
filename = paste0(output_dir, "dot_vafs")
save_plot(filename, dot_vafs, height = 2)
write_tsv(full_data, paste0(gsub("figures", "tables", filename), ".tsv"))
