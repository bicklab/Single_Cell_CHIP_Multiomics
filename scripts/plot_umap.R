
# UMAP plots ---------------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")


output_dir = "figures/supplemental/"

control_count = 0
dnmt3a_count = 0
tet2_count = 0

get_vaf = function(mutation, group) {
  results = str_extract_all(mutation, paste0(group, "[^)]*\\)"))
  percent_vector = lapply(results, str_extract_all, pattern = "(?<=\\().*%") %>% unlist()
  if (length(percent_vector) > 0) {
    max = substr(percent_vector, 1, nchar(percent_vector)-1) %>% as.numeric() %>% max()
    return(max)
  } else {
    return(0)
  }
}

get_vaf_identifier = function(vaf, group) {
  if (vaf > 0) {
    if (group == "DNMT3A") {
      dnmt3a_count <<- dnmt3a_count + 1 # modifies a global variable
      return(paste0(dnmt3a_count, ": ", vaf / 100))
    } else {
      tet2_count <<- tet2_count + 1 # modifies a global variable
      return(paste0(tet2_count, ": ", vaf / 100))
    }
  } else {
    control_count <<- control_count + 1 # modifies a global variable
    return(control_count)
  }
}

vaf_identifiers = singlet_data@meta.data %>%
  dplyr::select(orig.ident, MUTATION, MUTATION.GROUP) %>%
  unique() %>%
  rowwise() %>%
  mutate(VAF = get_vaf(MUTATION, MUTATION.GROUP)) %>%
  arrange(VAF) %>%
  mutate(vaf_identifier = paste(MUTATION.GROUP, get_vaf_identifier(VAF, MUTATION.GROUP)))

umap_patient_level_data = singlet_data@meta.data %>% 
  inner_join(vaf_identifiers)

(umap_patient_level_plot = ggplot(umap_patient_level_data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = `MUTATION.GROUP`), size = 0.4) +
  facet_wrap(~vaf_identifier) +
  format_no_axes_faceted)

filename = paste0(output_dir, "umap_patient_level")
save_plot(filename, umap_patient_level_plot, height = 8, width = 8.5, png = TRUE)
write_tsv(umap_patient_level_data, "tables/supplemental/umap_patient_level.tsv")


