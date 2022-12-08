
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")

options("digits" = 4)

clinical_data = read.xlsx("input_data/clinical data for brett v5.xlsx")

matched_mutations = singlet_data@meta.data %>%
  dplyr::select(orig.ident, MUTATION.GROUP) %>%
  unique()
matched_mutations$orig.ident[matched_mutations$orig.ident == "CH-21-001"] = "CH-20-001"

write_tsv(matched_mutations, "tables/supplemental/matched_mutations.tsv")

range_and_median = function(data) {
  return(paste0(median(data, na.rm = TRUE) %>% round(2), " (", min(data, na.rm = TRUE), " - ", max(data, na.rm = TRUE), ")"))
}

is_nondialysis = function(x) {
  if (x < 10) return(TRUE) else return(FALSE)
}

# standardize_ef = function(x) {
#   if (is.na(x)) {
#     return(x)
#   } else if (grepl("-", x)) {
#     values = str_split(x, "-")
#     values = lapply(values, as.numeric) %>% unlist()
#     return(mean(values))
#   } else if (grepl(">", x)) {
#     return(str_replace_all(x, ">", "") %>% as.numeric())
#   } else if (x %>% as.numeric() < 1) {
#     return(x %>% as.numeric() * 100)
#   } else {
#     return(x %>% as.numeric())
#   }
# }

grouped_clinical_data = inner_join(matched_mutations, clinical_data, by = c("orig.ident" = "patient_id")) %>%
  group_by(MUTATION.GROUP) %>%
  summarise(sex_male = paste0(sum(sex == "Male"), ":", sum(sex == "Female")),
            race_white = paste0(sum(race == "White", na.rm = TRUE), ":", sum(race != "White", na.rm = TRUE)),
            age = range_and_median(age_initial_sample_collection),
            
            hx_tobacco_use = paste0(sum(former_smoker == "yes", na.rm = TRUE), ":", sum(former_smoker == "no", na.rm = TRUE)),
            hx_CAD = paste0(sum(hx_cad == "yes", na.rm = TRUE), ":", sum(hx_cad == "no", na.rm = TRUE)),
            hx_cardiomyopathy = paste0(sum(hx_cardiomyopathy == "yes", na.rm = TRUE), ":", sum(hx_cardiomyopathy == "no", na.rm = TRUE)),
            hx_chf = paste0(sum(hx_CHF == "yes", na.rm = TRUE), ":", sum(hx_CHF == "no", na.rm = TRUE)),
            
            systolic_bp = range_and_median(systolic_bp),
            diastolic_bp = range_and_median(diastolic_bp),
            heart_rate = range_and_median(heart_rate),
            
            wbc = range_and_median(wbc),
            hgb = range_and_median(hgb),
            hct = range_and_median(hct),
            plt = range_and_median(plt),
            #neutrophil_pct = range_and_median(`neutrophils_%`),
            neutrophil_abs = range_and_median(`abs_neutrophils`),
            #lymphocyte_pct = range_and_median(`lymphocytes_%`),
            lymphocyte_abs = range_and_median(`abs_lymphs`),
            #monocyte_pct = range_and_median(`monocytes_%`),
            monocyte_abs = range_and_median(`abs_monocytes`),
            #eosinophil_pct = range_and_median(`eosinophils_%`),
            eosinophil_abs = range_and_median(`abs_eosinophils`),
            #basophil_pct = range_and_median(`basophils_%`),
            basophil_abs = range_and_median(`abs_basophils`),
            #granulocyte_pct = range_and_median(`imm_granulocyte_%`),
            granulocyte_abs = range_and_median(`abs_imm_granulocyte`),
            BUN = range_and_median(bun),
            creatinine = range_and_median(Filter(is_nondialysis, creatinine)),
            glucose = range_and_median(glucose),
            calcium = range_and_median(ca),
            cholesterol = range_and_median(total_cholesterol),
            HDL = range_and_median(hdl),
            LDL = range_and_median(ldl),
            triglycerides = range_and_median(triglycerides)) %>% #,
            #ef_pct = range_and_median(lapply(`EF_%`, standardize_ef) %>% unlist() %>% as.numeric())) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Category")

write_tsv(grouped_clinical_data, "tables/supplemental/grouped_clinical_data.tsv")
