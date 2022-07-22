
# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")

options("digits" = 4)

clinical_data = read.xlsx("input_data/clinical data for brett v3.xlsx")

matched_mutations = singlet_data@meta.data %>%
  select(orig.ident, MUTATION.GROUP) %>%
  unique()

range_and_median = function(data) {
  return(paste0(median(data, na.rm = TRUE) %>% round(2), " (", min(data, na.rm = TRUE), " - ", max(data, na.rm = TRUE), ")"))
}

is_nondialysis = function(x) {
  if (x < 10) return(TRUE) else return(FALSE)
}

standardize_ef = function(x) {
  if (is.na(x)) {
    return(x)
  } else if (grepl("-", x)) {
    values = str_split(x, "-")
    values = lapply(values, as.numeric) %>% unlist()
    return(mean(values))
  } else if (grepl(">", x)) {
    return(str_replace_all(x, ">", "") %>% as.numeric())
  } else if (x %>% as.numeric() < 1) {
    return(x %>% as.numeric() * 100)
  } else {
    return(x %>% as.numeric())
  }
}

grouped_clinical_data = inner_join(matched_mutations, clinical_data, by = c("orig.ident" = "patient_id")) %>%
  group_by(MUTATION.GROUP) %>%
  summarise(sex_male = paste0(sum(sex == "Male"), ":", sum(sex == "Female")),
            race_white = paste0(sum(race == "White", na.rm = TRUE), ":", sum(race != "White", na.rm = TRUE)),
            former_smoker = paste0(sum(former_smoker == "yes", na.rm = TRUE), ":", sum(former_smoker == "no", na.rm = TRUE)),
            hx_cad = paste0(sum(hx_cad == "yes", na.rm = TRUE), ":", sum(hx_cad == "no", na.rm = TRUE)),
            hx_cardiomyopathy = paste0(sum(hx_cardiomyopathy == "yes", na.rm = TRUE), ":", sum(hx_cardiomyopathy == "no", na.rm = TRUE)),
            hx_CHF = paste0(sum(hx_CHF == "yes", na.rm = TRUE), ":", sum(hx_CHF == "no", na.rm = TRUE)),
            median_age = range_and_median(age_initial_sample_collection),
            median_systolic_bp = range_and_median(systolic_bp),
            median_diastolic_bp = range_and_median(diastolic_bp),
            median_heart_rate = range_and_median(heart_rate),
            median_wbc = range_and_median(wbc),
            median_hgb = range_and_median(hgb),
            median_hct = range_and_median(hct),
            median_plt = range_and_median(plt),
            median_neutrophil_pct = range_and_median(`neutrophils_%`),
            median_neutrophil_abs = range_and_median(`abs_neutrophils`),
            median_lymphocyte_pct = range_and_median(`lymphocytes_%`),
            median_lymphocyte_abs = range_and_median(`abs_lymphs`),
            median_monocyte_pct = range_and_median(`monocytes_%`),
            median_monocyte_abs = range_and_median(`abs_monocytes`),
            median_eosinophil_pct = range_and_median(`eosinophils_%`),
            median_eosinophil_abs = range_and_median(`abs_eosinophils`),
            median_basophil_pct = range_and_median(`basophils_%`),
            median_basophil_abs = range_and_median(`abs_basophils`),
            median_granulocyte_pct = range_and_median(`imm_granulocyte_%`),
            median_granulocyte_abs = range_and_median(`abs_imm_granulocyte`),
            median_bun = range_and_median(bun),
            median_creatinine = range_and_median(Filter(is_nondialysis, creatinine)),
            median_glucose = range_and_median(glucose),
            median_ca = range_and_median(ca),
            median_cholesterol = range_and_median(total_cholesterol),
            median_hdl = range_and_median(hdl),
            median_ldl = range_and_median(ldl),
            median_triglycerides = range_and_median(triglycerides),
            median_ef_pct = range_and_median(lapply(`EF_%`, standardize_ef) %>% unlist() %>% as.numeric())) %>%
  
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Category")

write_tsv(grouped_clinical_data, "tables/supplemental/grouped_clinical_data.tsv")
