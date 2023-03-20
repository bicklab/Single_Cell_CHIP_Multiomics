install.packages("table1")

library(boot) 
library(table1)

melanoma2 <- melanoma

# Factor the basic variables that
# we're interested in
melanoma2$status <- 
  factor(melanoma2$status, 
         levels=c(2,1,3),
         labels=c("Alive", # Reference
                  "Melanoma death", 
                  "Non-melanoma death"))

table1(~ factor(sex) + age + factor(ulcer) + thickness | status, data=melanoma2)

# Dotplot plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/load_libraries.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")

clinical_data = openxlsx::read.xlsx("input_data/clinical data for brett v3.xlsx")

matched_mutations = read_tsv("tables/supplemental/matched_mutations.tsv")
grouped_clinical_data = inner_join(matched_mutations, clinical_data, by = c("orig.ident" = "patient_id"))

grouped_clinical_data$former_smoker[grouped_clinical_data$former_smoker == "yes"] = TRUE
grouped_clinical_data$former_smoker[grouped_clinical_data$former_smoker == "no"] = FALSE
grouped_clinical_data$former_smoker = as.logical(grouped_clinical_data$former_smoker)

grouped_clinical_data$hx_cad[grouped_clinical_data$hx_cad == "yes"] = TRUE
grouped_clinical_data$hx_cad[grouped_clinical_data$hx_cad == "no"] = FALSE
grouped_clinical_data$hx_cad = as.logical(grouped_clinical_data$hx_cad)

grouped_clinical_data$hx_CHF[grouped_clinical_data$hx_CHF == "yes"] = TRUE
grouped_clinical_data$hx_CHF[grouped_clinical_data$hx_CHF == "no"] = FALSE
grouped_clinical_data$hx_CHF = as.logical(grouped_clinical_data$hx_CHF)

continuous_mean <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("", "Mean (SD)"=sprintf("%s (%s)", MEAN, SD)))
}

rndr <- function(x, ...) {
  y <- render.default(x, ...)
  if (is.logical(x)) y[2] else y
}

(table1 = table1(~ Sex + Race + Age + `Hx tobacco use` + `Hx CAD` + `Hx CHF` +
                   `Systolic BP` + `Diastolic BP` + `Heart rate` | MUTATION.GROUP, data = grouped_clinical_data %>%
                   rename("Sex" = sex, "Race" = race, "Age" = age_initial_sample_collection, "Hx tobacco use" = former_smoker, "Hx CAD" = hx_cad, "Hx CHF" = hx_CHF, "Systolic BP" = systolic_bp, "Diastolic BP" = diastolic_bp, "Heart rate" = heart_rate), 
                 overall = NULL,
                 render.continuous = continuous_mean,
                 render = rndr))

# ULTIMATELY SAVED THE HTML FILE USING EXPORT FROM THE POPUP AND THEN CONVERTED TO PDF WITH ADOBE

blood_counts = grouped_clinical_data %>%
  select(c(orig.ident, wbc, hgb, hct, plt, abs_neutrophils, `neutrophils_%`, abs_lymphs, `lymphocytes_%`, abs_monocytes, `monocytes_%`, 
           abs_eosinophils, `eosinophils_%`, abs_basophils, `basophils_%`, abs_imm_granulocyte, `imm_granulocyte_%`))
write_tsv(blood_counts, "tables/supplemental/blood_counts.tsv")


# pdf("tables/supplemental/table1.pdf")
# rmarkdown::pandoc_convert("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/tables/supplemental/table1.html", to = "pdf")
# dev.off()
# 
# html_to_pdf <- function(html_file, pdf_file) {
#   cmd <- sprintf("pandoc %s -t latex -o %s", html_file, pdf_file)
#   system(cmd)
# }
# 
# html_to_pdf("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/tables/supplemental/table1.html", "tables/supplemental/table1.pdf")
# 
# 
# html_path = "/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/tables/supplemental/table1.html"
# pdf_path = "tables/supplemental/table1.pdf"
# 
# library(psycModel)
# html_to_pdf(html_path, pdf_path)
