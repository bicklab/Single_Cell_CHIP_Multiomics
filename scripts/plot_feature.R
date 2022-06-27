
# Feature plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

output_dir = "figures/supplemental/"

Idents(singlet_data) = "STIM"
unstim_data = subset(singlet_data, idents = "VEH")

table(singlet_data@meta.data$celltype.de)

features = c("IL7R", "LYZ", "CD14", "MS4A1", "CD8A", "GNLY", "NKG7", "CST3", "PPBP")
#features = c("FCGR3B", "JUND", "TXNIP", "FTH1", "GNLY", "SMIM24", "JAML", "CXCL8", "CD74", "TMSB10", "AIF1", "HIST1H4C")
feature_plot = FeaturePlot(singlet_data, features = features) &
  format_no_axes

filename = paste0(output_dir, "feature_canonical")
save_plot(filename, feature_plot)
