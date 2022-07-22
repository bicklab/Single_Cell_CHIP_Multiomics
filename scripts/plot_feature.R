
# Feature plots -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_singlet_data.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

output_dir = "figures/supplemental/"

Idents(singlet_data) = "STIM"
unstim_data = subset(singlet_data, idents = "VEH")

table(singlet_data@meta.data$celltype.de)
features = c("IL7R", "LYZ", "CD14", "MS4A1", "CD8A", "GNLY", "NKG7", "CST3", "PPBP")
plots = FeaturePlot(singlet_data, features = features, combine = FALSE)

for (i in 1:length(plots)) {
  if (i != 6) {
    plots[[i]] = plots[[i]] + NoLegend() + format_no_axes
  } else {
    plots[[i]] = plots[[i]] + format_no_axes
  }
}

feature_plot = cowplot::plot_grid(plotlist = plots)
filename = paste0(output_dir, "feature_canonical")
save_plot(filename, feature_plot, height = 8, width = 8)
