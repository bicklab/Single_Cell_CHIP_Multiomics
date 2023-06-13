library(RColorBrewer)
library(ggplot2)

# Image formatting --------------------------------------------------------


save_plot = function(plot, file, legend = TRUE, width = 4, height = 4, png = FALSE, transparent = FALSE) {
  if (png) { 
    ggplot2::ggsave(plot = plot, filename =  paste0(file, ".png"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in') 
  }
  ggplot2::ggsave(plot = plot, filename =  paste0(file, ".pdf"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in')
}

# Colors light to dark  ------------------------------------------------

format_colors = list(scale_color_gradient(low = brewer.pal(3, "Paired")[1], high = brewer.pal(3, "Paired")[2]))

# Colors dark to light  ------------------------------------------------

format_colors_desc = list(ggplot2::scale_color_gradient(low = brewer.pal(3, "Paired")[2], high = brewer.pal(3, "Paired")[1]))
format_colors_desc_tet2 = list(ggplot2::scale_color_gradient(low = brewer.pal(6, "Paired")[4], high = brewer.pal(6, "Paired")[3]))
format_colors_desc_dnmt3a = list(ggplot2::scale_color_gradient(low = brewer.pal(6, "Paired")[6], high = brewer.pal(6, "Paired")[5]))


# No axes -----------------------------------------------------------------

format_no_axes = list(scale_color_gradient(low = "#f0f0f0", high = brewer.pal(3, "Paired")[2], limits = c(0,6)),
                      theme(axis.line = element_blank(), 
                            axis.ticks = element_blank(), 
                            axis.text = element_blank(),
                            axis.title = element_blank()))

# No axes - faceted --------------------------------------------------------

format_no_axes_faceted = list(
  theme_bw(),
  guides(color = "none"),
  scale_color_manual(values = brewer.pal(3, "Paired")),
  theme(
    strip.text = element_text(size = rel(1.2)),
    strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank()),
  Seurat::NoAxes())

# Arrows -----------------------------------------------------------------

format_add_arrows = list(
  ggplot2::theme_bw(),
  ggplot2::guides(alpha = "none"),
  ggplot2::scale_fill_gradient(low = brewer.pal(3, "Paired")[1], high = brewer.pal(3, "Paired")[2]), # change the low to an alternate color for nonfaceting
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
    strip.background = ggplot2::element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(), 
    panel.border = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()),
  Seurat::NoAxes(),
  ggplot2::annotate(x=-16, xend=-16, y=-16, yend=-8, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate(x=-16, xend=-11, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate("text", label="UMAP 1", x=-17, -13, size = 4, angle = 90),
  ggplot2::annotate("text", label="UMAP 2", x=-14, y=-18, size = 4))

# Just arrows, no color
format_add_just_arrows = list(
  ggplot2::theme_bw(),
  ggplot2::guides(alpha = "none"),
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
    strip.background = ggplot2::element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(), 
    panel.border = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()),
  Seurat::NoAxes(),
  ggplot2::annotate(x=-16, xend=-16, y=-16, yend=-8, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate(x=-16, xend=-11, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate("text", label="UMAP 1", x=-17, -13, size = 4, angle = 90),
  ggplot2::annotate("text", label="UMAP 2", x=-14, y=-18, size = 4))


# Arrows 2 -----------------------------------------------------------------

format_add_arrows_2 = list(
  scale_fill_gradient(low = brewer.pal(3, "Paired")[1], high = brewer.pal(3, "Paired")[2]) # change the low to an alternate color for nonfaceting
)

# Arrows 3 -----------------------------------------------------------------

format_add_arrows_3 = list(
  scale_fill_gradient(low = brewer.pal(3, "Paired")[1], high = brewer.pal(3, "Paired")[2]) # change the low to an alternate color for nonfaceting
  
)

# Volcano -----------------------------------------------------------------

colors = c(brewer.pal(4, "Paired"), "#d9d9d9")
names(colors) = c("Inflammation", "Antigen presentation", "Cell adhesion", "Monocyte activation", "Other")
format_volcano = list(
  xlab("log2(fold change)"),
  ylab("-log10(p-value)"),
  labs(color = "Pathway"),
  theme_bw(),
  theme_classic(),
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12.5),
        title = element_text(size = 11),
        legend.key.size = unit(0, "lines")),
  scale_color_manual(values = colors)
)

# Volcano 2 -----------------------------------------------------------------

colors_2 = c(brewer.pal(4, "Paired"), "#d9d9d9")
names(colors_2) = c("T cell activation", "Intracellular signaling", "Immunomodulation", "Other")
format_volcano_2 = list(
  xlab("log2(fold change)"),
  ylab("-log10(p-value)"),
  labs(color = "Pathway"),
  theme_bw(),
  theme_classic(),
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12.5),
        title = element_text(size = 11),
        legend.key.size = unit(0, "lines")),
  scale_color_manual(values = colors_2)
)

# Volcano 3 -----------------------------------------------------------------

colors_3 = c(brewer.pal(4, "Paired"), "#d9d9d9") 
names(colors_3) = c("T cell activation/TCR signaling", "Intracellular signaling", "Immune recruitment", "Other")
format_volcano_3 = list(
  xlab("log2(fold change)"),
  ylab("-log10(p-value)"),
  labs(color = "Pathway"),
  theme_bw(),
  theme_classic(),
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 12.5),
        title = element_text(size = 11),
        legend.key.size = unit(0, "lines")),
  scale_color_manual(values = colors_3)
)


