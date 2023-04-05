library(wesanderson)

# Image formatting --------------------------------------------------------


save_plot = function(file, plot, legend = TRUE, width = 4, height = 4, png = FALSE, transparent = FALSE) {
  if (png) { 
    ggplot2::ggsave(plot = plot, filename =  paste0(file, ".png"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in') 
  }
  ggplot2::ggsave(plot = plot, filename =  paste0(file, ".pdf"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in')
}

# Colors light to dark red ------------------------------------------------

format_colors = list(scale_color_gradient(low = wes_palette("Royal1", 20, "continuous")[13], high = wes_palette("Royal1", 2)[2]))

# Colors dark to light red ------------------------------------------------

format_colors_desc = list(ggplot2::scale_color_gradient(low = wes_palette("Royal1", 2)[2], high = wes_palette("Royal1", 20, "continuous")[13]))

# No axes -----------------------------------------------------------------

format_no_axes = list(scale_color_gradient(low = "#f0f0f0", high = wes_palette("Royal1", 2)[2], limits = c(0,6)),
                      theme(axis.line = element_blank(), 
                            axis.ticks = element_blank(), 
                            axis.text = element_blank(),
                            axis.title = element_blank()))

# No axes - faceted --------------------------------------------------------

format_no_axes_faceted = list(
  theme_bw(),
  guides(color = "none"),
  scale_color_manual(values = wes_palette("Royal1", 20, "continuous")[c(1,7,18)]),
  theme(
    strip.text = element_text(size = rel(1.2)),
    strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank()),
  NoAxes())

# Arrows -----------------------------------------------------------------

format_add_arrows = list(
  ggplot2::theme_bw(),
  ggplot2::guides(alpha = "none"),
  ggplot2::scale_fill_gradient(low = wes_palette("Royal1", 20, "continuous")[15], high = wes_palette("Royal1", 8, "continuous")[3]), # change the low to an alternate color for nonfaceting
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
    strip.background = ggplot2::element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(), 
    panel.border = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()),
  NoAxes(),
  ggplot2::annotate(x=-16, xend=-16, y=-16, yend=-8, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate(x=-16, xend=-11, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))),
  ggplot2::annotate("text", label="UMAP 1", x=-17, -13, size = 4, angle = 90),
  ggplot2::annotate("text", label="UMAP 2", x=-14, y=-18, size = 4))

# Arrows 2 -----------------------------------------------------------------

format_add_arrows_2 = list(
  scale_fill_gradient(low = wes_palette("Royal1", 2)[1], high = wes_palette("Royal1", 20, "continuous")[18]) # change the low to an alternate color for nonfaceting
)

# Arrows 3 -----------------------------------------------------------------

format_add_arrows_3 = list(
  scale_fill_gradient(low = wes_palette("Royal1", 2)[1], high = wes_palette("Royal1", 8, "continuous")[8]) # change the low to an alternate color for nonfaceting
  
)

# Volcano -----------------------------------------------------------------

colors = c(wes_palette("Royal1", 8, "continuous")[c(1,2,8,3)], "#d9d9d9")
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

colors_2 = c(wes_palette("Royal1", 8, "continuous")[c(1,2,8)], "#d9d9d9") #c(wes_palette("Royal1", 20, "continuous")[c(1,7,18)], "#d9d9d9")
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

colors_3 = c(wes_palette("Royal1", 8, "continuous")[c(1,2,8)], "#d9d9d9") #c(wes_palette("Royal1", 20, "continuous")[c(1,7,18)], "#d9d9d9")
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

# Dotplot -----------------------------------------------------------------
dot_colors = c(wes_palette("Royal1", 8, "continuous")[c(1,2,8,3)])
names(dot_colors) = c("Inflammation", "Antigen presentation", "Cell adhesion", "Monocyte activation")

