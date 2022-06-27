
# Image formatting --------------------------------------------------------

save_plot = function(file, plot, legend = TRUE, width = 7, height = 7) {
  ggsave(plot = plot, filename =  paste0(file, ".png"), dpi = 600, width = ifelse(legend, width, width - 3), height = height, units = 'in')
  ggsave(plot = plot, filename =  paste0(file, ".pdf"), dpi = 600, width = ifelse(legend, width, width - 3), height = height, units = 'in')
}

# Colors light to dark red ------------------------------------------------

format_colors = list(scale_color_gradient(low = wes_palette("Royal1", 20, "continuous")[13], high = wes_palette("Royal1", 2)[2]))

# Colors dark to light red ------------------------------------------------

format_colors_desc = list(scale_color_gradient(low = wes_palette("Royal1", 2)[2], high = wes_palette("Royal1", 20, "continuous")[13]))

# No axes -----------------------------------------------------------------

format_no_axes = list(scale_color_gradient(low = "#f0f0f0", high = wes_palette("Royal1", 2)[2]),
                      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                            axis.title = element_blank()))

# No axes - faceted --------------------------------------------------------

format_no_axes_faceted = list(
  theme_bw(),
  guides(color = "none"),
  scale_color_manual(values = wes_palette("Royal1", 20, "continuous")[c(1,7,18)]),
  theme(
    strip.text = element_text(size = rel(1.2), face = "bold"),
    strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank()),
  NoAxes())

# Arrows -----------------------------------------------------------------

format_add_arrows = list(
  theme_bw(),
  guides(alpha = "none"),
  scale_fill_gradient(low = wes_palette("Royal1", 20, "continuous")[15], high = wes_palette("Royal1", 8, "continuous")[3]), # change the low to an alternate color for nonfaceting
  theme(
    strip.text = element_text(size = rel(1.5), face = "bold"),
    strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank()),
  NoAxes(),
  annotate(x=-16, xend=-16, y=-16, yend=-10, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))),
  annotate(x=-16, xend=-12, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))),
  annotate("text", label="UMAP 1", x=-17.5, -14, size = 4, angle = 90),
  annotate("text", label="UMAP 2", x=-15, y=-17, size = 4))

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
names(colors) = c("Inflammation", "Antigen Presentation", "Cell Adhesion", "Monocyte Activation", "Other")
format_volcano = list(
  xlab("log2(fold change)"),
  ylab("-log10(p-value)"),
  labs(color = "Pathway"),
  theme_bw(),
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()),
  scale_color_manual(values = colors)
)

# Volcano 2 -----------------------------------------------------------------

colors_2 = c(wes_palette("Royal1", 8, "continuous")[c(1,2,8)], "#d9d9d9") #c(wes_palette("Royal1", 20, "continuous")[c(1,7,18)], "#d9d9d9")
names(colors_2) = c("T Cell Activation", "Intracellular Signaling", "Immunomodulation", "Other")
format_volcano_2 = list(
  xlab("log2(fold change)"),
  ylab("-log10(p-value)"),
  labs(color = "Pathway"),
  theme_bw(),
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()),
  scale_color_manual(values = colors_2)
)
