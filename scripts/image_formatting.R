
# Image formatting --------------------------------------------------------

make_transparent = function(plots) {
  for (i in c(1:length(plots))) {
    plots[[i]] = plots[[i]] + theme(
      rect = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent', color = NA), #transparent legend bg
      legend.box.background = element_rect(fill='transparent', color = NA) #transparent legend panel
    )
  }
  return(plots)
}

save_plot = function(file, plot, legend = TRUE, width = 4, height = 4, png = FALSE, transparent = FALSE) {
  if (transparent) {
    t_plot = plot + theme(
      rect = element_rect(fill = "transparent"),
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent') #transparent legend panel
    ) + ggpubr::theme_transparent()
    ggsave(plot = t_plot, filename =  paste0(file, "_transparent.png"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in', bg = "transparent") 
  }
  if (png) { 
    ggsave(plot = plot, filename =  paste0(file, ".png"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in') 
  }
  ggsave(plot = plot, filename =  paste0(file, ".pdf"), dpi = 600, width = ifelse(legend, width, width - 1.5), height = height, units = 'in')
}

# Colors light to dark red ------------------------------------------------

format_colors = list(scale_color_gradient(low = wes_palette("Royal1", 20, "continuous")[13], high = wes_palette("Royal1", 2)[2]))

# Colors dark to light red ------------------------------------------------

format_colors_desc = list(scale_color_gradient(low = wes_palette("Royal1", 2)[2], high = wes_palette("Royal1", 20, "continuous")[13]))

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
  theme_bw(),
  guides(alpha = "none"),
  scale_fill_gradient(low = wes_palette("Royal1", 20, "continuous")[15], high = wes_palette("Royal1", 8, "continuous")[3]), # change the low to an alternate color for nonfaceting
  theme(
    strip.text = element_text(size = rel(1.2)),
    strip.background = element_rect(fill = "white", colour = "white", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank()),
  NoAxes(),
  annotate(x=-16, xend=-16, y=-16, yend=-8, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))),
  annotate(x=-16, xend=-11, y=-16, yend=-16, colour="black", lwd=1, geom="segment", arrow = arrow(length = unit(0.1, "inches"))),
  annotate("text", label="UMAP 1", x=-17, -13, size = 4, angle = 90),
  annotate("text", label="UMAP 2", x=-14, y=-18, size = 4))

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
