install.packages("tidyverse")

library(Seurat)
library(RColorBrewer)  
library(tidyverse)
#colors: 
colors <- brewer.pal(12, "Paired")

seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/seurat_with_clones_PB_10.16.23.rds")

obj.list <- SplitObject(seu, split.by = "STIM")
seu_veh <- obj.list$VEH
seu_veh <- NormalizeData(seu_veh)
seu_stim <- obj.list$STIM
seu_stim <- NormalizeData(seu_stim)

obj.list2 <- SplitObject(seu_veh, split.by = "GROUP")
seu_tet_veh <- obj.list2$TET2
seu_dnm_veh <- obj.list2$DNMT3A

Idents(seu_dnm_veh) <- "Clone"
Idents(seu_tet_veh) <- "Clone"

Idents(seu_veh) <- "MUTATION.GROUP"
Idents(seu_stim) <- "cell_type"

seu_veh@meta.data$MUTATION.GROUP[seu_veh@meta.data$MUTATION.GROUP == "none"] = "Control"
colnames(seu_veh@meta.data)
seu_tet_veh@meta.data$Clone[seu_tet_veh@meta.data$Clone == "Mutant"] = "TET2 Mutant"
seu_dnm_veh@meta.data$Clone[seu_dnm_veh@meta.data$Clone == "Mutant"] = "DNMT3A Mutant"

gg1 <- VlnPlot(subset(seu_veh, subset = cell_type == c("CD14 Monos")), 
        idents = c("TET2", "Control"),
        features = c("MIF", "CD44", "CD74"), 
        #group.by = "MUTATION.GROUP", 
        slot = "data", 
        add.noise = F,
        adjust = 2, 
        #cols = colors[c(6,2,4)],
        assay = "RNA") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank())
gg1
ggsave("~/Brett/CH_Multiomics/Analysis_out/violin_plots/FOSB_CD4_DNM.png", gg1, device = "png", dpi = 600)

#Get Average Expression Data
seu_veh@meta.data$combined <- paste0(seu_veh@meta.data$MUTATION.GROUP, "_", seu_veh@meta.data$Clone)
Idents(seu_veh) <- "combined"
AverageExpression(object = subset(seu_veh, subset = cell_type == c("CD14 Monos")),assays = "RNA", features = c("IL1B", "CXCL3", "CXCL1", "CCL4", "CCL2", "CCL7"))
