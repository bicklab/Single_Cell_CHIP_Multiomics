#Brett Heimlich 
#
#Use CellChat to compare between mutant and wiltype at baseline

library(RColorBrewer)
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(ggalluvial)
library(cowplot)
library(reticulate)
library("wesanderson")

library(RColorBrewer)  
#colors: 
colors <- brewer.pal(12, "Paired")

#read in 
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/seurat_with_clones_5.11.23.rds")

head(seu@meta.data)
#subset for non-determined genotyped cells 
#seu <- subset(seu, subset = Clone != "NA")

#check n per celltype
table(seu@meta.data$MUTATION.GROUP, seu@meta.data$cell_type)

#subset for relatively small n within celltype
seu <- subset(seu, subset = cell_type != "Erythroid-like cells")
seu <- subset(seu, subset = cell_type != "Platelets")
seu <- subset(seu, subset = cell_type != "Dendritic cells")
seu <- subset(seu, subset = cell_type != "CD16 Monos")

seu <- subset(seu, subset = cell_type_comb != "CD8+ T cells NA control")
seu <- subset(seu, subset = cell_type_comb != "CD4+ T cells NA control")
seu <- subset(seu, subset = cell_type_comb != "NK NA control")
seu <- subset(seu, subset = cell_type_comb != "B NA control")

seu <- subset(seu, subset = cell_type_comb != "CD8+ T cells NA chip")
seu <- subset(seu, subset = cell_type_comb !=  "NK NA chip")
seu <- subset(seu, subset = cell_type_comb != "CD4+ T cells NA chip")
seu <- subset(seu, subset = cell_type_comb != "B NA chip")
seu <- subset(seu, subset = cell_type_comb != "CD14 Monos NA chip")

#add annotation for downstream labels: 
seu@meta.data$cell_type_clone <- paste(seu@meta.data$cell_type, seu@meta.data$Clone)
seu@meta.data$cell_type_comb <- paste(seu@meta.data$cell_type_clone, seu@meta.data$GENOTYPE)
table(seu@meta.data$cell_type_clone)
unique(seu@meta.data$cell_type_comb)
#split
seu.stim <- SplitObject(seu, split.by = "STIM")
seu <- seu.stim$VEH
stim <- seu.stim$STIM

seu <- NormalizeData(seu)
stim <- NormalizeData(stim)

#recheck labels and numbers of cells: 
table(seu@meta.data$MUTATION.GROUP, seu@meta.data$cell_type_clone)
table(stim@meta.data$cell_type_clone)
######
# Generate cell chat objects
#######
cellchat.veh <- createCellChat(object = seu, meta = seu@meta.data, group.by = "cell_type_clone")
cellchat.stim <- createCellChat(object = stim, meta = stim@meta.data, group.by = "cell_type_clone")

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

cellchat.veh@DB <- CellChatDB.use
cellchat.stim@DB <- CellChatDB.use
#####################
### Preprocessing ###
#####################

cellchat.stim <- subsetData(cellchat.stim)
cellchat.stim <- identifyOverExpressedGenes(cellchat.stim)
cellchat.stim <- identifyOverExpressedInteractions(cellchat.stim)
cellchat.stim <- projectData(cellchat.stim, PPI.human)

cellchat.veh <- subsetData(cellchat.veh)
cellchat.veh <- identifyOverExpressedGenes(cellchat.veh)
cellchat.veh <- identifyOverExpressedInteractions(cellchat.veh)
cellchat.veh <- projectData(cellchat.veh, PPI.human)

#################
### Inference ###
#################

cellchat.stim <- computeCommunProb(cellchat.stim)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.stim <- filterCommunication(cellchat.stim, min.cells = 10)
cellchat.stim <- computeCommunProbPathway(cellchat.stim)
cellchat.stim <- aggregateNet(cellchat.stim)
cellchat.stim <- netAnalysis_computeCentrality(cellchat.stim, slot.name = "netP")

cellchat.veh <- computeCommunProb(cellchat.veh)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.veh <- filterCommunication(cellchat.veh, min.cells = 10)
cellchat.veh <- computeCommunProbPathway(cellchat.veh)
cellchat.veh <- aggregateNet(cellchat.veh)
cellchat.veh <- netAnalysis_computeCentrality(cellchat.veh, slot.name = "netP")

cellchat.stim <- updateCellChat(cellchat.stim)
cellchat.veh <- updateCellChat(cellchat.veh)

#Save all of this stuff if you want to use it again 
saveRDS(cellchat.stim, file = "~/Brett/CH_Multiomics/seurat_objects/cellchat/cellchat_clones_andcontrols_stim.rds")
saveRDS(cellchat.veh, file = "~/Brett/CH_Multiomics/seurat_objects/cellchat/cellchat_clones_andcontrols_veh.rds")

#cellchat.stim <- readRDS(paste0(dir,"cellchat_unstim_tet.rds"))
#cellchat.dnmt <- readRDS(paste0(dir,"cellchat_unstim_dnmt.rds"))
#cellchat.control <- readRDS(paste0(dir,"cellchat_unstim_control.rds"))
cellchat <- readRDS("~/Brett/CH_Multiomics/seurat_objects/cellchat/cellchat_clones_andcontrols_stim.rds")
cellchat <- updateCellChat(cellchat)
#####################
### Visualization ###
#####################

#Cell interaction network 
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#individual cell contributions
mat <- cellchat@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Celltype and strength of incoming/outgoing singnaling
gg1 <- netAnalysis_signalingRole_scatter(cellchat, title = "Incoming vs. Outgoing Signaling")
gg1
ggsave(filename= "014mut_stim_signalingRole_scatter.pdf", plot=gg1, width = 10, height = 8, units = 'in', dpi = 300)

cells <- levels(cellchat@idents)
#heatmaps
pdf(file = "heatmaps.pdf", width = 14, height = 14)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 8, height = 8, font.size =  10, color.heatmap = "Reds")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 8, height = 8, font.size =  10, color.heatmap = "Reds")
ht1 + ht2
dev.off()
?netAnalysis_signalingRole_heatmap
netAnalysis_signalingRole_heatmap(callchatter, pattern = "outgoing", width = 10, height = 20, signaling = pathways)

netMappingDEG(cellchat, features.name = "IL1B")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat@var.features$features.info
names(cellchat@var.features)
cd14 <- subsetCellChat(cellchat, idents.use = c("CD14 Mono Mutant", "CD14 Mono Wildtype"))
levels(cellchat@idents)

#Save progress 
cellchat <- readRDS("~/analysis_out/6722/mutant_analylsis/014_stim/014_stim_cellchat.rds")
saveRDS(cellchat, file = "046_stim_cellchat.rds")

setwd("~/Brett/CH_Multiomics/Analysis_out/")
#Supervised visualizations: 
# Hierarchy plots
# define `vertex.receiver` so that the left portion of the hierarchy plot shows signaling to cell of interest and the right portion shows signaling to other cell of interest
pathways <- cellchat@netP[["pathways"]]
cells <- levels(cellchat@idents)
vertex.receiver = c(3) # a numeric vector.  find these atL: levels(cellchat@idents)
for (j in 1:length(cells)) {
  ifelse((j%%2) == 1, vertex.receiver <- c(j,j+1), next)
  for (i in 1:length(pathways)) {
    pdf(file = paste0(cells[j], "_", pathways[i],"_hierarchy_clones_stim.pdf"), width = 16, height = 10)
    gg1 <- netVisual_aggregate(cellchat, signaling = pathways[i],  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.7)
    print(gg1)
    dev.off()
  }
}
