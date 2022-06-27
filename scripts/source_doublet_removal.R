
# Doublet removal -----------------------------------------------------------

source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_all_patients.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/source_individual_patients.R")
source("/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/scripts/image_formatting.R")

# Pilot strict quality control and doublet removal -----------------------
nrow(all_data@meta.data) # 109428
all_data = subset(all_data, subset = nFeature_RNA < 2500 & nCount_RNA < 7500)
nrow(all_data@meta.data) # 108569

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 2000)
all_data <- ScaleData(all_data)
all_data <- RunPCA(all_data)
all_data <- RunUMAP(all_data, dims = 1:10)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- all_data@meta.data$predicted.celltype.l2
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- all_data@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(all_data@meta.data))  ## Assuming 5% doublet formation rate - based on 10X standards
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
all_data <- doubletFinder_v3(all_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
all_data <- doubletFinder_v3(all_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_5428", sct = FALSE)

VlnPlot(all_data, features = c("DF.classifications_0.25_0.09_8143", "celltype.de"), group.by = "GENOTYPE")

barplot(table(all_data@meta.data$MUTATION.GROUP, all_data@meta.data$DF.classifications_0.25_0.09_4434))
barplot(prop.table(all_data@meta.data$DF.classifications_0.25_0.09_4434, all_data@meta.data$STIM))

doublet = all_data@meta.data
write_tsv(doublet %>% rownames_to_column(), file = "input_data/doublet_metadata.tsv") 

saveRDS(all_data, "rds_objects/doublet_labeled_data.rds")

Idents(all_data) = "DF.classifications_0.25_0.09_4434"
singlet_data = subset(all_data, idents = "Singlet")

nrow(singlet_data@meta.data)

saveRDS(singlet_data, "rds_objects/singlet_data.rds")


