#Brett Heimlich 
#Run Doublet finder on raw data (with metadata) 

#install doublet finder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(Seurat)

#import merged object 
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/raw/combined_rawseurat.rds")

#Filter
if (any(grepl(pattern = '^MT-', x = rownames(x = seu)))) {
  seu<- PercentageFeatureSet(
    object = seu,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

filtering <- function(x){
  subset(x = x, subset = (nCount_RNA <= 10000) & 
           (nFeature_RNA >= 200) & 
           (percent.mt < 20)) }

seu1 <- filtering(seu)

#separate and create list by lane
seu.list <- SplitObject(seu1, split.by = "LANE")
#perform loop/lapply to ID doublets

#make function to perform desired tasks related to doublet finder 
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:10)
})

seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- FindNeighbors(x, reduction = "pca", dims = 1:10)
  x <- FindClusters(x)
})

seu.list <- lapply(X = seu.list, FUN = function(x) {
  sweep.res.list <- paramSweep_v3(x, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  homotypic.prop <- modelHomotypic(x@meta.data$seurat_clusters)           ## ex: annotations <- x@meta.data$seurat_clusters
  nExp_poi <- round(0.075*nrow(x@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <- doubletFinder_v3(x, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
}) 

#put back into single object for downstream analysis
saveRDS(seu.list, "~/Brett/seurat_objects/seurat_list_afterDF.rds")
seu.list <- readRDS("~/Brett/CH_Multiomics/seurat_objects/old/Previous-versions/seurat_list_afterDF.rds")

#select only singlets
seu.list$`7079-1`$DF <- seu.list$`7079-1`$DF.classifications_0.25_0.09_1295
seu.list$`7079-2`$DF <- seu.list$`7079-2`$DF.classifications_0.25_0.09_1144
seu.list$`7079-3`$DF <- seu.list$`7079-3`$DF.classifications_0.25_0.09_1198
seu.list$`7079-4`$DF <- seu.list$`7079-4`$DF.classifications_0.25_0.09_929
seu.list$`6903-1`$DF <- seu.list$`6903-1`$DF.classifications_0.25_0.09_1289
seu.list$`7373-2`$DF <- seu.list$`7373-2`$DF.classifications_0.25_0.09_1173
seu.list$`7527-2`$DF <- seu.list$`7527-2`$DF.classifications_0.25_0.09_926
seu.list$`7527-1`$DF <- seu.list$`7527-1`$DF.classifications_0.25_0.09_1014

#subset for known missed doublets
idents(seu.list$`7079-3` )
head(seu.list$`7079-3`@meta.data)

Idents(seu.list$`7079-3`) <- "pANN_0.25_0.09_1198"

seu.list$`7079-3` <- subset(seu.list$`7079-3`, subset = pANN_0.25_0.09_1198 < 0.28)
#combine back into one object
seu <- merge(seu.list$`7079-1`, y = c(seu.list$`7079-2`, seu.list$`7079-3`, seu.list$`7079-4`, 
                                      seu.list$`6903-1`, seu.list$`7373-2`, seu.list$`7527-2`, 
                                      seu.list$`7527-1`), project = "MERGED", merge.data = TRUE)

head(seu@meta.data)

prop.table(table(seu@meta.data$DF, seu@meta.data$GROUP))

#save combined object before removal of doublets and integration
saveRDS(seu, "~/Brett/CH_Multiomics/seurat_objects/seurat_afterDF_refined5.10.23.rds")
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/seurat_afterDF_refined5.10.23.rds")

head(seu@meta.data)
table(seu@meta.data$DF)
#remove doublets
Idents(seu) <- "DF"
VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#Keep only singlets
seu <- subset(seu, idents = "Singlet")

#required pre-analysis for harmony 
seu <- Seurat::NormalizeData(seu, verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = T) %>% 
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = T)
seu <- RunUMAP(seu, dims = 1:20)

#Run Harmony: 
library(harmony)
seu.harmony <- seu %>% 
  RunHarmony(group.by.vars = c("STIM", "GROUP", "LANE"), plot_convergence = TRUE)

saveRDS(seu.harmony, "~/Brett/CH_Multiomics/seurat_objects/seurat_DF_Harmony.rds")
seu.harmony <- readRDS("~/Brett/seurat_objects/seurat_DF_Harmony_KEEP.rds")
#some visualizations

seu.harmony <- seu.harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5)

s1 <- DimPlot(seu.harmony, reduction = "harmony", label = F, pt.size = .1, group.by = "LANE", shuffle = T)
s2 <- DimPlot(seu.harmony, reduction = "harmony", label = F, pt.size = .1, group.by = "STIM", shuffle = T)
s3 <- DimPlot(seu.harmony, reduction = "harmony", label = F, pt.size = .1, group.by = "GROUP", shuffle = T)
s4 <- DimPlot(seu.harmony, reduction = "harmony", label = F, pt.size = .1, group.by = "GENOTYPE", shuffle = T)

s1 + s2 + s3 + s4

t1 <- DimPlot(seu, reduction = "umap", label = F, pt.size = .1, group.by = "LANE", shuffle = T)
t2 <- DimPlot(seu, reduction = "umap", label = F, pt.size = .1, group.by = "STIM", shuffle = T)
t3 <- DimPlot(seu, reduction = "umap", label = F, pt.size = .1, group.by = "GROUP", shuffle = T)
t4 <- DimPlot(seu, reduction = "umap", label = F, pt.size = .1, group.by = "GENOTYPE", shuffle = T)

t1 + t2 + t3 + t4
