#Brett Heimlich 
library(dplyr)
library(Seurat)
#Goal is to take individual seurat objects and combine them into a single object for downstream applications. 
# From there, will be filtered, normalized, and re-integrated. 

#Import required files 

#7079-1
CH_20_004 <- readRDS("~/Brett/seurat_objects/metadata/CH_20_004_seurat_metadata.rds")
CH_21_014 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_014_seurat_metadata.rds")

#7079-2
CH_21_034 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_034_seurat_metadata.rds")
CH_21_046 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_046_seurat_metadata.rds")

#7079-3
CH_21_033 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_033_seurat_metadata.rds")
CH_21_006 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_006_seurat_metadata.rds")

#7079-4
CH_21_013 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_013_seurat_metadata.rds")
CH_21_008 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_008_seurat_metadata.rds")

#6903
CH_21_037 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_037_seurat_metadata.rds")
CH_21_031 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_031_seurat_metadata.rds")

#7373
CH_21_020 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_020_seurat_metadata.rds")
CH_21_021 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_021_seurat_metadata.rds")

#7527
CH_21_029 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_029_seurat_metadata.rds")
CH_21_073 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_073_seurat_metadata.rds")
CH_21_074 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_074_seurat_metadata.rds")
CH_21_077 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_077_seurat_metadata.rds")
CH_21_079 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_079_seurat_metadata.rds")
CH_21_028 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_028_seurat_metadata.rds")
CH_20_001 <- readRDS("~/Brett/seurat_objects/metadata/CH_20_001_seurat_metadata.rds")
CH_20_002 <- readRDS("~/Brett/seurat_objects/metadata/CH_20_002_seurat_metadata.rds")
CH_20_005 <- readRDS("~/Brett/seurat_objects/metadata/CH_20_005_seurat_metadata.rds")
CH_21_036 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_036_seurat_metadata.rds")
CH_21_017 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_017_seurat_metadata.rds")
CH_21_002 <- readRDS("~/Brett/seurat_objects/metadata/CH_21_002_seurat_metadata.rds")


#install doublet finder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

#run doublet finder on each object 

#combine object 
#plan for downstream harmony integration. 

filtering <- function(x){
  subset(x = x, subset = (nCount_RNA <= 10000) & 
           (nFeature_RNA >= 200) & 
           (percent.mt < 20)) }

lapply(1:ncol(seu@meta.data), function(i) unique(seu@meta.data[[i]]))



#Combine 
#Join into one large seurat 
seu <- merge(CH_20_004, y = c(CH_21_014, CH_21_034, CH_21_046, CH_21_033, CH_21_006, CH_21_008, 
                                  CH_21_013, CH_21_037, CH_21_031, CH_21_020, CH_21_021, CH_21_029,
                              CH_21_073, CH_21_074, CH_21_077, CH_21_079, CH_21_028, CH_20_001,
                              CH_20_002, CH_20_005, CH_21_036, CH_21_017, CH_21_002), 
                 add.cell.ids = c("004", "014", "034", "046", "033", "006", 
                                  "008", "013", "037", "031", "020", "021", 
                                  "029", "073", "074", "077", "079", "028", 
                                  "001", "0002", "005", "036", "017", "002"), project = "MERGED")

#Check metadata
head(seu@meta.data)
lapply(1:ncol(seu@meta.data), function(i) unique(seu@meta.data[[i]]))
#b <- filter(seu@meta.data, is.na(seu@meta.data$STIM))

#save
setwd("~/Brett/seurat_objects/")
saveRDS(seu, "combined_rawseurat.rds")
seu <- readRDS("~/Brett/seurat_objects/raw/combined_rawseurat.rds")
#Filter
if (any(grepl(pattern = '^MT-', x = rownames(x = seu)))) {
  seu<- PercentageFeatureSet(
    object = seu,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

Idents(seu) <- "ProjectID"
VlnPlot(seu1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

filtering <- function(x){
  subset(x = x, subset = (nCount_RNA <= 10000) & 
           (nFeature_RNA >= 200) & 
           (percent.mt < 20)) }

seu1 <- filtering(seu)
nrow(seu1@meta.data)

saveRDS(seu1, "~/Brett/seurat_objects/filtered_seurat.rds")
