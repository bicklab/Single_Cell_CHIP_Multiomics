#Brett Heimlich 

#Goal is to Demultiplex and ID Samples and genoytes for downstream analysis. Product is labeled, unmodified Seurat object

#~~~~~~~~~~~~~~~~~~~~~~#
#### Pre-requisites ####
#~~~~~~~~~~~~~~~~~~~~~~#

options(stringsAsFactors = FALSE)
options(scipen = 999)
set.seed(100)

library(readxl)
library(Matrix)
library(Seurat)
library(tidyverse)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Read in and create seurat obj####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Read in raw data
data_dir <- '~/Brett/maester/7079/cellranger_output/Sample 1/'
list.files(data_dir)
counts1 <- Read10X(data.dir = data_dir)
umis1 <- counts1$"Gene Expression"
htos1 <- counts1$"Antibody Capture"
#adt1 <- counts1$"Antibody Capture"

#Select first 8 lines (HTO info) (or 4 below)
#htos1 <- htos1[1:8,]
#select next 17 lines for cell surface protein markers
#adt1 <- adt1[9:25,]

htos1 <- htos1[1:4,]
#adt1 <- adt1[7:146,]

joint.bcs <- intersect(colnames(umis1), colnames(htos1))
umis1 <- umis1[, joint.bcs]
htos1 <- as.matrix(htos1[, joint.bcs])
#adt1 <- as.matrix(adt1[, joint.bcs])
#rownames(adt1)
rownames(htos1)

#Create Seurat object and add protein and HTO Assays
seu <- CreateSeuratObject(counts = umis1)
seu[["HTO"]] <- CreateAssayObject(counts = htos1)
#seu[["ADT"]] <- CreateAssayObject(counts = adt1)

#Validate assays and assign default
seu
Assays(seu)
DefaultAssay(seu)

#change wd for where files should be saved
setwd("~/Brett/seurat_objects/")

#Export raw seurat object for later use if needed // Import if starting over 
#saveRDS(seu, "CH_21_016_raw_seurat.rds")
#seu <- readRDS("~/maester/7079/7079_analysis_output/7079_sample4_raw_seurat.rds")

#~~~~~~~~~~~~~~~~~~~#
#### Demultiplex ####
#~~~~~~~~~~~~~~~~~~~#

seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.95)
RidgePlot(seu, assay = "HTO", features = rownames(seu[["HTO"]])[1:4], ncol = 2)

#visualize heatmap 
#HTOHeatmap(seu, assay = "HTO", ncells = 5000)

# Visualize and Extract the singlets
table(seu$HTO_classification.global)

#change idents to find sing/doub/negs
Idents(seu) <- "HTO_classification.global"
VlnPlot(seu, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#Keep only singlets
seu.singlet <- subset(seu, idents = "Singlet")

#return Idents to originals 
#Idents(seu) <- "CellType"
#Idents(seu.singlet) <- "orig.ident"

#assign group names based on experimental conditions and genotype status
#################
#for 7527-2: 
#################
#Add chip status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("chip", "chip", "chip", "chip", "chip", "control")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVE ID 
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("CH-21-029", "CH-21-073", "CH-21-074", "CH-21-077", "CH-21-079", "CH-21-028")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#Add Project
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("7527", "7527", "7527", "7527", "7527", "7527")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#Add Lane
values <- c("7527-2", "7527-2", "7527-2", "7527-2", "7527-2", "7527-2")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#for 7527-1: 
#################
#Add chip status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("chip", "chip", "chip", "chip", "chip", "control")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVE ID 
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("CH-20-001", "CH-20-002", "CH-20-005", "CH-21-036", "CH-21-017", "CH-21-002")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#Add Project
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6")
values <- c("7527", "7527", "7527", "7527", "7527", "7527")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#Add Lane
values <- c("7527-1", "7527-1", "7527-1", "7527-1", "7527-1", "7527-1")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#For 7079-1: **Change this for each lane** 
#################
#Add stim status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4")
values <- c("chipstim", "chipstim", "chip", "chip")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVEID
values <- c("CH-20-004", "CH-21-014", "CH-20-004", "CH-21-014")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add chip status 
values <- c("chip", "chip", "chip", "chip")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("7079-1", "7079-1", "7079-1", "7079-1")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("7079", "7079", "7079", "7079")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#For 7079-2: 
#################

#Add stim status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4")
values <- c("chipstim", "chipstim", "chip", "chip")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVEID
values <- c("CH-21-046", "CH-21-034", "CH-21-046", "CH-21-034")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add chip status 
values <- c("chip", "chip", "chip", "chip")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("7079-2", "7079-2", "7079-2", "7079-2")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("7079", "7079", "7079", "7079")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#For 7079-3: 
#################

#Add stim status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4")
values <- c("chipstim", "chipstim", "chip", "chip")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVEID
values <- c("CH-21-006", "CH-21-033", "CH-21-006", "CH-21-033")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add chip status 
values <- c("chip", "chip", "chip", "chip")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("7079-3", "7079-3", "7079-3", "7079-3")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("7079", "7079", "7079", "7079")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#For 7079-4: 
#################

#Add stim status 
index <- c("sample-1", "sample-2", "sample-3", "sample-4")
values <- c("controlstim", "controlstim", "control", "control")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVEID
values <- c("CH-21-008", "CH-21-013", "CH-21-008", "CH-21-013")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add chip status 
values <- c("control", "control", "control", "control")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("7079-4", "7079-4", "7079-4", "7079-4")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("7079", "7079", "7079", "7079")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#This is for 6903:
#################
index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6", "sample-7", "sample-8")
values <- c("control", "control", "controlstim", "controlstim", "chip", "chip", "chipstim", "chipstim")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

index <- c("sample-1", "sample-2", "sample-3", "sample-4", "sample-5", "sample-6", "sample-7", "sample-8")
values <- c("CH-21-031", "CH-21-031", "CH-21-031", "CH-21-031", "CH-21-037", "CH-21-037", "CH-21-037", "CH-21-037")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

values <- c("control", "control", "control", "control", "chip", "chip", "chip", "chip")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("6903-1", "6903-1", "6903-1", "6903-1", "6903-1", "6903-1", "6903-1", "6903-1")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("6903", "6903", "6903", "6903", "6903", "6903", "6903", "6903")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#################
#This is for 7373-2:
#################

#Add stim status 
index <- c("sample-5", "sample-6", "sample-7", "sample-8")
values <- c("controlstim", "controlstim", "control", "control")
seu.singlet@meta.data$TYPE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add CHIVEID
values <- c("CH-21-021", "CH-21-020", "CH-21-021", "CH-21-020")
seu.singlet@meta.data$CHIVEID <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add chip status 
values <- c("control", "control", "control", "control")
seu.singlet@meta.data$CHIP <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add lane
values <- c("7373-2", "7373-2", "7373-2", "7373-2")
seu.singlet@meta.data$LANE <- values[match(seu.singlet@meta.data$hash.ID, index)]

#add project
values <- c("7373", "7373", "7373", "7373")
seu.singlet@meta.data$ProjectID <- values[match(seu.singlet@meta.data$hash.ID, index)]


#Check meta data
colnames(seu.singlet@meta.data)
seu.singlet@meta.data
VlnPlot(seu.singlet, features = "IL1B", group.by = "CHIVEID")

#Save seurat with meta-data IDs and groupings (if desired) 
#saveRDS(seu, "CH_21_016_ID_seurat.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Create Sample Specific Seurats ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

seu.singlet.list <- SplitObject(seu.singlet, split.by = "CHIVEID")
#Create individual seurat objects for each sample

#7079-1
CH_20_004 <- seu.singlet.list$`CH-20-004`
CH_21_014 <- seu.singlet.list$`CH-21-014`

saveRDS(CH_20_004, "CH_20_004_rawseurat.rds")
saveRDS(CH_21_014, "CH_21_014_rawseurat.rds")

#7079-2
CH_21_034 <- seu.singlet.list$`CH-21-034`
CH_21_046 <- seu.singlet.list$`CH-21-046`

saveRDS(CH_21_034, "CH_21_034_rawseurat.rds")
saveRDS(CH_21_046, "CH_21_046_rawseurat.rds")

#7079-3
CH_21_033 <- seu.singlet.list$`CH-21-033`
CH_21_006 <- seu.singlet.list$`CH-21-006`

saveRDS(CH_21_033, "CH_21_033_rawseurat.rds")
saveRDS(CH_21_006, "CH_21_006_rawseurat.rds")

#7079-4
CH_21_013 <- seu.singlet.list$`CH-21-013`
CH_21_008 <- seu.singlet.list$`CH-21-008`

saveRDS(CH_21_013, "CH_21_013_rawseurat.rds")
saveRDS(CH_21_008, "CH_21_008_rawseurat.rds")

CH_21_037 <- seu.singlet.list$`CH-21-037`
CH_21_031 <- seu.singlet.list$`CH-21-031`

saveRDS(CH_21_031, "CH_21_031_rawseurat.rds")
saveRDS(CH_21_037, "CH_21_037_rawseurat.rds")

#7373
CH_21_020 <- seu.singlet.list$`CH-21-020`
CH_21_021 <- seu.singlet.list$`CH-21-021`

saveRDS(CH_21_020, "CH_21_020_rawseurat.rds")
saveRDS(CH_21_021, "CH_21_021_rawseurat.rds")

#7527
CH_21_029 <- seu.singlet.list$`CH-21-029`
CH_21_073 <- seu.singlet.list$`CH-21-073`
CH_21_074 <- seu.singlet.list$`CH-21-074`
CH_21_077 <- seu.singlet.list$`CH-21-077`
CH_21_079 <- seu.singlet.list$`CH-21-079`
CH_21_028 <- seu.singlet.list$`CH-21-028`

CH_20_001 <- seu.singlet.list$`CH-20-001`
CH_20_002 <- seu.singlet.list$`CH-20-002`
CH_20_005 <- seu.singlet.list$`CH-20-005`
CH_21_036 <- seu.singlet.list$`CH-21-036`
CH_21_017 <- seu.singlet.list$`CH-21-017`
CH_21_002 <- seu.singlet.list$`CH-21-002`

#Save Files 

saveRDS(CH_21_029, "CH_21_029_rawseurat.rds")
saveRDS(CH_21_073, "CH_21_073_rawseurat.rds")
saveRDS(CH_21_074, "CH_21_074_rawseurat.rds")
saveRDS(CH_21_077, "CH_21_077_rawseurat.rds")
saveRDS(CH_21_079, "CH_21_079_rawseurat.rds")
saveRDS(CH_21_028, "CH_21_028_rawseurat.rds")

CH_21_029 <- seu.singlet.list$`CH-21-029`
CH_21_073 <- seu.singlet.list$`CH-21-073`
CH_21_074 <- seu.singlet.list$`CH-21-074`
CH_21_077 <- seu.singlet.list$`CH-21-077`
CH_21_079 <- seu.singlet.list$`CH-21-079`
CH_21_028 <- seu.singlet.list$`CH-21-028`

saveRDS(CH_20_001, "CH_20_001_rawseurat.rds")
saveRDS(CH_20_002, "CH_20_002_rawseurat.rds")
saveRDS(CH_20_005, "CH_20_005_rawseurat.rds")
saveRDS(CH_21_036, "CH_21_036_rawseurat.rds")
saveRDS(CH_21_017, "CH_21_017_rawseurat.rds")
saveRDS(CH_21_002, "CH_21_002_rawseurat.rds")

setwd("~/maester/7079/7079_analysis_output/sample_seurats")
saveRDS(CH_20_004, "CH_20_004_seurat.rds")
saveRDS(CH_21_014, "CH_21_014_seurat.rds")
saveRDS(CH_21_034, "CH_21_034_seurat.rds")
saveRDS(CH_21_046, "CH_21_046_seurat.rds")
saveRDS(CH_21_033, "CH_21_033_seurat.rds")
saveRDS(CH_21_006, "CH_21_006_seurat.rds")
saveRDS(CH_21_013, "CH_21_013_seurat.rds")
saveRDS(CH_21_008, "CH_21_008_seurat.rds")

saveRDS(CH_21_037, "CH_21_037_seurat.rds")
saveRDS(CH_21_031, "CH_21_031_seurat.rds")
#Read Files
