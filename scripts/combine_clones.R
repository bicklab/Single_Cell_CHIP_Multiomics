#Combine clone information by adding metadata to our final object. 

#Brett Heimlich 

library(Seurat)

#import main object 
seu <- readRDS("~/Brett/CH_Multiomics/seurat_objects/final_annotated_seurat_correctedDR_5.11.23.rds")
seu.list <- SplitObject(seu, split.by = "CHIVEID")

#import the af.dm (from all individual samples), make a smaller dataframe including the clone call and the cell barcode
files <- list.files(path = "~/Brett/CH_Multiomics/Maester_output/", pattern = "afmatrix_*", full.names = T, recursive = T)
#take this DF and then merge/add metadata to the main object. Repeat for all objects

varlist <- read_csv("~/Brett/CH_Multiomics/seurat_objects/Varlist.csv", col_names = T)
chivelist <- unique(varlist$CHIVEID)

for(i in 1:length(chivelist)){
  #i=1
  CHIVEID <- chivelist[i]
  lastthree <- substring(CHIVEID, 7,9) #trying to get the the last 3 of the CHIVEID
  if(lastthree == "002") {lastthree <- "0002"}
  x <- grep(CHIVEID, files)
  af.dm <- read.table(files[x]) 

  colnames(af.dm) <- paste0(lastthree, "_", colnames(af.dm))
  colnames(af.dm) <- gsub("\\.", "-", colnames(af.dm))
  
  mtvars <- varlist$mtvar[varlist$CHIVEID == CHIVEID]
  
  af_voi <- af.dm[mtvars,]
  af_voi <- t(af_voi)
  
  max <- rowMax(af_voi)
  names(max) <- rownames(af_voi)
  af_voi <- as.data.frame(af_voi)
  af_voi$max <- max 
  
  af_voi$Clone <- ifelse(af_voi$max > 1, "Mutant", "Wildtype")
  
  seu.list[[CHIVEID]] <- AddMetaData(seu.list[[CHIVEID]], select(af_voi, Clone))
  x <- NULL 
}

unstim_cd14monomarkers <- FindMarkers(seu_tet_veh, 
                                      group.by = "Clone",
                                      ident.1 = "Mutant",
                                      subset.ident = c("CD14 Monos", "CD16 Monos"),
                                      min.pct = 0.01,
                                      min.diff.pct = 0.01,
                                      assay = "RNA")

seu.list[[CHIVEID]]@meta.data

seu_clones <- merge(seu.list$`CH-20-002`, y = c(seu.list$`CH-20-001`, seu.list$`CH-21-002`, seu.list$`CH-20-004`,
                                                seu.list$`CH-20-005`, seu.list$`CH-21-006`, seu.list$`CH-21-008`, 
                                                seu.list$`CH-21-013`, seu.list$`CH-21-014`, seu.list$`CH-21-017`, 
                                                seu.list$`CH-21-020`, seu.list$`CH-21-021`, seu.list$`CH-21-028`, 
                                                seu.list$`CH-21-029`, seu.list$`CH-21-031`, seu.list$`CH-21-033`,
                                                seu.list$`CH-21-034`, seu.list$`CH-21-036`, seu.list$`CH-21-037`, 
                                                seu.list$`CH-21-046`, seu.list$`CH-21-073`, seu.list$`CH-21-074`, 
                                                seu.list$`CH-21-077`, seu.list$`CH-21-079`), project = "MERGED", merge.data = TRUE)
ch014 <- seu.singlet
ch014 <- seu.list$`CH-21-014`
split.ch014 <- SplitObject(ch014, split.by = "TYPE")
ch014veh <- split.ch014$chip

ch014veh <- NormalizeData(ch014veh)
ch014veh <- FindVariableFeatures(ch014veh)
ch014veh <- ScaleData(ch014veh)
Idents(ch014veh) <- "clonecells"

unstim_cd14monomarkers <- FindMarkers(ch014veh, 
                                      group.by = "Clone",
                                      ident.1 = "Mutant",
                                      #ident.2 = "Wildtype",
                                      subset.ident = c("CD14 Monos"),
                                      min.pct = 0.01,
                                      min.diff.pct = 0.01,
                                      assay = "RNA", 
                                      features = ch014veh@assays[["RNA"]]@var.features) 
                                      #slot = "counts")

saveRDS(seu_clones, file = "~/Brett/CH_Multiomics/seurat_objects/seurat_with_clones_5.18.23_AF2.rds")

clonecells <- ch014@meta.data$cell_type
clonecells <- as.data.frame(clonecells)
row.names(clonecells) <- row.names(ch014@meta.data)
clonecells$Clone <- ch014@meta.data$Clone

write_csv(clonecells, file = "~/Brett/CH_Multiomics/seurat_objects/ch014clonecells.csv", col_names = T)
