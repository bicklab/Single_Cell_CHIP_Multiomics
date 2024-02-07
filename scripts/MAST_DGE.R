#MAST DGE for mutant and WT comparisons

#libraries
library(ggplot2)
library(MAST)
library(Seurat)
library(SingleCellExperiment)
library(ggrepel)
library(data.table)

dir = "~/path"

#functions
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

#MAST DE custom function on subsetted populations/conditions of interest
##x = subsetted seurat object
##baseline = baseline factored level ("Wildtype")
##test_var = treatment condition for likelihood ratio test ("conditionMutant")
MAST_DE = function(x, baseline, test_var, comparison){
  
  DefaultAssay(x) = "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  
  message("Filtering variable features for genes expressed in at least 10% of cells")
  df <- PrctCellExpringGene(x, genes = VariableFeatures(x), group.by = "all")
  df_filtered <- df[df$Cell_proportion > 0.1,]
  x[['RNA']]@counts <- x[['RNA']]@counts[rownames(x[['RNA']]@counts) %in% df_filtered$Markers,]
  x[['RNA']]@data <- x[['RNA']]@data[rownames(x[['RNA']]@data) %in% df_filtered$Markers,]
  
  sca = as.SingleCellExperiment(x, assay = "RNA")
  
  #Calculate cellular detection rate
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  #factor levels, define contrasts, run linear model and lrt 
  cond<-factor(colData(sca)$Clone)
  cond<-relevel(cond, baseline)
  colData(sca)$condition<-cond
  sca_1 = SceToSingleCellAssay(sca, check_sanity = F)
  zlmCond <- zlm(~condition + cngeneson, sca_1)
  summaryCond <- summary(zlmCond, doLRT= test_var) 
  
  #organize results
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast==test_var & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==test_var & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>0.25], as.data.table(mcols(sca_1)), by='primerid')
  
  fcHurdle$label = fcHurdle$primerid %in% fcHurdleSig$primerid
  write.csv(fcHurdle, paste0('pooled_mutant_vs_wildtype/', comparison, ".csv"))

}



#Load seurat object
seu <- readRDS('seurat_with_clones.rds')


dim(seu)
counts <- seu[['RNA']]@counts
metadata <- seu@meta.data
seu <- CreateSeuratObject(counts, min.cells = 10, meta.data = metadata)
dim(seu)

#subset for samples of interest
DNMT3A <- c("CH-21-006", "CH-21-046")
TET2 <- c("CH-20-004", "CH-21-014", "CH-21-033", "CH-21-046", "CH-21-073", "CH-21-074")
ctrl <- c("CH-21-002", "CH-21-008", "CH-21-013", "CH-21-020", "CH-21-021", "CH-21-028", 'CH-21-031')
combined <-c(DNMT3A, TET2, ctrl)

Idents(seu) <- "CHIVEID"
seu <- subset(seu, idents = combined)

#ID healthy cells
Idents(seu) <- "Clone"
seu$Clone[is.na(seu$Clone)] <- "healthy"
table(seu$CHIVEID, seu$Clone)

#Subset all relevant cell types (done individually for reproducibility)
Idents(seu) <- "MUTATION.GROUP"
DNMT3A <- subset(seu, idents = "TET2", invert = T)
TET2 <- subset(seu, idents = "DNMT3A", invert = T)

Idents(DNMT3A) <- "cell_type"
DNMT3A.CD14 <- subset(DNMT3A, idents = "CD14 Monos")
DNMT3A.CD8 <- subset(DNMT3A, idents = "CD8+ T cells")
DNMT3A.CD4 <- subset(DNMT3A, idents = "CD4+ T cells")
DNMT3A.NK <- subset(DNMT3A, idents = "NK")
DNMT3A.B <- subset(DNMT3A, idents = "B")

Idents(TET2) <- "cell_type"
TET2.CD14 <-  subset(TET2, idents = "CD14 Monos")
TET2.CD8 <- subset(TET2, idents = "CD8+ T cells")
TET2.CD4 <- subset(TET2, idents = "CD4+ T cells")
TET2.NK <- subset(TET2, idents = "NK")
TET2.B <- subset(TET2, idents = "B")


MAST_DE(DNMT3A.CD14, baseline = "Wildtype", test_var = "conditionMutant", comparison = "DNMT3A CD14 Mutant vs Wildtype")
MAST_DE(DNMT3A.CD8, baseline = "Wildtype", test_var = "conditionMutant", comparison = "DNMT3A CD8 Mutant vs Wildtype")
MAST_DE(DNMT3A.CD4, baseline = "Wildtype", test_var = "conditionMutant", comparison = "DNMT3A CD4 Mutant vs Wildtype")
MAST_DE(DNMT3A.NK, baseline = "Wildtype", test_var = "conditionMutant", comparison = "DNMT3A NK Mutant vs Wildtype")
MAST_DE(DNMT3A.B, baseline = "Wildtype", test_var = "conditionMutant", comparison = "DNMT3A B Mutant vs Wildtype")

MAST_DE(TET2.CD14, baseline = "Wildtype", test_var = "conditionMutant", comparison = "TET2 CD14 Mutant vs Wildtype")
MAST_DE(TET2.CD8, baseline = "Wildtype", test_var = "conditionMutant", comparison = "TET2 CD8 Mutant vs Wildtype")
MAST_DE(TET2.CD4, baseline = "Wildtype", test_var = "conditionMutant", comparison = "TET2 CD4 Mutant vs Wildtype")
MAST_DE(TET2.NK, baseline = "Wildtype", test_var = "conditionMutant", comparison = "TET2 NK Mutant vs Wildtype")
MAST_DE(TET2.B, baseline = "Wildtype", test_var = "conditionMutant", comparison = "TET2 B Mutant vs Wildtype")


