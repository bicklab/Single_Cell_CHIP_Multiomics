#MFI_VAF correlation, venn diagram of overlapping genes, and genomic map of mitochondrial mutations

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   MFI_VAF correlation     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(ggplot2)
library(ggpubr)


VAF <- c(0.3,
         0.51,
         0,
         0,
         0,
         0,
         0,
         0.02,
         0.23,
         0.42
)


Median <- c(238.4431069,
            203.6374956,
            183.5920714,
            179.2788848,
            166.5154683,
            192.3716889,
            234.3239266,
            245.232115,
            190.3896784,
            272.7582826)

df <- data.frame(VAF, Median)
df

coefficient <- cor.test(df$VAF, df$Median)
coefficient$estimate

model <- lm(Median ~ VAF, data = df)
summary(model)
sigma(model)
summary(model)$r.squared

library(ggpubr)
ggscatter(df,x = "VAF", y = "Median", add = "reg.line", 
          add.params = list(color = "#ff7b7b", fill = "lightgray"),
          conf.int = T,  point = "FALSE") +
  #geom_point(alpha = 0.5)+
  geom_jitter(width = 0.03)+
  ylab("%pSTAT3(Y705)+/CD33+")+
  stat_cor(method = "pearson", geom = "label") +
  ggtitle("pSTAT3 activation correlates with CH VAF")+
  theme(plot.title = element_text(face = "bold"), panel.grid.major = element_blank(), panel.background = element_blank())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot of variant locations #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#


voi.ch <- read.csv('~/21_014_filtered_mito_variants.csv')$variant
message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))


# Wrangle variants
voi<- read.csv('~/coverage_table_for_loliplots_v2.csv')
voi_order <- voi[order(as.numeric(cutf(voi$mtvar, d = "_"))),]
voi_order.ch <- voi_order$mtvar
voi_order.cov <- voi_order$coverage
voi_along_chrM.num <- rep(0, 16569)
voi_along_chrM.num[as.numeric(cutf(voi_order.ch, d = "_"))] <- 1
labels.ch <- rep("", 16569)
labels.ch[as.numeric(cutf(voi_order.ch, d = "_"))] <- voi_order.ch

# Colors


library(RColorBrewer)
library(dplyr)

barcol.ch <- c("#1F78B4", "#E31A1C", "#33A02C", "#6A3D9A", "#FB9A99",  "#B2DF8A")
# Genes
GenePos.tib <- tibble(Names = c("MT.ATP6", "MT.ATP8", "MT.CO1", "MT.CO2", "MT.CO3", "MT.CYB", "MT.ND1", "MT.ND2", "MT.ND3",
                                "MT.ND4", "MT.ND4L", "MT.ND5", "MT.ND6", "MT.RNR1", "MT.RNR2"),
                      start = c(8527, 8366, 5904, 7586, 9207, 14747, 3307, 4470, 10059, 10760, 10470, 12337, 14149, 648, 1671), 
                      end = c(9207, 8572, 7445, 8269, 9990, 15887, 4262, 5511, 10404, 12137, 10766, 14148, 14673, 1601, 3229))
GenePos.tib <- GenePos.tib %>% arrange(start) %>%
  mutate(mid = round((end-start)/2+start,0), ycoord = rep(c(-0.05,-0.07), length.out = 15))

pdf("~/variant_locations_v3.pdf", height = 3, width = 10)
par(xpd=T)
barplot(voi_along_chrM.num, border = NA,
        yaxt = "n", xlab = "position", ylab = "coverage",
        cex.names = 0, space = 0, ylim = c(0,70))
axis(side = 1, at = c(0, 16569))
axis(side = 2, at = c(0,70))
points(x = as.numeric(cutf(voi_order.ch, d = "_")), y = voi_order.cov,
       pch = 16, col = barcol.ch)
segments(x0 = 0, x1 = 16569, y0 = 5, y1 = 5, lty = "dotted", col = "grey" )
legend(x = 0, y = 110, legend = voi_order$CHIVEID, bty = 'n', text.col = barcol.ch, cex = 0.7)
for (n in 1:nrow(GenePos.tib)) {
 arrows(x0 = GenePos.tib$end[n], x1 = GenePos.tib$start[n], y0 = GenePos.tib$ycoord[n], y1 = GenePos.tib$ycoord[n], code = 3, length = 0.05, col = "blue")
  text(x = GenePos.tib$mid[n], y = GenePos.tib$ycoord[n],
       labels = gsub("MT.", "", GenePos.tib$Names[n]), pos = 1, offset = runif(1), col = "blue", cex = 0.7)
}

textcol.ch <- na.omit(barcol.ch)
#attributes(textcol.ch) <- NULL
text(x = as.numeric(cutf(voi_order.ch, d = "_")), y = voi_order.cov, labels = voi_order.ch, pos = 4, srt = 45,
     col = textcol.ch, cex = 0.7)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Venn Diagram        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(VennDiagram)
library(tidyverse)
library(RColorBrewer)

mono_WT <- read.csv('~/MAST_DE_Mutant_vs_Wildtype/TET2 CD14 Mutant vs Wildtype.csv')
mono_ctrl <- read.csv('~/MAST_DE_Mutant_vs_Control/TET2 CD14 Mutant vs Control.csv')

colors <- brewer.pal(n=10,"Paired")
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"

mono_WT <- mono_WT %>% filter(label == "TRUE" & coef > 0)
length(mono_WT$primerid)


mono_ctrl <- mono_ctrl %>% filter(label == "TRUE" & coef > 0)

length(mono_ctrl$primerid)

venn.diagram(x = list(mono_WT$primerid, mono_ctrl$primerid),
             category.names = c(paste0("Wildtype comparison", " (",length(mono_WT$primerid), ")"),
                                paste0("Control comparison"," (", length(mono_ctrl$primerid), ")")),
             main = "Upregulated genes in TET2 Mutant CD14+ Comparisons",
             output= TRUE, 
             fill = c(alpha("#33A02C",0.3),alpha('#B2DF8A',0.3)),
             height = 5 , 
             width = 5 , 
             resolution = 600,
             units = "in",
             compression = "lzw",
             imagetype = 'svg',
             cat.col = c("#33A02C", '#B2DF8A'),
             cat.pos = c(180, 180),
             lwd = 0.5,
             filename = 'C:/Users/bhatp3/Downloads/wt_ctrl.svg')

overlap <- calculate.overlap(x = list("WT" = mono_WT$primerid, "Ctrl" = mono_ctrl$primerid))
write.csv(overlap$a3, 'tet2_CD14_overlap_siggenes.csv')
write.csv(overlap$a2, 'tet2_CD14_overlap_siggenes.csv')
write.csv(overlap$a1, 'tet2_CD14_overlap_siggenes.csv')


mono_WT <- read.csv('~/DNMT3A CD14 Mutant vs Wildtype.csv')
mono_ctrl <- read.csv('~/DNMT3A CD14 Mutant vs Control.csv')
mono_WT <- mono_WT %>% filter(label == "TRUE" & coef > 0)
mono_WT
mono_ctrl <- mono_ctrl %>% filter(label == "TRUE" & coef > 0)
mono_ctrl


venn.diagram(x = list(mono_WT$primerid, mono_ctrl$primerid),
             category.names = c(paste0("WT comparison", " (",length(mono_WT$primerid), ")"),
                                paste0("Ctrl comparison"," (", length(mono_ctrl$primerid), ")")),
             main = "Upregulated genes in DNMT3A Mutant CD14+ Comparisons",
             output= TRUE, 
             fill = c(alpha("#33A02C",0.3),alpha('#B2DF8A',0.3)),
             height = 5 , 
             width = 5 , 
             resolution = 600,
             units = "in",
             compression = "lzw",
             cat.col = c("#33A02C", '#B2DF8A'),
             cat.pos = c(180, 180),
             lwd = 0.5,
             filename = 'C:/Users/bhatp3/Downloads/dnm_wt_ctrl.png')

mono_tet2 <- read.table('~/CD14 Monos_TET2_vs_none.tsv', header = T)
mono_dnm <- read.table('~/CD14 Monos_DNMT3A_vs_none.tsv', header = T)
mono_tet2 <- mono_tet2 %>% filter(log2FoldChange > 0.25 & p.adjust(pvalue) > 0.05)
mono_tet2
mono_dnm <- mono_dnm %>% filter(log2FoldChange > 0.25 & p.adjust(pvalue) > 0.05)
mono_dnm


venn.diagram(x = list(mono_tet2$gene, mono_dnm$gene),
             category.names = c(paste0("TET2 comparison", " (",length(mono_tet2$gene), ")"),
                                paste0("DNMT3A comparison"," (", length(mono_dnm$gene), ")")),
             main = "Upregulated genes in DNMT3A vs TET2 CD14+ Comparisons",
             output= TRUE, 
             fill = c(alpha("#33A02C",0.3),alpha('#B2DF8A',0.3)),
             height = 5 , 
             width = 5 , 
             resolution = 600,
             units = "in",
             compression = "lzw",
             cat.col = c("#33A02C", '#B2DF8A'),
             cat.pos = c(180, 180),
             lwd = 0.5,
             filename = 'C:/Users/bhatp3/Downloads/dnm_tet2.png')
overlap <- calculate.overlap(x = list("tet2" = mono_tet2$gene, "dnmt3a" = mono_dnm$gene))
write.csv(overlap$a3, 'tet2_dnm_overlap_siggenes.csv')
write.csv(overlap$a2, 'tet2_CD14_overlap_siggenes.csv')
write.csv(overlap$a1, 'tet2_CD14_overlap_siggenes.csv')
