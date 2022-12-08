
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(wesanderson)
library(ggrepel)
library(xlsx)
library(ggradar)
library(future)
library(UpSetR)
library(openxlsx)
library(clusterProfiler)
library(DoubletFinder)
library(ggnewscale)
library(ggh4x)
library(ggpubr)


# Set directory -----------------------------------------------------------

directory = "/home/rstudio/bicklab-aps/edit/Single_Cell_CHIP/"
setwd(directory)
