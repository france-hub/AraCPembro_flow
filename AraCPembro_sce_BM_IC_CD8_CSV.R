rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('devtools')){
  install.packages("devtools")
}
if(!require('flowCore')){
  install.packages("flowCore")
}
if(!require('cytofCore')){
  devtools::install_github("nolanlab/cytofCore")
}
if(!require('cytofkit2')){
  #do not update hexbin
  devtools::install_github("JinmiaoChenLab/cytofkit2")
}
if(!require('FlowSOM')){
  install.packages("FlowSOM")
}
if(!require('cluster')){
  install.packages("cluster")
}
if(!require('Rtsne')){
  install.packages("Rtsne")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages("dplyr")
}
if(!require('ggthemes')){
  install.packages("ggthemes")
}
if(!require('RColorBrewer')){
  install.packages('RColorBrewer')
}
if(!require("uwot")){
  install.packages("uwot")
}
if(!require("CATALYST")){
  BiocManager::install("CATALYST")
}
if(!require("diffcyt")){
  BiocManager::install("diffcyt")
}
if(!require("stringr")){
  BiocManager::install("stringr")
}
if(!require("Rphenograph")){
  BiocManager::install("Rphenograph")
}
if(!require("scran")){
  BiocManager::install("scran")
}
if(!require("scater")){
  BiocManager::install("scater")
}
if(!require("ggcyto")){
  BiocManager::install("ggcyto")
}
if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
}
if(!require("flowWorkspace")){
  BiocManager::install("flowWorkspace")
}
if(!require("reshape2")){
  BiocManager::install("reshape2")
}
if(!require("ggrepel")){
  BiocManager::install("ggrepel")
}
if(!require("slingshot")){
  BiocManager::install("slingshot")
}
if(!require("knn.covertree")){
  devtools::install_github('flying-sheep/knn.covertree')
}
if(!require("theislab/destiny")){
  devtools::install_github('theislab/destiny')
}
if(!require("ggpubr")){
  install.packages('ggpubr')
}


###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
library(rstudioapi)
library(flowCore)
library(cytofCore)
library(cytofkit2)
library(FlowSOM)
library(cluster)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(RColorBrewer)
library(uwot)
library(CATALYST)
library(diffcyt)
library(stringr)
library(Rphenograph)
library(scran)
library(scater)
library(ggcyto)
library(SingleCellExperiment)
library(flowWorkspace)
library(reshape2)
library(ggrepel)
library(slingshot)
library(knn.covertree)
library(destiny)
library(readxl)
library(ggpubr)
})
############################
######## LOAD DATA #########
############################

### Load the (transformed, normalized) FCS files 

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

load("workspace_CD8_BM_IC.rds")

fs
############################
######BUILD SCE OBJECT######
############################

# Keyword ($CYT = "FACS")
names(keyword(fs[[1]]))[15] <- "$CYT"
ds <- keyword(fs[[1]])
l <- list(cyt = "\\$CYT$")
keep <- lapply(l, grep, names(ds))
ds[[keep$cyt]] <- "FACS"
keyword(fs[[1]])[[keep$cyt]] <- "FACS"

##Building panel dataframe
# Define channels of interest and marker_class
fcs_colname <- colnames(fs)[7:22]
marker_class <- c(rep("type", 4), rep("state",3), rep("type",7), "state", "type") 
antigen <- fcs_colname
length(marker_class) == length(fcs_colname)

#Panel
panel <- cbind(fcs_colname, antigen, marker_class)
panel <- as.data.frame(panel)
all(panel$fcs_colname %in% colnames(fs))

##Building metadata dataframe
condition <- FCSfiles
condition <- word(condition, 2,3, sep = "_")
condition[grepl("HD", condition)] <- "HD"

patient_id <- FCSfiles
patient_id <- word(patient_id, 1, sep = "_")

sample_id <- paste(patient_id, condition, sep = "_")

file_name <- FCSfiles

md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

##SCE object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, features = panel$fcs_colname)
sce@assays@data$exprs <- sce@assays@data$counts


#delete HD
sce <- filterSCE(sce, sample_id != "AA_HD")

## QC
# Density
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Counts
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Multidimensional scaling (MDS)
CATALYST::pbMDS(sce, color_by = "condition", label = NULL)

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

## Clustering
#FlowSOM
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)

delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta11", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = 'last')
plotExprHeatmap(sce, features = "type",
                by = "sample_id", k = "meta20", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last", row_anno = "condition")


#sce <- filterSCE(sce, cluster_id != 8, k = "meta8")
set.seed(1234)
n_cells <- 2000
exaggeration_factor <- 12.0
eta <- n_cells/exaggeration_factor
#sce <- runDR(sce, dr = "TSNE", cells = n_cells, features = "type", theta = 0.5, max_iter = 1000, 
             distMethod = "euclidean",
             PCA = TRUE, eta = eta, exaggeration_factor = 12.0)
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")
#sce <- runDR(sce, dr = "DiffusionMap", cells = n_cells, features = "type", assay = "exprs")

## Plots
# tSNE
plotDR(sce, "TSNE", color_by = "condition")
plotDR(sce, "TSNE", color_by = type_markers(sce))
plotDR(sce, "TSNE", color_by = "meta9") + facet_wrap("condition") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "TSNE", color_by = "meta9", facet_by = "sample_id") #by sample_id

# UMAP
plotDR(sce, "UMAP", color_by = "condition")
plotDR(sce, "UMAP", color_by = type_markers(sce)) #markers distribution
plotDR(sce, "UMAP", color_by = "meta11") + facet_wrap("condition") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "UMAP", color_by = "meta11", facet_by = "patient_id") #by sample_id

#Abundances
# Differential abundance
plotAbundances(sce, k = "meta6", by = "sample_id", group_by = "condition")
plotAbundances(sce, k = "meta11", by = "cluster_id", group_by = "condition")

#Diffusion Maps
plotDR(sce, "DiffusionMap", color_by = type_markers(sce), facet_by = "condition")

##Add annotations
#Read annotation file
annotation_table <- readxl::read_excel(file.choose())
annotation_table

# convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster, 
                                       levels = c("TEMRA CD57+", "C2", "T act-dys", "TDT", "T act", "T dys", 
                                                  "C7", "TEMRA", "C9", "T naive", "C11"))

save(list = ls(), file = "workspaceSCE_CD8_IC.rds")
load( "workspaceSCE_CD8_IC.rds")
# apply manual annotation
sce <- mergeClusters(sce, k = "meta11", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)
sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")

#Heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta11", scale = "last")

#Heatmap merging according to cluster_annotation
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "mean",
                scale = "last", bars = TRUE, perc = TRUE)

#Delete cluster < 1%
sce <- filterSCE(sce, cluster_id != "C11", k = "cluster_annotation")


plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "mean",
                scale = "last", bars = TRUE, perc = TRUE, bin_anno = TRUE)

annotation_table <- annotation_table[-11,]
annotation_table$new_cluster <- droplevels(annotation_table$new_cluster, exclude = "C11") 

sce <- mergeClusters(sce, k = "meta11", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

#UMAP with annotations
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#tSNE with annotations
plotDR(sce, "TSNE", color_by = "cluster_annotation") + facet_wrap("condition")

#Diffusion
plotDR(sce, "DiffusionMap", color_by = "cluster_annotation") + facet_wrap("condition")

# Plot abund + annotations
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)

# Statistical Analysis
ei <- sce@metadata$experiment_info

design <- createDesignMatrix(
  ei, cols_design = c("condition", "patient_id")
)

design <- design[,-c(20)] # not a full rank matrix, drop column 19

#DG NR vs DG CR (not sig)
contrast <- createContrast(c(0, 1, 0, 0, rep(0, 17)))
out_DA1 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA1, format_vals = TRUE) #check if significant

#post CR vs post NR (not sig)
contrast <- createContrast(c(0, 0, -1, 1, rep(0, 17)))
out_DA2 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA2, format_vals = TRUE)#check if significant

#DG CR post_CR (TEMRA sig)
contrast <- createContrast(c(0, 0, 1, 0, rep(0, 17)))
out_DA3 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA3, format_vals = TRUE) #check if significant
da_3 <- rowData(out_DA3$res)
#Add sig to boxplots
stat.test_3 <- as_tibble(da_3)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)
p.adj.signif <- c(rep("ns", 7), "4e-02", rep("ns",2))
y.position <- c(42, 10, 40, 15, 55, 20, 10, 30, 7,40)
group1 <- (rep("DG_CR",10))
group2 <- (rep("post_CR", 10))
stat.test_3 <- cbind(stat.test_3, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = "patient_id")
bxp <- bxp + stat_pvalue_manual(stat.test_3, label = "p.adj.signif", tip.length = 0.01, size = 2.5) 
bxp

#NR_DG vs NR_post
contrast <- createContrast(c(0, -1, 0, 1, rep(0, 17)))
out_DA4 <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA4, format_vals = TRUE)
da_4 <- rowData(out_DA4$res)
#Add sig to boxplots
stat.test_4 <- as_tibble(da_4)
p.adj.signif <- c("ns", "3.8e-04", "9.2e-03", "ns", "1.7e-02", rep("ns",4), "4.6e-02")
y.position <- c(42, 10, 40, 15, 55, 20, 10, 30, 7,40)
group1 <- (rep("DG_NR",10))
group2 <- (rep("post_NR", 10))
stat.test_4 <- cbind(stat.test_4, group1, group2, p.adj.signif, y.position)
bxp <- plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = "patient_id")
bxp <- bxp + stat_pvalue_manual(stat.test_4, label = "p.adj.signif", tip.length = 0.01, size = 2.5) 
bxp

# save workspace 
save(list = ls(), file = "workspaceSCE_CD8_IC.rds")
