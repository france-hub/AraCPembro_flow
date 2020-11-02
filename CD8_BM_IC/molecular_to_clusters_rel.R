suppressPackageStartupMessages({
  library(rstudioapi)
  library(devtools)
  library(flowCore)
  library(cytofCore)
  library(FlowSOM)
  library(cluster)
  library(ggplot2)
  library(dplyr)
  library(ggthemes)
  library(RColorBrewer)
  library(uwot)
  library(CATALYST)
  library(diffcyt)
  library(stringr)
  library(scran)
  library(scater)
  library(ggcyto)
  library(SingleCellExperiment)
  library(flowWorkspace)
  library(reshape2)
  library(ggrepel)
  library(slingshot)
  library(knn.covertree)
  library(readxl)
  library(flowStats)
  library(ggpubr)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

load("AraCPembro_step2_CD8_BM_IC.rds")

ind_apopt <- which(colData(sce)$patient_id %in% c('148', '151'))
ind_oncomet <- which(colData(sce)$patient_id %in% c('108', '140'))
ind_onc_sign <- which(colData(sce)$patient_id %in% c('121','142','113','125','140','113','130','152','145','108'))
ind_epi_dys <-  which(colData(sce)$patient_id %in% c('108','113','118','125','140'))
apopt <- rep("no_apopt", length(colData(sce)$patient_id))
apopt[ind_apopt] <- "yes_apopt"
oncomet <- rep("no_oncomet", length(colData(sce)$patient_id))
oncomet[ind_oncomet] <- "yes_oncomet"
onc_sign <- rep("no_onc_sign", length(colData(sce)$patient_id))
onc_sign[ind_onc_sign] <- "yes_onc_sign"
epi_dys <- rep("no_epi_dys", length(colData(sce)$patient_id))
epi_dys[ind_epi_dys] <- "yes_epi_dys"

colData(sce) <- cbind(colData(sce), apopt, oncomet, onc_sign, epi_dys)

names(colData(sce))

CATALYST::pbMDS(sce, color_by = "apopt", label = "patient_id")
CATALYST::pbMDS(sce, color_by = "oncomet", label = "patient_id")
CATALYST::pbMDS(sce, color_by = "onc_sign", label = "patient_id")
CATALYST::pbMDS(sce, color_by = "epi_dys", label = "patient_id")

pbMDS(sce, by = "both", k = "cluster_annotation", shape_by = "condition", size_by = TRUE)

# Plot UMAP
#apopt (pt 148, 151)
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("apopt") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
#oncomet (pt 108,140)
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("oncomet") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
#onc_sign (pt 121,142,113,125,140,113,130,152,145,108)
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("onc_sign") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("onc_sign")
plot + geom_density2d(binwidth = 0.002, colour = "black")
#UMAP molecular signatures-responses
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap(c("onc_sign", "condition")) + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions

#epi_dys (108,113,118,125,140)
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("epi_dys") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap(c("epi_dys", "condition")) + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions

plotDR(sce, "UMAP", color_by = "cluster_annotation", facet_by = "sample_id") #by sample_id

#Plot Abundancies
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "apopt", shape_by = NULL)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "oncomet", shape_by = NULL)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "onc_sign", shape_by = NULL)
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", group_by = "epi_dys", shape_by = NULL)

ind_onc_sign_exp <- which(sce@metadata$experiment_info$patient_id %in% c('121','142','113','125','140','113','130','152','145','108'))
onc_sign_exp <- rep("no_onc_sign", length(sce@metadata$experiment_info$patient_id))
onc_sign_exp[ind_onc_sign_exp] <- "yes_onc_sign"
sce@metadata$experiment_info <- cbind(sce@metadata$experiment_info, onc_sign_exp)
colnames(sce@metadata$experiment_info) <- c("sample_id", "condition", "patient_id", "n_cells", "onc_sign")

# Statistical Analysis
ei <- sce@metadata$experiment_info

design <- createDesignMatrix(
  ei, cols_design = c("onc_sign", "patient_id")
)

head(design) #check design matrix
ncol(design)
#which(colnames(design) == "patient_id152") #matrix not full rank: patient_id152 not estimable
design <- design[,-c(19)] # not a full rank matrix, drop column 19

#(not sig)
contrast <- createContrast(c(0, 1, rep(0, 18)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)

