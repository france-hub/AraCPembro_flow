rm(list = ls())

###############################
### LOAD REQUIRED PACKAGES ####
###############################
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

#load("AraCPembro_step1_CD4_BM_EC.rds")

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
antigen <- fcs_colname
marker_class <- c(rep("type", 4), "state", "type", "state", rep("type", 8), "state")
length(marker_class) == length(fcs_colname)

# Create panel
panel <- as.data.frame(cbind(fcs_colname, antigen, marker_class))
all(panel$fcs_colname %in% colnames(fs)) #check

##Building metadata dataframe
#Set conditions
condition <- FCSfiles
condition <- word(condition, 2,3, sep = "_") 

#Set patient_id
patient_id <- FCSfiles
patient_id <- word(patient_id, 1, sep = "_")

#Set sample_id
sample_id <- paste(patient_id, condition, sep = "_")

#Set file_name
file_name <- FCSfiles

#Create metadata dataframe
md <- cbind(file_name, sample_id, condition, patient_id)
md <- data.frame(md)

#Check if ids and md$file_name are the same
ids <- c(fsApply(fs, identifier))
ids%in%md$file_name

# create SingleCellExperiment object
sce <- CATALYST::prepData(fs, panel = panel, md = md, transform = FALSE, features = panel$fcs_colname)
sce@assays@data$exprs <- sce@assays@data$counts

## QC
# Density plots
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 4
p

# Check number of cells for each sample
n_cells(sce)
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Pseudo-bulk multidimensional scaling (pbMDS)
CATALYST::pbMDS(sce, color_by = "condition", label = "patient_id")

# Non Redundancy Score (NRS)
plotNRS(sce, features = type_markers(sce), color_by = "condition")

# Run FlowSOM and ConsensusClusterPlus
set.seed(1234)
sce <- cluster(sce,
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)
delta_area(sce)
plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "meta14", 
                bars = TRUE, perc = TRUE, col_anno = TRUE, scale = "last", bin_anno = TRUE)

#A lot of clusters < 1%, focus on most frequent clusters
sce <- filterSCE(sce, cluster_id %in% c(1,2,3,4,8,13,14), k = "meta14")

# Run dimensionality reduction - UMAP
n_cells <- 3000
n_events <- min(n_cells(sce))
if(!(n_cells < n_events))
  n_cells <- n_events
sce <- runDR(sce, dr =  "UMAP", cells = n_cells, features = "type")

# Plot UMAP
plotDR(sce, "UMAP", color_by = type_markers(sce)) #markers distribution
plotDR(sce, "UMAP", color_by = "meta11") + facet_wrap("condition") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3))) #different conditions
plotDR(sce, "UMAP", color_by = "meta11", facet_by = "sample_id") #by sample_id

##Add annotations
#Read annotation file
annotation_table <- readxl::read_excel("annotation_CD4_EC.xlsx")
annotation_table

# convert to factor with merged clusters in desired order
annotation_table$new_cluster <- factor(annotation_table$new_cluster)

# Apply manual annotation
sce <- mergeClusters(sce, k = "meta14", 
                     table = annotation_table, id = "cluster_annotation", overwrite = TRUE)

# Plot heatmap with annotations
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta14", m = "cluster_annotation", scale = "last", bin_anno = TRUE)

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "cluster_annotation", scale = "last", bin_anno = TRUE)

#Plot UMAP with annotations
dev.off()
plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot <- plotDR(sce, "UMAP", color_by = "cluster_annotation") + facet_wrap("condition")
plot + geom_density2d(binwidth = 0.006, colour = "black")

#Plot Abundancies
plotAbundances(sce, k = "cluster_annotation", by = "cluster_id", shape_by = NULL)

# Statistical Analysis
ei <- sce@metadata$experiment_info

design <- createDesignMatrix(
  ei, cols_design = c("condition", "patient_id")
)


head(design) #check design matrix

#which(colnames(design) == "patient_id151") #matrix not full rank: patient_id151 not estimable
design <- design[,-c(20)] # not a full rank matrix, drop column 20

#DG NR vs DG CR
contrast <- createContrast(c(0, 1, 0, 0, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)
da_1 <- rowData(out_DA$res)
sub_1 <- filterSCE(sce, condition == c("DG_CR", "DG_NR"))
CR_ind <- grep("CR", levels(sub_1$sample_id))
NR_ind <- grep("NR", levels(sub_1$sample_id))
levels(sub_1$sample_id) <- levels(sub_1$sample_id)[c(CR_ind, NR_ind)]
plotDiffHeatmap(sub_1, da_1, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")

#post CR vs post NR (not sig)
contrast <- createContrast(c(0, 0, -1, 1, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)
da_2 <- rowData(out_DA$res)
sub_2 <- filterSCE(sce, condition == c("post_CR", "post_NR"))
CR_ind <- grep("CR", levels(sub_1$sample_id))
NR_ind <- grep("NR", levels(sub_1$sample_id))
levels(sub_2$sample_id) <- levels(sub_2$sample_id)[c(CR_ind, NR_ind)]
plotDiffHeatmap(sub_2, da_1, all = TRUE, fdr = 0.05, sort_by = "padj", fun = "mean")


#DG CR post_CR (not sig)
contrast <- createContrast(c(0, 0, 1, 0, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)

#NR_DG vs NR_post (not sig)
contrast <- createContrast(c(0, -1, 0, 1, rep(0, 17)))
out_DA <- diffcyt(sce, design = design, contrast = contrast, method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "cluster_annotation", verbose = TRUE, subsampling = TRUE, transform = FALSE)
topTable(out_DA, format_vals = TRUE)

# save workspace 
save(list = ls(), file = "AraCPembro_step2_CD4_BM_EC.rds")
