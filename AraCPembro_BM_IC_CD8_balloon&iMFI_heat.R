rm(list = ls())

library(scales)
library(gplots)
library(reshape2)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory


load("workspace_CD8_APembro_IC.rds")


fcs_dir <- "fcs_IC"
if(!dir.exists(fcs_dir)){
  dir.create(fcs_dir)
  setwd(fcs_dir)
}

merged <- sce2fcs(sce, split_by = NULL, keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.FCS(merged, filename = "merged.fcs")
fs <- sce2fcs(sce, split_by = "cluster_annotation", keep_cd = TRUE, keep_dr = FALSE, assay = "counts")
write.flowSet(fs, outdir = fcs_dir)


df <- as.data.frame(colData(sce))
summary(df)

cond_cluster_df <- df %>% select(condition, cluster_annotation)
freq_table <- table(cond_cluster_df$condition, cond_cluster_df$cluster_annotation)
tot_CR_DG <- sum(freq_table["DG_CR", ])
tot_CR_post <- sum(freq_table["post_CR", ])
tot_NR_DG <-  sum(freq_table["DG_NR", ])
tot_NR_post <-  sum(freq_table["post_NR", ])

freq_table["DG_CR", ] <- freq_table["DG_CR", ]/tot_CR_DG*100
freq_table["post_CR", ] <- freq_table["post_CR", ]/tot_CR_post*100
freq_table["DG_NR", ] <- freq_table["DG_NR", ]/tot_NR_DG*100
freq_table["post_NR", ] <- freq_table["post_NR", ]/tot_NR_post*100


freq_table <- t(freq_table)
freq_table <- freq_table[c(3,2,1,6,8,5,9,4,10,7),]

dev.off()
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "mean",
                scale = "last", bars = TRUE, perc = TRUE)
freqdf <- cbind(rownames(freq_table), freq_table[, 1], freq_table[, 2], freq_table[, 3], freq_table[, 4])
colnames(freqdf) <- c("Cluster", "DG_CR", "DG_NR", "post_CR", "post_NR")
write.csv(freqdf, file = "freq_table.csv", row.names = FALSE)
balloon <- read.table(file.choose(), header = TRUE, sep = ",", stringsAsFactors = FALSE)
balloon$Cluster <- factor(balloon$Cluster, levels = rev(unique(balloon$Cluster)))
balloon_melted <- melt(balloon, sort = FALSE)

p <- ggplot(balloon_melted, aes(x = variable ,y = Cluster))

pp <- p + geom_point(aes(size = value), colour = "grey34") + theme(panel.background = element_blank()) +
  scale_size_area(max_size = 12)

pp

matrix_r <- read.table(file.choose(), sep = ",", stringsAsFactors = FALSE, header = TRUE)
matrix_r <- matrix_r[-c(4,12,13), ]
matrix_r <- matrix_r[,-c(2:5,18,19)]
matrix_r
colnames(matrix_r)
col_later <- gsub(".fcs", "", matrix_r$X)
matrix_r <- matrix_r[, c(2:25)]
colnames(matrix_r)
matrix_r[, grep("Freq", colnames(matrix_r))]
freqs <- matrix_r %>% select(grep("Freq", colnames(matrix_r)))
MFI <- matrix_r %>% select(grep("Mean", colnames(matrix_r)))
iMFI <- freqs*MFI
colnames(iMFI) <- c("CD25", "CD27", "CD28", "CD45RA", "CD57", "EOMES", "GRZB",
                    "ICOS", "KI67", "PD1", "TBET", "TCF1")
scaled_iMFI <- iMFI
for(i in c(1:ncol(iMFI))){
  scaled_iMFI[, i] <- rescale(iMFI[, i], to = c(0, 100))
}

heat.colors<-colorRampPalette(c("navy","blue4","blue","skyblue","khaki1","lightgoldenrod1","goldenrod1","orange"))(100)
##### Make a heatmap  ##change parameters if needed (type help(heatmap.2) to get a full description of the options)
matrix_r <- scaled_iMFI
scale.data <- as.matrix((t(matrix_r)-apply(t(matrix_r),1,mean))/apply(t(matrix_r),1,sd))
colnames(scale.data) <- col_later
##### Plot heatmap 
dev.off()
heatmap.2(as.matrix(t(scale.data)),
          dendrogram="both", scale="none",  na.color="grey",
          col = heat.colors, trace = "none", labRow = rownames(t(scale.data)), key = TRUE, keysize = 1, cexCol=1,
          density.info= "none", symkey=FALSE, margins=c(7,10), main="iMFI", 
          xlab="Markers", ylab= "CLUSTERS")

