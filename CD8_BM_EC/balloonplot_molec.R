rm(list = ls())

library(gplots)
library(scales)
library(reshape2)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

load("AraCPembro_step2_CD8_BM_EC.rds")

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


sce$cluster_annotation <- cluster_ids(sce, "cluster_annotation")
df <- data.frame(colData(sce))
head(df)
list <- list(apopt, oncomet, onc_sign, epi_dys)
all.frames <- lapply(seq_along(list), function(i) cbind(df, list[[i]]))
df <-  do.call(rbind,all.frames)
colnames(df)[which(colnames(df) == 'list[[i]]')] <- "condition_mol"

cond_cluster_df <- df %>% select(cluster_annotation, condition_mol)
freq_table <- table(cond_cluster_df$cluster_annotation, cond_cluster_df$condition_mol)

tot_cl1 <- sum(freq_table["Cluster 1", ])
tot_T_act <- sum(freq_table["T act", ])
tot_T_act_eff <-  sum(freq_table["T act-eff", ])
tot_T_senescent <-  sum(freq_table["T senescent", ])
tot_T_naive <- sum(freq_table["T naive", ])
tot_TEMRA_CD57posCD28pos <- sum(freq_table["TEMRA CD57+CD28+", ])
tot_T_CD57neg<-  sum(freq_table["TEMRA CD57-", ])
tot_TEMRA_CD57posCD28neg <-  sum(freq_table["TEMRA CD57+CD28-", ])

freq_table <- t(freq_table)

dev.off()
plotExprHeatmap(sce, features = type_markers(sce), k = "cluster_annotation", by = "cluster_id",  fun = "mean",
                scale = "last", bars = TRUE, perc = TRUE)
freqdf <- cbind(rownames(freq_table), freq_table[, 1], freq_table[, 2], freq_table[, 3], freq_table[, 4],
                freq_table[, 5], freq_table[, 6], freq_table[, 7], freq_table[, 8])
colnames(freqdf) <- c("Condition", colnames(freq_table))
write.csv(freqdf, file = "freq_table.csv", row.names = FALSE)
balloon <- read.table(file.choose(), header = TRUE, sep = ",", stringsAsFactors = FALSE)
balloon$Cluster <- factor(balloon$Cluster, levels = rev(unique(balloon$Cluster)))
balloon_melted <- melt(balloon, sort = FALSE)

p <- ggplot(balloon_melted, aes(x = variable ,y = Condition))

pp <- p + geom_point(aes(size = value), colour="black", fill="skyblue") + theme(panel.background = element_blank()) +
  scale_size_area(max_size = 12)

pp

p +geom_point(aes(size=value), shape=21, colour="black", fill="skyblue")+
  theme(panel.background=element_blank(), panel.border = element_rect(colour = "blue", fill=NA, size=1))+
  scale_size_area(max_size=20)
