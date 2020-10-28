rm(list = ls())

###############################
## INSTALL REQUIRED PACKAGES ##
###############################

if(!require('rstudioapi')) {
  install.packages('rstudioapi')
}
if(!require('Biobase')){
  install.packages('Biobase')
}
if(!require('flowCore')){
  install.packages("flowCore")
}
if(!require('FlowSOM')){
  install.packages("FlowSOM")
}
if(!require('ggplot2')){
  install.packages("ggplot2")
}
if(!require('dplyr')){
  install.packages('dplyr')
}
if(!require("flowVS")){
  install.packages(file.choose(), repos = NULL, type = "source")
}
if(!require('flowStats')){
  BiocManager::install('flowStats')
}
if(!require('flowSOMworkshop')){
  BiocManager::install('flowSOMworkshop')
}
if(!require('stringr')){
  BiocManager::install('stringr')
}

if(!require('cytofCore')){
  BiocManager::install('cytofCore')
}

###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2) 
  library(readxl) 
  library(flowCore) 
  library(flowDensity)
  library(flowAI) 
  library(FlowSOM) 
  library(FlowSOMworkshop) 
  library(flowStats)
  library(CytoNorm)
  library(PeacoQC) 
  library(cytofCore)
  library(flowVS)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define fcs_directory
csv <- "csv"
csvDirectory <- paste(PrimaryDirectory, csv, sep = "/")
dir.create(csvDirectory)

# Define workingDirectory
wdName <- "Working_Directory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#Define csv2fcsDirectory
csv2fcsDirectory <- "csv2fcsDirectory"
csv2fcsDirectory <- paste(PrimaryDirectory, csv2fcsDirectory, sep = "/")
dir.create(csv2fcsDirectory)

#List CSV files
CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

for(i in c(1:length(CSVfiles))){
  data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
  print(CSVfiles[i])
  print(fileName[i])
  cytofCore.write.FCS(as.matrix(data), 
                      filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                      what = "numeric")
}


# read flowSet
# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
fs <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)
colnames(fs)

QC_dir <- "QC"
preprocessed_dir <- "Preprocessed"
if(!dir.exists(QC_dir)){
  dir.create(QC_dir)
  dir.create(preprocessed_dir)
  dir.create(file.path(QC_dir, "flowAI"))
  dir.create(file.path(QC_dir, "PeacoQC"))
}

cellCount <- rep(NA, length(FCSfiles))
names(cellCount) <- FCSfiles

for (file in FCSfiles){
  ff <- read.FCS(file.path(csv2fcsDirectory, file))
  cellCount[file] <- nrow(ff)
  # Run flowAI on the samples
  resQC <- flow_auto_qc(fcsfiles = ff,
                        folder_results = file.path(QC_dir, "flowAI"),
                        output = 1)
  # Run PeacoQC
  resQC <- PeacoQC(ff = ff,
                   determine_good_cells = "all",
                   channels = c(7:22),
                   plot = TRUE,
                   output_folder = file.path(QC_dir, "PeacoQC"))
  ff <- ff[resQC$GoodCells,]
  saveRDS(resQC, file.path(QC_dir, "PeacoQC", gsub(".fcs", "_QC.RDS", file)))
  write.FCS(ff, 
            file.path(preprocessed_dir, file))
}

fs <- read.flowSet(list.files("Preprocessed",
                              pattern = ".fcs",
                              full.names = T))


colnames(fs)[7:22] <- c('TCF1', 'EOMES', 'CD27', 'CD45RA', 'CD4', 'FOXP3', 
                        'CD3', 'CD57', 'CD28', 'TBET', 'ICOS', 'GRZB', 
                        'PD1', 'CD25', 'L_D', 'KI67')

channels <- colnames(fs)
flowViz.par.set(theme =  trellis.par.get(), reset = TRUE)
as.vector(channels)
densityplot(~TCF1+EOMES+CD27+CD45RA+CD4+FOXP3+CD3+CD57+CD28+TBET+ICOS+GRZB+PD1+CD25+L_D+KI67, fs)

save(list = ls(), file =("workspace_CD8_BM_IC.rds"))
