rm(list = ls())

###############################
### LOAD REQUIRED PACKAGES ####
###############################

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2) 
  library(flowCore) 
  library(cytofCore)
  library(dplyr)
})

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# Define directory that contains csv files
csv <- "csv"
csvDirectory <- paste(PrimaryDirectory, csv, sep = "/")
dir.create(csvDirectory)

# Define working directory
wdName <- "WorkingDirectory"
workingDirectory <- paste(PrimaryDirectory, wdName, sep = "/")
dir.create(workingDirectory)

#Define csv to fcs directory
csv2fcsDirectory <- "csv2fcsDirectory"
csv2fcsDirectory <- paste(PrimaryDirectory, csv2fcsDirectory, sep = "/")
dir.create(csv2fcsDirectory)

#List CSV files and add .fcs extension
CSVfiles <- list.files(csvDirectory, pattern = ".csv$", full = FALSE)
fileName <- gsub(".csv", ".fcs", CSVfiles)

#Obtain fcs files from csv files
for(i in c(1:length(CSVfiles))){
  data <- read.csv(paste(csvDirectory, CSVfiles[i], sep = "/"))
  print(CSVfiles[i])
  print(fileName[i])
  cytofCore.write.FCS(as.matrix(data), 
                      filename = paste(csv2fcsDirectory, fileName[i], sep = "/"),
                      what = "numeric")
}

# Create flowSet from FCSfiles
FCSfiles <- list.files(csv2fcsDirectory, pattern = ".fcs$", full = FALSE)
fs <- read.flowSet(files = FCSfiles, path = csv2fcsDirectory, truncate_max_range = FALSE)

#Look at the colnames(fs) and rename those of interest with the specific markers of interest
colnames(fs) 

colnames(fs)[7:22] <- c('KLRG1', 'CD45RA', 'CD27', 'CD69', 'CD8', 'CD57', 
                        'L_D', 'TIGIT', 'CD28', 'DNAM1', 'CD56', 'CCR7', 
                        'PD1', 'CD25', 'Tim3', 'CD3')

save(list = ls(), file =("AraCPembro_step1_CD8_BM_EC.rds"))

