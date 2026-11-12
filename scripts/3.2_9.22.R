###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(Matrix)


##### 0.2. Define Parameters #####
scriptID <- "9.22."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameKm <- paste0(dataDir, "9.2.Km.RDS")
Km <- readRDS(file = fileNameKm)


##### 1.2. Modify Data #####
n <- nrow(Km)
labels <- rownames(Km)

env_suffixes <- c("_19RC", "_19RD", "_20RC", "_20RD")

new_labels <- as.vector(outer(labels, env_suffixes, paste0))

Km_all <- kronecker(matrix(1, 4, 4), Km)

rownames(Km_all) <- new_labels
colnames(Km_all) <- new_labels

fileNameKm <- paste0(saveDir, scriptID, "Km.all.RDS")
saveRDS(object = Km_all, file = fileNameKm)



