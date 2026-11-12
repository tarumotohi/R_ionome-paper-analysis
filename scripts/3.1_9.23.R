###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.23."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameIonPCA <- paste0(dataDir, "9.0.ionRILMean.PCA.RDS")
ion.pca <- readRDS(fileNameIonPCA)

dm.ionRC19 <- ion.pca$RC2019@completeObs
rownames(dm.ionRC19) <- gsub("-", "N", rownames(dm.ionRC19))
dm.ionRC19 <- as.data.frame(dm.ionRC19)

dm.ionRD19 <- ion.pca$RD2019@completeObs
rownames(dm.ionRD19) <- gsub("-", "N", rownames(dm.ionRD19))
dm.ionRD19 <- as.data.frame(dm.ionRD19)

dm.ionRC20 <- ion.pca$RC2020@completeObs
rownames(dm.ionRC20) <- gsub("-", "N", rownames(dm.ionRC20))
dm.ionRC20 <- as.data.frame(dm.ionRC20)

dm.ionRD20 <- ion.pca$RD2020@completeObs
rownames(dm.ionRD20) <- gsub("-", "N", rownames(dm.ionRD20))
dm.ionRD20 <- as.data.frame(dm.ionRD20)

rownames(dm.ionRC19) <- paste0(rownames(dm.ionRC19), "_19RC")
rownames(dm.ionRD19) <- paste0(rownames(dm.ionRD19), "_19RD")
rownames(dm.ionRC20) <- paste0(rownames(dm.ionRC20), "_20RC")
rownames(dm.ionRD20) <- paste0(rownames(dm.ionRD20), "_20RD")

dm.ionome_all <- bind_rows(dm.ionRC19, dm.ionRD19, dm.ionRC20, dm.ionRD20)

fileNamedm.ionome_all <- paste0(saveDir, scriptID, "ionome_all.RDS")
saveRDS(object = dm.ionome_all, file = fileNamedm.ionome_all)


##### 1.2. Kernel #####
Xi.all <-scale(dm.ionome_all)

Ki.all <- tcrossprod(Xi.all) / ncol(Xi.all)

h <- 2
d2 <- as.matrix(dist(Xi.all)^2)
Kig.all.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

fileNameKi.all <- paste0(saveDir, scriptID, "Ki.all.RDS")
saveRDS(object = Ki.all, file = fileNameKi.all)

fileNameKig.all.2 <- paste0(saveDir, scriptID, "Kig.all.2.RDS")
saveRDS(object = Kig.all.2, file = fileNameKig.all.2)



