###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####


##### 0.2. Define Parameters #####
scriptID <- "9.1."


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

dm.ionRD19 <- ion.pca$RD2019@completeObs
rownames(dm.ionRD19) <- gsub("-", "N", rownames(dm.ionRD19))

dm.ionRC20 <- ion.pca$RC2020@completeObs
rownames(dm.ionRC20) <- gsub("-", "N", rownames(dm.ionRC20))

dm.ionRD20 <- ion.pca$RD2020@completeObs
rownames(dm.ionRD20) <- gsub("-", "N", rownames(dm.ionRD20))

fileNameion.RC19 <- paste0(saveDir, scriptID, "ion.RC19.RDS")
saveRDS(object = dm.ionRC19, file = fileNameion.RC19)

fileNameion.RD19 <- paste0(saveDir, scriptID, "ion.RD19.RDS")
saveRDS(object = dm.ionRD19, file = fileNameion.RD19)

fileNameion.RC20 <- paste0(saveDir, scriptID, "ion.RC20.RDS")
saveRDS(object = dm.ionRC20, file = fileNameion.RC20)

fileNameion.RD20 <- paste0(saveDir, scriptID, "ion.RD20.RDS")
saveRDS(object = dm.ionRD20, file = fileNameion.RD20)


##### 1.2. Kernel #####
Xi.RC19 <-scale(dm.ionRC19)
Xi.RD19 <-scale(dm.ionRD19)
Xi.RC20 <-scale(dm.ionRC20)
Xi.RD20 <-scale(dm.ionRD20)

Ki.RC19 <- tcrossprod(Xi.RC19) / ncol(Xi.RC19)
Ki.RD19 <- tcrossprod(Xi.RD19) / ncol(Xi.RD19)
Ki.RC20 <- tcrossprod(Xi.RC20) / ncol(Xi.RC20)
Ki.RD20 <- tcrossprod(Xi.RD20) / ncol(Xi.RD20)

h <- 2
d2 <- as.matrix(dist(Xi.RC19)^2)
Kig.RC19.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

d2 <- as.matrix(dist(Xi.RD19)^2)
Kig.RD19.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

d2 <- as.matrix(dist(Xi.RC20)^2)
Kig.RC20.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

d2 <- as.matrix(dist(Xi.RD20)^2)
Kig.RD20.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

fileNameKi.RC19 <- paste0(saveDir, scriptID, "Ki.RC19.RDS")
saveRDS(object = Ki.RC19, file = fileNameKi.RC19)

fileNameKi.RD19 <- paste0(saveDir, scriptID, "Ki.RD19.RDS")
saveRDS(object = Ki.RD19, file = fileNameKi.RD19)

fileNameKi.RC20 <- paste0(saveDir, scriptID, "Ki.RC20.RDS")
saveRDS(object = Ki.RC20, file = fileNameKi.RC20)

fileNameKi.RD20 <- paste0(saveDir, scriptID, "Ki.RD20.RDS")
saveRDS(object = Ki.RD20, file = fileNameKi.RD20)

fileNameKig.RC19.2 <- paste0(saveDir, scriptID, "Kig.RC19.2.RDS")
saveRDS(object = Kig.RC19.2, file = fileNameKig.RC19.2)

fileNameKig.RD19.2 <- paste0(saveDir, scriptID, "Kig.RD19.2.RDS")
saveRDS(object = Kig.RD19.2, file = fileNameKig.RD19.2)

fileNameKig.RC20.2 <- paste0(saveDir, scriptID, "Kig.RC20.2.RDS")
saveRDS(object = Kig.RC20.2, file = fileNameKig.RC20.2)

fileNameKig.RD20.2 <- paste0(saveDir, scriptID, "Kig.RD20.2.RDS")
saveRDS(object = Kig.RD20.2, file = fileNameKig.RD20.2)



