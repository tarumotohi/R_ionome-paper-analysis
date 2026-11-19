###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(qtl)


##### 0.2. Define Parameters #####
scriptID <- "9.2."


##### 0.3. Define Directories #####
dataDir <- "data/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
fileNameGenoMat <- paste0(dataDir, "genoMat.csv")
genoMat <- read.csv(fileNameGenoMat,  row.names = 1)

geno <- genoMat[nchar(rownames(genoMat)) == 7 | nchar(rownames(genoMat)) == 24, ]

Xm_all <- as.matrix(geno)
Xm_all <- Xm_all - 1

Km_all <- tcrossprod(Xm_all) / ncol(Xm_all)

Km_RIL <- Km_all[nchar(rownames(Km_all)) == 7, nchar(colnames(Km_all)) == 7]

fileNameKm <- paste0(saveDir, scriptID, "Km.RDS")
saveRDS(object = Km_RIL, file = fileNameKm)



