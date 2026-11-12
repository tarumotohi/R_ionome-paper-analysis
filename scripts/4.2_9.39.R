###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####


##### 0.2. Define Parameters #####
scriptID <- "9.39."


##### 0.3. Define Directories #####
dataDir <- "data/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameGenoMat <- paste0(dataDir, "genoMat.csv")
genoMat <- read.csv(fileNameGenoMat,  row.names = 1)

geno <- genoMat[nchar(rownames(genoMat)) == 7 | nchar(rownames(genoMat)) == 24, ]

Xm <- as.matrix(geno)
Xm <- Xm - 1

Km <- tcrossprod(Xm) / ncol(Xm)

fileNameXm <- paste0(saveDir, scriptID, "Xm.RDS")
saveRDS(object = Xm, file = fileNameXm)

fileNameKm <- paste0(saveDir, scriptID, "Km.RDS")
saveRDS(object = Km, file = fileNameKm)



