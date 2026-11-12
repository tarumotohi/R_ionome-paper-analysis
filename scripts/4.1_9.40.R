###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(stringr)
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.40."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
#RIL
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

dm.ionome_RIL <- bind_rows(dm.ionRC19, dm.ionRD19, dm.ionRC20, dm.ionRD20)

#F2
fileNameionome_F2Raw <- paste0("data/", "2021_Tottori_RIL_IonomeXray.csv")
ionome_F2Raw <- read.csv(fileNameionome_F2Raw)

ionome_F2 <- ionome_F2Raw[nchar(ionome_F2Raw$variety) == 15, ]

substr(ionome_F2$variety, 4, 4) <- "N"
substr(ionome_F2$variety, 12, 12) <- "N"

ionome_F2$plot <- str_pad(ionome_F2$plot, width = 3, side = "left", pad = "0")

ionome_F2 <- ionome_F2[!rownames(ionome_F2) %in% c("732", "1136"), ]
rownames(ionome_F2) <- paste0(
  ionome_F2$variety, "_",
  ionome_F2$block, "-",
  ionome_F2$plot, "_",
  ionome_F2$ind
)

ionome_F2 <- ionome_F2[ , -(1:9)]

ionome_F2 <- ionome_F2[rowSums(!is.na(ionome_F2)) > 0, ]

dm.ionome_all <- rbind(dm.ionome_RIL, ionome_F2)


##### 1.2. Kernel #####
Xi.all <-scale(dm.ionome_all)

Ki.all <- tcrossprod(Xi.all) / ncol(Xi.all)

h <- 2
d2 <- as.matrix(dist(Xi.all)^2)
Kig.all.2 <- exp(- h * d2 / median(d2[upper.tri(d2)]))

fileNameKi.all <- paste0(saveDir, scriptID, "Ki.RIL_F2.RDS")
saveRDS(object = Ki.all, file = fileNameKi.all)

fileNameKig.all.2 <- paste0(saveDir, scriptID, "Kig.RIL_F2.2.RDS")
saveRDS(object = Kig.all.2, file = fileNameKig.all.2)



