###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(pcaMethods)


##### 0.2. Define Parameters #####
scriptID <- "9.0."


##### 0.3. Define Directories #####
dataDir <- "data/"
saveDir <- "data/ionome/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameIonome19 <- paste0(dataDir, "2019_Tottori_RIL_IonomeXray_20190909.csv")
ion19 <- read.csv(file = fileNameIonome19, header = TRUE)

fileNameIonome20 <- paste0(dataDir, "2020_Tottori_RIL_IonomeXray.csv")
ion20 <- read.csv(file = fileNameIonome20, header = TRUE)


##### 1.2. Modify Data #####
#### 2019 ####
ion19$plot <- as.factor(
  paste0(ion19$block, formatC(ion19$plot, width = 3, flag = "0")))

X <- ion19[,10:ncol(ion19)]
index <- ion19$plot
plot.mean19 <- apply(X, 2, tapply, index, mean, na.rm = T)

plot_ids <- rownames(plot.mean19)
variety_lookup <- ion19[!duplicated(ion19$plot), c("plot", "variety")]
ionMean19 <- data.frame(plot = plot_ids, plot.mean19)
ionMean19 <- merge(variety_lookup, ionMean19, by = "plot")

selector <- substr(ionMean19$variety, 1, 3) == "E71"
ionRILMean19 <- ionMean19[selector, ]

is.control <- substr(ionRILMean19$plot, 2, 2) == "C"
ionRILMean19RC <- ionRILMean19[is.control,] # Control
ionRILMean19RD <- ionRILMean19[!is.control,] # Drought

fileNameIonRILMean19RC <- paste0(saveDir, scriptID, "ionRILMean.19RC.csv")
write.csv(x = ionRILMean19RC, file = fileNameIonRILMean19RC)

fileNameIonRILMean19RD <- paste0(saveDir, scriptID, "ionRILMean.19RD.csv")
write.csv(x = ionRILMean19RD, file = fileNameIonRILMean19RD)


#### 2020 ####
ion20$plot <- as.factor(
  paste0(ion20$block, formatC(ion20$plot, width = 3, flag = "0")))

X <- ion20[,10:ncol(ion20)]
index <- ion20$plot
plot.mean20 <- apply(X, 2, tapply, index, mean, na.rm = T)

plot_ids <- rownames(plot.mean20)
variety_lookup <- ion20[!duplicated(ion20$plot), c("plot", "variety")]
ionMean20 <- data.frame(plot = plot_ids, plot.mean20)
ionMean20 <- merge(variety_lookup, ionMean20, by = "plot")

selector <- substr(ionMean20$variety, 1, 3) == "E71"
ionRILMean20 <- ionMean20[selector, ]

is.control <- substr(ionRILMean20$plot, 2, 2) == "C"
ionRILMean20RC <- ionRILMean20[is.control,]
ionRILMean20RD <- ionRILMean20[!is.control,]

fileNameIonRILMean20RC <- paste0(saveDir, scriptID, "ionRILMean.20RC.csv")
write.csv(x = ionRILMean20RC, file = fileNameIonRILMean20RC)

fileNameIonRILMean20RD <- paste0(saveDir, scriptID, "ionRILMean.20RD.csv")
write.csv(x = ionRILMean20RD, file = fileNameIonRILMean20RD)




###### 2. PCA ######
##### 2.1. PCA #####
#### 2019RC ####
is.ion19RC <- (substr(colnames(ionRILMean19RC), 1, 1) != "v") &
  (substr(colnames(ionRILMean19RC), 1, 1) != "p")

dmat19RC <- as.matrix(ionRILMean19RC[,is.ion19RC])

dmat19RC[is.nan(dmat19RC)] <- NA

rownames(dmat19RC) <- ionRILMean19RC$variety

all_na_rows19RC <- rowSums(!is.na(dmat19RC)) == 0
na_row_names19RC <- rownames(dmat19RC)[all_na_rows19RC]
print(na_row_names19RC)

dmat19RCSel <- dmat19RC[rowSums(!is.na(dmat19RC)) > 0, ]
na_count <- sum(is.na(dmat19RCSel))
print(na_count)

pca19RCSel <- pca(dmat19RCSel, method = "ppca", nPcs = 12, scale = "uv", completeObs = TRUE)

#### 2019RD ####
is.ion19RD <- (substr(colnames(ionRILMean19RD), 1, 1) != "v") &
  (substr(colnames(ionRILMean19RD), 1, 1) != "p")

dmat19RD <- as.matrix(ionRILMean19RD[,is.ion19RD])

dmat19RD[is.nan(dmat19RD)] <- NA

rownames(dmat19RD) <- ionRILMean19RD$variety

all_na_rows19RD <- rowSums(!is.na(dmat19RD)) == 0
na_row_names19RD <- rownames(dmat19RD)[all_na_rows19RD]
print(na_row_names19RD)

dmat19RDSel <- dmat19RD[rowSums(!is.na(dmat19RD)) > 0, ]
na_count <- sum(is.na(dmat19RDSel))
print(na_count)

pca19RDSel <- pca(dmat19RDSel, method = "ppca", nPcs = 12, scale = "uv", completeObs = TRUE)

#### 2020RC ####
is.ion20RC <- (substr(colnames(ionRILMean20RC), 1, 1) != "v") &
  (substr(colnames(ionRILMean20RC), 1, 1) != "p")

dmat20RC <- as.matrix(ionRILMean20RC[,is.ion20RC])

dmat20RC[is.nan(dmat20RC)] <- NA

rownames(dmat20RC) <- ionRILMean20RC$variety

all_na_rows20RC <- rowSums(!is.na(dmat20RC)) == 0
na_row_names20RC <- rownames(dmat20RC)[all_na_rows20RC]
print(na_row_names20RC)

dmat20RCSel <- dmat20RC[rowSums(!is.na(dmat20RC)) > 0, ]
na_count <- sum(is.na(dmat20RCSel))
print(na_count)

pca20RCSel <- pca(dmat20RCSel, method = "ppca", nPcs = 12, scale = "uv", completeObs = TRUE)


#### 2020RD ####
is.ion20RD <- (substr(colnames(ionRILMean20RD), 1, 1) != "v") &
  (substr(colnames(ionRILMean20RD), 1, 1) != "p")

dmat20RD <- as.matrix(ionRILMean20RD[,is.ion20RD])

dmat20RD[is.nan(dmat20RD)] <- NA

rownames(dmat20RD) <- ionRILMean20RD$variety

all_na_rows20RD <- rowSums(!is.na(dmat20RD)) == 0
na_row_names20RD <- rownames(dmat20RD)[all_na_rows20RD]
print(na_row_names20RD)

dmat20RDSel <- dmat20RD[rowSums(!is.na(dmat20RD)) > 0, ]
na_count <- sum(is.na(dmat20RDSel))
print(na_count)

pca20RDSel <- pca(dmat20RDSel, method = "ppca", nPcs = 12, scale = "uv", completeObs = TRUE)

ion.pca <- list(RC2019 = pca19RCSel, RD2019 = pca19RDSel,
                RC2020 = pca20RCSel, RD2020 = pca20RDSel)
fileNameIonPCA <- paste0(saveDir, scriptID, "ionRILMean.PCA.RDS")
saveRDS(ion.pca, file = fileNameIonPCA)



