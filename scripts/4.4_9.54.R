###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(stringr)
library(pbmcapply)
library(ggplot2)
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.54."

nIter <- 20000
burnIn <- 1000

mc.cores <- 6

seedFile <- "data/seedSet.rds"
if(file.exists(seedFile)) {
  seedSet <- readRDS(seedFile)
} else {
  seedSet <- sample(x = 1:1e4, size = 100)
  saveRDS(seedSet, file = seedFile)
}


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- paste0("results/ionome/RIL_F2/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}
BGLRDir <- paste0("results/BGLR/", scriptID, "/")
if(!dir.exists(BGLRDir)) dir.create(BGLRDir, recursive = TRUE)




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameKm <- paste0(dataDir, "9.39.Km.RDS")
Km <- readRDS(fileNameKm)

rows7 <- rownames(Km)[nchar(rownames(Km)) == 7]

suffixes <- c("_19RD","_19RC","_20RD","_20RC")

Km_new <- Km[rows7, rows7, drop=FALSE]

for(r in rows7){
  for(s in suffixes){
    new_name <- paste0(r, s)

    new_row <- Km_new[r, , drop=FALSE]

    Km_new <- rbind(Km_new, new_row)
    rownames(Km_new)[nrow(Km_new)] <- new_name

    new_col <- Km_new[, r, drop=FALSE]
    Km_new <- cbind(Km_new, new_col)
    colnames(Km_new)[ncol(Km_new)] <- new_name
  }
}

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

common_ids <- intersect(rownames(Km_new), rownames(dm.ionome_RIL))
Km_new <- Km_new[common_ids, common_ids, drop=FALSE] 
dm.ionome_RIL <- dm.ionome_RIL[common_ids, , drop=FALSE]




###### 2.  ######
##### 2.1. y(ion) = Y + u + e #####
treatments <- c("RD", "RC")
traits <- colnames(dm.ionome_RIL)

results1 <- list()

for(i in seq_along(treatments)) {
  
  treatment <- treatments[i]
  
  dm.ionome_RIL_treat <- dm.ionome_RIL[substr(rownames(dm.ionome_RIL), 11, 12) == treatment, ]
  
  common_RIL <- intersect(rownames(dm.ionome_RIL_treat), rownames(Km_new))
  Km_RIL <- Km_new[common_RIL, common_RIL]
  dm.ionome_RIL_treat <- dm.ionome_RIL_treat[common_RIL, , drop=FALSE] 
  
  print(all(rownames(dm.ionome_RIL_treat) == rownames(Km_RIL)))

  year <- substr(rownames(dm.ionome_RIL_treat), 9, 10)
  year <- factor(year) 

  eta <- list(
    list(X = model.matrix(~ year-1), model = "FIXED"),
    list(K = Km_RIL, model = "RKHS")
  )
  
  set.seed(seedSet[i])
  
  results1[[treatment]] <- list()

  e_hatResults <- list()
  hg2Results <- list()
  
  for(trait in traits) {
    
    modelDir <- paste0(BGLRDir, "model1_", treatment, "_", trait, "/")
    if(!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)

    fit <- BGLR(
      y = dm.ionome_RIL_treat[[trait]],
      ETA = eta,
      nIter = nIter,
      burnIn = burnIn,
      verbose = FALSE,
      saveAt = modelDir
    )

    u_hat_Km <- fit$ETA[[2]]$u  

    sigma_g2 <- fit$ETA[[2]]$varU
    sigma_e2 <- fit$varE
    
    hg2Results[[trait]] <- sigma_g2 / (sigma_g2 + sigma_e2)
    
    y_vec <- dm.ionome_RIL_treat[[trait]]
    
    res <- y_vec - fit$yHat
    
    names(res) <- rownames(dm.ionome_RIL_treat)

    e_hatResults[[trait]] <- res
  }

  fileNamehg2 <- paste0(saveDir, scriptID, "model1.", treatment, ".hg2.RDS")
  saveRDS(hg2Results, file = fileNamehg2)
  
  fileName_ehat <- paste0(saveDir, scriptID, "model1.", treatment, ".ehat.RDS")
  saveRDS(e_hatResults, file = fileName_ehat)
  
}



