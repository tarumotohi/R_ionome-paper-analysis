###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.45."

nIter <- 20000
burnIn <- 1000


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

fileNameKi <- paste0(dataDir, "9.40.Ki.RIL_F2.RDS")
Ki <- readRDS(file = fileNameKi)

#F2
fileNameShoot <- "data/shoot.csv"
shoot <- read.csv(file = fileNameShoot, row.names = 1)

fileNameYield <- "data/yield.csv"
yield <- read.csv(file = fileNameYield, row.names = 1)

shoot_F2 <- shoot[shoot$year == 2021 & nchar(shoot$variety) == 24,]
yield_F2 <- yield[yield$year == 2021 & nchar(yield$variety) == 24,]

yield_F2 <- mutate(yield_F2, SeedWeight_10_g
                   = round((SeedWeight_all_g / NumOfSeed) * 10, 3))

rownames(shoot_F2) <- shoot_F2$variety

yield_F2 <- yield_F2[!(is.na(yield_F2[,10]) & is.na(yield_F2[,11])), ]
rownames(yield_F2) <- yield_F2$variety

shoot_F2 <- shoot_F2[, -c(1:3, 5:9)]
yield_F2 <- yield_F2[, -c(1:3, 5:9)]

phenoF2 <- bind_rows(shoot_F2, yield_F2)

rownames(phenoF2) <- phenoF2$variety
phenoF2 <- phenoF2[, -c(1,3)]

colnames(phenoF2) <- gsub("FreshWeight_Shoot_g", "FSW", colnames(phenoF2))
colnames(phenoF2) <- gsub("GrowthStage", "GS", colnames(phenoF2))
colnames(phenoF2) <- gsub("NumberOfTillers", "NB", colnames(phenoF2))
colnames(phenoF2) <- gsub("NumberOfNodes", "NN", colnames(phenoF2))
colnames(phenoF2) <- gsub("TillerLength_cm", "BL", colnames(phenoF2))
colnames(phenoF2) <- gsub("SeedWeight_10_g", "TENSW", colnames(phenoF2))
colnames(phenoF2) <- gsub("SeedWeight_all_g", "TSW", colnames(phenoF2))
colnames(phenoF2) <- gsub("NumOfSeed", "NS", colnames(phenoF2))

phenoAll <-  phenoF2

dim(Km)
dim(Ki)
dim(phenoAll)

common_ids <- intersect(rownames(Ki), rownames(phenoAll))
Ki <- Ki[common_ids, common_ids, drop=FALSE]
phenoAll <- phenoAll[common_ids, , drop=FALSE]




###### 2.  ######
##### 2.1. Model 1 y = u + e #####
phenoCol <- colnames(phenoAll)

Blocks <- c("RC", "RD")

KmCountResults <- list()

allResults <- list()

for (Block in Blocks) {
  
  pheno_block <- phenoAll[substr(rownames(phenoAll), 17, 18) == Block, ]
  
  u_fitResults <- list()
  hg2Results <- list()
  KmCount <- list()
  
  for (phenoCC in phenoCol) {
    
    pheno_block_trait <- pheno_block[[phenoCC]]
    names(pheno_block_trait) <- rownames(pheno_block)
    
    pheno_block_trait_sel <- pheno_block_trait[!is.na(pheno_block_trait)]
    
    ind <- intersect(
      names(pheno_block_trait_sel),
      rownames(Km)
    )
    
    print(length(ind))
    KmCount[[phenoCC]] <- length(ind)
    
    y <- pheno_block_trait_sel[ind]
    
    eta <- list(
      list(K = Km[ind, ind], model = "RKHS")
    )
    
    modelDir <- paste0(BGLRDir, "model1_", Block, "_", phenoCC, "/")
    if(!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)
    
    model <- BGLR(
      y = y,
      ETA = eta,
      nIter = nIter,
      burnIn = burnIn,
      verbose = FALSE,
      saveAt = modelDir
    )
    
    u_hat_Km <- model$ETA[[1]]$u
    names(u_hat_Km) <- ind
    
    u_fitResults[[phenoCC]] <- u_hat_Km
    
    sigma_g2 <- model$ETA[[1]]$varU
    sigma_e2 <- model$varE
    
    hg2Results[[phenoCC]] <- sigma_g2 / (sigma_g2 + sigma_e2)
  }
  
  allResults[[Block]] <- u_fitResults
  
  KmCountResults[[Block]] <- KmCount
  
  fileNameuhat <- paste0(saveDir, scriptID, "model1.", Block, ".u_fit.RDS")
  saveRDS(u_fitResults, file = fileNameuhat)
  
  fileNamehg2 <- paste0(saveDir, scriptID, "model1.", Block, ".hg2.RDS")
  saveRDS(hg2Results, file = fileNamehg2)
}



