###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(stringr)
library(pbmcapply)


##### 0.2. Define Parameters #####
scriptID <- "9.32."

nIter <- 20000
burnIn <- 1000
thin <- 10

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
saveDir <- paste0("results/ionome/environment/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}
BGLRDir <- paste0("results/BGLR/", scriptID, "/")
if(!dir.exists(BGLRDir)) dir.create(BGLRDir, recursive = TRUE)




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameKm_unique <- paste0(dataDir, "9.2.Km.RDS")
Km_unique <- readRDS(fileNameKm_unique)

fileNameKm <- paste0(dataDir, "9.22.Km.all.RDS")
Km <- readRDS(file = fileNameKm)

fileNamePheno19RC <- paste0("data/", "Tottori_RIL.csv")
pheno19RC  <- read.csv(file = fileNamePheno19RC, row.names = 1)
pheno19RC  <- pheno19RC[, 2:16]

fileNamePheno19RD <- paste0("data/", "Tottori_RIL.csv")
pheno19RD <- read.csv(file = fileNamePheno19RD, row.names = 1)
pheno19RD <- pheno19RD[,18:32]
colnames(pheno19RD) <- gsub("\\.1$", "", colnames(pheno19RD))

fileNamePheno20RC <- paste0("data/", "Tottori.2020.RIL.csv")
pheno20RC <- read.csv(file = fileNamePheno20RC, row.names = 1)

pheno20RC <- pheno20RC[pheno20RC$block == "RC", ]
pheno20RC$block <- NULL   

rownames(pheno20RC) <- pheno20RC$variety
pheno20RC$variety <- NULL

fileNamePheno20RD <- paste0("data/", "Tottori.2020.RIL.csv")
pheno20RD <- read.csv(file = fileNamePheno20RD, row.names = 1)

pheno20RD <- pheno20RD[pheno20RD$block == "RD", ]
pheno20RD$block <- NULL   

rownames(pheno20RD) <- pheno20RD$variety
pheno20RD$variety <- NULL

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

fileNameKi <- paste0(dataDir, "9.23.Ki.all.RDS")
Ki <- readRDS(file = fileNameKi)


##### 1.2. Modify Data #####
commonIDs <- Reduce(intersect, list(
  rownames(Km_unique),
  rownames(pheno19RC),
  rownames(pheno19RD),
  rownames(pheno20RC),
  rownames(pheno20RD),
  rownames(dm.ionRC19),
  rownames(dm.ionRD19),
  rownames(dm.ionRC20),
  rownames(dm.ionRD20)
))

pheno19RC <- pheno19RC[commonIDs, , drop = FALSE]
pheno19RD <- pheno19RD[commonIDs, , drop = FALSE]
pheno20RC <- pheno20RC[commonIDs, , drop = FALSE]
pheno20RD <- pheno20RD[commonIDs, , drop = FALSE]

dm.ionRC19 <- dm.ionRC19[commonIDs, , drop = FALSE]
dm.ionRD19 <- dm.ionRD19[commonIDs, , drop = FALSE]
dm.ionRC20 <- dm.ionRC20[commonIDs, , drop = FALSE]
dm.ionRD20 <- dm.ionRD20[commonIDs, , drop = FALSE]

commonCols <- Reduce(intersect, list(
  colnames(pheno19RC),
  colnames(pheno19RD),
  colnames(pheno20RC),
  colnames(pheno20RD)
))

pheno19RC <- pheno19RC[, commonCols, drop = FALSE]
pheno19RD <- pheno19RD[, commonCols, drop = FALSE]
pheno20RC <- pheno20RC[, commonCols, drop = FALSE]
pheno20RD <- pheno20RD[, commonCols, drop = FALSE]

rownames(pheno19RC) <- paste0(rownames(pheno19RC), "_19RC")
rownames(pheno19RD) <- paste0(rownames(pheno19RD), "_19RD")
rownames(pheno20RC) <- paste0(rownames(pheno20RC), "_20RC")
rownames(pheno20RD) <- paste0(rownames(pheno20RD), "_20RD")

phenoAll <- rbind(
  pheno19RC, 
  pheno19RD, 
  pheno20RC, 
  pheno20RD
)

rownames(dm.ionRC19) <- paste0(rownames(dm.ionRC19), "_19RC")
rownames(dm.ionRD19) <- paste0(rownames(dm.ionRD19), "_19RD")
rownames(dm.ionRC20) <- paste0(rownames(dm.ionRC20), "_20RC")
rownames(dm.ionRD20) <- paste0(rownames(dm.ionRD20), "_20RD")

ionAll <- rbind(
  dm.ionRC19, 
  dm.ionRD19, 
  dm.ionRC20, 
  dm.ionRD20
)

Km <- Km[rownames(phenoAll), rownames(phenoAll)]

Ki <- Ki[rownames(phenoAll), rownames(phenoAll)]




###### 2.  ######
##### 2.2. Model 2 #####
year <- factor(str_sub(rownames(phenoAll), -4, -3))
treatment <- factor(str_sub(rownames(phenoAll), -2, -1))
phenoCol <- colnames(phenoAll)

Years  <- c("19", "20")
Blocks <- c("RC", "RD")

allResults <- list()

for (Year in Years) {
  for (Block in Blocks) {
    
    ind_train  <- which(!(year == Year & treatment == Block)) 
    ind_target <- which(year == Year & treatment == Block)
    
    pheno_train  <- phenoAll[ind_train, , drop = FALSE]
    pheno_target <- phenoAll[ind_target, , drop = FALSE]

    lineNames <- str_sub(rownames(phenoAll), 1, nchar(rownames(phenoAll)) - 5)
    uniqueLines <- unique(lineNames)

    cvResults <- list()

    for (phenoCC in phenoCol) {
      phenoNow <- phenoAll[[phenoCC]]
      names(phenoNow) <- rownames(phenoAll)

      seedResults <- list()
      
      for (seed in 1:10) {
        set.seed(seedSet[seed])
        
        folds <- sample(rep(1:10, length.out = length(uniqueLines)))
        names(folds) <- uniqueLines
        
        yPredVec <- rep(NA, length(phenoNow))
        
        cat(
          "\nYear:", Year,
          "Block:", Block,
          "Trait:", phenoCC,
          "Seed:", seed,
          "\n"
        )
        
        cvRes <- pbmclapply(1:10, function(i) {
          testLines <- uniqueLines[folds == i]
          
          trainInd_fold <- ind_train[!(lineNames[ind_train] %in% testLines)]
          
          y_train <- phenoNow
          y_train[-c(trainInd_fold)] <- NA  

          eta <- list(
            list(K = Km, model = "RKHS"),
            list(K = Ki, model = "RKHS")
          )
          
          model <- BGLR(
            y = y_train,
            ETA = eta,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            verbose = FALSE,
            saveAt = paste0(
              BGLRDir,
              "model2_", Year, Block, "_",
              phenoCC, "_seed", seed, "_fold", i, "_")
          )
          
          u_hat_Km <- model$ETA[[1]]$u 
          u_hat_Ki <- model$ETA[[2]]$u  
          yPred <- u_hat_Km + u_hat_Ki
          names(yPred) <- rownames(phenoAll)
          
          return(yPred)
        },
        mc.cores = mc.cores)
        
        yPredVec <- rep(NA, length(phenoNow))
        names(yPredVec) <- names(phenoNow)
        for (i in 1:10) {
          testLines <- uniqueLines[folds == i]
          testInd <- paste0(testLines, "_", Year, Block)
          yPredVec[testInd] <- cvRes[[i]][testInd]
        }
        
        pheno_targetNow <- pheno_target[[phenoCC]]
        names(pheno_targetNow ) <- rownames(pheno_target)
        
        plot(pheno_targetNow, yPredVec[names(pheno_targetNow)])
        cor(pheno_targetNow, yPredVec[names(pheno_targetNow)])
        
        common_labels <- intersect(names(pheno_targetNow), names(yPredVec))
        
        yPred_aligned   <- yPredVec[common_labels]
        pheno_aligned   <- pheno_targetNow[common_labels]
        
        correlation <- cor(yPred_aligned, pheno_aligned , use = "complete.obs")
        
        seedResults[[paste0("seed_", seed)]] <- list(
          predictions = yPredVec,
          correlation = correlation
        )
      }

      cvResults[[phenoCC]] <- seedResults
    }
    
    corAve <- sapply(cvResults, function(seedResults) {
      seedCorrelations <- sapply(seedResults, function(x) x$correlation)
      mean(seedCorrelations, na.rm = TRUE)
    })

    allCorrelations <- unlist(lapply(cvResults, function(seedResults) {
      sapply(seedResults, function(x) x$correlation)
    }))
    
    corAve["All"] <- mean(allCorrelations, na.rm = TRUE)
    
    fileNameCorAve <- paste0(saveDir, scriptID, "model2.", Year, Block, ".corAve.RDS")
    saveRDS(corAve, file = fileNameCorAve)
    
    allResults[[paste0(Year, Block)]] <- corAve
  }
}

fileNameallResults <- paste0(saveDir, scriptID, "model2.all.corAve.RDS")
saveRDS(allResults, file = fileNameallResults)


##### 2.4. Model 4 #####
year <- factor(str_sub(rownames(phenoAll), -4, -3))
treatment <- factor(str_sub(rownames(phenoAll), -2, -1))
phenoCol <- colnames(phenoAll)

yearTreatment <- factor(paste0(year, "_", treatment))

Years  <- c("19", "20")
Blocks <- c("RC", "RD")

allResults <- list()

for (Year in Years) {
  for (Block in Blocks) {
    
    ind_train  <- which(!(year == Year & treatment == Block)) 
    ind_target <- which(year == Year & treatment == Block)
    
    pheno_train  <- phenoAll[ind_train, , drop = FALSE]
    pheno_target <- phenoAll[ind_target, , drop = FALSE]
    
    lineNames <- str_sub(rownames(phenoAll), 1, nchar(rownames(phenoAll)) - 5)
    uniqueLines <- unique(lineNames)
    
    cvResults <- list()
    
    yPredResults <- list()
    
    for (phenoCC in phenoCol) {
      phenoNow <- phenoAll[[phenoCC]]
      names(phenoNow) <- rownames(phenoAll)
      
      seedResults <- list()
      
      for (seed in 1:10) {
        set.seed(seedSet[seed])
        
        folds <- sample(rep(1:10, length.out = length(uniqueLines)))
        names(folds) <- uniqueLines
        
        yPredVec <- rep(NA, length(phenoNow))
        
        cat(
          "\nYear:", Year,
          "Block:", Block,
          "Trait:", phenoCC,
          "Seed:", seed,
          "\n"
        )
        
        cvRes <- pbmclapply(1:10, function(i) {
          testLines <- uniqueLines[folds == i]
          
          trainInd_fold <- ind_train[!(lineNames[ind_train] %in% testLines)]
          
          y_train <- phenoNow
          y_train[-c(trainInd_fold)] <- NA 
          
          eta <- list(
            list(K = Km, model = "RKHS"),
            list(X = model.matrix(~ yearTreatment - 1), model = "FIXED")
          )
          
          model <- BGLR(
            y = y_train,
            ETA = eta,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            verbose = FALSE,
            saveAt = paste0(
              BGLRDir,
              "model4_", Year, Block, "_",
              phenoCC, "_seed", seed, "_fold", i, "_")
          )

          u_hat_Km <- model$ETA[[1]]$u 
        
          yPred <- u_hat_Km 
          names(yPred) <- rownames(phenoAll)
          
          return(yPred)
        },
        mc.cores = mc.cores)
        
        yPredVec <- rep(NA, length(phenoNow))
        names(yPredVec) <- names(phenoNow)
        for (i in 1:10) {
          testLines <- uniqueLines[folds == i]
          testInd <- paste0(testLines, "_", Year, Block)
          yPredVec[testInd] <- cvRes[[i]][testInd]
        }
        
        pheno_targetNow <- pheno_target[[phenoCC]]
        names(pheno_targetNow ) <- rownames(pheno_target)
        
        yPred_target <- yPredVec[names(pheno_targetNow)]
        if (is.null(yPredResults[[phenoCC]])) {
          yPredResults[[phenoCC]] <- list()
        }
        yPredResults[[phenoCC]][[paste0("seed_", seed)]] <- yPred_target
        
        plot(pheno_targetNow, yPredVec[names(pheno_targetNow)])
        cor(pheno_targetNow, yPredVec[names(pheno_targetNow)])
        
        common_labels <- intersect(names(pheno_targetNow), names(yPredVec))
        
        yPred_aligned   <- yPredVec[common_labels]
        pheno_aligned   <- pheno_targetNow[common_labels]
        
        correlation <- cor(yPred_aligned, pheno_aligned , use = "complete.obs")

        seedResults[[paste0("seed_", seed)]] <- list(
          predictions = yPredVec,
          correlation = correlation
        )
      }

      cvResults[[phenoCC]] <- seedResults
    }

    corAve <- sapply(cvResults, function(seedResults) {
      seedCorrelations <- sapply(seedResults, function(x) x$correlation)
      mean(seedCorrelations, na.rm = TRUE)
    })

    allCorrelations <- unlist(lapply(cvResults, function(seedResults) {
      sapply(seedResults, function(x) x$correlation)
    }))
    
    corAve["All"] <- mean(allCorrelations, na.rm = TRUE)

    fileNameCorAve <- paste0(saveDir, scriptID, "model4.", Year, Block, ".corAve.RDS")
    saveRDS(corAve, file = fileNameCorAve)

    fileNameyPred <- paste0(saveDir, scriptID, "model4.", Year, Block, ".u_pred.RDS")
    saveRDS(yPredResults, file = fileNameyPred)

    allResults[[paste0(Year, Block)]] <- corAve
  }
}

fileNameallResults <- paste0(saveDir, scriptID, "model4.all.corAve.RDS")
saveRDS(allResults, file = fileNameallResults)


##### 2.5. Model 5 #####
year <- factor(str_sub(rownames(phenoAll), -4, -3))
treatment <- factor(str_sub(rownames(phenoAll), -2, -1))
phenoCol <- colnames(phenoAll)

yearTreatment <- factor(paste0(year, "_", treatment))

Years  <- c("19", "20")
Blocks <- c("RC", "RD")

allResults <- list()

for (Year in Years) {
  for (Block in Blocks) {
    
    ind_train  <- which(!(year == Year & treatment == Block)) 
    ind_target <- which(year == Year & treatment == Block)
    
    pheno_train  <- phenoAll[ind_train, , drop = FALSE]
    pheno_target <- phenoAll[ind_target, , drop = FALSE]
    
    lineNames <- str_sub(rownames(phenoAll), 1, nchar(rownames(phenoAll)) - 5)
    uniqueLines <- unique(lineNames)

    cvResults <- list()
    
    for (phenoCC in phenoCol) {
      phenoNow <- phenoAll[[phenoCC]]
      names(phenoNow) <- rownames(phenoAll)
      
      seedResults <- list()

      for (seed in 1:10) {
        set.seed(seedSet[seed])
        
        folds <- sample(rep(1:10, length.out = length(uniqueLines)))
        names(folds) <- uniqueLines
        
        yPredVec <- rep(NA, length(phenoNow))
        
        cat(
          "\nYear:", Year,
          "Block:", Block,
          "Trait:", phenoCC,
          "Seed:", seed,
          "\n"
        )

        cvRes <- pbmclapply(1:10, function(i) {
          testLines <- uniqueLines[folds == i]

          trainInd_fold <- ind_train[!(lineNames[ind_train] %in% testLines)]
          
          y_train <- phenoNow
          y_train[-c(trainInd_fold)] <- NA 

          eta <- list(
            list(K = Km, model = "RKHS"),
            list(K = Ki, model = "RKHS"),
            list(X = model.matrix(~ yearTreatment - 1), model = "FIXED")
          )

          model <- BGLR(
            y = y_train,
            ETA = eta,
            nIter = nIter,
            burnIn = burnIn,
            thin = thin,
            verbose = FALSE,
            saveAt = paste0(
              BGLRDir,
              "model5_", Year, Block, "_",
              phenoCC, "_seed", seed, "_fold", i, "_")
          )
          
          u_hat_Km <- model$ETA[[1]]$u 
          u_hat_Ki <- model$ETA[[2]]$u  
          
          yPred <- u_hat_Km + u_hat_Ki
          names(yPred) <- rownames(phenoAll)
          
          return(yPred)
        },
        mc.cores = mc.cores)

        yPredVec <- rep(NA, length(phenoNow))
        names(yPredVec) <- names(phenoNow)
        for (i in 1:10) {
          testLines <- uniqueLines[folds == i]
          testInd <- paste0(testLines, "_", Year, Block)
          yPredVec[testInd] <- cvRes[[i]][testInd]
        }

        pheno_targetNow <- pheno_target[[phenoCC]]
        names(pheno_targetNow ) <- rownames(pheno_target)
        
        plot(pheno_targetNow, yPredVec[names(pheno_targetNow)])
        cor(pheno_targetNow, yPredVec[names(pheno_targetNow)])

        common_labels <- intersect(names(pheno_targetNow), names(yPredVec))
        
        yPred_aligned   <- yPredVec[common_labels]
        pheno_aligned   <- pheno_targetNow[common_labels]

        correlation <- cor(yPred_aligned, pheno_aligned , use = "complete.obs")

        seedResults[[paste0("seed_", seed)]] <- list(
          predictions = yPredVec,
          correlation = correlation
        )
      }

      cvResults[[phenoCC]] <- seedResults
    }

    corAve <- sapply(cvResults, function(seedResults) {
      seedCorrelations <- sapply(seedResults, function(x) x$correlation)
      mean(seedCorrelations, na.rm = TRUE)
    })

    allCorrelations <- unlist(lapply(cvResults, function(seedResults) {
      sapply(seedResults, function(x) x$correlation)
    }))
    
    corAve["All"] <- mean(allCorrelations, na.rm = TRUE)

    fileNameCorAve <- paste0(saveDir, scriptID, "model5.", Year, Block, ".corAve.RDS")
    saveRDS(corAve, file = fileNameCorAve)

    allResults[[paste0(Year, Block)]] <- corAve
  }
}

fileNameallResults <- paste0(saveDir, scriptID, "model5.all.corAve.RDS")
saveRDS(allResults, file = fileNameallResults)



