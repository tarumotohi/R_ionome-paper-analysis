###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(pbmcapply)


##### 0.2. Define Parameters #####
scriptID <- "9.7_2."

nIter <- 20000
burnIn <- 1000

mc.cores <- 6

yearBlock <- "20RD"


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- paste0("results/ionome/5Models/", yearBlock ,"/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}
BGLRDir <- paste0("results/BGLR/", scriptID, "/")
if(!dir.exists(BGLRDir)) dir.create(BGLRDir, recursive = TRUE)




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameKm <- paste0(dataDir, "9.2.Km.RDS")
Km <- readRDS(file = fileNameKm)

fileNameKig <- paste0(dataDir, "9.1.Kig.RD20.2.RDS")
Kig <- readRDS(file = fileNameKig)

fileNameKi <- paste0(dataDir, "9.1.Ki.RD20.RDS")
Ki <- readRDS(file = fileNameKi)

fileNamePheno <- paste0("data/", "Tottori.2020.RIL.csv")
pheno <- read.csv(file = fileNamePheno, row.names = 1)

phenoSel <- pheno[pheno$block == "RD", ]
phenoSel$block <- NULL   

rownames(phenoSel) <- phenoSel$variety
phenoSel$variety <- NULL

commonIDs <- Reduce(intersect, list(rownames(phenoSel), rownames(Km), rownames(Kig)))

phenoSel <- phenoSel[commonIDs, , drop = FALSE]
Km <- Km[commonIDs, commonIDs]
Kig <- Kig[commonIDs, commonIDs]
Ki <- Ki[commonIDs, commonIDs]



###### 2. GP ######
##### 2.1. Model 1 #####
phenoCol <- colnames(phenoSel)

cvResults <- list()

for (phenoCC in phenoCol) {
  phenoCCNow <- phenoSel[[phenoCC]]
  names(phenoCCNow) <- rownames(phenoSel)
  
  seedResults <- list()
  
  for (seed in 1:10) {
    set.seed(seed)
    folds <- sample(rep(1:10, length.out = length(phenoCCNow)))
    
    cvRes <- pbmclapply(1:10, function(i) {
      testInd <- which(folds == i)
      trainInd <- setdiff(1:length(phenoCCNow), testInd)
      
      y <- phenoCCNow
      y[testInd] <- NA
      
      eta <- list(list(K = Km, model = "RKHS"))
      
      model1 <- BGLR(
        y = y,
        ETA = eta,
        nIter = nIter,
        burnIn = burnIn,
        verbose = FALSE,
        saveAt = paste0(
          BGLRDir,
          "model1_", yearBlock, "_",
          phenoCC, "_seed", seed, "_fold", i, "_")
      )
      
      uHat <- model1$ETA[[1]]$u
      names(uHat) <- rownames(Km)
      uPred <- uHat[testInd]
      
      return(uPred)
    },
    mc.cores = mc.cores)
    
    yPredVec <- rep(NA, length(phenoCCNow))
    for (i in 1:10) {
      testInd <- which(folds == i)
      yPredVec[testInd] <- cvRes[[i]]
    }
    correlation <- cor(yPredVec, phenoCCNow, use = "complete.obs")
    
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

fileNameCorAve <- paste0(saveDir, scriptID, yearBlock, ".model1.corAve.RDS")
saveRDS(corAve, file = fileNameCorAve)


##### 2.1.5 Model 1.5 #####
phenoCol <- colnames(phenoSel)

cvResults <- list()

for (phenoCC in phenoCol) {
  phenoCCNow <- phenoSel[[phenoCC]]
  names(phenoCCNow) <- rownames(phenoSel)
  
  seedResults <- list()
  
  for (seed in 1:10) {
    set.seed(seed)
    folds <- sample(rep(1:10, length.out = length(phenoCCNow)))
    
    cvRes <- pbmclapply(1:10, function(i) {
      testInd <- which(folds == i)
      trainInd <- setdiff(1:length(phenoCCNow), testInd)
      
      y <- phenoCCNow
      y[testInd] <- NA
      
      eta <- list(list(K = Ki, model = "RKHS"))
      
      model1.5 <- BGLR(
        y = y,
        ETA = eta,
        nIter = nIter,
        burnIn = burnIn,
        verbose = FALSE,
        saveAt = paste0(
          BGLRDir,
          "model1.5_", yearBlock, "_",
          phenoCC, "_seed", seed, "_fold", i, "_")
      )
      
      uHat <- model1.5$ETA[[1]]$u
      names(uHat) <- rownames(Kig)
      uPred <- uHat[testInd]
      
      return(uPred)
    },
    mc.cores = mc.cores)
    
    yPredVec <- rep(NA, length(phenoCCNow))
    for (i in 1:10) {
      testInd <- which(folds == i)
      yPredVec[testInd] <- cvRes[[i]]
    }
    correlation <- cor(yPredVec, phenoCCNow, use = "complete.obs")
    
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

fileNameCorAve <- paste0(saveDir, scriptID, yearBlock, ".model1.5.corAve.RDS")
saveRDS(corAve, file = fileNameCorAve)


##### 2.2 Model 2 #####
phenoCol <- colnames(phenoSel)

cvResults <- list()

for (phenoCC in phenoCol) {
  phenoCCNow <- phenoSel[[phenoCC]]
  names(phenoCCNow) <- rownames(phenoSel)
  
  seedResults <- list()
  
  for (seed in 1:10) {
    set.seed(seed)
    folds <- sample(rep(1:10, length.out = length(phenoCCNow)))
    
    cvRes <- pbmclapply(1:10, function(i) {
      testInd <- which(folds == i)
      trainInd <- setdiff(1:length(phenoCCNow), testInd)
      
      y <- phenoCCNow
      y[testInd] <- NA
      
      eta <- list(list(K = Kig, model = "RKHS"))
      
      model2 <- BGLR(
        y = y,
        ETA = eta,
        nIter = nIter,
        burnIn = burnIn,
        verbose = FALSE,
        saveAt = paste0(
          BGLRDir,
          "model2_", yearBlock, "_",
          phenoCC, "_seed", seed, "_fold", i, "_")
      )
      
      uHat <- model2$ETA[[1]]$u
      names(uHat) <- rownames(Kig)
      uPred <- uHat[testInd]
      
      return(uPred)
    },
    mc.cores = mc.cores)
    
    yPredVec <- rep(NA, length(phenoCCNow))
    for (i in 1:10) {
      testInd <- which(folds == i)
      yPredVec[testInd] <- cvRes[[i]]
    }
    correlation <- cor(yPredVec, phenoCCNow, use = "complete.obs")
    
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

fileNameCorAve <- paste0(saveDir, scriptID, yearBlock, ".model2.corAve.RDS")
saveRDS(corAve, file = fileNameCorAve)


##### 2.2.5 Model 2.5 #####
phenoCol <- colnames(phenoSel)

cvResults <- list()

for (phenoCC in phenoCol) {
  phenoCCNow <- phenoSel[[phenoCC]]
  names(phenoCCNow) <- rownames(phenoSel)
  
  seedResults <- list()
  
  for (seed in 1:10) {
    set.seed(seed)
    folds <- sample(rep(1:10, length.out = length(phenoCCNow)))
    
    cvRes <- pbmclapply(1:10, function(i) {
      testInd <- which(folds == i)
      trainInd <- setdiff(1:length(phenoCCNow), testInd)
      
      y <- phenoCCNow
      y[testInd] <- NA
      
      eta <- list(list(K = Km, model = "RKHS"),
                  list(K = Ki, model = "RKHS"))
      
      model2.5 <- BGLR(
        y = y,
        ETA = eta,
        nIter = nIter,
        burnIn = burnIn,
        verbose = FALSE,
        saveAt = paste0(
          BGLRDir,
          "model2.5_", yearBlock, "_",
          phenoCC, "_seed", seed, "_fold", i, "_")
      )
      
      u_hat_Km <- model2.5$ETA[[1]]$u  
      u_hat_Kig <- model2.5$ETA[[2]]$u 
      
      uHat <- u_hat_Km + u_hat_Kig
      names(uHat) <- rownames(Km)
      uPred <- uHat[testInd]
      
      return(uPred)
    },
    mc.cores = mc.cores)
    
    yPredVec <- rep(NA, length(phenoCCNow))
    for (i in 1:10) {
      testInd <- which(folds == i)
      yPredVec[testInd] <- cvRes[[i]]
    }
    correlation <- cor(yPredVec, phenoCCNow, use = "complete.obs")
    
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

fileNameCorAve <- paste0(saveDir, scriptID, yearBlock, ".model2.5.corAve.RDS")
saveRDS(corAve, file = fileNameCorAve)


##### 2.3 Model 3 #####
phenoCol <- colnames(phenoSel)

cvResults <- list()

for (phenoCC in phenoCol) {
  phenoCCNow <- phenoSel[[phenoCC]]
  names(phenoCCNow) <- rownames(phenoSel)
  
  seedResults <- list()
  
  for (seed in 1:10) {
    set.seed(seed)
    folds <- sample(rep(1:10, length.out = length(phenoCCNow)))
    
    cvRes <- pbmclapply(1:10, function(i) {
      testInd <- which(folds == i)
      trainInd <- setdiff(1:length(phenoCCNow), testInd)
      
      y <- phenoCCNow
      y[testInd] <- NA
      
      eta <- list(list(K = Km, model = "RKHS"),
                  list(K = Kig, model = "RKHS"))
      
      model3 <- BGLR(
        y = y,
        ETA = eta,
        nIter = nIter,
        burnIn = burnIn,
        verbose = FALSE,
        saveAt = paste0(
          BGLRDir,
          "model3_", yearBlock, "_",
          phenoCC, "_seed", seed, "_fold", i, "_")
      )
      
      u_hat_Km <- model3$ETA[[1]]$u  
      u_hat_Kig <- model3$ETA[[2]]$u  
      
      uHat <- u_hat_Km + u_hat_Kig
      names(uHat) <- rownames(Km)
      uPred <- uHat[testInd]
      
      return(uPred)
    },
    mc.cores = mc.cores)
    
    yPredVec <- rep(NA, length(phenoCCNow))
    for (i in 1:10) {
      testInd <- which(folds == i)
      yPredVec[testInd] <- cvRes[[i]]
    }
    correlation <- cor(yPredVec, phenoCCNow, use = "complete.obs")
    
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

fileNameCorAve <- paste0(saveDir, scriptID, yearBlock, ".model3.corAve.RDS")
saveRDS(corAve, file = fileNameCorAve)


##### 2.4 5 Models #####
# ファイルパス
fileModel1 <- paste0(saveDir, scriptID, yearBlock, ".model1.corAve.RDS")
fileModel1.5 <- paste0(saveDir, scriptID, yearBlock, ".model1.5.corAve.RDS")
fileModel2 <- paste0(saveDir, scriptID, yearBlock, ".model2.corAve.RDS")
fileModel2.5 <- paste0(saveDir, scriptID, yearBlock, ".model2.5.corAve.RDS")
fileModel3 <- paste0(saveDir, scriptID, yearBlock, ".model3.corAve.RDS")

corModel1 <- readRDS(fileModel1)
corModel1.5 <- readRDS(fileModel1.5)
corModel2 <- readRDS(fileModel2)
corModel2.5 <- readRDS(fileModel2.5)
corModel3 <- readRDS(fileModel3)

traitNames <- Reduce(intersect, list(
  names(corModel1),
  names(corModel1.5),
  names(corModel2),
  names(corModel2.5),
  names(corModel3)
))

corSummary <- data.frame(
  Trait = traitNames,
  Model1 = unname(corModel1[traitNames]),
  Model1.5 = unname(corModel1.5[traitNames]),
  Model2 = unname(corModel2[traitNames]),
  Model2.5 = unname(corModel2.5[traitNames]),
  Model3 = unname(corModel3[traitNames])
)

if ("All" %in% names(corModel1)) {
  corSummary <- rbind(
    corSummary,
    data.frame(Trait = "All",
               Model1 = corModel1["All"],
               Model1.5 = corModel1.5["All"],
               Model2 = corModel2["All"],
               Model2.5 = corModel2.5["All"],
               Model3 = corModel3["All"])
  )
}

fileSummary <- paste0(saveDir, scriptID, yearBlock, ".corAve.allModels.RDS")
saveRDS(corSummary, file = fileSummary)

fileSummaryCsv <- paste0(saveDir, scriptID, yearBlock, ".corAve.allModels.csv")
write.csv(corSummary, file = fileSummaryCsv, row.names = FALSE)



