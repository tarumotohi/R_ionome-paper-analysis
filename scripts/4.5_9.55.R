###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(stringr)
library(ggplot2)
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.55."

nIter <- 20000
burnIn <- 1000

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

fileNameIonome_ehat_RC <- paste0(saveDir, "9.54.model1.RC.ehat.RDS")
Ionome_ehat_RC <- readRDS(fileNameIonome_ehat_RC)
E_hat_RC_mat <- do.call(cbind, Ionome_ehat_RC)

E_hat_mat_RC_scaled <- scale(E_hat_RC_mat)

Ki_RC <- E_hat_mat_RC_scaled %*% t(E_hat_mat_RC_scaled)
Ki_RC <- Ki_RC / ncol(E_hat_mat_RC_scaled)

fileNameIonome_ehat_RD <- paste0(saveDir, "9.54.model1.RD.ehat.RDS")
Ionome_ehat_RD <- readRDS(fileNameIonome_ehat_RD)
E_hat_RD_mat <- do.call(cbind, Ionome_ehat_RD)

E_hat_mat_RD_scaled <- scale(E_hat_RD_mat)

Ki_RD <- E_hat_mat_RD_scaled  %*% t(E_hat_mat_RD_scaled)
Ki_RD <- Ki_RD / ncol(E_hat_mat_RD_scaled)

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

phenoRIL <- rbind(pheno19RC, pheno19RD, pheno20RC, pheno20RD)

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

phenoAll <- rbind(phenoRIL, phenoF2)

dim(Km)
dim(phenoAll)




###### 2.  ######
##### 2.1. Model 1 y = Y + u + e #####
treatments <- c("RD", "RC")
traits <- colnames(phenoAll)

results1 <- list()

for(i in seq_along(treatments)) {
  
  treatment <- treatments[i]
  
  pheno_RIL <- phenoAll[substr(rownames(phenoAll), 11, 12) == treatment, ]
  pheno_F2 <- phenoAll[substr(rownames(phenoAll), 17, 18) == treatment, ]
  
  common_RIL <- intersect(rownames(pheno_RIL), rownames(Km_new))
  Km_RIL <- Km_new[common_RIL, common_RIL]
  pheno_RIL <- pheno_RIL[common_RIL, , drop=FALSE] 

  common_F2 <- intersect(rownames(pheno_F2), rownames(Km))

  col_idx <- match(substr(common_RIL, 1, 7), substr(colnames(Km), 1, 7))

  Km_F2 <- Km[common_F2, col_idx]

  colnames(Km_F2) <- common_RIL

  year <- substr(rownames(pheno_RIL), 9, 10)
  year <- factor(year) 

  eta <- list(
    list(X = model.matrix(~ year-1), model = "FIXED"),
    list(K = Km_RIL, model = "RKHS")
  )
  
  set.seed(seedSet[i])
  
  results1[[treatment]] <- list()

  u_fitResults <- list()
  hg2Results <- list()
  residualResults <- list()
  
  for(trait in traits) {
    
    modelDir <- paste0(BGLRDir, "model1_", treatment, "_", trait, "/")
    if(!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)

    fit <- BGLR(
      y = pheno_RIL[[trait]],
      ETA = eta,
      nIter = nIter,
      burnIn = burnIn,
      verbose = FALSE,
      saveAt = modelDir
    )

    u_hat_Km <- Km_F2 %*% fit$ETA[[2]]$u 

    common_names <- intersect(rownames(pheno_F2), rownames(u_hat_Km))
    
    pheno_vec <- pheno_F2[common_names, trait]
    u_hat_Km <- u_hat_Km[common_names, ]

    correlation <- cor(pheno_vec, u_hat_Km, use="complete.obs")
    print(correlation)
    
    results1[[treatment]][[trait]] <- list(
      pheno = pheno_vec,
      pred = u_hat_Km,
      cor = correlation
    )

    sigma_g2 <- fit$ETA[[2]]$varU
    sigma_e2 <- fit$varE
    
    hg2Results[[trait]] <- sigma_g2 / (sigma_g2 + sigma_e2)
    
    u_fitResults[[trait]] <- u_hat_Km

    residuals <- pheno_RIL[[trait]] - fit$yHat
    names(residuals) <- rownames(pheno_RIL)
    
    residualResults[[trait]] <- residuals
  }

  fileNameu_fit <- paste0(saveDir, scriptID, "model1.", treatment, ".u_fit.RDS")
  saveRDS(u_fitResults, file = fileNameu_fit)
  
  fileNamehg2 <- paste0(saveDir, scriptID, "model1.", treatment, ".hg2.RDS")
  saveRDS(hg2Results, file = fileNamehg2)
  
  fileNameResidual <- paste0(saveDir, scriptID,
                             "model1.", treatment, ".residual.RDS")
  saveRDS(residualResults, file = fileNameResidual)
  
}


##### 2.2. Model 2 y = Y + f + u + e #####
treatments <- c("RD", "RC")
traits <- colnames(phenoAll)

results2 <- list()

for(i in seq_along(treatments)) {
  
  treatment <- treatments[i]
  
  if(treatment == "RC") {
    Ki_current <- Ki_RC
  } else {
    Ki_current <- Ki_RD
  }
  
  pheno_RIL <- phenoAll[substr(rownames(phenoAll), 11, 12) == treatment, ]
  pheno_F2 <- phenoAll[substr(rownames(phenoAll), 17, 18) == treatment, ]
  
  common_RIL <- Reduce(intersect, list(rownames(pheno_RIL), rownames(Km_new), rownames(Ki_current)))
  Km_RIL <- Km_new[common_RIL, common_RIL]
  Ki_RIL <- Ki_current[common_RIL, common_RIL]
  pheno_RIL <- pheno_RIL[common_RIL, , drop=FALSE] 

  common_F2 <- Reduce(intersect, list(rownames(pheno_F2), rownames(Km)))

  col_idx <- match(substr(common_RIL, 1, 7), substr(colnames(Km), 1, 7))

  Km_F2 <- Km[common_F2, col_idx]

  year <- substr(rownames(pheno_RIL), 9, 10)
  year <- factor(year) 

  eta <- list(
    list(X = model.matrix(~ year-1), model = "FIXED"),
    list(K = Km_RIL, model = "RKHS"),
    list(K = Ki_RIL, model = "RKHS")
  )
  
  set.seed(seedSet[i])
  
  results2[[treatment]] <- list()

  u_fitResults <- list()
  hg2Results <- list()
  hf2Results <- list()
  f_hat_RILResults <- list()
  residualResults <- list()
  
  for(trait in traits) {
    
    modelDir <- paste0(BGLRDir, "model2_", treatment, "_", trait, "/")
    if(!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)

    fit <- BGLR(
      y = pheno_RIL[[trait]],
      ETA = eta,
      nIter = nIter,
      burnIn = burnIn,
      verbose = FALSE,
      saveAt = modelDir
    )

    u_hat_Km <- Km_F2 %*% fit$ETA[[2]]$u  
    
    pred <- u_hat_Km

    common_names <- intersect(rownames(pheno_F2), rownames(pred))
    
    u_hat_Km <- u_hat_Km[common_names, ]
    
    pheno_vec <- pheno_F2[common_names, trait]
    pred <- pred[common_names, ]

    correlation <- cor(pheno_vec, pred, use="complete.obs")
    print(correlation)
    
    results2[[treatment]][[trait]] <- list(
      pheno = pheno_vec,
      pred = pred,
      cor = correlation
    )
    
    u_fitResults[[trait]] <- u_hat_Km

    sigma_g2 <- fit$ETA[[2]]$varU
    sigma_f2 <- fit$ETA[[3]]$varU
    sigma_e2 <- fit$varE
    
    hg2Results[[trait]] <- sigma_g2 / (sigma_g2 + sigma_f2 + sigma_e2)
    hf2Results[[trait]] <- sigma_f2 / (sigma_g2 + sigma_f2 + sigma_e2)

    f_hat_RIL <- Ki_RIL %*% fit$ETA[[3]]$u
    f_hat_RIL <- f_hat_RIL - mean(f_hat_RIL)
    f_hat_RILResults[[trait]] <- f_hat_RIL

    residuals <- pheno_RIL[[trait]] - fit$yHat
    names(residuals) <- rownames(pheno_RIL)
    
    residualResults[[trait]] <- residuals
    
  }

  fileNameu_fit <- paste0(saveDir, scriptID, "model2.", treatment, ".u_fit.RDS")
  saveRDS(u_fitResults, file = fileNameu_fit)
  
  fileNamehg2 <- paste0(saveDir, scriptID, "model2.", treatment, ".hg2.RDS")
  saveRDS(hg2Results, file = fileNamehg2)
  
  fileNamehf2 <- paste0(saveDir, scriptID, "model2.", treatment, ".hf2.RDS")
  saveRDS(hf2Results, file = fileNamehf2)
  
  fileNamef_hat_RIL <- paste0(saveDir, scriptID, "model2.", treatment, ".f_hat_RIL.RDS")
  saveRDS(f_hat_RILResults, file = fileNamef_hat_RIL)
  
  fileNameResidual <- paste0(saveDir, scriptID,
                             "model2.", treatment, ".residual.RDS")
  saveRDS(residualResults, file = fileNameResidual)
  
}


##### 2.3.  #####
summary_list <- list()

for(treatment in treatments) {
  for(trait in traits) {
    cor1 <- results1[[treatment]][[trait]]$cor
    cor2 <- results2[[treatment]][[trait]]$cor
    
    summary_list <- append(summary_list, list(
      data.frame(
        Treatment = treatment,
        Trait = trait,  
        Model1 = cor1,
        Model2 = cor2
      )
    ))
  }
}

summary_df <- bind_rows(summary_list)

fileNamesummary_df <- paste0(saveDir, scriptID,  "prediction.1920RIL.21F2.csv")
write.csv(summary_df, file = fileNamesummary_df, row.names = FALSE)



