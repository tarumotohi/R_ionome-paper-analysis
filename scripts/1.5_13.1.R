###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(Hmisc)


##### 0.2. Define Parameters #####
scriptID <- "13.1."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "results/ionome/correlation/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}

yearBlock <- "19RD"




###### 1. Data Preparation ######
##### 1.1. Load Data #####
#RIL 19RD ionome
fileNameIonome <- paste0(dataDir, "9.0.ionRILMean.19RD.csv")
ionome <- read.csv(file = fileNameIonome, row.names = 1)

ionome[, 2] <- gsub("-", "N", ionome[, 2])
rownames(ionome) <- ionome[, 2]
ionome <- ionome[, -c(1,2)]

#RIL 19RD phenotype
fileNamePheno <- paste0("data/", "Tottori_RIL.csv")
pheno <- read.csv(file = fileNamePheno, row.names = 1)
phenoSel <- pheno[,18:32]
colnames(phenoSel) <- gsub("\\.1$", "", colnames(phenoSel))


##### 1.2. Modify Data #####
common <- intersect(rownames(ionome), rownames(phenoSel))
ionome_sub <- ionome[common, ]
pheno_sub <- phenoSel[common, ]




###### 2. Correlation ######
res <- rcorr(as.matrix(ionome_sub), as.matrix(pheno_sub))

cor_matrix <- res$r

fileNameCor_matrix <- paste0(dataDir, scriptID, "Tottori.RIL.", yearBlock, ".ion.pheno.cor.csv")
write.csv(x = cor_matrix, file = fileNameCor_matrix)

p_matrix <- res$P

fileNameP_matrix <- paste0(dataDir, scriptID, "Tottori.RIL.", yearBlock, ".ion.pheno.p.csv")
write.csv(x = p_matrix, file = fileNameP_matrix)



