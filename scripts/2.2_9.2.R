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
##### 1.1. Load Data #####
fileNameCross.0 <- paste0(dataDir, "geno_ril_bc.csv")
cross.0 <- read.cross(format = "csv", file = fileNameCross.0)
cross.0 <- jittermap(cross.0)

cross.0 <- convert2riself(cross.0)
cross.0

cross <- est.rf(cross.0)
cross <- replace.map(cross, est.map(cross, error.prob = 1e-04))

cross <- calc.genoprob(cross, step = 1, error.prob = 1e-04)

geno.progeny <- NULL
for(i in 1:length(cross$geno)) {
  print(i)
  link.pos <- names(cross$geno[[i]]$prob[1, ,1])
  is.pseudo <- substr(link.pos, 1, 3) == "loc"
  tmp <- cross$geno[[i]]$prob[ ,is.pseudo, 1]
  pos <- as.numeric(substr(colnames(tmp), 4, nchar(colnames(tmp))))
  colnames(tmp) <- paste0("C", formatC(i, width=2, flag="0"), 
                          "P", formatC(pos, width=3, flag="0"))
  geno.progeny <- cbind(geno.progeny, tmp)
}

genotype <- round(geno.progeny * 2) - 1

rownames(genotype) <- cross$pheno$id


##### 1.2. Kernel #####
Xm <- genotype
Km <- tcrossprod(Xm) / ncol(Xm)

fileNameKm <- paste0(saveDir, scriptID, "Km.RDS")
saveRDS(object = Km, file = fileNameKm)



