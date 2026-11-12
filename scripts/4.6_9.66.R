###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(readr)


##### 0.2. Define Parameters #####
scriptID <- "9.66."


##### 0.3. Define Directories #####
dataDir <- "data/"
saveDir <- paste0("data/ionome/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}
BGLRDir <- paste0("results/BGLR/", scriptID, "/")
if(!dir.exists(BGLRDir)) dir.create(BGLRDir, recursive = TRUE)




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fieldMap <- data.frame(read_csv(paste0(dataDir, "Tottori.RIL.PlotPosition.csv")))

rownames(fieldMap) <- paste0(
  "E71N",
  sprintf("%03d", fieldMap$plot)
)


###### 2.  ######
##### 2.1.  #####
D <- as.matrix(dist(fieldMap[,c("x.real", "y.real")]))
D2 <- D^2

h <- 1 / median(D2[upper.tri(D2)])
kPos <- exp(-h * D2)

fileNamekPos <- paste0(saveDir, scriptID, "kPos.RDS")
saveRDS(kPos, file = fileNamekPos)

kPos_2 <- exp(-2 * h * D2)

fileNamekPos_2 <- paste0(saveDir, scriptID, "kPos_2.RDS")
saveRDS(kPos_2, file = fileNamekPos_2)

kPos_0.5 <- exp(-0.5 * h * D2)

fileNamekPos_0.5 <- paste0(saveDir, scriptID, "kPos_0.5.RDS")
saveRDS(kPos_0.5, file = fileNamekPos_0.5)

h2 <- 1 / median(D[upper.tri(D)])
kPos_exp <- exp(-h2 * D)

fileNamekPos_exp <- paste0(saveDir, scriptID, "kPos_exp.RDS")
saveRDS(kPos_exp, file = fileNamekPos_exp)



