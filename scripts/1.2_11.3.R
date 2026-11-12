###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(pcaMethods)
library(ggplot2)
library(emmeans)


##### 0.2. Define Parameters #####
scriptID <- "11.3."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "results/ionome/PCA/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameIonome19RC <- paste0(dataDir, "9.0.ionRILMean.19RC.csv")
ion19RC <- read.csv(file = fileNameIonome19RC, header = TRUE)

fileNameIonome19RD <- paste0(dataDir, "9.0.ionRILMean.19RD.csv")
ion19RD <- read.csv(file = fileNameIonome19RD, header = TRUE)

fileNameIonome20RC <- paste0(dataDir, "9.0.ionRILMean.20RC.csv")
ion20RC <- read.csv(file = fileNameIonome20RC, header = TRUE)

fileNameIonome20RD <- paste0(dataDir, "9.0.ionRILMean.20RD.csv")
ion20RD <- read.csv(file = fileNameIonome20RD, header = TRUE)


##### 1.2. Modify Data #####
substr(ion19RC$variety, 4, 4) <- "N"
substr(ion19RD$variety, 4, 4) <- "N"
substr(ion20RC$variety, 4, 4) <- "N"
substr(ion20RD$variety, 4, 4) <- "N"

ion19RC <- ion19RC[, -c(1, 2)]
ion19RD <- ion19RD[, -c(1, 2)]
ion20RC <- ion20RC[, -c(1, 2)]
ion20RD <- ion20RD[, -c(1, 2)]




###### 2. BLUE  ######
##### 2.1. BLUE Calculation #####
#### RC (Control) ####
ion19RC$Year <- "2019"
ion20RC$Year <- "2020"

ionRC <- rbind(ion19RC, ion20RC)

ionRC$Year <- factor(ionRC$Year)
ionRC$variety <- factor(ionRC$variety)

traits <- setdiff(names(ionRC), c("variety", "Year"))

blue.list_RC <- list()

for(trait in traits){
  
  form <- as.formula(paste(trait, "~ Year + variety"))
  fit <- lm(form, data = ionRC)
  
  blue_RC <- as.data.frame(emmeans(fit, ~ variety))
  
  blue.list_RC[[trait]] <- blue_RC[, c("variety", "emmean")]
  names(blue.list_RC[[trait]])[2] <- trait
}

blue.data_RC <- Reduce(function(x, y) merge(x, y, by = "variety"), blue.list_RC)

#### RD (Drought) ####
ion19RD$Year <- "2019"
ion20RD$Year <- "2020"

ionRD <- rbind(ion19RD, ion20RD)

ionRD$Year <- factor(ionRD$Year)
ionRD$variety <- factor(ionRD$variety)

traits <- setdiff(names(ionRD), c("variety", "Year"))

blue.list_RD <- list()

for(trait in traits){
  
  form <- as.formula(paste(trait, "~ Year + variety"))
  fit <- lm(form, data = ionRD)
  
  blue_RD <- as.data.frame(emmeans(fit, ~ variety))
  
  blue.list_RD[[trait]] <- blue_RD[, c("variety", "emmean")]
  names(blue.list_RD[[trait]])[2] <- trait
}

blue.data_RD <- Reduce(function(x, y) merge(x, y, by = "variety"), blue.list_RD)


##### 2.2. PCA  #####
names(blue.data_RC)[-1] <- paste0(names(blue.data_RC)[-1], "_C")
names(blue.data_RD)[-1] <- paste0(names(blue.data_RD)[-1], "_D")

blue.data <- merge(blue.data_RC, blue.data_RD, by = "variety")
rownames(blue.data) <- blue.data[, 1]
blue.data <- blue.data[, -1]

pcaIon19.20 <- pcaMethods::pca(blue.data, method = "ppca", nPcs = 12, scale = "uv", completeObs = TRUE)

reference_var <- "As_C"
for (pc in 1:2) {
  if (pcaIon19.20@loadings[reference_var, pc] < 0) {
    pcaIon19.20@scores[, pc] <-
      -pcaIon19.20@scores[, pc]
    
    pcaIon19.20@loadings[, pc] <-
      -pcaIon19.20@loadings[, pc]
  }
}

slplot(pcaIon19.20, pcs = 1:2)

##### Scores ####
fileNameScoresIonpcaIon19.20 <- paste0(saveDir, scriptID, "pcaIon19.20.scores.csv")
write.csv(x = pcaIon19.20@scores, file = fileNameScoresIonpcaIon19.20)

#### Loadings ####
fileNameloadingsIonpcaIon19.20 <- paste0(saveDir, scriptID, "pcaIon19.20.loadings.ion_2.pdf")
pdf(file = fileNameloadingsIonpcaIon19.20, width = 15, height = 8)

load <- as.data.frame(loadings(pcaIon19.20)[, 1:2])

load$trait <- rownames(load)
load$group <- ifelse(grepl("_C$", load$trait), "C", "D")
load$label <- gsub("_C$|_D$", "", load$trait)

load_C <- load[load$group == "C", ]
load_D <- load[load$group == "D", ]

link <- merge(
  load_C,
  load_D,
  by = "label",
  suffixes = c("_C", "_D")
)

ggplot(load, aes(x = PC1, y = PC2, label = label, shape = group)) +
  geom_point(size = 6, color = "gray30") +
  geom_text(vjust = -1, size = 5, color = "black") +
  
  geom_segment(
    data = link,
    aes(
      x = PC1_C,
      y = PC2_C,
      xend = PC1_D,
      yend = PC2_D
    ),
    color = "gray60",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  scale_shape_manual(
    name = "Treatment",
    values = c("C" = 16,
               "D" = 1),
    labels = c(
      "C" = "Control",
      "D" = "Drought"
    ) 
  ) +
  theme_bw() +
  labs(
    title = "",
    x = "PC1 (17%)",
    y = "PC2 (10%)"
  ) +
  theme_minimal()+
  theme(
    text = element_text(size = 24),
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

dev.off()

#### Explained Variance Ratio ####
importance <- pcaIon19.20@R2

var_exp <- data.frame(
  PC = paste0("PC", 1:9),
  Variance = as.numeric(importance)[1:9] * 100
)

p_exp <-ggplot(var_exp, aes(x = PC, y = Variance)) +
  geom_col(fill = "#BBBBBB") +
  # geom_text(
  #   aes(label = sprintf("%.1f%%", Variance)),
  #   vjust = -0.3,
  #   size = 4
  # ) +
  labs(
    x = "Principal Component",
    y = "Variance explained (%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22)
  )

ggsave(
  filename = paste0(
    saveDir,
    scriptID,
    "pca.exp.pdf"
  ),
  plot = p_exp,
  width = 15,
  height = 8
)



