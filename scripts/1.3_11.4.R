###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(tidyverse)
library(emmeans)
library(patchwork)
library(Hmisc)

get_stars <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}


##### 0.2. Define Parameters #####
scriptID <- "11.4."


##### 0.3. Define Directories #####
dataDir <- "data/"
saveDir <- "results/ionome/PCA/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameScoresIonpcaIon19.20 <- paste0(saveDir, "11.3.pcaIon19.20.scores.csv")
ionBLUE.PCA.scores <- read.csv(file = fileNameScoresIonpcaIon19.20, header = TRUE)

#RIL 2019RC phenotype
fileNamePheno19RC <- paste0(dataDir, "Tottori_RIL.csv")
pheno19RC <- read.csv(file = fileNamePheno19RC, row.names = 1)
phenoSel19RC <- pheno19RC[, 2:16]

#RIL 19RD phenotype
fileNamePheno19RD <- paste0(dataDir, "Tottori_RIL.csv")
pheno19RD <- read.csv(file = fileNamePheno19RD, row.names = 1)
phenoSel19RD <- pheno19RD[,18:32]
colnames(phenoSel19RD) <- gsub("\\.1$", "", colnames(phenoSel19RD))

#RIL 20RC phenotype
fileNamePheno20RC <- paste0(dataDir, "Tottori.2020.RIL.csv")
pheno20RC <- read.csv(file = fileNamePheno20RC, row.names = 1)

phenoSel20RC <- pheno20RC[pheno20RC$block == "RC", ]
phenoSel20RC$block <- NULL   

#RIL 20RD phenotype
fileNamePheno20RD <- paste0(dataDir, "Tottori.2020.RIL.csv")
pheno20RD <- read.csv(file = fileNamePheno20RD, row.names = 1)

phenoSel20RD <- pheno20RD[pheno20RD$block == "RD", ]
phenoSel20RD$block <- NULL   


##### 1.2. Modify Data #####
all_traits <- c("BL", "FSW", "GS", "NB", "NN", "NS", "TENSW", "TSW")

phenoSel19RC <- phenoSel19RC[ ,all_traits]
phenoSel19RD <- phenoSel19RD[ ,all_traits]

phenoSel19RC$variety <- rownames(phenoSel19RC)
phenoSel19RD$variety <- rownames(phenoSel19RD)




###### 2. BLUEs of the eight aboveground traits  ######
#### RC ####
phenoSel19RC$Year <- "2019"
phenoSel20RC$Year <- "2020"

ionRC <- rbind(phenoSel19RC, phenoSel20RC)

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

blue.data_RC_pheno <- Reduce(function(x, y) merge(x, y, by = "variety"), blue.list_RC)

#### RD ####
phenoSel19RD$Year <- "2019"
phenoSel20RD$Year <- "2020"

ionRD <- rbind(phenoSel19RD, phenoSel20RD)

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

blue.data_RD_pheno <- Reduce(function(x, y) merge(x, y, by = "variety"), blue.list_RD)




###### 3. Correlation with Principal Component Scores  ######
ionBLUE.PCA.scores
colnames(ionBLUE.PCA.scores)[1] <- "variety"

common_variety <- Reduce(
  intersect,
  list(
    ionBLUE.PCA.scores$variety,
    blue.data_RC_pheno$variety,
    blue.data_RD_pheno$variety
  )
)

length(common_variety)

ionBLUE.PCA.scores_common <- ionBLUE.PCA.scores[
  match(common_variety, ionBLUE.PCA.scores$variety),
]

blue.data_RC_pheno_common <- blue.data_RC_pheno[
  match(common_variety, blue.data_RC_pheno$variety),
]

blue.data_RD_pheno_common <- blue.data_RD_pheno[
  match(common_variety, blue.data_RD_pheno$variety),
]

ionBLUE.PCA.scores_common <- ionBLUE.PCA.scores_common[, setdiff(names(ionBLUE.PCA.scores_common), "variety")]
blue.data_RD_pheno_common <- blue.data_RD_pheno_common[, setdiff(names(blue.data_RD_pheno_common), "variety")]
blue.data_RC_pheno_common  <- blue.data_RC_pheno_common[, setdiff(names(blue.data_RC_pheno_common), "variety")]

res_RC <- rcorr(
  as.matrix(ionBLUE.PCA.scores_common),
  as.matrix(blue.data_RC_pheno_common)
)

res_RD <- rcorr(
  as.matrix(ionBLUE.PCA.scores_common),
  as.matrix(blue.data_RD_pheno_common)
)

cor_RC <- res_RC$r[
  colnames(ionBLUE.PCA.scores_common),
  colnames(blue.data_RC_pheno_common),
  drop = FALSE
]

cor_RD <- res_RD$r[
  colnames(ionBLUE.PCA.scores_common),
  colnames(blue.data_RD_pheno_common),
  drop = FALSE
]

p_RC <- res_RC$P[
  colnames(ionBLUE.PCA.scores_common),
  colnames(blue.data_RC_pheno_common),
  drop = FALSE
]

p_RD <- res_RD$P[
  colnames(ionBLUE.PCA.scores_common),
  colnames(blue.data_RD_pheno_common),
  drop = FALSE
]

cor_CD <- rbind(
  cor_RC,
  cor_RD
)

rownames(cor_CD) <- c(
  paste0(rownames(cor_RC), "-Control"),
  paste0(rownames(cor_RD), "-Drought")
)

p_CD <- rbind(
  p_RC,
  p_RD
)

rownames(p_CD) <- c(
  paste0(rownames(p_RC), "-Control"),
  paste0(rownames(p_RD), "-Drought")
)

trait_labels <- c(
  "FSW" = "Fresh Shoot Weight",
  "GS" = "Growth Stage",
  "NN" = "Number of Nodes",
  "NB" = "Number of Branches",
  "BL" = "Main Stem Length ",
  "TENSW" = "10-Seed Weight ",
  "TSW" = "Total Seed Weight ",
  "NS" = "Number of Seeds"
)

cor_CD <- cor_CD[, names(trait_labels)]
colnames(cor_CD) <- trait_labels

p_CD <- p_CD[, names(trait_labels), drop = FALSE]
colnames(p_CD) <- trait_labels

rows <- c(
  "PC1-Control",
  "PC1-Drought",
  "PC2-Control",
  "PC2-Drought"
)

cor_CD_PC1_2 <- cor_CD[rows, , drop = FALSE]
cor_CD_PC1_2 <- round(cor_CD_PC1_2, 2)

p_CD_PC1_2 <- p_CD[
  rows,
  ,
  drop = FALSE
]

p_vec <- as.vector(p_CD_PC1_2)

p_adj_vec <- p.adjust(
  p_vec,
  method = "holm"
)

p_CD_PC1_2_holm <- matrix(
  p_adj_vec,
  nrow = nrow(p_CD_PC1_2),
  ncol = ncol(p_CD_PC1_2),
  dimnames = dimnames(p_CD_PC1_2)
)

stars_CD_PC1_2 <- matrix(
  sapply(
    p_CD_PC1_2_holm,
    get_stars
  ),
  nrow = nrow(p_CD_PC1_2_holm),
  ncol = ncol(p_CD_PC1_2_holm),
  dimnames = dimnames(p_CD_PC1_2_holm)
)

traitPC1_2_df <- cor_CD_PC1_2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Set") %>%
  pivot_longer(
    cols = -Set,
    names_to = "Trait",
    values_to = "Correlation"
  )

stars_df <- stars_CD_PC1_2 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Set") %>%
  pivot_longer(
    cols = -Set,
    names_to = "Trait",
    values_to = "stars"
  )

trait_order <- c(
  "Main Stem Length ",
  "Fresh Shoot Weight",
  "Growth Stage",
  "Number of Branches",
  "Number of Nodes",
  "Number of Seeds",
  "10-Seed Weight ",
  "Total Seed Weight "
)

traitPC1_2_df$Trait <- factor(
  traitPC1_2_df$Trait,
  levels = trait_order
)

set_order <- c(
  "PC1-Control",
  "PC1-Drought",
  "PC2-Control",
  "PC2-Drought"
)

traitPC1_2_df <- left_join(
  traitPC1_2_df,
  stars_df,
  by = c("Set", "Trait")
)

p <- ggplot(traitPC1_2_df, aes(x = Trait, y = Set, fill = Correlation)) +
  geom_tile(color = NA) +
  scale_x_discrete(limits = trait_order) +
  scale_y_discrete(limits = rev(set_order)) +
  scale_fill_gradient2(
    low = "#1F77B4",
    mid = "white",
    high = "#FF7F0E",
    midpoint = 0,
    limits = c(-1, 1),
    space = "Lab",
    name = ""
  ) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", Correlation), stars)
    ), size = 3.5
  ) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  ) +
  coord_fixed()

print(p)

ggsave(filename = paste0(saveDir, scriptID, "PCA.BLUE.ionPheno.correlation.png"),
       plot = p, width = 8, height = 4)




###### 4. Correlation coefficients among the BLUEs of the eight aboveground traits  ######
common_variety_pheno <- intersect(
  blue.data_RC_pheno$variety,
  blue.data_RD_pheno$variety
)

blue_RC_common <- blue.data_RC_pheno[
  match(common_variety_pheno, blue.data_RC_pheno$variety),
]

blue_RD_common <- blue.data_RD_pheno[
  match(common_variety_pheno, blue.data_RD_pheno$variety),
]

blue_RC_common_mat <- blue_RC_common[
  , setdiff(names(blue_RC_common), "variety")
]

blue_RD_common_mat <- blue_RD_common[
  , setdiff(names(blue_RD_common), "variety")
]

res_pheno_RC <- rcorr(
  as.matrix(blue_RC_common_mat),
  type = "pearson"
)

cor_pheno_RC <- res_pheno_RC$r
p_pheno_RC <- res_pheno_RC$P

res_pheno_RD <- rcorr(
  as.matrix(blue_RD_common_mat),
  type = "pearson"
)

cor_pheno_RD <- res_pheno_RD$r
p_pheno_RD <- res_pheno_RD$P

cor_pheno_RC <- cor_pheno_RC[
  names(trait_labels),
  names(trait_labels),
  drop = FALSE
]

cor_pheno_RD <- cor_pheno_RD[
  names(trait_labels),
  names(trait_labels),
  drop = FALSE
]

p_pheno_RC <- p_pheno_RC[
  names(trait_labels),
  names(trait_labels),
  drop = FALSE
]

p_pheno_RD <- p_pheno_RD[
  names(trait_labels),
  names(trait_labels),
  drop = FALSE
]

colnames(cor_pheno_RC) <- trait_labels
rownames(cor_pheno_RC) <- trait_labels

colnames(cor_pheno_RD) <- trait_labels
rownames(cor_pheno_RD) <- trait_labels

colnames(p_pheno_RC) <- trait_labels
rownames(p_pheno_RC) <- trait_labels

colnames(p_pheno_RD) <- trait_labels
rownames(p_pheno_RD) <- trait_labels

p_lower_RC <- p_pheno_RC[lower.tri(p_pheno_RC)]

p_adj_lower_RC <- p.adjust(
  p_lower_RC,
  method = "holm"
)

p_lower_RD <- p_pheno_RD[lower.tri(p_pheno_RD)]

p_adj_lower_RD <- p.adjust(
  p_lower_RD,
  method = "holm"
)

p_adj_RC <- matrix(NA, nrow = nrow(p_pheno_RC), ncol = ncol(p_pheno_RC), dimnames = dimnames(p_pheno_RC))
p_adj_RC[lower.tri(p_adj_RC)] <- p_adj_lower_RC

p_adj_RD <- matrix(NA, nrow = nrow(p_pheno_RD), ncol = ncol(p_pheno_RD), dimnames = dimnames(p_pheno_RD))
p_adj_RD[lower.tri(p_adj_RD)] <- p_adj_lower_RD

p_adj_RC[upper.tri(p_adj_RC)] <- t(p_adj_RC)[upper.tri(p_adj_RC)]
p_adj_RD[upper.tri(p_adj_RD)] <- t(p_adj_RD)[upper.tri(p_adj_RD)]

stars_pheno_RC <- matrix(
  sapply(p_adj_RC, function(x) if(is.na(x)) "" else get_stars(x)),
  nrow = nrow(p_adj_RC), ncol = ncol(p_adj_RC), dimnames = dimnames(p_adj_RC)
)

stars_pheno_RD <- matrix(
  sapply(p_adj_RD, function(x) if(is.na(x)) "" else get_stars(x)),
  nrow = nrow(p_adj_RD), ncol = ncol(p_adj_RD), dimnames = dimnames(p_adj_RD)
)

write.csv(cor_pheno_RC, file = paste0(saveDir, scriptID, "pheno.correlation.Control.csv"))
write.csv(cor_pheno_RD, file = paste0(saveDir, scriptID, "pheno.correlation.Drought.csv"))


##### RC #####
pheno_RC_df <- cor_pheno_RC %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(
    cols = -Var1,
    names_to = "Var2",
    values_to = "Correlation"
  )

stars_RC_df <- stars_pheno_RC %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(
    cols = -Var1,
    names_to = "Var2",
    values_to = "stars"
  )

pheno_RC_df <- left_join(pheno_RC_df, stars_RC_df, by = c("Var1", "Var2"))

trait_order_pheno <- trait_order 

p_rc <- ggplot(pheno_RC_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_x_discrete(limits = trait_order_pheno) +
  scale_y_discrete(limits = rev(trait_order_pheno)) +
  scale_fill_gradient2(
    low = "#1F77B4",
    mid = "white",
    high = "#FF7F0E",
    midpoint = 0,
    limits = c(-1, 1),
    space = "Lab",
    name = "Correlation"
  ) +
  geom_text(
    aes(label = ifelse(Correlation == 1, "1.00", paste0(sprintf("%.2f", Correlation), stars))),
    size = 3
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Control", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed()

print(p_rc)

##### RD #####
pheno_RD_df <- cor_pheno_RD %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(
    cols = -Var1,
    names_to = "Var2",
    values_to = "Correlation"
  )

stars_RD_df <- stars_pheno_RD %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(
    cols = -Var1,
    names_to = "Var2",
    values_to = "stars"
  )

pheno_RD_df <- left_join(pheno_RD_df, stars_RD_df, by = c("Var1", "Var2"))

p_rd <- ggplot(pheno_RD_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_x_discrete(limits = trait_order_pheno) +
  scale_y_discrete(limits = rev(trait_order_pheno)) +
  scale_fill_gradient2(
    low = "#1F77B4",
    mid = "white",
    high = "#FF7F0E",
    midpoint = 0,
    limits = c(-1, 1),
    space = "Lab",
    name = "Correlation"
  ) +
  geom_text(
    aes(label = ifelse(Correlation == 1, "1.00", paste0(sprintf("%.2f", Correlation), stars))),
    size = 3
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "Drought", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed()

print(p_rd)

p_combined <- p_rc + p_rd + plot_layout(guides = "collect") & 
  theme(legend.position = "right")

ggsave(
  filename = paste0(saveDir, scriptID, "pheno.correlation.Control_Drought.combined.png"),
  plot = p_combined, width = 12, height = 6
)



