###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(tidyverse)


##### 0.2. Define Parameters #####
scriptID <- "9.69."


##### 0.3. Define Directories #####
dataDir <- "results/ionome/RIL_F2/"
saveDir <- paste0("results/ionome/RIL_F2/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
#9.55 prediction
#Model 1 y = Y + u + e
fileNameModel1_RC_pred <- paste0(dataDir, "9.55.model1.RC.u_fit.RDS")
Model1_RC_pred <- readRDS(fileNameModel1_RC_pred)

fileNameModel1_RD_pred <- paste0(dataDir, "9.55.model1.RD.u_fit.RDS")
Model1_RD_pred <- readRDS(fileNameModel1_RD_pred)

#Model 2 y = Y + u + ionome + e
fileNameModel2_RC_pred <- paste0(dataDir, "9.55.model2.RC.u_fit.RDS")
Model2_RC_pred <- readRDS(fileNameModel2_RC_pred)

fileNameModel2_RD_pred <- paste0(dataDir, "9.55.model2.RD.u_fit.RDS")
Model2_RD_pred <- readRDS(fileNameModel2_RD_pred)

#9.45　2021
#Model 1 y = u + e 
fileNameModel1_RC_fit <- paste0(dataDir, "9.45.model1.RC.u_fit.RDS")
Model1_RC_fit <- readRDS(fileNameModel1_RC_fit)

fileNameModel1_RD_fit <- paste0(dataDir, "9.45.model1.RD.u_fit.RDS")
Model1_RD_fit <- readRDS(fileNameModel1_RD_fit)

# 9.67 prediction
#Model 3 y = u + distance + e
fileNameModel3_RC_pred <- paste0(dataDir, "9.67.model1.RC.u_fit.RDS")
Model3_RC_pred <- readRDS(fileNameModel3_RC_pred)

fileNameModel3_RD_pred <- paste0(dataDir, "9.67.model1.RD.u_fit.RDS")
Model3_RD_pred <- readRDS(fileNameModel3_RD_pred)

#Model 4 y = u + distance_2 + e
fileNameModel4_RC_pred <- paste0(dataDir, "9.67.model2.RC.u_fit.RDS")
Model4_RC_pred <- readRDS(fileNameModel4_RC_pred)

fileNameModel4_RD_pred <- paste0(dataDir, "9.67.model2.RD.u_fit.RDS")
Model4_RD_pred <- readRDS(fileNameModel4_RD_pred)

#Model 5 y = u + distance_0.5 + e
fileNameModel5_RC_pred <- paste0(dataDir, "9.67.model3.RC.u_fit.RDS")
Model5_RC_pred <- readRDS(fileNameModel5_RC_pred)

fileNameModel5_RD_pred <- paste0(dataDir, "9.67.model3.RD.u_fit.RDS")
Model5_RD_pred <- readRDS(fileNameModel5_RD_pred)

#Model 6 y = u + distance_exp + e
fileNameModel6_RC_pred <- paste0(dataDir, "9.67.model4.RC.u_fit.RDS")
Model6_RC_pred <- readRDS(fileNameModel6_RC_pred)

fileNameModel6_RD_pred <- paste0(dataDir, "9.67.model4.RD.u_fit.RDS")
Model6_RD_pred <- readRDS(fileNameModel6_RD_pred)


###### 2. Correlation ######
phenoCol <- names(Model1_RC_pred)

Blocks <- c("RC", "RD")

###### Model 1 ###### 
Model1_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model1_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model1_RD_pred
    Model_fit  <- Model1_RD_fit
  }

  cor_results <- list()
  
  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]
    
    common_names <- intersect(names(pred), names(fit))
    
    print(length(common_names))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]

    correlation <- cor(pred_common, fit_common, use = "complete.obs")

    cor_results[[trait]] <- correlation
  }

  Model1_results[[Block]] <- cor_results
  
}

###### Model2 ######
Model2_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model2_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model2_RD_pred
    Model_fit  <- Model1_RD_fit
  }

  cor_results <- list()

  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]

    common_names <- intersect(names(pred), names(fit))
    
    print(length(common_names))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]

    correlation <- cor(pred_common, fit_common, use = "complete.obs")

    cor_results[[trait]] <- correlation
  }
  
  Model2_results[[Block]] <- cor_results
  
}

###### Model3 ######
Model3_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model3_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model3_RD_pred
    Model_fit  <- Model1_RD_fit
  }
  
  cor_results <- list()
  
  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]
    
    common_names <- intersect(names(pred), names(fit))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]
    
    correlation <- cor(pred_common, fit_common,
                       use = "complete.obs")
    
    cor_results[[trait]] <- correlation
  }

  Model3_results[[Block]] <- cor_results
}

###### Model4 ######
Model4_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model4_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model4_RD_pred
    Model_fit  <- Model1_RD_fit
  }
  
  cor_results <- list()
  
  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]
    
    common_names <- intersect(names(pred), names(fit))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]
    
    correlation <- cor(pred_common, fit_common,
                       use = "complete.obs")
    
    cor_results[[trait]] <- correlation
  }

  Model4_results[[Block]] <- cor_results
}

###### Model5 ######
Model5_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model5_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model5_RD_pred
    Model_fit  <- Model1_RD_fit
  }
  
  cor_results <- list() 
  
  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]
    
    common_names <- intersect(names(pred), names(fit))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]
    
    correlation <- cor(pred_common, fit_common,
                       use = "complete.obs")
    
    cor_results[[trait]] <- correlation
  }

  Model5_results[[Block]] <- cor_results
}

###### Model 6 ######
Model6_results <- list()

for (Block in Blocks) {
  
  if (Block == "RC") {
    Model_pred <- Model6_RC_pred
    Model_fit  <- Model1_RC_fit
  } else if (Block == "RD") {
    Model_pred <- Model6_RD_pred
    Model_fit  <- Model1_RD_fit
  }
  
  cor_results <- list()
  
  for (trait in phenoCol) {
    
    pred <- Model_pred[[trait]]
    fit <- Model_fit[[trait]]
    
    common_names <- intersect(names(pred), names(fit))
    
    pred_common <- pred[common_names]
    fit_common <- fit[common_names]
    
    correlation <- cor(pred_common, fit_common,
                       use = "complete.obs")
    
    cor_results[[trait]] <- correlation
  }
  
  Model6_results[[Block]] <- cor_results
}


#all
all_results <- list(
  "Model1" = Model1_results,
  "Model2" = Model2_results,
  "Model3" = Model3_results,
  "Model4" = Model4_results,
  "Model5" = Model5_results,
  "Model6" = Model6_results
)




###### 3. plot ######
#accuracy
df_all <- lapply(names(all_results), function(model_name) {
  model_list <- all_results[[model_name]]
  
  lapply(names(model_list), function(block_name) {
    block_list <- model_list[[block_name]]
    data.frame(
      Model = model_name,
      Environment = block_name,
      Pheno = names(block_list),
      Cor = unlist(block_list)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

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

file_labels <- c(
  "RC" = "Control",
  "RD" = "Drought"
)

label_map <- c(
  "Model1" = "Model 3-1",
  "Model2" = "Model 3-2",
  "Model3" = "Model 3-3",
  "Model4" = "Model 3-4",
  "Model5" = "Model 3-5",
  "Model6" = "Model 3-6"
)

plotCor <- ggplot(df_all, aes(x = Pheno, y = Cor, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ Environment,  labeller = labeller(Environment = file_labels)) +
  theme_bw() +
  labs(y = "Accuracy", x = "Traits") +
  scale_x_discrete(labels = trait_labels) +
  scale_fill_manual(values = c(
    "Model1" = "#E69F00",
    "Model2" = "#56B4E9",
    "Model3" = "#C7E9C0",
    "Model4" = "#A1D99B",
    "Model5" = "#74C476",
    "Model6" = "#238B45"
  ),
  labels = label_map) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold")
  )

print(plotCor)

ggsave(
  filename = paste0(saveDir, scriptID, "correlation.1920RIL.21F2.uhat.png"),
  plot = plotCor,
  width = 15,
  height = 10
)



