###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(BGLR)
library(tidyverse)
library(dplyr)


##### 0.2. Define Parameters #####
scriptID <- "9.26_3."


##### 0.3. Define Directories #####
dataDir <- "results/ionome/5Models/"
saveDir <- "results/ionome/plot_5Models/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileName19RC <- paste0(dataDir, "19RC/", "9.3_2.19RC.corAve.allModels.csv")
CV.19RC <- read.csv(file = fileName19RC)

fileName19RD <- paste0(dataDir, "19RD/", "9.4_2.19RD.corAve.allModels.csv")
CV.19RD <- read.csv(file = fileName19RD)

fileName20RC <- paste0(dataDir, "20RC/", "9.6_2.20RC.corAve.allModels.csv")
CV.20RC <- read.csv(file = fileName20RC)

fileName20RD <- paste0(dataDir, "20RD/", "9.7_2.20RD.corAve.allModels.csv")
CV.20RD <- read.csv(file = fileName20RD)


##### 1.2. Modify Data #####
CV.19RC <- CV.19RC[-c(16,17), ]
CV.19RD <- CV.19RD[-c(16,17), ]
CV.20RC <- CV.20RC[-c(9,10), ]
CV.20RD <- CV.20RD[-c(9,10), ]

CV.list <- list(
  "2019RC"       = CV.19RC,
  "2019RD"       = CV.19RD,
  "2020RC"       = CV.20RC,
  "2020RD"       = CV.20RD
)

CV.df <- bind_rows(CV.list, .id = "source")

combined_df <- bind_rows(
  lapply(names(CV.list), function(name) {
    df <- CV.list[[name]]
    df_long <- pivot_longer(df, cols = starts_with("Model"),
                            names_to = "Model", values_to = "Value")
    df_long$File <- name
    df_long
  })
)

combined_df$File <- factor(combined_df$File, levels = names(CV.list))




###### 2. Plot ######
traits <- c("BL", "FSW", "GS", "NB", "NN", "NS", "TENSW", "TSW")

combined_df <- combined_df %>%
  dplyr::filter(Trait %in% traits)

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
  "2019RC" = "2019-\nControl",
  "2019RD" = "2019-\nDrought",
  "2020RC" = "2020-\nControl",
  "2020RD" = "2020-\nDrought"
)

label_map <- c(
  "Model1"   = "Model 1-1",
  "Model1.5" = "Model 1-2",
  "Model2"   = "Model 1-3",
  "Model2.5" = "Model 1-4",
  "Model3"   = "Model 1-5"
)

p_all <- ggplot(combined_df, aes(x = File, y = Value, fill = Model, alpha = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ Trait, scales = "fixed", ncol = 3,
             labeller = labeller(Trait = trait_labels)) +
  scale_x_discrete(labels = file_labels) +
  labs(title = "",
       x = "Year and treatment",
       y = "Accuracy",
       fill = "") +
  scale_fill_manual(values = c(
    "Model1" = "#E69F00", 
    "Model1.5" = "#56B4E9", 
    "Model2" = "#56B4E9",
    "Model2.5" = "#009E73",
    "Model3" = "#009E73"   
  ),
  labels = label_map) +
  scale_alpha_manual(values = c(
    "Model1" = 1,
    "Model1.5" = 1,
    "Model2" = 0.55,
    "Model2.5" = 1,
    "Model3" = 0.6
  ),
  labels = label_map) +
  guides(
    fill = guide_legend(
      title = "Model",
      override.aes = list(
        alpha = c(1, 1, 0.55, 1, 0.6)
      )
    )
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),   
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(filename = paste0(saveDir, scriptID, "CV.sameTreat.8trait.all.png"),
       plot = p_all, width = 15, height = 10)



