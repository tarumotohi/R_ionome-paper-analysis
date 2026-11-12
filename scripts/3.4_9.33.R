###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(tidyverse)
library(ggplot2)


##### 0.2. Define Parameters #####
scriptID <- "9.33."


##### 0.3. Define Directories #####
dataDir <- "results/ionome/environment/"
saveDir <- paste0("results/ionome/environment/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNameModel2 <- paste0(dataDir, "9.32.model2.all.corAve.RDS")
model2 <- readRDS(fileNameModel2)

fileNameModel4 <- paste0(dataDir, "9.32.model4.all.corAve.RDS")
model4 <- readRDS(fileNameModel4)

fileNameModel5 <- paste0(dataDir, "9.32.model5.all.corAve.RDS")
model5 <- readRDS(fileNameModel5)


##### 1.2. Modify Data #####
model2_df <- bind_rows(lapply(names(model2), function(env) {
  data.frame(
    Environment = env,
    Trait = names(model2[[env]]),
    Value = as.numeric(model2[[env]]),
    Model = "Model2"
  )
}))

model4_df <- bind_rows(lapply(names(model4), function(env) {
  data.frame(
    Environment = env,
    Trait = names(model4[[env]]),
    Value = as.numeric(model4[[env]]),
    Model = "Model4"
  )
}))

model5_df <- bind_rows(lapply(names(model5), function(env) {
  data.frame(
    Environment = env,
    Trait = names(model5[[env]]),
    Value = as.numeric(model5[[env]]),
    Model = "Model5"
  )
}))

df_245 <- rbind(model2_df, model4_df, model5_df)




###### 2. plot ######
##Model 2,4,5
df_245_modi <- df_245 %>%
  mutate(Model = factor(Model, levels = c("Model4", "Model2", "Model5")))

df_245_modi <- df_245_modi %>%
  filter(Trait != "All")

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
  "19RC" = "2019-Control",
  "19RD" = "2019-Drought",
  "20RC" = "2020-Control",
  "20RD" = "2020-Drought"
)

label_map <- c(
  "Model4" = "Model 2-1",
  "Model2" = "Model 2-2",
  "Model5" = "Model 2-3"
)

df_245_plot_modi <- ggplot(df_245_modi, aes(x = Trait, y = Value, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ Environment,  
             labeller = labeller(Environment = file_labels)) +
  theme_bw() +
  labs(title = "",
       y = "Accuracy", x = "Traits") +
  scale_x_discrete(labels = trait_labels) +
  scale_fill_manual(
    values = c(
      "Model2" = "#56B4E9", 
      "Model4" = "#E69F00",
      "Model5" = "#0072B2" 
    ),
    labels = label_map)+
  guides(
    fill = guide_legend(
      title = "Model"
    )
  ) +
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

ggsave(
  filename = paste0(saveDir, scriptID, "GP.environment.allTraits_245_modi.png"),
  plot = df_245_plot_modi,
  width = 15,
  height = 10
)



