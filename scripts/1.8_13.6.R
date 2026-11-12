###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(ggplot2)
library(patchwork)

get_stars <- function(p) {
  if (is.na(p)) return("")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}


##### 0.2. Define Parameters #####
scriptID <- "13.6."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- "results/ionome/correlation/"
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
fileNamecor19RC <- paste0(dataDir, "13.0.Tottori.RIL.19RC.ion.pheno.cor.csv")
cor19RC <- read.csv(file = fileNamecor19RC, row.names = 1)

fileNamecor19RD <- paste0(dataDir, "13.1.Tottori.RIL.19RD.ion.pheno.cor.csv")
cor19RD <- read.csv(file = fileNamecor19RD, row.names = 1)

fileNamecor20RC <- paste0(dataDir, "13.2.Tottori.RIL.20RC.ion.pheno.cor.csv")
cor20RC <- read.csv(file = fileNamecor20RC, row.names = 1)

fileNamecor20RD <- paste0(dataDir, "13.3.Tottori.RIL.20RD.ion.pheno.cor.csv")
cor20RD <- read.csv(file = fileNamecor20RD, row.names = 1)

fileNamep19RC <- paste0(dataDir, "13.0.Tottori.RIL.19RC.ion.pheno.p.csv")
p19RC <- read.csv(file = fileNamep19RC, row.names = 1)

fileNamep19RD <- paste0(dataDir, "13.1.Tottori.RIL.19RD.ion.pheno.p.csv")
p19RD <- read.csv(file = fileNamep19RD, row.names = 1)

fileNamep20RC <- paste0(dataDir, "13.2.Tottori.RIL.20RC.ion.pheno.p.csv")
p20RC <- read.csv(file = fileNamep20RC, row.names = 1)

fileNamep20RD <- paste0(dataDir, "13.3.Tottori.RIL.20RD.ion.pheno.p.csv")
p20RD <- read.csv(file = fileNamep20RD, row.names = 1)




###### 2. plot ######
##### 2.1. 3 traits #####
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

target_traits <- c("FSW", "GS", "TSW")

cor19RC_sel <- cor19RC[1:19, target_traits]
p19RC_sel <- p19RC[1:19, target_traits]

cor19RD_sel <- cor19RD[1:19, target_traits]
p19RD_sel <- p19RD[1:19, target_traits]

cor20RC_sel <- cor20RC[1:19, target_traits]
p20RC_sel <- p20RC[1:19, target_traits]

cor20RD_sel <- cor20RD[1:19, target_traits]
p20RD_sel <- p20RD[1:19, target_traits]

cor_list_sel <- list(cor19RC_sel, cor19RD_sel, cor20RC_sel, cor20RD_sel)
p_list_sel   <- list(p19RC_sel, p19RD_sel, p20RC_sel, p20RD_sel)
set_names <- c("2019-Control","2019-Drought","2020-Control","2020-Drought")

plot_list_3 <- list()

ionome_order <- c(
  "As", "Br", "Ca", "Cd", "Cl", "Co", "Cs", "Cu",
  "Fe", "K", "Mn", "Mo", "Ni", "P", "Rb", "S",
  "Si", "Sr", "Zn"
)

for(tr in target_traits){
  
  trait_df <- data.frame()
  
  p_all <- unlist(
    lapply(p_list_sel, function(x) x[, tr])
  )
  
  p_adj_all <- p.adjust(
    p_all,
    method = "holm"
  )

  start <- 1
  
  for(i in 1:4){
    
    cor_vec <- cor_list_sel[[i]][, tr]
    
    p_adj <- p_adj_all[start:(start + 18)]
    
    star <- sapply(p_adj, get_stars)
    
    df <- data.frame(
      Ionome = factor(
        rownames(cor_list_sel[[i]]),
        levels = ionome_order
      ),
      Correlation = cor_vec,
      stars = star,
      Set = set_names[i],
      stringsAsFactors = FALSE
    )
    
    trait_df <- rbind(trait_df, df)
    
    start <- start + 19
  }
  
  p <- ggplot(trait_df, aes(x = Ionome, y = Set, fill = Correlation)) +
    geom_tile(color = "NA") +
    scale_x_discrete(limits = ionome_order,
                     labels = ionome_order) +
    scale_y_discrete(limits = rev(set_names)) + 
    scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#FF7F0E",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "") +
    geom_text(aes(label = paste0(round(Correlation, 2), stars)), size = 3.5) +
    theme_minimal(base_size = 14) +
    labs(title = trait_labels[tr]) +
    theme(plot.title = element_text(face = "bold"), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)) +
    coord_fixed()
  
  print(p)
  
  plot_list_3[[tr]] <- p
}

combined_plot_3 <- wrap_plots(plot_list_3, ncol = 1)

png(
  paste0(saveDir, scriptID, "Tottori.RIL.3traits.ion.pheno.cor.holm.png"),
  width = 12,
  height = 10,
  units = "in",
  res = 300
)

print(combined_plot_3)

dev.off()


##### 2.2. all traits #####
all_traits <- c("BL", "FSW", "GS", "NB", "NN", "NS", "TENSW", "TSW")

cor19RC <- cor19RC[1:19, all_traits]
p19RC <- p19RC[1:19, all_traits]

cor19RD <- cor19RD[1:19, all_traits]
p19RD <- p19RD[1:19, all_traits]

cor20RC <- cor20RC[1:19, all_traits]
p20RC <- p20RC[1:19, all_traits]

cor20RD <- cor20RD[1:19, all_traits]
p20RD <- p20RD[1:19, all_traits]

cor_list <- list(cor19RC, cor19RD, cor20RC, cor20RD)
p_list <- list(p19RC, p19RD, p20RC, p20RD)
set_names <- c("2019-Control","2019-Drought","2020-Control","2020-Drought")

plot_list <- list()

for(tr in all_traits){
  
  trait_df <- data.frame()
  
  p_all <- unlist(
    lapply(p_list, function(x) x[, tr])
  )
  
  p_adj_all <- p.adjust(
    p_all,
    method = "holm"
  )
  
  start <- 1
  
  for(i in 1:4){
    
    cor_vec <- cor_list[[i]][, tr]
    
    p_adj <- p_adj_all[start:(start + 18)]
    
    star <- sapply(p_adj, get_stars)
    
    df <- data.frame(
      Ionome = factor(
        rownames(cor_list[[i]]),
        levels = ionome_order
      ),
      Correlation = cor_vec,
      stars = star,
      Set = set_names[i],
      stringsAsFactors = FALSE
    )
    
    trait_df <- rbind(trait_df, df)
    
    start <- start + 19
  }
  
  p <- ggplot(trait_df, aes(x = Ionome, y = Set, fill = Correlation)) +
    geom_tile(color = "NA") +
    scale_x_discrete(limits = ionome_order,
                     labels = ionome_order) +
    scale_y_discrete(limits = rev(set_names)) + 
    scale_fill_gradient2(low = "#1F77B4", mid = "white", high = "#FF7F0E",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "") +
    geom_text(aes(label = paste0(round(Correlation, 2), stars)), size = 3) +
    theme_minimal(base_size = 14) +
    labs(title = trait_labels[tr]) +
    theme(plot.title = element_text(face = "bold"), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA)) +
    coord_fixed()
  
  print(p)
  
  plot_list[[tr]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 1)

png(
  paste0(saveDir, scriptID, "Tottori.RIL.alltraits.ion.pheno.cor.holm.png"),
  width = 12,
  height = 27,
  units = "in",
  res = 300
)

print(combined_plot)

dev.off()



