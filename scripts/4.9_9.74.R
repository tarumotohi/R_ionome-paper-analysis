###### 0. Settings ######
rm(list = ls())


##### 0.1. Load Packages #####
library(ggplot2)
library(patchwork)
library(readr)


##### 0.2. Define Parameters #####
scriptID <- "9.74."


##### 0.3. Define Directories #####
dataDir <- "data/ionome/"
saveDir <- paste0("results/ionome/RIL_F2/")
if (!dir.exists(saveDir)) {
  dir.create(path = saveDir, recursive = TRUE)
}




###### 1. Data Preparation ######
##### 1.1. Load Data #####
#9.55 f_hat_RIL
fileNamef_hat_RIL_RC <- paste0(saveDir, "9.55.model2.RC.f_hat_RIL.RDS")
f_hat_RIL_RC <- readRDS(fileNamef_hat_RIL_RC)

fileNamef_hat_RIL_RD <- paste0(saveDir, "9.55.model2.RD.f_hat_RIL.RDS")
f_hat_RIL_RD <- readRDS(fileNamef_hat_RIL_RD)

#9.67 Pos_hat_RIL Model 2 y = Y + Kpos_2 + u + e
fileNameKpos_hat_RIL_RC <- paste0(saveDir, "9.67.model2.RC.Pos_hat_RIL.RDS")
Kpos_hat_RIL_RC <- readRDS(fileNameKpos_hat_RIL_RC)

fileNameKpos_hat_RIL_RD <- paste0(saveDir, "9.67.model2.RD.Pos_hat_RIL.RDS")
Kpos_hat_RIL_RD <- readRDS(fileNameKpos_hat_RIL_RD)

#map
fieldMap <- data.frame(read_csv(paste0("data/", "Tottori.RIL.PlotPosition.csv")))

rownames(fieldMap) <- paste0(
  "E71N",
  sprintf("%03d", fieldMap$plot)
)

fieldMap$id <- substr(rownames(fieldMap), 1, 7)




###### 2. plot ######
##### 2.1.  all #####
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

treatment_names <- c(
  RC = "Control",
  RD = "Drought"
)

effects <- list(
  Ionome = list(
    RC = f_hat_RIL_RC,
    RD = f_hat_RIL_RD
  ),
  Position = list(
    RC = Kpos_hat_RIL_RC,
    RD = Kpos_hat_RIL_RD
  )
)

model_names <- c(
  Ionome = "Ionome",
  Position = "Distance"
)

traits <- names(f_hat_RIL_RC)
years <- c("19", "20")

for(trait in traits){
  
  plots <- list()
  
  for(model in names(effects)){
    
    for(treatment in names(effects[[model]])){
      
      for(year in years){
        
        effect <- effects[[model]][[treatment]][[trait]]
        
        
        effect_year <- effect[
          substr(rownames(effect), 9, 10) == year,
          ,
          drop = FALSE
        ]
        
        
        id <- substr(
          rownames(effect_year),
          1,
          7
        )
        
        
        df <- data.frame(
          id = id,
          value = as.numeric(effect_year[,1])
        )
        
        
        df_merged <- merge(
          df,
          fieldMap,
          by = "id"
        )
        
        
        p <- ggplot(
          df_merged,
          aes(
            x = x.real,
            y = y.real,
            fill = value
          )
        ) +
          geom_tile(
            color = "grey50",
            linewidth = 0.15
          ) +
          scale_fill_gradient2(low = "#1F77B4", 
                               mid = "white", 
                               high = "#FF7F0E",
                               midpoint = 0
          ) +
          scale_x_reverse() +
          scale_y_reverse() +
          coord_fixed() +
          theme_minimal() +
          ggtitle(
            paste0(
              model_names[model],
              " / ",
              " 20",
              year,
              "-",
              treatment_names[treatment]
            )
          ) +
          theme(
            plot.title = element_text(
              hjust = 0.5,
              size = 15,
              face = "bold"
            ),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.3, "cm"),
            legend.text = element_text(size=8)
          )
        
        
        plots[[paste(
          model,
          treatment,
          year,
          sep="_"
        )]] <- p
        
      }
    }
  }
  
  trait_plot <- wrap_plots(
    plots,
    ncol = 2
  ) +
    plot_annotation(
      title = trait_labels[trait]
    ) +
    plot_layout(
      widths = c(1,1),
      heights = c(1,1,1,1)
    ) &
    theme(
      plot.margin = margin(5,5,5,5),
      plot.title = element_text(
        size = 15)
    ) 
  
  pdf(
    file = paste0(
      saveDir,
      scriptID,
      "Effect_map_",
      trait,
      ".pdf"
    ),
    width = 18,
    height = 8
  )
  
  print(trait_plot)
  
  dev.off()
  
}

##### 2.2.  2019 FSW #####
trait <- "FSW"
year <- "19"

plots <- list()

for(model in names(effects)){
  
  for(treatment in names(effects[[model]])){
    
    effect <- effects[[model]][[treatment]][[trait]]
    
    effect_year <- effect[
      substr(rownames(effect), 9, 10) == year,
      ,
      drop = FALSE
    ]
    
    id <- substr(
      rownames(effect_year),
      1,
      7
    )
    
    df <- data.frame(
      id = id,
      value = as.numeric(effect_year[,1])
    )
    
    df_merged <- merge(
      df,
      fieldMap,
      by = "id"
    )
    
    p <- ggplot(
      df_merged,
      aes(
        x = x.real,
        y = y.real,
        fill = value
      )
    ) +
      geom_tile(
        color = "grey50",
        linewidth = 0.15
      ) +
      scale_fill_gradient2(
        low = "#1F77B4",
        mid = "white",
        high = "#FF7F0E",
        midpoint = 0
      ) +
      scale_x_reverse() +
      scale_y_reverse() +
      coord_fixed() +
      theme_minimal() +
      ggtitle(
        paste0(
          model_names[model],
          " / ",
          "2019-",
          treatment_names[treatment]
        )
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          size = 15,
          face = "bold"
        ),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(size = 8)
      )
    
    plots[[paste(
      model,
      treatment,
      sep = "_"
    )]] <- p
    
  }
}

FSW_2019_plot <- wrap_plots(
  plots,
  ncol = 1
) +
  plot_annotation(
    title = trait_labels[trait]
  )

pdf(
  file = paste0(saveDir, scriptID, "Effect_map_FSW_2019.pdf"
  ),
  width = 10,
  height = 8
)

print(FSW_2019_plot)

dev.off()



