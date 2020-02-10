install.packages("readxl")
install.packages("tidyverse")
install.packages("gridExtra")
install.packages("RColorBrewer")

library(gridExtra)
library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


getwd()
setwd("/Users/arielleaiken/github/src/arielle")

DF1 <- read_csv("heatmapDF.csv")
colnames(DF1) <- c("row", "AAposition", "WT", "Mutation", "Radicicol_Score", "No_Radicicol_Score")

DF1["No_Radicicol_Score"] = DF1["No_Radicicol_Score"]*-1
DF1["Radicicol_Score"] = DF1["Radicicol_Score"]*-1

radicicol <- ggplot(data=DF1, aes(x=AAposition,y=Mutation,fill=Radicicol_Score)) +
         geom_tile() + theme(legend.position = "none") + ggtitle("With Radicicol") +
  scale_fill_distiller(type = "div", palette = "RdBu") 
  
#radicicol

no_radicicol <- ggplot(data=DF1, aes(x=AAposition,y=Mutation,fill=No_Radicicol_Score)) +
  geom_tile()  + ggtitle("Without Radicicol") + theme(legend.position = "none") +
  scale_fill_distiller(type = "div", palette = "RdBu")

#no_radicicol

grid.arrange(radicicol, no_radicicol, nrow = 1)

