
#------------------------------------------------------------------------
#open packages needed, set workind dir, read data into R:
#------------------------------------------------------------------------
library(gridExtra)
library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(grDevices)
library(nlme)
library(scales)
library(tidyr)

getwd()
#setwd("/Users/arielleaiken/github/src/arielle")


#read dataframe, omit NA values
DF1 <- read_csv("heatmapDF.csv")
DF1 <-na.omit(DF1)

#rename columns
colnames(DF1) <- c("row", "AAposition", "WT", "Mutation", "Radicicol_Score", "No_Radicicol_Score")

colnames(DF1)

#------------------------------------------------------------------------

#construct a LONG style data frame with both hsp90 active and hsp90 inactive
#scores combined

long_data <- gather(DF1, HSP90, Score, Radicicol_Score:No_Radicicol_Score, factor_key=TRUE)

#rename factor levels from Radicicol/noRadicicol to (HSP90) active/inactive 

levels(long_data$HSP90) <- c("inhibited","active")
long_data

#rescale the score column from -1 to +1
long_data$Score <- rescale(long_data$Score, to = c(-1, 1), 
                           from = range(long_data$Score, finite = TRUE))


#rename columns
colnames(long_data) <- c("row","position","wildtype","mutation","HSP90","score")
long_data
#plot heatmap

rescaled_together <- ggplot(data=long_data, 
                            aes(x=position,y=mutation,fill=score)) +
  geom_tile() + 
  theme(legend.position = "right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("HSP90 inhibited and active DMS") +  xlab("position") +
  scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
  scale_x_discrete(limits=seq(280,520,by=10)) +
  facet_wrap(HSP90~.,nrow=2)



#subtraction plot

active_col <- long_data[which(long_data$HSP90=='active'), ]$score

inhibited_col <- long_data[which(long_data$HSP90=='inhibited'), ]$score

sub_col <- active_col-inhibited_col

with_sub <- as.data.frame(
  cbind(
    DF1$AAposition,
    DF1$WT,
    DF1$Mutation,
    active_col,
    inhibited_col,
    sub_col
    
  )
)
colnames(with_sub) <- c("position","wildtype","mutation",
                        "HSP90_active","HSP90_inhibited","subtraction")
colnames(with_sub)

with_sub$position <- as.numeric(as.character(with_sub$position))
with_sub$subtraction <- as.numeric(as.character(with_sub$subtraction))


sub_plot<- ggplot(data=with_sub, 
                  aes(x=position,y=mutation,fill=subtraction)) +
  geom_tile() + 
  theme(legend.position = "right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("active-inhibited") +  xlab("position") +
  scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
  scale_x_discrete(limits=seq(280,520,by=10)) 

grid.arrange(rescaled_together,sub_plot,nrow=2)




