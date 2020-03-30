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
setwd("/Users/arielleaiken/github/src/arielle")


#read dataframe
DF1 <- read_csv("heatmapDF.csv")
DF1 <-na.omit(DF1)


#rename columns
colnames(DF1) <- c("row", "AAposition", "WT", "Mutation", "Radicicol_Score", "No_Radicicol_Score")

colnames(DF1)
DF1



#------------------------------------------------------------------------
#Rescaling method 1. Rescale both HSP90 active and HSP90 inactive columns 
#separately: 
#------------------------------------------------------------------------

#create new rescaled dataframe with range [-1, +1]
rescale(DF1$No_Radicicol_Score, to = c(-1, 1), 
        from = range(DF1$No_Radicicol_Score, finite = TRUE))



DF1_rescaled <- as.data.frame(cbind(DF1$row, 
                                    DF1$AAposition, 
                                    DF1$WT, 
                                    DF1$Mutation,
                                    DF1$Radicicol_Score, 
                                    DF1$No_Radicicol_Score,
                                    rescale(DF1$Radicicol_Score, to = c(-1, 1),
                                            from = range(DF1$Radicicol_Score, finite = TRUE)),
                                    rescale(DF1$No_Radicicol_Score, to = c(-1, 1),
                                            from = range(DF1$No_Radicicol_Score, finite = TRUE))))

DF1_rescaled
colnames(DF1_rescaled) <- c("row","position","wildtype","mutation",
                            "HSP90_inhibited","HSP90_active","HSP90_inhibited_rescaled",
                           "HSP90_active_rescaled")
colnames(DF1_rescaled)
DF1_rescaled$HSP90_active_rescaled

#force type to numeric
DF1_rescaled$HSP90_inhibited <- as.numeric(as.character(DF1_rescaled$HSP90_inhibited))
DF1_rescaled$HSP90_inhibited_rescaled <- as.numeric(as.character(DF1_rescaled$HSP90_inhibited_rescaled))
DF1_rescaled$HSP90_active <- as.numeric(as.character(DF1_rescaled$HSP90_active))
DF1_rescaled$HSP90_active_rescaled <- as.numeric(as.character(DF1_rescaled$HSP90_active_rescaled))
DF1_rescaled$position <- as.numeric(as.character(DF1_rescaled$position))

                                                 
#check that we no longer get the "levels" message (levels indicate non numeric values)                                                 
DF1_rescaled$HSP90_active_rescaled                                        
 

#create HI (HSP90 Inhibited) and HA (HSP90 Active) heatmaps                                          
HI <- ggplot(data=DF1_rescaled, 
       aes(x=position,y=mutation,fill=HSP90_inhibited_rescaled)) +
geom_tile() + 
theme(legend.position = "none", panel.grid.major = element_blank(), 
panel.grid.minor = element_blank()) +
ggtitle("HSP90 Inhibited") +  xlab("position") +
scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
scale_x_discrete(limits=seq(280,520,by=10))

HA <- ggplot(data=DF1_rescaled, 
              aes(x=position,y=mutation,fill=HSP90_active_rescaled)) +
  geom_tile() + 
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("HSP90 Active") +  xlab("position") +
  scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
  scale_x_discrete(limits=seq(280,520,by=10))

#arrange heatmaps onto one panel
grid.arrange(HI,HA,nrow=2)

#save csv here, unhash to run once, dont re-run. 
#write.csv(DF1_rescaled,"/Users/arielleaiken/github/src/arielle/rescaled.csv", row.names = FALSE)


#create subtracted (HA-HI) column for rescaled values and add to dataframe. 
#Subtracted values stored as 'S'
S <- DF1_rescaled$HSP90_active_rescaled - DF1_rescaled$HSP90_inhibited_rescaled

DF1_rescaled2 <- as.data.frame(cbind(DF1_rescaled, S))


#Create heatmap with subtracted scores. Heatmap stored as 'subtracted'
subtracted <- ggplot(data=DF1_rescaled, 
                     aes(x=position,y=mutation,fill=S)) +
  geom_tile() + 
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("Active-Inhibited") +  xlab("Amino Acid Position") +
  scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
  scale_x_discrete(limits=seq(280,520,by=10))
subtracted

#combine all three heatmaps
grid.arrange(HSI,HSA,subtracted, nrow=3)




#coeffcients for fitting a line for residual plot
coefs <- coef(lm(HSP90_inhibited_rescaled~HSP90_active_rescaled, 
        data= DF1_rescaled))
coefs
y=coefs[2]*DF1_rescaled$HSP90_active_rescaled + coefs[1]

ggplot(data=DF1_rescaled,aes(x=HSP90_active_rescaled,y=HSP90_inhibited_rescaled)) +
  geom_point(alpha=0.4,color=I("blue")) + 
  geom_smooth(color="green", size=.3) +
  xlab("HSP90 Active") +ylab("HSP90 Inactive") +ggtitle("Residual Plot") +
  geom_abline(slope=coefs[2],intercept=coefs[1], alpha=.4)

#--------------------------------------------------------------------------#
#Rescaling method 2. Combine hsp90 active and hsp90 inactive scores 
#into one column, THEN rescale together. 
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
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("active-inhibited") +  xlab("position") +
  scale_fill_gradient2(low = "red4", mid = "white",high = "royalblue4") +
  scale_x_discrete(limits=seq(280,520,by=10)) 

grid.arrange(rescaled_together,sub_plot,nrow=2)

