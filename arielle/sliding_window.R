############################################################################
############################################################################
###                                                                      ###
###                               20200410                               ###
###             SLIDING WINDOW CONSERVATION SCORE GENERATION             ###
###                                                                      ###
############################################################################
############################################################################


#set up workspace
setwd("/Users/arielleaiken/github/new_src_stuff/src/arielle")

library(gridExtra)
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(grDevices)
library(tidyr)
library(bio3d)
library(slider)
library(ggplot2)
library(evobiR)

DF <- read.csv("R_working_master_20200409.csv")



#################################################################
##                          Section1                           ##
##         generate sliding window conservation scores         ##
##            using evobiR function SlidingWindow()            ##
#################################################################

SW <-  SlidingWindow("mean", DF$conservation_EA, window=12, step=1)

#bind sliding window list to master data frame to use in ggplot. 
ggDF <- as.data.frame(cbind(DF[1:3487,], SW))

#rename new sliding window column
ggDF
colnames(ggDF) <- c(colnames(ggDF)[1:11], "conservation_sliding")
colnames(ggDF)


##################################################################
##                         Section 2                            ##
##              DMSO activity ~ conservation plots              ##
##################################################################


#activity~conservation using
#sliding window conservation scores 
#DMSO scores == HSP90_active 
DMSO_sw <- ggplot(data=ggDF, 
             aes(x=conservation_sliding,y=HSP90_active,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  facet_grid(structure~.) +
  theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-3,3))


#activity~conservation using
#DMSO scores == HSP90_active 
#compare to old (non sliding window) conservation scores
DMSO_orig <- ggplot(data=ggDF, 
               aes(x=conservation_EA,y=residuals,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  facet_grid(structure~.) +
  theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-3,3))




#compare the two versions using grid.arrage() 
#puts both plots on one panel
grid.arrange(DMSO_sw, DMSO_orig)





##################################################################
##                          Section 3                           ##
##             Residual scores ~ conservation plots             ##
##################################################################



#activity~conservation using
#sliding window conservation scores 
Residuals_sw <- ggplot(data=ggDF, 
       aes(x=conservation_sliding,y=residuals,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from residual scores") +
  facet_grid(structure~.) +
  theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-1,1))






#activity~conservation using
#compare to old (non sliding window) conservation scores
Residuals_orig <- ggplot(data=ggDF, 
       aes(x=conservation_EA,y=residuals,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from residual scores") +
  facet_grid(structure~.) +
  theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-1,1))

orig


#compare the two versions using grid.arrage() 
#puts both plots on one panel
grid.arrange(Residuals_sw, Residuals_orig)

