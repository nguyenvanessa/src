##################################################################
##                           20200409                           ##
##               Add secondary structure to plots               ##
##################################################################


#set working dir and  load packages
setwd("/Users/arielleaiken/github/new_src_stuff/src/arielle")
library(gridExtra)
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(grDevices)
library(tidyr)


#open csv file to be used for plotting
new_master <- read.csv("R_working_master_20200409.csv")


##################################################################
##                          Section 1:                          ##
##                 activity~SASA scores ggplots                 ##
##################################################################


#KINASE DOMAIN ACTIVE
SASA1_active <- ggplot(data=new_master,
                aes(x=CD_active_SASA,
                y=HSP90_active,
                color=structure)) +
          geom_point(alpha=.3) +
                ggtitle("activity~SASA kinase open",
                subtitle ="DMSO based activity scores") +
                xlab("SASA score for kinase domain in open conformation") +
                ylab("activity")



SASA1_active + facet_grid(structure~.)


#KINASE DOMAIN INACTIVE
SASA1_inactive <- ggplot(data=new_master,
                  aes(x=CD_inactive_SASA,
                  y=HSP90_active,
                  color=structure)) +
              geom_point(alpha=.3) +
                  ggtitle("activity~SASA kinase closed",
                  subtitle ="DMSO based activity scores") +
                  xlab("SASA score for kinase domain in closed conformation") +
                  ylab("activity")

SASA1_inactive +facet_grid(structure~.)




##################################################################
##                          Section 2:                          ##
##             activity~conservation scores ggplots             ##
##################################################################


#TREND LINE FOR RAW SCORES
coefs_DMSO <-coef(lm(HSP90_active~conservation_EA, data= new_master))
trend_line_DMSO<- coefs_DMSO[2]*new_master$conservation_EA + coefs_DMSO[1]

#RAW SCORES ACTIVITY ~ CONSERVATION SCORES 
ggplot(data=new_master, 
       aes(x=conservation_EA,y=HSP90_active,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  geom_smooth(data=new_master, color="black",
              aes(x=conservation_EA,y=trend_line_DMSO)) 

#FACETED VERSION OF RAW SCORES ACTIVITY ~ CONSERVATION SCORES 
ggplot(data=new_master, 
       aes(x=conservation_EA,y=HSP90_active,color=structure)) +
  geom_point(alpha=.4) +
      xlab("conservation score") +ylab("activity") +
      ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  facet_grid(structure~.) +
      theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-4,4))





#DIFFSEL BASED ACTIVITY SCORES ~ CONSERVATION SCORES
ggplot(data=new_master, 
       aes(x=conservation_EA,y=active_minus_inhib,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from diffsel scores") +
  facet_grid(structure~.) +
  theme(legend.position = "none") +
  coord_cartesian(
    ylim = c(-4,4))
 



#RESIDUAL BASED ACTIVITY SCORES ~ CONSERVATION SCORES
ggplot(data=new_master, 
       aes(x=conservation_EA,y=residuals,color=structure)) +
  geom_point(alpha=.4) +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from residual scores") +
  facet_grid(structure~.) +
  theme(legend.position = "none") 


#how many points are in each facet?
b<- new_master[which(new_master$structure=='beta') , ]
a<- new_master[which(new_master$structure=='helix') , ]
n<- new_master[which(new_master$structure=='none') , ]

length(b$position)
length(a$position)
length(n$position)

new_master <- arrange(new_master, position)
new_master


#write.csv(new_master,"/Users/arielleaiken/github/new_src_stuff/src/arielle/R_working_master_20200409.csv", row.names = FALSE)


colnames(new_master)
new_master[which(new_master$CD_active_SASA >200)  ,]

tail(arrange(new_master[which( new_master$structure=="sheet") ,] , 
        CD_active_SASA), 20)
          
                
tail(arrange(new_master[which( new_master$structure=="helix") ,] , 
             conservation_EA), 60)




#compare these two metrics SASA~Conservation:

coefs <-coef(lm(CD_active_SASA~conservation_EA, data=new_master))

coefs[1]
coefs[2]

line<- coefs[2]*new_master$conservation_EA + coefs[1]


ggplot(data=new_master, 
       aes(x=conservation_EA, y=CD_active_SASA, 
           alpha=I(.6), size=I(3), color=HSP90_active)) +
  geom_jitter() +
  geom_point() + 
    geom_smooth(data=new_master,
                              color="red",
                             aes(x=conservation_EA,
                                y=line,
                              size=.3)) +
  ggtitle("solvent accessibility ~ conservation ") +
  xlab("conservation score") + ylab("SASA scores") 



#Same SASA~conservation plot but faceted by secondary structure,
#colored by DMSO activity score (HSP90 active)
ggplot(data=new_master, 
       aes(x=conservation_EA, y=CD_active_SASA, 
           alpha=I(.6), size=I(3), color=HSP90_active)) +
  geom_jitter() +
  geom_point() + 
#  geom_smooth(data=new_master,
 #                            color="red",
  #                           aes(x=conservation_EA,
   #                              y=line,
    #                            size=.3)) +
  ggtitle("solvent accessibility ~ conservation ") +
  xlab("conservation score") + ylab("SASA scores") +
  facet_grid(structure~.)


