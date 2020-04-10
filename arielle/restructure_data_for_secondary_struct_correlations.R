setwd("/Users/arielleaiken/github/new_src_stuff/src/arielle")

library(gridExtra)
library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(grDevices)
library(tidyr)
library(bio3d)

##################################################################
##################################################################
##                          Section 1:                          ##
##          get secondary structure data from PDB               ##
##################################################################
##################################################################


#read in src PDB file
srcFULL<- read.pdb("2src.pdb")

#start and end positions for alpha helicies: 
helix <- as.data.frame(cbind(srcFULL$helix$start,
                    srcFULL$helix$end))
colnames(helix) <- c("start", "end")

#start and end positions for beta sheets:
sheet <- as.data.frame(cbind(srcFULL$sheet$start,
                    srcFULL$sheet$end))
colnames(sheet) <- c("start","end")

helix
sheet

####################################################################
##                           Section 2:                           ##
##          translate helix & sheet start/end data into           ##
##         factors that can go into a column for point muts data  ##
####################################################################



#combine helix and sheet dataframes
hs<-as.data.frame(rbind(helix,sheet))


#annotate secondary structure type 
type <- c(rep("helix",times=length(helix$start)),
rep("beta",times=length(sheet$start)))
struct_data <- as.data.frame(cbind(hs,type))

#order data by 'start' value 
struct_data <- arrange(struct_data,start)

#get length of each secondary structure
length <- struct_data$end - struct_data$start
struct_data<- as.data.frame(cbind(struct_data,length))

struct_data
#subset to just kinase domain positions
struct_data<- struct_data[13:32,]


#make a dataframe with start, end, type, and length of 'none' secondary structure
#combine with helix/sheet dataframe


struct_data$end

n_start <- struct_data$end[1:19]+1
n_end <- struct_data$start[2:20]-1


none_length <- n_end-n_start
nones <- as.data.frame(cbind(n_start,n_end, rep("none",length(n_start)),
                    none_length))

nones[9,] <- NA
nones <- na.omit(nones)
nones

colnames(nones) <- colnames(struct_data)
struct_data<- as.data.frame(rbind(struct_data,nones))
struct_data <- arrange(struct_data,start)
struct_data



#fix the begining and end of data frame so it accounts for position 270:519
#aka kinase domain 

struct_data[1,1] <- "270" #change start to 270
struct_data[1,4] <- "2" #change length of first row to 2

addin <- as.data.frame(cbind(518,519,"none",1))
colnames(addin) <- colnames(struct_data)
struct_data <- as.data.frame(rbind(struct_data,addin))

#create a list of types that should be equal to the length of kinase domain
just_types <- rep(struct_data$type,struct_data$length)

just_types
length(just_types)

struct_data

249-211

519-270 # equals 249
length(just_types) #why does this not equal 249 then? 

#because the lengths are off by 1. positions are inclusive.

struct_data$start <-as.numeric(as.character(struct_data$start))
struct_data$end <-as.numeric(as.character(struct_data$end))
struct_data$length <-as.numeric(as.character(struct_data$length))

struct_data$length <- struct_data$length+1

just_types <- rep(struct_data$type,struct_data$length)

just_types
length(just_types)

just_types <- as.data.frame(cbind(270:519,just_types))

long_struct <- as.data.frame(cbind(270:519,as.character(just_types)))
colnames(long_struct) <- c("positions","structure")

#################################################################
##                          Section 3:                         ##
##        combine structural data with master dataframe        ##
#################################################################

master <- read.csv("R_working_master_20200406.csv")



repetitive <- c()
position_in_master <- master$position
for(i in position_in_master){
  repetitive<-
    append(repetitive,
           as.character(
             long_struct[which(long_struct$position==i),]$structure))
}
length(repetitive)


new_master <- as.data.frame(cbind(master,repetitive))

colnames(new_master)[11] <- "structure"
head(new_master)


##################################################################
##                          Section 4:                          ##
##                            ggplots                           ##
##################################################################


#SASA plots


SASA1_active <- ggplot(data=new_master,
                       aes(x=CD_active_SASA,
                       y=HSP90_active,
                       color=structure)) +
                geom_point(alpha=.3) +
                ggtitle("activity~SASA kinase open",
subtitle ="DMSO based activity scores") +
                xlab("SASA score for kinase domain in open conformation") +
                ylab("activity")



SASA1_active

SASA1_inactive <- ggplot(data=new_master,
                       aes(x=CD_inactive_SASA,
                           y=HSP90_active,
                           color=structure)) +
  geom_point(alpha=.3) +
  ggtitle("activity~SASA kinase closed",
          subtitle ="DMSO based activity scores") +
  xlab("SASA score for kinase domain in closed conformation") +
  ylab("activity")
SASA1_inactive

grid.arrange(SASA1_active,SASA1_inactive)


#conservation plots

coefs_DMSO <-coef(lm(HSP90_active~conservation_EA, data= new_master))
trend_line_DMSO<- coefs_DMSO[2]*new_master$conservation_EA + coefs_DMSO[1]


ggplot(data=new_master, 
       aes(x=conservation_EA,y=HSP90_active,color=structure)) +
  geom_point() +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  geom_smooth(data=new_master, color="black",
            aes(x=conservation_EA,y=trend_line_DMSO)) 

#faceted version
ggplot(data=new_master, 
       aes(x=conservation_EA,y=HSP90_active,color=structure)) +
  geom_point() +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from raw yeast growth scores (no radicicol)") +
  #geom_smooth(data=new_master, color="black",
   #           aes(x=conservation_EA,y=trend_line_DMSO)) +
  facet_grid(structure~.) +
  theme(legend.position = "none") 
  



#diffsel
ggplot(data=new_master, 
       aes(x=conservation_EA,y=active_minus_inhib,color=structure)) +
  geom_point() +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from diffsel scores") +
  #geom_smooth(data=new_master, color="black",
  #           aes(x=conservation_EA,y=trend_line_DMSO)) +
  facet_grid(structure~.) +
  theme(legend.position = "none") 


#residual
ggplot(data=new_master, 
       aes(x=conservation_EA,y=residuals,color=structure)) +
  geom_point() +
  xlab("conservation score") +ylab("activity") +
  ggtitle("activity~conservation", 
          subtitle="activity from residual scores") +
  #geom_smooth(data=new_master, color="black",
  #           aes(x=conservation_EA,y=trend_line_DMSO)) +
  facet_grid(structure~.) +
  theme(legend.position = "none") 


#how many points are in each facet?
b<- new_master[which(new_master$structure=='sheet') , ]
a<- new_master[which(new_master$structure=='helix') , ]
n<- new_master[which(new_master$structure=='none') , ]

length(b$position)
length(a$position)
length(n$position)


tail(arrange(a, conservation_EA), 30)

