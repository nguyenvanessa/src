

###########################################################################
###########################################################################
###                                                                     ###
###                          APRIL 24 2020 AA                           ###
###                        CONSERVATION ANALYSIS                        ###
###                              A. AIKEN                               ###
###                                                                     ###
###########################################################################
###########################################################################

#################################################################
##                        load packages                        ##
#################################################################
library(ggseqlogo)
library(bio3d)
library(ggplot2)
library(ggExtra)
library(seqinr)
library(grDevices)
library(readr)
library(gridExtra)
library(dplyr)
##################################################################
##                read & structure required data                ##
##################################################################

working <- read.csv("R_working_master_20200422.csv")

#adjust factor levels
working$client_status <- 
  factor(working$client_status, 
         levels = c( "Strong_Inhibited","Moderate_Inhibited", 
                     "Weak_Inhibited","nonclient","Weak_Dependent",
                     "Moderate_Dependent", "Strong_Dependent"
         ))

#subset clients into separate dataframes 
SD<-filter(working, client_status=="Strong_Dependent")


############################################################################
############################################################################
###                                                                      ###
###                          SEQUENCE LOGO PLOT                          ###
###                                                                      ###
############################################################################
############################################################################


##################################################################
##              accepted input format is Matrices:              ##
##                  a position frequency matrix                 ##
##                  where the row is the letter                 ##
##                  and column is the position                  ##
##      Note: the matrix must be row named with the letter      ##
##################################################################


#################################################################
##  Replace SD with working dataframe --> includes all points  ##
#################################################################


SD <- working

levels(as.factor(a(SD$mutation)))
levels(as.factor(SD$position))
SD <- arrange(SD,mutation)

LogoMatrix <- matrix( 
  nrow=length(levels(as.factor(a(SD$mutation)))),
  ncol =length(levels(as.factor(SD$position))),
  dimnames=list(levels(as.factor(a(SD$mutation))))
)

LogoDF <- as.data.frame(LogoMatrix)
colnames(LogoDF) <- levels(as.factor(SD$position))
LogoMatrix <- as.matrix(LogoDF)
LogoMatrix

#set all values to 0
LogoDF[c(1:length(LogoDF))] <- 0

#change 3 letter AA code to 1 letter AA code
SD$mutation <- a(SD$mutation)
SD$wildtype <- a(SD$wildtype)

SD





#Fill in the matrix
for(i in (levels(as.factor(
  filter(SD, (mutation)=="A")$position
)))) {
  LogoDF["A",i] <- 
    SD[which(
      SD$mutation=='A' &
        SD$position==as.numeric(i)),"residual"]
}

for(i in (levels(as.factor(
  filter(SD, (mutation)=="C")$position
)))) {
  LogoDF["C",i] <- 
    SD[which(
      SD$mutation=='C' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="D")$position
)))) {
  LogoDF["D",i] <- 
    SD[which(
      SD$mutation=='D' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="E")$position
)))) {
  LogoDF["E",i] <- 
    SD[which(
      SD$mutation=='E' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="F")$position
)))) {
  LogoDF["F",i] <- 
    SD[which(
      SD$mutation=='F' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="H")$position
)))) {
  LogoDF["H",i] <- 
    SD[which(
      SD$mutation=='H' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="I")$position
)))) {
  LogoDF["I",i] <- 
    SD[which(
      SD$mutation=='I' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="K")$position
)))) {
  LogoDF["K",i] <- 
    SD[which(
      SD$mutation=='K' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="L")$position
)))) {
  LogoDF["L",i] <- 
    SD[which(
      SD$mutation=='L' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="M")$position
)))) {
  LogoDF["M",i] <- 
    SD[which(
      SD$mutation=='M' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="N")$position
)))) {
  LogoDF["N",i] <- 
    SD[which(
      SD$mutation=='N' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="P")$position
)))) {
  LogoDF["P",i] <- 
    SD[which(
      SD$mutation=='P' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="Q")$position
)))) {
  LogoDF["Q",i] <- 
    SD[which(
      SD$mutation=='Q' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="R")$position
)))) {
  LogoDF["R",i] <- 
    SD[which(
      SD$mutation=='R' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="S")$position
)))) {
  LogoDF["S",i] <- 
    SD[which(
      SD$mutation=='S' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="T")$position
)))) {
  LogoDF["T",i] <- 
    SD[which(
      SD$mutation=='T' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="V")$position
)))) {
  LogoDF["V",i] <- 
    SD[which(
      SD$mutation=='V' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="W")$position
)))) {
  LogoDF["W",i] <- 
    SD[which(
      SD$mutation=='W' &
        SD$position==as.numeric(i)),"residual"]
}
for(i in (levels(as.factor(
  filter(SD, (mutation)=="Y")$position
)))) {
  LogoDF["Y",i] <- 
    SD[which(
      SD$mutation=='Y' &
        SD$position==as.numeric(i)),"residual"]
}



LogoDF
LogoMatrix <- as.matrix(LogoDF)
LogoMatrix

ggseqlogo(abs(LogoMatrix), 
          method='prob', 
          col_scheme='auto') + 
  ylab('probability')  + facet_grid()
  


ggseqlogo(abs(LogoMatrix[,seq(from=1, to =249, by=10)]), 
          method='prob', 
          col_scheme='auto') + 
  ylab('probability?')  + ggtitle("every 10 positions",
                         subtitle = "letter height based on residual" )


LogoMatrix

##################################################################
##                           april 24                           ##
##################################################################

##################################################################
##                    conservation bar plots                    ##
##################################################################


#filter and store clients into separate variables
SI<-filter(working, client_status=="Strong_Inhibited")
MI<-filter(working, client_status=="Moderate_Inhibited")
WI<-filter(working, client_status=="Weak_Inhibited")
NC<-filter(working, client_status=="nonclient")
WD<-filter(working, client_status=="Weak_Dependent")
MD<-filter(working, client_status=="Moderate_Dependent")
SD<-filter(working, client_status=="Strong_Dependent")
NC_Cmean <- mean(NC$conservation) #non
SD_Cmean<- mean(SD$conservation) # strong dependent
SI_Cmean<- mean(SI$conservation) # strong inhibited
MD_Cmean <- mean(MD$conservation) #moderate dependent
MI_Cmean <- mean(MI$conservation) #moderate inhibited
WD_Cmean <- mean(WD$conservation) #weak depen
WI_Cmean <- mean(WI$conservation) #weak inhib

all_Cmean <- mean(working$conservation)

meancons <- as.data.frame(rbind(NC_Cmean,SD_Cmean,SI_Cmean,
                                MD_Cmean,MI_Cmean, WD_Cmean,
                                WI_Cmean, all_Cmean))


meancons[2] <- rownames(meancons)
colnames(meancons) <- c("mean_conservation", "client_status")
meancons


#plot
bars <- ggplot(meancons,  aes(x=client_status,y=mean_conservation,
                      fill=client_status)) +
  geom_col() + theme_bw() +
  scale_x_discrete(labels=c("all clients", "moderate dependent", "moderate inhibited",
                            "nonclients", "strong dependent", "strong inhibited",
                            "weak dependent", "weak inhibited")) +
  theme(axis.text.x = element_text(color="black", size=10, angle=90),
        legend.position = "none") +
  ggtitle("mean conservation of mutated positions for hsp90 clients")+
  scale_fill_viridis_d() 

#filtered again for just strong clients
meancons

bars_filtered <- ggplot(meancons[c(2,3),],  aes(x=client_status,y=mean_conservation,
                      fill=client_status)) +
  geom_col() + theme_bw() +
  scale_x_discrete(labels=c( "strong dependent", "strong inhibited")) +
  theme(axis.text.x = element_text(color="black", size=10),
        legend.position = "none") +
  ggtitle("mean conservation of mutated positions for hsp90 clients") +
  scale_fill_viridis_d(begin=.5, end=.7)


grid.arrange(bars,bars_filtered, nrow=2)






