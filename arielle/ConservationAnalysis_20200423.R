###########################################################################
###########################################################################
###                                                                     ###
###                          APRIL 23 2020 AA                           ###
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
SI<-filter(working, client_status=="Strong_Inhibited")
MI<-filter(working, client_status=="Moderate_Inhibited")
WI<-filter(working, client_status=="Weak_Inhibited")
NC<-filter(working, client_status=="nonclient")
WD<-filter(working, client_status=="Weak_Dependent")
MD<-filter(working, client_status=="Moderate_Dependent")
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


##################################################################
##   make an empty matrix (all entries are set to 0 to start)   ##
##################################################################

SD <- arrange(SD,position)
levels(as.factor(a(SD$mutation)))
levels(as.factor(SD$position))

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


##################################################################
##                      fill in the matrix                      ##
##################################################################
#This 'for' loops fills in each position in the matrix according 
#to the residual scores for each amino acid at positions in 
#'Strong Dependend' clients 

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

#convert the dataframe back to a matrix
LogoMatrix <- as.matrix(LogoDF)

##################################################################
##                      draw the logo plot                      ##
##################################################################


ggseqlogo(-LogoMatrix, 
            method='prob', 
              col_scheme='auto') + 
          ylab('probability') +
  ggtitle("Strong Dependent Clients")


