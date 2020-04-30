
#################################################################
##                        April 21 2020                        ##
##                    Conservation Analysis                    ##
##                    Using Residual Scores                    ##
#################################################################



#################################################################
##                        load packages                        ##
#################################################################
library(bio3d)
library(ggplot2)
library(ggExtra)
library(seqinr)
library(grDevices)
library(readr)
library(gridExtra)
library(dplyr)
#################################################################
##                          read data                          ##
#################################################################

#VN's residual scores
VNr <- read.csv("/Users/arielleaiken/github/src4/src/utils/residuals.csv")

#EA's conservation scores
EAc <- read.csv("/Users/arielleaiken/github/src4/src/arielle/data/EAconservation.csv")

#my data frame (containing my residual scores)
AAdf <- read.csv("/Users/arielleaiken/github/src4/src/arielle/data/R_working_master_20200415.csv")

############################################################################
############################################################################
###                                                                      ###
###                           RESTRUCTURE DATA                           ###
###                                                                      ###
############################################################################
############################################################################

#prelim adjustments
VNr <- VNr[3:length(VNr$X),]
length(VNr$X) == length(AAdf$position) 
split_this <- as.character(VNr[,1])
#split_this to be split into three columns: position, wildtype, and mutation

#################################################################
##                       position column                       ##
#################################################################
position <- c()
for(i in c(1:3499)) {
  position <- append(position, 
                     str_sub(split_this[i], 2, 4))
}
position <- as.numeric(as.character(position))

#################################################################
##                       wildtype column                       ##
#################################################################
wildtype <- c()
for(i in c(1:3499)) {
  wildtype <- append(wildtype, 
                     str_sub(split_this[i], 1, 1))
}
wildtype  <- aaa(wildtype)

#################################################################
##                       mutation column                       ##
#################################################################
mutation <- c()
for(i in c(1:3499)) {
  mutation <- append(mutation, 
                     str_sub(split_this[i], 5, 5))
}
mutation <- aaa(mutation)

#################################################################
##            combine position, wildtype, mutation,            ##
##            residual scores, conservation scores,            ##
##          and secondary structure into on dataframe          ##
##                      called 'working'                       ##
#################################################################

working <- as.data.frame(cbind(position,wildtype,mutation,VNr$envision_scaled_resids))


#this takes the positions from EA's conservation data
#and uses then to match conservation scores to the 
#positions in the 'working' dataframe:
repetitive <- c()
position_in_working <- working$position
for(i in position_in_working){
  repetitive<-
       append(repetitive,
               EAc[which(EAc$residue==i),]$normalized_conservation
    )
}

working <- as.data.frame(cbind(
  position,wildtype,mutation,
  VNr$envision_scaled_resids, 
  repetitive))

colnames(working) <- c("position","wildtype","mutation",
                       "residual", "conservation")

working <- arrange(working, conservation)

AAdf <- arrange(AAdf, conservation_EA)


working <- as.data.frame(cbind(
           working, AAdf$structure))

working <- arrange(working, position)

colnames(working) <- c("position","wildtype","mutation",
                       "residual", "conservation", "structure")

working$residual <- as.numeric(as.character(working$residual))
working$conservation <- as.numeric(as.character(working$conservation))
working$position <- as.numeric(as.character(working$position))

#save working dataframe: 
#write.csv(working, "/Users/arielleaiken/github/src4/src/arielle/R_working_master_20200421.csv", row.names = FALSE)



###########################################################################
###########################################################################
###                                                                     ###
###                                PLOTS                                ###
###                                                                     ###
###########################################################################
###########################################################################


##################################################################
##                 abs(residuals)~conservation                  ##
##################################################################

RESxCONS <- ggplot(data=working,
       aes(x=conservation,y=abs(residual), 
           color=structure)) +
  geom_point() 
RESxCONS 

RESxCONS + coord_cartesian(ylim=c(0,1))
RESxCONS + facet_grid(structure~.) + theme(legend.position = "none") 


##################################################################
##           get trend line for abs(residuals)~conservation     ##
##################################################################
coefs <-coef(lm(abs(residual)~conservation, data= working))
line<- coefs[2]*working$conservation + coefs[1]

RESxCONS +   coord_cartesian(ylim=c(0,1)) +
  geom_smooth(color="black",aes(x=working$conservation, y=line))

#################################################################
##                    residual~conservation                    ##
##                        ( no abs() )                         ##
#################################################################

resXcons <- ggplot(data=working,
                   aes(x=conservation,y=residual, 
                       color=structure)) +
            geom_point() 
resXcons 

resXcons + facet_grid(structure~.) + theme(legend.position = "none")

##################################################################
##                   uncat. notable positions                   ##
##################################################################

nota1 <- as.data.frame(rbind(
working[which(working$position==298) , ],
working[which(working$position==419) , ],
working[which(working$position==295) , ],
working[which(working$position==310) , ],
working[which(working$position==321) , ],
working[which(working$position==263) , ],
working[which(working$position==315) , ]
))

ggplot(nota1, aes(x=conservation,y=residual, color=structure)) +
  geom_text(data=nota1, label=nota1$position) +
  ggtitle("uncat noted positions")

#################################################################
##              Gonfloni et al notable positions               ##
##             that impair regulation when mutated             ##
#################################################################

Gonfloni <- as.data.frame(rbind(
working[which(working$position==289) , ],
working[which(working$position==290) , ],
working[which(working$position==307) , ],
working[which(working$position==327) , ]
))

Golfoni_plot <- ggplot(Gonfloni, aes(x=conservation,y=residual, color=client_status)) +
  geom_text(data=Gonfloni, label=Gonfloni$position) +
  ggtitle("Residues identified by Gonfloni et al 
that 'strongly impair regulation when mutated' ")

Golfoni_plot



##################################################################
##               Gonfloni et al notable positions               ##
##          that DO NOT impair regulation when mutated          ##
##################################################################

GonfloniX <- as.data.frame(rbind(
  working[which(working$position==302) , ],
  working[which(working$position==321) , ],
  working[which(working$position==323) , ],
  working[which(working$position==333) , ],
  working[which(working$position==334) , ],
  working[which(working$position==362) , ],
  working[which(working$position==363) , ],
  working[which(working$position==365) , ],
  working[which(working$position==399) , ],
  working[which(working$position==400) , ],
  working[which(working$position==497) , ],
  working[which(working$position==488) , ]
))

GolfoniX_plot <- ggplot(GonfloniX, aes(x=conservation,y=residual, color=structure)) +
  geom_text(data=GonfloniX, label=GonfloniX$position) +
  ggtitle("Residues identified by Gonfloni et al 
that 'strongly impair regulation when mutated' ")

GolfoniX_plot




##################################################################
##                combine impact and no impact                  ##
##               Gonfloni et al notable positions               ##
##################################################################


grid.arrange(Golfoni_plot,GolfoniX_plot)


Gonfloni[,7] <- "impact"
GonfloniX[,7] <- "no_impact"


#combined data
GonfloniC <-  as.data.frame(rbind(Gonfloni,GonfloniX))
colnames(GonfloniC)[7] <- "impact"

colnames(GonfloniC)

GonfloniC_plot <- ggplot(GonfloniC, aes(x=conservation,y=residual, color=client_status)) +
  geom_text(data=GonfloniC, label=GonfloniC$position) +
  ggtitle("Gonfloni et al positions")

GonfloniC_plot 
GonfloniC_plot + facet_grid(.~impact)
GonfloniC


##################################################################
##                        define clients                        ##
##################################################################

#positive resid= greater radicicol score = nonclient?  

sd(working$residual) #stand. dev = 0.2476952
mean(working$residual) #mean = essentially 0. 



workingtemp <- working
workingtemp<- as.data.frame(
  cbind(workingtemp,abs(workingtemp$residual)))
colnames(workingtemp) <- c(colnames(workingtemp)[1:6],"absres")
colnames(workingtemp)


strongclients <- workingtemp[which(workingtemp$absres >sd(workingtemp$absres)*2 ) , ]
clients <- workingtemp[which(workingtemp$absres > sd(workingtemp$absres) & workingtemp$absres < sd(workingtemp$absres)*2 ) , ] 
weakclients <- workingtemp[which(workingtemp$absres < sd(workingtemp$absres) & workingtemp$absres > sd(workingtemp$absres)/2  ) , ] 
nonclients<- workingtemp[which(workingtemp$absres < sd(workingtemp$absres)/2 ) , ]
#strong = greater than 2x sd(abs(residual))
#client = greater than sd(abs(residual))
#weak = less than sd(abs(residual))
#non = less than 1/2 sd(abs(residual))

strongclients[,8] <- "strong_client"
clients[,8] <- "client"
weakclients[,8] <- "weak_client"
nonclients[,8] <- "nonclient"

range(strongclients$absres)
range(clients$absres)
range(weakclients$absres)
range(nonclients$absres)


colnames(working)[8] <- "client_status"
working <-as.data.frame(rbind(strongclients,clients,weakclients,nonclients))
#overwrite working DF file: 
#write.csv(working, "/Users/arielleaiken/github/src4/src/arielle/R_working_master_20200421.csv", row.names = FALSE)


##################################################################
##                        back to plots,                        ##
##         now with client status information available         ##
##################################################################

working$client_status

VCS <- ggplot(working, aes(x=conservation, y=residual, color=client_status)) +
  geom_point() + ggtitle("visualized client status")





#################################################################
##          Lindquist et al + Gonfloni et al analysis          ##
#################################################################

#Lindquist 1999 mentioned Lys298Met reduced src activity 
working[which(working$position==298) , ]
# Lys --> Met  res:0.006287953    cons:0.6086163 
#label this as a lack "L"

Gonfloni$position #289,290,307,327,
Gonfloni$conservation
#label Gonfloni positions as G

VCS + geom_text(label="L", x=0.6086163,y=0.006287953, color="black") +
  geom_text(label="G", x=0.7008886,y=0,color="black")+
  geom_text(label="G", x=0.3033911,y=0,color="black")+
  geom_text(label="G", x=0.2812683,y=0,color="black")+
  geom_text(label="G", x=0.3381344,y=0,color="black")
  
  




