#################################################################
##                        April 22 2020                        ##
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

strongclients <- read.csv("data/strongclients.csv")
clients <- read.csv("data/clients.csv")
weakclients<- read.csv("data/weakclients.csv")
nonclients <- read.csv("data/nonclients") 


##################################################################
##            define inhibited and dependent clients            ##
##################################################################

StrongInhibited <- filter(strongclients, residual > 0 ) 
StrongDependent <-  filter(strongclients, residual < 0 ) 
ModerateInhibited <- filter(clients, residual > 0 ) 
ModerateDependent <- filter(clients, residual < 0 ) 
WeakInhibited <- filter(weakclients, residual > 0  )
WeakDependent <- filter(weakclients, residual < 0  )

StrongInhibited[8] <- "Strong_Inhibited"
StrongDependent[8] <- "Strong_Dependent"
ModerateInhibited[8] <-"Moderate_Inhibited"
ModerateDependent[8] <- "Moderate_Dependent"
WeakInhibited[8] <- "Weak_Inhibited"
WeakDependent[8] <- "Weak_Dependent"

#store as a dataframe 'working' and save. 
working <- as.data.frame(rbind(StrongDependent,StrongInhibited,
                               ModerateDependent,ModerateInhibited,
                               WeakDependent,WeakInhibited,
                               nonclients))
working <- arrange(working, position)
colnames(working)[8] <- "client_status"

#write.csv(working, "/Users/arielleaiken/github/src4/src/arielle/R_working_master_20200422.csv", row.names = FALSE)


#################################################################
##                            plots                            ##
#################################################################
working$client_status <- as.factor(working$client_status)
levels(working$client_status)

working$client_status <- 
  factor(working$client_status, 
         levels = c( "Strong_Dependent", "Moderate_Dependent", 
                    "Weak_Dependent","nonclient", "Weak_Inhibited",
                    "Moderate_Inhibited", "Strong_Inhibited"
                    ))

working$client_status <- 
  factor(working$client_status, 
         levels = rev(levels(working$client_status)))

ggplot(working, aes(x=conservation, y=residual, color=client_status)) +
  geom_point() + ggtitle("client status visualzation 2") +
  scale_color_manual(values = c("red", "orange", "green",
                                "black","blue","purple","hotpink"))

rp<-ggplot(filter(working, client_status!="nonclient"),
            aes(x=position, y=residual, color=client_status)) +
  geom_point() + ggtitle("client status by position") +
  scale_color_manual(values = c("red", "orange", "yellow",
                               "blue","purple","hotpink"))
rp
cp<-ggplot(
    filter(working, client_status!="nonclient"), 
      aes(x=position, y=conservation, color=client_status)) +
  geom_point() + ggtitle("conservation by position") +
  scale_color_manual(values = c("red", "orange", "yellow",
                                "blue","purple","hotpink"))
cp

cr <- ggplot(
  filter(working, client_status!="nonclient"), 
  aes(x=conservation, y=residual, color=client_status)) +
  geom_point() + ggtitle("client status conservation") +
  scale_color_manual(values = c("red", "orange", "yellow",
                                "blue","purple","hotpink"))
cr
grid.arrange(rp,cr, nrow=1)

ggplot(working, aes(x=position, y=conservation, color=structure)) +
  geom_point(alpha=.4)  + ggtitle("client status visualzation 2.2") +
  facet_grid(client_status~.) 


strong <- filter(working,client_status=="Strong_Inhibited" | client_status=="Strong_Dependent")

ggplot(strong, aes(x=position, y=conservation)) +
  geom_point(alpha=.4)  + ggtitle("client status visualzation 2.2") +
  facet_wrap(client_status~.) 


ggplot(strong, aes(x=conservation, fill=client_status)) + 
  geom_histogram() +
  facet_grid(client_status~.) + theme(legend.position = "none")

arrange(strong, desc(conservation))


strongnon <- filter(working,client_status=="Strong_Inhibited" | 
                    client_status=="Strong_Dependent" | 
                    client_status=="nonclient")


ggplot(strongnon, aes(x=conservation, fill=client_status)) + 
  geom_histogram() +
  facet_grid(client_status~.) + theme(legend.position = "none")


strongsamps <- as.data.frame(rbind(
                  StrongDependent %>% sample_n(100),
                  StrongInhibited %>% sample_n(100)))


p <- ggplot(strongsamps, aes(x=conservation, fill=V8)) + 
  geom_histogram() +
  facet_grid(V8~.) + theme(legend.position = "none")
p




