
###########################################################################
###########################################################################
###                                                                     ###
###                          APRIL 29 2020 AA                           ###
###                        CONSERVATION ANALYSIS                        ###
###                              A. AIKEN                               ###
###                                                                     ###
###########################################################################
###########################################################################



library(ggseqlogo)
library(bio3d)
library(ggplot2)
library(ggExtra)
library(seqinr)
library(grDevices)
library(readr)
library(gridExtra)
library(dplyr)
library(Bios2cor)  

working <- read.csv("data/R_working_master_20200422.csv")
src_fasta <- import.fasta("data/src_fasta.fasta")
pka_fasta <- import.fasta("data/pka_fasta.fasta")

pka_fasta$`sp|P17612|KAPCA_HUMAN`
src_fasta$`sp|P12931|SRC_HUMAN`



##################################################################
##                        strong clients                        ##
##################################################################

a(filter(working, position==274, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==274, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==302, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==302, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==329, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==329, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==363, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==363, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==367, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==367, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==369, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==369, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==373, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==373, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==375, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==375, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==379, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==379, client_status=="Strong_Inhibited")$mutation)

a(filter(working, position==456, client_status=="Strong_Dependent")$mutation)
a(filter(working, position==456, client_status=="Strong_Inhibited")$mutation)


strong <- c("NKGHH*S*QN", 
            "EYR******C", 
            "G********H", 
            "K********L",
            "P********M", 
            "S*********", 
            "W*********",
            "G*********", 
            "K*********",
            "P*********", 
            "S*********", 
            "W*********") 



strong
StrongLogo<-ggseqlogo(strong, method="prob") +
  xlab("position in c-src:
274  302  329  363  367  369  373  375  379  456
positions in PKA:
48   77   109  141  145  147  151  153  157  230") +
ggtitle("Strong clients")
StrongLogo


##################################################################
##                         weak clients                         ##
##################################################################

a(filter(working, position==274, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==274, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==302, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==302, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==329, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==329, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==363, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==363, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==367, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==367, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==369, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==369, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==373, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==373, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==375, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==375, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==379, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==379, client_status=="Weak_Inhibited")$mutation)

a(filter(working, position==456, client_status=="Weak_Dependent")$mutation)
a(filter(working, position==456, client_status=="Weak_Inhibited")$mutation)


weak <- c("IINAQFCFRD",
          "*AWRSPEGDI",
          "*DCE*TA*K*",
          "*M****R***",
          "*F****K***", 
          "******V***")

WeakLogo<-ggseqlogo(weak, method="prob") +
  xlab("position in c-src:
274  302  329  363  367  369  373  375  379  456
positions in PKA:
48   77   109  141  145  147  151  153  157  230") +
  ggtitle("Weak clients")
WeakLogo


#################################################################
##                         non clients                         ##
#################################################################

a(filter(working, position==274, client_status=="nonclient")$mutation)
a(filter(working, position==302, client_status=="nonclient")$mutation)
a(filter(working, position==329, client_status=="nonclient")$mutation)
a(filter(working, position==363, client_status=="nonclient")$mutation)
a(filter(working, position==367, client_status=="nonclient")$mutation)
a(filter(working, position==369, client_status=="nonclient")$mutation)
a(filter(working, position==373, client_status=="nonclient")$mutation)
a(filter(working, position==375, client_status=="nonclient")$mutation)
a(filter(working, position==379, client_status=="nonclient")$mutation)
a(filter(working, position==456, client_status=="nonclient")$mutation)


nonclient <- c("MCANNCDQNR",
               "*EDDD*LELQ",
               "*GEQE*PHME",
               "*LHGG**TFG",
               "*SFII***PK",
               "*TPKP***SP",
               "*WSM****TV",
               "*VTF****W*",
               "***P****V*",
               "***S******",
               "***T******",
               "***U******",
               "***V******",
               "VPYLVMISYT") #wildtype is also a nonclient

src_fasta$`sp|P12931|SRC_HUMAN`[c(274,302,329,363,367,369,373,375,379,456)] #wildtype 


nonclientLogo<-ggseqlogo(nonclient, method="prob") +
  xlab("position in c-src:
274  302  329  363  367  369  373  375  379  456
positions in PKA:
48   77   109  141  145  147  151  153  157  230") +
  ggtitle("Nonclients")
nonclientLogo

grid.arrange(StrongLogo,WeakLogo,nonclientLogo)



filter(as.data.frame(filter(working, position==274)$mutation), #7/13
        filter(working, position == 274)$mutation!="Stp")
7/13 #0.53

filter(as.data.frame(filter(working, position==302)$mutation), 
       filter(working, position == 302)$mutation!="Stp") #17
2/17 #0.12

filter(as.data.frame(filter(working, position==329)$mutation), 
       filter(working, position == 329)$mutation!="Stp") #17
2/17#0.12

filter(as.data.frame(filter(working, position==363)$mutation), 
       filter(working, position == 363)$mutation!="Stp") #18
1/18 #0.06

filter(as.data.frame(filter(working, position==367)$mutation), 
       filter(working, position == 367)$mutation!="Stp") #11
1/11 #0.09

filter(as.data.frame(filter(working, position==369)$mutation), 
       filter(working, position == 369)$mutation!="Stp") #6
0/6 #0

filter(as.data.frame(filter(working, position==373)$mutation), 
       filter(working, position == 373)$mutation!="Stp") #11
1/11 #0.09

filter(as.data.frame(filter(working, position==375)$mutation), 
       filter(working, position == 375)$mutation!="Stp") #15
0/15 #0

filter(as.data.frame(filter(working, position==379)$mutation), 
       filter(working, position == 379)$mutation!="Stp") #18
1/18 #0.06

filter(as.data.frame(filter(working, position==456)$mutation), 
       filter(working, position == 456)$mutation!="Stp") #16
5/16 #0.31











