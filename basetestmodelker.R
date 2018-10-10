library(FSA)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(itsadug) 
library(ggplot2)
library(nlme)
library(nls2)
library(knitr)

# GENERER UN JEU DE DONNEES 

ID<-c(1:50)
aget<-sample(c(1:12),50,replace=T)
rep<-rep(0,25)
samp<-sample(c(1:4),25,replace=T)
agem<-c(rep,samp)


a<-c(1,2,3,4,5,6,7,8,9,10,11,12)
mulen<-c(60, 100, 220, 260, 300, 340, 380, 405,430, 460, 485, 510)
taulen<-rep(20,12)

len<-rep(0,50)

for (i in 1:50){ # for each fish
  j<-aget[i]
  len[i]<-rnorm(1,mulen[j],taulen[j])
}

dataframe<-data.frame(ID,len,aget,agem)
RT<-rep(0,50)
RM<-rep(0,50)

# We have 50% TS and 50% TM at that point





      a1<-rnorm(50, 210, 50)
      a2<- rnorm(50, 490, 100)
      a3<-rnorm(50, 790, 100)
      a4<-rnorm(50, 1240, 100)
      a5<-rnorm(50, 1645, 100)
      a6<-rnorm(50, 2020,100)
      a7<-rnorm(50, 2350,100)
      a8<-rnorm(50, 2680,100)
      a9<-rnorm(50, 2910,100)
      a10<-rnorm(50, 3000,50)
      a11<-rnorm(50, 3100,50)
      a12<-rnorm(50, 3200, 50)
      
      dataframe<-cbind(dataframe, pheno, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
        
for (i in 1:length(dataframe$ID)){
  alpha<-5+as.numeric(dataframe$aget[i])

  RT[i]<-dataframe[i,alpha]
  
  if (agem[i]!=0){
      pheno[i]<-"TM"
      beta<-4+as.numeric(dataframe$agem[i])
    RM[i]<-dataframe[i,beta]}
  else 
    {pheno[i]="TS"
    RM[i]<-RT[i]}
    

    } # end of or j

  dataframe<-cbind(dataframe, RT, RM)

    
  # We convert the line matrix into a column matrix #
  datA<-gConvert(dataframe,in.pre="a",out.type="rad")
  dataframe <- gather(dataframe,agei,radi,a1:a12) 
  
  str_sub(dataframe$agei,start=1,end=1) <- "" # suppression of the 3 first character ne the column agei 
  #var$AgeT<-as.numeric(var$AgeT)
  dataframe %<>% mutate(agei=as.numeric(agei))%<>%
    filterD(!is.na(radi)) %>% # only the fish with growth increments are studied
    filterD(agei<=aget)%>%
    arrange(ID, agei)
  
  View(dataframe)
  
  # ADDING Environment of growth : 0 freshwater, 1 salt water
  
  milieu<-rep(0,length(dataframe$ID))
  for (i in 1:length(dataframe$ID)){
    if (dataframe$pheno[i]=="TS") {milieu[i] <- 0}
    else 
      if (dataframe$radi[i]<dataframe$RM[i]) {milieu[i]<-0}
      else {milieu[i]<-1}
  }

  dataframe<-cbind(dataframe,milieu)
  
  getwd() # "D:/THESE/ANALYSES/modele-complet/modelker"
  
  save(dataframe,file="dataframetestmodelker.RData")
  
# 