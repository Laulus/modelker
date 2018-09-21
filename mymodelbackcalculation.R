
## WORKING DIRECTORY

work.dir<-paste('H:/THESE/ANALYSES/modele-complet/modelker')

setwd(work.dir)
library(rjags)
load.module("glm")

#--------------- DATA --------------#
library(readr)
# var <- read_delim("vartest.csv",";", escape_double = FALSE, trim_ws = TRUE)
# var<-var[, c("ID","Riviere","Num_ec","Lecteur","Age_C","AgeT")]
# 
# str(var)
# 
# var<-na.omit(var)
# View(var)
# 
# N1=nrow(var)/8  # 60 fish in vartest => ok
# n.obs=nrow(var)


# data<-list(N=N1, n.obs=n.obs, ID=var$ID, Reader=as.factor(var$Lecteur), river=as.factor(var$River), n.age=max(var$AgeT), age=var$AgeT, idfish=unique(var$ID),idreader=unique(var$Lecteur), idriver=unique(var$River))
#______________________________________________________________________#

##-----------------------------MODEL ----------------------------------##

write(" model {
  
  # likelihood Backcalculation
  for (i in 1:N){
    len.total[i]~dnorm(mu.tot[i], tau.tot)
    mu.tot[i]<-a.retro + c.retro* exp(b.retro*log(len.scale[i]))
  } # end of loop i 
  
  # priors
  a.retro~dnorm(0,0.001) T(0,)
  b.retro~dnorm(0,0.001)
  c.retro~dnorm(0,0.001)
 
  tau.tot~dgamma(0.001,0.001)
  
  #predictions
  
  # WARNING RM=RT for sedentary individuals !!!!!!!!!!!!!!!!!
  
  for (i in 1:N) { # for each fish
    for (j in 1:age[i]) { # for each age, avec age[i] = total age of fish i
        # len[i,j] ~ dnorm(mu.t[i,j],tau.tot)
        # mupred[i,j,h] <- a[i]*(1-exp(-b[j]*(h-t0)))
        len[i,j] <- a.retro + len.total[i]*((rm[i]/rt[i])-a)*exp(b*log(r[i,j]/rm[i])) # attention à la matrice de départ pour les r[i,j]
      }
  }   # END OF MODEL,
  ", "mymodelbackcalculation.R")