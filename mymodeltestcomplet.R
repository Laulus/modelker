
## WORKING DIRECTORY

# work.dir<-paste('D:/THESE/ANALYSES/modele-complet/modelker')

work.dir<-paste('C:/Users/laulus/Desktop/THESE/ANALYSES/modele-complet/modelker')

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

data<-load("dataframetestmodelker.RData")
View(dataframe)

data<-list(N=nrow(dataframe),ID=as.factor(unique(dataframe$ID)),idfish=unique(dataframe$ID), len=unique(dataframe$len),length=dataframe$len, agei=dataframe$agei, radi=dataframe$radi, aget=dataframe$aget, agem=dataframe$agem, pheno=dataframe$pheno, RT=unique(dataframe$RT), rt=dataframe$RT, RM=dataframe$RM, milieu=dataframe$milieu)
#______________________________________________________________________#

##-----------------------------MODEL ----------------------------------##

write(" model {
  
  # likelihood Backcalculation
   for (i in 1:50){
     len[i]~dnorm(mu.tot[i], tau.tot)
     mu.tot[i]<-a.retro + c.retro* RT[i]^b.retro
  # len<-a.retro + c.retro* exp(b.retro*log(RT))
# mu.tot[i]<-a.retro + b.retro*RT[i]
  } # end of loop i
  
  # priors
  a.retro~dnorm(100,0.01) #T(0,)
  b.retro~dnorm(1,0.01)
  c.retro~dnorm(0.1,0.01)
 
  tau.tot~dgamma(0.001,0.001)
# sigma.tot<-1/tau.tot
  
  #predictions
  # WARNING RM=RT for sedentary individuals !!!!!!!!!!!!!!!!!

## DANS CETTE BOUCLE LA IL FAUT RAJOUTER DE LA VARIANCE SUR LES RADI,
# et ajouter de la variance sur les RT
# admettons qui augmente avec l'age j
# on aura donc une len retro i j moyenne et des intervalles de confiance
  
  for (i in 1:N) { # for each lines
        len.retro[i] <- a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*((radi[i]/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]
  }

  for (i in 1:N)  { #for each lines
    len.retro[i] ~ dnorm(mu[i],tau)
    mu[i] <- Linf*(1-exp(-k*(agei[i]-t0)))
    # Linf[i]<- a[river[i]] + c[ID[i]] # effet Linf par riviere et aussi effet random fish sur Linf
    # Linf[i]<- c[ID[i]] # on retire l'effet rivi?re sur Linf
    # k[i]<-b[river[i]] + d[ID[i]]
    # k[i]<-d[ID[i]]
    
  }

# # priors
#   for (i in idfish) {
#     # c[i]~dnorm(0,tauki)
#     d[i]~dnorm(0,tauli)       
#   }

k~dnorm(0,0.01)
Linf~dnorm(1000, 0.001)
tau~dgamma(0.001,0.001)
tauki~dgamma(0.001,0.001)
tauli~dgamma(0.001,0.001)
t0~ dunif(-10,0.5)

} # END OF MODEL
  ", "mymodelbackcalculation.R")



#____________________________________________________________________#

## -------------------------- ANALYSIS ------------------------------#

require(rjags)



nChains = 2 # Number of chains to run.

adaptSteps = 1000 # Number of steps to "tune" the samplers.

burnInSteps = 5000 # Number of steps to "burn-in" the samplers.

niter=50000 # Total number of steps in chains to save.





## PARAMETERS TO SAVE
parameters=c("a.retro","b.retro","c.retro","tau.tot",      #retrocalcul
             "Linf","k","tau","tauki","tauli","t0")
# parameters=c("a.retro","b.retro","tau.tot") 

#______________________________________________________________________#

## ------------------------------INITS---------------------------------#

inits<-function(){
  # inits "random" method 
  list(a.retro=100 , b.retro=1 , c.retro=0.1, tau.tot =0.01,
       Linf=600, t0=-0.2, k= 1, tau=0.1) #,tauki=0.3, tauli=0.3) # define inits for a b c and d to reestimate
  # list(a.retro=100 , b.retro=0.1 , tau.tot =0.1)
}

#________________________________________________________________________#
## Compile & adapt

### Start of the run ###

## Compile & adapt

# Create, initialize, and adapt the model:

jagsfit=jags.model(
  
  'mymodelbackcalculation.R',
  
  data,inits,
  
  n.chains=nChains,
  
  n.adapt = adaptSteps)



# Burn-in:

cat( "Burning in the MCMC chain.\n" )

#jagsfit$recompile()

update(jagsfit, n.iter=burnInSteps,progress.bar="text")



# The saved MCMC chain:

cat( "Sampling final MCMC chain.\n" )

fit.mcmc<-coda.samples(jagsfit,variable.names=parameters, n.iter=niter, thin=1)



## BACKUP

save(data,jagsfit,fit.mcmc,file=paste("Results_","mymodelbackcalculation",".RData",sep=""))


#______________________________________________________________________________#


# EXAMINE THE RESULTS



summary(fit.mcmc)
gelman.diag(fit.mcmc,multivariate=F) # test stat de m?lange chaine mcmc, convergence? multivariate psrf <1.1


library(mcmcplots)
traplot(fit.mcmc,"a.retro")
traplot(fit.mcmc,"b.retro")
traplot(fit.mcmc,"c.retro")
traplot(fit.mcmc, "t0")
traplot(fit.mcmc, "Linf")


