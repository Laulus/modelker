
## WORKING DIRECTORY

# work.dir<-paste('D:/THESE/ANALYSES/modele-complet/modelker')

work.dir<-paste('C:/Users/laulus/Desktop/THESE/ANALYSES/modele-complet/modelker')

setwd(work.dir)
library(rjags)
load.module("glm")

#--------------- DATA --------------#
library(readr)

# getwd()
data<-load("base1.RData") # projettut
# load("base2.RData") # projettut2

# data<-load("dataframetestmodelker.RData") 
# load the base you want to work on, dataframetestmodelker is a test
# base1 is old measures made by Eddy, Marc and others
# base 2 are recent readings made with the newest scalimetry methodology
# View(dataframe)
View(projettut)
summary(projettut)
# View(projettut2)


# for ( i in 1:length(projettut2$age.mer)){
#   if (is.na(projettut2$age.mer[i])==TRUE)
#   projettut2$age.mer[i]<-0
# }
# projettut2$age.mer<-as.numeric(projettut2$age.mer)


# ----AJOUT DE LA COLONNE MILIEU DANS LES BASES 1 et 2 (projettut et projettut2)

# projettut2<-subset(projettut2,!is.na(projettut2$RM)) # pour l'instant je retire les nA dans la colonne RM il faudra corriger
# projettut2$milieu<-rep(0,nrow(projettut2))
# for (i in 1:nrow(projettut2)){
# if ( projettut2$pheno[i]=="TM" && (projettut2$radi[i]>=projettut2$RM[i]) )
#  projettut2$milieu[i]<-1
# }
# 
# projettut$milieu<-rep(0,nrow(projettut))
# for (i in 1:nrow(projettut)){
#   if ( projettut$pheno[i]=="TM" && (projettut$radi[i]>=projettut$RM[i]) )
#     projettut$milieu[i]<-1
# }

# projettut2<-subset(projettut2,!is.na(projettut2$lf)) # pour l'instant je retire les nA dans la colonne RM il faudra corriger
# projettut2<-subset(projettut2,!is.na(projettut2$RT)) 
# 
# # dataframe<-dataframe
dataframe<-projettut
len<-rep(0,nlevels(as.factor(dataframe$id)))
RT<-rep(0,nlevels(as.factor(dataframe$id)))

for (i in 1:nlevels(as.factor(dataframe$id))){
  len[i]<-unique(dataframe$lf[which(dataframe$id==levels(as.factor(dataframe$id))[i])])
}

for (i in 1:nlevels(as.factor(dataframe$id))){
  RT[i]<-unique(dataframe$RT[which(dataframe$id==levels(as.factor(dataframe$id))[i])])
}
# dataframe<-projettut2
# colnames(dataframe)[11]<-"aget"
# colnames(dataframe)[12]<-"len"

# rajouter les différents autres facteurs dans le modèle , ie effet rivière, lecteur, cohorte ...

data<-list(N=nrow(dataframe),ID=as.factor(unique(dataframe$id)),idfish=unique(dataframe$id), len=len,length=dataframe$lf, agei=dataframe$agei, radi=dataframe$radi, aget=dataframe$age.tot, AGE=unique(dataframe$age.tot), pheno=dataframe$pheno, RT=RT, rt=dataframe$RT, RM=dataframe$RM, milieu=dataframe$milieu)
#______________________________________________________________________#

##-----------------------------MODEL ----------------------------------##

write("


model {

  # likelihood Backcalculation
   for (i in 1:length(idfish)){
     len[i]~dnorm(mu.tot[i], tau.tot)
     mu.tot[i]<-a.retro + c.retro* RT[i]^b.retro
  } # end of loop i

  # priors
  a.retro~dnorm(30,0.01) #T(0,)
  b.retro~dnorm(1,0.01)
  c.retro~dnorm(3,0.01)


# tau.tot <- 1/sqrt(sigma.tot)
# sigma.tot~dunif(0,10)

  tau.tot~dgamma(0.001,0.001)

  #predictions
  # WARNING RM=RT for sedentary individuals !!!!!!!!!!!!!!!!!

## DANS CETTE BOUCLE LA IL FAUT RAJOUTER DE LA VARIANCE SUR LES RADI,
# et ajouter de la variance sur les RT
# admettons qui augmente avec l'age j
# on aura donc une len retro i j moyenne et des intervalles de confiance



   for (i in 1:N) { # for each lines
         # len.retro[i] <- a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*(((radi[i]+eps[aget[i]])/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]
len.retro[i]~dnorm(mu.retro[i], tau.retro)
mu.retro[i]<-a.retro + ((length[i]*RM[i]/rt[i])-a.retro)*(((radi[i])/RM[i])^b.retro) # attention a la matrice de depart pour les r[i,j]  
   }

# tau.retro <- 1/sqrt(sigma.retro)
# sigma.retro~dunif(0,10)
tau.retro~dgamma(0.001,0.001)

# for (j in 1:max(aget[])){
# eps[j]~dnorm(0, 0.001)
# }


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
parameters=c("a.retro","b.retro","c.retro","tau.tot"
             ,"tau.retro"# "eps" ,      #retrocalcul
              )
# parameters=c("a.retro","b.retro","tau.tot") 

#______________________________________________________________________#

## ------------------------------INITS---------------------------------#

inits<-function(){
  # inits "random" method 
  list(a.retro=30 , b.retro=1 , c.retro=3, tau.tot =1, # eps=seq(0,100, length=max(data$aget)),
        tau.retro=1) 
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
# traplot(fit.mcmc, "t0")
# traplot(fit.mcmc, "Linf")


