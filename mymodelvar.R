
# N1 numbers of fish used to calculate variance on age
# fish for which we have 8 readings : 4 by readers
# we want to reattribute a probable age to each fish in the dataset N

#___________________________________________________________#

rm(list=ls())   # Clear memory



## WORKING DIRECTORY

work.dir<-paste('C:/Users/laulus.INRA/Desktop/THESE/ANALYSES/modele-complet')

setwd(work.dir)
library(rjags)
load.module("glm")

#--------------- DATA --------------#
library(readr)
var <- read_delim("vartest.csv",";", escape_double = FALSE, trim_ws = TRUE)
var<-var[, c("ID","River","Num_ec","Lecteur","Age_C","AgeT")]

str(var)

var<-na.omit(var)
View(var)

N1=nrow(var)/8  # 36 fish in vartest => ok

data<-list(N=N1, ID=var$ID, Reader=as.factor(var$Lecteur), river=as.factor(var$River), age=var$AgeT, idfish=unique(var$ID),idreader=unique(var$Lecteur), idriver=unique(var$River))
#______________________________________________________________________#

##-----------------------------MODEL ----------------------------------##

write("
      model {
      
      # likelihood
      
      ##  Multinomial logit link function         
      for (j in 1:n.obs) { 
      Age[j] ~dcat(PSI[ID[j], 1:8] )
      for (k in 1:n.age) { 
      PSI[j,k]  <- psi[j,k] / sum(psi[j,])   # constrain the transition such that their sum < 1 :   
      } # end loop k
      log(psi[j,9])<-0 #Reference category, set to zero
      for (k in 1 : (n.age+1)) { log(psi[j,k]) <- alpha[k] + eps[j]}    
      eps[j]<-  a[ID[j]]+ b[Reader[j]] +epsres[j]
      } # end loop j
      #for (i in 1:n.age) { 
      #Age.pred[i] ~dcat(PSI.pred[j,1:8] )
      #PSI.pred[j,k] <- alpha[k] + a[ID[j]]
      #}
      for (i in 1:N1)  {
      log(age[i])~dnorm(log(agetrue[i]),tauage)
      agetrue[i] <- age[i] * eps[i]
      eps[i]<-  a[ID[i]]+ b[Reader[i]] +epsres[i] # a is intra fish variance between 0 and 2 (8 scales)
      # b is inter readre variance between 0 and 2 # epsres variance not explained by intra nor inter
      }
      
      # priors
      for (i in idfish) {
      a[i]~dnorm(0,tauki)
      epsres[i]~dnorm(0,0.001) # residuals in age variance
      }
      for (i in idreader) {   
      b[i]~ dnorm(0,0.001)       
      }
      
      tauage~dgamma(0.001,0.001)
      tauki~dgamma(0.001,0.001)
      } # END OF MODEL,
      ", "mymodelvar.R")


#____________________________________________________________________#

# ANALYSIS

require(rjags)



nChains = 2 # Number of chains to run.

adaptSteps = 1000 # Number of steps to "tune" the samplers.

burnInSteps = 5000 # Number of steps to "burn-in" the samplers.

niter=50000 # Total number of steps in chains to save.





## PARAMETERS TO SAVE
parameters=c("a","b","c","tauage") 



#______________________________________________________________________#

## INITS

inits<-function(){
  # inits "random" method 
  list(a= , b= , # define inits for a b c and d to reestimate
}

#________________________________________________________________________#
## Compile & adapt

### Start of the run ###

## Compile & adapt

# Create, initialize, and adapt the model:

jagsfit=jags.model(
  
  'mymodel.R',
  
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

#save(data,jagsfit,fit.mcmc,file=paste("Results_",model.name,".RData",sep=""))


#______________________________________________________________________________#

#------------------------------------------------------------------------------

# EXAMINE THE RESULTS



summary(fit.mcmc)
gelman.diag(fit.mcmc,multivariate=F) # test stat de m?lange chaine mcmc, convergence? multivariate psrf <1.1


library(mcmcplots)
traplot(fit.mcmc,"a")
traplot(fit.mcmc,"b")
traplot(fit.mcmc,"c")
traplot(fit.mcmc,"d")
traplot(fit.mcmc,"eps.river")
traplot(fit.mcmc,"betar")
caterplot(fit.mcmc,c("beta0"),greek=TRUE, reorder=F)
caterplot(fit.mcmc,c("beta1"),greek=TRUE, reorder=F)
caterplot(fit.mcmc,c("eps.river"),greek=TRUE, reorder=F)
# caterplot(fit.mcmc,c("betar"),greek=TRUE, reorder=F)
#plot(fit.mcmc)

#------------ COEFFICIENTS --------------#
beta0.med<-summary(fit.mcmc)$quantiles[paste0("beta0[",1:10,"]"),"50%"] # on r?cup?re le vecteur des beta0
beta1.med<-summary(fit.mcmc)$quantiles[paste0("beta1[",1:10,"]"),"50%"] # on r?cup?re le vecteur des beta1
betar<-summary(fit.mcmc)$quantiles[paste0("betar[",1:3,"]"),"50%"]
eps.river<-summary(fit.mcmc)$quantiles[paste0("eps.river[",1:3,"]"),"50%"] # on r?cup?re le vecteur des eps.river

# deltaeps<-summary(fit.mcmc)$statistics[1:3,1:2]
#deltaeps<-summary(fit.mcmc)$statistics[40:42,1:2] # faux a modifier


#_______________________________________________________________#
#-------------PREDICT--------------#


#predictions
for (i in idfish) {
  for (j in idreader) {
    for (k in idriver) {
      agepred[i,j,k] ~ dnorm(agetruepred[i,j,k],tau)
      agetruepred[i,j,k] <-agetot[i] * (a[i]+b[i]+c[j]+d[k])
    }
  }   
}