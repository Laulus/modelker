## Not run: 
## remove (almost) everything in the working environment.
## You will get no warning, so don't do this unless you are really sure.
rm(list = ls())
## WORKING DIRECTORY

# work.dir<-paste('D:/THESE/ANALYSES/modele-complet/modelker')

#work.dir<-paste('C:/Users/laulus/Documents/Mirror/THESE/MODELE-COMPLET/aulus/model-new')

#setwd(work.dir)
library(rjags)

#--------------- DATA --------------#
library(readr)

#data<-read.jagsdata("data/data.txt")
#data<-read.bugsdata("data/data.txt")
data<-read.bugsdata("data/data500-1.txt") # data recent
#data<-read.bugsdata("data/data500-1.txt") # data all 
#data<-read.bugsdata("data/data500-2.txt") # data Eddy


# Compute vector with occasions of first capture
get.last <- function(x) max(which(!is.na(x)))
last <- apply(data$radius, 1, get.last) #

# # Max age in freshwater
# agefr=NULL
# for (i in 1:dim(data$radius)[1]) {
#   agefr[i] <- ifelse(last[i]>data$smoltAge[i],data$smoltAge[i], last[i] )
# }

# 
# dataToJags <- list(nScale=data$nScale
#                    ,age=data$age
#                    ,smoltAge=data$smoltAge
#                    ,bodySize=data$bodySize
#                    ,radius=data$radius
#                    ,fishMaxSize=data$fishMaxSize
#                    ,scaleMaxSize=data$scaleMaxSize
#                    ,smolt=data$smolt
# )

# on retire les individus capturés entre aout et novembre / car dernier annuli non marqué / sosu-estimation âge
toremove <- which(data$age < data$age.tot)
dataToJags <- list(nScale=data$nScale - length(toremove)
                  ,age=data$age[-toremove]
                 ,smoltAge=data$smoltAge[-toremove]
                 ,bodySize=data$bodySize[-toremove,]
                 ,radius=data$radius[-toremove,]
                 ,fishMaxSize=data$fishMaxSize[-toremove]
                 ,scaleMaxSize=data$scaleMaxSize[-toremove]
                 ,smolt=data$smolt[-toremove,]
                 )

# dataToJags <- list(nScale=data$nScale
#                    ,age=last#agefr #data$age
#                    ,smoltAge=data$smoltAge
#                    ,bodySize=data$bodySize
#                    ,radius=data$radius
#                    ,fishMaxSize=data$fishMaxSize
#                    ,scaleMaxSize=data$scaleMaxSize
# )


plot(NULL,xlim=c(0,max(data$radius,na.rm=T)),ylim=c(0,max(data$bodySize, na.rm=T)),xlab="radius",ylab="size")
points(data$radius,data$bodySize)


##-----------------------------MODEL ----------------------------------##

write("
      model {
      
     for (i in 1:nScale) { # 20 scales belonging to 10 individuals, unbalanced. 
     
           Y[i]<- ifelse(age[i]>smoltAge[i],2,1) # 2: at sea / 1: in freshwater
           
		for (j in 1:age[i]) { # maxAge of each scale
      
      ### likelihood
      # 1. retrocalcul
      bodySize[i,j] ~ dlnorm(allo[i,j],tau[1]) #T(,fishMaxSize[i]+1)
      allo[i,j] <- a[Y[i]] + X[i,j,1]* b[Y[i]] * (fishMaxSize[i]*X[i,j,2] -a[Y[i]])
      
      X[i,j,1] <- radius[i,j]/radius[i,smoltAge[i]]
      X[i,j,2] <- radius[i,smoltAge[i]]/scaleMaxSize[i]
      #X[i,j,1] <- ifelse(age[i]>smoltAge[i], radius[i,j]/radius[i,smoltAge[i]], radius[i,j]/radius[i,age[i]])
      #X[i,j,2] <- ifelse(age[i]>smoltAge[i], radius[i,smoltAge[i]]/scaleMaxSize[i], 1)

      } # END OF LOOP j
      } # END OF LOOP i
      
      ### priors
      #a[1]~dnorm(0, 0.01) 
      a[1]~dunif(0, 10) 
      #a[2]~dunif(0, 10)
      a[2] <- a[1] #+ delta[1]
      delta[1]~dnorm(0, 0.01) 

      #b[1]~dnorm(0, 0.01) 
      b[1]~dunif(0, 10) 
      b[2]~dunif(0, 10) 
      #b[2] <- b[1] #+ delta[2]
      delta[2]~dnorm(0, 0.01) 
      
      tau[1]<-1/(sigma[1]^2)
      sigma[1]~dunif(0,10)
      
      } # END OF MODEL
      
      ", "code/mymodel2-jags.R")

#____________________________________________________________________#

## -------------------------- ANALYSIS ------------------------------#

require(rjags)

nChains = 2 # Number of chains to run.
adaptSteps = 1000 # Number of steps to "tune" the samplers.
burnInSteps = 5000 # Number of steps to "burn-in" the samplers.
niter=10000 # Total number of steps in chains to save.


## PARAMETERS TO SAVE
parameters=c("a","b","sigma","delta")

#______________________________________________________________________#

## ------------------------------INITS---------------------------------#

inits<-function(){
  list(sigma=c(0.1)
       , a=c(1,NA)
       , b=c(0.1,NA)
       ,delta=c(0,0)
       )}


#________________________________________________________________________#
## Compile & adapt
### Start of the run ###
## Compile & adapt
# Create, initialize, and adapt the model:

jagsfit=jags.model(
  'code/mymodel2-jags.R',
  dataToJags,inits,
  n.chains=nChains,
  n.adapt = adaptSteps)

# Burn-in:
cat( "Burning in the MCMC chain.\n" )
#jagsfit$recompile()
update(jagsfit, n.iter=burnInSteps,progress.bar="text")

# The saved MCMC chain:
cat( "Sampling final MCMC chain.\n" )
fit.mcmc<-coda.samples(jagsfit,variable.names=parameters, n.iter=niter, thin=2)
# To save individuals len:
# len.mcmc <- coda.samples(jagsfit,
# variable.names = c("mu.retro"),
# n.iter = 1,thin = 1)

# max<-data$maxAge
size.mcmc<-coda.samples(jagsfit,variable.names=c("bodySize[1:20,1:12]"), n.iter=niter, thin=2)

## BACKUP
save(fit.mcmc,file=paste0("results/Results_data-1","model2-jags",".RData"))


#______________________________________________________________________________#
# EXAMINE THE RESULTS

#summary(fit.mcmc)
# # summary(len.mcmc)
#gelman.diag(fit.mcmc,multivariate=F) # test stat de m?lange chaine mcmc, convergence? multivariate psrf <1.1
# 
# 
library(mcmcplots)
# traplot(fit.mcmc,"a")
# traplot(fit.mcmc,"b")
# traplot(fit.mcmc,"sigmae")

traplot(fit.mcmc,parameters)
denplot(fit.mcmc,parameters)
caterplot(fit.mcmc,"a")
caterplot(fit.mcmc,"b")

mcmc<-as.matrix(fit.mcmc)
plot(mcmc[,"a[1]"],mcmc[,"b[1]"])
plot(mcmc[,"a"],mcmc[,"sigmae"])
plot(mcmc[,"b"],mcmc[,"sigmae"])



mcmc<-as.matrix(size.mcmc)
#size <- apply(mcmc, 2,quantile, probs=0.5)


plot(NULL,xlim=c(0,dim(data$bodySize)[2]),ylim=c(0,1000), xlab="Age",ylab="Body size")

last <- apply(dataToJags$bodySize, 1, get.last) #
for (i in 1:dataToJags$nScale){
points(last[i], data$bodySize[i,last[i]], col="lightgrey")
}

for (i in 1:20){
  
  points(dataToJags$smoltAge[i],dataToJags$bodySize[i,dataToJags$smoltAge[i]],pch=17,col=i)
  points(dataToJags$age[i],dataToJags$bodySize[i,dataToJags$age[i]],pch=16,col=i)
  points(1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"50%"],col=i, type="b")
  
  segments(1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"2.5%"],1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"97.5%"],col=i)

  }


smp <- sample(1:niter,1)
col <- c("steelblue1","royalblue3")
for (i in 1:20){
  
  points(dataToJags$smoltAge[i],dataToJags$bodySize[i,dataToJags$smoltAge[i]],pch=17,col=col[data$])
  points(dataToJags$age[i],dataToJags$bodySize[i,dataToJags$age[i]],pch=16,col=i)
  points(1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"50%"],col=i, type="b")
  
  segments(1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"2.5%"],1:dataToJags$age[i],summary(size.mcmc)[[2]][paste0("bodySize[",i,",",1:dataToJags$age[i],"]"),"97.5%"],col=i)
  
}



# # RECONSTRUCTION DES mu_retro # solution temporaire, pas propre
# 
# a.retro<-summary(fit.mcmc)$quantiles[paste0("a.retro"),"50%"]
# b.retro<-summary(fit.mcmc)$quantiles[paste0("b.retro"),"50%"]
# 
# # mu.retro<-(summary(len.mcmc)$empirical[paste0("mu.retro"),"Mean"])
# 
# 
# mu.retro<-rep(0,length(dataframe$lf))
# for ( i in 1:length(dataframe$lf)){
#   mu.retro[i]<-a.retro + ((dataframe$lf[i]*dataframe$RM[i]/dataframe$RT[i])-a.retro)*(((dataframe$radi[i])/dataframe$RM[i])^b.retro)
# }
# dataframe<-cbind(dataframe,mu.retro)
# save(dataframe,file="dataframeeddy.RData")
